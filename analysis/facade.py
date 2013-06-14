#==============================================================================
# facade.py
# ..
#==============================================================================
import __main__
import atexit
from multiprocessing import Manager, Queue
from plotting import PlottingProcess
from seren.analysis.SimBuffer import SimBuffer, BufferException

manager= Manager()

#TODO: in all the Python code, raise proper exceptions rather than a generic Exception
#TODO: the tests should not fail


#------------------------------------------------------------------------------
class Singletons:
    '''Container class for singletons object. They are:
        queue
            Queue for sending commands to the plotting process
        commands
            List of the commands shared with the plotting process.
            Caution: if you modify a command, you need to reassign it in the list
            to make the changes propagate to the other process
        completedqueue
            Queue used from the plotting process to signal the completion
            of a command
        globallimits
            Dict that for each quantity gives the limits
    '''
    queue = Queue()
    commands = manager.list()
    completedqueue = Queue()
    globallimits = manager.dict()

import commandsource as Commands
from data_fetcher import CreateUserQuantity, CreateTimeData
import signal
from time import sleep

#figure out if we are in interactive mode
try:
    __main__.__file__
    interactive=False
except AttributeError:
    interactive=True


#TODO: add function for resizing (programmatically) the figure
#------------------------------------------------------------------------------
def handle(e):
    '''This functions takes care of printing information about an error,
    if we are in interactive mode, or re-raising it, if we are in script mode
    (so that the execution of the script can stop if nobody catches the exception)
    '''
    if interactive:
        print str(e)
    else:
        raise e


#------------------------------------------------------------------------------
def loadsim(run_id, fileformat = 'ascii', buffer_flag = 'cache'):
    '''Given the run_id of a simulation, reads it from the disk'''
    SimBuffer.loadsim(run_id, fileformat=fileformat, buffer_flag=buffer_flag)
    return SimBuffer.get_current_sim()
    
def plot(x,y, type="default", snap="current", sim="current", overplot = False, autoscale = True, xunit="default", yunit="default", **kwargs):
    '''Plot particle data as a scatter plot.  Creates a new plotting window if
one does not already exist.

Required arguments:
    x          : Quantity on the x-axis. Must be a string.
    y          : Quantity on the y-axis. Must be a string.
        
Optional arguments:
    snap       : Number of the snapshot to plot. Defaults to 'current'.       
    sim        : Number of the simulation to plot. Defaults to 'current'.    
    overplot   : If True, overplots on the previous existing plot rather
                 than deleting it. Defaults to False.
    autoscale  : If True (default), the limits of the plot are set
                 automatically.  Can also be set to 'x' or 'y' to specify
                 that only one of the axis has to use autoscaling.
                 If False, autoscaling is not used. On an axis that does
                 not have autoscaling turned on, global limits are used
                 if defined for the plotted quantity.
    xunit      : Specify the unit to use for the plotting for the quantity
                 on the x-axis.
    yunit      : Specify the unit to use for the plotting for the quantity
                 on the y-axis.
    **kwargs   : Extra keyword arguments will be passed to the matplotlib 'plot'
                 function used to plot
    '''
    
    
    simno = get_sim_no(sim)
    #if we are plotting all species, call plot in turn
    if type=="all":
        sim = SimBuffer.get_sim_no(simno)
        snapobject = SimBuffer.get_snapshot_extended(sim, snap)
        nspecies = snapobject.GetNTypes()
        for ispecies in range(nspecies):
            plot(x,y,snapobject.GetSpecies(ispecies),snap,simno,(overplot or ispecies>0),autoscale,xunit,yunit,**kwargs)
        return
    command = Commands.ParticlePlotCommand(x, y, type, snap, simno, overplot, autoscale, xunit, yunit,**kwargs)
    data = command.prepareData(Singletons.globallimits)
    Singletons.queue.put([command, data])
    sleep(0.001)


#------------------------------------------------------------------------------
def plot_vs_time(y,sim="current",overplot=False,autoscale=True, xunit="default", yunit="default"):
    simno = get_sim_no(sim)
    command = Commands.PlotVsTime(y,simno,overplot,autoscale,xunit,yunit)
    data = command.prepareData(Singletons.globallimits)
    Singletons.queue.put([command, data])
    

#------------------------------------------------------------------------------
def render(x, y, render, snap="current", sim="current", overplot=False, autoscale=True, 
           autoscalerender=True, coordlimits=None, zslice=None,
           xunit="default", yunit="default", renderunit="default",
           res=64, interpolation='nearest'):
    '''Create a rendered plot from selected particle data.

Required arguments:
    x          : Quantity on the x-axis. Must be a string.
    y          : Quantity on the y-axis. Must be a string.
    renderdata : Quantity to be rendered. Must be a string.
        
Optional arguments:
    snap       : Number of the snapshot to plot. Defaults to 'current'.       
    sim        : Number of the simulation to plot. Defaults to 'current'.    
    overplot   : If True, overplots on the previous existing plot rather
                 than deleting it. Defaults to False.
    autoscale  : If True (default), the coordinate limits of the plot are
                 set automatically.  Can also be set to 'x' or 'y' to specify
                 that only one of the axis has to use autoscaling.
                 If False, autoscaling is not used. On an axis that does not
                 have autoscaling turned on, global limits are used if
                 defined for the plotted quantity.
    autoscalerender : Same as the autoscale, but for the rendered quantity.
    coordlimits : Specify the coordinate limits for the plot. In order of
                  precedence, the limits are set in this way:
                  - What this argument specifies. The value must be an
                    iterable of 4 elements: (xmin, xmax, ymin, ymax).
                  - If this argument is None (default), global settings for
                    the quantity are used.
                  - If global settings for the quantity are not defined,
                    the min and max of the data are used.
    zslice     : z-coordinate of the slice when doing a slice rendering.
                 Default is None, which produces a column-integrated plot.
                 If you set this variable, instead a slice rendering will
                 be done. 
    xunit      : Specify the unit to use for the plotting for the quantity
                 on the x-axis.
    yunit      : Specify the unit to use for the plotting for the quantity
                 on the y-axis.
    renderunit : Specify the unit to use for the plotting for the rendered
                 quantity.
    res        : Specify the resolution. Can be an integer number, in which
                 case the same resolution will be used on the two axes, or a
                 tuple (e.g., (xres, yres)) of two integer numbers, if you
                 want to specify different resolutions on the two axes.
    interpolation : Specify the interpolation to use. Default is nearest,
                    which will show the pixels of the rendering grid. If one
                    wants to smooth the image, bilinear or bicubic could be
                    used. See pyplot documentation for the full list of
                    possible values.
    '''
    if zslice is not None:
        zslice = float(zslice)
    simno = get_sim_no(sim)
    command = Commands.RenderPlotCommand(x, y, render, snap, simno, overplot, autoscale, autoscalerender, 
                                         coordlimits, zslice, xunit, yunit, renderunit, res, interpolation)
    data = command.prepareData(Singletons.globallimits)
    Singletons.queue.put([command, data])


#------------------------------------------------------------------------------
def renderslice(x, y, renderq, zslice, **kwargs):
    '''Thin wrapper around render that does slice rendering.

Required arguments:
    x          : Quantity on the x-axis. Must be a string.
    y          : Quantity on the y-axis. Must be a string.
    renderq    : Quantity to be rendered. Must be a string.
    zslice     : z-coordinate of the slice.

Optional arguments:
    See render function optional arguments
'''
    render(x, y, renderq, zslice=zslice, **kwargs)


#------------------------------------------------------------------------------
def addrenderslice(x, y, renderq, zslice, **kwargs):
    '''Thin wrapper around renderslice that sets overplot to True.  If autoscale is
not explicitly set, it will be set to False to preserve the existing settings.

Required arguments:
    x          : Quantity on the x-axis. Must be a string.
    y          : Quantity on the y-axis. Must be a string.
    renderq    : Quantity to be rendered. Must be a string.
    zslice     : z-coordinate of the slice.

Optional arguments:
    See render function optional arguments
'''


    try:
        kwargs['autoscale']
    except KeyError:
        kwargs['autoscale']=False
    render(x, y, renderq, zslice=zslice, overplot=True, **kwargs)


#------------------------------------------------------------------------------
def addrender(x, y, renderq, **kwargs):
    '''Thin wrapper around render that sets overplot to True.  If autoscale is not
explicitly set, it will be set to False to preserve the existing settings.

Required arguments:
    x          : Quantity on the x-axis. Must be a string.
    y          : Quantity on the y-axis. Must be a string.
    renderdata : Quantity to be rendered. Must be a string.
        
Optional arguments:
    See render function optional arguments

    '''
    try:
        kwargs['autoscale']
    except KeyError:
        kwargs['autoscale']=False
    render(x, y, renderq, overplot=True, **kwargs)


#------------------------------------------------------------------------------
def limit (quantity, min=None, max=None, auto=False, window='current', subfigure='current'):
    '''Set plot limits. Quantity is the quantity to limit.
    
Required arguments:
    quantity   : Set limits of this variable. Must be a string.

Optional arguments:
    min        : Minimum value of variable range.
    max        : Maximum value of variable range.
    auto       : If auto is set to True, then the limits for that quantity are 
                 set automatically. Otherwise, use the one given by x and y.
    window     : If window is set to 'global' is available, then any changes
                 will affect also future plots that do not have autoscaling 
                 turned on.
    subfigure  : If subfigure is set to 'all', the limits in all the figures or                 in all the subfigures of the current figure are set.
    '''
    if window=='all' and subfigure=='current':
        subfigure=='all'
    command = Commands.LimitCommand (quantity, min, max, auto, window, subfigure)
    Singletons.queue.put([command,None])
    if window=='global':
        okflag=Singletons.completedqueue.get()
        print okflag


#------------------------------------------------------------------------------
def addplot (x,y, **kwargs):
    '''Thin wrapper around plot that sets overplot to True.  All the other
arguments are the same. If autoscale is not explicitly set, it will be set
to False to preserve the existing settings.

Required arguments:
    x          : Quantity on the x-axis. Must be a string.
    y          : Quantity on the y-axis. Must be a string.
        
Optional arguments:
    See plot function optional arguments
    '''
    try:
        kwargs['autoscale']
    except KeyError:
        kwargs['autoscale']=False
    plot(x,y, overplot=True, **kwargs)


#------------------------------------------------------------------------------
def next():
    '''Advances the current snapshot of the current simulation'''
    try:
        snap(SimBuffer.get_no_next_snapshot())
    except BufferException as e:
        handle(e)


#------------------------------------------------------------------------------
def previous():
    '''Decrements the current snapshot of the current simulation'''
    try:
        snap(SimBuffer.get_no_previous_snapshot())
    except BufferException as e:
        handle(e)


#------------------------------------------------------------------------------
def snap(no):
    '''Jump to the given snapshot number of the current simulation.  Note that you
can use standard Numpy index notation (e.g., -1 is the last snapshot).

Required arguments:
    snapno     : Snapshot number
'''
    no = int(no)
    try:
        SimBuffer.set_current_snapshot_number(no)
    except BufferException as e:
        handle(e)
    
    update("current")


#------------------------------------------------------------------------------
def window(no = None):
    '''Changes the current window to the number specified. If the window
doesn\'t exist, recreate it.

Required arguments:
    winno      : Window number
    '''
    command = Commands.WindowCommand(no)
    data = None
    Singletons.queue.put([command,data])


#------------------------------------------------------------------------------
def subfigure(nx, ny, current):
    '''Creates a subplot in the current window.

Required arguments:
    nx         : x-grid size
    ny         : y-grid size
    current    : id of active sub-figure.  If sub-figure already exists,
                 then this sets the new active sub-figure.
    '''
    command = Commands.SubfigureCommand(nx, ny, current)
    data = None
    Singletons.queue.put([command,data])


#------------------------------------------------------------------------------
def newsim(paramfile=None, ndim=None):
    '''Create a new simulation object. Need to specify either the parameter file, either the number
    of dimensions. Note that it is not possible to change the number of dimensions afterwards.'''
    return SimBuffer.newsim(paramfile=paramfile, ndim=ndim)


#------------------------------------------------------------------------------
def setupsim():
    '''Set up the current simulation object. Note that after calling this function, no parameter change
    it\'s possible!'''
    sim = SimBuffer.get_current_sim()
    sim.SetupSimulation()


#------------------------------------------------------------------------------
def run(no=None):
    '''Run a simulation. If no argument is given, run the current one;
    otherwise queries the buffer for the given simulation number.
    If the simulation has not been setup, does it before running. 
    '''
    #gets the correct simulation object from the buffer
    try:
        if no is None:
            sim = SimBuffer.get_current_sim()
        else:
            sim = SimBuffer.get_sim_no(no)
    except BufferError as e:
        handle(e)
        
    #setup the simulation
    if not sim.setup:
        sim.SetupSimulation()
    SimBuffer.load_live_snapshot(sim)

    while sim.t < sim.tend and sim.Nsteps < sim.Nstepsmax:
        #TODO: maybe some of these operations could be done in another thread, so that the computation is
        #not slowed down when compared to the stand-alone c++ executable
        #But need to think carefully, because of the bloody GIL...
        snap_list = sim.InteractiveRun()
        for snap in snap_list:
            SimBuffer.add_snapshot(snap, sim)
        
        SimBuffer.load_live_snapshot(sim)
        update("live")


#------------------------------------------------------------------------------
def block():
    '''Stops the execution flow until the user presses enter.
    Useful in scripts, allowing to see a plot (which gets closed
    when the execution flow reaches the end of the script'''
    print "Press enter to quit..."
    raw_input()


#------------------------------------------------------------------------------
def update(type=None):
    '''Updates all the plots. You should never call directly this function,
    because all the plotting functions should call this function for you. 
    If you run into a situation when you need it, please contact the authors, because
    you probably just spotted a bug in the code.''' 
    #updates the plots
    for command in Singletons.commands:
        updateplot=False
        if type is None:
            updateplot=True
        else:
            try:
                if command.snap == type:
                    updateplot=True
            except AttributeError:
                updateplot=False
        if updateplot:
            data = command.prepareData(Singletons.globallimits)
            Singletons.queue.put([command, data])


#------------------------------------------------------------------------------
def savefig(name):
    '''Saves the current figure with the given name.
    Note that matplotlib figures out automatically the type of the file
    from the extension'''
    command = Commands.SaveFigCommand(name)
    data = None
    Singletons.queue.put([command,data])
    time.sleep(1e-3)


#------------------------------------------------------------------------------
def switch_nongui():
    '''Switches matplotlib backend, disabling interactive plotting.
    Useful in scripts where no interaction is required'''
    command = Commands.SwitchNonGui()
    data = None
    Singletons.queue.put([command,data])
    time.sleep(1e-3)


#------------------------------------------------------------------------------
def plotanalytical(x=None, y=None, snap = "current", sim = "current", overplot = True, 
                   autoscale = True, xunit="default", yunit="default"):
    '''Plots the analytical solution'''
        
    #TODO: figure out automatically the quantities to plot depending on current window    
    
    simno = get_sim_no(sim)
    command = Commands.AnalyticalPlotCommand(x, y, snap, simno, overplot, autoscale, xunit, yunit)
    data = command.prepareData(Singletons.globallimits)
    Singletons.queue.put([command, data])


#------------------------------------------------------------------------------
def rescale(quantity, unitname, window="current", subfig="current"):
    '''Rescales the specified quantity in the specified window to the specified unit'''
    command = Commands.RescaleCommand(quantity, unitname, window)
    Singletons.queue.put([command,None])
    okflag = Singletons.completedqueue.get()
    print okflag
    update()


#------------------------------------------------------------------------------
def sims():
    '''Print a list of the simulations'''
    print "These simulations are currently loaded into memory:"
    for num, sim in enumerate(SimBuffer.simlist):
        print str(num) + ' ' + sim.simparams.stringparams["run_id"]


#------------------------------------------------------------------------------
def snaps(simno):
    '''For the given simulation number, print a list of all the snapshots'''
    simno = int(simno)
    sim=SimBuffer.get_sim_no(simno)
    print "The run_id of the requested simulation is " + sim.simparams.stringparams["run_id"]
    print "These are the snapshots that we know about for this simulation:"
    for num, snap in enumerate(sim.snapshots):
        #TODO: snap.t is set correctly only the first time that the snapshot is read from the disc, should be fixed
        print str(num) + ' ' + snap.filename + " " + str(snap.t)
    try:
        live=None
        live=sim.live
    except AttributeError:
        pass
    if live is not None:
        print "In addition, there is a live snapshot in memory, at time " + str(live.t) 


#------------------------------------------------------------------------------
def set_current_sim(simno):
    '''Set the current simulation to the given number. Returns the newly set current simulation'''
    simno = int(simno)
    return SimBuffer.set_current_sim_no(simno)


#------------------------------------------------------------------------------
def get_sim_no(sim):
    if sim == "current":
        simno = SimBuffer.get_current_sim_no()
    else:
        simno = sim
    return simno


#------------------------------------------------------------------------------
def sigint(signum, frame):
    cleanup()


#------------------------------------------------------------------------------
def cleanup():
    Singletons.queue.put(["STOP",None])
    print "Waiting for background processes to finish..."
    plottingprocess.join()
    import sys
    sys.exit()


#------------------------------------------------------------------------------
def init():
    global plottingprocess
    plottingprocess = PlottingProcess(Singletons.queue, Singletons.commands, Singletons.completedqueue, Singletons.globallimits)
    plottingprocess.start()
    CreateUserQuantity('r','sqrt(x^2+y^2+z^2)')
    CreateUserQuantity('R','sqrt(x^2+y^2)')
    CreateUserQuantity('phi','arctan2(y,x)')
    CreateUserQuantity('theta','arccos(z/r)')
    CreateUserQuantity('vr','sin(theta)*cos(phi)*vx+sin(theta)*sin(phi)*vy+cos(theta)*vz')
    CreateUserQuantity('ar','sin(theta)*cos(phi)*ax+sin(theta)*sin(phi)*ay+cos(theta)*az')
    
    from compute import COM
    CreateTimeData('com_x',COM)
    CreateTimeData('com_y',COM,quantity='y')
    CreateTimeData('com_z',COM,quantity='z')

init()
 
signal.signal(signal.SIGINT, sigint)
signal.signal(signal.SIGTERM, sigint)
signal.signal(signal.SIGSEGV, sigint)
atexit.register(cleanup)


#------------------------------------------------------------------------------
if __name__=="__main__":
    loadsim('TEST')
    plot("x","rho")
    plotanalytical("x","rho")
    limit('x', -10.0, 10.0)
    snap(1)
    import time; time.sleep(2)
    next(); time.sleep(2)
    snap(8)
    limit('x', 'auto')
    print 'L1 error norm : ',L1errornorm("x","rho",1.0,8.0)

    block()
#    
#    loadsim('TEST')
#    plot("x","y", snap=0)
#    addplot("x", "y")
#    window()
#    plot("vx", "vy")
#    plot("vx", "x")
#    window()
#    plot("x","rho")
#    window()
#    subfigure(2,2,1)
#    plot("x", "y")
#    subfigure(2,2,2)
#    plot("vx", "vy")
#    subfigure(2,2,3)
#    plot("x", "rho")
#    subfigure(2,2,4)
#    plot("rho", "h")
#    addplot("rho", "m")
#    window(3)
#    addplot("rho", "h")
#    snap(99)
#    for i in range(10):
#        time.sleep(1)
#        previous()

