import __main__
import atexit
from multiprocessing import Manager, Queue
from plotting import PlottingProcess
from seren.analysis.SimBuffer import SimBuffer, BufferException

manager= Manager()

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
from scipy import interpolate
import signal
import numpy
from time import sleep

#figure out if we are in interactive mode
try:
    __main__.__file__
    interactive=False
except AttributeError:
    interactive=True

def handle(e):
    '''This functions takes care of printing information about an error,
    if we are in interactive mode, or re-raising it, if we are in script mode
    (so that the execution of the script can stop if nobody catches the exception)
    '''
    if interactive:
        print str(e)
    else:
        raise e
    
def loadsim(run_id, fileformat = 'ascii', buffer_flag = 'cache'):
    '''Given the run_id of a simulation, reads it from the disk'''
    SimBuffer.loadsim(run_id, fileformat=fileformat, buffer_flag=buffer_flag)
    return SimBuffer.get_current_sim()
    
def plot(x,y, snap="current", sim="current", overplot = False, autoscale = True, xunit="default", yunit="default"):
    '''Particle plot. Needed arguments:
    x
        Quantity on the x-axis. Must be a string.
    y
        Quantity on the x-axis. Must be a string.
        
    Optional arguments:
    snap
        Number of the snapshot to plot. Defaults to 'current'.       
    sim
        Number of the simulation to plot. Defaults to 'current'.    
    overplot
        If True, overplots on the previous existing plot rather than deleting it. Defaults to False.
    autoscale
        If True (default), the limits of the plot are set automatically.
        Can also be set to 'x' or 'y' to specify that only one of the axis has to use autoscaling.
        If False, autoscaling is not used. On an axis that does not have autoscaling turned on,
        global limits are used if defined for the plotted quantity.
    xunit
        Specify the unit to use for the plotting for the quantity on the x-axis.
    yunit
        Specify the unit to use for the plotting for the quantity on the y-axis.
    '''
    
    
    simno = get_sim_no(sim)
    command = Commands.ParticlePlotCommand(x, y, snap, simno, overplot, autoscale, xunit, yunit)
    data = command.prepareData(Singletons.globallimits)
    Singletons.queue.put([command, data])
    sleep(0.001)

def render(x, y, render, snap="current", sim="current", overplot=False, autoscale=True, 
           autoscalerender=True, coordlimits=None, zslice=None,
           xunit="default", yunit="default", renderunit="default",
           res=64, interpolation='nearest'):
    '''Render plot. Needed arguments:
    x
        Quantity on the x-axis. Must be a string.
    y
        Quantity on the x-axis. Must be a string.
    render
        Quantity to be rendered. Must be a string.
        
    Optional arguments:
    snap
        Number of the snapshot to plot. Defaults to 'current'.       
    sim
        Number of the simulation to plot. Defaults to 'current'.    
    overplot
        If True, overplots on the previous existing plot rather than deleting it. Defaults to False.
    autoscale
        If True (default), the coordinate limits of the plot are set automatically.
        Can also be set to 'x' or 'y' to specify that only one of the axis has to use autoscaling.
        If False, autoscaling is not used. On an axis that does not have autoscaling turned on,
        global limits are used if defined for the plotted quantity.
    autoscalerender
        Same as the previous one, but for the rendered quantity.
    coordlimits
        Specify the coordinate limits for the plot. In order of precedence, the limits are set in this way:
        -What this argument specifies. The value must be an iterable of 4 elements: (xmin, xmax, ymin, ymax).
        -If this argument is None (default), global settings for the quantity are used.
        -If global settings for the quantity are not defined, the min and max of the data are used.
    zslice
        z-coordinate of the slice when doing a slice rendering. Default is None, which produces a column-integrated
        plot. If you set this variable, instead a slice rendering will be done. 
    xunit
        Specify the unit to use for the plotting for the quantity on the x-axis.
    yunit
        Specify the unit to use for the plotting for the quantity on the y-axis.
    renderunit
        Specify the unit to use for the plotting for the rendered quantity.
    res
        Specify the resolution. Can be an integer number, in which case the same resolution will be used on the two axes,
        or a tuple (e.g., (xres, yres)) of two integer numbers, if you want to specify different resolutions on the two axes.
    interpolation
        Specify the interpolation to use. Default is nearest, which will show the pixels of the rendering grid. If one
        wants to smooth the image, bilinear or bicubic could be used. See pyplot documentation for the full list
        of possible values.
    '''
    simno = get_sim_no(sim)
    command = Commands.RenderPlotCommand(x, y, render, snap, simno, overplot, autoscale, autoscalerender, 
                                         coordlimits, zslice, xunit, yunit, renderunit, res, interpolation)
    data = command.prepareData(Singletons.globallimits)
    Singletons.queue.put([command, data])

def renderslice(x, y, renderq, zslice, **kwargs):
    '''Thin wrapper around render that does slice rendering.
    All the other arguments are the same. Needs the z-coordinate
    of the slice as a compulsory parameter'''
    render(x, y, renderq, zslice=zslice, **kwargs)
    
def addrenderslice(x, y, renderq, zslice, **kwargs):
    '''Thin wrapper around renderslice that sets overplot to True.
    All the other arguments are the same. If autoscale is not
    explicitly set, it will be set to False to preserve the
    existing settings.'''
    try:
        kwargs['autoscale']
    except KeyError:
        kwargs['autoscale']=False
    render(x, y, renderq, zslice=zslice, overplot=True, **kwargs)

def addrender(x, y, renderq, **kwargs):
    '''Thin wrapper around render that sets overplot to True.
    All the other arguments are the same. If autoscale is not
    explicitly set, it will be set to False to preserve the
    existing settings.'''
    try:
        kwargs['autoscale']
    except KeyError:
        kwargs['autoscale']=False
    render(x, y, renderq, overplot=True, **kwargs)
    
def limit (quantity, min=None, max=None, auto=False, window='current', subfigure='current'):
    '''Set plot limits. Quantity is the quantity to limit. If auto is set to True, then
    the limits for that quantity are set automatically. Otherwise, use the one given by x and y.
    By default, changes the limits only for the current subfigure of the current plot. One
    can specify the number for them or use the special keyword 'all' to change the limits in all
    the figures or in all the subfigures of the current figure. Also the keyword 'global' is available
    for window, which means that the change will affect also future plots that do not have autoscaling
    turned on. 
    '''
    if window=='all' and subfigure=='current':
        subfigure=='all'
    command = Commands.LimitCommand (quantity, min, max, auto, window, subfigure)
    Singletons.queue.put([command,None])
    if window=='global':
        okflag=Singletons.completedqueue.get()
        print okflag

def addplot (x,y, **kwargs):
    '''Thin wrapper around plot that sets overplot to True.
    All the other arguments are the same. If autoscale is not
    explicitly set, it will be set to False to preserve the
    existing settings.'''
    try:
        kwargs['autoscale']
    except KeyError:
        kwargs['autoscale']=False
    plot(x,y, overplot=True, **kwargs)
    
def next():
    '''Advances the current snapshot of the current simulation'''
    try:
        snap(SimBuffer.get_no_next_snapshot())
    except BufferException as e:
        handle(e)
        
def previous():
    '''Decrements the current snapshot of the current simulation'''
    try:
        snap(SimBuffer.get_no_previous_snapshot())
    except BufferException as e:
        handle(e)
        
def snap(no):
    '''Jump to the given snapshot number of the current simulation'''
    try:
        SimBuffer.set_current_snapshot_number(no)
    except BufferException as e:
        handle(e)
    
    update("current")
        
def window(no = None):
    '''Changes the current window to the number specified. If the
    window doesn't exist, recreate it.'''
    command = Commands.WindowCommand(no)
    data = None
    Singletons.queue.put([command,data])

def subfigure(nx, ny, current):
    '''Creates a subplot in the current window.
    The arguments nx and ny specify the grid size; the current
    arguments marks the subplot that will become the current one.
    If the plot already exists, just sets it as the current one.
    '''
    command = Commands.SubfigureCommand(nx, ny, current)
    data = None
    Singletons.queue.put([command,data])

def newsim(paramfile):
    '''Reads a parameter file and setups a simulation from it'''
    return SimBuffer.newsim(paramfile)

def newsimfromparams(paramfile):
    '''Creates a new simulation object and parses the parameter file,
    but does not do the setup, still allowing you to programmatically
    change some of the parameters'''
    return SimBuffer.newsimfromparams(paramfile)

def run(no=None):
    '''Run a simulation. If no argument is given, run the current one;
    otherwise queries the buffer for the given simulation number.
    If the simulation has not been setup (e.g., if you have used
    newsimfromparams), does it before running. 
    '''
    #gets the correct simulation object from the buffer
    try:
        if no is None:
            sim = SimBuffer.get_current_sim()
        else:
            sim = SimBuffer.get_sim_no(no)
    except BufferError as e:
        handle(e)
    #TODO: treat this as an exception
    if not sim.setup:
        print "The selected simulation was not set-up, we will do it for you"
        sim.SetupSimulation()
    
    while sim.t < sim.tend and sim.Nsteps < sim.Nstepsmax:
        #TODO: maybe some of these operations could be done in another thread, so that the computation is
        #not slowed down when compared to the stand-alone c++ executable
        #But need to think carefully, because of the bloody GIL...
        sim.InteractiveRun()
        SimBuffer.load_live_snapshot(sim)
        update("live")

def block():
    '''Stops the execution flow until the user presses enter.
    Useful in scripts, allowing to see a plot (which gets closed
    when the execution flow reaches the end of the script'''
    print "Press enter to quit..."
    raw_input()

def update(type=None):
    '''Updates all the plots. The user should never call directly this function.''' 
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

def savefig(name):
    '''Saves the current figure with the given name.
    Note that matplotlib figures out automatically the type of the file
    from the extension'''
    command = Commands.SaveFigCommand(name)
    data = None
    Singletons.queue.put([command,data])
    time.sleep(1e-3)

def switch_nongui():
    '''Switches matplotlib backend, disabling interactive plotting.
    Useful in scripts where no interaction is required'''
    command = Commands.SwitchNonGui()
    data = None
    Singletons.queue.put([command,data])
    time.sleep(1e-3)

def plotanalytical(x=None, y=None, snap = "current", sim = "current", overplot = True, 
                   autoscale = True, xunit="default", yunit="default"):
    '''Plots the analytical solution'''
        
    #TODO: figure out automatically the quantities to plot depending on current window    
    
    simno = get_sim_no(sim)
    command = Commands.AnalyticalPlotCommand(x, y, snap, simno, overplot, autoscale, xunit, yunit)
    data = command.prepareData(Singletons.globallimits)
    Singletons.queue.put([command, data])

def rescale(quantity, unitname, window="current", subfig="current"):
    '''Rescales the specified quantity in the specified window to the specified unit'''
    command = Commands.RescaleCommand(quantity, unitname, window)
    Singletons.queue.put([command,None])
    okflag = Singletons.completedqueue.get()
    print okflag
    update()

def L1errornorm(x=None, y=None, xmin=None, xmax=None, sim = "current", snap = "current"):
    '''Computes the L1 error norm from the simulation data relative to the analytical solution'''
    
    #get the simulation number from the buffer
    simno = get_sim_no(snap)
    
    #istantiate and setup the command object to retrieve analytical solution
    command1 = Commands.AnalyticalPlotCommand(x, y, snap, simno)
    adata = command1.prepareData(Singletons.globallimits)

    #istantiate and setup the 2nd command object to retrieve particle data
    command2 = Commands.ParticlePlotCommand(x, y, snap, simno)
    pdata = command2.prepareData(Singletons.globallimits)

    #cut arrays if limits are provided
    if xmin != None and xmax != None:
        aindex = numpy.logical_and( adata.x_data > xmin, adata.x_data < xmax)
        adata.x_data = adata.x_data[aindex]
        adata.y_data = adata.y_data[aindex]
        pindex = numpy.logical_and( pdata.x_data > adata.x_data.min(),
                                    pdata.x_data < adata.x_data.max() )
        pdata.x_data = pdata.x_data[pindex]
        pdata.y_data = pdata.y_data[pindex]

    
        
    #prepare interpolation function from analytical data
    #f = interpolate.interp1d(adata.x_data[::-1], adata.y_data[::-1], kind = 'linear', axis=0, bounds_error = False)
    f = interpolate.interp1d(adata.x_data[::], adata.y_data[::], kind = 'linear', axis=0, bounds_error = False)

    #compute error norm of particle data relative to analytical data
    L1 = numpy.linalg.norm((pdata.y_data - f(pdata.x_data)), ord=1)/pdata.x_data.size
    return L1

def get_sim_no(sim):
    if sim == "current":
        simno = SimBuffer.get_current_sim_no()
    else:
        simno = sim
    return simno

def sigint(signum, frame):
    cleanup()
    
def cleanup():
    Singletons.queue.put(["STOP",None])
    print "Waiting for background processes to finish..."
    plottingprocess.join()
    import sys
    sys.exit()

def init():
    global plottingprocess
    plottingprocess = PlottingProcess(Singletons.queue, Singletons.commands, Singletons.completedqueue, Singletons.globallimits)
    plottingprocess.start()

init()
 
signal.signal(signal.SIGINT, sigint)
signal.signal(signal.SIGTERM, sigint)
signal.signal(signal.SIGSEGV, sigint)
atexit.register(cleanup)

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

