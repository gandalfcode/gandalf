#==============================================================================
#  facade.py
#  Main gandalf library front-end when invoking gandalf from within python.
#
#  This file is part of GANDALF :
#  Graphical Astrophysics code for N-body Dynamics And Lagrangian Fluids
#  https://github.com/gandalfcode/gandalf
#  Contact : gandalfcode@gmail.com
#
#  Copyright (C) 2013  D. A. Hubber, G. Rosotti
#
#  GANDALF is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 2 of the License, or
#  (at your option) any later version.
#
#  GANDALF is distributed in the hope that it will be useful, but
#  WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  General Public License (http://www.gnu.org/licenses) for more details.
#==============================================================================
import __main__
import atexit
import time
import types
import defaults
from multiprocessing import Manager, Queue, Event
from plotting import PlottingProcess
from gandalf.analysis.SimBuffer import SimBuffer, BufferException

manager = Manager()

#TODO: in all the Python code, raise proper exceptions rather than a generic Exception
#TODO: the tests should not fail


#------------------------------------------------------------------------------
class Singletons_master:
    '''Container class for singletons object. They are:
queue          : Queue for sending commands to the plotting process
commands       : List of the commands shared with the plotting process.
                 Caution: if you modify a command, you must reassign it in the
                 list to make the changes propagate to the other process
completedqueue : Queue used from the plotting process to signal the completion
                 of a command
globallimits   : Dict that for each quantity gives the limits
'''
    queue = Queue()
    commands = manager.list()
    completedqueue = Queue()
    globallimits = manager.dict()
    free = Event()

class Singletons_serial(Singletons_master):
    @staticmethod
    def place_command(objects):
        Singletons.queue.put(objects)
        command, data = Singletons.queue.get()
        command.processCommand(plotting, data)

class Singletons_parallel(Singletons_master):
    @staticmethod
    def place_command(objects):
        Singletons.queue.put(objects)

if defaults.parallel:
    Singletons=Singletons_parallel
else:
    Singletons=Singletons_serial

import commandsource as Commands
from data_fetcher import CreateUserQuantity, CreateTimeData, UserQuantity
from data_fetcher import _KnownQuantities as KnownQuantities
import signal
from time import sleep
from statistics import structure_function
import subprocess
import tempfile
import glob
import os



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
def loadsim(run_id, fileformat=None, buffer_flag='cache'):
    '''Given the run_id of a simulation, reads it from the disk.
Returns the newly created simulation object.

Required arguments:
    run_id      : Simulation run identification string.

Optional qrguments:
    fileformat  : Format of all snapshot files of simulation.
    buffer_flag : Record snapshot data in simulation buffer.
'''
    SimBuffer.loadsim(run_id, fileformat=fileformat, buffer_flag=buffer_flag)
    return SimBuffer.get_current_sim()



class Plotting:
    def __init__(self):
        self.lastid=0
        import matplotlib.pyplot as plt
        self.plt=plt
        self.axesimages = {}
        self.commands=Singletons.commands
        self.commandsfigures = {}
        self.quantitiesfigures = {}
        self.globallimits = Singletons.globallimits


    def command_in_list (self, id):
        for command in Singletons.commands:
            if command.id==id:
                return True
        return False


#------------------------------------------------------------------------------
def plot(x, y, type="default", snap="current", sim="current",
         overplot=False, autoscale=False, xunit="default", yunit="default",
         xaxis="linear", yaxis="linear", **kwargs):
    '''Plot particle data as a scatter plot.  Creates a new plotting window if
one does not already exist.

Required arguments:
    x          : Quantity on the x-axis. Must be a string.
    y          : Quantity on the y-axis. Must be a string.

Optional arguments:
    type       : The type of the particles to plot (e.g. 'star' or 'sph').
    snap       : Number of the snapshot to plot. Defaults to 'current'.
    sim        : Number of the simulation to plot. Defaults to 'current'.
    overplot   : If True, overplots on the previous existing plot rather
                 than deleting it. Defaults to False.
    autoscale  : If True, the limits of the plot are set
                 automatically.  Can also be set to 'x' or 'y' to specify
                 that only one of the axis has to use autoscaling.
                 If False (default), autoscaling is not used. On an axis that does
                 not have autoscaling turned on, global limits are used
                 if defined for the plotted quantity.
    xunit      : Specify the unit to use for the plotting for the quantity
                 on the x-axis.
    yunit      : Specify the unit to use for the plotting for the quantity
                 on the y-axis.
    **kwargs   : Extra keyword arguments will be passed to matplotlib.
'''
    simno = get_sim_no(sim)
    overplot=to_bool(overplot)
    # If we are plotting all particle species, call plot in turn
    if type=="all":
        sim = SimBuffer.get_sim_no(simno)
        snapobject = SimBuffer.get_snapshot_extended(sim, snap)
        nspecies = snapobject.GetNTypes()
        for ispecies in range(nspecies):
            plot(x,y,snapobject.GetSpecies(ispecies),snap,simno,
                 (overplot or ispecies>0),autoscale,xunit,yunit,
                 xaxis,yaxis,**kwargs)
        return
    command = Commands.ParticlePlotCommand(x, y, type, snap, simno, overplot,
                                           autoscale, xunit, yunit,
                                           xaxis, yaxis, **kwargs)
    data = command.prepareData(Singletons.globallimits)
    Singletons.place_command([command, data])
    sleep(0.001)
    return data


#------------------------------------------------------------------------------
def time_plot(x, y, sim="current", overplot=False, autoscale=False,
              xunit="default", yunit="default", xaxis="linear",
              yaxis="linear", idx=None, idy=None, id=None,
              typex="default", typey="default", type="default", **kwargs):
    '''Plot two quantities as evolved in time one versus the another.  Creates
a new plotting window if one does not already exist.

Required arguments:
    x          : Quantity on x-axis. Must be a string. The quantity is looked
                 up in the quantities defined as a function of time. If it is
                 not found there, then we try to interpret it as a quantity
                 defined for a particle. In this case, the user needs to pass
                 either idx either id to specify which particle he wishes
                 to look-up.
    y          : Quantity on y-axis.  Must be a string. The interpretation is
                 like for the previous argument.

Optional arguments:
    sim        : Number of the simulation to plot. Defaults to 'current'.
    overplot   : If True, overplots on the previous existing plot rather
                 than deleting it. Defaults to False.
    autoscale  : If True, the limits of the plot are set
                 automatically.  Can also be set to 'x' or 'y' to specify
                 that only one of the axis has to use autoscaling.
                 If False (default), autoscaling is not used. On an axis that
                 does not have autoscaling turned on, global limits are used
                 if defined for the plotted quantity.
    xunit      : Specify the unit to use for the plotting for the quantity
                 on the x-axis.
    yunit      : Specify the unit to use for the plotting for the quantity
                 on the y-axis.
    idx        : id of the particle to plot on the x-axis. Ignored if the
                 quantity given (e.g., com_x) does not depend on the id.
    idy        : same as previous, on the y-axis.
    id         : same as the two previous ones. To be used when the id is the
                 same on both axes. If set, overwrites the passed idx and idy.
    typex      : type of particles on the x-axis. Ignored if the quantity
                 given does not depend on it
    typey      : as the previous one, on the y-axis.
    type       : as the previous ones, for both axis at the same time. If set,
                 overwrites typex and typey.
'''
    simno = get_sim_no(sim)
    overplot = to_bool(overplot)
    command = Commands.TimePlot(x, y,simno,overplot,autoscale,
                                  xunit,yunit,xaxis,yaxis,idx, idy, id,
                                  typex, typey, type, **kwargs)
    data = command.prepareData(Singletons.globallimits)

    Singletons.place_command([command, data])
    return data


#------------------------------------------------------------------------------
def render(x, y, render, snap="current", sim="current", overplot=False,
           autoscale=False, autoscalerender=False, coordlimits=None,
           zslice=None, xunit="default", yunit="default",
           renderunit="default", res=64, interpolation='nearest',lognorm=False,**kwargs):
    '''Create a rendered plot from selected particle data.

Required arguments:
    x          : Quantity on the x-axis. Must be a string.
    y          : Quantity on the y-axis. Must be a string.
    renderdata : Quantity to be rendered. Must be a string.

Optional arguments:
    snap       : Number of the snapshot to plot. Defaults to \'current\'.
    sim        : Number of the simulation to plot. Defaults to \'current\'.
    overplot   : If True, overplots on the previous existing plot rather
                 than deleting it. Defaults to False.
    autoscale  : If True, the coordinate limits of the plot are set
                 automatically.  Can also be set to 'x' or 'y' to specify
                 that only one of the axis has to use autoscaling.
                 If False (default), autoscaling is not used. On an axis that
                 does not have autoscaling turned on, global limits are used
                 if defined for the plotted quantity.
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
    lognorm    : Boolean flag specifying wheter the colour scale should be
                 logarithmic (default: linear). If you want to customise the
                 limits, use the vmin and vmax flags which are passed to
                 matplotlib
'''
    if zslice is not None:
        zslice = float(zslice)
    simno = get_sim_no(sim)
    overplot = to_bool(overplot)
    autoscalerender = to_bool(autoscalerender)
    if coordlimits is not None and isinstance(coordlimits, types.StringTypes):
        coordlimits = to_list (coordlimits, float)
    if isinstance(res, types.StringTypes):
        if res[0]=='[' or res[0]=='(':
            res = to_list(res,int)
        else:
            res = int(res)
    command = Commands.RenderPlotCommand(x, y, render, snap, simno, overplot,
                                         autoscale, autoscalerender,
                                         coordlimits, zslice, xunit, yunit,
                                         renderunit, res, interpolation,lognorm,**kwargs)
    data = command.prepareData(Singletons.globallimits)
    Singletons.place_command([command, data])
    return data


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
    data=render(x, y, renderq, zslice=zslice, **kwargs)
    return data


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
    data=render(x, y, renderq, zslice=zslice, overplot=True, **kwargs)
    return data


#------------------------------------------------------------------------------
def addrender(x, y, renderq, **kwargs):
    '''Thin wrapper around render that sets overplot to True.  If autoscale is
not explicitly set, it will be set to False to preserve the existing settings.

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
    data=render(x, y, renderq, overplot=True, **kwargs)
    return data



#------------------------------------------------------------------------------
def make_movie(filename, snapshots='all', window_no=0, fps=24):
    '''Generates movie for plots generated in given window'''

    # Remove all temporary files in the directory (in case they still exist)
    tmpfilelist = glob.glob('tmp.?????.png')
    for file in tmpfilelist:
        os.remove(file)

    sim = SimBuffer.get_current_sim()
    nframes = len(sim.snapshots)

    # Loop through all snapshots and create temporary images
    if snapshots == 'all':
        for isnap in range(len(sim.snapshots)):
            snap(isnap)
            tmpfile = 'tmp.' + str(isnap).zfill(5) + '.png'
            savefig(tmpfile)

    # Wait until all plotting processes have finished before making mp4 file
    if defaults.parallel:
        Singletons.free.wait()

    # Now join all temporary files together with ffmpeg
    subprocess.call(["ffmpeg","-y","-r",str(fps),"-i", "tmp.%05d.png", \
                     "-vcodec","mpeg4", "-qscale","5", "-r", str(fps), \
                     filename])

    # Now remove all temporary files just created to make movie
    tmpfilelist = glob.glob('tmp.?????.png')
    for file in tmpfilelist:
        os.remove(file)



#------------------------------------------------------------------------------
def limit(quantity, min=None, max=None, auto=False,
          window='current', subfigure='current'):
    '''Set plot limits. Quantity is the quantity to limit.

Required arguments:
    quantity   : Set limits of this variable. Must be a string.

Optional arguments:
    min        : Minimum value of variable range.
    max        : Maximum value of variable range.
    auto       : If auto is set to True, then the limits for that quantity are
                 set automatically. Otherwise, use the one given by max and min.
    window     : If window is set to 'global' is available, then any changes
                 will affect also future plots that do not have autoscaling
                 turned on.
    subfigure  : If subfigure is set to 'all', the limits in all the figures or
                 in all the subfigures of the current figure are set.
'''
    if min is not None:
        min = float(min)
    if max is not None:
        max = float(max)
    if not auto:
        auto=to_bool(auto)
    if window=='all' and subfigure=='current':
        subfigure=='all'
    command = Commands.LimitCommand(quantity, min, max, auto, window, subfigure)
    Singletons.place_command([command,None])
    if window=='global':
        okflag=Singletons.completedqueue.get()
        print okflag


#------------------------------------------------------------------------------
def addplot(x, y, **kwargs):
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
    data=plot(x, y, overplot=True, **kwargs)
    return data


#------------------------------------------------------------------------------
def next():
    '''Advances the current snapshot of the current simulation.
Return the new snapshot, or None if the call failed.'''
    try:
        snapshot=snap(SimBuffer.get_no_next_snapshot())
        return snapshot
    except BufferException as e:
        handle(e)


#------------------------------------------------------------------------------
def previous():
    '''Decrements the current snapshot of the current simulation.
Return the new snapshot, or None if the call failed.'''
    try:
        snapshot=snap(SimBuffer.get_no_previous_snapshot())
        return snapshot
    except BufferException as e:
        handle(e)


#------------------------------------------------------------------------------
def snap(no):
    '''Jump to the given snapshot number of the current simulation.  Note that
you can use standard Numpy index notation (e.g., -1 is the last snapshot).
Return the new snapshot, or None if the call failed.

Required arguments:
    snapno     : Snapshot number
'''
    no = int(no)
    snapshot=None
    try:
        snapshot=SimBuffer.set_current_snapshot_number(no)
    except BufferException as e:
        handle(e)
    if snapshot is not None:
        update("current")
    return snapshot


#------------------------------------------------------------------------------
def window(no = None):
    '''Changes the current window to the number specified. If the window
doesn\'t exist, recreate it.

Required arguments:
    winno      : Window number
'''
    if no is not None:
        no=int(no)
    command = Commands.WindowCommand(no)
    data = None
    Singletons.place_command([command,data])


#------------------------------------------------------------------------------
def subfigure(nx, ny, current):
    '''Creates a subplot in the current window.

Required arguments:
    nx         : x-grid size
    ny         : y-grid size
    current    : id of active sub-figure.  If sub-figure already exists,
                 then this sets the new active sub-figure.
'''
    nx = int(nx)
    ny = int(ny)
    current = int(current)
    command = Commands.SubfigureCommand(nx, ny, current)
    data = None
    Singletons.place_command([command,data])


#------------------------------------------------------------------------------
def newsim(paramfile=None, ndim=None, sim=None):
    '''Create a new simulation object. Need to specify either the parameter
file, or the number of dimensions and the simulation type. Note that it is not
possible to change the number of dimensions afterwards or simulation type
afterwards.
'''
    return SimBuffer.newsim(paramfile=paramfile, ndim=ndim, simtype=sim)


#------------------------------------------------------------------------------
def setupsim():
    '''Set up the current simulation object. Note that after calling this function,
no parameter change it\'s possible.
'''
    sim = SimBuffer.get_current_sim()
    sim.SetupSimulation()
    sim.simparams.RecordParametersToFile()


#------------------------------------------------------------------------------
def run(no=None):
    '''Run a simulation. If no argument is given, run the current one;
otherwise queries the buffer for the given simulation number.
If the simulation has not been setup, does it before running.

Optional arguments:
    no         : Simulation number
'''
    #gets the correct simulation object from the buffer
    try:
        if no is None:
            sim = SimBuffer.get_current_sim()
        else:
            no = int(no)
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
        #But need to think carefully, because of the GIL... (??)
        snap_list = sim.InteractiveRun()
        for snap in snap_list:
            SimBuffer.add_snapshot(snap, sim)

        SimBuffer.load_live_snapshot(sim)
        update("live")


#------------------------------------------------------------------------------
def block():
    '''Stops the execution flow until the user presses 'enter'.
Useful in scripts, allowing to see a plot (which gets closed
when the execution flow reaches the end of the script
'''
    print "Press enter to quit..."
    raw_input()


#------------------------------------------------------------------------------
def update(type=None):
    '''Updates all the plots. You should never call directly this function,
because all the plotting functions should call this function for you.
If you run into a situation when you need it, please contact the authors,
because you probably just spotted a bug in the code.
'''
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
            Singletons.place_command([command, data])


#------------------------------------------------------------------------------
def savefig(name):
    '''Saves the current figure with the given name.  Note that matplotlib
figures out automatically the type of the file from the extension.

Required arguments:
    name       : filename (including extension)

'''
    command = Commands.SaveFigCommand(name)
    data = None
    Singletons.place_command([command,data])
    time.sleep(1e-3)


#------------------------------------------------------------------------------
def switch_nongui():
    '''Switches matplotlib backend, disabling interactive plotting.
Useful in scripts where no interaction is required
'''
    command = Commands.SwitchNonGui()
    data = None
    Singletons.place_command([command,data])
    time.sleep(1e-3)


#------------------------------------------------------------------------------
def plotanalytical(x=None, y=None, ic="default", snap="current", sim="current",
                   overplot=True, autoscale=False, xunit="default",
                   yunit="default", time="snaptime"):
    '''Plots the analytical solution.  Reads the problem type from the \'ic\'
parameter and plots the appropriate solution if implemented.  If no solution
exists, then nothing is plotted.

Optional arguments:
    x          : Quantity on the x-axis. Must be a string.
    y          : Quantity on the y-axis. Must be a string.
    snap       : Number of the snapshot to plot. Defaults to 'current'.
    sim        : Number of the simulation to plot. Defaults to 'current'.
    overplot   : If True, overplots on the previous existing plot rather
                 than deleting it. Defaults to False.
    autoscale  : If True, the limits of the plot are set
                 automatically.  Can also be set to 'x' or 'y' to specify
                 that only one of the axis has to use autoscaling.
                 If False (default), autoscaling is not used. On an axis that does
                 not have autoscaling turned on, global limits are used
                 if defined for the plotted quantity.
    xunit      : Specify the unit to use for the plotting for the quantity
                 on the x-axis.
    yunit      : Specify the unit to use for the plotting for the quantity
                 on the y-axis.
    time       : Plots the analytical solution for the given time.
                 If not set, then reads the time from the sim or snapshot
'''
    #TODO: figure out automatically the quantities to plot depending on current window

    simno = get_sim_no(sim)
    overplot = to_bool(overplot)
    command = Commands.AnalyticalPlotCommand(x, y, ic, snap, simno, overplot,
                                             autoscale, xunit, yunit)

    data = command.prepareData(Singletons.globallimits, time)
    Singletons.place_command([command, data])
    return data


#------------------------------------------------------------------------------
def rescale(quantity, unitname, window="current", subfig="current"):
    '''Rescales the specified quantity in the specified window to the specified unit

Required arguments:
    quantity   : Quantity to be rescaled.  Must be a string.
    unitname   : Required unit for quantity.

Optional qrguments:
    window     : Window containing plot
    subfig     : Sub-figure in window containing plot
'''
    command = Commands.RescaleCommand(quantity, unitname, window)
    Singletons.place_command([command,None])
    okflag = Singletons.completedqueue.get()
    print okflag
    update()


#------------------------------------------------------------------------------
def sims():
    '''Print a list of the simulations to screen'''
    print "These simulations are currently loaded into memory:"
    for num, sim in enumerate(SimBuffer.simlist):
        print str(num) + ' ' + sim.simparams.stringparams["run_id"]


#------------------------------------------------------------------------------
def snaps(simno):
    '''For the given simulation number, print a list of all the snapshots

Required argument:
    simno      : Simulation number from which to print the snapshot list.
'''
    simno = int(simno)
    sim = SimBuffer.get_sim_no(simno)
    print "The run_id of the requested simulation is " + sim.simparams.stringparams["run_id"]
    print "These are the snapshots that we know about for this simulation:"
    for num, snap in enumerate(sim.snapshots):
        #TODO: snap.t is set correctly only the first time that the snapshot is read from the disc, should be fixed
        print str(num) + ' ' + snap.filename + " " + str(snap.t)
    try:
        live = None
        live = sim.live
    except AttributeError:
        pass
    if live is not None:
        print "In addition, there is a live snapshot in memory, at time " + str(live.t)


#------------------------------------------------------------------------------
def set_current_sim(simno):
    '''Set the current simulation to the given number.
Returns the newly set current simulation.

Required argument:
    simno      : Simulation number
'''
    simno = int(simno)
    return SimBuffer.set_current_sim_no(simno)


#------------------------------------------------------------------------------
def get_sim_no(sim):
    '''Returns the simulation id of the currently active simulation object

Required argument:
    sim        : Simulation
'''
    if sim == "current":
        simno = SimBuffer.get_current_sim_no()
    else:
        simno = int(sim)
    return simno

def get_data(quantity, snap="current",type="default",sim="current",unit="default" ):
    '''Returns the array with the data for the given quantity.
    The data is returned scaled to the specified unit
    
    Required argument:
        quantity        :The quantity required. Must be a string
        
    Optional arguments:
        type            :The type of the particles (e.g. 'star')
        snap            :Number of the snapshot. Defaults to 'current'
        sim             :Number of the simulation. Defaults to 'current'
        unit            :Specifies the unit to use to return the data
    '''
    simno = get_sim_no(sim)
    sim = SimBuffer.get_sim_no(simno)
    snapobject = SimBuffer.get_snapshot_extended(sim, snap)
    nspecies = snapobject.GetNTypes()
    if type=="all":
        raise Exception("You requested all particle types to get_data, but we can return only one array!")
    fetcher=UserQuantity(quantity)
    unitinfo,data,scaling,label=fetcher.fetch(type=type,snap=snapobject,unit=unit)
    return data*scaling

def get_render_data(x,y,quantity, sim="current",snap="current",
                    renderunit="default",
                    res=64,zslice=None,coordlimits=None):
    '''Return the rendered data for the given quantity. Useful when one needs
    to grid SPH data. The result is scaled to the specified unit. The options are 
    a subset of the options available to the 'render' function.
    
    Required arguments:
        x        : Quantity on the x-axis. Must be a string
        y        : Quantity on the y-axis. Must be a string
        quantity : Quantity to render.
        
    Optional arguments:
        snap     : Number of the snapshot to plot. Defaults to 'current'.
        sim      : Number of the simulation to plot. Defaults to 'current'
        renderunit: Unit to use for the rendered quantity
        res      : Resolution
        zslice   : z-coordinate of the slice when doing a slice rendering.
                   Default is None, which produces a column-integrated plot.
                   If you set this variable, a slice rendering will be
                   done instead.
        coordlimits: Limits of the coordinates on x and y. See documentation
                     of render.
                     
    Return:
        data     : The rendered data, scaled to the requested unit.
        
    
    '''
    if zslice is not None:
        zslice = float(zslice)
    simno=get_sim_no(sim)
    if coordlimits is not None and isinstance(coordlimits, types.StringTypes):
        coordlimits = to_list (coordlimits, float)
    if isinstance(res, types.StringTypes):
        if res[0]=='[' or res[0]=='(':
            res = to_list(res,int)
        else:
            res = int(res)
    command = Commands.RenderPlotCommand(x, y, quantity, snap, simno, True,
                                         True, True,
                                         coordlimits, zslice, "default", "default",
                                         renderunit, res, "nearest")
    data = command.prepareData(Singletons.globallimits)
    return data.render_data   
    


#------------------------------------------------------------------------------
def to_list(str_variable,type):
    '''Convert the input string to a list of the specified type'''
    parenthesis_open = ('[', '(')
    parenthesis_closed = (']',')')
    if str_variable[0] not in parenthesis_open or str_variable[-1] not in parenthesis_closed:
        raise ValueError('What you passed cannot be parsed as a tuple')

    splitted = str_variable[1:-1].split(',')
    return map(type,splitted)


#------------------------------------------------------------------------------
def to_bool(value):
    '''Parses the input string and convert it to a boolean. If the input is
    not a string, passes it to the built-in bool function (which means, that
    the result is False only if it is None or False).'''
    valid = {'true': True, 't': True, '1': True,
             'false': False, 'f': False, '0': False,
             }

    if not isinstance(value, types.StringTypes):
        return bool(value)

    lower_value = value.lower()
    if lower_value in valid:
        return valid[lower_value]
    else:
        raise ValueError('invalid literal for boolean: "%s"' % value)


#------------------------------------------------------------------------------
def sigint(signum, frame):
    cleanup()


#------------------------------------------------------------------------------
def cleanup():
    Singletons.place_command(["STOP",None])
    print "Waiting for background processes to finish..."
    plottingprocess.join()
    import sys
    sys.exit()

def ListFunctions():
    '''List the available functions defined in facade'''
    import gandalf_interpreter
    toexcludefunctions=gandalf_interpreter.toexcludefunctions
    functions = inspect.getmembers(facade, inspect.isfunction)
    functions=filter(lambda function: function not in toexcludefunctions, functions)
    print "The available functions in facade are: "
    for function in functions:
        print function.__name__

#------------------------------------------------------------------------------
def init():
    if defaults.parallel:
        global plottingprocess
        plottingprocess = PlottingProcess(Singletons.queue, Singletons.commands, Singletons.completedqueue, Singletons.globallimits, Singletons.free)
        plottingprocess.start()
    else:
        global plotting
        plotting=Plotting()
    CreateUserQuantity('r','sqrt(x^2+y^2+z^2)',scaling_factor='r', label='$r$')
    CreateUserQuantity('R','sqrt(x^2+y^2)',scaling_factor='r', label='$R$')
    CreateUserQuantity('phi','arctan2(y,x)', label='$\\phi$')
    CreateUserQuantity('theta','arccos(z/r)', label='$\\theta$')
    CreateUserQuantity('vr','sin(theta)*cos(phi)*vx+sin(theta)*sin(phi)*vy+cos(theta)*vz',scaling_factor='v',label='$v_r$')
    CreateUserQuantity('vR','sin(theta)*cos(phi)*vx+sin(theta)*sin(phi)*vy',scaling_factor='v',label='$v_R$')
    CreateUserQuantity('vphi','cos(phi)*vy-sin(phi)*vx',scaling_factor='v',label='$v_\\phi$')
    CreateUserQuantity('vtheta','cos(theta)*cos(phi)*vx+cos(theta)*sin(phi)*vy-sin(theta)*vz',scaling_factor='v', label='$v_\\theta$')
    CreateUserQuantity('ar','sin(theta)*cos(phi)*ax+sin(theta)*sin(phi)*ay+cos(theta)*az',scaling_factor='a',label='$a_r$')
    CreateUserQuantity('aR','sin(theta)*cos(phi)*ax+sin(theta)*sin(phi)*ay',scaling_factor='a',label='$a_R$')
    CreateUserQuantity('aphi','cos(phi)*vy-sin(phi)*vx',scaling_factor='a',label='$a_\\phi$')
    CreateUserQuantity('atheta','cos(theta)*cos(phi)*vx+cos(theta)*sin(phi)*vy-sin(theta)*vz',scaling_factor='a', label='$a_\\theta$')
    CreateUserQuantity('press','(gamma_eos - 1)*rho*u',scaling_factor='press',label='$P$')
    CreateUserQuantity('sound','sqrt(gamma_eos*(gamma_eos - 1)*u)',scaling_factor='v', label='$c_s$')
    CreateUserQuantity('temp','(gamma_eos - 1)*u*mu_bar',scaling_factor='temp',label='T')

    from data_fetcher import get_time_snapshot
    CreateTimeData('t',get_time_snapshot)

    from compute import COM
    CreateTimeData('com_x',COM)
    CreateTimeData('com_y',COM,quantity='y')
    CreateTimeData('com_z',COM,quantity='z')
    CreateTimeData('com_vx',COM,quantity='vx')
    CreateTimeData('com_vy',COM,quantity='vy')
    CreateTimeData('com_vz',COM,quantity='vz')



# Default code.  Run when facade.py is imported
#------------------------------------------------------------------------------
init()

if defaults.parallel:
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
