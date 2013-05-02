import analytical
from data import Data
from facade import SimBuffer
import numpy as np
from swig_generated.SphSim import RenderBase

'''This module contains all the source code for the commands that
are used to make the main process communicate with the plotting process.
Most of them correspond to functions defined in facade.
If you want to create your own command, remember to inherit from Command
and call the command constructor.
Each command must implement a processCommand function, which is the one called by
the plotting process to execute it.
'''


class Command:
    '''Base class from which all the other commands inherit.
    Each command has an unique id, which gets assigned by this class.
    Therefore, it's crucial to call the constructor of this class when
    creating a command.
    '''
    id = 0
    
    def __init__(self):
        #sets the unique identifier of the command
        Command.id += 1
        self.id = Command.id

class SwitchNonGui(Command):
    '''Command for switching to a non gui backend.
    The implementation uses a pyplot function.'''
    def __init__(self):
        Command.__init__(self)
        
    def processCommand(self, plotting, data):
        plotting.plt.switch_backend('Agg')
        
class SaveFigCommand(Command):
    '''Command for saving the current figure. The implementation
    uses a pyplot function'''
    def __init__(self, name):
        Command.__init__(self)
        self.name = name
        
    def processCommand (self, plotting, data):
        fig = plotting.plt.savefig(self.name)

class WindowCommand(Command):
    '''Thin wrapper around the figure function
    in pyplot. Also forces the figure to be drawn.'''
    def __init__(self, no = None):
        Command.__init__(self)
        self.no = no
        
    def processCommand (self, plotting, data):
        fig = plotting.plt.figure(self.no)
        fig.show()
        fig.canvas.draw()
        
class SubfigureCommand(Command):
    '''Thin wrapper around the subplot function
    in pyplot. Also forces the figure to be drawn.'''
    def __init__(self, nx, ny, current):
        Command.__init__(self)
        self.nx = nx
        self.ny = ny
        self.current = current
        
    def processCommand (self, plotting, data):
        ax = plotting.plt.subplot(self.nx, self.ny, self.current)
        fig = ax.figure 
        fig.show()
        fig.canvas.draw()
         
class PlotCommand(Command):
    '''Base class for all the plots.
    A plot must implement execute, which describes how to plot itself the first time,
    update, which describes how to redraw a plot, and prepare data, which describes
    how to extract the data from a simulation object.
    Most of the work for the bookkeeping of figures etc. is done inside the processCommand
    function, which calls execute or update depending on the plot already existing or not.
    The class also contains various helper functions that are needed to do the plots.
    '''
    
    quantitylabels = {'x': 'x', 'y': 'y', 'z': 'z', 'rho': '$\\rho$',
                      'vx': '$v_x$', 'vy': '$v_y$', 'vz': '$v_z$', 
                      'ax': '$a_x$', 'ay': '$a_y$', 'az': '$a_z$',
                      'm': 'm', 'h': 'h', 'u': 'u'}
    
    def __init__(self, xquantity, yquantity, snap, simno, 
                 overplot, autoscale, xunit="default", yunit="default"):
        Command.__init__(self)
        self.xquantity = xquantity
        self.yquantity = yquantity
        self.snap = snap
        self.sim = simno
        self.overplot = overplot
        self.autoscale = autoscale
        self.xunit = xunit
        self.yunit = yunit
        self.xlabel=""
        self.ylabel=""
        self.xunitname = ""
        self.yunitname = ""
        
    def processCommand(self, plotting, data):
        #work out if this is the first time or if the plot already exists          
        update = plotting.command_in_list(self.id)
        if update:
            #the plot already exist - retrieve the relevant objects...
            fig, ax, product = plotting.commandsfigures[self.id]
            #update the plot, the labels, the limits and redraw
            #remember: update is actually defined in the child class
            self.update(plotting, fig, ax, product, data)
            self.labels(ax)
            self.autolimits(ax)
            fig.canvas.draw()
        elif self.id > plotting.lastid:
            #firs time - get current figure and axis
            fig = plotting.plt.gcf()
            ax = fig.gca()
            if not self.overplot:
                #if we are not overplotting, clear what is there
                ax.clear()
                self.labels(ax)
            #call the function in the child class to do the plot
            product = self.execute(plotting, fig, ax, data)
            #set the autoscales on the axis
            if bool(self.autoscale):
                if self.autoscale=='x':
                    ax.set_autoscalex_on(True)
                    self.setlimits(plotting, ax, 'y')
                elif self.autoscale=='y':
                    ax.set_autoscaley_on(True)
                    self.setlimits(plotting, ax, 'x')
                else:
                    ax.set_autoscale_on(True)
            else:
                ax.set_autoscale_on(False)
                self.setlimits(plotting, ax, 'x')
                self.setlimits(plotting, ax, 'y')
            self.autolimits(ax)
            #the figure might have been shown already, but we need to redo it,
            #because (apparently) matplotlib does not provide any way of knowing
            #if that's the case
            fig.show()
            
            #associate the command with the figure, axis and product objects
            plotting.commands.append(self)
            plotting.commandsfigures[self.id]= fig, ax, product
            plotting.lastid = self.id
            
            #finally we need to draw, because the previous call to show does not update
            #the figure
            fig.canvas.draw()
            
            #register quantities and units
            plotting.quantitiesfigures[(fig, ax, 'x')] = self.xquantity, self.xunitname
            plotting.quantitiesfigures[(fig, ax, 'y')] = self.yquantity, self.yunitname
            
    def labels(self, ax):
        '''Write the labels on the x and y axes.
        Uses the quantitylabels dictionary for that.'''
        xlabel = self.quantitylabels[self.xquantity]
        if self.xlabel != "":
            xlabel += ' [$'+self.xlabel+'$]'
        ax.set_xlabel(xlabel)
        ylabel = self.quantitylabels[self.yquantity]
        if self.ylabel != "":
            ylabel += ' [$'+self.ylabel+'$]'
        ax.set_ylabel(ylabel)
        
    def autolimits(self, ax):
        '''Recomputes limits from the data and use them to update
        the limits, if some axis has autoscaling set. Note that if the
        axes have autoscaling off, nothing happens.'''
        ax.relim()
        ax.autoscale_view()
        
    def get_sim_and_snap(self):
        '''Retrieves from the buffer the desired sim and snap'''
        sim = SimBuffer.get_sim_no(self.sim)
        
        if not sim.setup:
            raise Exception("""Error: this simulation has not been set up! If you have set all the relevant parameters/
            initial conditions, please run the setupsim command to set it up.""")
        
        if self.snap == "current":
            snap = SimBuffer.get_current_snapshot_by_sim(sim)
            if sim.snapshots == []:
                self.snap = "live"
        elif self.snap == "live":
            snap = SimBuffer.get_live_snapshot_sim(sim)
        else:
            snap = SimBuffer.get_snapshot_number_sim(sim, self.snap)
        
        return sim, snap
    
    def get_array(self, axis, snap):
        '''Helper function to get the array of the quantity on the x, y or
        rendered axis.
        
        Inputs:
            axis
                String with the desidered axis ('x', 'y' or 'render')
            snap
                Snapshot object
                
        Output:
            data
                The array with the data, in code units
            scaling_factor
                The scaling factor for unit conversion (multiply data by this
                to get the requested unit)
        '''
        scaling_factor=1.
        quantity = getattr(self, axis+'quantity')
        unit = getattr(self, axis+'unit')
        
        self.check_requested_quantity(quantity, snap)
        
        data, scaling_factor = snap.ExtractArray(quantity, scaling_factor, unit)
        setattr(self,axis+'label',snap.label)
        setattr(self,axis+'unitname',snap.unitname)
        return data, scaling_factor
    
    def check_requested_quantity(self, quantity, snap):
        oned = ('x', 'vx', 'ax')
        twod = ('y', 'vy', 'ay')
        threed = ('z', 'vz', 'az')
        minus3 = quantity in threed and snap.ndim<3
        minus2 = quantity in twod and snap.ndim<2
        if minus3 or minus2:
            raise Exception("Error: you requested the quantity " + quantity + ", but the simulation is only in " + str(snap.ndim) + " dims")
        if self.snap != "live":
            if quantity in ('ax', 'ay', 'az'):
                raise Exception ("Error: accelerations are available only for live snapshots")
        
    def setlimits(self, plotting, ax, axis):
        '''Helper function to set the limits of a plot.
        Takes care of handling autoscaling on/off for different axes.
        '''
        try:
            min, max, auto = plotting.globallimits[getattr(self,axis+'quantity')]
        except KeyError:
            return
        if min is None or max is None:
            method = getattr(ax,'set_autoscale'+axis+'_on')
            method(True)
        else:
            method = getattr(ax, 'set_'+axis+'lim')
            method(min,max)
        

class ParticlePlotCommand (PlotCommand):
    '''Inherited class that does particle plotting.
    Uses the plot method of axis for plotting, saves the line generated
    and calls set_data on it for a fast redrawing.
    '''
    
    def __init__(self, xquantity, yquantity, snap, simno, overplot=False, 
                 autoscale=True, xunit="default", yunit="default"):
        PlotCommand.__init__(self, xquantity, yquantity, snap, simno, overplot, 
                             autoscale, xunit, yunit)
                
    def update(self, plotting, fig, ax, line, data):
        line.set_data(data.x_data,data.y_data)
        
    def execute(self, plotting, fig, ax, data) :
        line, = ax.plot(data.x_data, data.y_data, '.')
        return line
    
    def prepareData (self, globallimits):
        sim, snap = self.get_sim_and_snap()

        x_data, xscaling_factor = self.get_array('x', snap)
        y_data, yscaling_factor = self.get_array('y', snap)
        
        data = Data(x_data*xscaling_factor, y_data*yscaling_factor)
        return data
    
class AnalyticalPlotCommand (PlotCommand):
    '''Inherited class that does the plotting of analytical solutions.
    Like particle plotting, uses the plot method on the axis objects and
    the set_data method on the returned line object.'''
    
    def __init__(self, xquantity, yquantity, snap, simno, overplot=True, autoscale=True,
                 xunit="default", yunit="default"):
        PlotCommand.__init__(self, xquantity, yquantity, snap, simno, overplot, 
                             autoscale, xunit, yunit)
        
    def update(self, plotting, fig, ax, line, data):
        line.set_data(data.x_data, data.y_data)
        
    def execute(self, plotting, fig, ax, data):
        line, = ax.plot(data.x_data, data.y_data)
        return line
    
    def prepareData(self, globallimits):
        sim, snap = self.get_sim_and_snap()
        time = snap.t
        ictype = sim.simparams.stringparams["ic"]
        try:
            #try to find an analytical solution class
            #with the same name as the initial conditions
            #type
            analyticalclass = getattr(analytical,ictype)
        except AttributeError:
            raise CommandException, "We do not know an analytical solution for the requested initial condition"
        
        #instantiate the class responsible for computing the analytical solution
        #at the moment we immediately discard the computer object, but note
        #that in principle we could keep it and pass around
        computer = analyticalclass(sim, time)
        x_data, y_data = computer.compute(self.xquantity, self.yquantity)
        
        dummy, xscaling_factor = self.get_array('x', snap)
        dummy, yscaling_factor = self.get_array('y', snap)
        
        data = Data(x_data*xscaling_factor,y_data*yscaling_factor)
        return data

class RenderPlotCommand (PlotCommand):
    '''Child class that does the rendered plots.
    Uses imshow to do the plotting; to update, calls the set_array method on the returned image.'''
    #TODO: add colormap selection
    def __init__(self, xquantity, yquantity, renderquantity, snap, simno, overplot, autoscale,
                 autoscalerender, coordlimits, zslice=None, xunit="default", yunit="default", 
                 renderunit="default", res=64, interpolation='nearest'):
        PlotCommand.__init__(self, xquantity, yquantity, snap, simno, 
                             overplot, autoscale, xunit, yunit)
        self.renderquantity = renderquantity
        self.autoscalerender = autoscalerender
        self.coordlimits = coordlimits
        self.zslice = zslice
        self.zquantity = None
        self.zunit = "default"
        self.renderunit = renderunit
        self.renderunitname = ""
        self.res = res
        self.interpolation = interpolation
        
    def update(self, plotting, fig, ax, products, data):
        im, colorbar = products
        im.set_array(data.render_data)
        if self.autoscalerender:
            im.autoscale()
        else:
            #set limits
            try:
                min, max, auto = plotting.globallimits[self.renderquantity]
                if auto:
                    im.autoscale()
                else:
                    im.set_clim(min,max)
                    self.autoscalerender = False
            except KeyError:
                try:
                    min=self.zmin
                    max=self.zmax
                    im.set_clim(min,max)
                except AttributeError:
                    im.autoscale()
    
    def execute(self, plotting, fig, ax, data):
        im = ax.imshow(data.render_data, extent=(self.xmin, self.xmax, self.ymin, self.ymax), interpolation=self.interpolation)
        
        #set limits
        if not self.autoscalerender:
            try:
                min, max, auto = plotting.globallimits[self.renderquantity]
                if auto:
                    im.autoscale()
                else:
                    im.set_clim(min,max)
                    self.autoscalerender = False
            except KeyError:
                pass
            
        self.autolimits(ax)
        colorbar = fig.colorbar(im)
        products = (im, colorbar)
        plotting.axesimages[ax]=products
        plotting.quantitiesfigures[(fig, ax, 'render')] = self.renderquantity, self.renderunitname
        return products
    
    def prepareData(self, globallimits):
        sim, snap = self.get_sim_and_snap()
        
        x_data, xscaling_factor = self.get_array('x',snap)
        y_data, yscaling_factor = self.get_array('y', snap)
        #create the grid
        #set resolution
        try:
            xres = self.res[0]
            yres = self.res[1]
        except TypeError:
            xres = self.res
            yres = self.res

        #set limits
        if self.coordlimits is None:
            try:
                self.xmin, self.xmax, auto = globallimits[self.xquantity]
                if auto:
                    self.xmin = float(x_data.min())
                    self.xmax = float(x_data.max())
            except KeyError:
                self.xmin = float(x_data.min())
                self.xmax = float(x_data.max())
            try:
                self.ymin, self.ymax, auto = globallimits[self.yquantity]
                if auto:
                    self.ymin = float(y_data.min())
                    self.ymax = float(y_data.max())
            except KeyError:
                self.ymin = float(y_data.min())
                self.ymax = float(y_data.max())
        else:
            self.xmin=self.coordlimits[0]
            self.xmax=self.coordlimits[1]
            self.ymin=self.coordlimits[2]
            self.ymax=self.coordlimits[3]
        
        #create the rendering object
        rendering = RenderBase.RenderFactory(sim.ndims, sim)
        renderscaling_factor=1.
        #allocate the rendered array
        rendered = np.zeros(xres*yres, dtype=np.float32)
        
        #call column integrated or slice rendering, depending on dimensionality and parameters
        if sim.ndims < 3 or self.zslice is None:
            returncode, renderscaling_factor = rendering.CreateColumnRenderingGrid(xres, yres, self.xquantity, self.yquantity, self.renderquantity,
                                                 self.renderunit, self.xmin, self.xmax,
                                                 self.ymin, self.ymax, rendered, snap,renderscaling_factor)
        else:
            quantities = ['x','y','z']
            quantities.pop(quantities.index(self.xquantity))
            quantities.pop(quantities.index(self.yquantity))
            self.zquantity = quantities[0]
            z_data, z_scaling_factor = self.get_array('z', snap)
            returncode, renderscaling_factor = rendering.CreateSliceRenderingGrid(xres, yres, self.xquantity, self.yquantity, self.zquantity, self.renderquantity,
                                                 self.renderunit, self.xmin, self.xmax,
                                                 self.ymin, self.ymax, self.zslice, rendered, snap, renderscaling_factor)
        rendered = rendered.reshape(xres,yres)
#        data = Data(x*xscaling_factor, y*yscaling_factor, rendered*renderscaling_factor)
        data = Data(None, None, rendered*renderscaling_factor)
        
        return data

class LimitCommand(Command):
    '''This command sets the limits of a plot. It is closely connected to the limit function
    in the facade module. Refer to that to know what it does.'''
    def __init__(self, quantity, min, max, auto, window, subfigure):
        Command.__init__(self)
        self.quantity = quantity
        self.min = min
        self.max = max
        self.auto = auto
        self.window = window
        self.subfigure = subfigure

    def prepareData(self, globallimits):
        #we don't need to do anything
        pass

    def processCommand(self, plotting, data):
        #first retrieve the correct figs and ax objects
        try:
            figs, axs, line = plotting.commandsfigures[self.id]
        except KeyError:
            #this gets executed the first time the command is run
            
            #get relevant figures
            if self.window == 'current':
                fig = plotting.plt.gcf()
                figs=[fig]
            elif self.window == 'all' or self.window == 'global':
                #get all the figures from matplotlib
                managers = plotting.plt._pylab_helpers.Gcf.get_all_fig_managers()
                figs=map(lambda manager: manager.canvas.figure, managers)
            else:
                fig = plotting.plt.figure(self.window)
                figs=[fig]
            line = None
            
            #get all the relevant axes
            axs=[]
            for fig in figs:
                if self.subfigure == 'current':
                    ax = fig.gca()
                    axs.append(ax)
                elif self.subfigure == 'all':
                    axslist = fig.axes
                    for ax in axslist:
                        axs.append(ax)
                else:
                    ax = fig.add_subplot(self.subfigure)
                    axs.append(ax)
                    
            plotting.commands.append(self)
            plotting.commandsfigures[self.id]= figs, axs, line
            plotting.lastid = self.id
            if self.window=='global':
                plotting.globallimits[self.quantity] = self.min, self.max, self.auto
                plotting.completedqueue.put('limit completed')
                
        #loops over all the axes and looks if the desired quantity has been plotted there                
        for ax in axs:
            #this big loop takes care of the x and y axis
            for axis in ['x', 'y']:
                try:
                    quantity, unitname = plotting.quantitiesfigures[(ax.figure, ax, axis)]
                except KeyError:
                    continue
                if quantity == self.quantity:
                    #the quantity has been found - call the relevant matplotlib function that actually sets the limits
                    #we use reflection to get the correct name, depending if we are dealing with x or y axis
                    if self.auto:
                        methodname='set_autoscale'+axis+'_on'
                        method=getattr(ax,methodname)
                        method(self.auto)
                        ax.relim()
                        ax.autoscale_view()
                    else:
                        methodname='set_'+axis+'lim'
                        method = getattr(ax,methodname)
                        method(self.min,self.max)
                        
            #here we take care of the quantity when is rendered
            try:
                quantity, unitname = plotting.quantitiesfigures[(ax.figure, ax, 'render')]
            except KeyError:
                continue
            if quantity == self.quantity:
                #looks for the image
                im, colorbar=plotting.axesimages[ax]
                if self.auto:
                   im.autoscale()
                else:
                    im.set_clim(self.min,self.max) 
                    #retrieve the command
                    commandid = [key for key,value in plotting.commandsfigures.items() if value==(ax.figure, ax, (im, colorbar)) ][0]
                    for index,command in enumerate(plotting.commands):
                        if command.id == commandid:
                            command.autoscalerender=False
                            command.zmin=self.min
                            command.zmax=self.max
                            plotting.commands[index]=command
        
        #redraw the affected figures
        for fig in figs:
            fig.canvas.draw()

class RescaleCommand(Command):
    '''Rescaling command. Refer to the rescale function in facade.
    Need to be modified to work in the same way as limit.'''
    def __init__(self, quantity, unitname, window):
        Command.__init__(self)
        self.quantity = quantity
        self.unitname = unitname
        self.window = window
        
    def prepareData(self, globallimits):
        pass
    
    def processCommand(self, plotting, data):
        #only do something THE FIRST TIME
        if not plotting.command_in_list(self.id):
            #TODO: remove duplicated code with LimitCommand class
            try:
                fig, ax, line = plotting.commandsfigures[self.id]
            except KeyError:
                if self.window=='current':
                    fig = plotting.plt.gcf()
                    ax = plotting.plt.gca()
                else:
                    fig = plotting.plt.figure(self.window)
                    ax = fig.gca()
                line = None
                plotting.commands.append(self)
                plotting.commandsfigures[self.id]= fig, ax, line
                plotting.lastid = self.id
            for axis in ['x', 'y', 'render']:
                try:
                    quantity, unitname = plotting.quantitiesfigures[(fig, ax, axis)]
                except KeyError:
                    continue
                if quantity == self.quantity:
                    for index, command in enumerate(plotting.commands):
                        attrname = axis+'quantity'
                        try:
                            if getattr(command, attrname) == quantity:
                                setattr(command,axis+'unit',self.unitname)
                                plotting.commands[index]=command
                        except AttributeError:
                            pass
            plotting.completedqueue.put("rescale completed")
        
class CommandException(Exception):
    pass     
