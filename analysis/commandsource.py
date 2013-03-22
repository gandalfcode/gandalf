import analytical
from data import Data
from facade import SimBuffer
import numpy as np
from swig_generated.SphSim import Render

class Command:
    
    id = 0
    
    def __init__(self):
        #sets the unique identifier of the command
        Command.id += 1
        self.id = Command.id

class SwitchNonGui(Command):
    def __init__(self):
        Command.__init__(self)
        
    def processCommand(self, plotting, data):
        plotting.plt.switch_backend('Agg')
        
class SaveFigCommand(Command):
    def __init__(self, name):
        Command.__init__(self)
        self.name = name
        
    def processCommand (self, plotting, data):
        fig = plotting.plt.savefig(self.name)

class WindowCommand(Command):
    def __init__(self, no = None):
        Command.__init__(self)
        self.no = no
        
    def processCommand (self, plotting, data):
        fig = plotting.plt.figure(self.no)
        fig.show()
        fig.canvas.draw()
        
class SubfigureCommand(Command):
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
        update = plotting.command_in_list(self.id)
        if update:
            fig, ax, product = plotting.commandsfigures[self.id]
            self.update(plotting, fig, ax, product, data)
            self.labels(ax)
            self.autolimits(ax)
            fig.canvas.draw()
        elif self.id > plotting.lastid:
            fig = plotting.plt.gcf()
            ax = fig.gca()
            if not self.overplot:
                ax.clear()
                self.labels(ax)
            product = self.execute(plotting, fig, ax, data)
            #sets the autoscales on the axis
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
            plotting.commands.append(self)
            plotting.commandsfigures[self.id]= fig, ax, product
            plotting.lastid = self.id
            #finally we need to draw, because the previous call to show does not update
            #the figure
            fig.canvas.draw()
            #now we register quantities and units
            plotting.quantitiesfigures[(fig, ax, 'x')] = self.xquantity, self.xunitname
            plotting.quantitiesfigures[(fig, ax, 'y')] = self.yquantity, self.yunitname
            
    def labels(self, ax):
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
        scaling_factor=1.
        quantity = getattr(self, axis+'quantity')
        unit = getattr(self, axis+'unit')
        data, scaling_factor = snap.ExtractArray(quantity, scaling_factor, unit)
        setattr(self,axis+'label',snap.label)
        setattr(self,axis+'unitname',snap.unitname)
        return data, scaling_factor
    
    def setlimits(self, plotting, ax, axis):
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
            analyticalclass = getattr(analytical,ictype)
        except AttributeError:
            raise CommandException, "We do not know an analytical solution for the requested initial condition"
        computer = analyticalclass(sim, time)
        x_data, y_data = computer.compute(self.xquantity, self.yquantity)
        
        dummy, xscaling_factor = self.get_array('x', snap)
        dummy, yscaling_factor = self.get_array('y', snap)
        
        data = Data(x_data*xscaling_factor,y_data*yscaling_factor)
        return data

class RenderPlotCommand (PlotCommand):
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
        
        rendering = Render()
        renderscaling_factor=1.
        rendered = np.zeros(xres*yres, dtype=np.float32)
        if sim.ndim < 3 or self.zslice is None:
            returncode, renderscaling_factor = rendering.CreateColumnRenderingGrid(xres, yres, self.xquantity, self.yquantity, self.renderquantity,
                                                 self.renderunit, self.xmin, self.xmax,
                                                 self.ymin, self.ymax, rendered, snap, sim.sph, renderscaling_factor)
        else:
            quantities = ['x','y','z']
            quantities.pop(quantities.index(self.xquantity))
            quantities.pop(quantities.index(self.yquantity))
            self.zquantity = quantities[0]
            z_data, z_scaling_factor = self.get_array('z', snap)
            returncode, renderscaling_factor = rendering.CreateSliceRenderingGrid(xres, yres, self.xquantity, self.yquantity, self.zquantity, self.renderquantity,
                                                 self.renderunit, self.xmin, self.xmax,
                                                 self.ymin, self.ymax, self.zslice, rendered, snap, sim.sph, renderscaling_factor)
        rendered = rendered.reshape(xres,yres)
        np.set_printoptions(threshold='nan')
#        data = Data(x*xscaling_factor, y*yscaling_factor, rendered*renderscaling_factor)
        data = Data(None, None, rendered*renderscaling_factor)
        
        return data

class LimitCommand(Command):
    def __init__(self, quantity, min, max, auto, window, subfigure):
        Command.__init__(self)
        self.quantity = quantity
        self.min = min
        self.max = max
        self.auto = auto
        self.window = window
        self.subfigure = subfigure

    def prepareData(self, globallimits):
        pass

    def processCommand(self, plotting, data):
        #first retrieve the correct figs and ax objects
        try:
            figs, axs, line = plotting.commandsfigures[self.id]
        except KeyError:
            #this gets executed the first time the command is run
            if self.window == 'current':
                fig = plotting.plt.gcf()
                figs=[fig]
            elif self.window == 'all' or self.window == 'global':
                managers = plotting.plt._pylab_helpers.Gcf.get_all_fig_managers()
                figs=map(lambda manager: manager.canvas.figure, managers)
            else:
                fig = plotting.plt.figure(self.window)
                figs=[fig]
            line = None
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
        
        for fig in figs:
            fig.canvas.draw()

class RescaleCommand(Command):
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
