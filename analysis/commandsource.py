import analytical
from data import Data
from facade import SimBuffer
import numpy as np

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
                      'vx': '$v_x$', 'vy': '$v_y$', 'vz': '$v_z$', 'm': 'm',
                      'h': 'h', 'u': 'u'}
    
    def __init__(self, xquantity, yquantity, snap, simno, 
                 overplot, autoscale, xunit="default", yunit="default"):
        Command.__init__(self)
        self.xquantity = xquantity
        self.yquantity = yquantity
        self.snap = snap
        self.simno = simno
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
            if self.autoscale:
                self.autolimits(ax)
            fig.canvas.draw()
        elif self.id > plotting.lastid:
            fig = plotting.plt.gcf()
            ax = fig.gca()
            if not self.overplot:
                ax.clear()
                self.labels(ax)
            product = self.execute(plotting, fig, ax, data)
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
            plotting.quantitiesfigures[(fig, 'x')] = self.xquantity, self.xunitname
            plotting.quantitiesfigures[(fig, 'y')] = self.yquantity, self.yunitname
            
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
        ax.relim()
        ax.autoscale_view()
        
    def get_sim_and_snap(self):
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

class ParticlePlotCommand (PlotCommand):
    
    def __init__(self, xquantity, yquantity, snap, simno, overplot=False, 
                 autoscale=True, xunit="default", yunit="default"):
        PlotCommand.__init__(self, xquantity, yquantity, snap, simno, overplot, 
                             autoscale, xunit, yunit)
                
    def update(self, plotting, fig, ax, line, data):
        line.set_data(data.x_data,data.y_data)
        
    def execute(self, plotting, fig, ax, data) :
        line, = ax.plot(data.x_data, data.y_data, '.')
        if self.autoscale:
            self.autolimits(ax)
        return line
    
    def prepareData (self):
        sim, snap = self.get_sim_and_snap()
        
        xscaling_factor=1.
        x_data, xscaling_factor = snap.ExtractArray(self.xquantity, xscaling_factor, self.xunit)
        self.xlabel = snap.label
        self.xunitname = snap.unitname
        yscaling_factor=1.
        y_data, yscaling_factor = snap.ExtractArray(self.yquantity, yscaling_factor, self.yunit)
        self.ylabel = snap.label
        self.yunitname = snap.unitname
        
        data = Data(x_data*xscaling_factor, y_data*yscaling_factor)
        return data
    
class AnalyticalPlotCommand (PlotCommand):
    
    def __init__(self, xquantity, yquantity, snap, simno, overplot, autoscale):
        PlotCommand.__init__(self, xquantity, yquantity, snap, simno, overplot, autoscale)
        
    def update(self, plotting, fig, ax, line, data):
        line.set_data(data.x_data, data.y_data)
        
    def execute(self, plotting, fig, ax, data):
        line, = ax.plot(data.x_data, data.y_data)
        return line
    
    def prepareData(self):
        #TODO: remove duplicated code with the previous class
        sim, snap = self.get_sim_and_snap()
        time = snap.t
        ictype = sim.simparams.stringparams["ic"]
        try:
            analyticalclass = getattr(analytical,ictype)
        except AttributeError:
            raise CommandException, "We do not know an analytical solution for the requested initial condition"
        computer = analyticalclass(sim, time)
        x_data, y_data = computer.compute(self.xquantity, self.yquantity)
        data = Data(x_data,y_data)
        return data

class RenderPlotCommand (RenderCommand):
    
    def __init__(self, xquantity, yquantity, renderquantity, snap, simno, overplot, autoscale, 
                 xunit="default", yunit="default", renderunit="default", res=256):
        PlotCommand.__init__(self, xquantity, yquantity, snap, simno, 
                             overplot, autoscale, xunit, yunit)
        self.renderquantity = renderquantity
        self.renderunit = renderunit
        self.res = res
        
    def update(self, plotting, fig, ax, data):
        im.set_array(data.render_data)
    
    def execute(self, plotting, fig, ax, data):
        im, = ax.imshow(data.render_data)
        return im
    
    def prepareData(self):
        sim, snap = self.get_sim_and_snap()
        #TODO: refactor this code to have a function
        xscaling_factor=1.
        x_data, xscaling_factor = snap.ExtractArray(self.xquantity, xscaling_factor, self.xunit)
        self.xlabel = snap.label
        self.xunitname = snap.unitname
        yscaling_factor=1.
        y_data, yscaling_factor = snap.ExtractArray(self.yquantity, yscaling_factor, self.yunit)
        self.ylabel = snap.label
        self.yunitname = snap.unitname
        renderscaling_factor=1.
        render_data, renderscaling_factor = snap.ExtractArray(self.renderquantity, renderscaling_factor, self.renderunit)
        self.renderlabel = snap.label
        self.renderunitname=snap.unitname
        
        #creates the grid
        #TODO: use proper limits
        #how to get the limits? that's a bit tricky
        #for now, just use the min and max of the array
        try:
            xres = self.res[0]
            yres = self.res[1]
        except TypeError:
            xres = self.res
            yres = self.res
        x = np.linspace(x_data.min(), x_data.max(), xres)
        y = np.linspace(y_data.min(), y_data.max(), yres)
        xx, yy = np.meshgrid(x, y)
        
        #TODO: add in the wrapper; put correct arguments
        rendering = Render()
        rendered = rendering.CreateRenderingGrid()
        
        data = Data(x*xscaling_factor, y*yscaling_factor, rendered*renderscaling_factor)
        
        return data

class LimitCommand(Command):
    def __init__(self, quantity, min, max):
        Command.__init__(self)
        self.quantity = quantity
        self.min = min
        self.max = max

    def prepareData(self):
        pass

    def processCommand(self, plotting, data):
        try:
            fig, ax, line = plotting.commandsfigures[self.id]
        except KeyError:
            fig = plotting.plt.gcf()
            ax = plotting.plt.gca()
            line = None
            plotting.commands.append(self)
            plotting.commandsfigures[self.id]= fig, ax, line
            plotting.lastid = self.id
        if self.quantity == 'x':
            if self.min == 'auto':
                ax.relim()
                ax.autoscale()
            else:
                ax.set_xlim(self.min, self.max)
        elif self.quantity == 'y':
            if self.min == 'auto':
                ax.relim()
                ax.autoscale()
            else:
                ax.set_ylim(self.min, self.max)
        fig.canvas.draw()

class RescaleCommand(Command):
    def __init__(self, quantity, unitname, window):
        Command.__init__(self)
        self.quantity = quantity
        self.unitname = unitname
        self.window = window
        
    def prepareData(self):
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
                    quantity, unitname = plotting.quantitiesfigures[(fig, axis)]
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
            plotting.completedqueue.put("completed")
        
class CommandException(Exception):
    pass     