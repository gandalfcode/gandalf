import analytical
from data import Data
from facade import SimBuffer

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
    def __init__(self, xquantity, yquantity, autoscale):
        Command.__init__(self)
        self.xquantity = xquantity
        self.yquantity = yquantity
        self.overplot = False
        self.autoscale = autoscale
        
    def processCommand(self, plotting, data):             
        update = plotting.command_in_list(self.id)
        if update:
            fig, ax, product = plotting.commandsfigures[self.id]
            self.update(plotting, fig, ax, product, data)
            if self.autoscale:
                ax.relim()
                ax.autoscale_view()
            fig.canvas.draw()
        elif self.id > plotting.lastid:
            fig = plotting.plt.gcf()
            ax = fig.gca()
            if not self.overplot:
                ax.clear()
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

class ParticlePlotCommand (PlotCommand):
    
    def __init__(self, xquantity, yquantity, autoscale):
        PlotCommand.__init__(self, xquantity, yquantity, autoscale)
        
    def update(self, plotting, fig, ax, line, data):
        line.set_data(data.x_data,data.y_data)
        
    def execute(self, plotting, fig, ax, data) :
        line, = ax.plot(data.x_data, data.y_data, '.')
        return line
    
    def prepareData (self):
        sim = SimBuffer.get_sim_no(self.sim)
        if self.snap == "current":
            snap = SimBuffer.get_current_snapshot_by_sim(sim)
            if sim.snapshots == []:
                self.snap = "live"
        elif self.snap == "live":
            snap = SimBuffer.get_live_snapshot_sim(sim)
        else:
            snap = SimBuffer.get_snapshot_number_sim(sim, self.snap)
        x_data = snap.ExtractArray(self.xquantity)
        y_data = snap.ExtractArray(self.yquantity)
        data = Data(x_data, y_data)
        return data
    
class AnalyticalPlotCommand (PlotCommand):
    
    def __init__(self, xquantity, yquantity, autoscale):
        PlotCommand.__init__(self, xquantity, yquantity, autoscale)
        
    def update(self, plotting, fig, ax, line, data):
        line.set_data(data.x_data, data.y_data)
        
    def execute(self, plotting, fig, ax, data):
        line, = ax.plot(data.x_data, data.y_data)
        return line
    
    def prepareData(self):
        #TODO: remove duplicated code with the previous class
        #get the snapshot
        sim = SimBuffer.get_sim_no(self.sim)
        if self.snap == "current":
            snap = SimBuffer.get_current_snapshot_by_sim(sim)
            if sim.snapshots == []:
                self.snap = "live"
        elif self.snap == "live":
            snap = SimBuffer.get_live_snapshot_sim(sim)
        else:
            snap = SimBuffer.get_snapshot_number_sim(sim, self.snap)
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

class LimitCommand(Command):
    def __init__(self, quantity, min, max):
        Command.__init__(self)
        self.quantity = quantity
        self.min = min
        self.max = max
        self.snap = None

    def prepareData(self):
        pass

    def processCommand(self, plotting, data):
        self.snap = None
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
        
class CommandException(Exception):
    pass     