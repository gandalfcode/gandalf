from data import Data
from SimBuffer import SimBuffer

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
            fig.show()
            plotting.commands.append(self)
            plotting.commandsfigures[self.id]= fig, ax, product
            plotting.lastid = self.id
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
        sim = SimBuffer.simlist[self.sim]       
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