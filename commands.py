from data import Data
from SimBuffer import SimBuffer

class Command:
    
    id = 0
    
    def __init__(self):
        #sets the unique identifier of the command
        Command.id += 1
        self.id = Command.id

class WindowCommand(Command):
    def __init__(self):
        Command.__init__(self)
        
    def processCommand (self, plotting, data):
        fig = plotting.plt.figure()
        fig.show()
        
class PlotCommand(Command):
    def __init__(self, xquantity, yquantity, autoscale = True):
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
            if self.overplot:
                ax = fig.gca()
            else:
                fig.clear()
                ax = fig.add_subplot(111)
            product = self.execute(plotting, fig, ax, data)
            fig.show()
            plotting.commands.append(self)
            plotting.commandsfigures[self.id]= fig, ax, product
            plotting.lastid = self.id

class ParticlePlotCommand (PlotCommand):
    
    def __init__(self, xquantity, yquantity):
        PlotCommand.__init__(self, xquantity, yquantity)
        
    def update(self, plotting, fig, ax, line, data):
        line.set_data(data.x_data,data.y_data)
        
    def execute(self, plotting, fig, ax, data) :
        line, = ax.plot(data.x_data, data.y_data, '.')
        return line
    
    def prepareData (self):
        snap = SimBuffer.get_current_snapshot()
        x_data = snap.ExtractArray(self.xquantity)
        y_data = snap.ExtractArray(self.yquantity)
        data = Data(x_data, y_data)
        return data