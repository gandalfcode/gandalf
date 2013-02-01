from multiprocessing import Process
from Queue import Empty

class PlottingProcess (Process):
    def __init__(self, queue, commands):
        Process.__init__(self)
        self.queue = queue
        self.commands = commands
        self.commandsfigures = {}
        self.lastid = 0
    
    def run(self):
        import matplotlib.pyplot as plt
        self.plt = plt
        
        # Main loop
        while 1:
            self.mypause(0.01)
            
            #Reads data from the queue
            try:
                command, data = self.queue.get(block = False)
            except Empty:
                continue
            if command == "STOP":
                break
            
            self.remove_closed_figures()
            
            command.processCommand(self, data)
                
    def command_in_list (self, id):
        for command in self.commands:
            if id == command.id:
                return True
        return False
    
            
    def mypause (self, interval):
        plt = self.plt
        backend = plt.rcParams['backend']
        if backend in plt._interactive_bk:
            figManager = plt._pylab_helpers.Gcf.get_active()
            if figManager is not None:
                canvas = figManager.canvas
                canvas.start_event_loop(interval)
                return
    
        # No on-screen figure is active, so sleep() is all we need.
        import time
        time.sleep(interval)
    
    def remove_closed_figures(self):
        for index, commanddummy in enumerate(self.commands):
            fig, ax, line = self.commandsfigures[commanddummy.id]
            if not self.plt.fignum_exists(fig.number):
                self.commandsfigures.pop(commanddummy.id)
                self.commands.pop(index)
