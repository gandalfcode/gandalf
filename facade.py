import commands
from multiprocessing import Manager, Queue
from plotting import PlottingProcess
from SimBuffer import SimBuffer


class Singletons:
    queue = Queue()
    commands = Manager().list()
    
def loadsim(run_id):
    SimBuffer.loadsim(run_id)
    return SimBuffer.get_current_sim()
    
def plot(x,y, overplot = False):
    command = commands.ParticlePlotCommand(x, y)
    command.overplot = overplot
    data = command.prepareData()
    Singletons.queue.put([command, data])

def addplot (x,y):
    plot(x,y, True)
    
def next():
    snap(SimBuffer.get_no_next_snapshot())
        
def snap(no):
    SimBuffer.get_snapshot_number_current_sim(no)
    for command in Singletons.commands:
        data = command.prepareData()
        Singletons.queue.put([command, data])
        
def window():
    command = commands.WindowCommand()
    data = None
    Singletons.queue.put([command,data])
    
if __name__=="__main__":
    plottingprocess = PlottingProcess(Singletons.queue, Singletons.commands)
    plottingprocess.start()
    import time; time.sleep(1)
    loadsim('TEST')
    plot("x","y")
    window()
    plot("vx", "vy")
#    plot("vx", "x")
    window()
    plot("x","rho")
    addplot("rho", "h")
    for i in range(10):
        time.sleep(1)
        next()
    Singletons.queue.put(["STOP",None])
    plottingprocess.join()