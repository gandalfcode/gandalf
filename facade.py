import atexit
import commands
from multiprocessing import Manager, Queue
from plotting import PlottingProcess
from SimBuffer import SimBuffer, BufferException
import signal

class Singletons:
    queue = Queue()
    commands = Manager().list()
    
def loadsim(run_id):
    SimBuffer.loadsim(run_id)
    return SimBuffer.get_current_sim()
    
def plot(x,y, overplot = False, snap="current", autoscale = True):
    command = commands.ParticlePlotCommand(x, y, autoscale)
    command.overplot = overplot
    command.snap = snap
    data = command.prepareData()
    Singletons.queue.put([command, data])

def addplot (x,y, **kwargs):
    plot(x,y, True, **kwargs)
    
def next():
    try:
        snap(SimBuffer.get_no_next_snapshot())
    except BufferException as e:
        print str(e)
        
def previous():
    try:
        snap(SimBuffer.get_no_previous_snapshot())
    except BufferException as e:
        print str(e)
        
def snap(no):
    try:
        SimBuffer.set_current_snapshot_number(no)
    except BufferException as e:
        print str(e)
        return
    for command in Singletons.commands:
        if command.snap == "current": 
            data = command.prepareData()
            Singletons.queue.put([command, data])
        
def window(no = None):
    command = commands.WindowCommand(no)
    data = None
    Singletons.queue.put([command,data])

def subfigure(nx, ny, current):
    command = commands.SubfigureCommand(nx, ny, current)
    data = None
    Singletons.queue.put([command,data])

def init():
    global plottingprocess
    plottingprocess = PlottingProcess(Singletons.queue, Singletons.commands)
    plottingprocess.start()
    
init()

def sigint(signum, frame):
    cleanup()
    
def cleanup():
    Singletons.queue.put(["STOP",None])
    plottingprocess.join()
    import sys
    sys.exit()
 
signal.signal(signal.SIGINT, sigint)
signal.signal(signal.SIGTERM, sigint)
signal.signal(signal.SIGSEGV, sigint)
atexit.register(cleanup)

if __name__=="__main__":
    import time; time.sleep(1)
    loadsim('TEST')
    plot("x","y", snap=0)
    addplot("x", "y")
    window()
    plot("vx", "vy")
    plot("vx", "x")
    window()
    plot("x","rho")
    window()
    subfigure(2,2,1)
    plot("x", "y")
    subfigure(2,2,2)
    plot("vx", "vy")
    subfigure(2,2,3)
    plot("x", "rho")
    subfigure(2,2,4)
    plot("rho", "h")
    addplot("rho", "m")
    window(3)
    addplot("rho", "h")
    snap(200)
    for i in range(10):
        time.sleep(1)
        previous()
    Singletons.queue.put(["STOP", None])
    plottingprocess.join()

