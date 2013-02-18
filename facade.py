import atexit
import commandsource as Commands
from multiprocessing import Manager, Queue
from plotting import PlottingProcess
from SimBuffer import SimBuffer, BufferException
import signal
from time import sleep

class Singletons:
    queue = Queue()
    commands = Manager().list()
    
def loadsim(run_id):
    SimBuffer.loadsim(run_id)
    return SimBuffer.get_current_sim()
    
def plot(x,y, overplot = False, snap="current", autoscale = True, sim="current"):
    command = Commands.ParticlePlotCommand(x, y, autoscale)
    command.overplot = overplot
    command.snap = snap
    if sim == "current":
        simno = SimBuffer.get_current_sim_no()
    else:
        simno = sim
    #TODO: substitute the number of the simulation inside the plot command with the object
    command.sim = simno
    data = command.prepareData()
    Singletons.queue.put([command, data])
    sleep(0.001)

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
    
    update("current")
        
def window(no = None):
    command = Commands.WindowCommand(no)
    data = None
    Singletons.queue.put([command,data])

def subfigure(nx, ny, current):
    command = Commands.SubfigureCommand(nx, ny, current)
    data = None
    Singletons.queue.put([command,data])

def newsim(paramfile):
    SimBuffer.newsim(paramfile)

def run(no=None):
    #gets the correct simulation object from the buffer
    try:
        if no is None:
            sim = SimBuffer.get_current_sim()
        else:
            sim = SimBuffer.get_sim_no(no)
    except BufferError as e:
        print str(e)
        return
    if not sim.setup:
        print "The selected simulation has not been set-up. Please set it up before running"
        return
    
    sim.Run()
    
    SimBuffer.load_live_snapshot(sim)
    
    update("live")

def block():
    print "Press enter to quit..."
    raw_input()

def update(type=None):
    #updates the plots
    for command in Singletons.commands:
        updateplot=False
        if type is None:
            updateplot=True
        else:
            if command.snap == type:
                updateplot=True
        if updateplot:
            data = command.prepareData()
            Singletons.queue.put([command, data])

def savefig(name):
    command = Commands.SaveFigCommand(name)
    data = None
    Singletons.queue.put([command,data])
    time.sleep(1e-3)

def switch_nongui():
    command = Commands.SwitchNonGui()
    data = None
    Singletons.queue.put([command,data])
    time.sleep(1e-3)

def plotanalytical(x=None, y=None, overplot = True, sim = "current", snap = "current", autoscale = True):
    '''Plots the analytical solution'''
    
    #TODO: remove duplicated code from the plot function
    
    #TODO: figure out automatically the quantities to plot depending on current window    
    
    #get the simulation number from the buffer
    if sim == "current":
        simno = SimBuffer.get_current_sim_no()
    else:
        simno = sim
    
    #istantiate and setup the command object
    command = Commands.AnalyticalPlotCommand(x, y, autoscale)
    command.overplot = overplot
    command.sim = simno
    command.snap = snap
    data = command.prepareData()
    Singletons.queue.put([command, data])

def init():
    global plottingprocess
    plottingprocess = PlottingProcess(Singletons.queue, Singletons.commands)
    plottingprocess.start()
    
init()

def sigint(signum, frame):
    cleanup()
    
def cleanup():
    Singletons.queue.put(["STOP",None])
    print "Waiting for background processes to finish..."
    plottingprocess.join()
    import sys
    sys.exit()
 
signal.signal(signal.SIGINT, sigint)
signal.signal(signal.SIGTERM, sigint)
signal.signal(signal.SIGSEGV, sigint)
atexit.register(cleanup)

if __name__=="__main__":
    loadsim('TEST')
    plot("x","rho")
    plotanalytical("x","rho")
    import time; time.sleep(1)
    next(); time.sleep(1)
    snap(10)
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

