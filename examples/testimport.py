from seren.analysis.facade import *

import numpy as np

sim=newsim('testimport.dat')
sim.PreSetupForPython()
x=np.arange(10)
xx, yy = np.meshgrid(x,x)
sim.ImportArray(xx.flatten(),'x')
sim.ImportArray(yy.flatten(),'y')
sim.ImportArray(np.ones(100)*5,'m')
sim.ImportArray(np.ones(100),'u')
sim.SetupSimulation()
plot('x','y')
import time; time.sleep(1)
run()
window(2)
render("x","y","rho")
block()

