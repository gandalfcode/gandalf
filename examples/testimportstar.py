from seren.analysis.facade import *

import numpy as np

sim=newsim('testimport.dat')
sim.SetParam('Nstar',100)
sim.PreSetupForPython()
x=np.arange(10)
xx, yy = np.meshgrid(x,x)
sim.ImportArray(xx.flatten(),'x', 'star')
sim.ImportArray(yy.flatten(),'y', 'star')
sim.ImportArray(np.ones(100)*5,'m', 'star')
sim.SetupSimulation()
plot('x','y', type='star')
import time; time.sleep(1)
block()

