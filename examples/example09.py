#==============================================================================
# example09.py
# Create initial conditions for pure N-body simulation inside the python
# script, and then run the simulation to completion while plotting results.
#==============================================================================
from gandalf.analysis.facade import *
import numpy as np
import time

# Create empty numpy arrays for setting star initial conditions
Nstar = 3
x = np.zeros(Nstar)
y = np.zeros(Nstar)
vx = np.zeros(Nstar)
vy = np.zeros(Nstar)
m = np.zeros(Nstar)
h = 0.000001*np.ones(Nstar)

# Set values for each star individually (Note all velocities initially zero)
m[0] = 3.0;   x[0] = 1.0;   y[0] = 3.0
m[1] = 4.0;   x[1] = -2.0;  y[1] = -1.0
m[2] = 5.0;   x[2] = 1.0;   y[2] = -1.0

# Create new 1D simulation object and set parameters
sim = newsim(ndim=2,sim='nbody')
sim.SetParam('ic','python')
sim.SetParam('nbody','hermite4ts')
sim.SetParam('sub_systems',0)
sim.SetParam('Npec',3)
sim.SetParam('Nlevels',1)
sim.SetParam('Nstar',Nstar)
sim.SetParam('tend',80.0)
sim.SetParam('dt_snap',1.0)
sim.SetParam('noutputstep',128)
sim.SetParam('ndiagstep',2048)
sim.SetParam('dimensionless',1)
sim.SetParam('run_id','BURRAU1')
sim.SetParam('out_file_form','su')

# Call setup routines and import particle data
sim.PreSetupForPython()
sim.ImportArray(x,'x','star')
sim.ImportArray(y,'y','star')
sim.ImportArray(vx,'vx','star')
sim.ImportArray(vy,'vy','star')
sim.ImportArray(m,'m','star')
sim.ImportArray(h,'h','star')
sim.SetupSimulation()

# Plot the density of all particles near the shock
plot("x","y",type="star")
limit("x",-30.0,30.0,window="all")
limit("y",-20.0,40.0,window="all")

# Run simulation and save plot to file
run()
block()
