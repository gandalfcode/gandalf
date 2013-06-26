#==============================================================================
# pyicburrau.py
# Example of setting-up the initial conditions for the 3-body
# Burrau/pythagorean problem in python and then running the simulation and
# displaying the results.
#==============================================================================
from seren.analysis.facade import *
import numpy as np
import time

# Create a new (but unspecified) 2D simulation object
sim = newsim(ndim=2)

# Set all important code parameters
sim.SetParam('ic','python')
sim.SetParam('sim','nbody')
sim.SetParam('run_id','BURRAU1')
sim.SetParam('nbody','hermite4ts')
sim.SetParam('sub_systems',0)
sim.SetParam('Npec',5)
sim.SetParam('Nstar',3)
sim.SetParam('Nsph',0)
sim.SetParam('tend',80.0)
sim.SetParam('nbody_mult',0.025)
sim.SetParam('dt_snap',1.0)
sim.SetParam('nlevels',1)
sim.SetParam('x_boundary_lhs','open')
sim.SetParam('x_boundary_rhs','open')
sim.SetParam('y_boundary_lhs','open')
sim.SetParam('y_boundary_rhs','open')
sim.SetParam('dimensionless',1)

# Set number of particles and allocate local numpy arrays
Nstar = 3
x = np.zeros(Nstar)
y = np.zeros(Nstar)
vx = np.zeros(Nstar)
vy = np.zeros(Nstar)
m = np.zeros(Nstar)

# Set all (non-zero) values for star initial conditions
m[0] = 3.0;  x[0] = 1.0;   y[0] = 3.0
m[1] = 4.0;  x[1] = -2.0;  y[1] = -1.0
m[2] = 5.0;  x[2] = 1.0;   y[2] = -1.0

print m
print x
print y

# Now call setup routines to create objects and allocate memory 
# (N.B. parameters are now fixed and cannot be changed after this point)
sim.PreSetupForPython()

# Import all arrays to main memory
sim.ImportArray(x,'x','star')
sim.ImportArray(y,'y','star')
sim.ImportArray(vx,'vx','star')
sim.ImportArray(vy,'vy','star')
sim.ImportArray(m,'m','star')
sim.SetupSimulation()

# First, plot the initial conditions
plot("x","y",type='star')
limit('x',-10.0,10.0)
limit('y',-10.0,10.0)


# Now run the simualtion and plot the results in a new window
time.sleep(1)
run()
plot("x","y",type='star')
addplot("x","y",type='star',snap=0,color='blue')
block()
