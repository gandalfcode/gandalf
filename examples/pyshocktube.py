#==============================================================================
# pyshocktube.py
# Example of setting-up a 1D shocktube simulation inside a python script.
#==============================================================================
from seren.analysis.facade import *
import numpy as np
import time

# Create a new (but unspecified) 1D simulation object
sim = newsim(ndim=1)

# Set all important code parameters
sim.SetParam('ic','python')
sim.SetParam('sim','sph')
sim.SetParam('sph','gradh')
sim.SetParam('run_id','SHOCK1')
sim.SetParam('sph_integration','lfkdk')
sim.SetParam('gas_eos','energy_eqn')
sim.SetParam('energy_integration','PEC')
sim.SetParam('gamma_eos',1.6666666666666666666)
sim.SetParam('Npart',200)
sim.SetParam('tend',0.075)
sim.SetParam('dt_snap',1.0)
sim.SetParam('x_boundary_lhs','open')
sim.SetParam('x_boundary_rhs','open')
sim.SetParam('dimensionless',1)

# Set number of particles and allocate local numpy arrays
Npart = 200
x = np.linspace(-0.995,0.995,num=Npart)
vx = np.zeros(Npart)
m = np.ones(Npart)*2.0/Npart
u = np.ones(Npart)*1.5

# Loop over all particles and set properties, depending on whether the 
# particle is in the LHS shock region or the RHS shock region
for i in range(Npart):
    if x[i] < 0.0:
        vx[i] = 4.0
    else:
        vx[i] = -4.0

# Now call setup routines to create objects and allocate memory 
# (N.B. parameters are now fixed and cannot be changed after this point)
sim.PreSetupForPython()

# Import all arrays to main memory
sim.ImportArray(x,'x')
sim.ImportArray(vx,'vx')
sim.ImportArray(m,'m')
sim.ImportArray(u,'u')
sim.SetupSimulation()

# First, plot the initial conditions
plot("x","rho")

# Now run the simualtion and plot the results in a new window
run()
plot("x","rho")
block()
