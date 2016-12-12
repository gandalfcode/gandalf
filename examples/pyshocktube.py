#==============================================================================
#  pyshocktube.py
#  Example of setting-up a 1D shocktube simulation inside a python script.
#==============================================================================
from gandalf.analysis.facade import *
import numpy as np
import time

# Create a new (but unspecified) 1D simulation object
sim = newsim(ndim=1,sim='gradhsph')

# Set all important code parameters
sim.SetParam('ic','python')
sim.SetParam('sph','gradh')
sim.SetParam('run_id','SHOCK1')
sim.SetParam('sph_integration','lfkdk')
sim.SetParam('gas_eos','energy_eqn')
sim.SetParam('gamma_eos',1.6666666666666666)
sim.SetParam('Nhydro',200)
sim.SetParam('tend',0.075)
sim.SetParam('dt_snap',1.0)
sim.SetParam('boundary_lhs[0]','open')
sim.SetParam('boundary_rhs[0]','open')
sim.SetParam('dimensionless',1)

# Set main initial conditions parameters, plus set parameters needed 
# to plot analytical solutions.
Nhydro = 200
xmin = -2.0
xmax = 2.0
vfluid = 8.0
sim.SetParam('rhofluid1',1.0)
sim.SetParam('rhofluid2',1.0)
sim.SetParam('press1',1.0)
sim.SetParam('press2',1.0)
sim.SetParam('vfluid1[0]',vfluid)
sim.SetParam('vfluid2[0]',-vfluid)
sim.SetParam('boxmin[0]',xmin)
sim.SetParam('boxmax[0]',xmax)
deltax = (xmax - xmin) / Nhydro

# Allocate arrays for particles and set main values
x = np.linspace(xmin + 0.5*deltax,xmax - 0.5*deltax,num=Nhydro)
vx = np.zeros(Nhydro)
m = np.ones(Nhydro)*(xmax - xmin)/Nhydro
u = np.ones(Nhydro)*1.5

# Loop over all particles and set properties, depending on whether the 
# particle is in the LHS shock region or the RHS shock region
for i in range(Nhydro):
    if x[i] < 0.0:
        vx[i] = vfluid
    else:
        vx[i] = -vfluid

# Now call setup routines to create objects and allocate memory 
# (N.B. parameters are now fixed and cannot be changed after this point)
sim.PreSetupForPython()

# Import all arrays to main memory
sim.ImportArray(x,'x')
sim.ImportArray(vx,'vx')
sim.ImportArray(m,'m')
sim.ImportArray(u,'u')
sim.SetupSimulation()

# Plot the density with the analytical solution
subfigure(2,2,1)
plot("x","rho")
plotanalytical("x","rho",ic="shocktube")
limit("x",-0.45,0.45,window="all")

# Plot the x-velocity with the analytical solution
subfigure(2,2,2)
plot("x","vx")
plotanalytical("x","vx",ic="shocktube")
limit("x",-0.45,0.45,window="all")

# Plot the specific internal energy with the solution
subfigure(2,2,3)
plot("x","u")
plotanalytical("x","u",ic="shocktube")
limit("x",-0.45,0.45,window="all")

# Plot the smoothing length
subfigure(2,2,4)
plot("rho","h")

# 'Sleep hack' (to allow matplotlib to update the figure) 
# before running the simulation
time.sleep(1)
run()
block()
