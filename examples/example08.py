#==============================================================================
# example08.py
# Create initial conditions for SPH simulation inside the python script, and
# then run the simulation to completion while plotting results.
#==============================================================================
from gandalf.analysis.facade import *
import numpy as np
import time

# Set basic parameters for generating initial conditions
Nhydro = 200
vfluid = 4.0
xmin = -1.5
xmax = 1.5

# Set uniform line of Nhydro particles between the limits of xmin and xmax
# in local numpy arrays
deltax = (xmax - xmin) / Nhydro
x = np.linspace(xmin + 0.5*deltax,xmax - 0.5*deltax,num=Nhydro)
m = np.ones(Nhydro)*(xmax - xmin)/Nhydro

# Set velocities of shock-tube so v = vfluid for x < 0 and -vfluid for x > 0
vx = np.ones(Nhydro)*vfluid
vx[x > 0.0] = -vfluid

# Create new 1D simulation object and set parameters
sim = newsim(ndim=1,sim='sph')
sim.SetParam('ic','python')
sim.SetParam('gas_eos','isothermal')
sim.SetParam('Nhydro',Nhydro)
sim.SetParam('tend',0.2)
sim.SetParam('dt_snap',0.05)
sim.SetParam('dimensionless',1)
sim.SetParam('vfluid1[0]',vfluid)
sim.SetParam('vfluid2[0]',-vfluid)
sim.SetParam('boxmin[0]',xmin)
sim.SetParam('boxmax[0]',xmax)
sim.SetParam('run_id','SHOCKTUBE1')

# Call setup routines and import particle data
sim.PreSetupForPython()
sim.ImportArray(x,'x')
sim.ImportArray(vx,'vx')
sim.ImportArray(m,'m')
setupsim()
# Plot the density of all particles near the shock
plot("x","rho")
plotanalytical("x","rho",ic="shocktube")
limit("x",-0.17,0.17,window="all")
limit("rho",0,21.0,window="all")

# Run simulation and save plot to file
run()
savefig("shocktube.png")
block()
