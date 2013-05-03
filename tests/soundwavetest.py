#==============================================================================
# soundwavetest.py
# Run isothermal soundwave test using initial conditions file 'soundwave.dat'
#==============================================================================
from seren.analysis.facade import *
import time

# Create new simulation object and generate initial conditions
newsim('soundwave.dat')

# Plot density distribution along with analytical solution
plot("x","rho")
plotanalytical("x","rho")

# Use 'sleep' hack to refresh plot window.  Then run simulation.
time.sleep(1)
run()
block()
