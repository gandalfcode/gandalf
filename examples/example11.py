#==============================================================================
# example11.py
# Plot particle quantities in an alternative coordinate system.
#==============================================================================
from gandalf.analysis.facade import *

# Create simulation object from Kelvin-Helmholtz parameters file
sim = loadsim("SEDOV1")

# Plot the density as a function of radial position
plot("R","rho",snap=5)
block()
