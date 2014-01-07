#==============================================================================
# example12.py
# Plot particle quantities in an alternative coordinate system.
#==============================================================================
from gandalf.analysis.facade import *

# Create simulation object from Boss-Bodenheimer parameters file
sim = newsim("bossbodenheimer.dat")
sim.SetParam("tend",0.02)
setupsim()

# Run simulation and plot x-y positions of SPH particles in the default
# units specified in the `bossbodenheimer.dat' parameters file.
plot("x","y")
run()
block()

# After pressing return, re-plot last snapshot but in new specified units (au).
plot("x","y",xunit="au",yunit="au")
block()
