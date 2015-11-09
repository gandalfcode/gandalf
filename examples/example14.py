#==============================================================================
# example14.py
# Plot quantities of single particles as a function of time through all
# snapshots in the simulation.
#==============================================================================
from gandalf.analysis.facade import *

# Load in simulation from disk (from example 12)
sim = loadsim("BB1")

# Plot the density of particle `100' as a function of time
time_plot("t","rho",id=100,linestyle="-")
block()

# Now plot x-y track of particle 100 as it moves across computational domain
time_plot("x","y",id=100,linestyle="-")
block()
