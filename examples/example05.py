#==============================================================================
# example05.py
# Example to load in a previously run simulation from the disk, 
#==============================================================================
from gandalf.analysis.facade import *

# Load simulation with run_id `ADSOD1' (run in example 1)
loadsim("ADSOD1")

# Plot x vs. density for first snapshot and save to eps file
plot("x","rho",snap=0)
savefig("snap1.eps")

# Plot x vs. density for second snapshot and save to eps file
plot("x","rho",snap=1)
savefig("snap2.eps")

# Plot x vs. density for first and second snapshots and save to eps file
plot("x","rho",snap=0)
addplot("x","rho",snap=1)
savefig("snap12.eps")
block()
