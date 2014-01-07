#==============================================================================
# example06.py
# Load in an old simulation from the disk, run a new simulation and then
# plot both simulations in the same plot for comparison.
#==============================================================================
from gandalf.analysis.facade import *

# Load in Sod test simulation (run in example 1)
sim0 = loadsim("ADSOD1")

# Create new simulation object from `adsod.dat' parameters file, but modify
# to use no artificial viscosity and then run to completion.
sim1 = newsim("adsod.dat")
sim1.SetParam("run_id","ADSOD2")
sim1.SetParam("avisc","none")
setupsim()
run()

# Plot Sod test results, with (sim=0) and without (sim=1) artificial viscosity,
# on same figure and then save to eps file.
plot("x","rho",sim=0,snap=2)
addplot("x","rho",sim=1,snap=2)
savefig("sod12.eps")
block()
