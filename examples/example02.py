#==============================================================================
# example02.py
# Example to prepare a simulation from a parameters file, modify a parameter,
# then run the simulation to completion.
#==============================================================================
from gandalf.analysis.facade import *

# Create simulation from parameters file and modify a single parameter
sim = newsim("adsod.dat")
sim.SetParam("tend",1.5)

# Process modified parameters, set-up simulation and run to the end.
setupsim()
run()
