#==============================================================================
# example01.py
# Basic example to run a simulation from a parameters file.
#==============================================================================
from gandalf.analysis.facade import *

# Create simulation object from parameters file
sim = newsim("adsod.dat")

# Perform all set-up routines and then run simulation to completion
setupsim()
run()
