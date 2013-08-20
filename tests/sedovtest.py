#==============================================================================
# sedovtest.py
# Run the Sedov blastwave test using initial conditions specified in the 
# file 'sedov.dat'.
#==============================================================================
from gandalf.analysis.facade import *
import time

# Create new Sedov test simulation object from 'sedov.dat' file
newsim('sedov.dat')
setupsim()

# Run simulation and plot results to screen
plot("R","rho")
run()
block()


