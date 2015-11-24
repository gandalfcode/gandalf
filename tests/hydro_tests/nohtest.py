#==============================================================================
# nohtest.py
# Run Noh test using initial conditions file 'noh.dat'.
#==============================================================================
from gandalf.analysis.facade import *
import time

# Create new Noh test simulation object from 'noh.dat' file
newsim('noh.dat')
setupsim()

# Plot density distribution along with analytical solution
plot("R","rho")
plotanalytical("R","rho")
limit("R",0.0,0.4,window="all")
limit("rho",0.0,17.0,window="all")


run()
block()
