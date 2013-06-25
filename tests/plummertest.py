#==============================================================================
# plummertest.py
# Run the Plummer sphere test using initial conditions specified in the file 
# 'plummer.dat' and then plot various quantities.
#==============================================================================
from seren.analysis.facade import *
import time

# Create new shocktube simulation from 'adshock.dat' file
newsim('plummer.dat')
setupsim()

# Plot the density with the analytical solution
plot("x","y")
limit("x",-5.5,5.5,window="all")
limit("y",-5.5,5.5,window="all")

# 'Sleep hack' (to allow matplotlib to update the figure) 
# before running the simulation
time.sleep(1)
run()
plot("x","y")
block()
