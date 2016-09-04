#==============================================================================
# quadrupletest.py
# Run the binary star test using initial conditions specified in the file 
# 'quadruple.dat' and then plot the orbital quantities.
#==============================================================================
from gandalf.analysis.facade import *
import time

# Create new shocktube simulation from 'adshock.dat' file
newsim('quadruple.dat')
setupsim()

# Plot the density with the analytical solution
plot("x","y",type="star")
limit("x",-1.5,1.5,window="all")
limit("y",-1.5,1.5,window="all")

# 'Sleep hack' (to allow matplotlib to update the figure) 
# before running the simulation
time.sleep(1)
run()
block()
