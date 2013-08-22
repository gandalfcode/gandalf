#==============================================================================
# plummertest.py
# Run the Plummer sphere test using initial conditions specified in the file 
# 'plummer.dat' and then plot various quantities.
#==============================================================================
from gandalf.analysis.facade import *
import time

# Create new shocktube simulation from 'adshock.dat' file
newsim('plummer.dat')
setupsim()

# Plot the density with the analytical solution
subfigure(2,2,1)
plot("x","y",type='star')
#limit("x",-5.5,5.5,window="all")
#limit("y",-5.5,5.5,window="all")

subfigure(2,2,2)
plot("x","z",type='star')

subfigure(2,2,3)
plot("r","ar",type='star')

subfigure(2,2,4)
plot("r","vr",type='star')

# 'Sleep hack' (to allow matplotlib to update the figure) 
# before running the simulation
time.sleep(1)
run()
block()
