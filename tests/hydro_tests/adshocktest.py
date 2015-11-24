#==============================================================================
# adshocktest.py
# Run the adiabatic wall-shock test using initial conditions specified in
# the file 'adshock.dat' and then plotting important quantities together
# with the analytical solutions.
#==============================================================================
from gandalf.analysis.facade import *
import time

# Create new shocktube simulation from 'adshock.dat' file
newsim('adshock.dat')
setupsim()

# Plot the density with the analytical solution
subfigure(2,1,1)
plot("x","rho")
limit("rho",0.5,5.0)
plotanalytical("x","rho")
limit("x",-9.5,9.5,window="all")

# Plot the x-velocity with the analytical solution
subfigure(2,1,2)
plot("x","vx")
plotanalytical("x","vx")
limit("x",-9.5,9.5,window="all")

# 'Sleep hack' (to allow matplotlib to update the figure)
# before running the simulation
time.sleep(1)
run()
block()
