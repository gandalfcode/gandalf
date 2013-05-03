#==============================================================================
# adshocktest.py
# Run the adiabatic wall-shock test using initial conditions specified in
# the file 'adshock.dat' and then plotting important quantities together 
# with the analytical solutions.
#==============================================================================
from seren.analysis.facade import *
import time

# Create new shocktube simulation from 'adshock.dat' file
newsim('adshock.dat')

# Plot the density with the analytical solution
subfigure(2,2,1)
plot("x","rho")
plotanalytical("x","rho")
limit("x",-8.5,8.5,window="all")

# Plot the x-velocity with the analytical solution
subfigure(2,2,2)
plot("x","vx")
plotanalytical("x","vx")

# Plot the specific internal energy with the solution
subfigure(2,2,3)
plot("x","u")
plotanalytical("x","u")

# Plot the smoothing length
subfigure(2,2,4)
plot("x","h")

# 'Sleep hack' (to allow matplotlib to update the figure) 
# before running the simulation
time.sleep(1)
run()
block()
