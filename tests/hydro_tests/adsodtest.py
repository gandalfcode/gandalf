#==============================================================================
# adsodtest.py
# Run the adiabatic Sod shocktube test using initial conditions specified in
# the file 'adsod.dat' and then plotting important quantities together 
# with the analytical solutions.
#==============================================================================
from gandalf.analysis.facade import *
import time

# Create new shocktube simulation from 'adsod.dat' file
newsim('adsod.dat')
setupsim()

# Plot the density with the analytical solution
subfigure(2,2,1)
plot("x","rho")
plotanalytical("x","rho")
limit("x",-1.0,1.0,window="all")

# Plot the x-velocity with the analytical solution
subfigure(2,2,2)
plot("x","vx")
plotanalytical("x","vx")
limit("x",-1.0,1.0,window="all")

# Plot the specific internal energy with the solution
subfigure(2,2,3)
plot("x","u")
plotanalytical("x","u")
limit("x",-1.0,1.0,window="all")

# Plot the thermal pressure
subfigure(2,2,4)
plot("x","press")
plotanalytical("x","press")
limit("x",-1.0,1.0,window="all")

# 'Sleep hack' (to allow matplotlib to update the figure) 
# before running the simulation
time.sleep(2)
run()
block()
