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

# Plot results to screen

# Plot the density with the analytical solution
subfigure(2,2,1)
plot("R","rho")
plotanalytical("R","rho")
limit("R",0,1,window="all")


# Plot the velocity with the analytical solution
subfigure(2,2,2)
#plot("R","vr")
#plotanalytical("R","vr")
limit("R",0,1,window="all")

# Plot the specific internal energy with the solution
subfigure(2,2,3)
plot("R","u")
plotanalytical("R","u")
limit("R",0,1,window="all")
limit("u",0,10,window="all")

# Plot the thermal pressure
subfigure(2,2,4)
plot("R","press")
plotanalytical("R","press")
limit("R",0,1,window="all")

# Run simulation
run()


block()


