#==============================================================================
# jeanstest.py
# Run Noh test using initial conditions file 'noh.dat'.
#==============================================================================
from gandalf.analysis.facade import *
import time

# Create new Noh test simulation object from 'noh.dat' file
newsim('jeans.dat')
setupsim()

# Plot the density with the analytical solution
subfigure(3,1,1)
plot("x","rho")
plotanalytical("x","rho")
limit("x",0.0,1.0,window="all")

# Plot the x-velocity with the analytical solution
subfigure(3,1,2)
plot("x","vx")
plotanalytical("x","vx")
limit("x",0.0,1.0,window="all")

# Plot the acceleration with the solution
subfigure(3,1,3)
plot("x","ax")
plotanalytical("x","ax")
limit("x",0.0,1.0,window="all")

run()
block()
