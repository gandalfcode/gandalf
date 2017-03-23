#==============================================================================
# nbodyjeanstest.py
# Run N-body version of Jeans instability test (i.e. with no pressure support).
#==============================================================================
from gandalf.analysis.facade import *
import time

# Create new Noh test simulation object from 'noh.dat' file
newsim('nbodyjeans.dat')
setupsim()

# Plot the x-velocity with the analytical solution
subfigure(2,1,1)
plot("x","vx",type="star")
plotanalytical("x","vx")
limit("x",0.0,1.0,window="all")

# Plot the acceleration with the solution
subfigure(2,1,2)
plot("x","ax",type="star")
plotanalytical("x","ax")
limit("x",0.0,1.0,window="all")

run()
block()
