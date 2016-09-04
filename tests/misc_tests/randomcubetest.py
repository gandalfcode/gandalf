#==============================================================================
# randomcubetest.py
# Relax 2D cube of randomlly placed particles to a glass-like structure.
#==============================================================================
from gandalf.analysis.facade import *
import time

# Create new random cube simulation object from 'randomcube.dat' file
newsim('randomcube.dat')
setupsim()

# Plot various quantities of cube
subfigure(2,2,1)
plot("x","y")
subfigure(2,2,2)
plot("x","rho")
subfigure(2,2,3)
plot("x","vy")
subfigure(2,2,4)
plot("h","ax")

run()
block()


