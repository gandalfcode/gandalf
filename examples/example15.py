#==============================================================================
# example15.py
# Generate a movie from the snapshots of the simulation
#==============================================================================
from gandalf.analysis.facade import *

#Load in simulation from the disk (run in example01.py)
sim = loadsim('ADSOD1')

#Make a plot
plot('x','rho')

#Loop over all the snapshots and produce a movie
make_movie('MOVIE1.mp4')
block()
