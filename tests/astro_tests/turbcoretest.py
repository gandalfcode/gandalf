#==============================================================================
# turbcoretest.py
# Run the Boss-Bodenheimer test using initial conditions specified in the 
# file 'bossbodenheimer.dat'.
#==============================================================================
from gandalf.analysis.facade import *
from gandalf.analysis.data_fetcher import *
from gandalf.analysis.statistics import structure_function
import time

# Create and set-up simulation object from 'bossbodenheimer.dat' file
sim1 = loadsim('TURBCORE2')
#newsim('turbcore.dat')
#setupsim()


# Plot the particle positions
snap(5)
bins,vmean = structure_function(snap='current',nbin=12,rmin=0.001,rmax=1.0,npoints=8000)

#run()

block()
