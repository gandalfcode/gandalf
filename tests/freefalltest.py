#==============================================================================
# freefalltest.py
# Run the freefall collapse test using initial conditions specified in the
# file 'freefall.dat'.
#==============================================================================
from gandalf.analysis.facade import *
from gandalf.analysis.data_fetcher import *
from gandalf.analysis.compute import lagrangian_radii
import time

# Create new freefall collapse simulation from 'freefall.dat' file and then
# run the simulation
newsim('freefall.dat')
setupsim()
run()

# Create Lagrangian radii data for 10%, 50% and 90% mass radii.
CreateTimeData('lr1',lagrangian_radii,0.1)
CreateTimeData('lr2',lagrangian_radii,0.5)
CreateTimeData('lr3',lagrangian_radii,0.9)

# Plot Lagrangian radii as a function of time
time_plot("t","lr3",linestyle='-')
limit("lr3",0.0,1.05)
time_plot("t","lr2",overplot=True,linestyle='-')
time_plot("t","lr1",overplot=True,linestyle='-')

# Prevent program from closing before showing plot window
block()
