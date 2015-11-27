#==============================================================================
# bondiaccretiontest.py
# Run the Binary accretion test using initial conditions specified in the 
# file 'binaryaccretion.dat'.
#==============================================================================
from gandalf.analysis.facade import *
from gandalf.analysis.data_fetcher import *
import time

# Create and set-up simulation object from 'bossbodenheimer.dat' file
loadsim("BONDI1-R4.0")

# Plot the particle positions
subfigure(2,2,1)
plot("r","rho")

# Create rendered slice of density
subfigure(2,2,2)
plot("r","vr")

# Plot star/sink masses
subfigure(2,2,3)
time_plot("t","m",type='star',id=0)

# Plot sound and star speed
subfigure(2,2,4)
#plot("r","temp")
time_plot("t","r",type='star',id=0)

run()
block()
