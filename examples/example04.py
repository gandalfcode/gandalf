#==============================================================================
# example04.py
# Example to create a simple 2D random cube of particles and plot the x-y
# coordinates to screen as the simulation progresses.
#==============================================================================
from gandalf.analysis.facade import *

# Create simulation object from `glass.dat' parameters file
sim = newsim("glass.dat")

# Set the time interval where the matplotlib window is updated in seconds
sim.SetParam("dt_python",2.0)

# Set-up simulation, plot the x-y coordinates and run the simulation
setupsim()
plot("x","y")
run()

# Finally, save the figure to file
savefig("figure.eps")
block()
