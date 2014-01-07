#==============================================================================
# example07.py
# Overplot the analytical solution of a known problem while simulatneously
# running the simulation.
#==============================================================================
from gandalf.analysis.facade import *

# Create the new simulation object
sim = newsim("adsod.dat")
sim.SetParam("dt_python",2.0)
setupsim()

# Plot both the simulation results along with the analytical solution
# (N.B. the analytical solution is overplotted by default).
plot("x","rho")
plotanalytical("x","rho")

# Set the limits of the x-axis plot
limit("x",-1.1,1.1)

# Finally run the simulation to completion
run()
block()
