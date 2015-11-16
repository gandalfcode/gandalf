#==============================================================================
# example18.py
# Plot quantities as a function of time through all snapshots in the simulation.
# Extends on example14.py
#==============================================================================
from gandalf.analysis.facade import *
#In this case, we need also a function defined in compute
from gandalf.analysis.compute import lagrangian_radii

# Load in simulation from disk (from example 12)
sim = loadsim("BB1")

# Define the half-mass radius
CreateTimeData("half_r",lagrangian_radii,mfrac=0.5)

# Plot it
time_plot("t","half_r")
block()

# Define a function for computing the total mass
def ComputeMass(snap,type="default",unit="default"):
    m=get_data('m',snap=snap,type=type,unit=unit)
    return m.sum()

# Define the quantity
CreateTimeData("totm",ComputeMass)

# Plot it
time_plot("t","totm")
block()