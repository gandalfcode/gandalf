#==============================================================================
# example17.py
# Create a new user-defined quantity from a function and plot on window.
# Extends on example13.py
#==============================================================================
from gandalf.analysis.facade import *

# Load in simulation frmo disk
sim = loadsim("ADSOD1")

# Create new quantity, specific kinetic energy of particles, including the
# scaling factor (specific energy) and latex label.
CreateUserQuantity("ke","0.5*vx*vx",scaling_factor="u",
                   label="$\\frac{1}{2}v^2$")

# Plot defined quantity along with internal energy
plot("x","ke",snap=2)
limit("ke",0.0,2.7)
addplot("x","u",snap=2)
block()

# Define a function for computing kinetic energy
def ComputeKineticEnergy(snap,type="default",unit="default"):
    vx=get_data("vx",snap=snap,type=type,unit=unit)
    return 0.5*vx*vx
    
# Create another new quantity
CreateUserQuantity("ke2",ComputeKineticEnergy,scaling_factor="u",
                   label="$\\frac{1}{2}v^2$")
# Overplot the new quantity - see that it overlays perfectly on the previous one
addplot("x","ke2",snap=2)
block()