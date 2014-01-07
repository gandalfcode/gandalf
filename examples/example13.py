#==============================================================================
# example13.py
# Create a new user-defined quantity and plot on window.
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
