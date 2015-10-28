#==============================================================================
# example12.py
# Plot particle quantities in an alternative coordinate system.
#==============================================================================
from gandalf_dust.analysis.facade import *
from matplotlib.colors import LogNorm

# Create simulation object from Boss-Bodenheimer parameters file
sim = newsim("bossbodenheimer.dat")
sim.SetParam("tend",0.02)
setupsim()

# Run simulation and plot x-y positions of SPH particles in the default
# units specified in the `bossbodenheimer.dat' parameters file.

plot("x","y")
addplot("x","y",type="star")
limit("x",-0.007,0.007)
limit("y",-0.007,0.007)

window()
render("x","y","rho",res=256,#norm=LogNorm(),
       interpolation='bicubic')
limit("x",-0.007,0.007)
limit("y",-0.007,0.007)

run()
block()

# After pressing return, re-plot last snapshot but in new specified units (au).
plot("x","y",xunit="au",yunit="au")
window()
render("x","y","rho",res=256,#norm=LogNorm(),
       interpolation='bicubic')
limit("x",-0.007,0.007)
limit("y",-0.007,0.007)
block()
