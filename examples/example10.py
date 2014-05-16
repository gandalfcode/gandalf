#==============================================================================
# example10.py
# Create a rendered image during an interactive simulation.
#==============================================================================
from gandalf.analysis.facade import *

# Create simulation object from Kelvin-Helmholtz parameters file
sim = newsim("khi.dat")
setupsim()

# Generate rendered image of density field on the x-y plane with a
# resolution of 128 x 128 pixels.
render("x","y","rho",res=128)
limit("x",-0.5,0.5)
limit("y",-0.5,0.5)
limit("rho",1.0,2.0)
run()
block()
