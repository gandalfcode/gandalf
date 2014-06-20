#==============================================================================
# khitest.py
# Run Kelvin-Helmholtz test using initial conditions file 'khi.dat'.
#==============================================================================
from gandalf.analysis.facade import *
import time

# Create new simulation object and generate initial conditions
newsim("khi.dat")
setupsim()

# Create rendered density plot
render("x","y","rho",res=128)
limit("rho",1.0,2.0)
limit("x",-0.5,0.5)
limit("y",-0.5,0.5)

# Run the simulation and then block to retain the rendering window
run()
block()
