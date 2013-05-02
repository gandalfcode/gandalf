#==============================================================================
# khitest.py
# Run Kelvin-Helmholtz test using initial conditions file 'khi.dat'.
#==============================================================================
from seren.analysis.facade import *
import time

# Create new simulation object and generate initial conditions
newsim("khi.dat")

# Create rendered density plot
render("x","y","rho",res=128)

# Run the simulation and then block to retain the rendering window
run()
block()
