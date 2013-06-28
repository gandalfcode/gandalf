#==============================================================================
# plummertest.py
# Run the Plummer sphere test using initial conditions specified in the file 
# 'plummer.dat' and then plot various quantities.
#==============================================================================
from seren.analysis.facade import *
from seren.analysis.data_fetcher import *
from seren.analysis.compute import lagrangian_radii
import time

# Create new shocktube simulation from 'adshock.dat' file
newsim('hybridplummer.dat')
setupsim()

# Plot the density with the analytical solution
subfigure(2,2,1)
plot("x","y")
addplot("x","y",type='star')
limit("x",-5.5,5.5,window="all")
limit("y",-5.5,5.5,window="all")

subfigure(2,2,2)
plot("r","rho")

subfigure(2,2,3)
plot("r","ar")

subfigure(2,2,4)
plot("r","sound")
addplot("r","vr",type='star')

# 'Sleep hack' (to allow matplotlib to update the figure) 
# before running the simulation
time.sleep(1)
run()
block()
