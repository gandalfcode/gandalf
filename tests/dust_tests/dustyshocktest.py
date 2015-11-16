#==============================================================================
# adsodtest.py
# Run the adiabatic Sod shocktube test using initial conditions specified in
# the file 'adsod.dat' and then plotting important quantities together 
# with the analytical solutions.
#==============================================================================
from gandalf.analysis.facade import *
import time

# Create new shocktube simulation from 'adsod.dat' file
newsim('dustyshock.dat')
setupsim()

# Plot the density with the analytical solution
subfigure(2,2,1)
plot("x","rho", type='sph')
plot("x","rho", type='dust', overplot=True)
limit("x",-0.5,0.5,window="all")
limit("rho",0.2,1.1,window='all')

# Plot the x-velocity with the analytical solution
subfigure(2,2,2)
plot("x","vx", type='sph')
plot("x","vx", type='dust', overplot=True)
limit("x",-0.5,0.5,window="all")
limit("vx",-0.05,0.8,window='all')

# Plot the specific internal energy with the solution
subfigure(2,2,3)
plot("x","u", type='sph')
limit("x",-0.5,0.5,window="all")
limit("u",1.9,3.2,window='all')

# Plot the thermal pressure
subfigure(2,2,4)
plot("x","press", type='sph')
limit("x",-0.5,0.5,window="all")
limit("press",0.1,1.1,window='all')



# 'Sleep hack' (to allow matplotlib to update the figure) 
# before running the simulation
time.sleep(2)
run()
block()
