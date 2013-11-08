#==============================================================================
# binaryaccretiontest.py
# Run the Binary accretion test using initial conditions specified in the 
# file 'binaryaccretion.dat'.
#==============================================================================
from gandalf.analysis.facade import *
from gandalf.analysis.data_fetcher import *
import time

# Create and set-up simulation object from 'bossbodenheimer.dat' file
newsim('binaryaccretion.dat')
setupsim()

# Plot the particle positions
subfigure(2,2,1)
plot("x","y")
addplot("x","y",type="star")
#limit("x",-0.5,0.5)
#limit("y",-0.5,0.5)

# Create rendered slice of density
subfigure(2,2,2)
render("x","y","rho",res=128)
#limit("x",-0.5,0.5)
#limit("y",-0.5,0.5)

# Plot star/sink masses
subfigure(2,2,3)
plot("x","vy")
#plot_vs_time_quantity_particle("m",type='star',id=0)

# Plot sound and star speed
subfigure(2,2,4)
plot("x","vx")

run()
block()
