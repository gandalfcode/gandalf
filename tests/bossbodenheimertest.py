#==============================================================================
# bossbodenheimertest.py
# Run the Boss-Bodenheimer test using initial conditions specified in the 
# file 'bossbodenheimer.dat'.
#==============================================================================
from gandalf_dust.analysis.facade import *
from gandalf_dust.analysis.data_fetcher import *
import time
from matplotlib.colors import LogNorm

# Create and set-up simulation object from 'bossbodenheimer.dat' file
newsim('bossbodenheimer.dat')
setupsim()

# Plot the particle positions
subfigure(2,2,1)
plot("x","y")
addplot("x","y",type="star")
limit("x",-0.007,0.007)
limit("y",-0.007,0.007)

# Create rendered slice of density
subfigure(2,2,2)
render("x","y","rho",res=256,norm=LogNorm(),interpolation='bicubic')
limit("x",-0.007,0.007)
limit("y",-0.007,0.007)
#limit("rho",1.0e-16,1.0e-14)

# Plot EOS properties
subfigure(2,2,3)
plot("rho","temp",xaxis="log",yaxis="log")

# Plot sink properties
subfigure(2,2,4)
plot("x","vx")
limit("vx",-1.3,1.3)
addplot("x","vx",type="star")

run()
block()
