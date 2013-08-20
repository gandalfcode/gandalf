#==============================================================================
# bossbodenheimertest.py
# Run the Boss-Bodenheimer test using initial conditions specified in the 
# file 'bossbodenheimer.dat'.
#==============================================================================
from gandalf.analysis.facade import *
from gandalf.analysis.data_fetcher import *
import time

# Create and set-up simulation object from 'bossbodenheimer.dat' file
newsim('bossbodenheimer.dat')
setupsim()

# Plot the particle positions
subfigure(2,2,1)
plot("x","y")
addplot("x","y",type="star")
limit("x",-0.01,0.01)
limit("y",-0.01,0.01)

# Create rendered slice of density
subfigure(2,2,2)
render("x","y","rho",res=128,zslice=0.0)
limit("x",-0.01,0.01)
limit("y",-0.01,0.01)
limit("rho",1.0e-17,1.0e-13)

# Plot EOS properties
subfigure(2,2,3)
plot("rho","temp",xaxis="log",yaxis="log")

# Plot sink properties
subfigure(2,2,4)
plot("h","rho",xaxis="log",yaxis="log")

run()
block()
