from seren.analysis.facade import *
import time

newsim('noh.dat')
setupsim()

plot("x","rho")
plotanalytical("x","rho")
#plot("x","y")
run()
#plot("x","y")
plot("x","rho")
plotanalytical("x","rho")
block()


