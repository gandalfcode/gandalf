from seren.analysis.facade import *
import time

newsim('sedov.dat')
setupsim()

plot("R","rho")
#plot("x","y")
run()
#plot("x","y")
plot("R","rho")
block()


