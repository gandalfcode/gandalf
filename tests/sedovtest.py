from seren.analysis.facade import *
import time

newsim('sedov.dat')
setupsim()

plot("rcyl","rho")
#plot("x","y")
run()
#plot("x","y")
plot("rcyl","rho")
block()


