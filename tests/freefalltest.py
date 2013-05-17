from seren.analysis.facade import *
import time

newsim('freefall.dat')
setupsim()
plot("r","rho")
run()
plot("r","rho")
block()


