from seren.analysis.facade import *
import time

newsim('soundwave.dat')
plot("x","rho")
plotanalytical("x","rho")
time.sleep(1)
run()
block()
