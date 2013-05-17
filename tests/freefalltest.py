from seren.analysis.facade import *
from seren.analysis.compute import lagrangian_radii
import time

newsim('freefall.dat')
setupsim()
plot("r","rho")
print 'rlag : ',lagrangian_radii(0.5)
run()
print 'rlag : ',lagrangian_radii(0.5)
plot("r","rho")
block()
