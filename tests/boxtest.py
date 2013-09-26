from gandalf.analysis.facade import *
import time

newsim('box.dat')
setupsim()


#plot("h","rho")
plot("x","y")
run()
#limit("rho",0.0,1.5)
block()
