from facade import *
import time

newsim('adshock.dat')
subfigure(2,2,1)
plot("x","rho")
subfigure(2,2,2)
plot("x","vx")
subfigure(2,2,3)
plot("x","u")
subfigure(2,2,4)
plot("x","ax")
run()
block()


