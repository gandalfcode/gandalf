from seren.analysis.facade import *
from seren.analysis.data_fetcher import *
import time

newsim('bossbodenheimer.dat')
setupsim()
#render("x","y","rho",res=128,zslice=0.0)
#render("x","y","rho",res=32) #,xunit="pc",yunit="pc")
plot("x","m")
run()
#render("x","y","rho",res=128,zslice=0.0)
#plot("x","y") #,xunit="pc",yunit="pc")
block()
