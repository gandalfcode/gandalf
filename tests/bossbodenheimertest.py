from seren.analysis.facade import *
from seren.analysis.data_fetcher import *
import time

newsim('bossbodenheimer.dat')
setupsim()
render("x","y","rho")
#plot("x","y") #,xunit="pc",yunit="pc")
run()
render("x","y","rho")
#plot("x","y") #,xunit="pc",yunit="pc")
block()
