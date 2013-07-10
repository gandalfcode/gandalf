from seren.analysis.facade import *
from seren.analysis.data_fetcher import *
from seren.analysis.compute import lagrangian_radii
import time

#newsim('freefall.dat')
#setupsim()
#run()

loadsim('FREEFALL1')
CreateTimeData('lr1',lagrangian_radii,0.1)
CreateTimeData('lr2',lagrangian_radii,0.5)
CreateTimeData('lr3',lagrangian_radii,0.9)

plot_vs_time("lr2")
block()
