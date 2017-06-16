#==============================================================================
# freefalltest.py
# Run the freefall collapse test using initial conditions specified in the
# file 'freefall.dat'.
#==============================================================================
from gandalf.analysis.facade import *
from gandalf.analysis.data_fetcher import *
from gandalf.analysis.compute import particle_data
from gandalf.analysis.SimBuffer import SimBuffer, BufferException
import time
import matplotlib.pyplot as plt
import numpy as np
import math
from matplotlib import rc
from mpl_toolkits.axes_grid1 import AxesGrid



#--------------------------------------------------------------------------------------------------
rc('font', **{'family': 'normal', 'weight' : 'bold', 'size' : 16})
rc('text', usetex=True)


# Binary parameters
m1 = 0.5
m2 = 0.5
abin = 1.0
ebin = 0.5
etot0 = -0.5*m1*m2/abin
period = 2.0*math.pi*math.sqrt(abin*abin*abin/(m1 + m2))

xmin = -0.6
xmax = 2.1
ymin = -0.85
ymax = 0.85
xsize = xmax - xmin
ysize = ymax - ymin



CreateTimeData('x',particle_data,quantity='x')
CreateTimeData('y',particle_data,quantity='y')


# Leapfrog KDK
kdksim = newsim('binaryorbit.dat')
kdksim.SetParam('nbody','lfkdk')
setupsim()
run()
x_kdk = get_time_data("t","x")
y_kdk = get_time_data("t","y")


# Leapfrog DKD
dkdsim = newsim('binaryorbit.dat')
dkdsim.SetParam('nbody','lfdkd')
setupsim()
run()
x_dkd = get_time_data("t","x")
y_dkd = get_time_data("t","y")


# 4th-order Hermite
hermite4sim = newsim('binaryorbit.dat')
hermite4sim.SetParam('nbody','hermite4')
setupsim()
run()
x_hermite4 = get_time_data("t","x")
y_hermite4 = get_time_data("t","y")


# 4th-order Hermite TS
hermite4tssim = newsim('binaryorbit.dat')
hermite4tssim.SetParam('nbody','hermite4ts')
hermite4tssim.SetParam('Npec',5)
setupsim()
run()
x_4ts = get_time_data("t","x")
y_4ts = get_time_data("t","y")


# 6th-order Hermite
#hermite6tssim = newsim('binaryorbit.dat')
#hermite6tssim.SetParam('nbody','hermite6ts')
#hermite6tssim.SetParam('Npec',5)
#setupsim()
#run()
#x_6ts = get_time_data("t","x")
#y_6ts = get_time_data("t","y")



# Create matplotlib figure object with shared x-axis
#--------------------------------------------------------------------------------------------------
#fig, axarr = plt.subplots(2, 1, sharex='col', sharey='row', figsize=(10,4))
fig, axarr = plt.subplots(4, 1, figsize=(6,11), sharex='col', sharey='row')
fig.subplots_adjust(hspace=0.001, wspace=0.001)
fig.subplots_adjust(bottom=0.06, top=0.98, left=0.14, right=0.98)

axarr[0].set_ylabel(r"$y$")
axarr[0].set_ylim([ymin, ymax])
axarr[0].set_xlim([xmin, xmax])
axarr[0].plot(x_kdk.y_data, y_kdk.y_data, color="black", linestyle='-', label='Leapfrog KDK', lw=1.0)
axarr[0].text(xmin + 0.02*xsize, ymax - 0.1*ysize, "(a) Leapfrog-KDK", fontsize=12)

axarr[1].set_ylabel(r"$y$")
axarr[1].set_ylim([ymin, ymax])
axarr[1].plot(x_dkd.y_data, y_dkd.y_data, color="black", linestyle='-', label='Leapfrog DKD', lw=1.0)
axarr[1].text(xmin + 0.02*xsize, ymax - 0.1*ysize, "(b) Leapfrog-DKD", fontsize=12)

axarr[2].set_ylabel(r"$y$")
axarr[2].set_ylim([ymin, ymax])
axarr[2].plot(x_hermite4.y_data, y_hermite4.y_data, color="black", linestyle='-', label='4H', lw=1.0)
axarr[2].text(xmin + 0.02*xsize, ymax - 0.1*ysize, "(c) 4th-order Hermite", fontsize=12)

axarr[3].set_xlabel(r"$x$")
axarr[3].set_ylabel(r"$y$")
axarr[3].set_ylim([ymin, ymax])
axarr[3].plot(x_4ts.y_data, y_4ts.y_data, color="black", linestyle='-', label='4TS', lw=1.0)
axarr[3].text(xmin + 0.02*xsize, ymax - 0.1*ysize, "(d) 4th-order Hermite TS", fontsize=12)


plt.show()
fig.savefig('binaryorbit.pdf', dpi=50)


# Prevent program from closing before showing plot window
block()
