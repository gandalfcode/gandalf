#==============================================================================
# freefalltest.py
# Run the freefall collapse test using initial conditions specified in the
# file 'freefall.dat'.
#==============================================================================
from gandalf.analysis.facade import *
from gandalf.analysis.data_fetcher import *
from gandalf.analysis.compute import energy_error
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
ttot = 150.0
ttot /= period

errormin = 6.666666e-12
errormax = 3.3333e-1


# Energy error
CreateTimeData('energy_error',energy_error,etot0=etot0)

# Leapfrog KDK
kdksim = newsim('binary.dat')
kdksim.SetParam('nbody','lfkdk')
setupsim()
run()
kdkerror = get_time_data("t","energy_error",linestyle='-')

# Leapfrog DKD
dkdsim = newsim('binary.dat')
dkdsim.SetParam('nbody','lfdkd')
setupsim()
run()
dkderror = get_time_data("t","energy_error",linestyle='-')


# 4th-order Hermite
hermite4sim = newsim('binary.dat')
hermite4sim.SetParam('nbody','hermite4')
setupsim()
run()
hermite4error = get_time_data("t","energy_error",linestyle='-')


# 4th-order Hermite
hermite4tssim = newsim('binary.dat')
hermite4tssim.SetParam('nbody','hermite4ts')
hermite4tssim.SetParam('Npec',2)
setupsim()
run()
hermite4tserror = get_time_data("t","energy_error",linestyle='-')


# Normalise times
kdkerror.x_data /= period
dkderror.x_data /= period
hermite4error.x_data /= period
hermite4tserror.x_data /= period


# 6th-order Hermite
#hermite6tssim = newsim('binary.dat')
#hermite6tssim.SetParam('nbody','hermite6ts')
#hermite6tssim.SetParam('Npec',5)
#setupsim()
#run()
#hermite6tserror = get_time_data("t","energy_error",linestyle='-.')



# Create matplotlib figure object with shared x-axis
#--------------------------------------------------------------------------------------------------
#fig, axarr = plt.subplots(2, 1, sharex='col', sharey='row', figsize=(10,4))
fig, axarr = plt.subplots(1, 1, figsize=(14,4), sharex='row')
#fig.subplots_adjust(hspace=0.001, wspace=0.001)
fig.subplots_adjust(bottom=0.11, top=0.98, left=0.07, right=0.98)

axarr.set_ylabel(r"$\Delta E/E_0$")
axarr.set_xlabel(r"$t$")
axarr.set_yscale("log")
axarr.set_ylim([errormin, errormax])
axarr.set_xlim([0.0, ttot])
axarr.plot(kdkerror.x_data, kdkerror.y_data, color="red", linestyle=':', label='Leapfrog KDK', lw=1.0)
#axarr.plot(dkderror.x_data, dkderror.y_data, color="blue", linestyle='-', label='Leapfrog DKD', lw=1.0)
axarr.plot(hermite4error.x_data, hermite4error.y_data, color="black", linestyle='-', label='4th order Hermite', lw=1.0)
axarr.plot(hermite4tserror.x_data, hermite4tserror.y_data, color="blue", linestyle='--', label='4th order Hermite TS', lw=1.0)
#axarr.plot(hermite6tserror.x_data, hermite6tserror.y_data, color="blue", linestyle='--', label='6th order Hermite TS', lw=1.0)
#axarr.legend(fontsize=12)
legend = axarr.legend(loc='upper right', fontsize=12)


plt.show()
fig.savefig('binaryerror.pdf', dpi=50)


# Prevent program from closing before showing plot window
block()
