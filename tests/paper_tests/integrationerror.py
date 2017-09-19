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

errormin = 6.666666e-12
errormax = 6.666666e-1


# Energy error
CreateTimeData('energy_error',energy_error,etot0=etot0)

# Leapfrog KDK
kdksim = newsim('binary.dat')
kdksim.SetParam('nbody','lfkdk')
setupsim()
run()
kdkerror = get_time_data("t","energy_error",linestyle='-')

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
hermite4error.x_data /= period
hermite4tserror.x_data /= period



rplummer = 1.0
mplummer = 1.0
tmin     = 0.0
tmax     = 98.0
tsize    = tmax - tmin
tcross   = math.sqrt(6.0*rplummer*rplummer*rplummer/mplummer)
tmax     /= tcross
errormin = 6.66666e-11
errormax = 3.33333e-1

# Leapfrog KDK
kdksim2 = newsim('plummer.dat')
kdksim2.SetParam('nbody','lfkdk')
setupsim()
etot2 = kdksim2.GetInitialEnergy()
CreateTimeData('energy_error2',energy_error,etot0=etot2)
run()

# Energy error
kdkerror2 = get_time_data("t","energy_error2",linestyle='-')
kdkerror2.x_data /= tcross

# 4th-order Hermite
hermite4sim2 = newsim('plummer.dat')
hermite4sim2.SetParam('nbody','hermite4')
setupsim()
run()
hermite4error2 = get_time_data("t","energy_error2",linestyle='-')
hermite4error2.x_data /= tcross

# 4th-order Hermite
hermite4tssim2 = newsim('plummer.dat')
hermite4tssim2.SetParam('nbody','hermite4ts')
hermite4tssim2.SetParam('Npec',3)
setupsim()
run()
hermite4tserror2 = get_time_data("t","energy_error2",linestyle='-')
hermite4tserror2.x_data /= tcross


# Leapfrog KDK
kdksim3 = newsim('plummer.dat')
kdksim3.SetParam('nbody','lfkdk')
kdksim3.SetParam('Nlevels',5)
setupsim()
etot3 = kdksim3.GetInitialEnergy()
CreateTimeData('energy_error3',energy_error,etot0=etot3)
run()
kdkerror3 = get_time_data("t","energy_error2",linestyle='-')
kdkerror3.x_data /= tcross

# 4th-order Hermite
hermite4sim3 = newsim('plummer.dat')
hermite4sim3.SetParam('nbody','hermite4')
hermite4sim3.SetParam('Nlevels',5)
setupsim()
run()
hermite4error3 = get_time_data("t","energy_error3",linestyle='-')
hermite4error3.x_data /= tcross

# 4th-order Hermite
hermite4tssim3 = newsim('plummer.dat')
hermite4tssim3.SetParam('nbody','hermite4ts')
hermite4tssim3.SetParam('Npec',3)
hermite4tssim3.SetParam('Nlevels',5)
setupsim()
run()
hermite4tserror3 = get_time_data("t","energy_error3",linestyle='-')
hermite4tserror3.x_data /= tcross



# Create matplotlib figure object with shared x-axis
#--------------------------------------------------------------------------------------------------
fig, axarr = plt.subplots(3, 1, sharex='col', sharey='row', figsize=(8,18))
#fig, axarr = plt.subplots(2, 1, figsize=(14,8), sharex='row')
fig.subplots_adjust(hspace=0.0001, wspace=0.0001)
fig.subplots_adjust(bottom=0.06, top=0.99, left=0.11, right=0.98)


axarr[0].set_ylabel(r"$\Delta E/E_0$")
#axarr[0].set_xlabel(r"$t$")
axarr[0].set_yscale("log")
axarr[0].set_ylim([errormin, errormax])
axarr[0].set_xlim([0.0, tmax])
axarr[0].plot(kdkerror.x_data, kdkerror.y_data, color="red", linestyle=':', label='Leapfrog KDK', lw=1.0)
axarr[0].plot(hermite4error.x_data, hermite4error.y_data, color="black", linestyle='-', label='4th order Hermite', lw=1.0)
axarr[0].plot(hermite4tserror.x_data, hermite4tserror.y_data, color="blue", linestyle='--', label='4th order Hermite TS', lw=1.0)
axarr[0].text(tmin + 0.007*tsize, 0.085*errormax, "(a) Binary orbit, $e=0.1$", fontsize=16)
legend = axarr[0].legend(loc='upper right', fontsize=16)


axarr[1].set_ylabel(r"$\Delta E/E_0$")
axarr[1].set_xlabel(r"$t/t_{_{\rm CROSS}}$")
axarr[1].set_yscale("log")
axarr[1].set_ylim([errormin, errormax])
axarr[1].set_xlim([0.0, tmax])
axarr[1].plot(kdkerror2.x_data, kdkerror2.y_data, color="red", linestyle=':', label='Leapfrog KDK', lw=1.0)
axarr[1].plot(hermite4error2.x_data, hermite4error2.y_data, color="black", linestyle='-', label='4th order Hermite', lw=1.0)
axarr[1].plot(hermite4tserror2.x_data, hermite4tserror2.y_data, color="blue", linestyle='--', label='4th order Hermite TS', lw=1.0)
axarr[1].text(tmin + 0.007*tsize, 0.085*errormax, "(b) Plummer sphere, $N = 200$, $L = 1$", fontsize=16)


axarr[2].set_ylabel(r"$\Delta E/E_0$")
axarr[2].set_xlabel(r"$t/t_{_{\rm CROSS}}$")
axarr[2].set_yscale("log")
axarr[2].set_ylim([errormin, errormax])
axarr[2].set_xlim([0.0, tmax])
axarr[2].plot(kdkerror3.x_data, kdkerror3.y_data, color="red", linestyle=':', label='Leapfrog KDK', lw=1.0)
axarr[2].plot(hermite4error3.x_data, hermite4error3.y_data, color="black", linestyle='-', label='4th order Hermite', lw=1.0)
axarr[2].plot(hermite4tserror3.x_data, hermite4tserror3.y_data, color="blue", linestyle='--', label='4th order Hermite TS', lw=1.0)
axarr[2].text(tmin + 0.007*tsize, 0.085*errormax, "(c) Plummer sphere, $N = 200$, $L = 5$", fontsize=16)


plt.show()
fig.savefig('integrationerror.pdf', dpi=100)


# Prevent program from closing before showing plot window
block()
