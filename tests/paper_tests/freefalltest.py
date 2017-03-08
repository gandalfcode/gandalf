#==============================================================================
# freefalltest.py
# Run the freefall collapse test using initial conditions specified in the
# file 'freefall.dat'.
#==============================================================================
from gandalf.analysis.facade import *
from gandalf.analysis.data_fetcher import *
from gandalf.analysis.compute import lagrangian_radii
import time
import matplotlib.pyplot as plt
import numpy as np
import math
from matplotlib import rc
from mpl_toolkits.axes_grid1 import AxesGrid
import time


def FreefallSolution(rho, mfrac):
    iMAX = 1000
    tff = np.sqrt(3.0*3.1415/32.0/rho)
    r0 = math.pow(mfrac, 0.33333333333)

    r = np.range(0.0, r0, 0.001*r0)
    t = 2.0*tff*(math.acos(np.sqrt(r/r0)) + np.sqrt(r/r0)*np.sqrt(1.0 - r/r0))/math.pi
    t /= tff
    return t,r

    #r = np.zeros(iMAX)
    #t = np.zeros(iMAX)

    #i = 0
    #while i < iMAX:
    #    r[i] = r0*i/(iMAX - 1)
    #    t[i] = 2.0*tff*(math.acos(np.sqrt(r[i]/r0)) + np.sqrt(r[i]/r0)*np.sqrt(1.0 - r[i]/r0))/3.14157
    #    i = i + 1

    #t /= tff

    #return t,r



#rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('font', **{'family': 'normal', 'weight' : 'bold', 'size' : 16})
rc('text', usetex=True)

# Set all plot limits
tmin = 0.001
tmax = 1.0
rmin = 0.0
rmax = 1.1
rho  = 3.0/4.0/3.114157

CreateTimeData('lr1',lagrangian_radii,mfrac=0.1)
CreateTimeData('lr2',lagrangian_radii,mfrac=0.5)
CreateTimeData('lr3',lagrangian_radii,mfrac=0.9)


# Run the simulation
newsim('freefall.dat')
setupsim()
run()


# Load in simulation and create Lagrangian radii data for 10%, 50% and 90% mass radii.
r_data = get_data("r", sim=0, snap=0)
a_data = get_data("ar", sim=0, snap=0)
data10 = get_time_data("t","lr1",linestyle='-')
data50 = get_time_data("t","lr2",linestyle='-')
data90 = get_time_data("t","lr3",linestyle='-')

# Get analytical solutions for each mass fraction
t10, r10 = FreefallSolution(rho, 0.1)
t50, r50 = FreefallSolution(rho, 0.5)
t90, r90 = FreefallSolution(rho, 0.9)

# Normalise data to freefalltime
data10.x_data /= 1.1107
data50.x_data /= 1.1107
data90.x_data /= 1.1107

# Create matplotlib figure object with shared x-axis
#fig, axarr = plt.subplots(2, 1, sharex='col', sharey='row', figsize=(10,4))
fig, axarr = plt.subplots(1, 2, figsize=(10,4))
fig.subplots_adjust(hspace=0.001, wspace=0.001)
fig.subplots_adjust(bottom=0.11, top=0.97, left=0.11, right=0.98)


#axarr[0,0].set_xscale("log")
#axarr[0,0].set_yscale("log")
axarr[0].set_ylabel(r"$a_{_{\rm r}}$")
axarr[0].set_xlabel(r"$r$")
#axarr[0].set_xlim([tmin, tmax])
#axarr[0].set_ylim([rmin, rmax])
#axarr[0].plot(t10, r10, color="red", linestyle='-', lw=0.5)
axarr[0].scatter(r_data, a_data, color='black', marker=',', s=24.0, label='10\%')
axarr[0].legend(fontsize=12)

axarr[1].set_ylabel(r"$R_{_{\rm LAG}}$")
axarr[1].set_xlabel(r"$t/t_{_{\rm FF}}$")
axarr[1].set_xlim([tmin, tmax])
axarr[1].set_ylim([rmin, rmax])
axarr[1].plot(t10, r10, color="red", linestyle='-', lw=0.5)
axarr[1].plot(t50, r50, color="red", linestyle='-', lw=0.5)
axarr[1].plot(t90, r90, color="red", linestyle='-', lw=0.5)
axarr[1].scatter(data10.x_data, data10.y_data, color='black', marker=',', s=24.0, label='10\%')
axarr[1].scatter(data50.x_data, data50.y_data, color='black', marker='D', s=24.0, label='50\%')
axarr[1].scatter(data90.x_data, data90.y_data, color='black', marker='^', s=24.0, label='90\%')
axarr[1].legend(fontsize=12)

plt.show()
fig.savefig('freefall.pdf', dpi=50)

# Prevent program from closing before showing plot window
block()
