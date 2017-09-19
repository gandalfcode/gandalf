#==============================================================================
# freefalltest.py
# Run the freefall collapse test using initial conditions specified in the
# file 'freefall.dat'.
#==============================================================================
from gandalf.analysis.facade import *
from gandalf.analysis.data_fetcher import *
from gandalf.analysis.compute import lagrangian_radii
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

# Set all plot limits
tmin     = 0.00
tmax     = 98.0
rmin     = 0.13
rmax     = 11.0
tsize    = tmax - tmin
rsize    = rmax - rmin
rplummer = 1.0
mplummer = 1.0
tcross   = math.sqrt(6.0*rplummer*rplummer*rplummer/mplummer)
stride   = 8
sim_no   = 0

# Run the simulation
mainsim = newsim('plummer.dat')
setupsim()
run()

# Prepare 10%, 50% and 90% Lagrangian radii
CreateTimeData('lr1',lagrangian_radii,mfrac=0.1)
CreateTimeData('lr2',lagrangian_radii,mfrac=0.5)
CreateTimeData('lr3',lagrangian_radii,mfrac=0.9)
data10 = get_time_data("t","lr1",type='star')
data50 = get_time_data("t","lr2",type='star')
data90 = get_time_data("t","lr3",type='star')

# Normalise freefall data to units of the crossing time
data10.x_data /= tcross
data50.x_data /= tcross
data90.x_data /= tcross
tmax /= tcross
tsize /= tcross


# Run the simulation
hybridsim = newsim('hybridplummer.dat')
setupsim()
run()

CreateTimeData('gas1',lagrangian_radii,mfrac=0.1)
CreateTimeData('gas2',lagrangian_radii,mfrac=0.5)
CreateTimeData('gas3',lagrangian_radii,mfrac=0.9)
CreateTimeData('star1',lagrangian_radii,mfrac=0.1)
CreateTimeData('star2',lagrangian_radii,mfrac=0.5)
CreateTimeData('star3',lagrangian_radii,mfrac=0.9)
gas_data10 = get_time_data("t","gas1",type='sph')
gas_data50 = get_time_data("t","gas2",type='sph')
gas_data90 = get_time_data("t","gas3",type='sph')
star_data10 = get_time_data("t","star1",type='star')
star_data50 = get_time_data("t","star2",type='star')
star_data90 = get_time_data("t","star3",type='star')

# Normalise freefall data to units of the crossing time
gas_data10.x_data /= tcross
gas_data50.x_data /= tcross
gas_data90.x_data /= tcross
star_data10.x_data /= tcross
star_data50.x_data /= tcross
star_data90.x_data /= tcross


# Create matplotlib figure object with shared x-axis
#--------------------------------------------------------------------------------------------------
fig, axarr = plt.subplots(1, 1, figsize=(8,6), sharex='col', sharey='row')
fig.subplots_adjust(hspace=0.001, wspace=0.001)
fig.subplots_adjust(bottom=0.09, top=0.98, left=0.1, right=0.98)

#axarr[0].set_ylabel(r"$R_{_{\rm LAG}}$")
##axarr[0].set_xlabel(r"$t/t_{_{\rm CROSS}}$")
#axarr[0].set_yscale("log")
#axarr[0].set_xlim([tmin, tmax])
#axarr[0].set_ylim([rmin, rmax])
#axarr[0].plot(data10.x_data, data10.y_data, linestyle='-', color='red', lw=1.0, label='Stars')
#axarr[0].plot(data50.x_data, data50.y_data, linestyle='-', color='red', lw=1.0)
#axarr[0].plot(data90.x_data, data90.y_data, linestyle='-', color='red', lw=1.0)
#axarr[0].text(tmin + 0.03*tsize, 0.7, "$10\%$", fontsize=12)
#axarr[0].text(tmin + 0.03*tsize, 2.0, "$50\%$", fontsize=12)
#axarr[0].text(tmin + 0.03*tsize, 5.0, "$90\%$", fontsize=12)
#axarr[0].text(tmin + 0.03*tsize, rmax - 0.3*rsize, "(a) $N_s = 500$", fontsize=12)
#legend = axarr[0].legend(loc='upper right', fontsize=12)

axarr.set_ylabel(r"$R_{_{\rm LAG}}$")
axarr.set_xlabel(r"$t/t_{_{\rm CROSS}}$")
axarr.set_yscale("log")
axarr.set_xlim([tmin, tmax])
axarr.set_ylim([rmin, rmax])
axarr.plot(star_data10.x_data, star_data10.y_data, linestyle='-', color='red', lw=1.0, label='Stars')
axarr.plot(star_data50.x_data, star_data50.y_data, linestyle='-', color='red', lw=1.0)
axarr.plot(star_data90.x_data, star_data90.y_data, linestyle='-', color='red', lw=1.0)
axarr.plot(gas_data10.x_data, gas_data10.y_data, linestyle='--', color='black', lw=1.0, label='Gas')
axarr.plot(gas_data50.x_data, gas_data50.y_data, linestyle='--', color='black', lw=1.0)
axarr.plot(gas_data90.x_data, gas_data90.y_data, linestyle='--', color='black', lw=1.0)
axarr.text(tmin + 0.02*tsize, 0.7, "$10\%$", fontsize=16)
axarr.text(tmin + 0.02*tsize, 1.5, "$50\%$", fontsize=16)
axarr.text(tmin + 0.02*tsize, 4.2, "$90\%$", fontsize=16)
#axarr.text(tmin + 0.02*tsize, rmax - 0.2*rsize, "$N_s = 500$, $N_g = 5000$", fontsize=16)
legend = axarr.legend(loc='upper right', fontsize=16)

plt.show()
fig.savefig('hybridplummer.pdf', dpi=50)

# Prevent program from closing before showing plot window
block()
