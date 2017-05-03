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
def FreefallSolution(rho, mfrac):
    tff = np.sqrt(3.0*3.1415/32.0/rho)
    r0 = math.pow(mfrac, 0.33333333333)
    r = np.arange(0.0, 0.99999, 0.0001)
    t = np.arccos(np.sqrt(r/r0)) + np.sqrt(r/r0)*np.sqrt(1.0 - r/r0)
    t *= 2.0/3.14157
    return t,r



#--------------------------------------------------------------------------------------------------
rc('font', **{'family': 'normal', 'weight' : 'bold', 'size' : 16})
rc('text', usetex=True)

# Set all plot limits
tmin   = 0.00
tmax   = 1.0
rmin   = 0.005
rmax   = 1.05
rho    = 3.0/4.0/3.114157
radius = 1.0
stride = 8
sim_no = 0

# Analytical solutions for grav. acceleration and potential
r_acc = np.arange(rmin, 0.99999*rmax, 0.001)
a_acc = -r_acc
gpot_acc = 0.5*(r_acc*r_acc - 3.0*radius*radius)

# Run the simulation
mainsim = newsim('freefall.dat')
setupsim()
run()

# Get grav. data for simulaton for plotting
r_data    = get_data("r", sim=sim_no, snap=0)
a_data    = get_data("ar", sim=sim_no, snap=0)
gpot_data = get_data("gpot", sim=sim_no, snap=0)

# Prepare 10%, 50% and 90% Lagrangian radii
CreateTimeData('lr1',lagrangian_radii,mfrac=0.05)
CreateTimeData('lr2',lagrangian_radii,mfrac=0.2)
CreateTimeData('lr3',lagrangian_radii,mfrac=0.5)
CreateTimeData('lr4',lagrangian_radii,mfrac=1.0)
data05  = get_time_data("t","lr1")
data20  = get_time_data("t","lr2")
data50  = get_time_data("t","lr3")
data100 = get_time_data("t","lr4")

# Get analytical solutions for each mass fraction
t05, r05   = FreefallSolution(rho, 0.05)
t20, r20   = FreefallSolution(rho, 0.2)
t50, r50   = FreefallSolution(rho, 0.5)
t100, r100 = FreefallSolution(rho, 1.0)

# Normalise freefall data to units of freefall time
data05.x_data  /= 1.1107
data20.x_data  /= 1.1107
data50.x_data  /= 1.1107
data100.x_data /= 1.1107


# Create matplotlib figure object with shared x-axis
#--------------------------------------------------------------------------------------------------
#fig, axarr = plt.subplots(2, 1, sharex='col', sharey='row', figsize=(10,4))
fig, axarr = plt.subplots(2, 1, figsize=(7,8), sharex='col')
fig.subplots_adjust(hspace=0.0001, wspace=0.0001)
fig.subplots_adjust(bottom=0.07, top=0.97, left=0.12, right=0.98)

axarr[0].set_ylabel(r"$a_{_{\rm r}}$")
#axarr[0].set_xlabel(r"$r$")
axarr[0].set_xlim([0.0001, 1.04])
axarr[0].set_ylim(-1.1,0.1)
axarr[0].plot(r_acc, a_acc, color="red", linestyle='-', lw=0.5)
axarr[0].scatter(r_data[::stride], a_data[::stride], color='black', marker='.', s=4.0)
#axarr[0].legend(fontsize=12)
axarr[1].set_ylabel(r"$\phi$")
axarr[1].set_xlabel(r"$r$")
axarr[1].set_xlim([rmin, rmax])
axarr[1].set_ylim(-1.55,-0.85)
axarr[1].plot(r_acc, gpot_acc, color="red", linestyle='-', lw=0.5)
axarr[1].scatter(r_data[::stride], -gpot_data[::stride], color='black', marker='.', s=4.0)
#axarr[1].legend(fontsize=12)

plt.show()
fig.savefig('sphereaccel.pdf', dpi=50)



# Create matplotlib figure object with shared x-axis
#--------------------------------------------------------------------------------------------------
fig2, axarr2 = plt.subplots(1, 1, figsize=(7,5))
fig2.subplots_adjust(hspace=0.001, wspace=0.001)
fig2.subplots_adjust(bottom=0.1, top=0.99, left=0.1, right=0.98)

axarr2.set_ylabel(r"$R_{_{\rm LAG}}$")
axarr2.set_xlabel(r"$t/t_{_{\rm FF}}$")
axarr2.set_xlim([tmin, tmax])
axarr2.set_ylim([rmin, rmax])
axarr2.plot(t05, r05, color="red", linestyle='-', lw=0.5)
axarr2.plot(t20, r20, color="red", linestyle='-', lw=0.5)
axarr2.plot(t50, r50, color="red", linestyle='-', lw=0.5)
axarr2.plot(t100, r100, color="red", linestyle='-', lw=0.5)
axarr2.scatter(data05.x_data, data05.y_data, color='black', marker='.', s=24.0, label='$5\%$')
axarr2.scatter(data20.x_data, data20.y_data, color='black', marker=',', s=24.0, label='$20\%$')
axarr2.scatter(data50.x_data, data50.y_data, color='black', marker='D', s=24.0, label='$50\%$')
axarr2.scatter(data100.x_data, data100.y_data, color='black', marker='^', s=24.0, label='$100\%$')
axarr2.legend(fontsize=12)

plt.show()
fig2.savefig('freefall.pdf', dpi=50)

# Prevent program from closing before showing plot window
block()
