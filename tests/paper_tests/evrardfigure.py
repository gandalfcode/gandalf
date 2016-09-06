#==============================================================================
# evrardfigure.py
# ...
#==============================================================================
from gandalf.analysis.facade import *
import matplotlib.pyplot as plt
import numpy as np
import math
from matplotlib import rc
from mpl_toolkits.axes_grid1 import AxesGrid
import time

#rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('font', **{'family': 'normal', 'weight' : 'bold', 'size' : 16})
rc('text', usetex=True)


# Set all plot limits
xmin   = 0.008
xmax   = 0.8
rhomin = 0.03333
rhomax = 666.666
vxmin  = -1.8
vxmax  = 0.2
umin = 0.03333
umax = 3.3333
pmin = 0.03333
pmax = 333.333
entropymin = 0.001
entropymax = 0.22
pi = 3.1415
richardSolution = True
stride = 8


if richardSolution == True:

    solution_file = 'evrard-solution-t0.8.dat'
    # Read in the header:
    with open(solution_file, 'r') as f:
        header = f.readline()[1].split()

    # Get the rest of the data:
    r, rho, u, v = np.genfromtxt(solution_file).T
    press = (2/3.) * rho*u
    entropy = press / rho**(5/3.)

else:

    # Scaling factors from 'profile.idl' file
    rost = 0.75 #/pi
    est = 1.054811e-1
    pst = rost*est
    vst = math.sqrt(est)
    rst = 1.257607


    # Read in solution file (Non-pythonic way!!)
    with open('ppm1oaa', 'r') as f:
        firstline = f.readline()
        header = firstline.split()
        Ntot = 350
        t = float(header[1])*math.sqrt(6.67e-8) #*8.0/pi/pi

        r = np.zeros(Ntot)
        v = np.zeros(Ntot)
        rho = np.zeros(Ntot)
        press = np.zeros(Ntot)
        entropy = np.zeros(Ntot)

        print 'N : ',Ntot,'    time : ',t

        read_data = f.readlines()
        i = 0

        for line in read_data:
            data = line.split()
            r[i]     = float(data[1])/rst*1e-11
            rho[i]   = float(data[2])/rost
            v[i]     = float(data[4])/vst*1e-8
            press[i] = float(data[3])/pst*1e-16
            entropy[i] = press[i]/math.pow(rho[i], 1.666666);
            entropy[i] = 0.66666666*0.05*math.pow(2.0*pi*r[i], 0.6666666)
            i += 1


# Extract data from Grad-h SPH simulation
snapno = 4
loadsim('EVRARD-GRADH')
x0   = get_data('r', sim=0, snap=snapno)
rho0 = get_data('rho', sim=0, snap=snapno)
v0   = get_data('vr', sim=0, snap=snapno)
u0   = get_data('u', sim=0, snap=snapno)
p0   = get_data('press', sim=0, snap=snapno)
#entropy0 = get_data('entropy', sim=0, snap=snapno)
entropy0 = p0 / rho0**(5/3.)

loadsim('EVRARD-MFM-MOVING')
x1   = get_data('r', sim=1, snap=snapno)
rho1 = get_data('rho', sim=1, snap=snapno)
v1   = get_data('vr', sim=1, snap=snapno)
u1   = get_data('u', sim=1, snap=snapno)
p1   = get_data('press', sim=1, snap=snapno)
#entropy1 = get_data('entropy', sim=1, snap=snapno)
entropy1 = p1 / rho1**(5/3.)

# Create matplotlib figure object with shared x-axis
fig, axarr = plt.subplots(3, 2, sharex='col', sharey='row', figsize=(8,8))
fig.subplots_adjust(hspace=0.001, wspace=0.001)
fig.subplots_adjust(bottom=0.08, top=0.99, left=0.12, right=0.97)


print 'r   : ',np.min(r),np.max(r)
print 'rho : ',np.min(rho),np.max(rho)
print 'v   : ',np.min(v),np.max(v)

axarr[0,0].set_xscale("log")
axarr[0,0].set_yscale("log")
axarr[0,0].set_ylabel(r"$\rho$")
axarr[0,0].set_xlim([xmin, xmax])
axarr[0,0].set_ylim([rhomin, rhomax])
axarr[0,0].scatter(x0[::stride], rho0[::stride], color='black', marker='.', s=1.0, label='Gradh-SPH')
axarr[0,0].plot(r, rho, color="red", label='PPM solution')
axarr[0,0].legend(fontsize=12)

axarr[0,1].set_xscale("log")
axarr[0,1].set_yscale("log")
axarr[0,1].set_xlim([xmin, xmax])
axarr[0,1].set_ylim([rhomin, rhomax])
axarr[0,1].scatter(x1[::stride], rho1[::stride], color='black', marker='.', s=1.0, label='MFM moving')
axarr[0,1].plot(r, rho, color="red", label='PPM solution')
axarr[0,1].legend(fontsize=12)

axarr[1,0].set_xscale("log")
axarr[1,0].set_ylabel(r"$v$")
axarr[1,0].set_xlim([xmin, xmax])
axarr[1,0].set_ylim([vxmin, vxmax])
axarr[1,0].scatter(x0[::stride], v0[::stride], color='black', marker='.', s=1.0)
axarr[1,0].plot(r, v, color="red", label='PPM solution')

axarr[1,1].set_xscale("log")
axarr[1,1].set_xlim([xmin, xmax])
axarr[1,1].set_ylim([vxmin, vxmax])
axarr[1,1].scatter(x1[::stride], v1[::stride], color='black', marker='.', s=1.0, label='MFM moving')
axarr[1,1].plot(r, v, color="red", label='PPM solution')

axarr[2,0].set_xscale("log")
axarr[2,0].set_ylabel(r"$v$")
axarr[2,0].set_xlim([xmin, xmax])
axarr[2,0].set_ylim([entropymin, entropymax])
axarr[2,0].scatter(x0[::stride], entropy0[::stride], color='black', marker='.', s=1.0)
axarr[2,0].plot(r, entropy, color="red", label='PPM solution')
axarr[1,0].legend(fontsize=12)

axarr[2,1].set_xscale("log")
axarr[2,1].set_xlim([xmin, xmax])
axarr[2,1].set_ylim([entropymin, entropymax])
axarr[2,1].scatter(x1[::stride], entropy1[::stride], color='black', marker='.', s=1.0, label='MFM moving')
axarr[2,1].plot(r, entropy, color="red", label='PPM solution')


#fig.tight_layout()
plt.show()
fig.savefig('evrard.pdf', dpi=50)
