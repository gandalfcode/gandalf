#==============================================================================
# adsodtest.py
#==============================================================================
from gandalf.analysis.facade import *
import matplotlib.pyplot as plt
from matplotlib import rc
from mpl_toolkits.axes_grid1 import AxesGrid
import time


#rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('font', **{'family': 'normal', 'weight' : 'bold', 'size' : 11})
rc('text', usetex=True)

# Set all plot limits
xmin   = 0.002
xmax   = 0.39
rhomin = 0.01
rhomax = 5.4
vmin  = 0.001
vmax  = 2.2
umin   = 0.02
umax   = 2000.0
stride = 8

# Extract data from Grad-h SPH simulation
loadsim('SEDOV2D-GRADHSPH')
x0   = get_data('R', sim=0, snap=1)
rho0 = get_data('rho', sim=0, snap=1)
v0   = get_data('vR2d', sim=0, snap=1)

# Extract data from static MFV simulation
loadsim('SEDOV2D-MFV-STATIC')
x1   = get_data('R', sim=1, snap=1)
rho1 = get_data('rho', sim=1, snap=1)
v1   = get_data('vR2d', sim=1, snap=1)

# Extract data from moving MFV simulation
loadsim('SEDOV2D-MFV-MOVING')
x2   = get_data('R', sim=2, snap=1)
rho2 = get_data('rho', sim=2, snap=1)
v2   = get_data('vR2d', sim=2, snap=1)

# Extract data from moving MFM simulation
loadsim('SEDOV2D-MFM-MOVING')
x3   = get_data('R', sim=3, snap=1)
rho3 = get_data('rho', sim=3, snap=1)
v3   = get_data('vR2d', sim=3, snap=1)


# Extract data for analytical solution
rhodata = get_analytical_data("R","rho",sim=0,snap=1)
vxdata = get_analytical_data("R","vR2d",sim=0,snap=1)
udata = get_analytical_data("R","u",sim=0,snap=1)

# Create matplotlib figure object with shared x-axis
fig, axarr = plt.subplots(2, 4, sharex='col', sharey='row', figsize=(12,5.5))
fig.subplots_adjust(hspace=0.02, wspace=0.0001)
fig.subplots_adjust(bottom=0.1, top=0.98, left=0.06, right=0.99)


# Grad-h SPH simulation
axarr[0,0].set_ylabel(r"$\rho$", fontsize=20)
axarr[0,0].set_xlim([xmin, xmax])
axarr[0,0].set_ylim([rhomin, rhomax])
axarr[0,0].scatter(x0[::stride], rho0[::stride], color='black', marker='.', s=1.0, label='Gradh-SPH')
axarr[0,0].plot(rhodata.x_data, rhodata.y_data, color="red", label='Solution')
axarr[0,0].legend(fontsize=12)

axarr[1,0].set_ylabel(r"$v_{_R}$", fontsize=20)
axarr[1,0].set_xlabel(r"$x$", fontsize=20)
axarr[1,0].set_ylim([vmin, vmax])
axarr[1,0].scatter(x0[::stride], v0[::stride], color='black', marker='.', s=1.0)
axarr[1,0].plot(vxdata.x_data, vxdata.y_data, color="red")


# Static MFV simulation
axarr[0,1].set_xlim([xmin, xmax])
axarr[0,1].set_ylim([rhomin, rhomax])
axarr[0,1].scatter(x1[::stride], rho1[::stride], color='black', marker='.', s=1.0, label='MFV-static')
axarr[0,1].plot(rhodata.x_data, rhodata.y_data, color="red", label='Solution')
axarr[0,1].legend(fontsize=12)

axarr[1,1].set_xlabel(r"$x$", fontsize=20)
axarr[1,1].scatter(x1[::stride], v1[::stride], color='black', marker='.', s=1.0)
axarr[1,1].plot(vxdata.x_data, vxdata.y_data, color="red")


# Moving MFV simulation
axarr[0,2].set_xlim([xmin, xmax])
axarr[0,2].set_ylim([rhomin, rhomax])
axarr[0,2].scatter(x2[::stride], rho2[::stride], color='black', marker='.', s=1.0, label='MFV-moving')
axarr[0,2].plot(rhodata.x_data, rhodata.y_data, color="red", label='Solution')
axarr[0,2].legend(fontsize=12)

axarr[1,2].set_xlabel(r"$x$", fontsize=20)
axarr[1,2].scatter(x2[::stride], v2[::stride], color='black', marker='.', s=1.0)
axarr[1,2].plot(vxdata.x_data, vxdata.y_data, color="red")


# Moving MFM simulation
axarr[0,3].set_xlim([xmin, xmax])
axarr[0,3].set_ylim([rhomin, rhomax])
axarr[0,3].scatter(x3[::stride], rho3[::stride], color='black', marker='.', s=1.0, label='MFM-moving')
axarr[0,3].plot(rhodata.x_data, rhodata.y_data, color="red", label='Solution')
axarr[0,3].legend(fontsize=12)

axarr[1,3].set_xlabel(r"$x$", fontsize=20)
axarr[1,3].scatter(x3[::stride], v3[::stride], color='black', marker='.', s=1.0)
axarr[1,3].plot(vxdata.x_data, vxdata.y_data, color="red")


#fig.tight_layout()
plt.show()
fig.savefig('sedov2d.pdf', dpi=50)
