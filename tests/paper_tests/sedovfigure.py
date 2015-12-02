#==============================================================================
# adsodtest.py
#==============================================================================
from gandalf.analysis.facade import *
import matplotlib.pyplot as plt
from matplotlib import rc
from mpl_toolkits.axes_grid1 import AxesGrid
import time


#rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('font', **{'family': 'normal', 'weight' : 'bold', 'size' : 12})
rc('text', usetex=True)

# Set all plot limits
xmin   = 0.002
xmax   = 0.43
rhomin = 0.01
rhomax = 5.4
vmin  = 0.001
vmax  = 2.2
umin   = 0.02
umax   = 2000.0

# Extract data from Grad-h SPH simulation
loadsim('SEDOV2D-GRADHSPH-L5')
x0   = get_data('R', sim=0, snap=6)
rho0 = get_data('rho', sim=0, snap=6)
vx0  = get_data('vR2d', sim=0, snap=6)
u0   = get_data('u', sim=0, snap=6)

# Extract data from MFV simulation
loadsim('SEDOV2D-MFV-MOVING-L5')
x1   = get_data('R', sim=1, snap=6)
rho1 = get_data('rho', sim=1, snap=6)
vx1  = get_data('vR2d', sim=1, snap=6)
u1   = get_data('u', sim=1, snap=6)

# Extract data from MFV simulation
loadsim('SEDOV2D-MFV-STATIC-L5')
x2   = get_data('R', sim=2, snap=6)
rho2 = get_data('rho', sim=2, snap=6)
vx2  = get_data('vR2d', sim=2, snap=6)
u2   = get_data('u', sim=2, snap=6)

# Extract data for analytical solution
rhodata = get_analytical_data("R","rho",sim=0,snap=6)
vxdata = get_analytical_data("R","vR2d",sim=0,snap=6)
udata = get_analytical_data("R","u",sim=0,snap=6)

# Create matplotlib figure object with shared x-axis
fig, axarr = plt.subplots(2, 3, sharex='col', sharey='row', figsize=(11,7))
fig.subplots_adjust(hspace=0.001, wspace=0.001)
fig.subplots_adjust(bottom=0.08, top=0.98, left=0.07, right=0.99)


axarr[0,0].set_ylabel(r"$\rho$", fontsize=20)
axarr[0,0].set_xlim([xmin, xmax])
axarr[0,0].set_ylim([rhomin, rhomax])
axarr[0,0].scatter(x0, rho0, color='black', marker='.', s=1.0, label='Gradh-SPH')
axarr[0,0].plot(rhodata.x_data, rhodata.y_data, color="red", label='Solution')
axarr[0,0].legend(fontsize=12)

axarr[1,0].set_ylabel(r"$v_{_R}$", fontsize=20)
axarr[1,0].set_xlabel(r"$x$", fontsize=20)
axarr[1,0].set_ylim([vmin, vmax])
#axarr[1,0].set_ylim([vxmin, vxmax])
axarr[1,0].scatter(x0, vx0, color='black', marker='.', s=1.0, label='Gradh-SPH')
axarr[1,0].plot(vxdata.x_data, vxdata.y_data, color="red", label='Solution')

#axarr[2,0].set_xlabel(r"$x$")
#axarr[2,0].set_ylabel(r"$u$")
#axarr[2,0].set_ylim([umin, umax])
#axarr[2,0].set_yscale("log")
##axarr[2,0].set_ylim([umin, umax])
#axarr[2,0].scatter(x0, u0, color='black', marker='+', label='Gradh-SPH')
#axarr[2,0].plot(udata.x_data, udata.y_data, color="red", label='Riemann solution')


axarr[0,1].set_xlim([xmin, xmax])
axarr[0,1].set_ylim([rhomin, rhomax])
axarr[0,1].scatter(x1, rho1, color='black', marker='.', s=1.0, label='MFV-moving')
axarr[0,1].plot(rhodata.x_data, rhodata.y_data, color="red", label='Solution')
axarr[0,1].legend(fontsize=12)

#axarr[1,1].set_ylim([vxmin, vxmax])
axarr[1,1].set_xlabel(r"$x$", fontsize=20)
axarr[1,1].scatter(x1, vx1, color='black', marker='.', s=1.0, label='MFV')
axarr[1,1].plot(vxdata.x_data, vxdata.y_data, color="red", label='Solution')

#axarr[2,1].set_xlabel(r"$x$")
##axarr[2,1].set_yscale("log")
##axarr[2,1].set_ylim([umin, umax])
##axarr[2,1].scatter(x1, u1, color='black', marker='+', label='MFV')
#axarr[2,1].plot(udata.x_data, udata.y_data, color="red", label='Riemann solution')


axarr[0,2].set_xlim([xmin, xmax])
axarr[0,2].set_ylim([rhomin, rhomax])
axarr[0,2].scatter(x2, rho2, color='black', marker='.', s=1.0, label='MFV-static')
axarr[0,2].plot(rhodata.x_data, rhodata.y_data, color="red", label='Solution')
axarr[0,2].legend(fontsize=12)

#axarr[1,2].set_ylim([vxmin, vxmax])
axarr[1,2].set_xlabel(r"$x$", fontsize=20)
axarr[1,2].scatter(x2, vx2, color='black', marker='.', s=1.0, label='MFV')
axarr[1,2].plot(vxdata.x_data, vxdata.y_data, color="red", label='Solution')

#axarr[2,2].set_xlabel(r"$x$")
#axarr[2,2].set_yscale("log")
##axarr[2,2].set_ylim([umin, umax])
#axarr[2,2].scatter(x2, u2, color='black', marker='+', label='MFV')
#axarr[2,2].plot(udata.x_data, udata.y_data, color="red", label='Riemann solution')

#fig.tight_layout()
plt.show()
fig.savefig('sedov.pdf')
