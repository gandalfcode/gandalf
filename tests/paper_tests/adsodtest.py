#==============================================================================
# adsodtest.py
#==============================================================================
from gandalf.analysis.facade import *
import matplotlib.pyplot as plt
from matplotlib import rc
from mpl_toolkits.axes_grid1 import AxesGrid
import time


#rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('font', **{'family': 'normal', 'weight' : 'bold', 'size' : 16})
rc('text', usetex=True)

# Set all plot limits
xmin   = -9.5
xmax   = 9.5
rhomin = 0.15
rhomax = 1.15
vxmin  = -0.15
vxmax  = 0.95
umin   = 1.7
umax   = 2.7

# Extract data from Grad-h SPH simulation
loadsim('ADSOD-GRADHSPH')
x0   = get_data('x', sim=0, snap=10)
rho0 = get_data('rho', sim=0, snap=10)
vx0  = get_data('vx', sim=0, snap=10)
u0   = get_data('u', sim=0, snap=10)

# Extract data from MFV simulation
loadsim('ADSOD-MFV-MOVING')
x1   = get_data('x', sim=1, snap=10)
rho1 = get_data('rho', sim=1, snap=10)
vx1  = get_data('vx', sim=1, snap=10)
u1   = get_data('u', sim=1, snap=10)

# Extract data from MFV simulation
loadsim('ADSOD-MFV-STATIC')
x2   = get_data('x', sim=2, snap=10)
rho2 = get_data('rho', sim=2, snap=10)
vx2  = get_data('vx', sim=2, snap=10)
u2   = get_data('u', sim=2, snap=10)

# Extract data for analytical solution
rhodata = get_analytical_data("x","rho",sim=0,snap=10)
vxdata = get_analytical_data("x","vx",sim=0,snap=10)
udata = get_analytical_data("x","u",sim=0,snap=10)

# Create matplotlib figure object with shared x-axis
fig, axarr = plt.subplots(3, 3, sharex='col', sharey='row', figsize=(10,7))
fig.subplots_adjust(hspace=0.001, wspace=0.001)
fig.subplots_adjust(bottom=0.08, top=0.99, left=0.06, right=0.99)

axarr[0,0].set_ylabel(r"$\rho$")
axarr[0,0].set_xlim([xmin, xmax])
axarr[0,0].set_ylim([rhomin, rhomax])
axarr[0,0].scatter(x0, rho0, color='black', marker='+', label='Gradh-SPH')
#axarr[0,0].scatter(x1, rho1, color='blue', marker='+', label='MFV-moving')
#axarr[0,0].scatter(x2, rho2, color='green', marker='x', label='MFV-static')
axarr[0,0].plot(rhodata.x_data, rhodata.y_data, color="red", label='Riemann solution')
axarr[0,0].legend(fontsize=12)

axarr[1,0].set_ylabel(r"$v_x$")
axarr[1,0].set_ylim([vxmin, vxmax])
axarr[1,0].scatter(x0, vx0, color='black', marker='+', label='Gradh-SPH')
#axarr[1,0].scatter(x1, vx1, color='blue', marker='+', label='MFV')
axarr[1,0].plot(vxdata.x_data, vxdata.y_data, color="red", label='Riemann solution')

axarr[2,0].set_xlabel(r"$x$")
axarr[2,0].set_ylabel(r"$u$")
axarr[2,0].set_ylim([umin, umax])
axarr[2,0].scatter(x0, u0, color='black', marker='+', label='Gradh-SPH')
#axarr[2,0].scatter(x1, u1, color='blue', marker='+', label='MFV')
axarr[2,0].plot(udata.x_data, udata.y_data, color="red", label='Riemann solution')

axarr[0,1].set_xlim([xmin, xmax])
axarr[0,1].set_ylim([rhomin, rhomax])
#axarr[0,1].scatter(x0, rho0, color='black', marker='*', label='Gradh-SPH')
axarr[0,1].scatter(x1, rho1, color='black', marker='+', label='MFV-moving')
#axarr[0,1].scatter(x2, rho2, color='green', marker='x', label='MFV-static')
axarr[0,1].plot(rhodata.x_data, rhodata.y_data, color="red", label='Riemann solution')
axarr[0,1].legend(fontsize=12)

axarr[1,1].set_ylim([vxmin, vxmax])
#axarr[1,1].scatter(x0, vx0, color='black', marker='*', label='Gradh-SPH')
axarr[1,1].scatter(x1, vx1, color='black', marker='+', label='MFV')
axarr[1,1].plot(vxdata.x_data, vxdata.y_data, color="red", label='Riemann solution')

axarr[2,1].set_xlabel(r"$x$")
axarr[2,1].set_ylim([umin, umax])
#axarr[2,1].scatter(x0, u0, color='black', marker='*', label='Gradh-SPH')
axarr[2,1].scatter(x1, u1, color='black', marker='+', label='MFV')
axarr[2,1].plot(udata.x_data, udata.y_data, color="red", label='Riemann solution')

axarr[0,2].set_xlim([xmin, xmax])
axarr[0,2].set_ylim([rhomin, rhomax])
#axarr[0,1].scatter(x0, rho0, color='black', marker='*', label='Gradh-SPH')
axarr[0,2].scatter(x2, rho2, color='black', marker='+', label='MFV-static')
#axarr[0,1].scatter(x2, rho2, color='green', marker='x', label='MFV-static')
axarr[0,2].plot(rhodata.x_data, rhodata.y_data, color="red", label='Riemann solution')
axarr[0,2].legend(fontsize=12)

axarr[1,2].set_ylim([vxmin, vxmax])
#axarr[1,1].scatter(x0, vx0, color='black', marker='*', label='Gradh-SPH')
axarr[1,2].scatter(x2, vx2, color='black', marker='+', label='MFV')
axarr[1,2].plot(vxdata.x_data, vxdata.y_data, color="red", label='Riemann solution')

axarr[2,2].set_xlabel(r"$x$")
axarr[2,2].set_ylim([umin, umax])
#axarr[2,1].scatter(x0, u0, color='black', marker='*', label='Gradh-SPH')
axarr[2,2].scatter(x2, u2, color='black', marker='+', label='MFV')
axarr[2,2].plot(udata.x_data, udata.y_data, color="red", label='Riemann solution')

#fig.tight_layout()
plt.show()
fig.savefig('adsod.pdf')
