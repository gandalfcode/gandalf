#==============================================================================
# adsodtest.py
#==============================================================================
from gandalf.analysis.facade import *
import matplotlib.pyplot as plt
from matplotlib import rc
from mpl_toolkits.axes_grid1 import AxesGrid
import time


rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)


# Set all plot limits
xmin   = -9.5
xmax   = 9.5
rhomin = 0.2
rhomax = 1.15
vxmin  = -0.2
vxmax  = 0.95
umin   = 1.7
umax   = 2.8

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
fig, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, figsize=(6,12))

ax1.set_ylabel(r"$\rho$")
ax1.set_xlim([xmin, xmax])
ax1.set_ylim([rhomin, rhomax])
ax1.scatter(x0, rho0, color='black', marker='*', label='Gradh-SPH')
ax1.scatter(x1, rho1, color='blue', marker='+', label='MFV-moving')
ax1.scatter(x2, rho2, color='green', marker='x', label='MFV-static')
ax1.plot(rhodata.x_data, rhodata.y_data, color="red", label='Riemann solution')
ax1.legend()

ax2.set_ylabel(r"$v_x$")
ax2.set_ylim([vxmin, vxmax])
ax2.scatter(x0, vx0, color='black', marker='*', label='Gradh-SPH')
ax2.scatter(x1, vx1, color='blue', marker='+', label='MFV')
ax2.plot(vxdata.x_data, vxdata.y_data, color="red", label='Riemann solution')

ax3.set_xlabel(r"x")
ax3.set_ylabel(r"$u$")
ax3.set_ylim([umin, umax])
ax3.scatter(x0, u0, color='black', marker='*', label='Gradh-SPH')
ax3.scatter(x1, u1, color='blue', marker='+', label='MFV')
ax3.plot(udata.x_data, udata.y_data, color="red", label='Riemann solution')


plt.tight_layout()
plt.show()
fig.savefig('adsod-rho.pdf')


block()
