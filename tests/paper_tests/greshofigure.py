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
xmin = 0.002
xmax = 0.69
vmin = -0.1
vmax = 1.1


# Extract data from Grad-h SPH simulation
loadsim('GRESHO-GRADHSPH')
x0   = get_data('R', sim=0, snap=1)
v0   = get_data('vphi', sim=0, snap=1)

# Extract data from static MFV simulation
loadsim('GRESHO-MFV-STATIC')
x1   = get_data('R', sim=1, snap=1)
v1   = get_data('vphi', sim=1, snap=1)

# Extract data from moving MFV simulation
loadsim('GRESHO-MFV-MOVING')
x2   = get_data('R', sim=2, snap=1)
v2   = get_data('vphi', sim=2, snap=1)

# Extract data from moving MFM simulation
loadsim('GRESHO-MFM-MOVING')
x3   = get_data('R', sim=3, snap=1)
v3   = get_data('vphi', sim=3, snap=1)


# Extract data for analytical solution
vxdata = get_analytical_data("R","vphi",sim=0,snap=1)
pressdata = get_analytical_data("R","press",sim=0,snap=1)

# Create matplotlib figure object with shared x-axis
fig, axarr = plt.subplots(1, 4, sharex='col', sharey='row', figsize=(14.5,3.5))
fig.subplots_adjust(hspace=0.02, wspace=0.0001)
fig.subplots_adjust(bottom=0.15, top=0.98, left=0.04, right=0.99)


# Grad-h SPH simulation
axarr[0].set_ylabel(r"$v_{\phi}$", fontsize=14)
axarr[0].set_xlabel(r"$R$", fontsize=14)
axarr[0].set_xlim([xmin, xmax])
axarr[0].set_ylim([vmin, vmax])
axarr[0].scatter(x0, v0, color='black', marker='.', s=1.0, label='Gradh-SPH')
axarr[0].plot(vxdata.x_data, vxdata.y_data, color="red", label='Solution')
axarr[0].legend(fontsize=12)


# Static MFV simulation
axarr[1].set_xlabel(r"$R$", fontsize=16)
axarr[1].set_xlim([xmin, xmax])
axarr[1].set_ylim([vmin, vmax])
axarr[1].scatter(x1, v1, color='black', marker='.', s=1.0, label='MFV-static')
axarr[1].plot(vxdata.x_data, vxdata.y_data, color="red", label='Solution')
axarr[1].legend(fontsize=12)


# Moving MFV simulation
axarr[2].set_xlabel(r"$R$", fontsize=16)
axarr[2].set_xlim([xmin, xmax])
axarr[2].set_ylim([vmin, vmax])
axarr[2].scatter(x2, v2, color='black', marker='.', s=1.0, label='MFV-moving')
axarr[2].plot(vxdata.x_data, vxdata.y_data, color="red", label='Solution')
axarr[2].legend(fontsize=12)


# Moving MFM simulation
axarr[3].set_xlabel(r"$R$", fontsize=16)
axarr[3].set_xlim([xmin, xmax])
axarr[3].set_ylim([vmin, vmax])
axarr[3].scatter(x3, v3, color='black', marker='.', s=1.0, label='MFM-moving')
axarr[3].plot(vxdata.x_data, vxdata.y_data, color="red", label='Solution')
axarr[3].legend(fontsize=12)


#fig.tight_layout()
plt.show()
fig.savefig('gresho.pdf', dpi=50)
