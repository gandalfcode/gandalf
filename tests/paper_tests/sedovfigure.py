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
xmax   = 0.37
rhomin = 0.01
rhomax = 5.9
vmin   = -0.3
vmax   = 2.2
umin   = 0.02
umax   = 2000.0
stride = 4

# Extract data from MFV simulation
mfvsim_l1 = newsim("sedov2d-mfv-moving.dat")
mfvsim_l1.SetParam('Nlevels',1)
mfvsim_l1.SetParam('time_step_limiter','none')
setupsim()
run()
R0   = get_data('R')
rho0 = get_data('rho')
v0   = get_data('vR2d')

# Extract data from MFV simulation
mfvsim_l10 = newsim("sedov2d-mfv-moving.dat")
mfvsim_l10.SetParam('time_step_limiter','none')
mfvsim_l10.SetParam('Nlevels',10)
setupsim()
run()
R1   = get_data('R')
rho1 = get_data('rho')
v1   = get_data('vR2d')

# Extract data from MFV simulation
mfvsim_l10_sm = newsim("sedov2d-mfv-moving.dat")
mfvsim_l10_sm.SetParam('time_step_limiter','simple')
mfvsim_l10_sm.SetParam('Nlevels',10)
mfvsim_l10_sm.SetParam('level_diff_max',1)
setupsim()
run()
R2   = get_data('R')
rho2 = get_data('rho')
v2   = get_data('vR2d')

# Extract data from MFV simulation
gradhsphsim_l10 = newsim("sedov2d-gradh.dat")
gradhsphsim_l10.SetParam('time_step_limiter','none')
gradhsphsim_l10.SetParam('Nlevels',10)
gradhsphsim_l10.SetParam('level_diff_max',20)
setupsim()
run()
R3   = get_data('R')
rho3 = get_data('rho')
v3   = get_data('vR2d')

# Extract data from MFV simulation
gradhsphsim_l10_sm = newsim("sedov2d-gradh.dat")
gradhsphsim_l10_sm.SetParam('time_step_limiter','simple')
gradhsphsim_l10_sm.SetParam('Nlevels',10)
gradhsphsim_l10_sm.SetParam('level_diff_max',1)
setupsim()
run()
R4   = get_data('R')
rho4 = get_data('rho')
v4   = get_data('vR2d')


# Extract data for analytical solution
rhodata = get_analytical_data("R","rho")
vxdata = get_analytical_data("R","vR2d")
udata = get_analytical_data("R","u")

# Create matplotlib figure object with shared x-axis
fig, axarr = plt.subplots(2, 5, sharex='col', sharey='row', figsize=(14,5.5))
fig.subplots_adjust(hspace=0.0001, wspace=0.0001)
fig.subplots_adjust(bottom=0.1, top=0.98, left=0.05, right=0.99)


# MFV with 1 timestep level
axarr[0,0].set_ylabel(r"$\rho$", fontsize=11)
axarr[0,0].set_xlim([xmin, xmax])
axarr[0,0].set_ylim([rhomin, rhomax])
axarr[0,0].scatter(R0[::stride], rho0[::stride], color='black', marker='.', s=1.0, label='MFV, L=1')
axarr[0,0].plot(rhodata.x_data, rhodata.y_data, color="red", label='Solution')
axarr[0,0].legend(fontsize=11)
axarr[1,0].set_ylabel(r"$v_{_R}$", fontsize=11)
axarr[1,0].set_xlabel(r"$x$", fontsize=11)
axarr[1,0].set_ylim([vmin, vmax])
axarr[1,0].scatter(R0[::stride], v0[::stride], color='black', marker='.', s=1.0)
axarr[1,0].plot(vxdata.x_data, vxdata.y_data, color="red")

# MFV with 10 timestep levels
axarr[0,1].set_xlim([xmin, xmax])
axarr[0,1].set_ylim([rhomin, rhomax])
axarr[0,1].scatter(R1[::stride], rho1[::stride], color='black', marker='.', s=1.0, label='MFV, L=10')
axarr[0,1].plot(rhodata.x_data, rhodata.y_data, color="red", label='Solution')
axarr[0,1].legend(fontsize=11)
axarr[1,1].set_xlabel(r"$x$", fontsize=11)
axarr[1,1].scatter(R1[::stride], v1[::stride], color='black', marker='.', s=1.0)
axarr[1,1].plot(vxdata.x_data, vxdata.y_data, color="red")

# MFV with 10 timestep levels + SM limiter
axarr[0,2].set_xlim([xmin, xmax])
axarr[0,2].set_ylim([rhomin, rhomax])
axarr[0,2].scatter(R2[::stride], rho2[::stride], color='black', marker='.', s=1.0, label='MFV, L=10 + SM')
axarr[0,2].plot(rhodata.x_data, rhodata.y_data, color="red", label='Solution')
axarr[0,2].legend(fontsize=11)
axarr[1,2].set_xlabel(r"$x$", fontsize=11)
axarr[1,2].scatter(R2[::stride], v2[::stride], color='black', marker='.', s=1.0)
axarr[1,2].plot(vxdata.x_data, vxdata.y_data, color="red")

# Gradh-SPH with 10 timestep levels
axarr[0,3].set_xlim([xmin, xmax])
axarr[0,3].set_ylim([rhomin, rhomax])
axarr[0,3].scatter(R3[::stride], rho3[::stride], color='black', marker='.', s=1.0, label='Gradh, L=10')
axarr[0,3].plot(rhodata.x_data, rhodata.y_data, color="red", label='Solution')
axarr[0,3].legend(fontsize=11)
axarr[1,3].set_xlabel(r"$x$", fontsize=11)
axarr[1,3].scatter(R3[::stride], v3[::stride], color='black', marker='.', s=1.0)
axarr[1,3].plot(vxdata.x_data, vxdata.y_data, color="red")

# Gradh-SPH with 10 timestep levels + SM limiter
axarr[0,4].set_xlim([xmin, xmax])
axarr[0,4].set_ylim([rhomin, rhomax])
axarr[0,4].scatter(R4[::stride], rho4[::stride], color='black', marker='.', s=1.0, label='Gradh, L=10 + SM')
axarr[0,4].plot(rhodata.x_data, rhodata.y_data, color="red", label='Solution')
axarr[0,4].legend(fontsize=11)
axarr[1,4].set_xlabel(r"$x$", fontsize=11)
axarr[1,4].scatter(R4[::stride], v4[::stride], color='black', marker='.', s=1.0)
axarr[1,4].plot(vxdata.x_data, vxdata.y_data, color="red")


plt.show()
fig.savefig('sedov.pdf', dpi=50)
