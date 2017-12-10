#==============================================================================
# adsodtest.py
#==============================================================================
from gandalf.analysis.facade import *
import matplotlib.pyplot as plt
from matplotlib import rc
from mpl_toolkits.axes_grid1 import AxesGrid
import time


#rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('font', **{'family': 'normal', 'weight' : 'bold', 'size' : 14})
rc('text', usetex=True)

# Set all plot limits
xmin   = -9.5
xmax   = 9.5
rhomin = 0.15
rhomax = 1.115
vxmin  = -0.15
vxmax  = 0.95
umin   = 1.7
umax   = 2.7
Nleft  = 240
Nright = 60
Nunequal = 60


# Extract data from Grad-h SPH simulation
gradhsphsim = newsim("adsod-gradhsph.dat")
gradhsphsim.SetParam('gas_eos',"isothermal")
gradhsphsim.SetParam('gamma_eos',1.11)
gradhsphsim.SetParam('Nlattice1[0]',Nleft)
gradhsphsim.SetParam('Nlattice2[0]',Nright)
setupsim()
run()
x0   = get_data('x') #, sim=0, snap=10)
rho0 = get_data('rho') #, sim=0, snap=10)
vx0  = get_data('vx') #, sim=0, snap=10)

# Extract data from Grad-h SPH simulation
gradhsphsim_unequal = newsim("adsod-gradhsph.dat")
gradhsphsim_unequal.SetParam('gas_eos',"isothermal")
gradhsphsim_unequal.SetParam('gamma_eos',1.11)
gradhsphsim_unequal.SetParam('Nlattice1[0]',Nunequal)
gradhsphsim_unequal.SetParam('Nlattice2[0]',Nunequal)
setupsim()
run()
x1   = get_data('x') #, sim=0, snap=10)
rho1 = get_data('rho') #, sim=0, snap=10)
vx1  = get_data('vx') #, sim=0, snap=10)

# Extract data from MFV simulation
mfvsim = newsim("adsod-mfv-moving.dat")
mfvsim.SetParam('gas_eos',"isothermal")
mfvsim.SetParam('gamma_eos',1.11)
mfvsim.SetParam('Nlattice1[0]',Nleft)
mfvsim.SetParam('Nlattice2[0]',Nright)
setupsim()
run()
x2   = get_data('x') #, sim=1, snap=10)
rho2 = get_data('rho') #, sim=1, snap=10)
vx2  = get_data('vx') #, sim=1, snap=10)

# Extract data from MFV simulation
mfvsim_unequal = newsim("adsod-mfv-moving.dat")
mfvsim_unequal.SetParam('gas_eos',"isothermal")
mfvsim_unequal.SetParam('gamma_eos',1.11)
mfvsim_unequal.SetParam('Nlattice1[0]',Nunequal)
mfvsim_unequal.SetParam('Nlattice2[0]',Nunequal)
setupsim()
run()
x3   = get_data('x') #, sim=1, snap=10)
rho3 = get_data('rho') #, sim=1, snap=10)
vx3  = get_data('vx') #, sim=1, snap=10)

# Extract data from MFV simulation
mfmsim = newsim("adsod-mfm-moving.dat")
mfmsim.SetParam('gas_eos',"isothermal")
mfmsim.SetParam('gamma_eos',1.11)
mfmsim.SetParam("riemann_solver", "exact")
setupsim()
run()
x4   = get_data('x') #, sim=2, snap=10)
rho4 = get_data('rho') #, sim=2, snap=10)
vx4  = get_data('vx') #, sim=2, snap=10)

# Extract data from MFV simulation
mfmsim_unequal = newsim("adsod-mfm-moving.dat")
mfmsim_unequal.SetParam('gas_eos',"isothermal")
mfmsim_unequal.SetParam('gamma_eos',1.11)
mfmsim_unequal.SetParam("riemann_solver", "exact")
mfmsim_unequal.SetParam('Nlattice1[0]',Nunequal)
mfmsim_unequal.SetParam('Nlattice2[0]',Nunequal)
setupsim()
run()
x5   = get_data('x') #, sim=1, snap=10)
rho5 = get_data('rho') #, sim=1, snap=10)
vx5  = get_data('vx') #, sim=1, snap=10)

# Extract data from MFV simulation
mfmsim_hllc = newsim("adsod-mfm-moving.dat")
mfmsim_hllc.SetParam('gas_eos',"isothermal")
mfmsim_hllc.SetParam('gamma_eos',1.11)
mfmsim_hllc.SetParam("riemann_solver", "hllc")
setupsim()
run()
x6   = get_data('x') #, sim=2, snap=10)
rho6 = get_data('rho') #, sim=2, snap=10)
vx6  = get_data('vx') #, sim=2, snap=10)

# Extract data from MFV simulation
mfmsim_unequal_hllc = newsim("adsod-mfm-moving.dat")
mfmsim_unequal_hllc.SetParam('gas_eos',"isothermal")
mfmsim_unequal_hllc.SetParam('gamma_eos',1.11)
mfmsim_unequal_hllc.SetParam("riemann_solver", "hllc")
mfmsim_unequal_hllc.SetParam('Nlattice1[0]',Nunequal)
mfmsim_unequal_hllc.SetParam('Nlattice2[0]',Nunequal)
setupsim()
run()
x7   = get_data('x') #, sim=1, snap=10)
rho7 = get_data('rho') #, sim=1, snap=10)
vx7  = get_data('vx') #, sim=1, snap=10)


# Extract data for analytical solution
rhodata = get_analytical_data("x","rho") #,sim=0,snap=10)
vxdata = get_analytical_data("x","vx") #,sim=0,snap=10)
udata = get_analytical_data("x","u") #,sim=0,snap=10)


# Create matplotlib figure object with shared x-axis
fig, axarr = plt.subplots(1, 3, sharex='col', sharey='row', figsize=(13,3.5))
fig.subplots_adjust(hspace=0.0001, wspace=0.0001)
fig.subplots_adjust(bottom=0.14, top=0.99, left=0.045, right=0.99)

axarr[0].set_ylabel(r"$\rho$")
axarr[0].set_xlabel(r"$x$")
axarr[0].set_xlim([xmin, xmax])
axarr[0].set_ylim([rhomin, rhomax])
axarr[0].plot(rhodata.x_data, rhodata.y_data, linestyle='-', color="red", label='Exact solution')
axarr[0].scatter(x0, rho0, marker='o', facecolors='none', edgecolors='blue', s=10, label='Gradh-SPH, equal-mass')
axarr[0].scatter(x1, rho1, color='black', marker='+', s=32, label='Gradh-SPH, unequal-mass')
axarr[0].legend(fontsize=10)

axarr[1].set_xlabel(r"$x$")
axarr[1].set_xlim([xmin, xmax])
axarr[1].set_ylim([rhomin, rhomax])
axarr[1].plot(rhodata.x_data, rhodata.y_data, linestyle='-', color="red", label='Exact solution')
axarr[1].scatter(x2, rho2, marker='o', facecolors='none', edgecolors='blue', s=10, label='MFV, equal-mass')
axarr[1].scatter(x3, rho3, color='black', marker='+', s=32, label='MFV, unequal-mass')
axarr[1].legend(fontsize=10)

axarr[2].set_xlabel(r"$x$")
axarr[2].set_xlim([xmin, xmax])
axarr[2].set_ylim([rhomin, rhomax])
axarr[2].plot(rhodata.x_data, rhodata.y_data, linestyle='-', color="red", label='Exact solution')
axarr[2].scatter(x4, rho4, marker='o', facecolors='none', edgecolors='blue', s=10, label='MFM, equal-mass')
axarr[2].scatter(x5, rho5, color='black', marker='+', s=32, label='MFM, unequal-mass')
axarr[2].legend(fontsize=10)

#axarr[3].set_xlabel(r"$x$")
#axarr[3].set_xlim([xmin, xmax])
#axarr[3].set_ylim([rhomin, rhomax])
#axarr[3].plot(rhodata.x_data, rhodata.y_data, linestyle='-', color="red", label='Exact solution')
#axarr[3].scatter(x6, rho6, color='black', marker='+', label='MFM, equal-mass')
#axarr[3].scatter(x7, rho7, marker='o', facecolors='none', edgecolors='blue', label='MFM, unequal-mass')
#axarr[3].legend(fontsize=10)


plt.show()
fig.savefig('isosod.pdf', dpi=50)
