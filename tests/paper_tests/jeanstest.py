#==============================================================================
# jeanstest.py
# Run Noh test using initial conditions file 'noh.dat'.
#==============================================================================
from gandalf.analysis.facade import *
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
xmin   = 0.0
xmax   = 1.0
rhomin = 0.95
rhomax = 1.05
stride = 7
sim_no = 0

# Create new Noh test simulation object from 'noh.dat' file
sim1 = newsim('jeans.dat')
#sim1.SetParam('Nstepsmax', 1)
setupsim()
run()

# Get grav. data for simulaton for plotting
x_data   = get_data("x", sim=sim_no)
vx_data  = get_data("vx", sim=sim_no)
rho_data = get_data("rho", sim=sim_no)

rho_analytical = get_analytical_data("x", "rho", sim=sim_no)
vx_analytical  = get_analytical_data("x", "vx", sim=sim_no)


# Create matplotlib figure object with shared x-axis
#--------------------------------------------------------------------------------------------------
fig, axarr = plt.subplots(2, 1, figsize=(7,10), sharex='row')
fig.subplots_adjust(hspace=0.001, wspace=0.001)
fig.subplots_adjust(bottom=0.08, top=0.97, left=0.13, right=0.98)

axarr[0].set_ylabel(r"$\rho$")
axarr[0].set_xlim([xmin, xmax])
axarr[0].plot(rho_analytical.x_data, rho_analytical.y_data, color="red", linestyle='-', lw=0.5)
axarr[0].scatter(x_data[::stride], rho_data[::stride], color='black', marker='.', s=4.0)
axarr[1].set_ylabel(r"$v_{x}$")
axarr[1].set_xlabel(r"$x$")
axarr[1].set_xlim([xmin, xmax])
axarr[1].plot(vx_analytical.x_data, vx_analytical.y_data, color="red", linestyle='-', lw=0.5)
axarr[1].scatter(x_data[::stride], vx_data[::stride], color='black', marker='.', s=4.0)

plt.show()
fig.savefig('jeanstest.eps', dpi=50)


run()
block()
