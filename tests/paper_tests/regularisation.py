#==============================================================================
# adsodtest.py
#==============================================================================
from gandalf.analysis.facade import *
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc
from mpl_toolkits.axes_grid1 import AxesGrid
import time

#rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('font', **{'family': 'normal', 'weight' : 'bold', 'size' : 14})
rc('text', usetex=True)

# Set all plot limits
xmin    = -0.5
xmax    = 0.5
ymin    = -0.5
ymax    = 0.5
rhomin  = 0.25
rhomax  = 2.25
xsize   = xmax - xmin
rhosize = rhomax - rhomin
Npart   = 512
Nregmax = 100


x = np.arange(xmin, xmax, (xmax - xmin)/1000)
rho = 1.0 + 0.5*np.sin(2.0*np.pi*x)


# MC sampling with no regularisation
noregsim = newsim("basicic.dat")
noregsim.SetParam('regularise_particle_ics',0)
setupsim()
run()
x0   = get_data('x')
y0   = get_data('y')
rho0 = get_data('rho')

# MC sampling with no regularisation
antregsim = newsim("basicic.dat")
antregsim.SetParam('regularise_particle_ics',1)
antregsim.SetParam('alpha_reg',0.0)
antregsim.SetParam('rho_reg',0.9)
antregsim.SetParam('Nreg',Nregmax)
antregsim.SetParam('Nhydro',Npart)
setupsim()
run()
x1   = get_data('x')
y1   = get_data('y')
rho1 = get_data('rho')

# MC sampling with no regularisation
gn2011regsim = newsim("basicic.dat")
gn2011regsim.SetParam('regularise_particle_ics',1)
gn2011regsim.SetParam('alpha_reg',0.1)
gn2011regsim.SetParam('rho_reg',0.0)
gn2011regsim.SetParam('Nreg',Nregmax)
gn2011regsim.SetParam('Nhydro',Npart)
setupsim()
run()
x2   = get_data('x')
y2   = get_data('y')
rho2 = get_data('rho')

# MC sampling with no regularisation
regsim = newsim("basicic.dat")
regsim.SetParam('regularise_particle_ics',1)
regsim.SetParam('alpha_reg',0.1)
regsim.SetParam('rho_reg',0.9)
regsim.SetParam('Nreg',Nregmax)
regsim.SetParam('Nhydro',Npart)
setupsim()
run()
x3   = get_data('x')
y3   = get_data('y')
rho3 = get_data('rho')




# Create matplotlib figure object with shared x-axis
fig, axarr = plt.subplots(2, 4, sharex='col', sharey='row', figsize=(13,6.5))
fig.subplots_adjust(hspace=0.0001, wspace=0.0001)
fig.subplots_adjust(bottom=0.072, top=0.99, left=0.05, right=0.99)

axarr[0,0].set_ylabel(r"$\rho$")
axarr[0,0].set_xlim([xmin, xmax])
axarr[0,0].set_ylim([rhomin, rhomax])
axarr[0,0].plot(x, rho, linestyle='-', color="red")
axarr[0,0].scatter(x0, rho0, marker='.', s=16)
axarr[0,0].text(xmin + 0.03*xsize, rhomax - 0.1*rhosize, '(a) Unregularised', color='black', size=14)
axarr[1,0].set_ylabel(r"$y$")
axarr[1,0].set_xlabel(r"$x$")
axarr[1,0].set_ylim([ymin, ymax])
axarr[1,0].scatter(x0, y0, color='black', marker='.', s=16)

axarr[0,1].set_xlim([xmin, xmax])
axarr[0,1].set_ylim([rhomin, rhomax])
axarr[0,1].plot(x, rho, linestyle='-', color="red")
axarr[0,1].scatter(x2, rho2, marker='.', s=16)
axarr[0,1].text(xmin + 0.03*xsize, rhomax - 0.1*rhosize, '(b) GN2011', color='black', size=14)
axarr[1,1].set_xlabel(r"$x$")
axarr[1,1].set_ylim([ymin, ymax])
axarr[1,1].scatter(x2, y2, color='black', marker='.', s=16)

axarr[0,2].set_xlim([xmin, xmax])
axarr[0,2].set_ylim([rhomin, rhomax])
axarr[0,2].plot(x, rho, linestyle='-', color="red")
axarr[0,2].scatter(x1, rho1, marker='.', s=16)
axarr[0,2].text(xmin + 0.03*xsize, rhomax - 0.1*rhosize, '(c) W95', color='black', size=14)
axarr[1,2].set_xlabel(r"$x$")
axarr[1,2].set_ylim([ymin, ymax])
axarr[1,2].scatter(x1, y1, color='black', marker='.', s=16)

axarr[0,3].set_xlim([xmin, xmax])
axarr[0,3].set_ylim([rhomin, rhomax])
axarr[0,3].plot(x, rho, linestyle='-', color="red")
axarr[0,3].scatter(x3, rho3, marker='.', s=16)
axarr[0,3].text(xmin + 0.03*xsize, rhomax - 0.1*rhosize, '(d) GANDALF (Combined)', color='black', size=14)
axarr[1,3].set_xlabel(r"$x$")
axarr[1,3].set_ylim([ymin, ymax])
axarr[1,3].scatter(x3, y3, color='black', marker='.', s=16)


plt.show()
fig.savefig('regularisation.pdf', dpi=50)
