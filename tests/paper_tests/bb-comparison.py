#==============================================================================
# bossbodenheimertest.py
# ...
#==============================================================================
from gandalf.analysis.facade import *
from gandalf.analysis.data_fetcher import *
import time
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import math
from matplotlib import rc
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import AxesGrid
from scipy.optimize import curve_fit


#--------------------------------------------------------------------------------------------------
rc('font', **{'family': 'normal', 'weight' : 'bold', 'size' : 16})
rc('text', usetex=True)


Npart = 16000
rhomin = 1.0e-18
rhomax = 1.0e-12
xmin = -0.007
xmax = 0.007
ymin = -0.007
ymax = 0.007
xsize = xmax - xmin
ysize = ymax - ymin



# Lowest resolution simulation (16 particles)
gradhsim = newsim("bb-gradh.dat")
gradhsim.SetParam('Nhydro',Npart)
gradhsim.SetParam('rho_sink',2.0e-13)
gradhsim.SetParam('tend',0.037)
setupsim()
run()

gradh_data_0 = get_render_data('x', 'y', 'rho', sim=0, snap=12, res=256, coordlimits=(xmin, xmax, ymin, ymax))
gradh_data_1 = get_render_data('x', 'y', 'rho', sim=0, snap=20, res=256, coordlimits=(xmin, xmax, ymin, ymax))
gradh_data_2 = get_render_data('x', 'y', 'rho', sim=0, snap=28, res=256, coordlimits=(xmin, xmax, ymin, ymax))
gradh_data_3 = get_render_data('x', 'y', 'rho', sim=0, snap=36, res=256, coordlimits=(xmin, xmax, ymin, ymax))
gradh_data_4 = get_render_data('x', 'y', 'rho', sim=0, snap=44, res=256, coordlimits=(xmin, xmax, ymin, ymax))
gradh_data_5 = get_render_data('x', 'y', 'rho', sim=0, snap=52, res=256, coordlimits=(xmin, xmax, ymin, ymax))



# Lowest resolution simulation (16 particles)
mfmsim = newsim("bb-mfm.dat")
mfmsim.SetParam('Nhydro',Npart)
mfmsim.SetParam('rho_sink',2.0e-13)
mfmsim.SetParam('tend',0.037)
setupsim()
run()

mfm_data_0 = get_render_data('x', 'y', 'rho', sim=1, snap=12, res=256, coordlimits=(xmin, xmax, ymin, ymax))
mfm_data_1 = get_render_data('x', 'y', 'rho', sim=1, snap=20, res=256, coordlimits=(xmin, xmax, ymin, ymax))
mfm_data_2 = get_render_data('x', 'y', 'rho', sim=1, snap=28, res=256, coordlimits=(xmin, xmax, ymin, ymax))
mfm_data_3 = get_render_data('x', 'y', 'rho', sim=1, snap=36, res=256, coordlimits=(xmin, xmax, ymin, ymax))
mfm_data_4 = get_render_data('x', 'y', 'rho', sim=1, snap=44, res=256, coordlimits=(xmin, xmax, ymin, ymax))
mfm_data_5 = get_render_data('x', 'y', 'rho', sim=1, snap=52, res=256, coordlimits=(xmin, xmax, ymin, ymax))



fig, ax = plt.subplots(2, 6, sharey='row', figsize=(12,4))
fig.subplots_adjust(hspace=0.01, wspace=0.01)
fig.subplots_adjust(bottom=0.01, top=0.99, left=0.01, right=0.99)

ax[0][0].set_axis_off()
ax[0][0].imshow(gradh_data_0, interpolation='nearest', cmap=cm.jet, norm=LogNorm(vmin=rhomin, vmax=rhomax), extent=(xmin, xmax, ymin, ymax))
ax[0][0].text(xmin + 0.02*xsize, ymax - 0.08*ysize, '(a) Grad-h', color='white', size=10)
ax[0][0].text(xmin + 0.56*xsize, ymax - 0.08*ysize, r'$t = 0.016\,{\rm Myr}$', color='white', size=10)
ax[0][1].set_axis_off()
ax[0][1].imshow(gradh_data_1, interpolation='nearest', cmap=cm.jet, norm=LogNorm(vmin=rhomin, vmax=rhomax), extent=(xmin, xmax, ymin, ymax))
ax[0][1].text(xmin + 0.56*xsize, ymax - 0.08*ysize, r'$t = 0.02\,{\rm Myr}$', color='white', size=10)
ax[0][2].set_axis_off()
ax[0][2].imshow(gradh_data_2, interpolation='nearest', cmap=cm.jet, norm=LogNorm(vmin=rhomin, vmax=rhomax), extent=(xmin, xmax, ymin, ymax))
ax[0][2].text(xmin + 0.56*xsize, ymax - 0.08*ysize, r'$t = 0.024\,{\rm Myr}$', color='white', size=10)
ax[0][3].set_axis_off()
ax[0][3].imshow(gradh_data_3, interpolation='nearest', cmap=cm.jet, norm=LogNorm(vmin=rhomin, vmax=rhomax), extent=(xmin, xmax, ymin, ymax))
ax[0][3].text(xmin + 0.56*xsize, ymax - 0.08*ysize, r'$t = 0.028\,{\rm Myr}$', color='white', size=10)
ax[0][4].set_axis_off()
ax[0][4].imshow(gradh_data_4, interpolation='nearest', cmap=cm.jet, norm=LogNorm(vmin=rhomin, vmax=rhomax), extent=(xmin, xmax, ymin, ymax))
ax[0][4].text(xmin + 0.56*xsize, ymax - 0.08*ysize, r'$t = 0.032\,{\rm Myr}$', color='white', size=10)
ax[0][5].set_axis_off()
ax[0][5].imshow(gradh_data_5, interpolation='nearest', cmap=cm.jet, norm=LogNorm(vmin=rhomin, vmax=rhomax), extent=(xmin, xmax, ymin, ymax))
ax[0][5].text(xmin + 0.56*xsize, ymax - 0.08*ysize, r'$t = 0.036\,{\rm Myr}$', color='white', size=10)

ax[1][0].set_axis_off()
ax[1][0].imshow(mfm_data_0, interpolation='nearest', cmap=cm.jet, norm=LogNorm(vmin=rhomin, vmax=rhomax), extent=(xmin, xmax, ymin, ymax))
ax[1][0].text(xmin + 0.02*xsize, ymax - 0.08*ysize, '(b) MFM', color='white', size=10)
ax[1][1].set_axis_off()
ax[1][1].imshow(mfm_data_1, interpolation='nearest', cmap=cm.jet, norm=LogNorm(vmin=rhomin, vmax=rhomax), extent=(xmin, xmax, ymin, ymax))
ax[1][2].set_axis_off()
ax[1][2].imshow(mfm_data_2, interpolation='nearest', cmap=cm.jet, norm=LogNorm(vmin=rhomin, vmax=rhomax), extent=(xmin, xmax, ymin, ymax))
ax[1][3].set_axis_off()
ax[1][3].imshow(mfm_data_3, interpolation='nearest', cmap=cm.jet, norm=LogNorm(vmin=rhomin, vmax=rhomax), extent=(xmin, xmax, ymin, ymax))
ax[1][4].set_axis_off()
ax[1][4].imshow(mfm_data_4, interpolation='nearest', cmap=cm.jet, norm=LogNorm(vmin=rhomin, vmax=rhomax), extent=(xmin, xmax, ymin, ymax))
ax[1][5].set_axis_off()
ax[1][5].imshow(mfm_data_5, interpolation='nearest', cmap=cm.jet, norm=LogNorm(vmin=rhomin, vmax=rhomax), extent=(xmin, xmax, ymin, ymax))



plt.show()
fig.savefig('bb-comparison.pdf', dpi=50)


block()
