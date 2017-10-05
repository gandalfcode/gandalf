#==============================================================================
# khitest.py
#==============================================================================
from gandalf.analysis.facade import *
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import rc
import time


rc('font', **{'family': 'serif', 'size' : 14})

# Set all plot limits
xmin = -0.5
xmax = 0.5
ymin = -0.5
ymax = 0.5
rhomin = 0.9
rhomax = 2.1

#loadsim('KHI-GRADH')
khi_gradh_sim = newsim('khi-gradh.dat')
setupsim()
run()
data0_0 = get_render_data('x', 'y', 'rho', sim=0, snap=5, res=256)
data0_1 = get_render_data('x', 'y', 'rho', sim=0, snap=10, res=256)
data0_2 = get_render_data('x', 'y', 'rho', sim=0, snap=15, res=256)
data0_3 = get_render_data('x', 'y', 'rho', sim=0, snap=20, res=256)
data0_4 = get_render_data('x', 'y', 'rho', sim=0, snap=25, res=256)

#loadsim('KHI-MFV-MOVING')
khi_mfv_sim = newsim('khi-mfv-moving.dat')
setupsim()
run()
data1_0 = get_render_data('x', 'y', 'rho', sim=1, snap=5, res=256)
data1_1 = get_render_data('x', 'y', 'rho', sim=1, snap=10, res=256)
data1_2 = get_render_data('x', 'y', 'rho', sim=1, snap=15, res=256)
data1_3 = get_render_data('x', 'y', 'rho', sim=1, snap=20, res=256)
data1_4 = get_render_data('x', 'y', 'rho', sim=1, snap=25, res=256)

#loadsim('KHI-MFV-MOVING')
khi_mfm_sim = newsim('khi-mfv-moving.dat')
khi_mfm_sim.SetParam('zero_mass_flux',1)
khi_mfm_sim.SetParam('run_id','KHI-MFM-MOVING')
setupsim()
run()
data2_0 = get_render_data('x', 'y', 'rho', sim=2, snap=5, res=256)
data2_1 = get_render_data('x', 'y', 'rho', sim=2, snap=10, res=256)
data2_2 = get_render_data('x', 'y', 'rho', sim=2, snap=15, res=256)
data2_3 = get_render_data('x', 'y', 'rho', sim=2, snap=20, res=256)
data2_4 = get_render_data('x', 'y', 'rho', sim=2, snap=25, res=256)


fig, ax = plt.subplots(3, 5, sharey='row', figsize=(15,9))
fig.subplots_adjust(hspace=0.01, wspace=0.01)
fig.subplots_adjust(bottom=0.0, top=1.0, left=0.0, right=1.0)


#fig = plt.figure(figsize=(10,3))
#fig.subplots_adjust(wspace=0.001,hspace=0.001)

#ax1.set_xlim([xmin, xmax])
#ax2.set_xlim([xmin, xmax])
#ax3.set_xlim([xmin, xmax])
#ax1.set_ylim([ymin, ymax])

#ax1 = fig.add_subplot(1,3,1,aspect='equal',xlim=[xmin,xmax],ylim=[ymin,ymax])
ax[0][0].text(14, 22, '(a) Grad-h SPH', color='white', size=14)
ax[0][0].set_axis_off()
ax[0][0].imshow(data0_0, interpolation='nearest', cmap=cm.jet, vmin=rhomin, vmax=rhomax)
ax[0][1].set_axis_off()
ax[0][1].imshow(data0_1, interpolation='nearest', cmap=cm.jet, vmin=rhomin, vmax=rhomax)
ax[0][2].set_axis_off()
ax[0][2].imshow(data0_2, interpolation='nearest', cmap=cm.jet, vmin=rhomin, vmax=rhomax)
ax[0][3].set_axis_off()
ax[0][3].imshow(data0_3, interpolation='nearest', cmap=cm.jet, vmin=rhomin, vmax=rhomax)
ax[0][4].set_axis_off()
ax[0][4].imshow(data0_4, interpolation='nearest', cmap=cm.jet, vmin=rhomin, vmax=rhomax)

ax[1][0].text(14, 22, '(b) MFV', color='white', size=14)
ax[1][0].set_axis_off()
ax[1][0].imshow(data1_0, interpolation='nearest', cmap=cm.jet, vmin=rhomin, vmax=rhomax)
ax[1][1].set_axis_off()
ax[1][1].imshow(data1_1, interpolation='nearest', cmap=cm.jet, vmin=rhomin, vmax=rhomax)
ax[1][2].set_axis_off()
ax[1][2].imshow(data1_2, interpolation='nearest', cmap=cm.jet, vmin=rhomin, vmax=rhomax)
ax[1][3].set_axis_off()
ax[1][3].imshow(data1_3, interpolation='nearest', cmap=cm.jet, vmin=rhomin, vmax=rhomax)
ax[1][4].set_axis_off()
ax[1][4].imshow(data1_4, interpolation='nearest', cmap=cm.jet, vmin=rhomin, vmax=rhomax)

ax[2][0].text(14, 22, '(c) MFM', color='white', size=14)
ax[2][0].set_axis_off()
ax[2][0].imshow(data2_0, interpolation='nearest', cmap=cm.jet, vmin=rhomin, vmax=rhomax)
ax[2][1].set_axis_off()
ax[2][1].imshow(data2_1, interpolation='nearest', cmap=cm.jet, vmin=rhomin, vmax=rhomax)
ax[2][2].set_axis_off()
ax[2][2].imshow(data2_2, interpolation='nearest', cmap=cm.jet, vmin=rhomin, vmax=rhomax)
ax[2][3].set_axis_off()
ax[2][3].imshow(data2_3, interpolation='nearest', cmap=cm.jet, vmin=rhomin, vmax=rhomax)
ax[2][4].set_axis_off()
ax[2][4].imshow(data2_4, interpolation='nearest', cmap=cm.jet, vmin=rhomin, vmax=rhomax)

#ax2 = fig.add_subplot(1,3,2,aspect='equal',sharey=ax1,xlim=[xmin,xmax],ylim=[ymin,ymax])
#ax[1].set_axis_off()
#ax[1].imshow(data1, interpolation='nearest', cmap=cm.jet, vmin=rhomin, vmax=rhomax)
#ax[1].text(10, 22, '(b) MFV', color='white', size=14)

#ax3 = fig.add_subplot(1,3,3,aspect='equal',sharey=ax1,xlim=[xmin,xmax],ylim=[ymin,ymax])
#ax[2].set_axis_off()
#ax[2].imshow(data2, interpolation='nearest', cmap=cm.jet, vmin=rhomin, vmax=rhomax)
#ax[2].text(10, 22, '(c) MFM', color='white', size=14)

#fig.tight_layout()
plt.show()
fig.savefig('khi.pdf', dpi=100)
