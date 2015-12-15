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
rhomin = 1.0
rhomax = 2.0


loadsim('KHI-GRADH')
data0 = get_render_data('x', 'y', 'rho', sim=0, snap=10, res=256)

loadsim('KHI-MFV-MOVING')
data1 = get_render_data('x', 'y', 'rho', sim=1, snap=10, res=256)

loadsim('KHI-MFV-STATIC')
data2 = get_render_data('x', 'y', 'rho', sim=2, snap=10, res=256)


fig, (ax1, ax2, ax3) = plt.subplots(1, 3, sharey='row', figsize=(9,3))
fig.subplots_adjust(hspace=0.02, wspace=0.02)
fig.subplots_adjust(bottom=0.0, top=1.0, left=0.0, right=1.0)
ax1.set_axis_off()
ax2.set_axis_off()
ax3.set_axis_off()

#fig = plt.figure(figsize=(10,3))
#fig.subplots_adjust(wspace=0.001,hspace=0.001)

#ax1.set_xlim([xmin, xmax])
#ax2.set_xlim([xmin, xmax])
#ax3.set_xlim([xmin, xmax])
#ax1.set_ylim([ymin, ymax])

#ax1 = fig.add_subplot(1,3,1,aspect='equal',xlim=[xmin,xmax],ylim=[ymin,ymax])
ax1.imshow(data0, interpolation='nearest', cmap=cm.jet, vmin = rhomin, vmax = rhomax)
ax1.text(14, 22, '(a) Grad-h SPH', color='white', size=14)

#ax2 = fig.add_subplot(1,3,2,aspect='equal',sharey=ax1,xlim=[xmin,xmax],ylim=[ymin,ymax])
ax2.imshow(data1, interpolation='nearest', cmap=cm.jet, vmin = rhomin, vmax = rhomax)
ax2.text(10, 22, '(b) MFV-moving', color='white', size=14)

#ax3 = fig.add_subplot(1,3,3,aspect='equal',sharey=ax1,xlim=[xmin,xmax],ylim=[ymin,ymax])
ax3.imshow(data2, interpolation='nearest', cmap=cm.jet, vmin = rhomin, vmax = rhomax)
ax3.text(10, 22, '(c) MFV-static', color='white', size=14)

#fig.tight_layout()
plt.show()
fig.savefig('khi.eps', dpi=50)
