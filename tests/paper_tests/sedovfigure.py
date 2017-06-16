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
xmin    = 0.002
xmax    = 0.37
rhomin  = 0.01
rhomax  = 4.9
vmin    = -0.3
vmax    = 2.2
umin    = 0.02
umax    = 2000.0
xsize   = xmax - xmin
rhosize = rhomax - rhomin
vsize   = vmax - vmin
usize   = umax - umin
stride  = 4


# Extract data from MFV simulation
mfvsim_l1 = newsim("sedov2d-mfv-moving.dat")
mfvsim_l1.SetParam('Nlevels',1)
mfvsim_l1.SetParam('time_step_limiter','none')
setupsim()
start = time.time()
run()
end = time.time()
print 'TIME : ',end - start
R0   = get_data('R')
rho0 = get_data('rho')
v0   = get_data('vR2d')

# Extract data from MFV simulation
mfvsim_l10 = newsim("sedov2d-mfv-moving.dat")
mfvsim_l10.SetParam('time_step_limiter','none')
mfvsim_l10.SetParam('Nlevels',5)
setupsim()
start = time.time()
run()
end = time.time()
print 'TIME : ',end - start
R1   = get_data('R')
rho1 = get_data('rho')
v1   = get_data('vR2d')

# Extract data from MFV simulation
mfvsim_l10_sm = newsim("sedov2d-mfv-moving.dat")
mfvsim_l10_sm.SetParam('time_step_limiter','simple')
mfvsim_l10_sm.SetParam('Nlevels',10)
mfvsim_l10_sm.SetParam('level_diff_max',1)
setupsim()
start = time.time()
run()
end = time.time()
print 'TIME : ',end - start
R2   = get_data('R')
rho2 = get_data('rho')
v2   = get_data('vR2d')

# Extract data from MFV simulation
mfvsim_l10_arepo = newsim("sedov2d-mfv-moving.dat")
mfvsim_l10_arepo.SetParam('time_step_limiter','conservative')
mfvsim_l10_arepo.SetParam('Nlevels',10)
mfvsim_l10_arepo.SetParam('level_diff_max',1)
setupsim()
start = time.time()
run()
end = time.time()
print 'TIME : ',end - start
R3   = get_data('R')
rho3 = get_data('rho')
v3   = get_data('vR2d')

# Extract data from MFV simulation
gradhsphsim_l1 = newsim("sedov2d-gradh.dat")
gradhsphsim_l1.SetParam('time_step_limiter','none')
gradhsphsim_l1.SetParam('Nlevels',1)
setupsim()
start = time.time()
run()
end = time.time()
print 'TIME : ',end - start
R4   = get_data('R')
rho4 = get_data('rho')
v4   = get_data('vR2d')

# Extract data from MFV simulation
gradhsphsim_l10 = newsim("sedov2d-gradh.dat")
gradhsphsim_l10.SetParam('time_step_limiter','none')
gradhsphsim_l10.SetParam('Nlevels',10)
gradhsphsim_l10.SetParam('level_diff_max',20)
setupsim()
start = time.time()
run()
end = time.time()
print 'TIME : ',end - start
R5   = get_data('R')
rho5 = get_data('rho')
v5   = get_data('vR2d')

# Extract data from MFV simulation
gradhsphsim_l10_sm = newsim("sedov2d-gradh.dat")
gradhsphsim_l10_sm.SetParam('time_step_limiter','simple')
gradhsphsim_l10_sm.SetParam('Nlevels',10)
gradhsphsim_l10_sm.SetParam('level_diff_max',1)
setupsim()
start = time.time()
run()
end = time.time()
print 'TIME : ',end - start
R6   = get_data('R')
rho6 = get_data('rho')
v6   = get_data('vR2d')

# Extract data from MFV simulation
gradhsphsim_l10_arepo = newsim("sedov2d-gradh.dat")
gradhsphsim_l10_arepo.SetParam('time_step_limiter','simple')
gradhsphsim_l10_arepo.SetParam('Nlevels',10)
gradhsphsim_l10_arepo.SetParam('level_diff_max',1)
setupsim()
start = time.time()
run()
end = time.time()
print 'TIME : ',end - start
R7   = get_data('R')
rho7 = get_data('rho')
v7   = get_data('vR2d')


# Extract data for analytical solution
rhodata = get_analytical_data("R","rho")
vxdata = get_analytical_data("R","vR2d")
udata = get_analytical_data("R","u")

# Create matplotlib figure object with shared x-axis
fig, axarr = plt.subplots(2, 3, sharex='col', sharey='row', figsize=(12,7))
fig.subplots_adjust(hspace=0.0001, wspace=0.0001)
fig.subplots_adjust(bottom=0.075, top=0.99, left=0.035, right=0.99)


# Grad-h SPH/MFV with 1 timestep level
axarr[0,0].set_ylabel(r"$\rho$", fontsize=14)
axarr[0,0].set_xlim([xmin, xmax])
axarr[0,0].set_ylim([rhomin, rhomax])
axarr[0,0].scatter(R4[::stride], rho4[::stride], color='black', marker='.', s=1.0, label='Grad-h, L=1')
axarr[0,0].plot(rhodata.x_data, rhodata.y_data, color="red")
axarr[0,0].text(xmin + 0.02*xsize, rhomax - 0.08*rhosize, '(a) Grad-h', color='black', size=14)
axarr[0,0].text(xmin + 0.82*xsize, rhomax - 0.08*rhosize, r'$L = 1$', color='black', size=14)
axarr[1,0].set_ylabel(r"$\rho$", fontsize=14)
axarr[1,0].set_xlabel(r"$x$", fontsize=14)
axarr[1,0].set_xlim([xmin, xmax])
axarr[1,0].set_ylim([rhomin, rhomax])
axarr[1,0].scatter(R0[::stride], rho0[::stride], color='black', marker='.', s=1.0, label='MFV, L=1')
axarr[1,0].plot(rhodata.x_data, rhodata.y_data, color="red")
axarr[1,0].text(xmin + 0.02*xsize, rhomax - 0.08*rhosize, '(b) MFV', color='black', size=14)
axarr[1,0].text(xmin + 0.82*xsize, rhomax - 0.08*rhosize, r'$L = 1$', color='black', size=14)


# MFV with 10 timestep levels
axarr[0,1].set_xlim([xmin, xmax])
axarr[0,1].set_ylim([rhomin, rhomax])
axarr[0,1].scatter(R5[::stride], rho5[::stride], color='black', marker='.', s=1.0, label='Grad-h, L=10')
axarr[0,1].plot(rhodata.x_data, rhodata.y_data, color="red")
axarr[0,1].text(xmin + 0.8*xsize, rhomax - 0.08*rhosize, r'$L = 10$', color='black', size=14)
#axarr[1,1].set_xlabel(r"$x$", fontsize=14)
#axarr[1,1].set_xlim([xmin, xmax])
#axarr[1,1].set_ylim([rhomin, rhomax])
#axarr[1,1].scatter(R1[::stride], rho1[::stride], color='black', marker='.', s=1.0, label='MFV, L=10')
#axarr[1,1].plot(rhodata.x_data, rhodata.y_data, color="red")


# MFV with 10 timestep levels + SM limiter
axarr[0,2].set_xlim([xmin, xmax])
axarr[0,2].set_ylim([rhomin, rhomax])
axarr[0,2].scatter(R6[::stride], rho6[::stride], color='black', marker='.', s=1.0, label='Grad-h, L=10 + SM09')
axarr[0,2].plot(rhodata.x_data, rhodata.y_data, color="red")
axarr[0,2].text(xmin + 0.6*xsize, rhomax - 0.08*rhosize, r'$L = 10 + {\rm SM09}$', color='black', size=14)
axarr[1,2].set_xlabel(r"$x$", fontsize=14)
axarr[1,2].set_xlim([xmin, xmax])
axarr[1,2].set_ylim([rhomin, rhomax])
axarr[1,2].scatter(R2[::stride], rho2[::stride], color='black', marker='.', s=1.0, label='MFV, L=10 + SM09')
axarr[1,2].text(xmin + 0.6*xsize, rhomax - 0.08*rhosize, r'$L = 10 + {\rm SM09}$', color='black', size=14)
axarr[1,2].plot(rhodata.x_data, rhodata.y_data, color="red")


# Gradh-SPH with 10 timestep levels
#axarr[0,2].set_xlim([xmin, xmax])
#axarr[0,2].set_ylim([rhomin, rhomax])
#axarr[0,2].scatter(R7[::stride], rho7[::stride], color='black', marker='.', s=1.0, label='Grad-h, L=10 + S10')
#axarr[0,2].plot(rhodata.x_data, rhodata.y_data, color="red")
#axarr[0,2].text(xmin + 0.65*xsize, rhomax - 0.08*rhosize, r'$L = 10 + {\rm S10}$', color='black', size=14)
axarr[1,1].set_xlabel(r"$x$", fontsize=14)
axarr[1,1].set_xlim([xmin, xmax])
axarr[1,1].set_ylim([rhomin, rhomax])
axarr[1,1].scatter(R3[::stride], rho3[::stride], color='black', marker='.', s=1.0, label='MFV, L=10 + S10')
axarr[1,1].plot(rhodata.x_data, rhodata.y_data, color="red")
axarr[1,1].text(xmin + 0.65*xsize, rhomax - 0.08*rhosize, r'$L = 10 + {\rm S10}$', color='black', size=14)


plt.show()
fig.savefig('sedov.pdf', dpi=50)
