#==============================================================================
# bossbodenheimertest.py
# ...
#==============================================================================
from gandalf.analysis.facade import *
from gandalf.analysis.data_fetcher import *
import time
import matplotlib.pyplot as plt
import numpy as np
import math
from matplotlib import rc
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import AxesGrid
from scipy.optimize import curve_fit


#--------------------------------------------------------------------------------------------------
rc('font', **{'family': 'normal', 'weight' : 'bold', 'size' : 16})
rc('text', usetex=True)


numMinRes = 4000
numSimMax = 7
nmin = numMinRes
nmax = numMinRes*math.pow(2,numSimMax-1)
mmin = 0.04
mmax = 0.21

gradh_Nres = []
gradh_msink = []
gradh_smooth_msink = []

mfm_Nres = []
mfm_msink = []
mfm_smooth_msink = []


# Lowest resolution simulation (16 particles)
for i in range(numSimMax):
    Npart = int(numMinRes*math.pow(2,i))
    sim = newsim("bb-gradh.dat")
    sim.SetParam('Nhydro',Npart)
    sim.SetParam('rho_sink',1.0e-13)
    sim.SetParam('smooth_accretion',0)
    sim.SetParam('tend',0.03)
    setupsim()
    run()
    sinkMassData = get_data("m", type='star')
    gradh_Nres.append(Npart)
    gradh_msink.append(0.5*(sinkMassData[0] + sinkMassData[1]))


# Lowest resolution simulation (16 particles)
for i in range(numSimMax):
    Npart = int(numMinRes*math.pow(2,i))
    sim = newsim("bb-gradh.dat")
    sim.SetParam('Nhydro',Npart)
    sim.SetParam('rho_sink',1.0e-13)
    sim.SetParam('smooth_accretion',1)
    sim.SetParam('tend',0.03)
    setupsim()
    run()
    sinkMassData = get_data("m", type='star')
    gradh_smooth_msink.append(0.5*(sinkMassData[0] + sinkMassData[1]))


# Lowest resolution simulation (16 particles)
for i in range(numSimMax):
    Npart = int(numMinRes*math.pow(2,i))
    sim = newsim("bb-mfm.dat")
    sim.SetParam('Nhydro',Npart)
    sim.SetParam('rho_sink',1.0e-13)
    sim.SetParam('smooth_accretion',0)
    sim.SetParam('tend',0.03)
    setupsim()
    run()
    sinkMassData = get_data("m", type='star')
    mfm_Nres.append(Npart)
    mfm_msink.append(0.5*(sinkMassData[0] + sinkMassData[1]))


# Lowest resolution simulation (16 particles)
for i in range(numSimMax):
    Npart = int(numMinRes*math.pow(2,i))
    sim = newsim("bb-mfm.dat")
    sim.SetParam('Nhydro',Npart)
    sim.SetParam('rho_sink',1.0e-13)
    sim.SetParam('smooth_accretion',1)
    sim.SetParam('tend',0.03)
    setupsim()
    run()
    sinkMassData = get_data("m", type='star')
    mfm_smooth_msink.append(0.5*(sinkMassData[0] + sinkMassData[1]))



fig, ax = plt.subplots(figsize=(8.5,5.5))
fig.subplots_adjust(hspace=0.02, wspace=0.0001)
fig.subplots_adjust(bottom=0.1, top=0.98, left=0.12, right=0.98)


ax.set_xscale("log")
#ax.set_yscale("log")
ax.set_xlim([nmin, nmax])
ax.set_ylim([mmin, mmax])
ax.set_xlabel(r"$N$", fontsize=12)
ax.set_ylabel(r"$m$", fontsize=12)
ax.scatter(gradh_Nres, gradh_msink, color='black', marker='+', s=32.0, label='Gradh')
ax.plot(gradh_Nres, gradh_msink, linestyle=':', color='black')
ax.scatter(gradh_Nres, gradh_smooth_msink, color='black', marker='+', s=32.0, label='SPH, smooth')
ax.plot(gradh_Nres, gradh_smooth_msink, linestyle=':', color='black')
ax.scatter(mfm_Nres, mfm_msink, color='red', marker='^', s=32.0, label='MFM')
ax.plot(mfm_Nres, mfm_msink, linestyle=':', color='red')
ax.scatter(mfm_Nres, mfm_smooth_msink, color='red', marker='^', s=32.0, label='MFM, smooth')
ax.plot(mfm_Nres, mfm_smooth_msink, linestyle=':', color='red')
legend = ax.legend(loc='upper right', fontsize=12)


plt.show()
fig.savefig('bossbodenheimer.pdf', dpi=50)


block()
