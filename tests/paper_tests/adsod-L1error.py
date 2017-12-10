#==============================================================================
# adsod-L1error.py
# Computes L1 error norms for various resolutions of the adiabatic Sod test.
# Used to demonstrate scaling of error with resolution.
#==============================================================================
from gandalf.analysis.facade import *
from gandalf.analysis.compute import L1errornorm
import matplotlib.pyplot as plt
import numpy as np
import time
import math


xmin = 33.0
xmax = 3333.0
ymin = 3.3e-4
ymax = 1.666666e-1
numSims = 7


# Make simple N^-2 line
Nmax = 1000
x_ideal = np.arange(xmin, xmax, 2.0)
y_ideal = 0.5*ymax*xmin/x_ideal


# Set empty lists to store results from each resolution
gradh_Nres = []
gradh_L1values = []

mfv_Nres = []
mfv_L1values = []

mfm_Nres = []
mfm_L1values = []

mfv_static_values = []
mfv_static_L1values = []


# Lowest resolution simulation (16 particles)
for i in range(numSims):
    Nlhs = int(32*math.pow(2,i))
    Nrhs = int(8*math.pow(2,i))
    sim = newsim("adsod-gradhsph.dat")
    sim.SetParam('Nlattice1[0]',Nlhs)
    sim.SetParam('Nlattice2[0]',Nrhs)
    sim.SetParam('riemann_solver','hllc')
    setupsim()
    run()
    gradh_Nres.append(Nlhs+Nrhs)
    gradh_L1values.append(L1errornorm("shocktube","x","rho",-10.0,10.0))
    #plot("x","rho")
    #plotanalytical("x","rho")


# Lowest resolution simulation (16 particles)
for i in range(numSims):
    Nlhs = int(20*math.pow(2,i))
    Nrhs = int(20*math.pow(2,i))
    sim = newsim("adsod-mfv-moving.dat")
    sim.SetParam('Nlattice1[0]',Nlhs)
    sim.SetParam('Nlattice2[0]',Nrhs)
    sim.SetParam('riemann_solver','hllc')
    setupsim()
    run()
    mfv_Nres.append(Nlhs+Nrhs)
    mfv_L1values.append(L1errornorm("shocktube","x","rho",-10.0,10.0))
    #plot("x","rho")
    #plotanalytical("x","rho")


# Lowest resolution simulation (16 particles)
for i in range(numSims):
    Nlhs = int(20*math.pow(2,i))
    Nrhs = int(20*math.pow(2,i))
    sim = newsim("adsod-mfm-moving.dat")
    sim.SetParam('Nlattice1[0]',Nlhs)
    sim.SetParam('Nlattice2[0]',Nrhs)
    sim.SetParam('riemann_solver','hllc')
    setupsim()
    run()
    mfm_Nres.append(Nlhs+Nrhs)
    mfm_L1values.append(L1errornorm("shocktube","x","rho",-10.0,10.0))
    #plot("x","rho")
    #plotanalytical("x","rho")



fig, ax = plt.subplots(figsize=(7.5,6.5))
fig.subplots_adjust(hspace=0.02, wspace=0.0001)
fig.subplots_adjust(bottom=0.1, top=0.98, left=0.12, right=0.98)


ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim([xmin, xmax])
ax.set_ylim([ymin, ymax])
ax.set_xlabel(r"$N$", fontsize=12)
ax.set_ylabel(r"$L1$", fontsize=12)
ax.scatter(gradh_Nres, gradh_L1values, color='black', marker='+', s=32.0, label='Gradh-SPH')
ax.scatter(mfm_Nres, mfm_L1values, color='red', marker='^', s=32.0, label='MFM')
ax.scatter(mfv_Nres, mfv_L1values, color='blue', marker='x', s=32.0, label='MFV')
#ax.scatter(mfv_static_Nres, mfv_static_L1values, color='green', marker='*', s=16.0, label='MFV-static')
ax.plot(x_ideal, y_ideal, linestyle=':', color='red', label='$\propto N^{-1}$')
legend = ax.legend(loc='upper right', fontsize=12)

print gradh_Nres
print gradh_L1values


plt.show()
fig.savefig('adsod-L1error.pdf', dpi=50)
