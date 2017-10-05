#==============================================================================
# soundwave-L1error.py
# Computes L1 error norms for various resolutions of the adiabatic Sod test.
# Used to demonstrate scaling of error with resolution.
#==============================================================================
from gandalf.analysis.facade import *
from gandalf.analysis.compute import L1errornorm
import matplotlib.pyplot as plt
import numpy as np
import time
import math


xmin = 8.0
xmax = 4096.0
ymin = 3.3333e-12
ymax = 3.3333e-6
numSimMax = 8


# Make simple N^-2 line
Nmax = 1000
x_ideal = np.arange(xmin, xmax, 2.0)
y_ideal = 0.5*ymax*xmin*xmin/x_ideal/x_ideal


# Set empty lists to store results from each resolution
gradh_Nres = []
gradh_L1values = []

# Set empty lists to store results from each resolution
gradh_unnorm_Nres = []
gradh_unnorm_L1values = []

mfv_Nres = []
mfv_L1values = []

mfm_Nres = []
mfm_L1values = []

#mfv_static_Nres = []
#mfv_static_L1values = []


# Grad-h SPH
for i in range(numSimMax):
    Nhydro = int(16*math.pow(2,i))
    sim = newsim("soundwave-gradhsph.dat")
    sim.SetParam('Nhydro',Nhydro)
    setupsim()
    run()
    L1error = L1errornorm("soundwave", x="x", y="rho", xmin=0.0, xmax=1.0, normalise=1.0)
    gradh_Nres.append(Nhydro)
    gradh_L1values.append(L1error)
    L1error = L1errornorm("soundwave", x="x", y="rho", xmin=0.0, xmax=1.0)
    gradh_unnorm_Nres.append(Nhydro)
    gradh_unnorm_L1values.append(L1error)
    #plot("x","rho")
    #plotanalytical("x","rho")


# MFV
for i in range(numSimMax):
    Nhydro = int(16*math.pow(2,i))
    sim = newsim("soundwave-mfv-moving.dat")
    sim.SetParam('Nhydro',Nhydro)
    setupsim()
    run()
    L1error = L1errornorm("soundwave", x="x", y="rho", xmin=0.0, xmax=1.0, normalise=1.0)
    mfv_Nres.append(Nhydro)
    mfv_L1values.append(L1error)
    #plot("x","rho")
    #plotanalytical("x","rho")

# MFV static
#for i in range(numSimMax):
#    Nhydro = int(16*math.pow(2,i))
#    sim = newsim("soundwave-mfv-moving.dat")
#    sim.SetParam('static_particles',1)
#    sim.SetParam('Nhydro',Nhydro)
#    setupsim()
#    run()
#    L1error = L1errornorm("soundwave", x="x", y="rho",
#                          xmin=0.0, xmax=1.0, normalise=1.0)
#    mfv_static_Nres.append(Nhydro)
#    mfv_static_L1values.append(L1error)
#    plot("x","rho")
#    plotanalytical("x","rho")

# MFM
for i in range(numSimMax):
    Nhydro = int(16*math.pow(2,i))
    sim = newsim("soundwave-mfv-moving.dat")
    sim.SetParam('zero_mass_flux',1)
    sim.SetParam('Nhydro',Nhydro)
    setupsim()
    run()
    L1error = L1errornorm("soundwave", x="x", y="rho", xmin=0.0, xmax=1.0, normalise=1.0)
    mfm_Nres.append(Nhydro)
    mfm_L1values.append(L1error)
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
ax.scatter(gradh_Nres, gradh_L1values, color='black', marker='+', s=48.0, label='Gradh-SPH (normalised)')
#ax.scatter(gradh_unnorm_Nres, gradh_unnorm_L1values, color='black', marker='o', s=48.0, label='Gradh-SPH')
ax.scatter(mfm_Nres, mfm_L1values, color='red', marker='^', s=48.0, label='MFM')
ax.scatter(mfv_Nres, mfv_L1values, color='blue', marker='x', s=48.0, label='MFV')
#ax.scatter(mfv_static_Nres, mfv_static_L1values, color='green', marker='*', s=16.0, label='MFV-static')
ax.plot(x_ideal, y_ideal, linestyle=':', color='red', label='$\propto N^{-2}$')
legend = ax.legend(loc='upper right', fontsize=12)


print gradh_Nres
print gradh_L1values
print mfv_L1values

plt.show()
fig.savefig('soundwave-L1error.pdf', dpi=50)


block()
