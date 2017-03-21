#==============================================================================
# soundwave-L1error.py
# Computes L1 error norms for various resolutions of the adiabatic Sod test.
# Used to demonstrate scaling of error with resolution.
#==============================================================================
from gandalf.analysis.facade import *
from gandalf.analysis.compute import L1errornorm
import matplotlib.pyplot as plt
import time
import math


xmin = 8.0
xmax = 8192.0
ymin = 3.3e-12
ymax = 3.3e-6
numSimMax = 9


# Set empty lists to store results from each resolution
gradh_Nres = []
gradh_L1values = []

mfv_Nres = []
mfv_L1values = []



# Lowest resolution simulation (16 particles)
for i in range(8):
    Nhydro = int(8*math.pow(2,i))
    sim = newsim("soundwave-gradhsph.dat")
    sim.SetParam('Nhydro',Nhydro)
    setupsim()
    run()
    L1error = L1errornorm("soundwave", x="x", y="rho",
                          xmin=0.0, xmax=1.0, normalise=1.0)
    gradh_Nres.append(Nhydro)
    gradh_L1values.append(L1error)
    plot("x","rho")
    plotanalytical("x","rho")


# Lowest resolution simulation (16 particles)
for i in range(8):
    Nhydro = int(8*math.pow(2,i))
    sim = newsim("soundwave-mfv-moving.dat")
    sim.SetParam('Nhydro',Nhydro)
    setupsim()
    run()
    L1error = L1errornorm("soundwave", x="x", y="rho",
                          xmin=0.0, xmax=1.0, normalise=1.0)
    mfv_Nres.append(Nhydro)
    mfv_L1values.append(L1error)
    plot("x","rho")
    plotanalytical("x","rho")


fig, ax = plt.subplots(figsize=(6.5,5.5))
#fig.subplots_adjust(hspace=0.02, wspace=0.0001)
#fig.subplots_adjust(bottom=0.14, top=0.98, left=0.16, right=0.98)

ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim([xmin, xmax])
ax.set_ylim([ymin, ymax])
ax.set_xlabel(r"$N$", fontsize=12)
ax.set_ylabel(r"$L1$", fontsize=12)
ax.scatter(gradh_Nres, gradh_L1values, color='black', marker='+', s=8.0, label='Gradh-SPH')
ax.scatter(mfv_Nres, mfv_L1values, color='blue', marker='x', s=8.0, label='MFV')
ax.plot([xmin, xmax], [2.0e-4, 2.0e-4/(xmax-xmin)/(xmax-xmin)], linestyle=':', color='red')
ax.legend(fontsize=12)



print gradh_Nres
print gradh_L1values
print mfv_L1values

plt.show()
fig.savefig('soundwave-L1error.eps', dpi=50)


block()
