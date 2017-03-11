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


xmin = 128.0
xmax = 16384.0
ymin = 3.3e-12
ymax = 3.3e-6
numSimMax = 9


# Make simple N^-2 line
Nmax = 1000
x = np.arange(xmin, xmax, 2.0)
ydirect = 0.00000025*x*x
ytree = 0.000005*x*np.log(x)

# Set empty lists to store results from each resolution
Nvalues = []
bruteforce_time = []
kdtree_time = []



# Brute-force computation
for i in range(numSimMax):
    Nhydro = int(128*math.pow(2,i))
    sim = newsim("freefall.dat")
    sim.SetParam('Nhydro',Nhydro)
    sim.SetParam('tsnapfirst',2.0)
    sim.SetParam("Nstepsmax",4)
    sim.SetParam("neib_search","bruteforce")
    setupsim()

    start = time.time()
    run()
    end = time.time()

    Nvalues.append(Nhydro)
    bruteforce_time.append(end - start)



# KD-tree computation
for i in range(numSimMax):
    Nhydro = int(128*math.pow(2,i))
    sim = newsim("freefall.dat")
    sim.SetParam('Nhydro',Nhydro)
    sim.SetParam('tsnapfirst',2.0)
    sim.SetParam("Nstepsmax",4)
    sim.SetParam("neib_search","kdtree")
    setupsim()

    start = time.time()
    run()
    end = time.time()

    kdtree_time.append(end - start)


fig, ax = plt.subplots(figsize=(7,5))
#fig.subplots_adjust(hspace=0.02, wspace=0.0001)
#fig.subplots_adjust(bottom=0.14, top=0.98, left=0.16, right=0.98)

ax.set_xscale("log")
ax.set_yscale("log")
#ax.set_xlim([xmin, xmax])
#ax.set_ylim([ymin, ymax])
ax.set_xlabel(r"$N$", fontsize=12)
ax.set_ylabel(r"$t_{_{\rm CPU}}$", fontsize=12)
ax.scatter(Nvalues, bruteforce_time, color='black', marker='+', s=8.0, label='Direct sum')
ax.scatter(Nvalues, kdtree_time, color='black', marker='^', s=8.0, label='KD-tree')
ax.plot(x, ydirect, linestyle=':', color='red')
ax.plot(x, ytree, linestyle=':', color='red')
ax.legend(fontsize=12)



plt.show()
fig.savefig('treescaling.eps', dpi=50)


block()
