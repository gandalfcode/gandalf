#==============================================================================
# soundwave-L1error.py
# Computes L1 error norms for various resolutions of the adiabatic Sod test.
# Used to demonstrate scaling of error with resolution.
#==============================================================================
from gandalf.analysis.facade import *
from gandalf.analysis.compute import L1errornorm
from matplotlib import rc
import matplotlib.pyplot as plt
import numpy as np
import time
import math


#rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('font', **{'family': 'normal', 'weight' : 'bold', 'size' : 16})
rc('text', usetex=True)


numBruteforceSimMax = 10
numTreeSimMax = 17
numSimMax = max(numBruteforceSimMax, numTreeSimMax)
Nstepsmax = 1
thetamaxsqd = 0.15
Nhydromin = 128
xmin = 0.999*Nhydromin
xmax = 1.001*Nhydromin*math.pow(2,numSimMax-1)
ymin = 6.666e-4
ymax = 3.333e3


# Make simple N, N^2 and NlogN lines
Nmax = 1000
x = np.arange(xmin, xmax, 2.0)
ydirect = 0.0000001*x*x
ytree = 0.000002*x*np.log(x)
yN = 0.000009*x


# Set empty lists to store results from each resolution
bruteforce_N = []
bruteforce_time = []
kdtree_N = []
kdtree_time = []



# Brute-force computation
for i in range(numBruteforceSimMax):
    Nhydro = int(Nhydromin*math.pow(2,i))
    print i,numSimMax,Nhydro
    sim = newsim("freefall.dat")
    sim.SetParam('Nhydro',Nhydro)
    sim.SetParam('tsnapfirst',2.0)
    sim.SetParam('ntreebuildstep',1)
    sim.SetParam("Nstepsmax",Nstepsmax)
    sim.SetParam("neib_search","bruteforce")
    setupsim()

    start = time.time()
    run()
    end = time.time()
    dt = end - start  #sim.GetBlockTime("SPH_ALL_FORCES")

    bruteforce_N.append(Nhydro)
    bruteforce_time.append(dt)



# KD-tree computation
for i in range(numTreeSimMax):
    Nhydro = int(Nhydromin*math.pow(2,i))
    sim = newsim("freefall.dat")
    sim.SetParam('Nhydro',Nhydro)
    sim.SetParam('tsnapfirst',2.0)
    sim.SetParam('ntreebuildstep',1)
    sim.SetParam('thetamaxsqd',thetamaxsqd)
    sim.SetParam('multipole','fast_monopole')
    sim.SetParam('Nleafmax',3)
    sim.SetParam("Nstepsmax",Nstepsmax)
    sim.SetParam("neib_search","kdtree")
    setupsim()

    start = time.time()
    run()
    end = time.time()
    dt = end - start  #sim.GetBlockTime("SPH_ALL_FORCES")

    kdtree_N.append(Nhydro)
    kdtree_time.append(dt)



fig, ax = plt.subplots(figsize=(7,5))
#fig.subplots_adjust(hspace=0.02, wspace=0.0001)
fig.subplots_adjust(bottom=0.11, top=0.98, left=0.115, right=0.98)

ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim([xmin, xmax])
ax.set_ylim([ymin, ymax])
ax.set_xlabel(r"$N$", fontsize=16)
ax.set_ylabel(r"$t_{_{\rm CPU}}$", fontsize=16)
ax.scatter(bruteforce_N, bruteforce_time, color='black', marker='+', s=16.0, label='Direct sum')
ax.scatter(kdtree_N, kdtree_time, color='black', marker='^', s=16.0, label='KD-tree')
ax.plot(x, yN, linestyle='-.', color='red', label='${\cal O}(N)$')
ax.plot(x, ydirect, linestyle=':', color='red', label='${\cal O}(N^2)$')
ax.plot(x, ytree, linestyle='--', color='red', label='${\cal O}(N\,\log{N})$')
legend = ax.legend(loc='upper left', fontsize=16)



plt.show()
fig.savefig('treescaling.pdf', dpi=50)


block()
