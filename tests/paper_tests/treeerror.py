#==============================================================================
# treeerror.py
# Run the freefall collapse test using initial conditions specified in the
# file 'freefall.dat'.
#==============================================================================
from __future__ import division
from gandalf.analysis.facade import *
from gandalf.analysis.data_fetcher import *
from gandalf.analysis.compute import lagrangian_radii
from gandalf.analysis.SimBuffer import SimBuffer, BufferException
import time
import matplotlib.pyplot as plt
import numpy as np
import math
from matplotlib import rc
from mpl_toolkits.axes_grid1 import AxesGrid


#--------------------------------------------------------------------------------------------------
def ForceError(N, ax_bf, ay_bf, az_bf, ax_tree, ay_tree, az_tree):
    error = 0.0
    i = 0
    while i < N:
        da = [ax_tree[i] - ax_bf[i], ay_tree[i] - ay_bf[i], az_tree[i] - az_bf[i]]
        amagsqd = ax_bf[i]*ax_bf[i] + ay_bf[i]*ay_bf[i] + az_bf[i]*az_bf[i]
        errori = (da[0]*da[0] + da[1]*da[1] + da[2]*da[2]) #/amagsqd
        error = error + errori
        i = i + 1
        #print i, error, errori, da[0]*da[0] + da[1]*da[1] + da[2]*da[2], amagsqd
    error = error/N
    error = math.sqrt(error)
    print 'ERROR : ',N,error
    return error


#--------------------------------------------------------------------------------------------------
#rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('font', **{'family': 'normal', 'weight' : 'bold', 'size' : 16})
rc('text', usetex=True)

# Set all plot limits
tmin      = 0.001
tmax      = 1.0
rmin      = 0.0
rmax      = 1.0
rho       = 3.0/4.0/3.114157
stride    = 8
numSimMax = 6
thetamin  = 0.1
thetamax  = 1.0
sim_no    = 0
theta_values               = []
monopole_error_values      = []
monopole_runtimes          = []
fast_monopole_error_values = []
fast_monopole_runtimes     = []
quadrupole_error_values    = []
quadrupole_runtimes        = []


# Run brute-force test to get accurate forces
#--------------------------------------------------------------------------------------------------
bruteforce_sim = newsim("freefall.dat")
bruteforce_sim.SetParam('neib_search','bruteforce')
bruteforce_sim.SetParam('Nstepsmax',1)
bruteforce_sim.SetParam('run_id','FREEFALL1-BF')
setupsim()
run()

ax_bf = get_data("ax")#, sim=sim_no, snap=0)
ay_bf = get_data("ay")#, sim=sim_no, snap=0)
az_bf = get_data("az")#, sim=sim_no, snap=0)


# Fast monpole errors as a function of opening angle
#--------------------------------------------------------------------------------------------------
for i in range(numSimMax):
    sim_no = sim_no + 1
    exponent = 1.0*i/(numSimMax - 1)
    theta = thetamin*math.pow(thetamax/thetamin,exponent)
    thetamaxsqd = theta*theta
    print 'THETAMAX : ',theta,exponent
    sim = newsim("freefall.dat")
    sim.SetParam('thetamaxsqd',thetamaxsqd)
    sim.SetParam('Nstepsmax',1)
    setupsim()
    start = time.time()
    run()
    end = time.time()
    ax_tree = get_data("ax") #, sim=sim_no, snap=0)
    ay_tree = get_data("ay") #, sim=sim_no, snap=0)
    az_tree = get_data("az") #, sim=sim_no, snap=0)
    N = ax_tree.size
    error = ForceError(N, ax_bf, ay_bf, az_bf, ax_tree, ay_tree, az_tree)
    theta_values.append(theta)
    fast_monopole_error_values.append(error)
    fast_monopole_runtimes.append(end - start)


# Monopole errors as a function of opening angle
#--------------------------------------------------------------------------------------------------
for i in range(numSimMax):
    sim_no = sim_no + 1
    exponent = 1.0*i/(numSimMax - 1)
    theta = thetamin*math.pow(thetamax/thetamin,exponent)
    thetamaxsqd = theta*theta
    sim = newsim("freefall.dat")
    sim.SetParam('thetamaxsqd',thetamaxsqd)
    sim.SetParam('multipole',"monopole")
    sim.SetParam('Nstepsmax',1)
    setupsim()
    start = time.time()
    run()
    end = time.time()
    ax_tree = get_data("ax") #, sim=sim_no, snap=0)
    ay_tree = get_data("ay") #, sim=sim_no, snap=0)
    az_tree = get_data("az") #, sim=sim_no, snap=0)
    N = ax_tree.size
    error = ForceError(N, ax_bf, ay_bf, az_bf, ax_tree, ay_tree, az_tree)
    monopole_error_values.append(error)
    monopole_runtimes.append(end - start)



# Create matplotlib figure object with shared x-axis
#--------------------------------------------------------------------------------------------------
fig, axarr = plt.subplots(2, 1, figsize=(7,10))
#fig.subplots_adjust(hspace=0.001, wspace=0.1)
fig.subplots_adjust(bottom=0.09, top=0.97, left=0.11, right=0.96)

axarr[0].set_xscale("log")
axarr[0].set_yscale("log")
axarr[0].set_ylabel(r"$t_{_{\rm CPU}}$")
axarr[0].set_xlabel(r"$\delta\,{\bf a}$")
axarr[0].plot(fast_monopole_error_values, fast_monopole_runtimes, color="blue", linestyle='-', label='Fast monopole')
axarr[0].plot(monopole_error_values, monopole_runtimes, color="red", linestyle='-', label='Monopole')
axarr[0].legend(fontsize=12)

axarr[1].set_xscale("log")
axarr[1].set_yscale("log")
axarr[1].set_ylabel(r"$\delta\,{\bf a}$")
axarr[1].set_xlabel(r"$\theta$")
axarr[1].set_xlim([thetamin, thetamax])
axarr[1].plot(theta_values, fast_monopole_error_values, color="blue", linestyle='-', label='Fast monopole')
axarr[1].plot(theta_values, monopole_error_values, color="red", linestyle='-', label='Monopole')
axarr[1].legend(fontsize=12)

plt.show()
fig.savefig('treeerror.pdf', dpi=50)

# Prevent program from closing before showing plot window
block()
