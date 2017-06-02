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
    print 'FORCE ERROR : ',N,error
    return error


#--------------------------------------------------------------------------------------------------
def PotentialError(N, gpot_bf, gpot_tree):
    error = 0.0
    i = 0
    while i < N:
        dgpot = gpot_tree[i] - gpot_bf[i]
        gpotmagsqd = gpot_bf[i]*gpot_bf[i]
        errori = dgpot*dgpot
        error = error + errori
        i = i + 1
    error = error/N
    error = math.sqrt(error)
    print 'GPOT ERROR : ',N,error
    return error



#--------------------------------------------------------------------------------------------------
#rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('font', **{'family': 'normal', 'weight' : 'bold', 'size' : 16})
rc('text', usetex=True)


# Set all plot limits
Nhydro      = 16000
Nstepsmax   = 1
tmin        = 0.001
tmax        = 1.0
rmin        = 0.0
rmax        = 1.0
rho         = 3.0/4.0/3.14157
stride      = 8
numSimMax   = 8
thetamin    = 0.1
thetamax    = 1.0
macerrormin = 3.3333e-10
macerrormax = 3.3333e-3
errormin    = 3.3333e-8
errormax    = 3.3333e-2
sim_no      = 0

theta_values                        = []
monopole_geo_error_values           = []
monopole_geo_poterror_values        = []
monopole_geo_runtimes               = []
fast_monopole_geo_error_values      = []
fast_monopole_geo_poterror_values   = []
fast_monopole_geo_runtimes          = []
quadrupole_geo_error_values         = []
quadrupole_geo_poterror_values      = []
quadrupole_geo_runtimes             = []
fast_quadrupole_geo_error_values    = []
fast_quadrupole_geo_poterror_values = []
fast_quadrupole_geo_runtimes        = []

macerror_values = []
monopole_gadget_error_values           = []
monopole_gadget_poterror_values        = []
monopole_gadget_runtimes               = []
fast_monopole_gadget_error_values      = []
fast_monopole_gadget_poterror_values   = []
fast_monopole_gadget_runtimes          = []
quadrupole_gadget_error_values         = []
quadrupole_gadget_poterror_values      = []
quadrupole_gadget_runtimes             = []
fast_quadrupole_gadget_error_values    = []
fast_quadrupole_gadget_poterror_values = []
fast_quadrupole_gadget_runtimes        = []



# Run brute-force test to get accurate forces
#--------------------------------------------------------------------------------------------------
bruteforce_sim = newsim("freefall.dat")
bruteforce_sim.SetParam('neib_search','bruteforce')
bruteforce_sim.SetParam('Nstepsmax',Nstepsmax)
bruteforce_sim.SetParam('Nhydro',Nhydro)
bruteforce_sim.SetParam('run_id','FREEFALL1-BF')
setupsim()
run()

ax_bf   = get_data("ax")
ay_bf   = get_data("ay")
az_bf   = get_data("az")
gpot_bf = get_data("gpot")


# GEOMETRIC MAC
# Fast monopole errors as a function of opening angle
#--------------------------------------------------------------------------------------------------
for i in range(numSimMax):
    sim_no = sim_no + 1
    exponent = 1.0*i/(numSimMax - 1)
    theta = thetamin*math.pow(thetamax/thetamin,exponent)
    thetamaxsqd = theta*theta
    print 'THETAMAX : ',theta,exponent
    sim = newsim("freefall.dat")
    sim.SetParam('gravity_mac',"geometric")
    sim.SetParam('thetamaxsqd',thetamaxsqd)
    sim.SetParam('multipole',"fast_monopole")
    sim.SetParam('Nstepsmax',Nstepsmax)
    sim.SetParam('Nhydro',Nhydro)
    setupsim()
    start = time.time()
    run()
    end = time.time()
    ax_tree   = get_data("ax")
    ay_tree   = get_data("ay")
    az_tree   = get_data("az")
    gpot_tree = get_data("gpot")
    N         = ax_tree.size
    error     = ForceError(N, ax_bf, ay_bf, az_bf, ax_tree, ay_tree, az_tree)
    poterror  = PotentialError(N, gpot_bf, gpot_tree)
    theta_values.append(theta)
    fast_monopole_geo_error_values.append(error)
    fast_monopole_geo_poterror_values.append(poterror)
    fast_monopole_geo_runtimes.append(end - start)


# Monopole errors as a function of opening angle
#--------------------------------------------------------------------------------------------------
for i in range(numSimMax):
    sim_no = sim_no + 1
    exponent = 1.0*i/(numSimMax - 1)
    theta = thetamin*math.pow(thetamax/thetamin,exponent)
    thetamaxsqd = theta*theta
    sim = newsim("freefall.dat")
    sim.SetParam('gravity_mac',"geometric")
    sim.SetParam('thetamaxsqd',thetamaxsqd)
    sim.SetParam('multipole',"monopole")
    sim.SetParam('Nstepsmax',Nstepsmax)
    sim.SetParam('Nhydro',Nhydro)
    setupsim()
    start = time.time()
    run()
    end = time.time()
    ax_tree   = get_data("ax")
    ay_tree   = get_data("ay")
    az_tree   = get_data("az")
    gpot_tree = get_data("gpot")
    N         = ax_tree.size
    error     = ForceError(N, ax_bf, ay_bf, az_bf, ax_tree, ay_tree, az_tree)
    poterror  = PotentialError(N, gpot_bf, gpot_tree)
    monopole_geo_error_values.append(error)
    monopole_geo_poterror_values.append(poterror)
    monopole_geo_runtimes.append(end - start)


# Quadrupole errors as a function of opening angle
#--------------------------------------------------------------------------------------------------
for i in range(numSimMax):
    sim_no = sim_no + 1
    exponent = 1.0*i/(numSimMax - 1)
    theta = thetamin*math.pow(thetamax/thetamin,exponent)
    thetamaxsqd = theta*theta
    sim = newsim("freefall.dat")
    sim.SetParam('gravity_mac',"geometric")
    sim.SetParam('thetamaxsqd',thetamaxsqd)
    sim.SetParam('multipole',"quadrupole")
    sim.SetParam('Nstepsmax',Nstepsmax)
    sim.SetParam('Nhydro',Nhydro)
    setupsim()
    start = time.time()
    run()
    end = time.time()
    ax_tree   = get_data("ax")
    ay_tree   = get_data("ay")
    az_tree   = get_data("az")
    gpot_tree = get_data("gpot")
    N         = ax_tree.size
    error     = ForceError(N, ax_bf, ay_bf, az_bf, ax_tree, ay_tree, az_tree)
    poterror  = PotentialError(N, gpot_bf, gpot_tree)
    quadrupole_geo_error_values.append(error)
    quadrupole_geo_poterror_values.append(poterror)
    quadrupole_geo_runtimes.append(end - start)


# Fast quadrupole errors as a function of opening angle
#--------------------------------------------------------------------------------------------------
for i in range(numSimMax):
    sim_no = sim_no + 1
    exponent = 1.0*i/(numSimMax - 1)
    theta = thetamin*math.pow(thetamax/thetamin,exponent)
    thetamaxsqd = theta*theta
    sim = newsim("freefall.dat")
    sim.SetParam('gravity_mac',"geometric")
    sim.SetParam('thetamaxsqd',thetamaxsqd)
    sim.SetParam('multipole',"fast_quadrupole")
    sim.SetParam('Nstepsmax',Nstepsmax)
    sim.SetParam('Nhydro',Nhydro)
    setupsim()
    start = time.time()
    run()
    end = time.time()
    ax_tree   = get_data("ax")
    ay_tree   = get_data("ay")
    az_tree   = get_data("az")
    gpot_tree = get_data("gpot")
    N         = ax_tree.size
    error     = ForceError(N, ax_bf, ay_bf, az_bf, ax_tree, ay_tree, az_tree)
    poterror  = PotentialError(N, gpot_bf, gpot_tree)
    fast_quadrupole_geo_error_values.append(error)
    fast_quadrupole_geo_poterror_values.append(poterror)
    fast_quadrupole_geo_runtimes.append(end - start)


# GADGET2 MAC
# Fast monopole errors as a function of opening angle
#--------------------------------------------------------------------------------------------------
for i in range(numSimMax):
    sim_no = sim_no + 1
    exponent = 1.0*i/(numSimMax - 1)
    macerror = macerrormin*math.pow(macerrormax/macerrormin, exponent)
    print 'MACERROR : ',macerror
    sim = newsim("freefall.dat")
    sim.SetParam('gravity_mac',"gadget2")
    sim.SetParam('macerror',macerror)
    sim.SetParam('multipole',"fast_monopole")
    sim.SetParam('Nstepsmax',Nstepsmax)
    sim.SetParam('Nhydro',Nhydro)
    setupsim()
    start = time.time()
    run()
    end = time.time()
    ax_tree   = get_data("ax")
    ay_tree   = get_data("ay")
    az_tree   = get_data("az")
    gpot_tree = get_data("gpot")
    N         = ax_tree.size
    error     = ForceError(N, ax_bf, ay_bf, az_bf, ax_tree, ay_tree, az_tree)
    poterror  = PotentialError(N, gpot_bf, gpot_tree)
    macerror_values.append(macerror)
    fast_monopole_gadget_error_values.append(error)
    fast_monopole_gadget_poterror_values.append(poterror)
    fast_monopole_gadget_runtimes.append(end - start)


# Monopole errors as a function of opening angle
#--------------------------------------------------------------------------------------------------
for i in range(numSimMax):
    sim_no = sim_no + 1
    exponent = 1.0*i/(numSimMax - 1)
    macerror = macerrormin*math.pow(macerrormax/macerrormin, exponent)
    print 'MACERROR : ',macerror
    sim = newsim("freefall.dat")
    sim.SetParam('gravity_mac',"gadget2")
    sim.SetParam('macerror',macerror)
    sim.SetParam('multipole',"monopole")
    sim.SetParam('Nstepsmax',Nstepsmax)
    sim.SetParam('Nhydro',Nhydro)
    setupsim()
    start = time.time()
    run()
    end = time.time()
    ax_tree   = get_data("ax")
    ay_tree   = get_data("ay")
    az_tree   = get_data("az")
    gpot_tree = get_data("gpot")
    N         = ax_tree.size
    error     = ForceError(N, ax_bf, ay_bf, az_bf, ax_tree, ay_tree, az_tree)
    poterror  = PotentialError(N, gpot_bf, gpot_tree)
    monopole_gadget_error_values.append(error)
    monopole_gadget_poterror_values.append(poterror)
    monopole_gadget_runtimes.append(end - start)


# Quadrupole errors as a function of opening angle
#--------------------------------------------------------------------------------------------------
for i in range(numSimMax):
    sim_no = sim_no + 1
    exponent = 1.0*i/(numSimMax - 1)
    macerror = macerrormin*math.pow(macerrormax/macerrormin, exponent)
    print 'MACERROR : ',macerror
    sim = newsim("freefall.dat")
    sim.SetParam('gravity_mac',"gadget2")
    sim.SetParam('macerror',macerror)
    sim.SetParam('multipole',"quadrupole")
    sim.SetParam('Nstepsmax',Nstepsmax)
    sim.SetParam('Nhydro',Nhydro)
    setupsim()
    start = time.time()
    run()
    end = time.time()
    ax_tree   = get_data("ax")
    ay_tree   = get_data("ay")
    az_tree   = get_data("az")
    gpot_tree = get_data("gpot")
    N         = ax_tree.size
    error     = ForceError(N, ax_bf, ay_bf, az_bf, ax_tree, ay_tree, az_tree)
    poterror  = PotentialError(N, gpot_bf, gpot_tree)
    quadrupole_gadget_error_values.append(error)
    quadrupole_gadget_poterror_values.append(poterror)
    quadrupole_gadget_runtimes.append(end - start)


# Fast quadrupole errors as a function of opening angle
#--------------------------------------------------------------------------------------------------
for i in range(numSimMax):
    sim_no = sim_no + 1
    exponent = 1.0*i/(numSimMax - 1)
    macerror = macerrormin*math.pow(macerrormax/macerrormin, exponent)
    print 'MACERROR : ',macerror
    sim = newsim("freefall.dat")
    sim.SetParam('gravity_mac',"gadget2")
    sim.SetParam('macerror',macerror)
    sim.SetParam('multipole',"fast_quadrupole")
    sim.SetParam('Nstepsmax',Nstepsmax)
    sim.SetParam('Nhydro',Nhydro)
    setupsim()
    start = time.time()
    run()
    end = time.time()
    ax_tree   = get_data("ax")
    ay_tree   = get_data("ay")
    az_tree   = get_data("az")
    gpot_tree = get_data("gpot")
    N         = ax_tree.size
    error     = ForceError(N, ax_bf, ay_bf, az_bf, ax_tree, ay_tree, az_tree)
    poterror  = PotentialError(N, gpot_bf, gpot_tree)
    fast_quadrupole_gadget_error_values.append(error)
    fast_quadrupole_gadget_poterror_values.append(poterror)
    fast_quadrupole_gadget_runtimes.append(end - start)





# Create matplotlib figure object with shared x-axis
#--------------------------------------------------------------------------------------------------
fig, axarr = plt.subplots(1, 2, figsize=(14,7), sharey='col')
fig.subplots_adjust(hspace=0.001, wspace=0.0001)
fig.subplots_adjust(bottom=0.09, top=0.97, left=0.09, right=0.96)

#axarr[0].set_xscale("log")
#axarr[0].set_yscale("log")
#axarr[0].set_ylabel(r"$t_{_{\rm CPU}}$")
#axarr[0].set_xlabel(r"$\delta\,{\bf a}$")
#axarr[0].plot(monopole_geo_poterror_values, monopole_geo_runtimes, color="red", linestyle='-', label='Monopole')
#axarr[0].plot(fast_monoqpole_geo_poterror_values, fast_monopole_geo_runtimes, color="blue", linestyle='-', label='Fast monopole')
#axarr[0].plot(quadrupole_geo_poterror_values, quadrupole_geo_runtimes, color="black", linestyle='-', label='Quadrupole')
#axarr[0].plot(fast_quadrupole_geo_poterror_values, fast_quadrupole_geo_runtimes, color="green", linestyle='-', label='Fast quadrupole')
#axarr[0].legend(fontsize=12)

axarr[0].set_xscale("log")
axarr[0].set_yscale("log")
axarr[0].set_ylabel(r"$\delta\,{\bf a}$")
axarr[0].set_xlabel(r"$\theta$")
axarr[0].set_xlim([thetamin, thetamax])
axarr[0].set_ylim([errormin, errormax])
axarr[0].plot(theta_values, monopole_geo_error_values, color="red", linestyle='-', label='Monopole')
axarr[0].plot(theta_values, fast_monopole_geo_error_values, color="blue", linestyle=':', label='Fast monopole')
axarr[0].plot(theta_values, quadrupole_geo_error_values, color="black", linestyle='--', label='Quadrupole')
#axarr[0].plot(theta_values, fast_quadrupole_geo_error_values, color="green", linestyle='-', label='Fast quadrupole')
axarr[0].scatter(theta_values, monopole_geo_error_values, color="red")
axarr[0].scatter(theta_values, fast_monopole_geo_error_values, color="blue")
axarr[0].scatter(theta_values, quadrupole_geo_error_values, color="black")
#axarr[0].scatter(theta_values, fast_quadrupole_geo_error_values, color="green")
legend0 = axarr[0].legend(loc='upper left', fontsize=12)

axarr[1].set_xscale("log")
axarr[1].set_yscale("log")
#axarr[1].set_ylabel(r"$\delta\,\phi$")
axarr[1].set_xlabel(r"$\alpha_{_{\rm MAC}}$")
axarr[1].set_xlim([macerrormin, macerrormax])
axarr[1].set_ylim([errormin, errormax])
axarr[1].plot(macerror_values, monopole_gadget_error_values, color="red", linestyle='-', label='Monopole')
axarr[1].plot(macerror_values, fast_monopole_gadget_error_values, color="blue", linestyle='-', label='Fast monopole')
axarr[1].plot(macerror_values, quadrupole_gadget_error_values, color="black", linestyle='-', label='Quadrupole')
#axarr[1].plot(macerror_values, fast_quadrupole_gadget_error_values, color="green", linestyle='-', label='Fast quadrupole')
axarr[1].scatter(macerror_values, monopole_gadget_error_values, color="red")
axarr[1].scatter(macerror_values, fast_monopole_gadget_error_values, color="blue")
axarr[1].scatter(macerror_values, quadrupole_gadget_error_values, color="black")
#axarr[1].scatter(macerror_values, fast_quadrupole_gadget_error_values, color="green")
legend1 = axarr[0].legend(loc='upper left', fontsize=12)

plt.show()
fig.savefig('treeerror.pdf', dpi=50)

# Prevent program from closing before showing plot window
block()
