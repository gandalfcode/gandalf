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
        errori = (da[0]*da[0] + da[1]*da[1] + da[2]*da[2]) /amagsqd
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
Nhydro       = 32000
Nstepsmax    = 1
tmin         = 0.003333
tmax         = 3.33333
rmin         = 0.0
rmax         = 1.0
rho          = 3.0/4.0/3.14157
stride       = 8
numSimMax    = 9
Nleafmax     = 6
thetamin     = 0.05
thetamax     = 1.1
thetasize    = thetamax - thetamin
macerrormin  = 1.5e-7
macerrormax  = 8.5e-4
macerrorsize = macerrormax - macerrormin
errormin     = 3.3333e-8
errormax     = 6.6666e-2
errorsize    = errormax - errormin
timemin      = 3.3333e-3
timemax      = 0.8 #2.0
timesize     = timemax - timemin
sim_no       = 0
thetamac     = 0.6
particle_distribution = "random"

xerror    = errormin*math.pow(10.0, 0.025*math.log10(errormax/errormin))
yerror    = errormin*math.pow(10.0, 0.925*math.log10(errormax/errormin))
xtheta    = thetamin*math.pow(10.0, 0.025*math.log10(thetamax/thetamin))
xmacerror = macerrormin*math.pow(10.0, 0.025*math.log10(macerrormax/macerrormin))
ytime     = timemin*math.pow(10.0, 0.925*math.log10(timemax/timemin))

theta_values                           = []
monopole_geo_error_values              = []
monopole_geo_poterror_values           = []
monopole_geo_runtimes                  = []
fast_monopole_geo_error_values         = []
fast_monopole_geo_poterror_values      = []
fast_monopole_geo_runtimes             = []
quadrupole_geo_error_values            = []
quadrupole_geo_poterror_values         = []
quadrupole_geo_runtimes                = []
fast_quadrupole_geo_error_values       = []
fast_quadrupole_geo_poterror_values    = []
fast_quadrupole_geo_runtimes           = []

macerror_values                        = []
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

quadrupole_eigen_error_values          = []
quadrupole_eigen_poterror_values       = []
quadrupole_eigen_runtimes              = []
fast_quadrupole_eigen_error_values     = []
fast_quadrupole_eigen_poterror_values  = []
fast_quadrupole_eigen_runtimes         = []



# Run brute-force test to get accurate forces
#--------------------------------------------------------------------------------------------------
bruteforce_sim = newsim("freefall.dat")
bruteforce_sim.SetParam("particle_distribution",particle_distribution)
bruteforce_sim.SetParam('neib_search','bruteforce')
bruteforce_sim.SetParam('Nstepsmax',Nstepsmax)
bruteforce_sim.SetParam('Nhydro',Nhydro)
bruteforce_sim.SetParam('run_id','FREEFALL1-BF')
setupsim()
ax_bf     = get_data("ax")
ay_bf     = get_data("ay")
az_bf     = get_data("az")
gpot_bf   = get_data("gpot")

start = time.time() #bruteforce_sim.GetBlockTime("SPH_ALL_FORCES")
run()
end = time.time() #bruteforce_sim.GetBlockTime("SPH_ALL_FORCES")
dt_direct = end - start


# GEOMETRIC MAC
# Cell monopole errors as a function of opening angle
#--------------------------------------------------------------------------------------------------
for i in range(numSimMax):
    sim_no = sim_no + 1
    exponent = 1.0*i/(numSimMax - 1)
    theta = thetamin*math.pow(thetamax/thetamin,exponent)
    thetamaxsqd = theta*theta
    sim = newsim("freefall.dat")
    sim.SetParam("particle_distribution",particle_distribution)
    sim.SetParam('gravity_mac',"geometric")
    sim.SetParam('thetamaxsqd',thetamaxsqd)
    sim.SetParam('multipole',"fast_monopole")
    sim.SetParam('Nleafmax',Nleafmax)
    sim.SetParam('Nstepsmax',Nstepsmax)
    sim.SetParam('Nhydro',Nhydro)
    setupsim()
    ax_tree   = get_data("ax")
    ay_tree   = get_data("ay")
    az_tree   = get_data("az")
    gpot_tree = get_data("gpot")
    N         = ax_tree.size
    error     = ForceError(N, ax_bf, ay_bf, az_bf, ax_tree, ay_tree, az_tree)
    poterror  = PotentialError(N, gpot_bf, gpot_tree)
    start = time.time() #sim.GetBlockTime("SPH_ALL_FORCES")
    run()
    end = time.time() #sim.GetBlockTime("SPH_ALL_FORCES")
    dt = end - start
    theta_values.append(theta)
    fast_monopole_geo_error_values.append(error)
    fast_monopole_geo_poterror_values.append(poterror)
    fast_monopole_geo_runtimes.append(dt/dt_direct) #(end - start)


# Monopole errors as a function of opening angle
#--------------------------------------------------------------------------------------------------
for i in range(numSimMax):
    sim_no = sim_no + 1
    exponent = 1.0*i/(numSimMax - 1)
    theta = thetamin*math.pow(thetamax/thetamin,exponent)
    thetamaxsqd = theta*theta
    sim = newsim("freefall.dat")
    sim.SetParam("particle_distribution",particle_distribution)
    sim.SetParam('gravity_mac',"geometric")
    sim.SetParam('thetamaxsqd',thetamaxsqd)
    sim.SetParam('multipole',"monopole")
    sim.SetParam('Nleafmax',Nleafmax)
    sim.SetParam('Nstepsmax',Nstepsmax)
    sim.SetParam('Nhydro',Nhydro)
    setupsim()
    ax_tree   = get_data("ax")
    ay_tree   = get_data("ay")
    az_tree   = get_data("az")
    gpot_tree = get_data("gpot")
    N         = ax_tree.size
    error     = ForceError(N, ax_bf, ay_bf, az_bf, ax_tree, ay_tree, az_tree)
    poterror  = PotentialError(N, gpot_bf, gpot_tree)
    start = time.time() #sim.GetBlockTime("SPH_ALL_FORCES")
    run()
    end = time.time() #sim.GetBlockTime("SPH_ALL_FORCES")
    dt = end - start
    monopole_geo_error_values.append(error)
    monopole_geo_poterror_values.append(poterror)
    monopole_geo_runtimes.append(dt/dt_direct) #(end - start)


# Quadrupole errors as a function of opening angle
#--------------------------------------------------------------------------------------------------
for i in range(numSimMax):
    sim_no = sim_no + 1
    exponent = 1.0*i/(numSimMax - 1)
    theta = thetamin*math.pow(thetamax/thetamin,exponent)
    thetamaxsqd = theta*theta
    sim = newsim("freefall.dat")
    sim.SetParam("particle_distribution",particle_distribution)
    sim.SetParam('gravity_mac',"geometric")
    sim.SetParam('thetamaxsqd',thetamaxsqd)
    sim.SetParam('multipole',"quadrupole")
    sim.SetParam('Nleafmax',Nleafmax)
    sim.SetParam('Nstepsmax',Nstepsmax)
    sim.SetParam('Nhydro',Nhydro)
    setupsim()
    ax_tree   = get_data("ax")
    ay_tree   = get_data("ay")
    az_tree   = get_data("az")
    gpot_tree = get_data("gpot")
    N         = ax_tree.size
    error     = ForceError(N, ax_bf, ay_bf, az_bf, ax_tree, ay_tree, az_tree)
    poterror  = PotentialError(N, gpot_bf, gpot_tree)
    start = time.time() #sim.GetBlockTime("SPH_ALL_FORCES")
    run()
    end = time.time() #sim.GetBlockTime("SPH_ALL_FORCES")
    dt = end - start
    quadrupole_geo_error_values.append(error)
    quadrupole_geo_poterror_values.append(poterror)
    quadrupole_geo_runtimes.append(dt/dt_direct) #(end - start)


# Cell quadrupole errors as a function of opening angle
#--------------------------------------------------------------------------------------------------
for i in range(numSimMax):
    sim_no = sim_no + 1
    exponent = 1.0*i/(numSimMax - 1)
    theta = thetamin*math.pow(thetamax/thetamin,exponent)
    thetamaxsqd = theta*theta
    sim = newsim("freefall.dat")
    sim.SetParam("particle_distribution",particle_distribution)
    sim.SetParam('gravity_mac',"geometric")
    sim.SetParam('thetamaxsqd',thetamaxsqd)
    sim.SetParam('multipole',"fast_quadrupole")
    sim.SetParam('Nleafmax',Nleafmax)
    sim.SetParam('Nstepsmax',Nstepsmax)
    sim.SetParam('Nhydro',Nhydro)
    setupsim()
    ax_tree   = get_data("ax")
    ay_tree   = get_data("ay")
    az_tree   = get_data("az")
    gpot_tree = get_data("gpot")
    N         = ax_tree.size
    error     = ForceError(N, ax_bf, ay_bf, az_bf, ax_tree, ay_tree, az_tree)
    poterror  = PotentialError(N, gpot_bf, gpot_tree)
    start = time.time() #sim.GetBlockTime("SPH_ALL_FORCES")
    run()
    end = time.time() #sim.GetBlockTime("SPH_ALL_FORCES")
    dt = end - start
    fast_quadrupole_geo_error_values.append(error)
    fast_quadrupole_geo_poterror_values.append(poterror)
    fast_quadrupole_geo_runtimes.append(dt/dt_direct)  #(end - start)


# GADGET2 MAC
# Cell monopole errors as a function of opening angle
#--------------------------------------------------------------------------------------------------
for i in range(numSimMax):
    sim_no = sim_no + 1
    exponent = 1.0*i/(numSimMax - 1)
    macerror = macerrormin*math.pow(macerrormax/macerrormin, exponent)
    sim = newsim("freefall.dat")
    sim.SetParam("particle_distribution",particle_distribution)
    sim.SetParam('gravity_mac',"gadget2")
    sim.SetParam('macerror',macerror)
    sim.SetParam('thetamaxsqd',thetamac*thetamac)
    sim.SetParam('multipole',"fast_monopole")
    sim.SetParam('Nleafmax',Nleafmax)
    sim.SetParam('Nstepsmax',Nstepsmax)
    sim.SetParam('Nhydro',Nhydro)
    setupsim()
    ax_tree   = get_data("ax")
    ay_tree   = get_data("ay")
    az_tree   = get_data("az")
    gpot_tree = get_data("gpot")
    N         = ax_tree.size
    error     = ForceError(N, ax_bf, ay_bf, az_bf, ax_tree, ay_tree, az_tree)
    poterror  = PotentialError(N, gpot_bf, gpot_tree)
    macerror_values.append(macerror)
    start = time.time() #sim.GetBlockTime("SPH_ALL_FORCES")
    run()
    end = time.time() #sim.GetBlockTime("SPH_ALL_FORCES")
    dt = end - start
    fast_monopole_gadget_error_values.append(error)
    fast_monopole_gadget_poterror_values.append(poterror)
    fast_monopole_gadget_runtimes.append(dt/dt_direct)  #(end - start)


# Monopole errors as a function of opening angle
#--------------------------------------------------------------------------------------------------
for i in range(numSimMax):
    sim_no = sim_no + 1
    exponent = 1.0*i/(numSimMax - 1)
    macerror = macerrormin*math.pow(macerrormax/macerrormin, exponent)
    sim = newsim("freefall.dat")
    sim.SetParam("particle_distribution",particle_distribution)
    sim.SetParam('gravity_mac',"gadget2")
    sim.SetParam('macerror',macerror)
    sim.SetParam('thetamaxsqd',thetamac*thetamac)
    sim.SetParam('multipole',"monopole")
    sim.SetParam('Nleafmax',Nleafmax)
    sim.SetParam('Nstepsmax',Nstepsmax)
    sim.SetParam('Nhydro',Nhydro)
    setupsim()
    ax_tree   = get_data("ax")
    ay_tree   = get_data("ay")
    az_tree   = get_data("az")
    gpot_tree = get_data("gpot")
    N         = ax_tree.size
    error     = ForceError(N, ax_bf, ay_bf, az_bf, ax_tree, ay_tree, az_tree)
    poterror  = PotentialError(N, gpot_bf, gpot_tree)
    start = time.time() #sim.GetBlockTime("SPH_ALL_FORCES")
    run()
    end = time.time() #sim.GetBlockTime("SPH_ALL_FORCES")
    dt = end - start
    monopole_gadget_error_values.append(error)
    monopole_gadget_poterror_values.append(poterror)
    monopole_gadget_runtimes.append(dt/dt_direct)  #(end - start)


# Quadrupole errors as a function of opening angle
#--------------------------------------------------------------------------------------------------
for i in range(numSimMax):
    sim_no = sim_no + 1
    exponent = 1.0*i/(numSimMax - 1)
    macerror = macerrormin*math.pow(macerrormax/macerrormin, exponent)
    sim = newsim("freefall.dat")
    sim.SetParam("particle_distribution",particle_distribution)
    sim.SetParam('gravity_mac',"gadget2")
    sim.SetParam('macerror',macerror)
    sim.SetParam('thetamaxsqd',thetamac*thetamac)
    sim.SetParam('multipole',"quadrupole")
    sim.SetParam('Nleafmax',Nleafmax)
    sim.SetParam('Nstepsmax',Nstepsmax)
    sim.SetParam('Nhydro',Nhydro)
    setupsim()
    ax_tree   = get_data("ax")
    ay_tree   = get_data("ay")
    az_tree   = get_data("az")
    gpot_tree = get_data("gpot")
    N         = ax_tree.size
    error     = ForceError(N, ax_bf, ay_bf, az_bf, ax_tree, ay_tree, az_tree)
    poterror  = PotentialError(N, gpot_bf, gpot_tree)
    start = time.time() #sim.GetBlockTime("SPH_ALL_FORCES")
    run()
    end = time.time() #sim.GetBlockTime("SPH_ALL_FORCES")
    dt = end - start
    quadrupole_gadget_error_values.append(error)
    quadrupole_gadget_poterror_values.append(poterror)
    quadrupole_gadget_runtimes.append(dt/dt_direct)  #(end - start)


# Cell quadrupole errors as a function of opening angle
#--------------------------------------------------------------------------------------------------
for i in range(numSimMax):
    sim_no = sim_no + 1
    exponent = 1.0*i/(numSimMax - 1)
    macerror = macerrormin*math.pow(macerrormax/macerrormin, exponent)
    sim = newsim("freefall.dat")
    sim.SetParam("particle_distribution",particle_distribution)
    sim.SetParam('gravity_mac',"gadget2")
    sim.SetParam('macerror',macerror)
    sim.SetParam('thetamaxsqd',thetamac*thetamac)
    sim.SetParam('multipole',"fast_quadrupole")
    sim.SetParam('Nleafmax',Nleafmax)
    sim.SetParam('Nstepsmax',Nstepsmax)
    sim.SetParam('Nhydro',Nhydro)
    setupsim()
    ax_tree   = get_data("ax")
    ay_tree   = get_data("ay")
    az_tree   = get_data("az")
    gpot_tree = get_data("gpot")
    N         = ax_tree.size
    error     = ForceError(N, ax_bf, ay_bf, az_bf, ax_tree, ay_tree, az_tree)
    poterror  = PotentialError(N, gpot_bf, gpot_tree)
    start = time.time() #sim.GetBlockTime("SPH_ALL_FORCES")
    run()
    end = time.time() #sim.GetBlockTime("SPH_ALL_FORCES")
    dt = end - start
    fast_quadrupole_gadget_error_values.append(error)
    fast_quadrupole_gadget_poterror_values.append(poterror)
    fast_quadrupole_gadget_runtimes.append(dt/dt_direct)  #(end - start)




# EIGENVALUE MAC
# Quadrupole errors as a function of opening angle
#--------------------------------------------------------------------------------------------------
for i in range(numSimMax):
    sim_no = sim_no + 1
    exponent = 1.0*i/(numSimMax - 1)
    macerror = macerrormin*math.pow(macerrormax/macerrormin, exponent)
    sim = newsim("freefall.dat")
    sim.SetParam("particle_distribution",particle_distribution)
    sim.SetParam('gravity_mac',"eigenmac")
    sim.SetParam('macerror',macerror)
    sim.SetParam('thetamaxsqd',thetamac*thetamac)
    sim.SetParam('multipole',"quadrupole")
    sim.SetParam('Nleafmax',Nleafmax)
    sim.SetParam('Nstepsmax',Nstepsmax)
    sim.SetParam('Nhydro',Nhydro)
    setupsim()
    ax_tree   = get_data("ax")
    ay_tree   = get_data("ay")
    az_tree   = get_data("az")
    gpot_tree = get_data("gpot")
    N         = ax_tree.size
    error     = ForceError(N, ax_bf, ay_bf, az_bf, ax_tree, ay_tree, az_tree)
    poterror  = PotentialError(N, gpot_bf, gpot_tree)
    start = time.time() #sim.GetBlockTime("SPH_ALL_FORCES")
    run()
    end = time.time() #sim.GetBlockTime("SPH_ALL_FORCES")
    dt = end - start
    quadrupole_eigen_error_values.append(error)
    quadrupole_eigen_poterror_values.append(poterror)
    quadrupole_eigen_runtimes.append(dt/dt_direct)  #(end - start)


# Cell quadrupole errors as a function of opening angle
#--------------------------------------------------------------------------------------------------
for i in range(numSimMax):
    sim_no = sim_no + 1
    exponent = 1.0*i/(numSimMax - 1)
    macerror = macerrormin*math.pow(macerrormax/macerrormin, exponent)
    sim = newsim("freefall.dat")
    sim.SetParam("particle_distribution",particle_distribution)
    sim.SetParam('gravity_mac',"eigenmac")
    sim.SetParam('macerror',macerror)
    sim.SetParam('thetamaxsqd',thetamac*thetamac)
    sim.SetParam('multipole',"fast_quadrupole")
    sim.SetParam('Nleafmax',Nleafmax)
    sim.SetParam('Nstepsmax',Nstepsmax)
    sim.SetParam('Nhydro',Nhydro)
    setupsim()
    ax_tree   = get_data("ax")
    ay_tree   = get_data("ay")
    az_tree   = get_data("az")
    gpot_tree = get_data("gpot")
    N         = ax_tree.size
    error     = ForceError(N, ax_bf, ay_bf, az_bf, ax_tree, ay_tree, az_tree)
    poterror  = PotentialError(N, gpot_bf, gpot_tree)
    start = time.time() #sim.GetBlockTime("SPH_ALL_FORCES")
    run()
    end = time.time() #sim.GetBlockTime("SPH_ALL_FORCES")
    dt = end - start
    fast_quadrupole_eigen_error_values.append(error)
    fast_quadrupole_eigen_poterror_values.append(poterror)
    fast_quadrupole_eigen_runtimes.append(dt/dt_direct)  #(end - start)



print "THETA     : ",theta_values
print "MONOERROR : ",monopole_geo_error_values
print "QUADERROR : ",quadrupole_geo_error_values
print "MONOTIMES : ",monopole_geo_runtimes
print "QUADTIMES : ",quadrupole_geo_runtimes

# Create matplotlib figure object with shared x-axis
#--------------------------------------------------------------------------------------------------
fig, axarr = plt.subplots(1, 3, sharey='col', figsize=(14,5))
fig.subplots_adjust(hspace=0.001, wspace=0.0001)
fig.subplots_adjust(bottom=0.1, top=0.98, left=0.06, right=0.99)

axarr[0].set_xscale("log")
axarr[0].set_yscale("log")
axarr[0].set_ylabel(r"$t_{_{\rm CPU}}$")
axarr[0].set_xlabel(r"$\delta\,{\bf a}$")
axarr[0].set_xlim([errormin, errormax])
axarr[0].set_ylim([timemin, timemax])
axarr[0].plot(monopole_geo_error_values, monopole_geo_runtimes, color="red", linestyle='-', label='Monopole')
axarr[0].plot(fast_monopole_geo_error_values, fast_monopole_geo_runtimes, color="blue", linestyle=':', label='Cell monopole')
axarr[0].plot(quadrupole_geo_error_values, quadrupole_geo_runtimes, color="black", linestyle='--', label='Quadrupole')
axarr[0].plot(fast_quadrupole_geo_error_values, fast_quadrupole_geo_runtimes, color="green", linestyle='-.', label='Cell quadrupole')
axarr[0].scatter(monopole_geo_error_values, monopole_geo_runtimes, color="red")
axarr[0].scatter(fast_monopole_geo_error_values, fast_monopole_geo_runtimes, color="blue")
axarr[0].scatter(quadrupole_geo_error_values, quadrupole_geo_runtimes, color="black")
axarr[0].scatter(fast_quadrupole_geo_error_values, fast_quadrupole_geo_runtimes, color="green")
axarr[0].text(xerror, ytime, '(a) Geometric MAC', color='black', size=14)
legend0 = axarr[0].legend(loc='lower left', fontsize=12)

axarr[1].set_xscale("log")
axarr[1].set_yscale("log")
#axarr[1].set_ylabel(r"$t_{_{\rm CPU}}$")
axarr[1].set_xlabel(r"$\delta\,{\bf a}$")
axarr[1].set_xlim([errormin, errormax])
axarr[1].set_ylim([timemin, timemax])
axarr[1].plot(monopole_gadget_error_values, monopole_gadget_runtimes, color="red", linestyle='-', label='Monopole')
axarr[1].plot(fast_monopole_gadget_error_values, fast_monopole_gadget_runtimes, color="blue", linestyle=':', label='Cell monopole')
axarr[1].plot(quadrupole_gadget_error_values, quadrupole_gadget_runtimes, color="black", linestyle='--', label='Quadrupole')
axarr[1].plot(fast_quadrupole_gadget_error_values, fast_quadrupole_gadget_runtimes, color="green", linestyle='-.', label='Cell quadrupole')
axarr[1].scatter(monopole_gadget_error_values, monopole_gadget_runtimes, color="red")
axarr[1].scatter(fast_monopole_gadget_error_values, fast_monopole_gadget_runtimes, color="blue")
axarr[1].scatter(quadrupole_gadget_error_values, quadrupole_gadget_runtimes, color="black")
axarr[1].scatter(fast_quadrupole_gadget_error_values, fast_quadrupole_gadget_runtimes, color="green")
axarr[1].text(xerror, ytime, '(b) GADGET-MAC', color='black', size=14)
axarr[1].set_yticks([])
#legend1 = axarr[1].legend(loc='upper right', fontsize=12)

axarr[2].set_xscale("log")
axarr[2].set_yscale("log")
#axarr[2].set_ylabel(r"$t_{_{\rm CPU}}$")
axarr[2].set_xlabel(r"$\delta\,{\bf a}$")
axarr[2].set_xlim([errormin, errormax])
axarr[2].set_ylim([timemin, timemax])
axarr[2].plot(quadrupole_eigen_error_values, quadrupole_eigen_runtimes, color="black", linestyle='--', label='Quadrupole')
axarr[2].plot(fast_quadrupole_eigen_error_values, fast_quadrupole_eigen_runtimes, color="green", linestyle='-.', label='Cell quadrupole')
axarr[2].scatter(quadrupole_eigen_error_values, quadrupole_eigen_runtimes, color="black")
axarr[2].scatter(fast_quadrupole_eigen_error_values, fast_quadrupole_eigen_runtimes, color="green")
axarr[2].text(xerror, ytime, '(c) Eigen-MAC', color='black', size=14)
axarr[2].set_yticks([])
#legend2 = axarr[2].legend(loc='upper right', fontsize=12)

plt.show()
fig.savefig('treeerrortime.pdf', dpi=50)



# Create matplotlib figure object with shared x-axis
#--------------------------------------------------------------------------------------------------
fig, axarr = plt.subplots(1, 3, sharey='col', figsize=(14,5))
fig.subplots_adjust(hspace=0.001, wspace=0.0001)
fig.subplots_adjust(bottom=0.1, top=0.98, left=0.06, right=0.99)


axarr[0].set_xscale("log")
axarr[0].set_yscale("log")
axarr[0].set_ylabel(r"$\delta\,{\bf a}$")
axarr[0].set_xlabel(r"$\theta$")
axarr[0].set_xlim([thetamin, thetamax])
axarr[0].set_ylim([errormin, errormax])
axarr[0].plot(theta_values, monopole_geo_error_values, color="red", linestyle='-', label='Monopole')
axarr[0].plot(theta_values, fast_monopole_geo_error_values, color="blue", linestyle=':', label='Cell monopole')
axarr[0].plot(theta_values, quadrupole_geo_error_values, color="black", linestyle='--', label='Quadrupole')
axarr[0].plot(theta_values, fast_quadrupole_geo_error_values, color="green", linestyle='-.', label='Cell quadrupole')
axarr[0].scatter(theta_values, monopole_geo_error_values, color="red")
axarr[0].scatter(theta_values, fast_monopole_geo_error_values, color="blue")
axarr[0].scatter(theta_values, quadrupole_geo_error_values, color="black")
axarr[0].scatter(theta_values, fast_quadrupole_geo_error_values, color="green")
axarr[0].text(xtheta, yerror, '(a) Geometric MAC', color='black', size=14)
legend0 = axarr[0].legend(loc='lower right', fontsize=16)

axarr[1].set_xscale("log")
axarr[1].set_yscale("log")
#axarr[1].set_ylabel(r"$\delta\,\phi$")
axarr[1].set_xlabel(r"$\alpha_{_{\rm MAC}}$")
axarr[1].set_xlim([macerrormin, macerrormax])
axarr[1].set_ylim([errormin, errormax])
axarr[1].plot(macerror_values, monopole_gadget_error_values, color="red", linestyle='-', label='Monopole')
axarr[1].plot(macerror_values, fast_monopole_gadget_error_values, color="blue", linestyle=':', label='Cell monopole')
axarr[1].plot(macerror_values, quadrupole_gadget_error_values, color="black", linestyle='--', label='Quadrupole')
axarr[1].plot(macerror_values, fast_quadrupole_gadget_error_values, color="green", linestyle='-.', label='Cell quadrupole')
axarr[1].scatter(macerror_values, monopole_gadget_error_values, color="red")
axarr[1].scatter(macerror_values, fast_monopole_gadget_error_values, color="blue")
axarr[1].scatter(macerror_values, quadrupole_gadget_error_values, color="black")
axarr[1].scatter(macerror_values, fast_quadrupole_gadget_error_values, color="green")
axarr[1].text(xmacerror, yerror, '(b) GADGET-MAC', color='black', size=16)
axarr[1].set_yticks([])
#legend1 = axarr[1].legend(loc='lower right', fontsize=12)

axarr[2].set_xscale("log")
axarr[2].set_yscale("log")
#axarr[2].set_ylabel(r"$\delta\,\phi$")
axarr[2].set_xlabel(r"$\alpha_{_{\rm MAC}}$")
axarr[2].set_xlim([macerrormin, macerrormax])
axarr[2].set_ylim([errormin, errormax])
axarr[2].plot(macerror_values, quadrupole_eigen_error_values, color="black", linestyle='--', label='Quadrupole')
axarr[2].plot(macerror_values, fast_quadrupole_eigen_error_values, color="green", linestyle='-.', label='Cell quadrupole')
axarr[2].scatter(macerror_values, quadrupole_eigen_error_values, color="black")
axarr[2].scatter(macerror_values, fast_quadrupole_eigen_error_values, color="green")
axarr[2].text(xmacerror, yerror, '(c) Eigen-MAC', color='black', size=16)
axarr[2].set_yticks([])
#legend1 = axarr[2].legend(loc='lower right', fontsize=12)

plt.show()
fig.savefig('treeerror.pdf', dpi=50)

# Prevent program from closing before showing plot window
block()
