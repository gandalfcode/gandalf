#==============================================================================
# jeanstest.py
# Run Noh test using initial conditions file 'noh.dat'.
#==============================================================================
from gandalf.analysis.facade import *
import time
import matplotlib.pyplot as plt
import numpy as np
import math
from matplotlib import rc
from mpl_toolkits.axes_grid1 import AxesGrid
from scipy.optimize import curve_fit




#--------------------------------------------------------------------------------------------------
rc('font', **{'family': 'normal', 'weight' : 'bold', 'size' : 16})
rc('text', usetex=True)


# Set all plot limits
xmin       = 0.0
xmax       = 1.0
rhomin     = 0.95
rhomax     = 1.05
stride     = 7
sim_no     = 0
#tend       = 0.05
amp        = 0.025
rho0       = 1.0
press0     = 1.0
wavelength = 1.0
mu_bar     = 1.0
lmin       = 0.05*wavelength
lmax       = 1.95*wavelength
ratiomin   = 0.0
ratiomax   = 2.0
numSimMax  = 24

stable_lambdas   = []
stable_periods   = []
unstable_lambdas = []
unstable_periods = []

stable_wavelengths   = np.arange(lmin, 0.999*wavelength, 0.001)
unstable_wavelengths = np.arange(1.001*wavelength, lmax, 0.001)

x_solution = np.arange(0.0, wavelength, 0.001*wavelength)


# Analytical solutions for stable and unstable modes (for function fitting)
#--------------------------------------------------------------------------------------------------
#def jeans_stable_solution(x, omega, rhofit):
#    return rhofit*(1.0 + amp*np.sin(2.0*math.pi*x/wavelength)*np.cos(omega*tend))

#def jeans_unstable_solution(x, omega, rhofit):
#    return rho0*(1.0 + amp*np.sin(2.0*math.pi*x/wavelength)*np.cosh(omega*tend))

def oscillation_period(lambda_jeans):
    return np.sqrt(np.pi/rho0)*wavelength/np.sqrt(lambda_jeans*lambda_jeans - wavelength*wavelength)

def growth_timescale(lambda_jeans):
    return np.sqrt(1/4.0/np.pi/rho0)*wavelength/np.sqrt(wavelength*wavelength - lambda_jeans*lambda_jeans)



# Perform small suite of simulations to compare to expected with dispersion relation.
#--------------------------------------------------------------------------------------------------
for i in range(numSimMax):

    # Chose jeans length
    ratio  = ratiomin + (i + 0.5)*(ratiomax - ratiomin)/numSimMax
    if ratio > 0.99 and ratio < 1.01: continue

    lambda_jeans = wavelength/ratio
    csound = lambda_jeans*math.sqrt(rho0/math.pi)
    temp0  = mu_bar*csound*csound
    if wavelength > lambda_jeans:
        tsim = growth_timescale(lambda_jeans)
    else:
        tsim = 0.3*oscillation_period(lambda_jeans)

    print 'LAMBDA_JEANS : ',lambda_jeans,tsim

    sim = newsim('jeans.dat')
    sim.SetParam('tend', tsim)
    sim.SetParam('tsnapfirst', tsim)
    sim.SetParam('temp0', temp0)
    setupsim()
    run()

    # Get grav. data for simulaton for plotting
    x_data   = get_data("x", sim=sim_no)
    vx_data  = get_data("vx", sim=sim_no)
    rho_data = get_data("rho", sim=sim_no)
    rho_analytical = get_analytical_data("x", "rho", sim=sim_no)
    vx_analytical  = get_analytical_data("x", "vx", sim=sim_no)
    sim_no = sim_no + 1

    # Compute best-fit for simulation
    if wavelength > lambda_jeans:
        omega_analytical = 1.0/growth_timescale(lambda_jeans)
        print "OMEGA_ANA : ",omega_analytical
        def jeans_unstable_solution(x, omega, rhofit):
            return rho0*(1.0 + amp*np.sin(2.0*math.pi*x/wavelength)*np.cosh(omega*tsim))
        popt, pcov = curve_fit(jeans_unstable_solution, x_data, rho_data,
                               bounds=([0.5*omega_analytical, 0.9*rho0], [2.0*omega_analytical, 1.1*rho0]))
        unstable_lambdas.append(wavelength/lambda_jeans)
        unstable_periods.append(1.0/popt[0])

        print 'UNSTABLE SOLUTION : ',popt[0],popt[1]

        fig, axarr = plt.subplots(1, 1, figsize=(7,10), sharex='row')
        axarr.set_ylabel(r"$\rho$")
        axarr.set_xlim([xmin, xmax])
        axarr.plot(rho_analytical.x_data, rho_analytical.y_data, color="red", linestyle='-', lw=0.5)
        axarr.plot(x_solution, jeans_unstable_solution(x_solution, *popt), color='blue')
        axarr.scatter(x_data[::stride], rho_data[::stride], color='black', marker='.', s=4.0)
        #plt.show()
        #block()

    else:
        omega_analytical = 2.0*math.pi/oscillation_period(lambda_jeans)
        print "OMEGA_ANA : ",omega_analytical
        def jeans_stable_solution(x, omega, rhofit):
            return rhofit*(1.0 + amp*np.sin(2.0*math.pi*x/wavelength)*np.cos(omega*tsim))
        popt, pcov = curve_fit(jeans_stable_solution, x_data, rho_data,
                               bounds=([0.7*omega_analytical, 0.9*rho0], [1.5*omega_analytical, 1.1*rho0]))
        stable_lambdas.append(wavelength/lambda_jeans)
        stable_periods.append(2.0*math.pi/popt[0])

        print 'STABLE SOLUTION : ',popt[0],popt[1]

        fig, axarr = plt.subplots(1, 1, figsize=(7,10), sharex='row')
        axarr.set_ylabel(r"$\rho$")
        axarr.set_xlim([xmin, xmax])
        axarr.plot(rho_analytical.x_data, rho_analytical.y_data, color="red", linestyle='-', lw=0.5)
        axarr.plot(x_solution, jeans_stable_solution(x_solution, *popt), color='blue')
        axarr.scatter(x_data[::stride], rho_data[::stride], color='black', marker='.', s=4.0)
        #plt.show()
        #block()


# Create matplotlib figure object with shared x-axis
#--------------------------------------------------------------------------------------------------
fig, axarr = plt.subplots(2, 1, figsize=(7,10), sharex='row')
fig.subplots_adjust(hspace=0.001, wspace=0.001)
fig.subplots_adjust(bottom=0.08, top=0.97, left=0.13, right=0.98)

axarr[0].set_ylabel(r"$\rho$")
axarr[0].set_xlim([xmin, xmax])
axarr[0].plot(rho_analytical.x_data, rho_analytical.y_data, color="red", linestyle='-', lw=0.5)
#axarr[0].plot(x_solution, jeans_unstable_solution(x_solution, *popt), color='blue')
axarr[0].scatter(x_data[::stride], rho_data[::stride], color='black', marker='.', s=4.0)
axarr[1].set_ylabel(r"$v_{x}$")
axarr[1].set_xlabel(r"$x$")
axarr[1].set_xlim([xmin, xmax])
axarr[1].plot(vx_analytical.x_data, vx_analytical.y_data, color="red", linestyle='-', lw=0.5)
axarr[1].scatter(x_data[::stride], vx_data[::stride], color='black', marker='.', s=4.0)

plt.show()
fig.savefig('jeanstest.eps', dpi=50)



# Dispersion relation figure
#--------------------------------------------------------------------------------------------------
fig2, axarr2 = plt.subplots(1, 1, figsize=(7,5), sharex='row')
fig2.subplots_adjust(hspace=0.001, wspace=0.001)
fig2.subplots_adjust(bottom=0.11, top=0.97, left=0.1, right=0.98)


x_stable_aux = wavelength/stable_wavelengths
y_stable = oscillation_period(x_stable_aux)

x_unstable_aux = wavelength/unstable_wavelengths
y_unstable = growth_timescale(x_unstable_aux)


axarr2.set_xlabel(r"$\lambda/\lambda_{_{\rm JEANS}}$")
axarr2.set_ylabel(r"$T$")
#axarr2.set_yscale("log")
axarr2.set_xlim([lmin, lmax])
axarr2.set_ylim(0.0, 3.0)
#axarr[0].plot(rho_analytical.x_data, rho_analytical.y_data, color="red", linestyle='-', lw=0.5)
axarr2.plot(stable_wavelengths, y_stable, color='red')
axarr2.plot(unstable_wavelengths, y_unstable, color='red')
axarr2.plot([1.0, 1.0], [0.0, 3.0], color='blue', linestyle='-.', linewidth=0.5)
#axarr2.scatter(ratio_values, omega_values, color='black', marker='*', s=8.0)
axarr2.scatter(stable_lambdas, stable_periods, marker='o', facecolors='none', edgecolors='black', s=24.0, label='Oscillating')
axarr2.scatter(unstable_lambdas, unstable_periods, facecolors='none', edgecolors='black', marker='*', s=24.0, label='Collapsing')
legend = axarr2.legend(loc='upper right', fontsize=12)

plt.show()
fig2.savefig('jeansdispersion.eps', dpi=50)


run()
block()
