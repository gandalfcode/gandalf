#==============================================================================
# freefalltest.py
# Run the freefall collapse test using initial conditions specified in the
# file 'freefall.dat'.
#==============================================================================
from gandalf.analysis.facade import *
from gandalf.analysis.data_fetcher import *
from gandalf.analysis.compute import lagrangian_radii
import matplotlib.pyplot as plt
from scipy.optimize import brentq
from math import pi,sqrt,acos
import time

def r_inverted(r,t):
    r05=sqrt(r)
    return 2./pi*(acos(r05)+r05*sqrt(1-r))-t

def freefall_analytical_radius(t):
    return brentq(r_inverted, 0, 1,args=t)

def timeratiofreefall(snap,type=None,unit="default"):
    t_ff = np.pi/2*np.sqrt(0.5)
    return snap.t/t_ff
    sim=snap.sim



if __name__=="__main__":
    # Create new freefall collapse simulation from 'freefall.dat' file and then
    # run the simulation
    import sys
    if len(sys.argv) > 1:
        newsim(sys.argv[1])
    else:
        newsim('freefall.dat')
    setupsim()
    run()
    
    # Create Lagrangian radii data for 10%, 50% and 90% mass radii.
    CreateTimeData('lr1',lagrangian_radii,mfrac=0.1)
    r_5=CreateTimeData('lr2',lagrangian_radii,mfrac=0.5).fetch()[1]
    CreateTimeData('lr3',lagrangian_radii,mfrac=0.9,label='r/R$_0$')
    time=CreateTimeData('tr',timeratiofreefall).fetch()[1]
    initial_radius=r_5[0]
    r_5/=r_5[0]
    
    analytical_5 = np.empty_like(r_5)
    for i,t in enumerate(time):
        analytical_5[i]=freefall_analytical_radius(t)
    
    print 'Error norm:', 
    print np.linalg.norm((analytical_5 - r_5)*initial_radius, ord=1)/time.size

    # Plot Lagrangian radii as a function of time
    time_plot("tr","lr3",linestyle='-')
    limit("lr3",0.0,1.05)
    time_plot("tr","lr2",overplot=True,linestyle='-')
    time_plot("tr","lr1",overplot=True,linestyle='-')
    plt.plot(time,analytical_5*initial_radius)
    plt.gca().set_ylabel('r/R$_0$')
    plt.gca().set_xlabel('t/t$_\mathrm{ff}$')
    # plt.gca().set_xscale('log')
    # plt.gca().set_yscale('log')
    # plt.gca().set_ylim(0.01,1)
    # plt.gca().set_xlim(1,0.001)
    plt.draw()
    
    # Prevent program from closing before showing plot window
    block()
