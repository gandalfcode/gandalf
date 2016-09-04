#==============================================================================
# freefalltest.py
# Run the freefall collapse test using initial conditions specified in the
# file 'freefall.dat'.
#==============================================================================
from gandalf.analysis.facade import *
from gandalf.analysis.data_fetcher import *
from gandalf.analysis.compute import lagrangian_radii
import matplotlib.pyplot as plt
import time

def timeratiofreefall(snap,type=None,unit="default"):
    t_ff = np.pi/2*np.sqrt(0.5)
    return snap.t/t_ff
    sim=snap.sim
#     time_unit_obj = sim.simunits.t
#     if unit=="default":
#         unit = time_unit_obj.outunit
#     unitinfo = UnitInfo()
#     unitinfo.name = unit
#     unitinfo.label = time_unit_obj.LatexLabel(unit)
#     scaling_factor = time_unit_obj.OutputScale(unit)
#     label = 't'
#     return unitinfo, snap.t, scaling_factor, label



# Create new freefall collapse simulation from 'freefall.dat' file and then
# run the simulation
newsim('freefall.dat')
setupsim()
run()

# Create Lagrangian radii data for 10%, 50% and 90% mass radii.
CreateTimeData('lr1',lagrangian_radii,mfrac=0.1)
CreateTimeData('lr2',lagrangian_radii,mfrac=0.5)
CreateTimeData('lr3',lagrangian_radii,mfrac=0.9,label='r/R$_0$')
CreateTimeData('t',timeratiofreefall)



# Plot Lagrangian radii as a function of time
time_plot("t","lr3",linestyle='-')
limit("lr3",0.0,1.05)
time_plot("t","lr2",overplot=True,linestyle='-')
time_plot("t","lr1",overplot=True,linestyle='-')
plt.gca().set_ylabel('r/R$_0$')
plt.gca().set_xlabel('t/t$_\mathrm{ff}$')
# plt.gca().set_xscale('log')
# plt.gca().set_yscale('log')
# plt.gca().set_ylim(0.01,1)
# plt.gca().set_xlim(1,0.001)
plt.draw()

# Prevent program from closing before showing plot window
block()
