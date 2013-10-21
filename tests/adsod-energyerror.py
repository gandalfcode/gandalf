#==============================================================================
# adsod-energyerror.py
# Run the adiabatic Sod shocktube test for various values of the timestep
# multiplier to determine the energy error scaling relative to the timestep.
#==============================================================================
from gandalf.analysis.facade import *
import time


dtmult = []
energyvalues = []


# Simulation 1, tmult = 0.01
# -----------------------------------------------
sim1 = newsim("adsod.dat")
sim1.SetParam("courant_mult",0.01)
sim1.SetParam("accel_mult",0.01)
setupsim()
run()
dtmult.append(0.01)
energyvalues.append(sim1.diag.Eerror)


# Simulation 1, tmult = 0.02
# -----------------------------------------------
sim2 = newsim("adsod.dat")
sim2.SetParam("courant_mult",0.02)
sim2.SetParam("accel_mult",0.02)
setupsim()
run()
dtmult.append(0.02)
energyvalues.append(sim2.diag.Eerror)


# Simulation 1, tmult = 0.04
# -----------------------------------------------
sim3 = newsim("adsod.dat")
sim3.SetParam("courant_mult",0.04)
sim3.SetParam("accel_mult",0.04)
setupsim()
run()
dtmult.append(0.04)
energyvalues.append(sim3.diag.Eerror)


# Simulation 1, tmult = 0.08
# -----------------------------------------------
sim4 = newsim("adsod.dat")
sim4.SetParam("courant_mult",0.08)
sim4.SetParam("accel_mult",0.08)
setupsim()
run()
dtmult.append(0.08)
energyvalues.append(sim4.diag.Eerror)


# Simulation 1, tmult = 0.16
# -----------------------------------------------
sim5 = newsim("adsod.dat")
sim5.SetParam("courant_mult",0.16)
sim5.SetParam("accel_mult",0.16)
setupsim()
run()
dtmult.append(0.16)
energyvalues.append(sim5.diag.Eerror)


# Simulation 1, tmult = 0.32
# -----------------------------------------------
sim6 = newsim("adsod.dat")
sim6.SetParam("courant_mult",0.32)
sim6.SetParam("accel_mult",0.32)
setupsim()
run()
dtmult.append(0.32)
energyvalues.append(sim6.diag.Eerror)


# For now, print out values to screen
print dtmult
print energyvalues
block()
