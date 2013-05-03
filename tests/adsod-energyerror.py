#==============================================================================
# adsod-energyerror.py
# Run the adiabatic Sod shocktube test for various values of the timestep
# multiplier to determine the energy error scaling relative to the timestep.
#==============================================================================
from seren.analysis.facade import *
import time


dtmult = []
energyvalues = []


# Simulation 1, tmult = 0.01
# -----------------------------------------------
sim1 = newsimfromparams("adsod.dat")
sim1.simparams.floatparams["courant_mult"] = 0.01
sim1.simparams.floatparams["accel_mult"] = 0.01
sim1.SetupSimulation()
run()
dtmult.append(0.01)
energyvalues.append(sim1.diag.Eerror)


# Simulation 1, tmult = 0.02
# -----------------------------------------------
sim2 = newsimfromparams("adsod.dat")
sim2.simparams.floatparams["courant_mult"] = 0.02
sim2.simparams.floatparams["accel_mult"] = 0.02
sim2.SetupSimulation()
run()
dtmult.append(0.02)
energyvalues.append(sim2.diag.Eerror)


# Simulation 1, tmult = 0.04
# -----------------------------------------------
sim3 = newsimfromparams("adsod.dat")
sim3.simparams.floatparams["courant_mult"] = 0.04
sim3.simparams.floatparams["accel_mult"] = 0.04
sim3.SetupSimulation()
run()
dtmult.append(0.04)
energyvalues.append(sim3.diag.Eerror)


# Simulation 1, tmult = 0.08
# -----------------------------------------------
sim4 = newsimfromparams("adsod.dat")
sim4.simparams.floatparams["courant_mult"] = 0.08
sim4.simparams.floatparams["accel_mult"] = 0.08
sim4.SetupSimulation()
run()
dtmult.append(0.08)
energyvalues.append(sim4.diag.Eerror)


# Simulation 1, tmult = 0.16
# -----------------------------------------------
sim5 = newsimfromparams("adsod.dat")
sim5.simparams.floatparams["courant_mult"] = 0.16
sim5.simparams.floatparams["accel_mult"] = 0.16
sim5.SetupSimulation()
run()
dtmult.append(0.16)
energyvalues.append(sim5.diag.Eerror)


# Simulation 1, tmult = 0.32
# -----------------------------------------------
sim6 = newsimfromparams("adsod.dat")
sim6.simparams.floatparams["courant_mult"] = 0.32
sim6.simparams.floatparams["accel_mult"] = 0.32
sim6.SetupSimulation()
run()
dtmult.append(0.32)
energyvalues.append(sim6.diag.Eerror)


# For now, print out values to screen
print dtmult
print energyvalues
block()
