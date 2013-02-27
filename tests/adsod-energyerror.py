from seren.analysis.facade import *
import time


dtmult = []
energyvalues = []


sim0 = newsimfromparams("adsod.dat")
sim0.simparams.floatparams["courant_mult"] = 0.01
sim0.simparams.floatparams["accel_mult"] = 0.01
sim0.SetupSimulation()
run()
dtmult.append(0.01)
energyvalues.append(sim0.diag.Eerror)

sim1 = newsimfromparams("adsod.dat")
sim1.simparams.floatparams["courant_mult"] = 0.02
sim1.simparams.floatparams["accel_mult"] = 0.02
sim1.SetupSimulation()
run()
dtmult.append(0.02)
energyvalues.append(sim1.diag.Eerror)

sim2 = newsimfromparams("adsod.dat")
sim2.simparams.floatparams["courant_mult"] = 0.04
sim2.simparams.floatparams["accel_mult"] = 0.04
sim2.SetupSimulation()
run()
dtmult.append(0.04)
energyvalues.append(sim2.diag.Eerror)

sim3 = newsimfromparams("adsod.dat")
sim3.simparams.floatparams["courant_mult"] = 0.08
sim3.simparams.floatparams["accel_mult"] = 0.08
sim3.SetupSimulation()
run()
dtmult.append(0.08)
energyvalues.append(sim3.diag.Eerror)

sim4 = newsimfromparams("adsod.dat")
sim4.simparams.floatparams["courant_mult"] = 0.16
sim4.simparams.floatparams["accel_mult"] = 0.16
sim4.SetupSimulation()
run()
dtmult.append(0.16)
energyvalues.append(sim4.diag.Eerror)

sim5 = newsimfromparams("adsod.dat")
sim5.simparams.floatparams["courant_mult"] = 0.32
sim5.simparams.floatparams["accel_mult"] = 0.32
sim5.SetupSimulation()
run()
dtmult.append(0.32)
energyvalues.append(sim5.diag.Eerror)


print dtmult
print energyvalues

block()


