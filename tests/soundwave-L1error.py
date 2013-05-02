#==============================================================================
# soundwave-L1error.py
# Run the soundwave test at various resolutions calculating the L1 error norm
# for each case.  For a second-order scheme (as SPH is in smooth-flow), 
#==============================================================================
from seren.analysis.facade import *
from seren.analysis.compute import L1errornorm
import time


# Create empty arrays to contain results (resolution vs. L1 error norms)
Nres = []
L1values = []


# Simulation 1, N = 8
# -----------------------------------------------
sim1 = newsimfromparams("soundwave.dat")
sim1.simparams.intparams["Npart"] = 8
sim1.SetupSimulation()
run()
Nres.append(8)
L1values.append(L1errornorm("x","rho",0.01,0.99))
plot("x","rho")
plotanalytical("x","rho")


# Simulation 2, N = 16
# -----------------------------------------------
sim2 = newsimfromparams("soundwave.dat")
sim2.simparams.intparams["Npart"] = 16
sim2.SetupSimulation()
run()
Nres.append(16)
L1values.append(L1errornorm("x","rho",0.01,0.99))
addplot("x","rho")


# Simulation 3, N = 32
# -----------------------------------------------
sim3 = newsimfromparams("soundwave.dat")
sim3.simparams.intparams["Npart"] = 32
sim3.SetupSimulation()
run()
Nres.append(32)
L1values.append(L1errornorm("x","rho",0.01,0.99))
addplot("x","rho")


# Simulation 4, N = 64
# -----------------------------------------------
sim4 = newsimfromparams("soundwave.dat")
sim4.simparams.intparams["Npart"] = 64
sim4.SetupSimulation()
run()
Nres.append(64)
L1values.append(L1errornorm("x","rho",0.01,0.99))
addplot("x","rho")


# Simulation 5, N = 128
# -----------------------------------------------
sim5 = newsimfromparams("soundwave.dat")
sim5.simparams.intparams["Npart"] = 128
sim5.SetupSimulation()
run()
Nres.append(128)
L1values.append(L1errornorm("x","rho",0.01,0.99))
addplot("x","rho")


# Simulation 6, N = 256
# -----------------------------------------------
sim6 = newsimfromparams("soundwave.dat")
sim6.simparams.intparams["Npart"] = 256
sim6.SetupSimulation()
run()
Nres.append(256)
L1values.append(L1errornorm("x","rho",0.01,0.99))
addplot("x","rho")


# For now, print out values to screen
print Nres
print L1values
block()


