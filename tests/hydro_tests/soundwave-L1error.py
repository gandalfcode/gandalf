#==============================================================================
# soundwave-L1error.py
# Run the soundwave test at various resolutions calculating the L1 error norm
# for each case.  For a second-order scheme (as SPH is in smooth-flow), 
#==============================================================================
from gandalf.analysis.facade import *
from gandalf.analysis.compute import L1errornorm
import time


# Create empty arrays to contain results (resolution vs. L1 error norms)
Nres = []
L1values = []


# Simulation 1, N = 8
# -----------------------------------------------
sim1 = newsim("soundwave.dat")
sim1.SetParam("Nhydro",8)
setupsim()
run()
Nres.append(8)
L1values.append(L1errornorm("soundwave","x","rho",0.01,0.99))
plot("x","rho")
plotanalytical("x","rho")


# Simulation 2, N = 16
# -----------------------------------------------
sim2 = newsim("soundwave.dat")
sim2.SetParam("Nhydro",16)
setupsim()
run()
Nres.append(16)
L1values.append(L1errornorm("soundwave","x","rho",0.01,0.99))
addplot("x","rho")


# Simulation 3, N = 32
# -----------------------------------------------
sim3 = newsim("soundwave.dat")
sim3.SetParam("Nhydro",32)
setupsim()
run()
Nres.append(32)
L1values.append(L1errornorm("soundwave","x","rho",0.01,0.99))
addplot("x","rho")


# Simulation 4, N = 64
# -----------------------------------------------
sim4 = newsim("soundwave.dat")
sim4.SetParam("Nhydro",64)
setupsim()
run()
Nres.append(64)
L1values.append(L1errornorm("soundwave","x","rho",0.01,0.99))
addplot("x","rho")


# Simulation 5, N = 128
# -----------------------------------------------
sim5 = newsim("soundwave.dat")
sim5.SetParam("Nhydro",128)
setupsim()
run()
Nres.append(128)
L1values.append(L1errornorm("soundwave","x","rho",0.01,0.99))
addplot("x","rho")


# Simulation 6, N = 256
# -----------------------------------------------
sim6 = newsim("soundwave.dat")
sim6.SetParam("Nhydro",256)
sim6.SetupSimulation()
run()
Nres.append(256)
L1values.append(L1errornorm("soundwave","x","rho",0.01,0.99))
addplot("x","rho")


# For now, print out values to screen
print Nres
print L1values
block()


