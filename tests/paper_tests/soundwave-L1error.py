#==============================================================================
# soundwave-L1error.py
# Computes L1 error norms for various resolutions of the adiabatic Sod test.
# Used to demonstrate scaling of error with resolution.
#==============================================================================
from gandalf.analysis.facade import *
from gandalf.analysis.compute import L1errornorm
import time
import math


# Set empty lists to store results from each resolution
gradh_Nres = []
gradh_L1values = []


# Lowest resolution simulation (16 particles)
for i in range(5):
    Nhydro = int(8*math.pow(2,i))
    sim = newsim("soundwave.dat")
    sim.SetParam('Nhydro',Nhydro)
    setupsim()
    run()
    gradh_Nres.append(Nhydro)
    gradh_L1values.append(L1errornorm("soundwave","x","rho",0.0,1.0))
    plot("x","rho")
    plotanalytical("x","rho")



print gradh_Nres
print gradh_L1values


block()
