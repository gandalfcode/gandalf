#==============================================================================
# adsod-L1error.py
# Computes L1 error norms for various resolutions of the adiabatic Sod test.
# Used to demonstrate scaling of error with resolution.
#==============================================================================
from gandalf.analysis.facade import *
from gandalf.analysis.compute import L1errornorm
import time


# Set empty lists to store results from each resolution
Nres = []
L1values = []


# Lowest resolution simulation (LHS 128 particles, RHS 16 particles)
sim1 = newsim("adsod.dat")
sim1.SetParam('Nlattice1[0]',128)
sim1.SetParam('Nlattice2[0]',16)
setupsim()
run()
Nres.append(128+16)
L1values.append(L1errornorm("shocktube","x","rho",-2.0,2.0))
plot("x","rho")
plotanalytical("x","rho")


# LHS 256 particles, RHS 32 particles
sim2 = newsim("adsod.dat")
sim2.SetParam('Nlattice1[0]',256)
sim2.SetParam('Nlattice2[0]',32)
setupsim()
run()
Nres.append(256+32)
L1values.append(L1errornorm("shocktube","x","rho",-2.0,2.0))
addplot("x","rho")


# LHS 512 particles, RHS 64 particles
sim3 = newsim("adsod.dat")
sim3.SetParam('Nlattice1[0]',512)
sim3.SetParam('Nlattice2[0]',64)
setupsim()
run()
Nres.append(512+64)
L1values.append(L1errornorm("shocktube","x","rho",-2.0,2.0))
addplot("x","rho")


# LHS 1024 particles, RHS 128 particles
sim4 = newsim("adsod.dat")
sim4.SetParam('Nlattice1[0]',1024)
sim4.SetParam('Nlattice2[0]',128)
setupsim()
run()
Nres.append(1024+128)
L1values.append(L1errornorm("shocktube","x","rho",-2.0,2.0))
addplot("x","rho")


# LHS 2048 particles, RHS 256 particles
sim5 = newsim("adsod.dat")
sim5.SetParam('Nlattice1[0]',2048)
sim5.SetParam('Nlattice2[0]',256)
setupsim()
run()
Nres.append(2048+256)
L1values.append(L1errornorm("shocktube","x","rho",-2.0,2.0))
addplot("x","rho")


print Nres
print L1values

block()



