from seren.analysis.facade import *
from seren.analysis.compute import L1errornorm
import time


Nres = []
L1values = []


sim0 = newsimfromparams("soundwave.dat")
sim0.simparams.intparams["Npart"] = 8
sim0.SetupSimulation()
run()
Nres.append(16)
L1values.append(L1errornorm("x","rho",0.01,0.99))
plot("x","rho")
plotanalytical("x","rho")



sim1 = newsimfromparams("soundwave.dat")
sim1.simparams.intparams["Npart"] = 16
sim1.SetupSimulation()
run()
Nres.append(16)
L1values.append(L1errornorm("x","rho",0.01,0.99))
addplot("x","rho")



sim2 = newsimfromparams("soundwave.dat")
sim2.simparams.intparams["Npart"] = 32
sim2.SetupSimulation()
run()
Nres.append(32)
L1values.append(L1errornorm("x","rho",0.01,0.99))
addplot("x","rho")



sim3 = newsimfromparams("soundwave.dat")
sim3.simparams.intparams["Npart"] = 64
sim3.SetupSimulation()
run()
Nres.append(64)
L1values.append(L1errornorm("x","rho",0.01,0.99))
addplot("x","rho")



sim4 = newsimfromparams("soundwave.dat")
sim4.simparams.intparams["Npart"] = 128
sim4.SetupSimulation()
run()
Nres.append(128)
L1values.append(L1errornorm("x","rho",0.01,0.99))
addplot("x","rho")



sim5 = newsimfromparams("soundwave.dat")
sim5.simparams.intparams["Npart"] = 256
sim5.SetupSimulation()
run()
Nres.append(128)
L1values.append(L1errornorm("x","rho",0.01,0.99))
addplot("x","rho")



print Nres
print L1values

block()


