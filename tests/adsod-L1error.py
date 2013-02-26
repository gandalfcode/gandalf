from seren.analysis.facade import *
import time


Nres = []
L1values = []


sim1 = newsim("adsod1024.dat")
run()
Nres.append(1024+128)
L1values.append(L1errornorm("x","rho",-2.5,2.5))
plot("x","rho")
plotanalytical("x","rho")


sim2 = newsim("adsod512.dat")
run()
Nres.append(512+64)
L1values.append(L1errornorm("x","rho",-2.5,2.5))
addplot("x","rho")


sim3 = newsim("adsod256.dat")
run()
Nres.append(256+32)
L1values.append(L1errornorm("x","rho",-2.5,2.5))
addplot("x","rho")


sim4 = newsim("adsod128.dat")
run()
Nres.append(128+16)
L1values.append(L1errornorm("x","rho",-2.5,2.5))
addplot("x","rho")


print Nres
print L1values

block()



