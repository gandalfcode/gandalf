from gandalf.analysis.facade import *
import numpy as np
from matplotlib import rcParams
import matplotlib.pyplot as plt
rcParams['font.size'] = 16
print rcParams.keys()
rcParams['font.family'] = 'serif'


loadsim('ADSOD1')
snap(-2)

x   = get_data('x')
rho = get_data('rho')
vx  = get_data('vx')

Ngas = len(x) / 2


plt.plot(x[:Ngas], rho[:Ngas], 'r+', label='gas')
plt.plot(x[Ngas:], rho[Ngas:], 'ko', label='dust')
plt.xlabel('$x$')
plt.ylabel(r'$\rho$')
#plt.xlim(-1,1)
plt.legend(loc='best')
plt.figure()
plt.plot(x[:Ngas], vx[:Ngas], 'r+', label='gas')
plt.plot(x[Ngas:], vx[Ngas:], 'ko', label='dust')
#plt.xlim(-1,1)
#plt.ylim(-0.1,0.7)
plt.xlabel('$x$')
plt.ylabel('$v_x$')
plt.legend(loc='best')

plt.show()
