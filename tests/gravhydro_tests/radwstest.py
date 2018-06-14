#==============================================================================
# jeanstest.py
# Run Noh test using initial conditions file 'noh.dat'.
#==============================================================================
from gandalf.analysis.facade import *
import time

# Create new cloud collapse simulation object
newsim('radws_test.dat')
setupsim()
run()



# Now compute the central tempertature of the last snapshot
tmax = snap(-1).t
rho = get_data('rho')
U   = get_data('u')

args = rho.argsort()[-10:]
rho_c = rho[args].mean()
U_c   = U[args].mean()

print 'Time:', tmax
print 'Central density / Specific energy:', rho_c, U_c
