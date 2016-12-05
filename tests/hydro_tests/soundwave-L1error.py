#==============================================================================
# soundwave-L1error.py
# Run the soundwave test at various resolutions calculating the L1 error norm
# for each case.  For a second-order scheme (as SPH is in smooth-flow), 
#==============================================================================
from gandalf.analysis.facade import *
from gandalf.analysis.compute import L1errornorm
import time


# Create empty arrays to contain results (resolution vs. L1 error norms)
resolutions = [64,128]
L1values = []

def run_and_plot_soundwave(res):
    sim = newsim("soundwave.dat")
    sim.SetParam("Nhydro",res)
    setupsim()
    run()
    errnorm=L1errornorm("soundwave","x","rho",0.01,0.99)
    addplot("x","rho")
    return errnorm
  
for res in resolutions:
    L1values.append(run_and_plot_soundwave(res))
plotanalytical("x","rho",overplot=True)

# For now, print out values to screen
print resolutions
print L1values
block()


