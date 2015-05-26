#==============================================================================
# example03.py
# Example to create a `blank' simulation object, set all important parameters,
# then run the simulation to completion.
#==============================================================================
from gandalf.analysis.facade import *

# Create `blank' simulation object (2-dimensional SPH)
sim = newsim(ndim=2,sim="sph")

# Set all important simulation parameters in order to create a grad-h SPH
# simulation of a 2D Sedov blast test initially from a 64x64 lattice.
sim.SetParam("ic","sedov")
sim.SetParam("run_id","SEDOV1")
sim.SetParam("Nlattice1[0]",64)
sim.SetParam("Nlattice1[1]",64)
sim.SetParam("boxmin[0]",-1.0)
sim.SetParam("boxmin[1]",-1.0)
sim.SetParam("boxmax[0]",1.0)
sim.SetParam("boxmax[1]",1.0)
sim.SetParam("boundary_lhs[0]","periodic")
sim.SetParam("boundary_rhs[0]","periodic")
sim.SetParam("boundary_lhs[1]","periodic")
sim.SetParam("boundary_rhs[1]","periodic")
sim.SetParam("dimensionless",1)
sim.SetParam("Nlevels",10)
sim.SetParam("tend",0.5)
sim.SetParam("tsnapfirst",0.0)
sim.SetParam("dt_snap",0.1)

# Now perform all set-up routines and run simulation
setupsim()
run()
