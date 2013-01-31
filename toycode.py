from numpy import *
from pylab import *
import copy
#from SimBuffer import *
import matplotlib.pyplot as plt
import SphSim




# ============================================================================
# Main program
# ============================================================================


# Create simulation 1
sim1 = SphSim.SphSimulation()

sim1.Setup()
#buf.add_simulation(sim1)

snap1 = SphSim.SphSnapshot()

snap1.CopyDataFromSimulation(2,sim1.sph.Nsph,sim1.sph.sphdata)
x1 = copy.deepcopy(snap1.ExtractArray("x"))
y1 = copy.deepcopy(snap1.ExtractArray("rho"))
ion()
plt.scatter(x1,y1)
plt.show()


sim1.Run(2000)


snap1.CopyDataFromSimulation(2,sim1.sph.Nsph,sim1.sph.sphdata)
x1 = copy.deepcopy(snap1.ExtractArray("x"))
y1 = copy.deepcopy(snap1.ExtractArray("rho"))
plt.scatter(x1,y1,color="green")
ioff()
plt.show()


exit()
