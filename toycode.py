from numpy import *
from pylab import *
import copy
from SimBuffer import *
import matplotlib.pyplot as plt
import SphSim
import SphSnap




# ============================================================================
# Main program
# ============================================================================


# Create main simulation buffer object
buf = SimBuffer()

# Define and initialise global variables
state = 0
debug = 0
script = 'none'


# Create simulation 1
sim1 = SphSim.SphSimulation()

sim1.Setup()
buf.add_simulation(sim1)

snap1 = SphSnap.SphSnapshot()

snap1.CopyDataFromSimulation(3,sim1.sph.Nsph,sim1.sph.sphdata)
x1 = copy.deepcopy(snap1.ExtractArray("x"))
y1 = copy.deepcopy(snap1.ExtractArray("rho"))
ion()
plt.x
plt.scatter(x1,y1)
plt.show()


for i in range(1000):
    sim1.MainLoop()


snap1.CopyDataFromSimulation(3,sim1.sph.Nsph,sim1.sph.sphdata)
x1 = copy.deepcopy(snap1.ExtractArray("x"))
y1 = copy.deepcopy(snap1.ExtractArray("rho"))
plt.scatter(x1,y1,color="green")
ioff()
plt.show()


exit()
