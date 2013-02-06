import numpy as np
import shocktub

class AnalyticalSolution():
    def __init__(self):
        pass
    
    def compute(self, parameters):
        pass
    
class shocktube (AnalyticalSolution):
    def __init__(self, sim, time):
        AnalyticalSolution.__init__(self)
        self.time = time
        #extract the parameters
        simfloatparams = sim.simparams.floatparams
        self.RHOinL = simfloatparams["rhofluid1"]
        self.RHOinR = simfloatparams["rhofluid2"]
        self.UinL = simfloatparams["vfluid1[0]"]
        self.UinR = simfloatparams["vfluid2[0]"]
        self.PinL = simfloatparams["press1"]
        self.PinR = simfloatparams["press2"]
        self.xL = simfloatparams["boxmin[0]"]
        self.x0 = 0.
        self.xR = simfloatparams["boxmax[0]"]
        self.time = time
        self.iMAX = 1000
        if sim.simparams.stringparams["gas_eos"] == "isothermal":
            self.gamma = 1+1e-5
        else:
            self.gamma = simfloatparams["gamma_eos"]
        
    def compute(self, x, y):
        shocktub.shocktub(self.RHOinL, self.RHOinR, self.UinL, self.UinR,
                 self.PinL, self.PinR, self.xL, self.x0, self.xR,
                 self.time, self.iMAX, self.gamma)
        data=np.genfromtxt('sod.out',names=['x','rho','vx','pressure'])
        return data[x], data[y]