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
    
    
class soundwave (AnalyticalSolution):
    def __init__(self, sim, time):
        AnalyticalSolution.__init__(self)
        self.time = time
        #extract the parameters
        simfloatparams = sim.simparams.floatparams
        self.rho = simfloatparams["rhofluid1"]
        self.press = simfloatparams["press1"]
        self.temp0 = simfloatparams["temp0"]
        self.mu_bar = simfloatparams["mu_bar"]
        self.amp = simfloatparams["amp"]
        self.xL = simfloatparams["boxmin[0]"]
        self.xR = simfloatparams["boxmax[0]"]
        if sim.simparams.stringparams["gas_eos"] == "isothermal":
            self.gamma = 1.0
            self.csound = np.sqrt(self.temp0/self.mu_bar)
        else:
            self.gamma = simfloatparams["gamma_eos"]
            self.csound = np.sqrt(self.gamma*self.press/self.rho)
        self.wlambda = self.xR - self.xL
        self.kwave = 2.0*3.1415/self.wlambda
        self.omega = 2.0*3.1415*self.csound/self.wlambda
        self.time = time
        self.iMAX = 1000
        if sim.simparams.stringparams["gas_eos"] == "isothermal":
            self.gamma = 1+1e-5
        else:
            self.gamma = simfloatparams["gamma_eos"]
        
    def compute(self, ix, iy):
        x = np.arange(self.xL,self.xR,1.0/self.iMAX)
        rho = self.rho*(1.0 + self.amp*np.sin(self.kwave*x - self.omega*self.time))
        vx = self.csound*self.amp*np.sin(self.kwave*x - self.omega*self.time)
        if ix == "x" and iy == "rho": return x,rho
        if ix == "x" and iy == "vx": return x,vx