#==============================================================================
# analytical.py
# ..
#==============================================================================
import numpy as np
import shocktub

'''This module contains the classes responsible for computing the analytical solution, in the case we know it.
There is a template empty class, which shows that the class must expose a compute function.
'''


#------------------------------------------------------------------------------
class AnalyticalSolution():
    '''
    Template for an analytical solution. Exposes a compute function,
    which is the one called by the function in facade. Takes two strings,
    x and y, which are the names of the two quantities that get returned
    in a tuple.
    All the parameters should be passed to __init__, so that compute
    doesn\'t need any other information to do its job (in this way, the
    object can be passed around after creation if needed).
    '''
    def __init__(self):
        pass
    
    def compute(self, x, y):
        pass


#------------------------------------------------------------------------------
class freefall (AnalyticalSolution):
    '''Analytical solution for the Noh problem test.
    When instantiated, it gets passed the sim object and
    the time for which the solution is requested. The values
    of the relevant variables are pulled from the simulation
    object and saved inside the object.
    '''
    def __init__(self, sim, time):
        AnalyticalSolution.__init__(self)
        self.time = time
        #extract the parameters
        simfloatparams = sim.simparams.floatparams
        self.radius = simfloatparams["radius"]
        self.rho = simfloatparams["rhofluid1"]
        self.time = time
        self.ndim = sim.ndims
        self.iMAX = 1000
        
    def compute(self, x, y):
        '''Computes the exact solution of the Noh problem for
        1, 2 and 3 dimensions'''
        r = np.arange(0.0,self.radius,1.0/self.iMAX)
        rho = self.rho*np.ones(self.iMAX)
        if self.time > 0:
            bound = 0.3333333333333*self.time
            i = 0
            while i < self.iMAX:
                if self.ndim == 1:
                    if x[i] < bound: rho[i] = 4.0*self.rho
                    else: rho[i] = self.rho
                elif self.ndim == 2:
                    if x[i] < bound: rho[i] = 16.0*self.rho
                    else: rho[i] = self.rho*(1.0 + self.time/x[i])
                elif self.ndim == 1:
                    if x[i] < bound: rho[i] = 64.0*self.rho
                    else: rho[i] = self.rho*pow(1.0 + self.time/x[i],2.0)
                i = i + 1
        return x,rho


#------------------------------------------------------------------------------
class noh (AnalyticalSolution):
    '''Analytical solution for the Noh problem test.
    When instantiated, it gets passed the sim object and
    the time for which the solution is requested. The values
    of the relevant variables are pulled from the simulation
    object and saved inside the object.
    '''
    def __init__(self, sim, time):
        AnalyticalSolution.__init__(self)
        self.time = time
        #extract the parameters
        simfloatparams = sim.simparams.floatparams
        self.radius = simfloatparams["radius"]
        self.rho = simfloatparams["rhofluid1"]
        self.time = time
        self.ndim = sim.ndims
        self.iMAX = 1000
        self.gamma = simfloatparams["gamma_eos"]
        
    def compute(self, x, y):
        '''Computes the exact solution of the Noh problem for
        1, 2 and 3 dimensions'''
        x = np.arange(0.0,self.radius,1.0/self.iMAX)
        rho = self.rho*np.ones(self.iMAX)
        if self.time > 0:
            bound = 0.3333333333333*self.time
            i = 0
            while i < self.iMAX:
                if self.ndim == 1:
                    if x[i] < bound: rho[i] = 4.0*self.rho
                    else: rho[i] = self.rho
                elif self.ndim == 2:
                    if x[i] < bound: rho[i] = 16.0*self.rho
                    else: rho[i] = self.rho*(1.0 + self.time/x[i])
                elif self.ndim == 1:
                    if x[i] < bound: rho[i] = 64.0*self.rho
                    else: rho[i] = self.rho*pow(1.0 + self.time/x[i],2.0)
                i = i + 1
        return x,rho


#------------------------------------------------------------------------------
class shocktube (AnalyticalSolution):
    '''Analytical solution for the shocktube test.
    When instantiated, it gets passed the sim object and
    the time for which the solution is requested. The values
    of the relevant variables are pulled from the simulation
    object and saved inside the object.
    '''
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
        '''Computes the exact solution of the Riemann problem.
        Gets passed two strings with the quantities that are needed
        (e.g., 'rho' and 'pressure').'''
        #calls the fortran module that computes the state
        shocktub.shocktub(self.RHOinL, self.RHOinR, self.UinL, self.UinR,
                 self.PinL, self.PinR, self.xL, self.x0, self.xR,
                 self.time, self.iMAX, self.gamma)
        #reads the data from the text file produced
        data=np.genfromtxt('sod.out',names=['x','rho','vx','pressure','u'])
        return data[x], data[y]
    

#------------------------------------------------------------------------------
class soundwave (AnalyticalSolution):
    '''Analytical solution for the soundwave.
    The initialization works in the same way as the shocktube.
    '''
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
        self.kwave = 2.0*np.pi/self.wlambda
        self.omega = 2.0*np.pi*self.csound/self.wlambda
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
        ax = -self.csound*self.csound*self.kwave*self.rho*self.amp*np.cos(self.kwave*x - self.omega*self.time)
        if ix == "x" and iy == "rho": 
            return x,rho
        elif ix == "x" and iy == "vx": 
            return x,vx
        elif ix == "x" and iy == "ax": 
            return x,ax
        else:
            raise KeyError("There were errors in the quantity you requested")
