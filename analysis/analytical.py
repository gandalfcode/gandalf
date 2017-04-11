#==============================================================================
#  analytical.py
#  Contains all python classes for computing analytical solutions for key
#  test problems for plotting or analysis.
#
#  This file is part of GANDALF :
#  Graphical Astrophysics code for N-body Dynamics And Lagrangian Fluids
#  https://github.com/gandalfcode/gandalf
#  Contact : gandalfcode@gmail.com
#
#  Copyright (C) 2013  D. A. Hubber, G. Rosotti
#
#  GANDALF is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 2 of the License, or
#  (at your option) any later version.
#
#  GANDALF is distributed in the hope that it will be useful, but
#  WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  General Public License (http://www.gnu.org/licenses) for more details.
#==============================================================================
import numpy as np
from swig_generated.SphSim import ShocktubeSolution
from scipy.interpolate import interp1d
#from scipy.special import gamma as GammaF

'''This module contains the classes responsible for computing the analytical solution, in the case we know it.
There is a template empty class, which shows that the class must expose a compute function.
'''


#------------------------------------------------------------------------------
class AnalyticalSolution():
    '''Template for an analytical solution. Exposes a compute function, which
is the one called by the function in facade. Takes two strings, x and y, which
are the names of the two quantities that get returned in a tuple.  All the
parameters should be passed to __init__, so that compute doesn\'t need any
other information to do its job (in this way, the object can be passed around
after creation if needed).
    '''
    def __init__(self):
        pass

    def compute(self, x, y):
        pass


#------------------------------------------------------------------------------
class freefall (AnalyticalSolution):
    '''Analytical solution for the freefall collapse test.'''

    def __init__(self, sim, time):
        AnalyticalSolution.__init__(self)
        self.time = time

        # Extract the parameters
        simfloatparams = sim.simparams.floatparams
        self.radius = simfloatparams["radius"]
        self.rho    = simfloatparams["rhofluid1"]
        self.time   = time
        self.ndim   = sim.ndims
        self.iMAX   = 1000
        self.tff    = np.sqrt(3.0*3.1415/32.0/self.rho)
        self.mtot   = 4.0*3.14157*self.rho*math.pow(self.radius,3)/3

    # TODO : This is not the freefall solution.  Need to update!!
    def compute(self, x, y):
        '''Computes exact solution for freefall collapse problem'''
        r = np.zeros(self.iMAX)
        a = np.zeros(self.iMAX)
        gpot = np.zeros(self.iMAX)

        i = 0
        while i < self.iMAX:
            r[i] = self.radius*i/(self.iMAX - 1)
            a[i] = r[i]
            gpot[i] = 0.0;
            i = i + 1

        if ix == "t" and iy == "r":
            return t,r
        else:
            raise KeyError("There were errors in the quantity you requested")



#------------------------------------------------------------------------------
class gresho (AnalyticalSolution):
    '''Analytical solution for the Gresho-Chan vortex problem test.'''

    def __init__(self, sim, time):
        AnalyticalSolution.__init__(self)
        self.time = time

        # Extract the parameters
        simfloatparams = sim.simparams.floatparams
        self.radius = simfloatparams["boxmax[0]"] - simfloatparams["boxmin[0]"]
        self.ndim   = sim.ndims
        self.iMAX   = 1000
        self.gamma  = simfloatparams["gamma_eos"]
        self.rho    = 1.0

    def compute(self, ix, iy):
        '''Computes the exact solution of the Gresho vortex problem.'''
        R = np.arange(0.0,self.radius,1.0/self.iMAX)
        press = self.rho*np.ones(self.iMAX)
        vphi = np.zeros(self.iMAX)

        i = 0
        while i < self.iMAX:
            if R[i] < 0.2:
                vphi[i] = 5.0*R[i]
                press[i] = 5.0 + 12.5*R[i]*R[i]
            elif R[i] < 0.4:
                vphi[i] = 2.0 - 5.0*R[i]
                press[i] = 9.0 + 12.5*R[i]*R[i] - 20.0*R[i] + 4.0*np.log10(R[i]/0.2)
            else:
                vphi[i] = 0.0
                press[i] = 3.0 + 4.0*np.log10(2.0)
            i = i + 1

        if ix == "R" and iy == "press":
            return R,press
        elif ix == "R" and iy == "vphi":
            return R,vphi
        else:
            raise KeyError("There were errors in the quantity you requested")


#------------------------------------------------------------------------------
class jeans (AnalyticalSolution):
    '''Analytical solution for the Jeans instability test with 1D sinusoidal perturbation.'''

    def __init__(self, sim, time):
        AnalyticalSolution.__init__(self)
        self.time = time

        # Extract the parameters
        simfloatparams = sim.simparams.floatparams
        self.rho     = simfloatparams["rhofluid1"]
        self.press   = simfloatparams["press1"]
        self.temp0   = simfloatparams["temp0"]
        self.mu_bar  = simfloatparams["mu_bar"]
        self.time    = time
        self.ndim    = sim.ndims
        self.iMAX    = 1000
        self.amp     = simfloatparams["amp"]
        self.xL      = simfloatparams["boxmin[0]"]
        self.xR      = simfloatparams["boxmax[0]"]
        self.wlambda = self.xR - self.xL
        self.kwave   = 2.0*3.14159/self.wlambda

        if sim.simparams.stringparams["sim"] == "nbody":
            self.gamma  = 1.0
            self.csound = 0.0
            self.omega  = np.sqrt(4.0*3.14159*self.rho)
            self.wlambdaJeans = 0.0

        else:

            if sim.simparams.stringparams["gas_eos"] == "isothermal":
                self.gamma = 1.0
                self.csound = np.sqrt(self.temp0/self.mu_bar)
            else:
                self.gamma = simfloatparams["gamma_eos"]
                self.csound = np.sqrt(self.gamma*self.press/self.rho)

            self.wlambdaJeans = np.sqrt(3.14159*self.csound*self.csound/self.rho)
            if self.wlambda < self.wlambdaJeans:
                self.omega = 2.0*3.14159*self.csound*np.sqrt(1.0/self.wlambda/self.wlambda - 1.0/self.wlambdaJeans/self.wlambdaJeans)
            elif self.wlambdaJeans < self.wlambda:
                self.omega = 2.0*3.14159*self.csound*np.sqrt(1.0/self.wlambdaJeans/self.wlambdaJeans - 1.0/self.wlambda/self.wlambda)
            else:
                self.omega = 0.0

    def compute(self, ix, iy):
        x = np.arange(self.xL,self.xR,1.0/self.iMAX)
        if self.wlambda < self.wlambdaJeans:
            rho = self.rho*(1.0 + self.amp*np.sin(self.kwave*x)*np.cos(self.omega*self.time))
            vx = -self.amp*self.omega*np.cos(self.kwave*x)*np.sin(self.omega*self.time)/self.kwave
            ax = -self.amp*self.omega*self.omega*np.cos(self.kwave*x)*np.cos(self.omega*self.time)/self.kwave
        elif self.wlambdaJeans < self.wlambda:
            rho = self.rho*(1.0 + self.amp*np.sin(self.kwave*x)*np.cosh(self.omega*self.time))
            vx = self.amp*self.omega*np.cos(self.kwave*x)*np.sinh(self.omega*self.time)/self.kwave
            ax = self.amp*self.omega*self.omega*np.cos(self.kwave*x)*np.cosh(self.omega*self.time)/self.kwave
        else:
            rho = self.amp*self.rho*np.cos(self.kwave*x)
            vx = 0
            ax = 0

        if ix == "x" and iy == "rho":
            return x,rho
        elif ix == "x" and iy == "vx":
            return x,vx
        elif ix == "x" and iy == "ax":
            return x,ax
        else:
            raise KeyError("There were errors in the quantity you requested")


#------------------------------------------------------------------------------
class noh (AnalyticalSolution):
    '''Analytical solution for the Noh problem test.'''

    def __init__(self, sim, time):
        AnalyticalSolution.__init__(self)
        self.time = time

        # Extract the parameters
        simfloatparams = sim.simparams.floatparams
        self.radius = simfloatparams["radius"]
        self.rho    = simfloatparams["rhofluid1"]
        self.time   = time
        self.ndim   = sim.ndims
        self.iMAX   = 1000
        self.gamma  = simfloatparams["gamma_eos"]

    def compute(self, x, y):
        '''Computes the exact solution of the Noh problem for
        1, 2 and 3 dimensions'''
        r = np.linspace(0.0,self.radius,self.iMAX)
        rho = self.rho*np.ones(self.iMAX)
        if self.time > 0:
            bound = 0.3333333333333*self.time
            inside = r<bound
            outside = r>bound
            if self.ndim ==1:
                rho[inside]=4.0*self.rho
                rho[outside]=self.rho
            elif self.ndim ==2:
                rho[inside]=16.0*self.rho
                rho[outside]=self.rho*(1.0 + self.time/r[outside])
            elif self.ndim ==3:
                rho[inside]=64.0*self.rho
                rho[outside]=self.rho*(1.0 + self.time/r[outside])**2.0
        return r,rho


#------------------------------------------------------------------------------
class shocktube (AnalyticalSolution):
    '''Analytical solution for the 1D shock tube test.'''

    def __init__(self, sim, time):
        AnalyticalSolution.__init__(self)
        self.time = time

        # Extract the parameters
        simfloatparams = sim.simparams.floatparams
        self.RHOinL = simfloatparams["rhofluid1"]
        self.RHOinR = simfloatparams["rhofluid2"]
        self.UinL   = simfloatparams["vfluid1[0]"]
        self.UinR   = simfloatparams["vfluid2[0]"]
        self.PinL   = simfloatparams["press1"]
        self.PinR   = simfloatparams["press2"]
        self.xL     = simfloatparams["boxmin[0]"]
        self.xR     = simfloatparams["boxmax[0]"]
        self.x0     = 0.5*(self.xL + self.xR)
        self.time   = time
        self.iMAX   = 16384
        if sim.simparams.stringparams["gas_eos"] == "isothermal":
            self.gamma = 1.0 + 1e-5

            cs2 = simfloatparams['temp0'] / simfloatparams['mu_bar']

            self.PinL   = self.RHOinL * cs2
            self.PinR   = self.RHOinR * cs2
        else:
            self.gamma = simfloatparams["gamma_eos"]

    def compute(self, x, y):
        '''Computes the exact solution of the Riemann problem.
        Gets passed two strings with the quantities that are needed
        (e.g., \'rho\' and \'pressure\').'''

        shocktube = ShocktubeSolution(self.RHOinL, self.RHOinR, self.UinL, self.UinR,
                                      self.PinL, self.PinR, self.xL, self.x0, self.xR,
                                      self.time, self.iMAX, self.gamma)

        if shocktube.single:
            prec_type=np.float32
        else:
            prec_type=np.float64
        xdata = np.zeros(self.iMAX, dtype=prec_type)
        ydata = np.zeros(self.iMAX, dtype=prec_type)

        xdata = shocktube.ComputeShocktubeSolution(x, self.iMAX)
        ydata = shocktube.ComputeShocktubeSolution(y, self.iMAX)

        return xdata, ydata



#------------------------------------------------------------------------------
class soundwave (AnalyticalSolution):
    '''Analytical solution for the 1D sound wave perturbation test.'''

    def __init__(self, sim, time):
        AnalyticalSolution.__init__(self)
        self.time = time

        # Extract the parameters
        simfloatparams = sim.simparams.floatparams
        self.rho    = simfloatparams["rhofluid1"]
        self.press  = simfloatparams["press1"]
        self.temp0  = simfloatparams["temp0"]
        self.mu_bar = simfloatparams["mu_bar"]
        self.amp    = simfloatparams["amp"]
        self.xL     = simfloatparams["boxmin[0]"]
        self.xR     = simfloatparams["boxmax[0]"]
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


#-------------------------------------------------------------------------------
class SedovSolution(object):
    """Class for simple approximations to the Sedov (1959) solution
    for point explosion.

    Solutions are taken from Book (1991), who present solutions by KorobeYnikov
    et al. (1961)

    args:
        E : float
            explosion energy
        rho : float
              density of the medium (initially constant everywhere)
        gamma : float, default=1.4
                ratio of specific heats
        nu : int, default=3
               number of dimensions which the explosion occurs in
        w : float, default=0
               power in the density law, rho=rho0 r^-w
    """
    def __init__(self, E, rho, gamma=1.4, nu=3,w=0., tiny=1e-50):
        self._tiny = tiny
        self._E = E
        self._gamma = gamma

        self._rho0 = rho
        self._rho1 = ((gamma + 1.)/(gamma - 1.))*rho

        self._nDim = nu
        self._w = w

        # Constants for the parametic equations:
        w1 = (3*nu - 2 + gamma*(2-nu))/(gamma + 1.)
        w2 = (2.*(gamma-1) + nu)/gamma
        w3 = nu*(2.-gamma)

        b0 = 1./(nu*gamma - nu + 2)
        b2 = (gamma-1.)/(gamma*(w2-w))
        b3 = (nu-w)/(float(gamma)*(w2-w))
        b5 = (2.*nu-w*(gamma+1))/(w3-w)
        b6 = 2./(nu+2-w)
        b1 = b2 + (gamma+1.)*b0 - b6
        b4 = b1*(nu-w)*(nu+2.-w)/(w3-w)
        b7 = w*b6
        b8 = nu*b6

        # C0 is surface area function in ndim. Use approx:
        if any(nu == np.array([1,2,3])):
            C0 = 2*(nu-1)*np.pi + (nu-2)*(nu-3)
        else: # General form requires Gamma function
            from scipy.special import gamma as GammaF
            C0 = (2**nu)*(np.pi**((nu-1)/2.))*GammaF((nu+1)/2.)/GammaF(nu)
        C5 = 2./(gamma - 1)
        C6 = (gamma + 1)/2.
        C1 = C5*gamma
        C2 = C6/gamma
        C3 = (nu*gamma - nu + 2.)/((w1-w)*C6)
        C4 = (nu + 2. - w)*b0*C6

        # Lambdas for setting up the interpolating functions:
        ETA = lambda F: (F**-b6)*((C1*(F-C2))** b2          )*((C3*(C4-F))**(  -b1  ))
        D   = lambda F: (F**-b7)*((C1*(F-C2))**(b3-  w  *b2))*((C3*(C4-F))**(b4+w*b1))*((C5*(C6-F))**-b5)
        P   = lambda F: (F** b8)*((C3*(C4-F))**(b4+(w-2)*b1))*((C5*(C6-F))**(1 -  b5))
        V   = lambda F: ETA(F)*F

        # Characterize the solution
        if w1 > w:
            Fmin = C2
        else:
            Fmin = C6

        F = np.logspace(np.log10(Fmin),0,1e5)

        # Sort the etas for our interpolation function
        eta = ETA(F)
        F = F[eta.argsort()]
        eta.sort()

        d = D(F)
        p = P(F)
        v = V(F)

        # If min(eta) != 0 then all values for eta < min(eta) = 0
        if eta[0] > 0:
            e01 = [0., eta[0]*(1-1e-10)]
            d01 = [0.,0]
            p01 = [0.,0]
            v01 = [0.,0]

            eta = np.concatenate([np.array(e01),eta])
            d   = np.concatenate([np.array(d01),  d])
            p   = np.concatenate([np.array(p01),  p])
            v   = np.concatenate([np.array(v01),  v])

        # Set up our interpolation functions
        self._d = interp1d(eta,d,bounds_error=False, fill_value=1./self._rho1)
        self._p = interp1d(eta,p,bounds_error=False, fill_value=0.)
        self._v = interp1d(eta,v,bounds_error=False, fill_value=0.)

        # Finally Calculate the normalization of R_s:
        I = eta**(nu-1)*(d*v**2 + p)

        I    = 0.5*(I[1: ]  +   I[:-1])
        deta =     (eta[1:] - eta[:-1])

        alpha = (I*deta).sum()*(8*C0)/((gamma**2-1.)*(nu+2.-w)**2)
        self._C = (1./alpha)**(1./(nu+2-w))

    # Shock properties
    def R_s(self,t):
        """Outer radius at time t"""
        t = np.maximum(t, self._tiny)
        return self._C *(self.E*t**2/self.rho0)**(1./(self._nDim + 2-self._w))

    def V_s(self,t):
        """Velocity of the shock wave"""
        t = np.maximum(t, self._tiny)
        return (2./(self._nDim+2-self._w)) * self.R_s(t) / t

    def P_s(self,t):
        """Post shock pressure"""
        return (2./(self.gamma+1))*self.rho0*self.V_s(t)**2

    @property
    def Rho_s(self,t=0):
        """Post shock density"""
        return self._rho1

    # Position dependent variables
    def rho(self,r,t):
        """Density at radius, r, and time, t."""
        eta = r/self.R_s(t) ;
        return self.Rho_s*self._d(eta)

    def P(self,r,t):
        """pressure at radius, r, and time, t."""
        eta = r/self.R_s(t) ;
        return self.P_s(t)*self._p(eta)

        # Position dependent variable
    def v(self,r,t):
        """Density at radius, r, and time, t."""
        eta = r/self.R_s(t) ;
        return self._v(eta)*(2/(self.gamma+1))*self.V_s(t)

    def u(self, r,t):
        """Internal energy at radius, r, and time, t."""
        return self.P(r,t)/(self.rho(r,t)*(self.gamma-1))

    def Entropy(self,r,t):
        """Entropy at radius, r, and time, t."""
        return self.P(r,t)/self.rho(r,t)**self.gamma

    # Other properties
    @property
    def E(self):
        """Total energy"""
        return self._E

    @property
    def gamma(self):
        """Ratio of specific heats"""
        return self._gamma

    @property
    def rho0(self):
        """Background density"""
        return self._rho0



class sedov (AnalyticalSolution):
    '''Analytical solution for the sedov blast wave problem'''

    def __init__(self, sim, time):
        fparams = sim.simparams.floatparams
        iparams = sim.simparams.intparams

        rho0 = fparams['rhofluid1']
        E0 = 1
        gamma = fparams['gamma_eos']
        w = 0 # Power law index
        ndim = iparams['ndim']

        self._sol = SedovSolution(E0, rho0, gamma=gamma,w=w,nu=ndim)

        # Save domain
        Rmax = 0;
        for i in range(ndim):
            Rmax += (0.5*(fparams['boxmax['+str(i)+']'] -
                          fparams['boxmin['+str(i)+']']))**2
        Rmax = Rmax**0.5

        self._r = np.linspace(0, Rmax, 1001)[1:]
        self._t = time

    def compute(self, x, y):
        return map(self._get_data, [x,y])

    def _get_data(self, x):
        if x == 'R':
            return self._r
        elif x == 'vR2d':
            return self._sol.v(self._r, self._t)
        elif x == 'vr':
            return self._sol.v(self._r, self._t)
        elif x == 'rho':
            return self._sol.rho(self._r, self._t)
        elif x == 'press':
            return self._sol.P(self._r, self._t)
        elif x == 'u':
            return self._sol.u(self._r, self._t)
        else:
            raise AttributeError("Sedov solution for variable %s not known"%x)
