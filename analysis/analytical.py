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
import shocktub

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

    # TODO : This is not the freefall solution.  Need to update!!
    def compute(self, x, y):
        '''Computes exact solution for freefall collapse problem'''
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
        r = np.arange(0.0,self.radius,1.0/self.iMAX)
        rho = self.rho*np.ones(self.iMAX)
        if self.time > 0:
            bound = 0.3333333333333*self.time
            i = 0
            while i < self.iMAX:
                if self.ndim == 1:
                    if r[i] < bound: rho[i] = 4.0*self.rho
                    else: rho[i] = self.rho
                elif self.ndim == 2:
                    if r[i] < bound: rho[i] = 16.0*self.rho
                    else: rho[i] = self.rho*(1.0 + self.time/r[i])
                elif self.ndim == 3:
                    if r[i] < bound: rho[i] = 64.0*self.rho
                    else: rho[i] = self.rho*pow(1.0 + self.time/r[i],2.0)
                i = i + 1
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
        self.iMAX   = 50000
        if sim.simparams.stringparams["gas_eos"] == "isothermal":
            self.gamma = 1.0 + 1e-5
        else:
            self.gamma = simfloatparams["gamma_eos"]

    def compute(self, x, y):
        '''Computes the exact solution of the Riemann problem.
        Gets passed two strings with the quantities that are needed
        (e.g., \'rho\' and \'pressure\').'''

        # TODO : Needs to be updated; perhaps re-write in C++
        # Calls the fortran module that computes the state
        shocktub.shocktub(self.RHOinL, self.RHOinR, self.UinL, self.UinR,
                          self.PinL, self.PinR, self.xL, self.x0, self.xR,
                          self.time, self.iMAX, self.gamma)

        # Reads the data from the text file produced
        data = np.genfromtxt('sod.out',names=['x','rho','vx','press','u'])
        return data[x], data[y]


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
        if sim.simparams.stringparams["gas_eos"] == "isothermal":
            self.gamma = 1.0 + 1e-5
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
