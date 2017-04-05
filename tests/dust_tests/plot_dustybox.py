#==============================================================================
#  plot_dustybox.py
#  Plots the time evolution of a the dust and gas velocity for a single particle
#  in the DUSTYBOX problem.
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
from gandalf.analysis.facade import *
import numpy as np
import matplotlib.pyplot as plt
import sys

class DriftVelocitySolution(object):
    def __init__(self, K, vg, vd, eps=None):
        self._K =  K
        self._vg0 = vg 
        self._vd0 = vd
        self._dv0 = vd - vg 

        if eps:
            self._eps = eps
        else:
            self._eps = 0
        
    def vcom(self, t):
        eps = self._eps
        return ((1-eps)*self._vg0 + eps*self._vd0)

    def dv(self,t):
        x = t * self._K        
        eps = self._eps
        return self._dv0 * np.exp(-x)

    def vg(self, t):
        eps = self._eps
        return self.vcom(t) - eps * self.dv(t)

    def vd(self, t):
        eps = self._eps
        return self.vcom(t) + (1-eps) * self.dv(t)
    
def DriftVelocitySolutionFactory():
    params = SimBuffer.get_current_sim().simparams
    K  = params.floatparams['drag_coeff']
    vg = params.floatparams['vfluid1[0]']
    vd = params.floatparams['vfluid2[0]']
    if params.stringparams['dust_forces'] == 'test_particle':
        eps = 0 
    elif params.stringparams['dust_forces'] == 'full_twofluid':
        rhod = params.floatparams['dust_mass_factor']
        eps = rhod / (1. + rhod)
    else:
        raise AttributeError('Dust simulation type not recognises')
    sol = DriftVelocitySolution(K, vg, vd, eps=eps)
    return sol
    
if __name__=="__main__":
    # Load the sim
    if len(sys.argv) > 1:
        sim = loadsim(sys.argv[1])
    else:
        sim = loadsim('DUSTYBOX')

    # Set up the analytical solution
    sol=DriftVelocitySolutionFactory()   
    
    # Find the maximum time
    tmax = snap(-1).t
    t = np.linspace(0, tmax, 10**3)
    
    # Plot the gas
    time_plot('t', 'vx', id=0, type='sph')
    plt.plot(t, sol.vg(t), 'k')
    
    # Plot the dust
    time_plot('t', 'vx', id=0, type='dust', overplot=True)
    plt.plot(t, sol.vd(t), 'k--')
    
    plt.xlim(0, t[-1])
    limit('vx',min(sol.vd(t).min(), sol.vg(t).min()), max(sol.vd(t).max(), sol.vg(t).max()), window='all')

    # Check energy conservation
    t  = []
    Ek = []
    U  = []
    s = snap(0)
    while True:
        t.append(s.t)

        m = get_data('m', type='sph')
        U.append((m*get_data('u', type='sph')).sum())

        _E = 0
        for tp in ['sph', 'dust']:
            m = get_data('m', type=tp)
            for j in range(sim.simparams.intparams['ndim']):
                _E += (0.5*m*get_data('v' + chr(ord('x')+j),type=tp)**2).sum()
        Ek.append(_E)

        try:
            s = next()
        except:
            break


    Etot = np.array(Ek) + np.array(U)
    
    plt.figure()
    plt.subplot(211)
    plt.plot(t, 1 - Etot/Etot[0])
    plt.xlabel('t')
    plt.ylabel('$\Delta E/E$')

    plt.subplot(212)
    plt.plot(t, Ek-Ek[0], label='$E_k$')
    plt.plot(t, U[0] - U,  label='$-U$')
    plt.xlabel('t')
    plt.ylabel('$\Delta E$')
    plt.legend(loc='best')

    plt.show()
