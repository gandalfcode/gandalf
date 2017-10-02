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
import matplotlib as mpl
mpl.rcParams['figure.figsize'] = 6.64, 4.98 # twice the mnras fig size
#                                             for a 1-column figure
mpl.rcParams['font.size'] = 16
mpl.rcParams['lines.linewidth'] = 2
mpl.rcParams['axes.linewidth'] = 1
mpl.rcParams['legend.numpoints'] = 1
mpl.rcParams['legend.handlelength'] = 2
mpl.rcParams['legend.frameon'] = 0
mpl.rcParams['legend.fontsize'] = 16
mpl.rcParams['lines.markeredgewidth'] = 2
mpl.rcParams['lines.markersize'] = 6
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42


class DriftVelocitySolution(object):
    def __init__(self, K, vg, vd, eps=None, rho=1):
        self._K =  K
        self._vg0 = vg 
        self._vd0 = vd
        self._dv0 = vd - vg 
        self._rho = rho
        
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

    def Ek(self, t):
        eps = self._eps
        return 0.5*self._rho*(self.vg(t)**2 * (1-eps) + self.vd(t)**2 * eps)
    
def DriftVelocitySolutionFactory():
    params = SimBuffer.get_current_sim().simparams
    K    = params.floatparams['drag_coeff']
    vg   = params.floatparams['vfluid1[0]']
    vd   = params.floatparams['vfluid2[0]']
    rhog = params.floatparams['rhofluid1']
    if params.stringparams['dust_forces'] == 'test_particle':
        eps = 0 
    elif params.stringparams['dust_forces'] == 'full_twofluid':
        rhod = params.floatparams['dust_mass_factor']
        eps = rhod / (1. + rhod)
    else:
        raise AttributeError('Dust simulation type not recognises')
    sol = DriftVelocitySolution(K, vg, vd, eps=eps, rho=rhog/(1-eps))
    return sol

def volume():
    params = SimBuffer.get_current_sim().simparams
    ndim = params.intparams['ndim']
    fp = params.floatparams
    V = 1.
    for i in range(ndim):
        V *= fp['boxmax[%d]'%i] - fp['boxmin[%d]'%i]

    return V
    
def Compute_Energy_Evolution():
    t  = []
    Ek = []
    U  = []
    s = snap(0)
    V = volume()
    while True:
        t.append(s.t)

        m = get_data('m', type='sph')
        U.append((m*get_data('u', type='sph')).sum()/V)

        E = 0
        for tp in ['sph', 'dust']:
            m = get_data('m', type=tp)
            for j in range(sim.simparams.intparams['ndim']):
                E += (0.5*m*get_data('v'+chr(ord('x')+j),type=tp)**2).sum()/V
        Ek.append(E)

        try:
            s = next()
        except:
            break

    return map(np.array, [t, Ek, U])


# Set up sims with both SPH and MFV

Ks = [0.01, 0.1, 1., 10., 100.]
schemes = ['sph', 'mfvmuscl']
colours = ['ko', 'cx']

tmax = 1.0
ti = np.linspace(0, tmax, 10**3)

for c, s in zip(colours, schemes):
    for K in Ks:
        print s, K
        sim = newsim("dustybox.dat", sim=s)
        sim.SetParam("drag_coeff", K)
        setupsim()
        run()

        sol=DriftVelocitySolutionFactory()
        t, Ek, U = Compute_Energy_Evolution()

        if s == 'sph': plt.plot(ti, sol.Ek(ti), ':', c='0.5')
        plt.plot(t, Ek, c, markerfacecolor='none')
        plt.xlabel('$t$')
        plt.ylabel('$E_k$')

plt.subplots_adjust(top=0.98, bottom=0.11,
                    left=0.15, right=0.98)
plt.savefig('dustybox.pdf')
plt.show()
