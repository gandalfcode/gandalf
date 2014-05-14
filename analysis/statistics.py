#==============================================================================
#  statistics.py
#  Contains useful helper routines for calculating important statistical
#  quantities for individual snapshot files.
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
from gandalf.analysis.facade import Singletons, SimBuffer, BufferException
import commandsource as Commands
import numpy as np
from math import *
import random
from scipy import interpolate
from data_fetcher import UserQuantity

'''This module collects helper functions to compute useful quantities'''



#------------------------------------------------------------------------------
def structure_function(snap, type="default", nbin=8, npoints=1000,
                       rmin=0.001, rmax=10.0):
    '''Calculate the structure function for a given snapshot'''

    # Return all relevant particle data (positions and velocity)
    x  = UserQuantity("x").fetch(type, snap)[1]
    y  = UserQuantity("y").fetch(type, snap)[1]
    z  = UserQuantity("x").fetch(type, snap)[1]
    vx = UserQuantity("vx").fetch(type, snap)[1]
    vy = UserQuantity("vy").fetch(type, snap)[1]
    vz = UserQuantity("vz").fetch(type, snap)[1]
    n  = x.size
    r  = np.zeros(n)
    v  = np.zeros(n)
    

    # Create logarithmic bins based on inputted range
    bins = np.linspace(log10(rmin),log10(rmax),nbin+1)
    print bins
    vmean = np.zeros(nbin+2)
    npart = np.zeros(nbin+2)
    

    # Loop through a random selection of points
    for j in range(npoints):
    #for i in range(n):
        i = random.randrange(0,n-1)
        x0 = x[i]; y0 = y[i]; z0 = z[i]
        vx0 = vx[i]; vy0 = vy[i]; vz0 = vz[i]

        r[:] = (x[:] - x0)**2 + (y[:] - y0)**2 + (z[:] - z0)**2
        v[:] = (vx[:] - vx0)**2 + (vy[:] - vy0)**2 + (vz[:] - vz0)**2
        r = np.log10(np.sqrt(r))
        #v = np.sqrt(v)
        #v = np.log10(np.sqrt(v))

        # Now discretise values into bins and sum up values
        binpos = np.digitize(r,bins,right=True)
        for k in range(n):
            vmean[binpos[k]] += v[k]
            npart[binpos[k]] += 1

        #print r
        #print v
        #print pos


    # Normalise logarithmic bins
    for j in range(nbin):
        #vmean[j] = vmean[j]/(10**bins[j+1] - 10**bins[j])
        if npart[j] > 0: vmean[j] = vmean[j]/npart[j]

    vmean = np.log10(vmean)


    # Print out results when finished (for now)
    print bins
    print vmean
    print npart

    return bins[0:nbin],vmean[0:nbin]



#------------------------------------------------------------------------------
def density_pdf(snap, type="default", nbin=16, rhomin="auto", rhomax="auto"):
    '''Calculate the probability density function of the density field'''

    rho = UserQuantity("rho").fetch(type, snap)[1]

    if rhomin == "auto": rhomin = np.amin(rho)
    if rhomax == "auto": rhomax = np.amax(rho)

    bins = np.linspace(log10(rhomin),log10(rhomax),nbin+1)
    rhopdf = np.zeros(nbin+2)

    binpos = np.digitize(rho,bins,right=True)
    for i in range(n):
        rhopdf[binpos[i]] += 1

    # Normalise logarithmic bins
    for j in range(nbin):
        rhopdf[j] = rhopdf[j]/(10**bins[j+1] - 10**bins[j])

    print bins
    print rhopdf

    return bins[0:nbin],rhopdf[0:nbin]
