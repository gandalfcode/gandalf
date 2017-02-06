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
import random
from math import *
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
    z  = UserQuantity("z").fetch(type, snap)[1]
    vx = UserQuantity("vx").fetch(type, snap)[1]
    vy = UserQuantity("vy").fetch(type, snap)[1]
    vz = UserQuantity("vz").fetch(type, snap)[1]
    n  = x.size
    r  = np.zeros(n)
    vsqd = np.zeros(n)
    

    # Create logarithmic bins based on inputted range
    bins = np.linspace(log10(rmin),log10(rmax),nbin+1)
    vsqdmean = np.zeros(nbin+2)
    npart = np.zeros(nbin+2)
    

    # Loop through a random selection of points
    for j in range(npoints):
        i = random.randrange(0,n-1)
        x0 = x[i]; y0 = y[i]; z0 = z[i]
        vx0 = vx[i]; vy0 = vy[i]; vz0 = vz[i]

        r[:] = (x[:] - x0)**2
        r[:] += (y[:] - y0)**2
        r[:] += (z[:] - z0)**2
        r = np.log10(np.sqrt(r))

        vsqd[:] = (vx[:] - vx0)**2
        vsqd[:] += (vy[:] - vy0)**2
        vsqd[:] += (vz[:] - vz0)**2

        # Now discretise values into bins and sum up values
        binpos = np.digitize(r,bins,right=True)
        for jj in range(npoints):
            ii = random.randrange(0,n-1)
            vsqdmean[binpos[ii]] += vsqd[ii]
            npart[binpos[ii]] += 1


    # Normalise velocity bins to arithmetic mean values
    for j in range(nbin):
        if npart[j] > 0: vsqdmean[j] = vsqdmean[j]/npart[j]

    vsqdmean = np.log10(vsqdmean)


    # Print out results when finished (for now)
    print bins
    print vsqdmean
    print npart

    return bins[0:nbin],vsqdmean[0:nbin]



#------------------------------------------------------------------------------
def density_pdf(snap, type="default", nbin=32, rhomin="auto", rhomax="auto"):
    '''Calculate the probability density function of the density field'''

    rho = UserQuantity("rho").fetch(type, snap)[1]

    if rhomin == "auto": rhomin = np.amin(rho)
    if rhomax == "auto": rhomax = np.amax(rho)
    n = rho.size
    rho = np.log10(rho)

    bins = np.linspace(log10(rhomin),log10(rhomax),nbin+1)
    #rhopdf = np.zeros(nbin+2)

    rhopdf = np.histogram(rho,bins=bins)[0]
    #binpos = np.digitize(rho,bins,right=True)
    #for i in range(n):
    #    rhopdf[binpos[i]] += 1

    prectype = rho.dtype
    rhopdf = rhopdf.astype(prectype)

    # Normalise logarithmic bins
    for j in range(nbin):
        rhopdf[j] = rhopdf[j]/(10**bins[j+1] - 10**bins[j])

    print bins
    print rhopdf

    return bins[0:nbin],rhopdf
