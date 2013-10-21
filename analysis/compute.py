#==============================================================================
#  compute.py
#  ..
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
from gandalf.analysis.facade import get_sim_no, Singletons
import commandsource as Commands
import numpy as np
from scipy import interpolate
from data_fetcher import UserQuantity

'''This module collects helper functions to compute useful quantities'''


#------------------------------------------------------------------------------
def particle_data(snap, quantity, type="default", id=None):
    '''Return for a given snapshot, a given quantity of a given type. If id is
not specified, return an array with the quantity for each particle. Otherwise,
return a scalar quantity only for the given particle.
'''
    values = UserQuantity(quantity).fetch(type, snap)[1]
    if id == None:
        return values
    else:
        return values[id]


#------------------------------------------------------------------------------
def COM(snap, quantity='x', type="default"):
    x = UserQuantity(quantity).fetch(type, snap)[1]
    m = UserQuantity('m').fetch(type, snap)[1]
    
    return (x*m).sum()/m.sum()



#------------------------------------------------------------------------------
def L1errornorm(ic, x=None, y=None, xmin=None, xmax=None, 
                sim="current", snap="current"):
    '''Computes the L1 error norm from the simulation data relative to the analytical solution'''
    
    # Get the simulation number from the buffer
    simno = get_sim_no(snap)
    
    # Instantiate and setup the command object to retrieve analytical solution
    command1 = Commands.AnalyticalPlotCommand(x, y, ic, snap, simno)
    adata = command1.prepareData(Singletons.globallimits)

    # Instantiate and setup the 2nd command object to retrieve particle data
    command2 = Commands.ParticlePlotCommand(x, y, "sph", snap, simno)
    pdata = command2.prepareData(Singletons.globallimits)

    # Cut arrays if limits are provided
    if xmin != None and xmax != None:
        aindex = np.logical_and( adata.x_data > xmin, adata.x_data < xmax)
        adata.x_data = adata.x_data[aindex]
        adata.y_data = adata.y_data[aindex]
        pindex = np.logical_and( pdata.x_data > adata.x_data.min(),
                                    pdata.x_data < adata.x_data.max() )
        pdata.x_data = pdata.x_data[pindex]
        pdata.y_data = pdata.y_data[pindex]
   
    # Prepare interpolation function from analytical data
    #f = interpolate.interp1d(adata.x_data[::-1], adata.y_data[::-1], kind = 'linear', axis=0, bounds_error = False)
    f = interpolate.interp1d(adata.x_data[::], adata.y_data[::], kind = 'linear', axis=0, bounds_error = False)

    # Compute error norm of particle data relative to analytical data
    L1 = np.linalg.norm((pdata.y_data - f(pdata.x_data)), ord=1)/pdata.x_data.size
    return L1


#------------------------------------------------------------------------------
def lagrangian_radii(snap, mfrac=0.5, type="default"):
    '''Computes the Lagrangian radii from all particles in simulation'''
    
    r = UserQuantity('r').fetch(type, snap)[1]
    m = UserQuantity('m').fetch(type, snap)[1]

    # Find particle ids in order of increasing radial distance
    porder = np.argsort(r)
    m_ordered = m[porder]
    mcumulative = np.cumsum(m_ordered)
    mtotal = mcumulative[-1]
    mlag = mfrac*mtotal
    index = np.searchsorted(mcumulative,mlag)
    return 0.5*(r[porder[index-1]] + r[porder[index]])


