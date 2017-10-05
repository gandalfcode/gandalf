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
from gandalf.analysis.facade import get_sim_no, Singletons, SimBuffer, BufferException
import commandsource as Commands
import numpy as np
from scipy import interpolate
from data_fetcher import UserQuantity
from swig_generated.SphSim import UnitInfo


'''This module collects helper functions to compute useful quantities'''


#------------------------------------------------------------------------------
def particle_data(snap, quantity, type="default", unit="default", id=None):
    '''Return for a given snapshot, a given quantity of a given type. If id is
not specified, return an array with the quantity for each particle. Otherwise,
return a scalar quantity only for the given particle.
'''
    unitinfo, values, scaling_factor, label = UserQuantity(quantity).fetch(type, snap, unit=unit)
    if values.size == 0:
        values_to_return = np.nan
    else:
        if id == None:
            values_to_return = values
        else:
            values_to_return = values[id]
    return unitinfo, values_to_return, scaling_factor, label


#------------------------------------------------------------------------------
def time_derivative(snap, quantity, type="default", unit="default", id=None):
    '''Return for a given snapshot, the time derivative of a given quantity
    of a given type. If possible, use central difference. Otherwise, use
    forward/backward difference. If id is not specified, return an array with
    the quantity for each particle. Otherwise, return a scalar quantity only
    for the given particle.
    '''
    from data_fetcher import get_time_snapshot

    if unit != "default":
        raise NotImplementedError("""time_derivative implemented only with default units
        at the moment!""")

    # Find previous and next snapshots.  If first or last snapshot, then
    # return the current snapshot to compute a value
    try:
        snap1 = SimBuffer.get_previous_snapshot_from_object(snap)
    except BufferException:
        snap1 = snap
    try:
        snap2 = SimBuffer.get_next_snapshot_from_object(snap)
    except BufferException:
        snap2 = snap

    # Return array of values of quantity.  If either are empty, return nan
    quantityunitinfo, values1, quantityscaling_factor, label = UserQuantity(quantity).fetch(type, snap1)
    values2 = UserQuantity(quantity).fetch(type, snap2)[1]

    timeunitinfo, time, timescaling_factor, tlabel = get_time_snapshot(snap)
    scaling_factor = quantityscaling_factor/timescaling_factor
    unitinfo = UnitInfo()
    unitinfo.name= quantityunitinfo.name + "_" + timeunitinfo.name
    unitinfo.label= quantityunitinfo.label + "\\ " + timeunitinfo.label + "^{-1}"

    # Calculate the time derivative with central difference and return value
    tdiff = snap2.t - snap1.t
    if values1.size == 0 or values2.size == 0:
        timederiv = np.nan
    else:
        if id == None:
            timederiv = (values2 - values1)/tdiff
        else:
            timederiv = (values2[id] - values1[id])/tdiff
    return unitinfo, timederiv, scaling_factor, label+"_t"


#------------------------------------------------------------------------------
def COM(snap, quantity='x', type="default", unit="default"):
    ''' Computes the centre-of-mass value of a given vector component'''
    xunitinfo, x, xscaling_factor, xlabel = UserQuantity(quantity).fetch(type, snap, unit)
    m = UserQuantity('m').fetch(type, snap)[1]

    com = (x*m).sum()/m.sum()

    return xunitinfo, com, xscaling_factor, xlabel+'_COM'


#------------------------------------------------------------------------------
def L1errornorm(ic, x=None, y=None, xmin=None, xmax=None, normalise=None,
                sim="current", snap="current", type="sph"):
    '''Computes the L1 error norm from the simulation data relative to the analytical solution'''

    # Get the simulation number from the buffer
    simno = get_sim_no(snap)

    # Instantiate and setup the command object to retrieve analytical solution
    command1 = Commands.AnalyticalPlotCommand(x, y, ic, snap, simno)
    adata = command1.prepareData(Singletons.globallimits)

    # Instantiate and setup the 2nd command object to retrieve particle data
    command2 = Commands.ParticlePlotCommand(x, y, type, snap, simno)
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

    # Normalise quantity to given average value (if one is given)
    if normalise != None:
        av = sum(pdata.y_data)/pdata.y_data.size
        pdata.y_data = pdata.y_data/av/normalise

    # Prepare interpolation function from analytical data
    #f = interpolate.interp1d(adata.x_data[::-1], adata.y_data[::-1], kind = 'linear', axis=0, bounds_error = False)
    f = interpolate.interp1d(adata.x_data[::], adata.y_data[::], kind = 'linear', axis=0, bounds_error = False)

    # Compute error norm of particle data relative to analytical data
    L1 = np.linalg.norm((pdata.y_data - f(pdata.x_data)), ord=1)/pdata.x_data.size
    return L1


#------------------------------------------------------------------------------
def lagrangian_radii(snap, mfrac=0.5, type="default", unit="default"):
    '''Computes the Lagrangian radii from all particles in simulation'''
    runitinfo, r, rscaling_factor, rlabel = UserQuantity('r').fetch(type, snap, unit)
    m = UserQuantity('m').fetch(type, snap, unit)[1]

    # Find particle ids in order of increasing radial distance
    porder      = np.argsort(r)
    m_ordered   = m[porder]
    mcumulative = np.cumsum(m_ordered)
    mtotal      = mcumulative[-1]
    mlag        = mfrac*mtotal
    index       = np.searchsorted(mcumulative,mlag)
    lag_radius  = 0.5*(r[porder[index-1]] + r[porder[index]])
    return runitinfo, lag_radius, rscaling_factor, 'lag_radius_' + str(mfrac)


#------------------------------------------------------------------------------
def energy_error(snap, etot0, type="default", unit="default"):
    '''Computes the energy error of all particles in the simulation'''
    vx   = UserQuantity('vx').fetch(type, snap, unit)[1]
    vy   = UserQuantity('vy').fetch(type, snap, unit)[1]
    vz   = UserQuantity('vz').fetch(type, snap, unit)[1]
    m    = UserQuantity('m').fetch(type, snap, unit)[1]
    gpot = UserQuantity('gpot').fetch(type, snap, unit)[1]
    N    = m.size

    # Loop over all particles and compute
    energy = 0.5*m*(vx*vx + vy*vy + vz*vz) - 0.5*m*gpot
    etot   = np.sum(energy)
    error  = abs((etot - etot0)/etot0)

    return error
