#==============================================================================
#  disc.py
#  Contains class definitions for identifying gravitationally bound clumps
#  of gas to stars, as a simple method for identifying discs in simulations.
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
from data_fetcher import UserQuantity, SimBuffer
import pyximport; pyximport.install(setup_args={'include_dirs':[np.get_include()]})
import extract_disc_cython


#------------------------------------------------------------------------------
class Blob:
    '''Main class for defining a clump of gas in a simultion'''


    #--------------------------------------------------------------------------
    def __init__(self, ids, snap, type):
        self.ids = ids
        self.snap = snap
        self.type = type


    #--------------------------------------------------------------------------
    def n_particles(self):
        return self.ids.sum()


    #--------------------------------------------------------------------------
    def mass(self, unit='default'):
        '''Returns total mass of gas inside blob'''
        
        mass_info = UserQuantity('m').fetch(self.type, self.snap, unit=unit)
        mass_type = mass_info[1]
        scaling_factor = mass_info[2]
        return mass_type[self.ids].sum()*scaling_factor


    #--------------------------------------------------------------------------
    def SPH_positions_sim_frame(self,unit='default'):
        '''Returns array of positions of particles comprising gas blob'''
        
        ndim = self.snap.ndim
        positions = np.zeros((self.n_particles(),ndim))
        first_coordinate_info = UserQuantity('x').fetch(self.type,self.snap,unit=unit)
        scaling_factor = first_coordinate_info[2]
        positions[:,0] = first_coordinate_info[1][self.ids]
        if ndim>1:
            positions[:,1] = UserQuantity('y').fetch(self.type,self.snap,unit=unit)[1][self.ids]
        if ndim>2:
            positions[:,2] = UserQuantity('z').fetch(self.type,self.snap,unit=unit)[1][self.ids]
        return positions*scaling_factor


    #--------------------------------------------------------------------------
    def SPH_velocities_sim_frame(self,unit='default'):
        '''Returns array of velocities of particles comprising gas blob'''
        
        ndim = self.snap.ndim
        velocities = np.zeros((self.n_particles(),ndim))
        first_coordinate_info = UserQuantity('vx').fetch(self.type,self.snap,unit=unit)
        scaling_factor = first_coordinate_info[2]
        velocities[:,0] = first_coordinate_info[1][self.ids]
        if ndim>1:
            velocities[:,1] = UserQuantity('vy').fetch(self.type,self.snap,unit=unit)[1][self.ids]
        if ndim>2:
            velocities[:,2] = UserQuantity('vz').fetch(self.type,self.snap,unit=unit)[1][self.ids]
        return velocities*scaling_factor


#------------------------------------------------------------------------------
class Ambient_gas(Blob):
    pass


#------------------------------------------------------------------------------
class Disc(Blob):
    
    def __init__(self, star, ids, snap, type):
        Blob.__init__(self, ids, snap, type)
        self.star=star
    
    def radius(self, lagradius=0.5):
        raise NotImplementedError
    
    def angular_momentum(self):
        raise NotImplementedError
    
    def rotation_axis(self):
        raise NotImplementedError
    
    def surface_density(self):
        raise NotImplementedError
    
    def SPH_positions_star_frame(self):
        raise NotImplementedError
        
    def SPH_velocities_star_frame(self):
        raise NotImplementedError
    

#------------------------------------------------------------------------------
def extract_discs (snapno, sim, type='default', eccenlimit=0.9,
                   distancelimit=1., limiteigenvalues=0.2):
    '''This function takes a snapshot and a simulation number ("current" is
    also fine) and looks for which particles are bound to the stars. It
    returns a tuple, consisting of an Ambient_gas object (representing the
    gas that is not bound to any star) and of a list of Disc objects, one for
    each star.
    '''
    
    snap = SimBuffer.get_snapshot_extended(sim,snapno)

    parameters = dict(
                  eccenlimit = eccenlimit,
                  distancelimit = distancelimit,
                  limiteigenvalues = limiteigenvalues,
                  )

    # Query the number of dimensions
    ndim = snap.ndim
    
    # Fetch the data - we use code units
    
    # First extract coordinates, velocities and masses for the given type
    x_type = UserQuantity('x').fetch(type, snap)[1]
    vx_type = UserQuantity('vx').fetch(type, snap)[1]
    x_star = UserQuantity('x').fetch('star',snap)[1]
    vx_star = UserQuantity('vx').fetch('star',snap)[1]
    
    if ndim>1:
        y_type = UserQuantity('y').fetch(type, snap)[1]
        vy_type = UserQuantity('vy').fetch(type, snap)[1]
        y_star = UserQuantity('y').fetch('star',snap)[1]
        vy_star = UserQuantity('vy').fetch('star',snap)[1]
    
    if ndim>2:
        z_type = UserQuantity('z').fetch(type, snap)[1]
        vz_type = UserQuantity('vz').fetch(type, snap)[1]
        z_star = UserQuantity('z').fetch('star',snap)[1]
        vz_star = UserQuantity('vz').fetch('star',snap)[1]
    
    m_type = UserQuantity('m').fetch(type, snap)[1]
    m_star = UserQuantity('m').fetch('star',snap)[1]
    
    n_star = snap.GetNparticlesType('star')
    
    if ndim==2:
        owner = extract_disc_cython.flag_owner2d(x_type,y_type,vx_type,vy_type,m_type,x_star,y_star,vx_star,vy_star,m_star,parameters)
    elif ndim==3:
        owner = extract_disc_cython.flag_owner3d(x_type,y_type,z_type,vx_type,vy_type,vz_type,m_type,x_star,y_star,z_star,vx_star,vy_star,vz_star,m_star,parameters)
    
    # Loops over all stars and creates the disc objects
    disclist=[]
    for istar in range(n_star):
        ids = ( owner==istar )
        disc = Disc(istar, ids, snap, type)
        disclist.append(disc)
        print 'mass disc number', istar, ':',disc.mass()
        
    # Create the ambient gas object
    ids = (owner==-1)
    ambient = Ambient_gas(ids, snap, type)
    print 'ambient gas mass:', ambient.mass()
    
    # Return a tuple with the ambient gas and a list of discs
    return (ambient, disclist)
