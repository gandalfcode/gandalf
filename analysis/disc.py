import numpy as np
from data_fetcher import UserQuantity
import pyximport; pyximport.install(setup_args={'include_dirs':[np.get_include()]})
import extract_disc_cython

class Blob:
    def __init__(self, ids, snap, type):
        self.ids=ids
        self.snap=snap
        self.type=type
    
    def n_particles(self):
        return self.ids.sum()
    
    def mass(self, unit='default'):
        mass_info=UserQuantity('m').fetch(self.type, self.snap, unit=unit)
        mass_type=mass_info[1]
        scaling_factor=mass_info[2]
        return mass_type[self.ids].sum()*scaling_factor
    
    def SPH_positions_sim_frame(self,unit='default'):
        ndim=self.snap.ndim
        positions = np.zeros((self.n_particles(),ndim))
        first_coordinate_info=UserQuantity('x').fetch(self.type,self.snap,unit=unit)
        scaling_factor=first_coordinate_info[2]
        positions[:,0]=first_coordinate_info[1][self.ids]
        if ndim>1:
            positions[:,1]=UserQuantity('y').fetch(self.type,self.snap,unit=unit)[1][self.ids]
        if ndim>2:
            positions[:,2]=UserQuantity('z').fetch(self.type,self.snap,unit=unit)[1][self.ids]
        return positions*scaling_factor
    
    def SPH_velocities_sim_frame(self,unit='default'):
        ndim=self.snap.ndim
        velocities = np.zeros((self.n_particles(),ndim))
        first_coordinate_info=UserQuantity('vx').fetch(self.type,self.snap,unit=unit)
        scaling_factor=first_coordinate_info[2]
        velocities[:,0]=first_coordinate_info[1][self.ids]
        if ndim>1:
            velocities[:,1]=UserQuantity('vy').fetch(self.type,self.snap,unit=unit)[1][self.ids]
        if ndim>2:
            velocities[:,2]=UserQuantity('vz').fetch(self.type,self.snap,unit=unit)[1][self.ids]
        return velocities*scaling_factor

class Ambient_gas(Blob):
    pass

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
    


    
def extract_discs (snap, type='default', eccenlimit=0.9, distancelimit=1., limiteigenvalues=0.2):
    '''This function takes a snapshot and looks for which particles are bound
    to the stars. It returns a tuple, consisting of an Ambient_gas object
    (representing the gas that is not bound to any star) and of a list
    of Disc objects, one for each star.
    '''
    
    parameters = dict(
                  eccenlimit = eccenlimit,
                  distancelimit = distancelimit,
                  limiteigenvalues = limiteigenvalues,
                  )
    
    #query the number of dimensions
    ndim=snap.ndim
    
    #fetch the data - we use code units
    
    #first extract coordinates, velocities and masses for the given type
    x_type = UserQuantity('x').fetch(type, snap)[1]
    vx_type = UserQuantity('vx').fetch(type, snap)[1]
    x_star=UserQuantity('x').fetch('star',snap)[1]
    vx_star=UserQuantity('vx').fetch('star',snap)[1]
    
    if ndim>1:
        y_type = UserQuantity('y').fetch(type, snap)[1]
        vy_type = UserQuantity('vy').fetch(type, snap)[1]
        y_star=UserQuantity('y').fetch('star',snap)[1]
        vy_star=UserQuantity('vy').fetch('star',snap)[1]
    
    if ndim>2:
        z_type = UserQuantity('z').fetch(type, snap)[1]
        vz_type = UserQuantity('vz').fetch(type, snap)[1]
        z_star=UserQuantity('z').fetch('star',snap)[1]
        vz_star=UserQuantity('vz').fetch('star',snap)[1]
    
    m_type = UserQuantity('m').fetch(type, snap)[1]
    m_star=UserQuantity('m').fetch('star',snap)[1]
    
    n_star = snap.GetNparticlesType('star')
    
    if ndim==2:
        owner = extract_disc_cython.flag_owner2d(x_type,y_type,vx_type,vy_type,m_type,x_star,y_star,vx_star,vy_star,m_star,parameters)
    elif ndim==3:
        owner = extract_disc_cython.flag_owner3d(x_type,y_type,z_type,vx_type,vy_type,vz_type,m_type,x_star,y_star,z_star,vx_star,vy_star,vz_star,m_star,parameters)
    
    #now loops over the stars and create the discs
    disclist=[]
    for istar in range(n_star):
        ids= ( owner==istar )
        disc=Disc(istar, ids, snap, type)
        disclist.append(disc)
        print 'mass disc number', istar, ':',disc.mass()
        
    #create the ambient gas object
    ids= (owner==-1)
    ambient=Ambient_gas(ids, snap, type)
    print 'ambient gas mass:', ambient.mass()
    
    #return a tuple with the ambient gas and a list of discs
    return (ambient, disclist)