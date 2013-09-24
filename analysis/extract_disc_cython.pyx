import cython
import numpy as np
cimport numpy as np
from libc.math cimport sqrt

cdef double G=1


cpdef flag_owner(np.ndarray[double, ndim=2] datasph, np.ndarray[double, ndim=2] datastar, parameters):
    """
    Finds the star owner of each sph particle. Does not check for binaries.
    Returns an array containing the id of the owner for each particle, or -1 if
    it is unbound.
    
    The owner star is simply flagged as the one with the smallest binding energy. There is also a
    limit on the allowed eccentricity and distance from the star, that are read from the parameters dictionary.
    """
    cdef double xp, yp, zp
    cdef double xrel, yrel, zrel
    cdef int iparticle, istar
    cdef int nstars, nparticles
    cdef double vxp, vyp, vzp
    cdef double distance, relvel
    cdef double xrelvel, yrelvel, zrelvel
    cdef double mstar, en, enmin
    #dimensionless units, so G=1
    cdef double costheta, sintheta, eccen
    cdef double eccenlimit = parameters["eccenlimit"]
    cdef double distancelimit = parameters["distancelimit"]

    nstars=datastar[:,0].size
    nparticles=datasph[:,0].size

    
    #for each particle, finds the star to which it is mostly bound
    owner = - np.ones(nparticles, dtype='int32')
    for iparticle in range(nparticles):
        enmin=0.
        xp=datasph[iparticle,0]
        yp=datasph[iparticle,1]
        zp=datasph[iparticle,2]
        vxp=datasph[iparticle,3]
        vyp=datasph[iparticle,4]
        vzp=datasph[iparticle,5]
        for istar in range(nstars):
            xrel=xp-datastar[istar,0]
            yrel=yp-datastar[istar,1]
            zrel=zp-datastar[istar,2]
            distance=sqrt(xrel**2+yrel**2+zrel**2)
            xrelvel=vxp-datastar[istar,3]
            yrelvel=vyp-datastar[istar,4]
            zrelvel=vzp-datastar[istar,5]
            relvel=sqrt(xrelvel**2+yrelvel**2+zrelvel**2)
            mstar = datastar[istar,6]
            en = 0.5 * relvel**2 - G * mstar / distance
            costheta = (xrel*xrelvel+yrel*yrelvel+zrel*zrelvel)/distance/relvel
            sintheta = sqrt(1-costheta**2)
            eccen = distance * relvel * sintheta
            #energy criterion
            if (en < enmin and eccen < eccenlimit and distance<distancelimit):
                enmin = en
                owner[iparticle]=istar
                
    return owner

@cython.wraparound(False)
@cython.boundscheck(False)
def get_disc_array(int istar, np.ndarray[int, ndim=1] owner, np.ndarray[double, ndim=2] datasph, np.ndarray[double, ndim=2] datastar):
    '''
    Scans the datasph array and retrieves the particles belonging to the istar star. Note that this can be used also when the gas particles
    have been classified as ambient gas. If the particles belong to a real star, shifts the coordinates and the velocity relatively to the star.
    Returns an array with their data.
    '''
    cdef int ndisc
    cdef int iparticle
    cdef int idisc
    
    ndisc=0
    for iparticle in range(owner.size):
        if owner[iparticle] == istar:
            ndisc += 1 
    
    discdata = np.zeros((ndisc,7))
    idisc = 0
    for iparticle in range(owner.size):
        if owner[iparticle] == istar:
            discdata[idisc,:] = datasph[iparticle, :]
            if istar != -1:
                discdata[idisc, :6] -= datastar[istar, :6]
            idisc += 1
            
    return discdata

cpdef compute_energy_to_cluster(np.ndarray[double, ndim=2] dataparticles, np.ndarray[double, ndim=2] datastar):
    '''Computes the energy with respect to the cluster for all the particles in array dataparticles.
    Returns an array containing all the energies.'''
    cdef int iparticle, istar
    cdef int nparticles, nstars
    cdef double en, xp, yp, zp, vxp, vyp, vzp
    cdef double xrel, distance, xrelvel, yrelvel, zrelvel
    cdef double yrel, zrel
    cdef double relvel, mstar
    
    nparticles=dataparticles[:,0].size
    nstars=datastar[:,0].size
    energy=np.zeros(nparticles)
    
    unboundparticles=[]
    
    for iparticle in range(nparticles):
        en=0
        xp=dataparticles[iparticle,0]
        yp=dataparticles[iparticle,1]
        zp=dataparticles[iparticle,2]
        vxp=dataparticles[iparticle,3]
        vyp=dataparticles[iparticle,4]
        vzp=dataparticles[iparticle,5]
        for istar in range(nstars):
            xrel=xp-datastar[istar,0]
            yrel=yp-datastar[istar,1]
            zrel=zp-datastar[istar,2]
            distance=sqrt(xrel**2+yrel**2+zrel**2)
            xrelvel=vxp-datastar[istar,3]
            yrelvel=vyp-datastar[istar,4]
            zrelvel=vzp-datastar[istar,5]
            relvel=sqrt(xrelvel**2+yrelvel**2+zrelvel**2)
            mstar = datastar[istar,6]
            en += 0.5 * relvel**2 - G * mstar / distance
        energy[iparticle] = en
            
    return energy
