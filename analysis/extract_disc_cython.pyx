import cython
import numpy as np
cimport numpy as np
from libc.math cimport sqrt

cdef double G=1


cpdef flag_owner3d(float[:] x_type, float[:] y_type, float[:] z_type,
                   float[:] vx_type, float[:] vy_type, float[:] vz_type, float[:] m_type,
                   float[:] x_star, float[:] y_star, float[:] z_star,
                   float[:] vx_star, float[:] vy_star, float[:] vz_star,
                   float[:] m_star, parameters):
    """
    Finds the star owner of each sph particle. Does not check for binaries.
    Returns an array containing the id of the owner for each particle, or -1 if
    it is unbound.
    
    The owner star is simply flagged as the one with the smallest binding energy. There is also a
    limit on the allowed eccentricity and distance from the star, that are read from the parameters dictionary.
    """
    cdef float xp, yp, zp
    cdef float xrel, yrel, zrel
    cdef int iparticle, istar
    cdef int nstars, nparticles
    cdef float vxp, vyp, vzp
    cdef float distance, relvel
    cdef float xrelvel, yrelvel, zrelvel
    cdef float mstar, en, enmin
    #dimensionless units, so G=1
    cdef float costheta, sintheta, eccen
    cdef float eccenlimit = parameters["eccenlimit"]
    cdef float distancelimit = parameters["distancelimit"]

    nparticles=x_type.size
    nstars=x_star.size
    
    #for each particle, finds the star to which it is mostly bound
    owner = - np.ones(nparticles, dtype='int32')
    for iparticle in range(nparticles):
        enmin=0.
        xp=x_type[iparticle]
        yp=y_type[iparticle]
        zp=z_type[iparticle]
        vxp=vx_type[iparticle]
        vyp=vy_type[iparticle]
        vzp=vz_type[iparticle]
        for istar in range(nstars):
            xrel=xp-x_star[istar]
            yrel=yp-y_star[istar]
            zrel=zp-z_star[istar]
            distance=sqrt(xrel**2+yrel**2+zrel**2)
            xrelvel=vxp-vx_star[istar]
            yrelvel=vyp-vy_star[istar]
            zrelvel=vzp-vz_star[istar]
            relvel=sqrt(xrelvel**2+yrelvel**2+zrelvel**2)
            mstar = m_star[istar]
            en = 0.5 * relvel**2 - G * mstar / distance
            costheta = (xrel*xrelvel+yrel*yrelvel+zrel*zrelvel)/distance/relvel
            sintheta = sqrt(1-costheta**2)
            eccen = distance * relvel * sintheta
            #energy criterion
            if (en < enmin and eccen < eccenlimit and distance<distancelimit):
                enmin = en
                owner[iparticle]=istar
                
    return owner


cpdef flag_owner2d(float[:] x_type, float[:] y_type, float[:] vx_type, float[:] vy_type, float[:] m_type,
                   float[:] x_star, float[:] y_star, float[:] vx_star, float[:] vy_star, float[:] m_star,
                   parameters):
    """
    Finds the star owner of each sph particle. Does not check for binaries.
    Returns an array containing the id of the owner for each particle, or -1 if
    it is unbound.
    
    The owner star is simply flagged as the one with the smallest binding energy. There is also a
    limit on the allowed eccentricity and distance from the star, that are read from the parameters dictionary.
    """
    cdef float xp, yp
    cdef float xrel, yrel
    cdef int iparticle, istar
    cdef int nstars, nparticles
    cdef float vxp, vyp
    cdef float distance, relvel
    cdef float xrelvel, yrelvel
    cdef float mstar, en, enmin
    #dimensionless units, so G=1
    cdef float costheta, sintheta, eccen
    cdef float eccenlimit = parameters["eccenlimit"]
    cdef float distancelimit = parameters["distancelimit"]

    nstars=x_type.size
    nparticles=x_star.size

    
    #for each particle, finds the star to which it is mostly bound
    owner = - np.ones(nparticles, dtype='int32')
    for iparticle in range(nparticles):
        enmin=0.
        xp=x_type[iparticle]
        yp=y_type[iparticle]
        vxp=vx_type[iparticle]
        vyp=vy_type[iparticle]
        for istar in range(nstars):
            xrel=xp-x_star[istar]
            yrel=yp-y_star[istar]
            distance=sqrt(xrel**2+yrel**2)
            xrelvel=vxp-vx_star[istar]
            yrelvel=vyp-vy_star[istar]
            relvel=sqrt(xrelvel**2+yrelvel**2)
            mstar = m_star[istar]
            en = 0.5 * relvel**2 - G * mstar / distance
            costheta = (xrel*xrelvel+yrel*yrelvel)/distance/relvel
            sintheta = sqrt(1-costheta**2)
            eccen = distance * relvel * sintheta
            #energy criterion
            if (en < enmin and eccen < eccenlimit and distance<distancelimit):
                enmin = en
                owner[iparticle]=istar
                
    return owner


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
