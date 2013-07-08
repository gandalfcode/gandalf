//=============================================================================
//  SphParticle.h
//  Main SPH particle data structure
//=============================================================================


#ifndef _SPH_PARTICLE_H_
#define _SPH_PARTICLE_H_


#include "Precision.h"
#include "Constants.h"


//=============================================================================
//  Structure SphParticle
/// \brief  SPH particle data structure.
/// \author D. A. Hubber, G. Rosotti
/// \date   03/04/2013
//=============================================================================
template <int ndim>
struct SphParticle {
  bool active;                      ///< Flag if active (i.e. recompute step)
  bool potmin;                      ///< Is particle at a potential minima?
  int iorig;                        ///< Original particle i.d.
  int itype;                        ///< SPH particle type
  int level;                        ///< Current timestep level of particle
  int nstep;                        ///< Integer step-size of particle
  FLOAT r[ndim];                    ///< Position
  FLOAT v[ndim];                    ///< Velocity
  FLOAT a[ndim];                    ///< Total acceleration
  FLOAT r0[ndim];                   ///< Position at beginning of step
  FLOAT v0[ndim];                   ///< Velocity at beginning of step
  FLOAT a0[ndim];                   ///< Acceleration at beginning of step
  FLOAT agrav[ndim];                ///< Gravitational acceleration
  FLOAT u;                          ///< Specific internal energy
  FLOAT u0;                         ///< u at beginning of step
  FLOAT dudt;                       ///< Compressional heating rate
  FLOAT dudt0;                      ///< dudt at beginning of step
  FLOAT m;                          ///< Particle mass
  FLOAT h;                          ///< SPH smoothing length
  FLOAT invh;                       ///< 1 / h
  FLOAT hfactor;                    ///< invh^(ndim + 1)
  FLOAT rho;                        ///< SPH density
  FLOAT invrho;                     ///< 1 / rho
  FLOAT press;                      ///< Thermal pressure
  FLOAT pfactor;                    ///< Pressure factor in SPH EOM
  FLOAT div_v;                      ///< Velocity divergence
  FLOAT invomega;                   ///< grad-h omega/f correction term
  FLOAT zeta;                       ///< grad-h gravity correction term
  FLOAT chi;                        ///< grad-h star-gravity correction term
  FLOAT q;                          ///< Internal energy density
  FLOAT invq;                       ///< 1 / q
  FLOAT sound;                      ///< Sound speed
  FLOAT gpot;                       ///< Gravitational potential
  FLOAT gpe;                        ///< Gravitational potential energy
  FLOAT gradP[ndim];                ///< Pressure gradient
  FLOAT gradrho[ndim];              ///< Density gradient
  FLOAT gradv[ndim][ndim];          ///< Velocity gradient matrix
  FLOAT rhomax;                     ///< Maximum neighbour density
  FLOAT rhomin;                     ///< Minimum neighbour density
  FLOAT pressmax;                   ///< Maximum neighbour pressure
  FLOAT pressmin;                   ///< Minimum neighbour pressure
  FLOAT vmax[ndim];                 ///< Maximum neighbour velocity
  FLOAT vmin[ndim];                 ///< Minimum neighbour velocity
  DOUBLE dt;                        ///< Particle timestep


  // SPH particle constructor to initialise all values
  // --------------------------------------------------------------------------
  SphParticle()
  {
    active = false;
    potmin = false;
    iorig = -1;
    itype = -1;
    level = 0;
    for (int k=0; k<ndim; k++) r[k] = (FLOAT) 0.0;
    for (int k=0; k<ndim; k++) v[k] = (FLOAT) 0.0;
    for (int k=0; k<ndim; k++) a[k] = (FLOAT) 0.0;
    for (int k=0; k<ndim; k++) r0[k] = (FLOAT) 0.0;
    for (int k=0; k<ndim; k++) v0[k] = (FLOAT) 0.0;
    for (int k=0; k<ndim; k++) a0[k] = (FLOAT) 0.0;
    for (int k=0; k<ndim; k++) agrav[k] = (FLOAT) 0.0;
    u = (FLOAT) 0.0;
    u0 = (FLOAT) 0.0;
    dudt = (FLOAT) 0.0;
    dudt0 = (FLOAT) 0.0;
    m = (FLOAT) 0.0;
    h = (FLOAT) 0.0;
    invh = (FLOAT) 0.0;
    hfactor = (FLOAT) 0.0;
    rho = (FLOAT) 0.0;
    invrho = (FLOAT) 0.0;
    press = (FLOAT) 0.0;
    pfactor = (FLOAT) 0.0;
    invomega = (FLOAT) 0.0;
    zeta = (FLOAT) 0.0;
    chi = (FLOAT) 0.0;
    q = (FLOAT) 0.0;
    invq = (FLOAT) 0.0;
    gpot = (FLOAT) 0.0;
    gpe = (FLOAT) 0.0;
    sound = (FLOAT) 0.0;
    for (int k=0; k<ndim; k++) gradP[k] = (FLOAT) 0.0;
    for (int k=0; k<ndim; k++) gradrho[k] = (FLOAT) 0.0;
    for (int k=0; k<ndim; k++)
      for (int kk=0; kk<ndim; kk++) gradv[k][kk] = (FLOAT) 0.0;
    dt = (DOUBLE) 0.0;
    rhomin = (FLOAT) 0.0;
    rhomax = (FLOAT) 0.0;
    pressmin = (FLOAT) 0.0;
    pressmax = (FLOAT) 0.0;
    for (int k=0; k<ndim; k++) vmax[k] = 0.0;
    for (int k=0; k<ndim; k++) vmin[k] = 0.0;
  }

};
#endif
