//=============================================================================
//  NbodyParticle.h
//  Main star particle data structure
//=============================================================================


#ifndef _NBODY_PARTICLE_H_
#define _NBODY_PARTICLE_H_


#include "Precision.h"
#include "Constants.h"


//=============================================================================
//  Class NbodyParticle
/// \brief  N-body particle data structure
/// \author D. A. Hubber
/// \date   15/04/2013
//=============================================================================
template <int ndim>
class NbodyParticle
{
 public:

  bool active;                      ///< Flag if active (i.e. recompute step)
  int level;                        ///< Current timestep level    
  DOUBLE r[ndim];                   ///< Position
  DOUBLE v[ndim];                   ///< Velocity
  DOUBLE a[ndim];                   ///< Acceleration
  DOUBLE adot[ndim];                ///< Time derivative of acceleration (jerk)
  DOUBLE adot2[ndim];               ///< 2nd time derivative of acceleration
  DOUBLE adot3[ndim];               ///< 3rd time derivative of acceleration
  DOUBLE r0[ndim];                  ///< Position at beginning of step
  DOUBLE v0[ndim];                  ///< Velocity at beginning of step
  DOUBLE a0[ndim];                  ///< Acceleration at beginning of step
  DOUBLE adot0[ndim];               ///< Jerk at beginning of step
  DOUBLE m;                         ///< Star mass
  DOUBLE h;                         ///< Smoothing length
  DOUBLE invh;                      ///< 1 / h
  DOUBLE hfactor;                   ///< invh^(ndim + 1)
  DOUBLE gpot;                      ///< Gravitational potential


  // Star particle constructor to initialise all values
  // --------------------------------------------------------------------------
  NbodyParticle()
  {
    active = false;
    for (int k=0; k<ndim; k++) r[k] = 0.0;
    for (int k=0; k<ndim; k++) v[k] = 0.0;
    for (int k=0; k<ndim; k++) a[k] = 0.0;
    for (int k=0; k<ndim; k++) adot[k] = 0.0;
    for (int k=0; k<ndim; k++) r0[k] = 0.0;
    for (int k=0; k<ndim; k++) v0[k] = 0.0;
    for (int k=0; k<ndim; k++) a0[k] = 0.0;
    for (int k=0; k<ndim; k++) adot0[k] = 0.0;
    m = 0;
    h = 0;
    invh = 0.0;
    hfactor = 0.0;
  } 

};
#endif
