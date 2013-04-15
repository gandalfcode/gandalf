//=============================================================================
// SphLeapfrogKDK.cpp
// Contains functions for integrating SPH particle positions and velocities 
// using the leapfrog kick-drift-kick (KDK) scheme.
//=============================================================================


#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <math.h>
#include "Sph.h"
#include "SphKernel.h"
#include "SphIntegration.h"
#include "SphParticle.h"
#include "EOS.h"
#include "Debug.h"
using namespace std;




//=============================================================================
//  SphLeapfrogKDK::SphLeapfrogKDK
/// SphLeapfrogKDK class constructor
//=============================================================================
template <int ndim>
SphLeapfrogKDK<ndim>::SphLeapfrogKDK( DOUBLE accel_mult_aux,
			       DOUBLE courant_mult_aux) :
  SphIntegration<ndim>(accel_mult_aux, courant_mult_aux)
{
}



//============================================================================
//  SphLeapfrogKDK::~SphLeapfrog()
/// SphLeapfrogKDK class destructor
//=============================================================================
template <int ndim>
SphLeapfrogKDK<ndim>::~SphLeapfrogKDK()
{
}



//=============================================================================
//  SphLeapfrogKDK::AdvanceParticles
/// Integrate particle positions to 2nd order, and particle velocities to 1st
/// order from the beginning of the step to the current simulation time, i.e. 
/// $r(t+dt) = r(t) + v(t)*dt + 0.5*a(t)*dt^2$, 
/// $v(t+dt) = v(t) + a(t)*dt$.
/// Also set particles at the end of step as 'active' in order to compute 
/// the end-of-step force computation.
//=============================================================================
template <int ndim>
void SphLeapfrogKDK<ndim>::AdvanceParticles
(int n,                             ///< [in] Integer time in block time struct
 int level_step,                    ///< [in] Current block level
 int Nsph,                          ///< [in] No. of SPH particles
 SphParticle<ndim> *sph,            ///< [inout] SPH particle data array
 FLOAT timestep)                    ///< [in] Base timestep value
{
  int i;                            // Particle counter
  int k;                            // Dimension counter
  int nstep;                        // Particle (integer) step size
  FLOAT dt;                         // Timestep since start of step

  debug2("[SphLeapfrogKDK::AdvanceParticles]");

  // Advance positions and velocities of all SPH particles
  // --------------------------------------------------------------------------
#pragma omp parallel for default(shared) private(dt,k,nstep)
  for (i=0; i<Nsph; i++) {

    // Compute time since beginning of current step
    nstep = pow(2,level_step - sph[i].level);
    if (n%nstep == 0) dt = timestep*(FLOAT) nstep;
    else dt = timestep*(FLOAT) (n%nstep);

    // Advance particle positions and velocities
    for (k=0; k<ndim; k++) sph[i].r[k] = sph[i].r0[k] + 
      (sph[i].v0[k] + 0.5*sph[i].a[k]*sph[i].dt)*dt;
    for (k=0; k<ndim; k++) sph[i].v[k] = sph[i].v0[k] + sph[i].a0[k]*dt;

    // Set particle as active at end of step
    if (n%nstep == 0) sph[i].active = true;
    else sph[i].active = false;
  }
  // --------------------------------------------------------------------------

  return;
}
 


//=============================================================================
//  SphLeapfrogKDK::CorrectionTerms
/// Compute velocity integration to second order at the end of the step by 
/// adding a second order correction term.  The full integration becomes
/// v(t+dt) = v(t) + 0.5*(a(t) + a(t+dt))*dt 
//=============================================================================
template <int ndim>
void SphLeapfrogKDK<ndim>::CorrectionTerms
(int n,                             ///< [in] Integer time in block time struct
 int level_step,                    ///< Current block level
 int Nsph,                          ///< No. of SPH particles
 SphParticle<ndim> *sph,            ///< SPH particle data array
 FLOAT timestep)                    ///< Base timestep value
{
  int i;                            // Particle counter
  int k;                            // Dimension counter
  int nstep;                        // Particle (integer) step size

  debug2("[SphLeapfrogKDK::CorrectionTerms]");

  // --------------------------------------------------------------------------
#pragma omp parallel for default(shared) private(k,nstep)
  for (i=0; i<Nsph; i++) {
    nstep = pow(2,level_step - sph[i].level);
    if (n%nstep == 0)
      for (k=0; k<ndim; k++) sph[i].v[k] += 
	     (FLOAT) 0.5*(sph[i].a[k] - sph[i].a0[k])*timestep*(FLOAT) nstep;
  }
  // --------------------------------------------------------------------------

  return;
}



//=============================================================================
//  SphLeapfrogKDK::EndTimestep
/// Record all important SPH particle quantities at the end of the step for  
/// the start of the new timestep.
// ============================================================================
template <int ndim>
void SphLeapfrogKDK<ndim>::EndTimestep
(int n,                             ///< [in] Integer time in block time struct
 int level_step,                    ///< Current block level
 int Nsph,                          ///< No. of SPH particles
 SphParticle<ndim> *sph)            ///< SPH particle data array
{
  int i;                            // Particle counter
  int k;                            // Dimension counter
  int nstep;                        // Particle (integer) step size

  debug2("[SphLeapfrogKDK::EndTimestep]");

  // --------------------------------------------------------------------------
#pragma omp parallel for default(shared) private(k,nstep)
  for (i=0; i<Nsph; i++) {
    nstep = pow(2,level_step - sph[i].level);
    if (n%nstep == 0) {
      for (k=0; k<ndim; k++) sph[i].r0[k] = sph[i].r[k];
      for (k=0; k<ndim; k++) sph[i].v0[k] = sph[i].v[k];
      for (k=0; k<ndim; k++) sph[i].a0[k] = sph[i].a[k];
      //sph[i].active = false;
      sph[i].active = true;
    }
  }
  // --------------------------------------------------------------------------

  return;
}



// Template class instances for each dimensionality value (1, 2 and 3)
template class SphLeapfrogKDK<1>;
template class SphLeapfrogKDK<2>;
template class SphLeapfrogKDK<3>;
