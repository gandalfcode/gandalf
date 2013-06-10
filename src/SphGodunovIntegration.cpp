//=============================================================================
//  SphGodunovIntegration.cpp
//  Contains functions for integrating SPH particle positions and velocities 
//  using the conservative integration scheme for Godunov SPH described by 
//  Inutsuka (2002).
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
//  SphGodunovIntegration::SphGodunovIntegration
/// SphGodunovIntegration class constructor
//=============================================================================
template <int ndim>
SphGodunovIntegration<ndim>::SphGodunovIntegration(
			       DOUBLE accel_mult_aux, 
			       DOUBLE courant_mult_aux) :
  SphIntegration<ndim>(accel_mult_aux, courant_mult_aux)
{
}



//=============================================================================
//  SphGodunovIntegration::~SphGodunovIntegration()
/// SphGodunovIntegration class destructor
//=============================================================================
template <int ndim>
SphGodunovIntegration<ndim>::~SphGodunovIntegration()
{
}



//=============================================================================
//  SphGodunovIntegration::AdvanceParticles
/// Integrate particle positions to 2nd order, and particle velocities to 1st
/// order from the beginning of the step to the current simulation time, i.e. 
/// $r(t+dt) = r(t) + v(t)*dt + 0.5*a(t)*dt^2$, 
/// $v(t+dt) = v(t) + a(t)*dt$.
/// Also set particles at the end of step as 'active' in order to compute 
/// the end-of-step force computation.
//=============================================================================
template <int ndim>
void SphGodunovIntegration<ndim>::AdvanceParticles
(int n,                             ///< [in] Current integer time
 int Nsph,                          ///< [in] No. of SPH particles
 SphParticle<ndim> *sphdata,        ///< [inout] SPH particle data array
 FLOAT timestep)                    ///< [in] Minimum timestep level
{
  int i;                            // Particle counter
  int k;                            // Dimension counter
  int nstep;                        // Particle (integer) step size
  FLOAT dt;                         // Timestep since start of step

  debug2("[SphGodunovIntegration::AdvanceParticles]");

  // --------------------------------------------------------------------------
#pragma omp parallel for default(shared) private(dt,k,nstep)
  for (i=0; i<Nsph; i++) {

    // Compute time since beginning of current step
    nstep = sphdata[i].nstep;
    if (n%nstep == 0) dt = timestep*(FLOAT) nstep;
    else dt = timestep*(FLOAT) (n%nstep);

    // Advance particle positions and velocities
    for (k=0; k<ndim; k++) sphdata[i].r[k] = sphdata[i].r0[k] + 
      (sphdata[i].v0[k] + 0.5*sphdata[i].a[k]*sphdata[i].dt)*dt;
    for (k=0; k<vdim; k++) sphdata[i].v[k] = sphdata[i].v0[k] + 
      sphdata[i].a0[k]*dt;

    // Set particle as active at end of step
    if (n%nstep == 0) sphdata[i].active = true;
    else sphdata[i].active = false;
  }
  // --------------------------------------------------------------------------

  return;
}
 


//=============================================================================
//  SphGodunovIntegration::CorrectionTerms
/// Empty definition (No correction terms to apply).
//=============================================================================
template <int ndim>
void SphGodunovIntegration<ndim>::CorrectionTerms
(int n, int Nsph, SphParticle<ndim> *sph, FLOAT timestep)
{
  return;
}



//=============================================================================
//  SphGodunovIntegration::EndTimestep
/// Record all important SPH particle quantities at the end of the step for 
/// the start of the new timestep.
//=============================================================================
template <int ndim>
void SphGodunovIntegration<ndim>::EndTimestep
(int n,                             ///< [in] Integer time in block time struct
 int Nsph,                          ///< [in] No. of SPH particles
 SphParticle<ndim> *sphdata)        ///< [inout] SPH particle data array
{
  int i;                            // Particle counter
  int k;                            // Dimension counter
  int nstep;                        // Particle (integer) step size

  debug2("[SphGodunovIntegration::EndTimestep]");

  // --------------------------------------------------------------------------
#pragma omp parallel for default(shared) private(k,nstep)
  for (i=0; i<Nsph; i++) {
    nstep = sphdata[i].nstep;
    if (n%nstep == 0) {
      for (k=0; k<ndim; k++) sphdata[i].r0[k] = sphdata[i].r[k];
      for (k=0; k<ndim; k++) sphdata[i].v0[k] = sphdata[i].v[k];
      for (k=0; k<ndim; k++) sphdata[i].a0[k] = sphdata[i].a[k];
      //sphdata[i].active = false;
      sphdata[i].active = true;
    }
  }
  // --------------------------------------------------------------------------

  return;
}



//=============================================================================
//  SphIntegration::Timestep
/// Default timestep size for SPH particles.  Takes the minimum of : 
/// (i)  const*h/(sound_speed + h*|div_v|)    (Courant condition)
/// (ii) const*sqrt(h/|a|)                    (Acceleration condition)
//=============================================================================
template <int ndim>
DOUBLE SphGodunovIntegration<ndim>::Timestep
(SphParticle<ndim> &part,           ///< Reference to SPH particle
 int hydro_forces)                  ///< Computing hydro forces or not
{
  DOUBLE timestep;                  // Variable to record/compare timesteps
  //DOUBLE amag;

  // Courant condition
  timestep = this->courant_mult*part.h/(part.sound + small_number_dp);

  // Local convergence/divergence condition
  timestep = min(timestep,this->courant_mult/
                 (fabs(part.div_v) + small_number_dp));

  //Acceleration condition
  //amag = sqrt(DotProduct(part.a,part.a,ndim));
  //timestep = min(timestep, accel_mult*sqrt(part.h/(amag + small_number_dp)));

  return timestep;
}



template class SphGodunovIntegration<1>;
template class SphGodunovIntegration<2>;
template class SphGodunovIntegration<3>;
