// ============================================================================
// SphLeapfrogKDK.cpp
// Contains functions for integrating SPH particle positions and velocities 
// using the leapfrog kick-drift-kick (KDK) scheme.
// ============================================================================


#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <math.h>
#include "Dimensions.h"
#include "Sph.h"
#include "SphKernel.h"
#include "SphIntegration.h"
#include "SphParticle.h"
#include "EOS.h"
#include "Debug.h"
using namespace std;



// ============================================================================
// SphLeapfrogKDK::SphLeapfrogKDK()
// ============================================================================
SphLeapfrogKDK::SphLeapfrogKDK(int ndimaux, int vdimaux, 
			       DOUBLE accel_mult_aux, 
			       DOUBLE courant_mult_aux) :
  SphIntegration(ndimaux, vdimaux, accel_mult_aux, courant_mult_aux)
{
}



// ============================================================================
// SphLeapfrogKDK::~SphLeapfrog()
// ============================================================================
SphLeapfrogKDK::~SphLeapfrogKDK()
{
}



// ============================================================================
// SphLeapfrogKDK::AdvanceParticles
// Integrate particle positions to 2nd order, and particle velocities to 1st
// order from the beginning of the step to the current simulation time, i.e. 
// r(t+dt) = r(t) + v(t)*dt + 0.5*a(t)*dt^2, 
// v(t+dt) = v(t) + a(t)*dt.
// Also set particles at the end of step as 'active' in order to compute 
// the end-of-step force computation.
// ============================================================================
void SphLeapfrogKDK::AdvanceParticles(int n, int level_step, int Nsph,
				      SphParticle *sph, FLOAT timestep)
{
  int i;                                // Particle counter
  int k;                                // Dimension counter
  int nstep;                            // Particle (integer) step size
  FLOAT dt;                             // Timestep since start of step

  debug2("[SphLeapfrogKDK::AdvanceParticles]");

  for (i=0; i<Nsph; i++) {
    nstep = pow(2,level_step - sph[i].level);
    if (n%nstep == 0) dt = timestep*(FLOAT) nstep;
    else dt = timestep*(FLOAT) (n%nstep);
    for (k=0; k<ndim; k++) sph[i].r[k] = sph[i].r0[k] + sph[i].v0[k]*dt
      + (FLOAT) 0.5*sph[i].a0[k]*dt*dt;
    for (k=0; k<vdim; k++) sph[i].v[k] = sph[i].v0[k] + sph[i].a0[k]*dt;
    if (n%nstep == 0) sph[i].active = true;
  }

  return;
}
 


// ============================================================================
// SphLeapfrogKDK::CorrectionTerms
// Compute velocity integration to second order at the end of the step by 
// adding a second order correction term.  The full integration becomes
// v(t+dt) = v(t) + 0.5*(a(t) + a(t+dt))*dt 
// ============================================================================
void SphLeapfrogKDK::CorrectionTerms(int n, int level_step, int Nsph,
				     SphParticle *sph, FLOAT timestep)
{
  int i;                                // Particle counter
  int k;                                // Dimension counter
  int nstep;                            // Particle (integer) step size

  debug2("[SphLeapfrogKDK::CorrectionTerms]");

  for (i=0; i<Nsph; i++) {
    nstep = pow(2,level_step - sph[i].level);
    if (n%nstep == 0)
      for (k=0; k<ndim; k++) sph[i].v[k] += 
	(FLOAT) 0.5*(sph[i].a[k] - sph[i].a0[k])*timestep*(FLOAT) nstep;
  }

  return;
}



// ============================================================================
// SphLeapfrogKDK::EndTimestep
// Record all important SPH particle quantities at the end of the step for the 
// start of the new timestep.
// ============================================================================
void SphLeapfrogKDK::EndTimestep(int n, int level_step, 
				 int Nsph, SphParticle *sph)
{
  int i;                                // Particle counter
  int k;                                // Dimension counter
  int nstep;                            // Particle (integer) step size

  debug2("[SphLeapfrogKDK::EndTimestep]");

  for (i=0; i<Nsph; i++) {
    nstep = pow(2,level_step - sph[i].level);
    if (n%nstep == 0) {
      for (k=0; k<ndim; k++) sph[i].r0[k] = sph[i].r[k];
      for (k=0; k<ndim; k++) sph[i].v0[k] = sph[i].v[k];
      for (k=0; k<ndim; k++) sph[i].a0[k] = sph[i].a[k];
      sph[i].active = false;
    }
  }

  return;
}
