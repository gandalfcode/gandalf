// ============================================================================
// NbodyLeapfrogKDK.cpp
// Contains functions for integrating star particle positions and velocities 
// using the leapfrog kick-drift-kick (KDK) scheme.
// ============================================================================


#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <math.h>
#include "Dimensions.h"
#include "Sph.h"
#include "NbodyIntegration.h"
#include "StarParticle.h"
#include "Debug.h"
using namespace std;



// ============================================================================
// NbodyLeapfrogKDK::NbodyLeapfrogKDK()
// ============================================================================
NbodyLeapfrogKDK::NbodyLeapfrogKDK(int ndimaux, int vdimaux, 
				   DOUBLE nbody_mult_aux) : 
  NbodyIntegration(ndimaux, vdimaux, nbody_mult_aux)
{
}



// ============================================================================
// NbodyLeapfrogKDK::~NbodyLeapfrog()
// ============================================================================
NbodyLeapfrogKDK::~NbodyLeapfrogKDK()
{
}



// ============================================================================
// NbodyLeapfrogKDK::AdvanceParticles
// Integrate particle positions to 2nd order, and particle velocities to 1st
// order from the beginning of the step to the current simulation time, i.e. 
// r(t+dt) = r(t) + v(t)*dt + 0.5*a(t)*dt^2, 
// v(t+dt) = v(t) + a(t)*dt.
// Also set particles at the end of step as 'active' in order to compute 
// the end-of-step force computation.
// ============================================================================
void NbodyLeapfrogKDK::AdvanceParticles(int n, int level_step, int Nsystem,
					StarParticle *sph, Double timestep)
{
  int i;                                // Particle counter
  int k;                                // Dimension counter
  int nstep;                            // Particle (integer) step size
  DOUBLE dt;                            // Timestep since start of step

  debug2("[NbodyLeapfrogKDK::AdvanceParticles]");

  for (i=0; i<Nsystem; i++) {
    nstep = pow(2,level_step - system[i].level);
    if (n%nstep == 0) dt = timestep*(DOUBLE) nstep;
    else dt = timestep*(DOUBLE) (n%nstep);
    for (k=0; k<ndim; k++) system[i].r[k] = system[i].r0[k] + 
			     system[i].v0[k]*dt + 0.5*sph[i].a0[k]*dt*dt;
    for (k=0; k<vdim; k++)
      system[i].v[k] = system[i].v0[k] + system[i].a0[k]*dt;
    if (n%nstep == 0) system[i].active = true;
  }

  return;
}
 


// ============================================================================
// NbodyLeapfrogKDK::CorrectionTerms
// Compute velocity integration to second order at the end of the step by 
// adding a second order correction term.  The full integration becomes
// v(t+dt) = v(t) + 0.5*(a(t) + a(t+dt))*dt 
// ============================================================================
void NbodyLeapfrogKDK::CorrectionTerms(int n, int level_step, int Nsystem,
				       StarParticle *system, DOUBLE timestep)
{
  int i;                                // Particle counter
  int k;                                // Dimension counter
  int nstep;                            // Particle (integer) step size

  debug2("[NbodyLeapfrogKDK::CorrectionTerms]");

  for (i=0; i<Nsystem; i++) {
    nstep = pow(2,level_step - system[i].level);
    if (n%nstep == 0)
      for (k=0; k<ndim; k++) system[i].v[k] += 
			       0.5*(system[i].a[k] - system[i].a0[k])*
			       timestep*(DOUBLE) nstep;
  }

  return;
}



// ============================================================================
// NbodyLeapfrogKDK::EndTimestep
// Record all important star particle quantities at the end of the step 
// for the start of the new timestep.
// ============================================================================
void NbodyLeapfrogKDK::EndTimestep(int n, int level_step, 
				   int Nsystem, StarParticle *system)
{
  int i;                                // Particle counter
  int k;                                // Dimension counter
  int nstep;                            // Particle (integer) step size

  debug2("[NbodyLeapfrogKDK::EndTimestep]");

  for (i=0; i<Nsystem; i++) {
    nstep = pow(2,level_step - system[i].level);
    if (n%nstep == 0) {
      for (k=0; k<ndim; k++) system[i].r0[k] = system[i].r[k];
      for (k=0; k<ndim; k++) system[i].v0[k] = system[i].v[k];
      for (k=0; k<ndim; k++) system[i].a0[k] = system[i].a[k];
      system[i].active = false;
    }
  }

  return;
}
