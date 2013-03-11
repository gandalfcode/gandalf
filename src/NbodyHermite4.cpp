// ============================================================================
// NbodyHermite4.cpp
// Contains functions for integrating star particle positions and velocities 
// using the 4th-order Hermite scheme (Makino & Aarseth 1992).
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
// NbodyHermite4::NbodyHermite4()
// ============================================================================
NbodyHermite4::NbodyHermite4(int ndimaux, int vdimaux, DOUBLE nbody_mult_aux) :
  NbodyIntegration(ndimaux, vdimaux, nbody_mult_aux)
{
}



// ============================================================================
// NbodyHermite4::~NbodyLeapfrog()
// ============================================================================
NbodyHermite4::~NbodyHermite4()
{
}



// ============================================================================
// NbodyHermite4::AdvanceParticles
// Integrate particle positions to 3nd order, and particle velocities to 2nd
// order from the beginning of the step to the current simulation time, i.e. 
// r(t+dt) = r(t) + v(t)*dt + a(t)*dt^2/2 + adot(t)*dt^3/6, 
// v(t+dt) = v(t) + a(t)*dt + adot(t)*dt^2/2.
// Also set particles at the end of step as 'active' in order to compute 
// the end-of-step force computation.
// ============================================================================
void NbodyHermite4::AdvanceParticles(int n, int level_step, int Nsystem,
				     StarParticle *system, Double timestep)
{
  int i;                                // Particle counter
  int k;                                // Dimension counter
  int nstep;                            // Particle (integer) step size
  DOUBLE dt;                            // Timestep since start of step

  debug2("[NbodyHermite4::AdvanceParticles]");

  // Advance all systems
  // --------------------------------------------------------------------------
  for (i=0; i<Nsystem; i++) {

    // Compute time since beginning of step
    nstep = pow(2,level_step - system[i].level);
    if (n%nstep == 0) dt = timestep*(DOUBLE) nstep;
    else dt = timestep*(DOUBLE) (n%nstep);

    // Advance positions to third order and velocities to second order
    for (k=0; k<ndim; k++) system[i].r[k] = system[i].r0[k] + 
      system[i].v0[k]*dt + 0.5*system[i].a0[k]*dt*dt + 
      onesixth*system[i].adot0[k]*dt*dt*dt;
    for (k=0; k<vdim; k++) system[i].v[k] = system[i].v0[k] + 
      system[i].a0[k]*dt + 0.5*system[i].adot0[k]*dt*dt;

    // If at end of step, set system particle as active
    if (n%nstep == 0) system[i].active = true;
  }
  // --------------------------------------------------------------------------

  return;
}
 


// ============================================================================
// NbodyHermite4::CorrectionTerms
// Compute 2nd and 3rd time derivatives of the acceleration using Hermite 
// interpolation.  Finally correct positions to 5th order and velocities to 
// 4th order using higher-order derivatives.
// ============================================================================
void NbodyHermite4::CorrectionTerms(int n, int level_step, int Nsystem,
				    StarParticle *system, DOUBLE timestep)
{
  int i;                                // Particle counter
  int k;                                // Dimension counter
  int nstep;                            // Particle (integer) step size

  debug2("[NbodyHermite4::CorrectionTerms]");

  // Loop over all system particles
  // --------------------------------------------------------------------------
  for (i=0; i<Nsystem; i++) {
    nstep = pow(2,level_step - system[i].level);

    // If at end of step, compute correction terms
    if (n%nstep == 0) {
      for (k=0; k<ndim; k++) {

	system[i].adot2[k] = 0.0;
	system[i].adot3[k] = 0.0;

	system[i].r[k] += 0.0;
	system[i].v[k] += 0.0;

      }
    }

  }
  // --------------------------------------------------------------------------

  return;
}



// ============================================================================
// NbodyHermite4::EndTimestep
// Record all important star particle quantities at the end of the step 
// for the start of the new timestep.
// ============================================================================
void NbodyHermite4::EndTimestep(int n, int level_step, 
				int Nsystem, StarParticle *system)
{
  int i;                                // Particle counter
  int k;                                // Dimension counter
  int nstep;                            // Particle (integer) step size

  debug2("[NbodyHermite4::EndTimestep]");

  // Loop over all system particles
  // --------------------------------------------------------------------------
  for (i=0; i<Nsystem; i++) {
    nstep = pow(2,level_step - system[i].level);

    // If at end of the current step, set quantites for start of new step
    if (n%nstep == 0) {
      for (k=0; k<ndim; k++) system[i].r0[k] = system[i].r[k];
      for (k=0; k<ndim; k++) system[i].v0[k] = system[i].v[k];
      for (k=0; k<ndim; k++) system[i].a0[k] = system[i].a[k];
      for (k=0; k<ndim; k++) system[i].adot0[k] = system[i].adot[k];
      system[i].active = false;
    }

  }
  // --------------------------------------------------------------------------

  return;
}
