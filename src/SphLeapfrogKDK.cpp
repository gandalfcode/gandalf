//=============================================================================
//  SphLeapfrogKDK.cpp
//  Contains functions for integrating SPH particle positions and velocities 
//  using the leapfrog kick-drift-kick (KDK) scheme.
//
//  This file is part of GANDALF :
//  Graphical Astrophysics code for N-body Dynamics And Lagrangian Fluids
//  https://github.com/gandalfcode/gandalf
//  Contact : gandalfcode@gmail.com
//
//  Copyright (C) 2013  D. A. Hubber, G. Rosotti
//
//  GANDALF is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 2 of the License, or
//  (at your option) any later version.
//
//  GANDALF is distributed in the hope that it will be useful, but
//  WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  General Public License (http://www.gnu.org/licenses) for more details.
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
SphLeapfrogKDK<ndim>::SphLeapfrogKDK(DOUBLE accel_mult_aux, 
                                     DOUBLE courant_mult_aux,
                                     DOUBLE energy_mult_aux,
				     eosenum gas_eos_aux,
				     tdaviscenum tdavisc_aux) :
  SphIntegration<ndim>(accel_mult_aux, courant_mult_aux, energy_mult_aux,
		       gas_eos_aux, tdavisc_aux)
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
 FLOAT timestep,                    ///< [in] Base timestep value
 Sph<ndim> *sph)                    ///< [inout] Pointer to SPH object
{
  int dn;                           // Integer time since beginning of step
  int i;                            // Particle counter
  int k;                            // Dimension counter
  int nstep;                        // Particle (integer) step size
  FLOAT dt;                         // Timestep since start of step
  SphParticle<ndim> *part;          // SPH particle pointer

  debug2("[SphLeapfrogKDK::AdvanceParticles]");
  timing->StartTimingSection("LFKDK_ADVANCE_PARTICLES",2);

  // Advance positions and velocities of all SPH particles
  //---------------------------------------------------------------------------
#pragma omp parallel for default(none) private(dn,dt,i,k,nstep,part)\
  shared(n,sph,timestep)
  for (i=0; i<sph->Nsph; i++) {

    // Compute time since beginning of current step
    nstep = sph->sphintdata[i].nstep;
    dn = n - sph->sphintdata[i].nlast;
    dt = timestep*(FLOAT) dn;
    part = sph->sphintdata[i].part;

    // Advance particle positions and velocities
    for (k=0; k<ndim; k++) part->r[k] = sph->sphintdata[i].r0[k] +
      sph->sphintdata[i].v0[k]*dt + 0.5*sph->sphintdata[i].a0[k]*dt*dt;
    for (k=0; k<ndim; k++) part->v[k] =
      sph->sphintdata[i].v0[k] + sph->sphintdata[i].a0[k]*dt;

    // Integrate time-dependent viscosity
    if (tdavisc != notdav) part->alpha += part->dalphadt*timestep;

    // Integrate explicit energy equation
    if (gas_eos == energy_eqn) 
      part->u = sph->sphintdata[i].u0 + sph->sphintdata[i].dudt0*dt;

    // Set particle as active at end of step
    if (dn == nstep) part->active = true;
    else part->active = false;
  }
  //---------------------------------------------------------------------------

  timing->EndTimingSection("LFKDK_ADVANCE_PARTICLES");

  return;
}
 


//=============================================================================
//  SphLeapfrogKDK::CorrectionTerms
/// Compute velocity integration to second order at the end of the step by 
/// adding a second order correction term.  The full integration becomes
/// $v(t+dt) = v(t) + 0.5*(a(t) + a(t+dt))*dt$. 
//=============================================================================
template <int ndim>
void SphLeapfrogKDK<ndim>::CorrectionTerms
(int n,                             ///< [in] Integer time in block time struct
 FLOAT timestep,                    ///< [in] Base timestep value
 Sph<ndim> *sph)                    ///< [inout] Pointer to SPH object
{
  int dn;                           // Integer time since beginning of step
  int i;                            // Particle counter
  int k;                            // Dimension counter
  int nstep;                        // Particle (integer) step size
  SphParticle<ndim> *part;          // SPH particle pointer

  debug2("[SphLeapfrogKDK::CorrectionTerms]");
  timing->StartTimingSection("LFKDK_CORRECTION_TERMS",2);

  //---------------------------------------------------------------------------
#pragma omp parallel for default(none) private(dn,i,k,nstep,part)\
  shared(n,sph,timestep)
  for (i=0; i<sph->Nsph; i++) {
    dn = n - sph->sphintdata[i].nlast;
    nstep = sph->sphintdata[i].nstep;
    part = sph->sphintdata[i].part;
    if (dn == nstep)
      for (k=0; k<ndim; k++) part->v[k] += timestep*(FLOAT) nstep*
        (FLOAT) 0.5*(part->a[k] - sph->sphintdata[i].a0[k]);
    if (dn == nstep && gas_eos == energy_eqn) part->u +=
      0.5*(part->dudt - sph->sphintdata[i].dudt0)*timestep*(FLOAT) nstep;

  }
  //---------------------------------------------------------------------------

  timing->EndTimingSection("LFKDK_CORRECTION_TERMS");

  return;
}



//=============================================================================
//  SphLeapfrogKDK::EndTimestep
/// Record all important SPH particle quantities at the end of the step for  
/// the start of the new timestep.
//=============================================================================
template <int ndim>
void SphLeapfrogKDK<ndim>::EndTimestep
(int n,                             ///< [in] Integer time in block time struct
 FLOAT timestep,                    ///< [in] Base timestep value
 Sph<ndim> *sph)                    ///< [inout] Pointer to SPH object
{
  int dn;                           // Integer time since beginning of step
  int i;                            // Particle counter
  int k;                            // Dimension counter
  int nstep;                        // Particle (integer) step size
  SphParticle<ndim> *part;          // SPH particle pointer

  debug2("[SphLeapfrogKDK::EndTimestep]");
  timing->StartTimingSection("LFKDK_END_TIMESTEP",2);

  //---------------------------------------------------------------------------
#pragma omp parallel for default(none) private(dn,i,k,nstep,part) shared(n,sph,timestep)
  for (i=0; i<sph->Nsph; i++) {
    dn = n - sph->sphintdata[i].nlast;
    nstep = sph->sphintdata[i].nstep;
    part = sph->sphintdata[i].part;

    if (dn == nstep) {
      for (k=0; k<ndim; k++) part->v[k] += timestep*(FLOAT) nstep*
        (FLOAT) 0.5*(part->a[k] - sph->sphintdata[i].a0[k]);
      for (k=0; k<ndim; k++) sph->sphintdata[i].r0[k] = part->r[k];
      for (k=0; k<ndim; k++) sph->sphintdata[i].v0[k] = part->v[k];
      for (k=0; k<ndim; k++) sph->sphintdata[i].a0[k] = part->a[k];
      if (gas_eos == energy_eqn) {
	part->u += 0.5*(part->dudt - sph->sphintdata[i].dudt0)*
	  timestep*(FLOAT) nstep;
	sph->sphintdata[i].u0 = part->u;
	sph->sphintdata[i].dudt0 = part->dudt;
      }
      sph->sphintdata[i].nlast = n;
      part->active = false;
    }
  }
  //---------------------------------------------------------------------------

  timing->EndTimingSection("LFKDK_END_TIMESTEP");

  return;
}



//=============================================================================
//  SphLeapfrogKDK::CheckTimesteps
/// Check through all SPH particles to see if the maximum neighbour timestep 
/// exceeds the maximum allowed difference (level_diff_max).  If so, then 
/// check if we can prematurely finish the current timestep and move the 
/// to a lower timestep level.  Returns the number of particles whose timesteps
/// have been reduced.
//=============================================================================
template <int ndim>
int SphLeapfrogKDK<ndim>::CheckTimesteps
(int level_diff_max,                ///< [in] Max. allowed SPH neib dt diff
 int n,                             ///< [in] Integer time in block time struct
 Sph<ndim> *sph)                    ///< [inout] Pointer to SPH object
{
  int activecount = 0;              // No. of ptcls with new active timesteps
  int dn;                           // Integer time since beginning of step
  int i;                            // Particle counter
  int level_new;                    // New level of particle
  int nnewstep;                     // New step size of particle
  SphParticle<ndim> *part;          // SPH particle pointer

  debug2("[SphLeapfrogKDK::CheckTimesteps]");
  timing->StartTimingSection("LFKDK_CHECK_TIMESTEPS",2);

  //---------------------------------------------------------------------------
#pragma omp parallel for default(none) private(dn,level_new,nnewstep,part)\
  shared(level_diff_max,n,sph) reduction(+:activecount)
  for (i=0; i<sph->Nsph; i++) {
    dn = n - sph->sphintdata[i].nlast;
    part = sph->sphintdata[i].part;
    if (dn == sph->sphintdata[i].nstep) continue;

    // Check if neighbour timesteps are too small.  If so, then reduce 
    // timestep if possible
    if (part->levelneib - part->level > level_diff_max) {
      level_new = part->levelneib - level_diff_max;
      nnewstep = sph->sphintdata[i].nstep/pow(2,level_new - part->level);

      // If new level is correctly synchronised, then change all quantities
      if (dn%nnewstep == 0) {
        part->level = level_new;
        if (dn > 0) sph->sphintdata[i].nstep = dn;
        part->active = true;
        activecount++;
      }
    }
  }
  //---------------------------------------------------------------------------

  timing->EndTimingSection("LFKDK_CHECK_TIMESTEPS");

  return activecount;
}



// Template class instances for each dimensionality value (1, 2 and 3)
template class SphLeapfrogKDK<1>;
template class SphLeapfrogKDK<2>;
template class SphLeapfrogKDK<3>;
