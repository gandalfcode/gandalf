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
template <int ndim, template <int> class ParticleType>
SphLeapfrogKDK<ndim, ParticleType>::SphLeapfrogKDK(DOUBLE accel_mult_aux,
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
template <int ndim, template <int> class ParticleType>
SphLeapfrogKDK<ndim, ParticleType>::~SphLeapfrogKDK()
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
template <int ndim, template <int> class ParticleType>
void SphLeapfrogKDK<ndim, ParticleType >::AdvanceParticles
(int n,                             ///< [in] Integer time in block time struct
 FLOAT timestep,                    ///< [in] Base timestep value
 int Npart,                         ///< [in] Number of particles
 SphParticle<ndim>* sph_gen)        ///< [inout] Pointer to SPH particle array
{
  int dn;                           // Integer time since beginning of step
  int i;                            // Particle counter
  int k;                            // Dimension counter
  int nstep;                        // Particle (integer) step size
  FLOAT dt;                         // Timestep since start of step
  ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sph_gen);

  debug2("[SphLeapfrogKDK::AdvanceParticles]");
  timing->StartTimingSection("SPH_LFKDK_ADVANCE_PARTICLES",2);


  // Advance positions and velocities of all SPH particles
  //---------------------------------------------------------------------------
#pragma omp parallel for default(none) private(dn,dt,i,k,nstep)\
  shared(cout,n,Npart,sphdata,timestep)
  for (i=0; i<Npart; i++) {
    SphParticle<ndim>& part = sphdata[i];
    if (part.itype == dead) continue;

    // Compute time since beginning of current step
    nstep = part.nstep;
    dn = n - part.nlast;
    dt = timestep*(FLOAT) dn;

    // Advance particle positions and velocities
    for (k=0; k<ndim; k++) part.r[k] = part.r0[k] +
      part.v0[k]*dt + 0.5*part.a0[k]*dt*dt;
    for (k=0; k<ndim; k++) part.v[k] =
      part.v0[k] + part.a0[k]*dt;

    // Integrate time-dependent viscosity
    if (tdavisc != notdav) part.alpha += part.dalphadt*timestep;

    // Integrate explicit energy equation
    if (gas_eos == energy_eqn) 
      part.u = part.u0 + part.dudt0*dt;

    // Set particle as active at end of step
    if (dn == nstep) part.active = true;
    else part.active = false;

    // Flag all dead particles as inactive here
    if (part.itype == dead) part.active = false;
  }
  //---------------------------------------------------------------------------

  timing->EndTimingSection("SPH_LFKDK_ADVANCE_PARTICLES");

  return;
}
 


//=============================================================================
//  SphLeapfrogKDK::CorrectionTerms
/// Compute velocity integration to second order at the end of the step by 
/// adding a second order correction term.  The full integration becomes
/// $v(t+dt) = v(t) + 0.5*(a(t) + a(t+dt))*dt$. 
//=============================================================================
template <int ndim, template <int> class ParticleType>
void SphLeapfrogKDK<ndim, ParticleType>::CorrectionTerms
(int n,                             ///< [in] Integer time in block time struct
 FLOAT timestep,                    ///< [in] Base timestep value
 int Npart,                         ///< [in] Number of particles
 SphParticle<ndim>* sph_gen)        ///< [inout] Pointer to SPH particle array
{
  int dn;                           // Integer time since beginning of step
  int i;                            // Particle counter
  int k;                            // Dimension counter
  int nstep;                        // Particle (integer) step size
  ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sph_gen);

  debug2("[SphLeapfrogKDK::CorrectionTerms]");
  timing->StartTimingSection("SPH_LFKDK_CORRECTION_TERMS",2);


  //---------------------------------------------------------------------------
#pragma omp parallel for default(none) private(dn,i,k,nstep)\
  shared(n,Npart,sphdata,timestep)
  for (i=0; i<Npart; i++) {
    SphParticle<ndim>& part = sphdata[i];
    if (part.itype == dead) continue;

    // Compute time since beginning of current step
    dn = n - part.nlast;
    nstep = part.nstep;

    if (dn == nstep)
      for (k=0; k<ndim; k++) part.v[k] += timestep*(FLOAT) nstep*
        (FLOAT) 0.5*(part.a[k] - part.a0[k]);
    if (dn == nstep && gas_eos == energy_eqn) part.u +=
      0.5*(part.dudt - part.dudt0)*timestep*(FLOAT) nstep;

  }
  //---------------------------------------------------------------------------

  timing->EndTimingSection("SPH_LFKDK_CORRECTION_TERMS");

  return;
}



//=============================================================================
//  SphLeapfrogKDK::EndTimestep
/// Record all important SPH particle quantities at the end of the step for  
/// the start of the new timestep.
//=============================================================================
template <int ndim, template <int> class ParticleType>
void SphLeapfrogKDK<ndim, ParticleType>::EndTimestep
(int n,                             ///< [in] Integer time in block time struct
 FLOAT timestep,                    ///< [in] Base timestep value
 int Npart,                         ///< [in] Number of particles
 SphParticle<ndim>* sph_gen)        ///< [inout] Pointer to SPH particle array
{
  int dn;                           // Integer time since beginning of step
  int i;                            // Particle counter
  int k;                            // Dimension counter
  int nstep;                        // Particle (integer) step size
  ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sph_gen);

  debug2("[SphLeapfrogKDK::EndTimestep]");
  timing->StartTimingSection("SPH_LFKDK_END_TIMESTEP",2);


  //---------------------------------------------------------------------------
#pragma omp parallel for default(none) private(dn,i,k,nstep) \
  shared(n,Npart,sphdata,timestep)
  for (i=0; i<Npart; i++) {
    SphParticle<ndim>& part = sphdata[i];
    if (part.itype == dead) continue;

    dn = n - part.nlast;
    nstep = part.nstep;

    if (dn == nstep) {
      for (k=0; k<ndim; k++) part.v[k] += timestep*(FLOAT) nstep*
        (FLOAT) 0.5*(part.a[k] - part.a0[k]);
      for (k=0; k<ndim; k++) part.r0[k] = part.r[k];
      for (k=0; k<ndim; k++) part.v0[k] = part.v[k];
      for (k=0; k<ndim; k++) part.a0[k] = part.a[k];
      if (gas_eos == energy_eqn) {
        part.u += 0.5*(part.dudt - part.dudt0)*timestep*(FLOAT) nstep;
        part.u0 = part.u;
        part.dudt0 = part.dudt;
      }
      part.nlast = n;
      part.active = false;
    }
  }
  //---------------------------------------------------------------------------

  timing->EndTimingSection("SPH_LFKDK_END_TIMESTEP");

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
template <int ndim, template <int> class ParticleType>
int SphLeapfrogKDK<ndim, ParticleType>::CheckTimesteps
(int level_diff_max,                ///< [in] Max. allowed SPH neib dt diff
 int level_step,                    ///< [in] Level of base timestep
 int n,                             ///< [in] Integer time in block time struct
 int Npart,                         ///< [in] Number of particles
 SphParticle<ndim>* sph_gen)        ///< [inout] Pointer to SPH particle array
{
  int activecount = 0;              // No. of ptcls with new active timesteps
  int dn;                           // Integer time since beginning of step
  int i;                            // Particle counter
  int level_new;                    // New level of particle
  int nnewstep;                     // New step size of particle
  ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sph_gen);

  debug2("[SphLeapfrogKDK::CheckTimesteps]");
  timing->StartTimingSection("SPH_LFKDK_CHECK_TIMESTEPS",2);


  //---------------------------------------------------------------------------
#pragma omp parallel for default(none) private(dn,i,level_new,nnewstep)	\
  shared(level_diff_max,level_step,n,Npart,sphdata) reduction(+:activecount)
  for (i=0; i<Npart; i++) {
    SphParticle<ndim>& part = sphdata[i];
    if (part.itype == dead) continue;

    dn = n - part.nlast;
    if (dn == part.nstep) continue;

    // Check if neighbour timesteps are too small.  If so, then reduce 
    // timestep if possible
    if (part.levelneib - part.level > level_diff_max) {
      level_new = part.levelneib - level_diff_max;
      nnewstep = pow(2,level_step - level_new);

      // If new level is correctly synchronised, then change all quantities
      if (dn%nnewstep == 0) {
        if (dn > 0) part.nstep = dn;
        part.level = level_new;
        part.active = true;
        activecount++;
      }
    }
  }
  //---------------------------------------------------------------------------

  timing->EndTimingSection("SPH_LFKDK_CHECK_TIMESTEPS");

  return activecount;
}



// Template class instances for each dimensionality value (1, 2 and 3) and particle type
template class SphLeapfrogKDK<1, GradhSphParticle>;
template class SphLeapfrogKDK<2, GradhSphParticle>;
template class SphLeapfrogKDK<3, GradhSphParticle>;
template class SphLeapfrogKDK<1, SM2012SphParticle>;
template class SphLeapfrogKDK<2, SM2012SphParticle>;
template class SphLeapfrogKDK<3, SM2012SphParticle>;
