//=================================================================================================
//  SphLeapfrogDKD.cpp
//  Contains functions for integrating SPH particle positions and velocities
//  using the leapfrog drift-kick-drift (DKD) scheme.
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
//=================================================================================================


#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <math.h>

#include "../Headers/Integration.h"
#include "Sph.h"
#include "SmoothingKernel.h"
#include "Particle.h"
#include "EOS.h"
#include "Debug.h"
using namespace std;



//=================================================================================================
//  SphLeapfrogDKD::SphLeapfrogDKD
/// SphLeapfrogDKD class constructor
//=================================================================================================
template <int ndim, template <int> class ParticleType>
SphLeapfrogDKD<ndim, ParticleType>::SphLeapfrogDKD(DOUBLE accel_mult_aux, DOUBLE courant_mult_aux,
                                                   DOUBLE energy_mult_aux, bool energy_integration,
                                                   tdaviscenum tdavisc_aux) :
  SphIntegration<ndim>(accel_mult_aux, courant_mult_aux, energy_mult_aux, energy_integration,
                       tdavisc_aux)
{
}



//=================================================================================================
//  SphLeapfrogDKD::~SphLeapfrog()
/// SphLeapfrogDKD class destructor
//=================================================================================================
template <int ndim, template <int> class ParticleType>
SphLeapfrogDKD<ndim, ParticleType>::~SphLeapfrogDKD()
{
}



//=================================================================================================
//  SphLeapfrogDKD::AdvanceParticles
/// Integrate both particle positions and velocities to 1st order from the
/// beginning of the step to the current simulation time, i.e.
/// $r(t+dt) = r(t) + v(t)*dt$,
/// $v(t+dt) = v(t) + a(t)*dt$.
/// If particle has reached the half-step, then set as active to compute
/// the 'kick' acceleration terms.
//=================================================================================================
template <int ndim, template <int> class ParticleType>
void SphLeapfrogDKD<ndim, ParticleType>::AdvanceParticles
 (const int n,                         ///< [in] Integer time in block time struct
  const FLOAT t,                       ///< [in] Current simulation time
  const FLOAT timestep,                ///< [in] Base timestep value
  Hydrodynamics<ndim>* hydro)
{
  int dn;                              // Integer time since beginning of step
  int nstep;                           // Particle (integer) step size
  int i;                               // Particle counter
  int k;                               // Dimension counter
  FLOAT dt;                            // Timestep since start of step

  debug2("[SphLeapfrogDKD::AdvanceParticles]");
  CodeTiming::BlockTimer timer = timing->StartNewTimer("SPH_LFDKD_ADVANCE_PARTICLES");

  Sph<ndim>* sph = reinterpret_cast<Sph<ndim>*>(hydro);
  ParticleType<ndim>* sphdata = reinterpret_cast<ParticleType<ndim>*>(sph->GetSphParticleArray());

  // Advance positions and velocities of all SPH particles
  //-----------------------------------------------------------------------------------------------
#pragma omp parallel for default(none) private(dn,dt,i,k,nstep) shared(n,sphdata,sph,t,timestep)
  for (i=0; i<sph->Nhydro; i++) {
    SphParticle<ndim>& part = sphdata[i];
    if (part.flags.is_dead()) continue;

    // Compute time since beginning of current step
    nstep = part.nstep;
    dn = n - part.nlast;
    //dt = timestep*(FLOAT) dn;
    dt = t - part.tlast;

    // Advance particle positions and velocities depending on if we're before
    // or after the half-step.
    if (dn < nstep) {
      for (k=0; k<ndim; k++) part.r[k] = part.r0[k] + part.v0[k]*dt;
      for (k=0; k<ndim; k++) part.v[k] = part.v0[k] + part.a[k]*dt;
    }
    else {
      for (k=0; k<ndim; k++) part.v[k] = part.v0[k] + part.a[k]*dt;
      for (k=0; k<ndim; k++) part.r[k] = part.r0[k] + (FLOAT) 0.5*(part.v0[k] + part.v[k])*dt;
    }

    // Integrate time-dependent viscosity
    if (tdavisc != notdav) part.alpha += part.dalphadt*timestep;

    // Integrate explicit energy equation
    if (energy_integration) part.u = part.u0 + part.dudt*dt;

    // Set particle as active at end of step
    if (dn == nstep/2) part.flags.set(active);
    else part.flags.unset(active);

    // Flag all dead particles as inactive here
    if (part.flags.is_dead()) part.flags.unset(active);

  }
  //-----------------------------------------------------------------------------------------------

  return;
}
//=================================================================================================
//  SphLeapfrogDKD::SetActiveParticles
/// Set or unset the active flag for all particles based upon whther the particles need a force
/// calculation this timestep.
//=================================================================================================
template <int ndim, template <int> class ParticleType>
void SphLeapfrogDKD<ndim, ParticleType>::SetActiveParticles
(const int n,                         ///< [in] Current timestep number
 Hydrodynamics<ndim>* hydro)
{
  Sph<ndim>* sph = reinterpret_cast<Sph<ndim>*>(hydro);
  ParticleType<ndim>* sphdata = reinterpret_cast<ParticleType<ndim>*>(sph->GetSphParticleArray());

#pragma omp parallel for default(none) shared(n,sphdata,sph)
  for (int i=0; i<sph->Nhydro; i++) {
    SphParticle<ndim>& part = sphdata[i];
    int dn = n - part.nlast;

    // Force calculation is at mid-point of step
    if (dn == (part.nstep/2))
      part.flags.set(active);
    else
      part.flags.unset(active);
  }
}



//=================================================================================================
//  SphLeapfrogDKD::EndTimestep
/// Record all important SPH particle quantities at the end of the step for
/// the start of the new timestep.
//=================================================================================================
template <int ndim, template <int> class ParticleType>
void SphLeapfrogDKD<ndim, ParticleType>::EndTimestep
 (const int n,                         ///< [in] Integer time in block time struct
  const FLOAT t,                       ///< [in] Current simulation time
  const FLOAT timestep,                ///< [in] Base timestep value
  Hydrodynamics<ndim>* hydro)
{
  int dn;                              // Integer time since beginning of step
  int nstep;                           // Particle (integer) step size
  int i;                               // Particle counter
  int k;                               // Dimension counter

  debug2("[SphLeapfrogDKD::EndTimestep]");
  CodeTiming::BlockTimer timer = timing->StartNewTimer("SPH_LFDKD_END_TIMESTEP");

  Sph<ndim>* sph = reinterpret_cast<Sph<ndim>*>(hydro);
  ParticleType<ndim>* sphdata = reinterpret_cast<ParticleType<ndim>*>(sph->GetSphParticleArray());

  //-----------------------------------------------------------------------------------------------
#pragma omp parallel for default(none) private(dn,i,k,nstep) shared(n,sphdata,sph,t,timestep)
  for (i=0; i<sph->Nhydro; i++) {
    SphParticle<ndim>& part = sphdata[i];
    if (part.flags.is_dead()) continue;

    if (part.flags.check(end_timestep)) {
      for (k=0; k<ndim; k++) part.r0[k] = part.r[k];
      for (k=0; k<ndim; k++) part.v0[k] = part.v[k];
      for (k=0; k<ndim; k++) part.a0[k] = part.a[k];

      if (energy_integration) {
        part.u0 = part.u;
      }

      part.nlast   = n;
      part.tlast   = t;
      part.dt      = part.dt_next;
      part.dt_next = 0;
      part.flags.unset(active);
      part.flags.unset(end_timestep);
    }
  }
  //-----------------------------------------------------------------------------------------------

  return;
}



//=================================================================================================
//  SphLeapfrogDKD::CheckTimesteps
/// Record all important SPH particle quantities at the end of the step for
/// the start of the new timestep.
//=================================================================================================
template <int ndim, template <int> class ParticleType>
int SphLeapfrogDKD<ndim, ParticleType>::CheckTimesteps
 (const int level_diff_max,            ///< [in] Max. allowed SPH neib dt diff
  const int level_step,                ///< [in] Level of base timestep
  const int n,                         ///< [in] Integer time in block time struct
  const FLOAT timestep,                ///< [in] Current time-step
  Hydrodynamics<ndim>* hydro)
{
  int dn;                              // Integer time since beginning of step
  int level_new;                       // New timestep level
  int nnewstep;                        // New integer timestep
  int activecount = 0;                 // No. of newly active particles
  int i;                               // Particle counter

  debug2("[SphLeapfrogDKD::CheckTimesteps]");
  CodeTiming::BlockTimer timer = timing->StartNewTimer("SPH_LFDKD_CHECK_TIMESTEPS");

  Sph<ndim>* sph = reinterpret_cast<Sph<ndim>*>(hydro);
  ParticleType<ndim>* sphdata = reinterpret_cast<ParticleType<ndim>*>(sph->GetSphParticleArray());

  //-----------------------------------------------------------------------------------------------
#pragma omp parallel for default(none) private(dn,level_new,nnewstep) \
  shared(level_diff_max,n,level_step,sphdata,sph) reduction(+:activecount)
  for (i=0; i<sph->Nhydro; i++) {
    SphParticle<ndim>& part = sphdata[i];
    if (part.flags.is_dead()) continue;

    dn = n - part.nlast;
    if (dn == part.nstep) continue;

    // Check if neighbour timesteps are too small.  If so, then reduce timestep if possible
    if (part.levelneib - part.level > level_diff_max) {
      level_new = part.levelneib - level_diff_max;
      nnewstep = pow(2,level_step - level_new);

      // If new level is correctly synchronised at the half-step of the
      // new-step (where acceleration is computed), then change all quantities
      if ((2*dn)%nnewstep == 0) {
        if (dn > 0) part.nstep = dn;
        part.level = level_new;
        part.flags.set(active);
        activecount++;
      }
    }
  }
  //-----------------------------------------------------------------------------------------------

  return activecount;
}



// Template class instances for each dimensionality value (1, 2 and 3)
template class SphLeapfrogDKD<1, GradhSphParticle>;
template class SphLeapfrogDKD<2, GradhSphParticle>;
template class SphLeapfrogDKD<3, GradhSphParticle>;
template class SphLeapfrogDKD<1, SM2012SphParticle>;
template class SphLeapfrogDKD<2, SM2012SphParticle>;
template class SphLeapfrogDKD<3, SM2012SphParticle>;
