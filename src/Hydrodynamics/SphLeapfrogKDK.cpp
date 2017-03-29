//=================================================================================================
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
//=================================================================================================


#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <math.h>

#include "../Headers/Integration.h"
#include "Debug.h"
#include "EOS.h"
#include "Particle.h"
#include "SmoothingKernel.h"
#include "Sph.h"
using namespace std;



//=================================================================================================
//  SphLeapfrogKDK::SphLeapfrogKDK
/// SphLeapfrogKDK class constructor
//=================================================================================================
template <int ndim, template <int> class ParticleType>
SphLeapfrogKDK<ndim, ParticleType>::SphLeapfrogKDK
 (DOUBLE _accel_mult, DOUBLE _courant_mult, DOUBLE _energy_mult,
  eosenum _gas_eos, tdaviscenum _tdavisc) :
  SphIntegration<ndim>(_accel_mult, _courant_mult, _energy_mult, _gas_eos, _tdavisc)
{
}



//=================================================================================================
//  SphLeapfrogKDK::~SphLeapfrog()
/// SphLeapfrogKDK class destructor
//=================================================================================================
template <int ndim, template <int> class ParticleType>
SphLeapfrogKDK<ndim, ParticleType>::~SphLeapfrogKDK()
{
}



//=================================================================================================
//  SphLeapfrogKDK::AdvanceParticles
/// Integrate particle positions to 2nd order, and particle velocities to 1st
/// order from the beginning of the step to the current simulation time, i.e.
/// $r(t+dt) = r(t) + v(t)*dt + 0.5*a(t)*dt^2$,
/// $v(t+dt) = v(t) + a(t)*dt$.
/// Also set particles at the end of step as 'active' in order to compute
/// the end-of-step force computation.
//=================================================================================================
template <int ndim, template <int> class ParticleType>
void SphLeapfrogKDK<ndim, ParticleType >::AdvanceParticles
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

  debug2("[SphLeapfrogKDK::AdvanceParticles]");
  CodeTiming::BlockTimer timer = timing->StartNewTimer("SPH_LFKDK_ADVANCE_PARTICLES");

  Sph<ndim>* sph = reinterpret_cast<Sph<ndim>*>(hydro);
  ParticleType<ndim>* sphdata = reinterpret_cast<ParticleType<ndim>*>(sph->GetSphParticleArray());

  // Advance positions and velocities of all SPH particles
  //-----------------------------------------------------------------------------------------------
#pragma omp parallel for default(none) private(dn,dt,i,k,nstep) shared(sphdata, sph)
  for (i=0; i<sph->Nhydro; i++) {
    SphParticle<ndim>& part = sph->GetSphParticlePointer(i);
    if (part.flags.is_dead()) continue;

    // Compute time since beginning of current step
    nstep = part.nstep;
    dn = n - part.nlast;
    //dt = timestep*(FLOAT) dn;
    dt = t - part.tlast;

    // Advance particle positions and velocities
    for (k=0; k<ndim; k++) part.r[k] = part.r0[k] + part.v0[k]*dt + 0.5*part.a0[k]*dt*dt;
    for (k=0; k<ndim; k++) part.v[k] = part.v0[k] + part.a0[k]*dt;

    // Integrate time-dependent viscosity
    if (tdavisc != notdav) part.alpha += part.dalphadt*timestep;

    // Integrate explicit energy equation
    if (gas_eos == energy_eqn) part.u = part.u0 + part.dudt0*dt;

    // Set particle as active at end of step
    if (dn == nstep) part.flags.set(active);
    else part.flags.unset(active);

    // Flag all dead particles as inactive here
    if (part.flags.is_dead()) part.flags.unset(active);

  }
  //-----------------------------------------------------------------------------------------------

  return;
}



//=================================================================================================
//  SphLeapfrogKDK::CorrectionTerms
/// Compute velocity integration to second order at the end of the step by
/// adding a second order correction term.  The full integration becomes
/// $v(t+dt) = v(t) + 0.5*(a(t) + a(t+dt))*dt$.
//=================================================================================================
template <int ndim, template <int> class ParticleType>
void SphLeapfrogKDK<ndim, ParticleType>::CorrectionTerms
 (const int n,                         ///< [in] Integer time in block time struct
  const FLOAT t,                       ///< [in] Current simulation time
  const FLOAT timestep,                ///< [in] Base timestep value
  Hydrodynamics<ndim>* hydro)
{
  int dn;                              // Integer time since beginning of step
  int nstep;                           // Particle (integer) step size
  int i;                               // Particle counter
  int k;                               // Dimension counter

  debug2("[SphLeapfrogKDK::CorrectionTerms]");
  CodeTiming::BlockTimer timer = timing->StartNewTimer("SPH_LFKDK_CORRECTION_TERMS");

  Sph<ndim>* sph = reinterpret_cast<Sph<ndim>*>(hydro);
  ParticleType<ndim>* sphdata = reinterpret_cast<ParticleType<ndim>*>(sph->GetSphParticleArray());

  //-----------------------------------------------------------------------------------------------
#pragma omp parallel for default(none) private(dn,i,k,nstep)  shared(sphdata, sph)
  for (i=0; i<sph->Nhydro; i++) {
    SphParticle<ndim>& part = sphdata[i];
    if (part.flags.is_dead()) continue;

    // Compute time since beginning of current step
    dn = n - part.nlast;
    nstep = part.nstep;

    if (dn == nstep) {
      for (k=0; k<ndim; k++) part.v[k] += (t - part.tlast)* //timestep*(FLOAT) nstep*
        (FLOAT) 0.5*(part.a[k] - part.a0[k]);
    }
    if (dn == nstep && gas_eos == energy_eqn) {
      part.u += 0.5*(part.dudt - part.dudt0)*(t - part.tlast); //timestep*(FLOAT) nstep;

      // In spurious cases where correction term can lead to negative energies, simply use
      // 1st-order integration (which should guarantee positive energy due to the CFL condition)
      if (part.u <= (FLOAT) 0.0) part.u = part.u0 + part.dudt0*(t - part.tlast);
    }

  }
  //-----------------------------------------------------------------------------------------------

  return;
}


//=================================================================================================
//  SphLeapfrogKDK::SetActiveParticles
/// Set or unset the active flag for all particles based upon whther the particles need a force
/// calculation this timestep.
//=================================================================================================
template <int ndim, template <int> class ParticleType>
void SphLeapfrogKDK<ndim, ParticleType>::SetActiveParticles
(const int n,                         ///< [in] Current timestep number
 Hydrodynamics<ndim>* hydro)
{
  Sph<ndim>* sph = reinterpret_cast<Sph<ndim>*>(hydro);
  ParticleType<ndim>* sphdata = reinterpret_cast<ParticleType<ndim>*>(sph->GetSphParticleArray());

#pragma omp parallel for default(none) shared(sphdata, sph)
  for (int i=0; i<sph->Nhydro; i++) {
    SphParticle<ndim>& part = sphdata[i];

    int dn = n - part.nlast;

    // Force calculation is at end of step
    if (part.flags.check(end_timestep) || dn == part.nstep)
      part.flags.set(active);
    else
      part.flags.unset(active);
  }
}

//=================================================================================================
//  SphLeapfrogKDK::EndTimestep
/// Record all important SPH particle quantities at the end of the step for
/// the start of the new timestep.
//=================================================================================================
template <int ndim, template <int> class ParticleType>
void SphLeapfrogKDK<ndim, ParticleType>::EndTimestep
 (const int n,                         ///< [in] Integer time in block time struct
  const FLOAT t,                       ///< [in] Current simulation time
  const FLOAT timestep,                ///< [in] Base timestep value
  Hydrodynamics<ndim>* hydro)
{
  int dn;                              // Integer time since beginning of step
  int nstep;                           // Particle (integer) step size
  int i;                               // Particle counter
  int k;                               // Dimension counter

  debug2("[SphLeapfrogKDK::EndTimestep]");
  CodeTiming::BlockTimer timer = timing->StartNewTimer("SPH_LFKDK_END_TIMESTEP");

  Sph<ndim>* sph = reinterpret_cast<Sph<ndim>*>(hydro);
  ParticleType<ndim>* sphdata = reinterpret_cast<ParticleType<ndim>*>(sph->GetSphParticleArray());

  //-----------------------------------------------------------------------------------------------
#pragma omp parallel for default(none) private(dn,i,k,nstep) shared(sphdata, sph)
  for (i=0; i<sph->Nhydro; i++) {
    SphParticle<ndim>& part = sphdata[i];
    if (part.flags.is_dead()) continue;


    if (part.flags.check(end_timestep)) {
      for (k=0; k<ndim; k++) part.v[k] += 0.5 * (t - part.tlast) * (part.a[k] - part.a0[k]);
      for (k=0; k<ndim; k++) part.r0[k] = part.r[k];
      for (k=0; k<ndim; k++) part.v0[k] = part.v[k];
      for (k=0; k<ndim; k++) part.a0[k] = part.a[k];

      // If using an adiabatic energy equation, then explictly integrate the internal energy
      if (gas_eos == energy_eqn) {
        part.u     += 0.5*(part.dudt - part.dudt0)*(t - part.tlast); //timestep*(FLOAT) nstep;

        // In spurious cases where correction term can lead to negative energies, simply use
        // 1st-order integration (which should guarantee positive energy due to the CFL condition)
        if (part.u <= (FLOAT) 0.0) part.u = part.u0 + part.dudt0*(t - part.tlast);

        part.u0    = part.u;
        part.dudt0 = part.dudt;
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
//  SphLeapfrogKDK::CheckTimesteps
/// Check through all SPH particles to see if the maximum neighbour timestep exceeds the maximum
/// allowed difference (level_diff_max).  If so, then check if we can prematurely finish the
/// current timestep and move the to a lower timestep level.  Returns the number of particles
/// whose timesteps have been reduced.
//=================================================================================================
template <int ndim, template <int> class ParticleType>
int SphLeapfrogKDK<ndim, ParticleType>::CheckTimesteps
 (const int level_diff_max,            ///< [in] Max. allowed SPH neib dt diff
  const int level_step,                ///< [in] Level of base timestep
  const int n,                         ///< [in] Integer time in block time struct
  const FLOAT timetep,                 ///< [in] Current time-step
  Hydrodynamics<ndim>* hydro)
{
  int dn;                              // Integer time since beginning of step
  int level_new;                       // New timestep level
  int nnewstep;                        // New integer timestep
  int activecount = 0;                 // No. of newly active particles
  int i;                               // Particle counter

  debug2("[SphLeapfrogKDK::CheckTimesteps]");
  CodeTiming::BlockTimer timer = timing->StartNewTimer("SPH_LFKDK_CHECK_TIMESTEPS");

  Sph<ndim>* sph = reinterpret_cast<Sph<ndim>*>(hydro);
  ParticleType<ndim>* sphdata = reinterpret_cast<ParticleType<ndim>*>(sph->GetSphParticleArray());
  //-----------------------------------------------------------------------------------------------
#pragma omp parallel for default(none) private(dn,i,level_new,nnewstep) \
  shared(sphdata, sph) reduction(+:activecount)
  for (i=0; i<sph->Nhydro; i++) {
    SphParticle<ndim>& part = sphdata[i];
    if (part.flags.is_dead()) continue;

    dn = n - part.nlast;
    if (dn == part.nstep) continue;

    // Check if neighbour timesteps are too small.  If so, then reduce timestep if possible
    if (part.levelneib - part.level > level_diff_max) {
      level_new = part.levelneib - level_diff_max;
      nnewstep  = pow(2,level_step - level_new);

      // If new level is correctly synchronised, then change all quantities
      if (dn%nnewstep == 0) {
        if (dn > 0) part.nstep = dn;
        part.level  = level_new;
        part.flags.set(active);
        activecount++;
      }
    }
  }
  //-----------------------------------------------------------------------------------------------

  return activecount;
}



// Template class instances for each dimensionality value (1, 2 and 3) and particle type
template class SphLeapfrogKDK<1, GradhSphParticle>;
template class SphLeapfrogKDK<2, GradhSphParticle>;
template class SphLeapfrogKDK<3, GradhSphParticle>;
template class SphLeapfrogKDK<1, SM2012SphParticle>;
template class SphLeapfrogKDK<2, SM2012SphParticle>;
template class SphLeapfrogKDK<3, SM2012SphParticle>;
