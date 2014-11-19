//=================================================================================================
//  EnergyPEC.cpp
//  Contains functions for energy equation integration using a
//  Predict-Evaluate-Correct (PEC) scheme.
//  N.B. this PEC scheme is the same as integrating the particle velocities
//  in the Leapfrog KDK scheme.
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
#include "Sph.h"
#include "SphKernel.h"
#include "SphIntegration.h"
#include "SphParticle.h"
#include "EOS.h"
#include "EnergyEquation.h"
#include "Debug.h"
using namespace std;



//=================================================================================================
//  EnergyEquation::EnergyEquation()
/// EnergyEquation constructor
//=================================================================================================
template <int ndim>
EnergyEquation<ndim>::EnergyEquation(DOUBLE energy_mult_aux) :
  energy_mult(energy_mult_aux)
{
}



//=================================================================================================
//  EnergyEquation::~EnergyEquation()
/// EnergyEquation destructor
//=================================================================================================
template <int ndim>
EnergyEquation<ndim>::~EnergyEquation()
{
}



// Class instances for each dimensionality (1, 2 and 3)
template class EnergyEquation<1>;
template class EnergyEquation<2>;
template class EnergyEquation<3>;
template class NullEnergy<1>;
template class NullEnergy<2>;
template class NullEnergy<3>;



//=================================================================================================
//  EnergyPEC::EnergyPEC()
/// EnergyPEC class constructor
//=================================================================================================
template <int ndim, template <int> class ParticleType>
EnergyPEC<ndim,ParticleType>::EnergyPEC(DOUBLE energy_mult_aux) :
  EnergyEquation<ndim>(energy_mult_aux)
{
}



//=================================================================================================
//  EnergyPEC::~EnergyPEC()
/// EnergyPEC class destructor
//=================================================================================================
template <int ndim, template <int> class ParticleType>
EnergyPEC<ndim,ParticleType>::~EnergyPEC()
{
}



//=================================================================================================
//  EnergyPEC::EnergyIntegration
/// Integrate internal energy to first order from the beginning of the step to
/// the current simulation time, i.e. u(t+dt) = u(t) + dudt(t)*dt
//=================================================================================================
template <int ndim, template <int> class ParticleType>
void EnergyPEC<ndim,ParticleType>::EnergyIntegration
 (const int n,                         ///< [in] Integer time in block time struct
  const int Npart,                     ///< [in] Number of particles
  const FLOAT t,                       ///< [in] Current simulation time
  const FLOAT timestep,                ///< [in] Base timestep value
  SphParticle<ndim>* sph_gen)          ///< [inout] Pointer to SPH particle array
{
  int i;                               // Particle counter
  FLOAT dt;                            // Timestep since start of step
  ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sph_gen);

  debug2("[EnergyPEC::EnergyIntegration]");
  timing->StartTimingSection("ENERGY_PEC_INTEGRATION",2);

  //-----------------------------------------------------------------------------------------------
#pragma omp parallel for default(none) private(dt,i) shared(Npart,sphdata)
  for (i=0; i<Npart; i++) {
    SphParticle<ndim>& part = sphdata[i];
    if (part.itype == dead) continue;

    // Compute time since beginning of current step
    dt = t - part.tlast;
    part.u = part.u0 + part.dudt0*dt;
  }
  //-----------------------------------------------------------------------------------------------

  timing->EndTimingSection("ENERGY_PEC_INTEGRATION");

  return;
}



//=================================================================================================
//  EnergyPEC::CorrectionTerms
/// Compute energy integration to second order at the end of the step by
/// adding a second order correction term.  The full integration becomes
/// $u(t+dt) = u(t) + 0.5*(dudt(t) + dudt(t+dt))*dt$.
//=================================================================================================
template <int ndim, template <int> class ParticleType>
void EnergyPEC<ndim,ParticleType>::EnergyCorrectionTerms
 (const int n,                         ///< [in] Integer time in block time struct
  const int Npart,                     ///< [in] Number of particles
  const FLOAT t,                       ///< [in] Current simulation time
  const FLOAT timestep,                ///< [in] Base timestep value
  SphParticle<ndim>* sph_gen)          ///< [inout] Pointer to SPH particle array
{
  int dn;                           // Integer time since beginning of step
  int i;                            // Particle counter
  ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sph_gen);

  debug2("[EnergyPEC::EnergyCorrectionTerms]");
  timing->StartTimingSection("ENERGY_PEC_CORRECTION_TERMS",2);

  //-----------------------------------------------------------------------------------------------
#pragma omp parallel for default(none) private(dn,i) shared(n,Npart,sphdata,timestep)
  for (i=0; i<Npart; i++) {
    SphParticle<ndim>& part = sphdata[i];
    if (part.itype == dead) continue;

    // Compute time since beginning of current step
    dn = n - part.nlast;

    if (dn == part.nstep) {
      part.u += 0.5*(part.dudt - part.dudt0)*(t - part.tlast); //timestep*(FLOAT) nstep;
    }

  }
  //-----------------------------------------------------------------------------------------------

  timing->EndTimingSection("ENERGY_PEC_CORRECTION_TERMS");

  return;
}



//=================================================================================================
//  EnergyPEC::EndTimestep
/// Record all important thermal quantities at the end of the step for start of the new timestep.
//=================================================================================================
template <int ndim, template <int> class ParticleType>
void EnergyPEC<ndim,ParticleType>::EndTimestep
 (const int n,                         ///< [in] Integer time in block time struct
  const int Npart,                     ///< [in] Number of particles
  const FLOAT t,                       ///< [in] Current simulation time
  const FLOAT timestep,                ///< [in] Base timestep value
  SphParticle<ndim>* sph_gen)          ///< [inout] Pointer to SPH particle array
{
  int dn;                           // Integer time since beginning of step
  int i;                            // Particle counter
  ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sph_gen);

  debug2("[EnergyPEC::EndTimestep]");
  timing->StartTimingSection("ENERGY_PEC_END_TIMESTEP",2);

  //-----------------------------------------------------------------------------------------------
#pragma omp parallel for default(none) private(dn,i) shared(n,Npart,sphdata,timestep)
  for (i=0; i<Npart; i++) {
    SphParticle<ndim>& part = sphdata[i];
    if (part.itype == dead) continue;
    dn = n - part.nlast;

    if (dn == part.nstep) {
      part.u     += 0.5*(part.dudt - part.dudt0)*(t - part.tlast); //timestep*(FLOAT) nstep;
      part.u0    = part.u;
      part.dudt0 = part.dudt;
    }

  }
  //-----------------------------------------------------------------------------------------------

  timing->EndTimingSection("ENERGY_PEC_END_TIMESTEP");

  return;
}



//=================================================================================================
//  EnergyPEC::Timestep
/// Compute explicit timestep such that u cannot change by a large fraction in one step,
/// i.e. dt = const*u/|dudt + epsilon| where epsilon is to prevent the denominator becoming zero.
//=================================================================================================
template <int ndim, template <int> class ParticleType>
DOUBLE EnergyPEC<ndim,ParticleType>::Timestep
(SphParticle<ndim> &part)           ///< [inout] SPH particle reference
{
  return this->energy_mult*(DOUBLE) (part.u/(fabs(part.dudt) + small_number));
}



template class EnergyPEC<1, GradhSphParticle>;
template class EnergyPEC<2, GradhSphParticle>;
template class EnergyPEC<3, GradhSphParticle>;
template class EnergyPEC<1, SM2012SphParticle>;
template class EnergyPEC<2, SM2012SphParticle>;
template class EnergyPEC<3, SM2012SphParticle>;
