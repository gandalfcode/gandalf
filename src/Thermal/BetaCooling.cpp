//=================================================================================================
//  BetaCooling
//  Simple beta-cooling model for self-gravitating discs
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

#include <math.h>
#include <iostream>
#include <vector>

#include "EnergyEquation.h"
#include "Parameters.h"
#include "SimUnits.h"
#include "Nbody.h"


template<int ndim, template <int> class ParticleType>
BetaCooling<ndim,ParticleType>::BetaCooling
(Parameters* params,
 SimUnits* units)
 : EnergyEquation<ndim>(params->floatparams["energy_mult"]),
   Beta(params->floatparams["BetaCoolingParam"]),
   nbody(NULL)
{  };



//=================================================================================================
//  BetaCooling::EnergyIntegration
/// Integrate internal energy to first order from the beginning of the step to
/// the current simulation time, i.e. u(t+dt) = u(t) + dudt(t)*dt .
//=================================================================================================
template <int ndim, template <int> class ParticleType>
void BetaCooling<ndim,ParticleType>::EnergyIntegration
 (const int n,                         ///< [in] Integer time in block time struct
  const FLOAT t,                       ///< [in] Current simulation time
  const FLOAT timestep,                ///< [in] Base timestep value
  Hydrodynamics<ndim>* hydro)
{
  int i;                               // Particle counter
  ParticleType<ndim>* partdata = hydro->template GetParticleArray<ParticleType>();
  StarParticle<ndim>* star = nbody->stardata;


  debug2("[BetaCooling::EnergyIntegration]");
  CodeTiming::BlockTimer timer = timing->StartNewTimer("BETACOOLING_INTEGRATION");


  //-----------------------------------------------------------------------------------------------
#pragma omp parallel for default(none) private(i) shared(partdata, hydro, star)
  for (i=0; i<hydro->Nhydro; i++) {
    ParticleType<ndim>& part = partdata[i];
    if (part.flags.is_dead()) continue;

    FLOAT dr[ndim] ;
    for (int j=0; j<ndim;j++) dr[j] = part.r[j]-star[0].r[j];
    FLOAT r2 = DotProduct(dr,dr,ndim);
    FLOAT r = sqrt(r2);

    FLOAT Omega = sqrt(star[0].m / (r*r2)) ;

    // Compute time since beginning of current step
    FLOAT dt = timestep*(n - part.nlast);

    // Use an approximate form that gives the correct limits
    part.u = (part.u0 + part.dudt0*dt) / (1 + dt * Omega / Beta) ;
  }
  //-----------------------------------------------------------------------------------------------

  return;
}


//=================================================================================================
//  BetaCooling::EndTimestep
/// Compute the temperature at the end of the step
//=================================================================================================
template <int ndim, template <int> class ParticleType>
void BetaCooling<ndim,ParticleType>::EndTimestep
 (const int n,                         ///< [in] Integer time in block time struct
  const FLOAT t,                       ///< [in] Current simulation time
  const FLOAT timestep,                ///< [in] Base timestep value
  Hydrodynamics<ndim>* hydro)
{
  int i;                               // Particle counter
  ParticleType<ndim>* partdata = hydro->template GetParticleArray<ParticleType>();
  StarParticle<ndim>* star = nbody->stardata;


  debug2("[BetaCooling::EnergyIntegration]");
  CodeTiming::BlockTimer timer = timing->StartNewTimer("BETACOOLING_INTEGRATION");


  //-----------------------------------------------------------------------------------------------
#pragma omp parallel for default(none) private(i) shared(partdata, hydro, star)
  for (i=0; i<hydro->Nhydro; i++) {
    ParticleType<ndim>& part = partdata[i];
    if (part.flags.is_dead()) continue;

    if (part.flags.check(end_timestep)) {
      FLOAT dr[ndim] ;
      for (int j=0; j<ndim;j++) dr[j] = part.r[j]-star[0].r[j];
      FLOAT r2 = DotProduct(dr,dr,ndim);
      FLOAT r = sqrt(r2);

      FLOAT Omega = sqrt(star[0].m / (r*r2)) ;

      // Compute time since beginning of current step
      FLOAT dt = timestep*(n - part.nlast);

      // Use an approximate form that gives the correct limits
      part.u = (part.u0 + 0.5*(part.dudt + part.dudt0)*dt) ;
      if (part.u < 0) part.u = part.u0 + part.dudt0*dt ;
      part.u /= (1 + dt  * Omega / Beta) ;
    }
  }
  //-----------------------------------------------------------------------------------------------

  return;
}

template class BetaCooling<1, GradhSphParticle>;
template class BetaCooling<2, GradhSphParticle>;
template class BetaCooling<3, GradhSphParticle>;
template class BetaCooling<1, SM2012SphParticle>;
template class BetaCooling<2, SM2012SphParticle>;
template class BetaCooling<3, SM2012SphParticle>;

