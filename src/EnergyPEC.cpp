//=============================================================================
//  EnergyPEC.cpp
//  Contains functions for energy equation integration using a 
//  Predict-Evalulate-Correct (PEC) scheme.
//  N.B. this PEC scheme is the same as integrating the particle velocities 
//  in the Leapfrog KDK scheme.
//
//  This file is part of GANDALF :
//  Graphical Astrophysics code for N-body Dynamics and Lagrangian Fluids
//  https://github.com/gandalfcode/gandalf
//  Contact : gandalfcode@gmail.com
//
//  Copyright (C) 2013  D. A. Hubber, G Rosotti
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
#include "EnergyEquation.h"
#include "Debug.h"
using namespace std;




//=============================================================================
//  EnergyEquation::EnergyEquation()
/// EnergyEquation constructor
//=============================================================================
template <int ndim>
EnergyEquation<ndim>::EnergyEquation(DOUBLE energy_mult_aux) :
  energy_mult(energy_mult_aux)
{
}



//=============================================================================
//  EnergyEquation::~EnergyEquation()
/// EnergyEquation destructor
//=============================================================================
template <int ndim>
EnergyEquation<ndim>::~EnergyEquation()
{
}



// Class instances for each dimensionality (1, 2 and 3)
template class EnergyEquation<1>;
template class EnergyEquation<2>;
template class EnergyEquation<3>;



//=============================================================================
//  EnergyPEC::EnergyPEC()
/// EnergyPEC class constructor
//=============================================================================
template <int ndim>
EnergyPEC<ndim>::EnergyPEC(DOUBLE energy_mult_aux) :
  EnergyEquation<ndim>(energy_mult_aux)
{
}



//=============================================================================
//  EnergyPEC::~EnergyPEC()
/// EnergyPEC class destructor
//=============================================================================
template <int ndim>
EnergyPEC<ndim>::~EnergyPEC()
{
}



//=============================================================================
//  EnergyPEC::EnergyIntegration
/// Integrate internal energy to first order from the beginning of the step to 
/// the current simulation time, i.e. u(t+dt) = u(t) + dudt(t)*dt
//=============================================================================
template <int ndim>
void EnergyPEC<ndim>::EnergyIntegration
(int n,                             ///< [in] Integer time in block time struct
 int Nsph,                          ///< [in] No. of SPH particles
 SphParticle<ndim> *sphdata,        ///< [inout] SPH particle data array
 FLOAT timestep)                    ///< [in] Base timestep value
{
  int dn;                           // Integer time since beginning of step
  int i;                            // Particle counter
  int nstep;                        // Particle (integer) step size
  FLOAT dt;                         // Timestep since start of step

  debug2("[EnergyPEC::EnergyIntegration]");

  // --------------------------------------------------------------------------
#pragma omp parallel for default(none) private(dt,nstep,i,dn) \
     shared(sphdata, n, timestep, Nsph)
  for (i=0; i<Nsph; i++) {
    nstep = sphdata[i].nstep;
    dn = n - sphdata[i].nlast;
    dt = timestep*(FLOAT) dn;
    sphdata[i].u = sphdata[i].u0 + sphdata[i].dudt0*dt;
  }
  // --------------------------------------------------------------------------

  return;
}
 


//=============================================================================
//  EnergyPEC::CorrectionTerms
/// Compute energy integration to second order at the end of the step by 
/// adding a second order correction term.  The full integration becomes
/// $u(t+dt) = u(t) + 0.5*(dudt(t) + dudt(t+dt))*dt$.
//=============================================================================
template <int ndim>
void EnergyPEC<ndim>::EnergyCorrectionTerms
(int n,                             ///< [in] Integer time in block time struct
 int Nsph,                          ///< [in] No. of SPH particles
 SphParticle<ndim> *sphdata,        ///< [inout] SPH particle data array
 FLOAT timestep)                    ///< [in] Base timestep value
{
  int dn;                           // Integer time since beginning of step
  int i;                            // Particle counter
  int nstep;                        // Particle (integer) step size

  debug2("[EnergyPEC::EnergyCorrectionTerms]");

  // --------------------------------------------------------------------------
#pragma omp parallel for default(none) private(nstep,dn,i) \
     shared(n,Nsph,sphdata,timestep)
  for (i=0; i<Nsph; i++) {
    dn = n - sphdata[i].nlast;
    nstep = sphdata[i].nstep;
    if (dn == nstep) sphdata[i].u += 
      0.5*(sphdata[i].dudt - sphdata[i].dudt0)*timestep*(FLOAT) nstep;
  }
  // --------------------------------------------------------------------------

  return;
}



//=============================================================================
//  EnergyPEC::EndTimestep
/// Record all important thermal quantities at the end of the step for the 
/// start of the new timestep.
//=============================================================================
template <int ndim>
void EnergyPEC<ndim>::EndTimestep
(int n,                             ///< [in] Integer time in block time struct
 int Nsph,                          ///< [in] No. of SPH particles
 SphParticle<ndim> *sphdata)        ///< [inout] SPH particle data array
{
  int dn;                           // Integer time since beginning of step
  int i;                            // Particle counter
  int nstep;                        // Particle (integer) step size

  debug2("[EnergyPEC::EndTimestep]");

  // --------------------------------------------------------------------------
#pragma omp parallel for default(none) private(dn,i,nstep) \
      shared(n,Nsph,sphdata)
  for (i=0; i<Nsph; i++) {
    dn = n - sphdata[i].nlast;
    nstep = sphdata[i].nstep;
    if (dn == nstep) {
      sphdata[i].u0 = sphdata[i].u;
      sphdata[i].dudt0 = sphdata[i].dudt;
    }
  }
  // --------------------------------------------------------------------------

  return;
}



//=============================================================================
//  EnergyPEC::Timestep
/// Compute explicit timestep such that u cannot change by a large fraction 
/// in one step, i.e. dt = const*u/|dudt + epsilon| 
/// where epsilon is to prevent the denominator becoming zero.
//=============================================================================
template <int ndim>
DOUBLE EnergyPEC<ndim>::Timestep
(SphParticle<ndim> &part            ///< [inout] SPH particle reference
)
{
  return this->energy_mult*(DOUBLE) (part.u/(fabs(part.dudt) + small_number));
}



template class EnergyPEC<1>;
template class EnergyPEC<2>;
template class EnergyPEC<3>;



