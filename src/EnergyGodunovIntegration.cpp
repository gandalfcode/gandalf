//=============================================================================
//  EnergyGodunovIntegration.cpp
//  Energy integration scheme for Inutsuka (2002) Godunov SPH algorithm.
//  Conserves total energy (thermal + kinetic) to machine precision for
//  direct sum forces and global timesteps.
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
#include "SmoothingKernel.h"
#include "SphIntegration.h"
#include "Particle.h"
#include "EOS.h"
#include "EnergyEquation.h"
#include "Exception.h"
#include "Debug.h"
using namespace std;




//=============================================================================
//  EnergyGodunovIntegration::EnergyGodunovIntegration()
/// EnergyGodunovIntegration class constructor
//=============================================================================
template <int ndim>
EnergyGodunovIntegration<ndim>::EnergyGodunovIntegration(DOUBLE energy_mult_aux) :
  EnergyEquation<ndim>(energy_mult_aux)
{
}



//=============================================================================
//  EnergyGodunovIntegration::~EnergyGodunovIntegration()
/// EnergyGodunovIntegration class destructor
//=============================================================================
template <int ndim>
EnergyGodunovIntegration<ndim>::~EnergyGodunovIntegration()
{
}



//=============================================================================
//  EnergyGodunovIntegration::EnergyIntegration
/// Integrate internal energy to first order from the beginning of the step to
/// the current simulation time, i.e. $u(t+dt) = u(t) + dudt(t)*dt$
//=============================================================================
template <int ndim>
void EnergyGodunovIntegration<ndim>::EnergyIntegration
(int n,                             ///< [in] Integer time in block time struct
 int Nhydro,                          ///< [in] No. of SPH particles
 SphParticle<ndim> *sphdata_gen,  ///< [inout] SPH particle data array
 FLOAT timestep)                    ///< [in] Base timestep value
{
  int dn;                           // Integer time since beginning of step
  int i;                            // Particle counter
  int nstep;                        // Particle (integer) step size
  FLOAT dt;                         // Timestep since start of step

  debug2("[EnergyGodunovIntegration::EnergyIntegration]");

  GodunovSphParticle<ndim>* sphdata = static_cast<GodunovSphParticle<ndim>* > (sphdata_gen);

  //---------------------------------------------------------------------------
#pragma omp parallel for default(none) private(dn,dt,i,nstep) \
  shared(n,Nhydro,timestep,cout,sphdata)
  for (i=0; i<Nhydro; i++) {
    SphParticle<ndim>& part = sphdata[i];

    nstep = sphdata[i].nstep;
    dn = n - sphdata[i].nlast;
    dt = timestep*(FLOAT) dn;
    part.u = sphdata[i].u0 + part.dudt*dt;

    //if (i == 0) {
    //  cout << "Energy int : " << i << "   " << sphdata[i].u0 << "    " << part.u << "    " << part.dudt << "    " << dt << endl;
    //}

    if (part.u != part.u) {
      cout << "Something wrong with energy integration (NaN) : " << endl;
      cout << part.u << "   " << sphdata[i].u0 << "  " << part.dudt
	   << "   " << dt << "   " << nstep << "    " << timestep << endl;
      exit(0);
    }
    if (part.u < small_number) {
      cout << "Something wrong with energy integration (0) : " << endl;
      cout << part.u << "   " << sphdata[i].u0 << "  " << part.dudt
	   << "   " << dt << "   " << nstep << "    "
	   << sphdata[i].u0/part.dudt << endl;
      cout << " dt_courant : " << part.h/part.sound << "    "
	   << sphdata[i].u0/(part.dudt + small_number) << "    "
	   << timestep << endl;
      string message = "Problem with energy integration (0)";
      ExceptionHandler::getIstance().raise(message);
    }
  }
  //---------------------------------------------------------------------------

  return;
}




//=============================================================================
//  EnergyGodunovIntegration::EndTimestep
/// Record all important thermal quantities at the end of the step for the
/// start of the new timestep.
//=============================================================================
template <int ndim>
void EnergyGodunovIntegration<ndim>::EndTimestep
(int n,                             ///< [in] Integer time in block time struct
 int Nhydro,                          ///< [in] No. of SPH particles
 FLOAT timestep,                    ///< [in] Base timestep value
 SphParticle<ndim> *sphdata)  ///< [inout] SPH particle data array
{
  int dn;                           // Integer time since beginning of step
  int i;                            // Particle counter
  int nstep;                        // Particle (integer) step size

  debug2("[EnergyGodunovIntegration::EndTimestep]");

  //---------------------------------------------------------------------------
#pragma omp parallel for default(none) private(dn,i,nstep)	\
  shared(n,Nhydro,sphdata)
  for (i=0; i<Nhydro; i++) {
    SphParticle<ndim>& part = sphdata[i];
    dn = n - sphdata[i].nlast;
    nstep = sphdata[i].nstep;
    if (n == sphdata[i].nlast) {
      //if (n == 0 || n - sphdata[i].nlast == sphdata[i].nstep) {
      sphdata[i].u0 = part.u;
      sphdata[i].dudt0 = part.dudt;
    }
  }
  //---------------------------------------------------------------------------

  return;
}



//=============================================================================
//  EnergyGodunovIntegration::Timestep
/// Compute explicit timestep such that u cannot change by a large fraction
/// in one step, i.e. dt = const*u/|dudt + epsilon| where epsilon is a small
/// number to prevent the denominator becoming zero.
//=============================================================================
template <int ndim>
DOUBLE EnergyGodunovIntegration<ndim>::Timestep(SphParticle<ndim> &part)
{
  return this->energy_mult*(DOUBLE) (part.u/(fabs(part.dudt) + small_number));
}



template class EnergyGodunovIntegration<1>;
template class EnergyGodunovIntegration<2>;
template class EnergyGodunovIntegration<3>;
