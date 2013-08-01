//=============================================================================
//  SphGodunovIntegration.cpp
//  Contains functions for integrating SPH particle positions and velocities 
//  using the conservative integration scheme for Godunov SPH described by 
//  Inutsuka (2002).
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
#include "Debug.h"
using namespace std;



//=============================================================================
//  SphGodunovIntegration::SphGodunovIntegration
/// SphGodunovIntegration class constructor
//=============================================================================
template <int ndim>
SphGodunovIntegration<ndim>::SphGodunovIntegration(
			       DOUBLE accel_mult_aux, 
			       DOUBLE courant_mult_aux) :
  SphIntegration<ndim>(accel_mult_aux, courant_mult_aux)
{
}



//=============================================================================
//  SphGodunovIntegration::~SphGodunovIntegration()
/// SphGodunovIntegration class destructor
//=============================================================================
template <int ndim>
SphGodunovIntegration<ndim>::~SphGodunovIntegration()
{
}



//=============================================================================
//  SphGodunovIntegration::AdvanceParticles
/// Integrate particle positions to 2nd order, and particle velocities to 1st
/// order from the beginning of the step to the current simulation time, i.e. 
/// $r(t+dt) = r(t) + v(t)*dt + 0.5*a(t)*dt^2$, 
/// $v(t+dt) = v(t) + a(t)*dt$.
/// Also set particles at the end of step as 'active' in order to compute 
/// the end-of-step force computation.
//=============================================================================
template <int ndim>
void SphGodunovIntegration<ndim>::AdvanceParticles
(int n,                             ///< [in] Current integer time
 int Nsph,                          ///< [in] No. of SPH particles
 SphParticle<ndim> *sphdata,        ///< [inout] SPH particle data array
 FLOAT timestep)                    ///< [in] Minimum timestep level
{
  int dn;                           // Integer time since beginning of step
  int i;                            // Particle counter
  int k;                            // Dimension counter
  int nstep;                        // Particle (integer) step size
  FLOAT dt;                         // Timestep since start of step

  debug2("[SphGodunovIntegration::AdvanceParticles]");

  // --------------------------------------------------------------------------
#pragma omp parallel for default(none) private(dt,k,nstep,i,dn)\
  shared(n,Nsph,sphdata,timestep)
  for (i=0; i<Nsph; i++) {

    // Compute time since beginning of current step
    nstep = sphdata[i].nstep;
    dn = n - sphdata[i].nlast;
    dt = timestep*(FLOAT) dn;

    // Advance particle positions and velocities
    for (k=0; k<ndim; k++) sphdata[i].r[k] = sphdata[i].r0[k] + 
      (sphdata[i].v0[k] + 0.5*sphdata[i].a[k]*sphdata[i].dt)*dt;
    for (k=0; k<vdim; k++) sphdata[i].v[k] = sphdata[i].v0[k] + 
      sphdata[i].a0[k]*dt;

    // Set particle as active at end of step
    if (dn == nstep) sphdata[i].active = true;
    else sphdata[i].active = false;
  }
  // --------------------------------------------------------------------------

  return;
}
 


//=============================================================================
//  SphGodunovIntegration::CorrectionTerms
/// Empty definition (No correction terms to apply).
//=============================================================================
template <int ndim>
void SphGodunovIntegration<ndim>::CorrectionTerms
(int n, int Nsph, SphParticle<ndim> *sph, FLOAT timestep)
{
  return;
}



//=============================================================================
//  SphGodunovIntegration::EndTimestep
/// Record all important SPH particle quantities at the end of the step for 
/// the start of the new timestep.
//=============================================================================
template <int ndim>
void SphGodunovIntegration<ndim>::EndTimestep
(int n,                             ///< [in] Integer time in block time struct
 int Nsph,                          ///< [in] No. of SPH particles
 SphParticle<ndim> *sphdata)        ///< [inout] SPH particle data array
{
  int dn;                           // Integer time since beginning of step
  int i;                            // Particle counter
  int k;                            // Dimension counter
  int nstep;                        // Particle (integer) step size

  debug2("[SphGodunovIntegration::EndTimestep]");

  // --------------------------------------------------------------------------
#pragma omp parallel for default(none) private(k,nstep,i,dn)\
  shared(n,Nsph,sphdata)
  for (i=0; i<Nsph; i++) {
    dn = n - sphdata[i].nlast;
    nstep = sphdata[i].nstep;

    if (n == sphdata[i].nlast) {
      for (k=0; k<ndim; k++) sphdata[i].r0[k] = sphdata[i].r[k];
      for (k=0; k<ndim; k++) sphdata[i].v0[k] = sphdata[i].v[k];
      for (k=0; k<ndim; k++) sphdata[i].a0[k] = sphdata[i].a[k];
      sphdata[i].active = false;
      sphdata[i].nlast = n;
    }
  }
  // --------------------------------------------------------------------------

  return;
}



//=============================================================================
//  SphGodunovIntegration::CheckTimesteps
/// Record all important SPH particle quantities at the end of the step for  
/// the start of the new timestep.
// ============================================================================
template <int ndim>
int SphGodunovIntegration<ndim>::CheckTimesteps
(int level_diff_max, int n,                             ///< [in] Integer time in block time struct
 int Nsph,                          ///< [in] No. of SPH particles
 SphParticle<ndim> *sphdata)        ///< [inout] SPH particle data array
{
  int activecount = 0;              // ..
  int dn;                           // Integer time since beginning of step
  int i;                            // Particle counter
  int k;                            // Dimension counter
  int level_new;                    // ..
  int nnewstep;                     // ..
  int nstep;                        // Particle (integer) step size

  debug2("[SphLeapfrogKDK::CheckTimesteps]");

  // --------------------------------------------------------------------------
#pragma omp parallel for default(none) private(dn,k,level_new,nnewstep,nstep)\
  shared(level_diff_max,n,Nsph,sphdata) reduction(+:activecount)
  for (i=0; i<Nsph; i++) {
    dn = n - sphdata[i].nlast;
    nstep = sphdata[i].nstep;

    // Check if neighbour timesteps are too small.  If so, then reduce 
    // timestep if possible
    if (sphdata[i].levelneib - sphdata[i].level > level_diff_max) {
      level_new = sphdata[i].levelneib - level_diff_max;
      nnewstep = sphdata[i].nstep/pow(2,level_new - sphdata[i].level);

      // If new level is correctly synchronised, then change all quantities
      if (n%nnewstep == 0) {
	nstep = dn;
	sphdata[i].level = level_new;
	if (dn > 0) sphdata[i].nstep = nstep;
	sphdata[i].active = true;
	activecount++;
      }
    }
  }
  // --------------------------------------------------------------------------

  return activecount;
}



//=============================================================================
//  SphIntegration::Timestep
/// Default timestep size for SPH particles.  Takes the minimum of : 
/// (i)  const*h/(sound_speed + h*|div_v|)    (Courant condition)
/// (ii) const*sqrt(h/|a|)                    (Acceleration condition)
//=============================================================================
template <int ndim>
DOUBLE SphGodunovIntegration<ndim>::Timestep
(SphParticle<ndim> &part,           ///< Reference to SPH particle
 int hydro_forces)                  ///< Computing hydro forces or not
{
  DOUBLE timestep;                  // Variable to record/compare timesteps
  //DOUBLE amag;

  // Courant condition
  timestep = this->courant_mult*part.h/(part.sound + small_number_dp);

  // Local convergence/divergence condition
  timestep = min(timestep,this->courant_mult/
                 (fabs(part.div_v) + small_number_dp));

  //Acceleration condition
  //amag = sqrt(DotProduct(part.a,part.a,ndim));
  //timestep = min(timestep, accel_mult*sqrt(part.h/(amag + small_number_dp)));

  return timestep;
}



template class SphGodunovIntegration<1>;
template class SphGodunovIntegration<2>;
template class SphGodunovIntegration<3>;
