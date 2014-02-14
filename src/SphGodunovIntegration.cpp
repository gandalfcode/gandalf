//=============================================================================
//  SphGodunovIntegration.cpp
//  Contains functions for integrating SPH particle positions and velocities 
//  using the conservative integration scheme for Godunov SPH described by 
//  Inutsuka (2002).
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
 FLOAT timestep,                    ///< [in] Minimum timestep level
 Sph<ndim> *sph)                    ///< [inout] Pointer to SPH object
{
  int dn;                           // Integer time since beginning of step
  int i;                            // Particle counter
  int k;                            // Dimension counter
  int nstep;                        // Particle (integer) step size
  FLOAT dt;                         // Timestep since start of step
  SphParticle<ndim> *part;          // SPH particle pointer

  debug2("[SphGodunovIntegration::AdvanceParticles]");

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
      (sph->sphintdata[i].v0[k] + 0.5*part->a[k]*part->dt)*dt;
    for (k=0; k<vdim; k++) part->v[k] = sph->sphintdata[i].v0[k] +
      sph->sphintdata[i].a0[k]*dt;

    // Set particle as active at end of step
    if (dn == nstep) part->active = true;
    else part->active = false;
  }
  //---------------------------------------------------------------------------

  return;
}
 


//=============================================================================
//  SphGodunovIntegration::CorrectionTerms
/// Empty definition (No correction terms to apply).
//=============================================================================
template <int ndim>
void SphGodunovIntegration<ndim>::CorrectionTerms
(int n,                             ///< [in] Integer time in block time struct
 FLOAT timestep,                    ///< [in] Base timestep value
 Sph<ndim> *sph)                    ///< [inout] Pointer to SPH object
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
 FLOAT timestep,                    ///< [in] Base timestep value
 Sph<ndim> *sph)                    ///< [inout] Pointer to SPH object
{
  int dn;                           // Integer time since beginning of step
  int i;                            // Particle counter
  int k;                            // Dimension counter
  int nstep;                        // Particle (integer) step size
  SphParticle<ndim> *part;          // SPH particle pointer

  debug2("[SphGodunovIntegration::EndTimestep]");

  //---------------------------------------------------------------------------
#pragma omp parallel for default(none) private(dn,i,k,nstep,part) shared(n,sph)
  for (i=0; i<sph->Nsph; i++) {
    dn = n - sph->sphintdata[i].nlast;
    nstep = sph->sphintdata[i].nstep;
    part = sph->sphintdata[i].part;

    if (n == sph->sphintdata[i].nlast) {
      for (k=0; k<ndim; k++) sph->sphintdata[i].r0[k] = part->r[k];
      for (k=0; k<ndim; k++) sph->sphintdata[i].v0[k] = part->v[k];
      for (k=0; k<ndim; k++) sph->sphintdata[i].a0[k] = part->a[k];
      sph->sphintdata[i].nlast = n;
      part->active = false;
      //part->active = true;
    }
  }
  //---------------------------------------------------------------------------

  return;
}



//=============================================================================
//  SphGodunovIntegration::CheckTimesteps
/// Record all important SPH particle quantities at the end of the step for  
/// the start of the new timestep.
//=============================================================================
template <int ndim>
int SphGodunovIntegration<ndim>::CheckTimesteps
(int level_diff_max,                ///< [in] Max. allowed SPH neib dt diff
 int n,                             ///< [in] Integer time in block time struct
 Sph<ndim> *sph)                    ///< [inout] Pointer to SPH object
{
  int activecount = 0;              // ..
  int dn;                           // Integer time since beginning of step
  int i;                            // Particle counter
  int k;                            // Dimension counter
  int level_new;                    // New level of particle
  int nnewstep;                     // New step size of particle
  int nstep;                        // Particle (integer) step size
  SphParticle<ndim> *part;          // SPH particle pointer

  debug2("[SphLeapfrogKDK::CheckTimesteps]");

  //---------------------------------------------------------------------------
#pragma omp parallel for default(none) shared(level_diff_max,n,sph)\
  private(dn,i,k,level_new,nnewstep,nstep,part) reduction(+:activecount)
  for (i=0; i<sph->Nsph; i++) {
    dn = n - sph->sphintdata[i].nlast;
    nstep = sph->sphintdata[i].nstep;
    if (dn == nstep) continue;

    // Check if neighbour timesteps are too small.  If so, then reduce 
    // timestep if possible
    if (part->levelneib - part->level > level_diff_max) {
      level_new = part->levelneib - level_diff_max;
      nnewstep = sph->sphintdata[i].nstep/pow(2,level_new - part->level);

      // If new level is correctly synchronised, then change all quantities
      if (n%nnewstep == 0) {
        nstep = dn;
        part->level = level_new;
        if (dn > 0) sph->sphintdata[i].nstep = dn; //nstep;
        part->active = true;
        activecount++;
      }
    }
  }
  //---------------------------------------------------------------------------

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
 Sph<ndim> *sph)                    ///< Pointer to SPH object
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
