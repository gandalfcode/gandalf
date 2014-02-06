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
                                     DOUBLE courant_mult_aux) :
  SphIntegration<ndim>(accel_mult_aux, courant_mult_aux)
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
 int Nsph,                          ///< [in] No. of SPH particles
 SphIntParticle<ndim> *sphintdata,  ///< [inout] SPH particle integration data
 FLOAT timestep)                    ///< [in] Base timestep value
{
  int dn;                           // Integer time since beginning of step
  int i;                            // Particle counter
  int k;                            // Dimension counter
  int nstep;                        // Particle (integer) step size
  FLOAT dt;                         // Timestep since start of step
  SphParticle<ndim> *part;          // SPH particle pointer

  debug2("[SphLeapfrogKDK::AdvanceParticles]");

  // Advance positions and velocities of all SPH particles
  //---------------------------------------------------------------------------
#pragma omp parallel for default(none) private(dn,dt,i,k,nstep,part)\
  shared(n,Nsph,sphintdata,timestep)
  for (i=0; i<Nsph; i++) {

    // Compute time since beginning of current step
    nstep = sphintdata[i].nstep;
    dn = n - sphintdata[i].nlast;
    dt = timestep*(FLOAT) dn;
    part = sphintdata[i].part;

    // Advance particle positions and velocities
    for (k=0; k<ndim; k++) part->r[k] = sphintdata[i].r0[k] +
      sphintdata[i].v0[k]*dt + 0.5*sphintdata[i].a0[k]*dt*dt;
    for (k=0; k<ndim; k++) part->v[k] =
      sphintdata[i].v0[k] + sphintdata[i].a0[k]*dt;

    // Set particle as active at end of step
    if (dn == nstep) part->active = true;
    else part->active = false;
  }
  //---------------------------------------------------------------------------

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
 int Nsph,                          ///< [in] No. of SPH particles
 SphIntParticle<ndim> *sphintdata,  ///< [inout] SPH particle integration data
 FLOAT timestep)                    ///< [in] Base timestep value
{
  int dn;                           // Integer time since beginning of step
  int i;                            // Particle counter
  int k;                            // Dimension counter
  int nstep;                        // Particle (integer) step size
  SphParticle<ndim> *part;          // SPH particle pointer

  debug2("[SphLeapfrogKDK::CorrectionTerms]");

  //---------------------------------------------------------------------------
#pragma omp parallel for default(none) private(dn,i,k,nstep,part)\
  shared(n,Nsph,sphintdata,timestep)
  for (i=0; i<Nsph; i++) {
    dn = n - sphintdata[i].nlast;
    nstep = sphintdata[i].nstep;
    part = sphintdata[i].part;
    if (dn == nstep)
      for (k=0; k<ndim; k++) part->v[k] += timestep*(FLOAT) nstep*
        (FLOAT) 0.5*(part->a[k] - sphintdata[i].a0[k]);
  }
  //---------------------------------------------------------------------------

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
 int Nsph,                          ///< [in] No. of SPH particles
 SphIntParticle<ndim> *sphintdata)  ///< [inout] SPH particle integration data
{
  int dn;                           // Integer time since beginning of step
  int i;                            // Particle counter
  int k;                            // Dimension counter
  int nstep;                        // Particle (integer) step size
  SphParticle<ndim> *part;          // SPH particle pointer

  debug2("[SphLeapfrogKDK::EndTimestep]");

  //---------------------------------------------------------------------------
#pragma omp parallel for default(none) private(dn,i,k,nstep,part)\
  shared(n,Nsph,sphintdata)
  for (i=0; i<Nsph; i++) {
    dn = n - sphintdata[i].nlast;
    nstep = sphintdata[i].nstep;
    part = sphintdata[i].part;

    if (dn == nstep) {
      for (k=0; k<ndim; k++) sphintdata[i].r0[k] = part->r[k];
      for (k=0; k<ndim; k++) sphintdata[i].v0[k] = part->v[k];
      for (k=0; k<ndim; k++) sphintdata[i].a0[k] = part->a[k];
      part->active = false;
      sphintdata[i].nlast = n;
    }
  }
  //---------------------------------------------------------------------------

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
 int Nsph,                          ///< [in] No. of SPH particles
 SphIntParticle<ndim> *sphintdata)  ///< [inout] SPH particle integration data
{
  int activecount = 0;              // No. of ptcls with new active timesteps
  int dn;                           // Integer time since beginning of step
  int i;                            // Particle counter
  int level_new;                    // New level of particle
  int nnewstep;                     // New step size of particle
  SphParticle<ndim> *part;          // SPH particle pointer

  debug2("[SphLeapfrogKDK::CheckTimesteps]");

  //---------------------------------------------------------------------------
#pragma omp parallel for default(none) private(dn,level_new,nnewstep,part)\
  shared(level_diff_max,n,Nsph,sphintdata) reduction(+:activecount)
  for (i=0; i<Nsph; i++) {
    dn = n - sphintdata[i].nlast;
    part = sphintdata[i].part;
    if (dn == sphintdata[i].nstep) continue;

    // Check if neighbour timesteps are too small.  If so, then reduce 
    // timestep if possible
    if (part->levelneib - part->level > level_diff_max) {
      level_new = part->levelneib - level_diff_max;
      nnewstep = sphintdata[i].nstep/pow(2,level_new - part->level);

      // If new level is correctly synchronised, then change all quantities
      if (dn%nnewstep == 0) {
        part->level = level_new;
        if (dn > 0) sphintdata[i].nstep = dn;
        part->active = true;
        activecount++;
      }
    }
  }
  //---------------------------------------------------------------------------

  return activecount;
}



// Template class instances for each dimensionality value (1, 2 and 3)
template class SphLeapfrogKDK<1>;
template class SphLeapfrogKDK<2>;
template class SphLeapfrogKDK<3>;
