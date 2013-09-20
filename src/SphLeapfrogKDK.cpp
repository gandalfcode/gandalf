//=============================================================================
//  SphLeapfrogKDK.cpp
//  Contains functions for integrating SPH particle positions and velocities 
//  using the leapfrog kick-drift-kick (KDK) scheme.
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
 SphParticle<ndim> *sphdata,        ///< [inout] SPH particle data array
 FLOAT timestep)                    ///< [in] Base timestep value
{
  int dn;                           // Integer time since beginning of step
  int i;                            // Particle counter
  int k;                            // Dimension counter
  int nstep;                        // Particle (integer) step size
  FLOAT dt;                         // Timestep since start of step

  debug2("[SphLeapfrogKDK::AdvanceParticles]");

  // Advance positions and velocities of all SPH particles
  // --------------------------------------------------------------------------
#pragma omp parallel for default(none) private(dt,k,nstep,i, dn)\
  shared(sphdata,Nsph,n,timestep)
  for (i=0; i<Nsph; i++) {

    // Compute time since beginning of current step
    nstep = sphdata[i].nstep;
    dn = n - sphdata[i].nlast;
    dt = timestep*(FLOAT) dn;

    // Advance particle positions and velocities
    for (k=0; k<ndim; k++) sphdata[i].r[k] = sphdata[i].r0[k] + 
      sphdata[i].v0[k]*dt + 0.5*sphdata[i].a0[k]*dt*dt;
    for (k=0; k<ndim; k++) sphdata[i].v[k] = 
      sphdata[i].v0[k] + sphdata[i].a0[k]*dt;

    // Set particle as active at end of step
    if (dn == nstep) sphdata[i].active = true;
    else sphdata[i].active = false;
  }
  // --------------------------------------------------------------------------

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
 SphParticle<ndim> *sphdata,        ///< [inout] SPH particle data array
 FLOAT timestep)                    ///< [in] Base timestep value
{
  int dn;                           // Integer time since beginning of step
  int i;                            // Particle counter
  int k;                            // Dimension counter
  int nstep;                        // Particle (integer) step size

  debug2("[SphLeapfrogKDK::CorrectionTerms]");

  // --------------------------------------------------------------------------
#pragma omp parallel for default(none) private(dn,k,nstep,i)\
  shared(n,Nsph,sphdata,timestep)
  for (i=0; i<Nsph; i++) {
    dn = n - sphdata[i].nlast;
    nstep = sphdata[i].nstep;
    if (dn == nstep)
      for (k=0; k<ndim; k++) sphdata[i].v[k] += timestep*(FLOAT) nstep*
        (FLOAT) 0.5*(sphdata[i].a[k] - sphdata[i].a0[k]);
  }
  // --------------------------------------------------------------------------

  return;
}



//=============================================================================
//  SphLeapfrogKDK::EndTimestep
/// Record all important SPH particle quantities at the end of the step for  
/// the start of the new timestep.
// ============================================================================
template <int ndim>
void SphLeapfrogKDK<ndim>::EndTimestep
(int n,                             ///< [in] Integer time in block time struct
 int Nsph,                          ///< [in] No. of SPH particles
 SphParticle<ndim> *sphdata)        ///< [inout] SPH particle data array
{
  int dn;                           // Integer time since beginning of step
  int i;                            // Particle counter
  int k;                            // Dimension counter
  int nstep;                        // Particle (integer) step size

  debug2("[SphLeapfrogKDK::EndTimestep]");

  // --------------------------------------------------------------------------
#pragma omp parallel for default(none) private(dn,k,nstep,i)\
  shared(n,Nsph,sphdata)
  for (i=0; i<Nsph; i++) {
    dn = n - sphdata[i].nlast;
    nstep = sphdata[i].nstep;

    if (dn == nstep) {
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
//  SphLeapfrogKDK::CheckTimesteps
/// Record all important SPH particle quantities at the end of the step for  
/// the start of the new timestep.
// ============================================================================
template <int ndim>
int SphLeapfrogKDK<ndim>::CheckTimesteps
(int level_diff_max,                ///< [in] Max. allowed SPH neib dt diff
 int n,                             ///< [in] Integer time in block time struct
 int Nsph,                          ///< [in] No. of SPH particles
 SphParticle<ndim> *sphdata)        ///< [inout] SPH particle data array
{
  int activecount = 0;              // ..
  int dn;                           // Integer time since beginning of step
  int i;                            // Particle counter
  int k;                            // Dimension counter
  int level_new;                    // New level of particle
  int nnewstep;                     // New step size of particle
  int nstep;                        // Particle (integer) step size

  debug2("[SphLeapfrogKDK::CheckTimesteps]");

  // --------------------------------------------------------------------------
#pragma omp parallel for default(none) private(dn,k,level_new,nnewstep,nstep)\
  shared(level_diff_max,n,Nsph,sphdata) reduction(+:activecount)
  for (i=0; i<Nsph; i++) {
    dn = n - sphdata[i].nlast;
    nstep = sphdata[i].nstep;
    if (dn == nstep) continue;

    // Check if neighbour timesteps are too small.  If so, then reduce 
    // timestep if possible
    if (sphdata[i].levelneib - sphdata[i].level > level_diff_max) {
      level_new = sphdata[i].levelneib - level_diff_max;
      nnewstep = sphdata[i].nstep/pow(2,level_new - sphdata[i].level);

      // If new level is correctly synchronised, then change all quantities
      if (n%nnewstep == 0) {
        nstep = dn;
        sphdata[i].level = level_new;
        if (dn > 0) sphdata[i].nstep = dn; //nstep;
        sphdata[i].active = true;
        activecount++;
      }
    }
  }
  // --------------------------------------------------------------------------

  return activecount;
}



// Template class instances for each dimensionality value (1, 2 and 3)
template class SphLeapfrogKDK<1>;
template class SphLeapfrogKDK<2>;
template class SphLeapfrogKDK<3>;
