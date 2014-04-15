//=============================================================================
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
//  SphLeapfrogDKD::SphLeapfrogDKD
/// SphLeapfrogDKD class constructor
//=============================================================================
template <int ndim, template <int> class ParticleType>
SphLeapfrogDKD<ndim, ParticleType>::SphLeapfrogDKD(DOUBLE accel_mult_aux,
                                     DOUBLE courant_mult_aux,
                                     DOUBLE energy_mult_aux,
				     eosenum gas_eos_aux,
				     tdaviscenum tdavisc_aux) :
  SphIntegration<ndim>(accel_mult_aux, courant_mult_aux, energy_mult_aux,
		       gas_eos_aux, tdavisc_aux)
{
}



//============================================================================
//  SphLeapfrogDKD::~SphLeapfrog()
/// SphLeapfrogDKD class destructor
//=============================================================================
template <int ndim, template <int> class ParticleType>
SphLeapfrogDKD<ndim, ParticleType>::~SphLeapfrogDKD()
{
}



//=============================================================================
//  SphLeapfrogDKD::AdvanceParticles
/// Integrate both particle positions and velocities to 1st order from the 
/// beginning of the step to the current simulation time, i.e. 
/// $r(t+dt) = r(t) + v(t)*dt$,
/// $v(t+dt) = v(t) + a(t)*dt$.
/// If particle has reached the half-step, then set as active to compute 
/// the 'kick' acceleration terms.
//=============================================================================
template <int ndim, template <int> class ParticleType>
void SphLeapfrogDKD<ndim, ParticleType>::AdvanceParticles
(int n,                             ///< [in] Integer time in block time struct
 FLOAT timestep,                    ///< [in] Base timestep value
 int npart,                         ///< [in] Number of particles
 SphParticle<ndim>* sph_gen)        ///< [inout] Pointer to SPH particle array
{
  int dn;                           // Integer time since beginning of step
  int i;                            // Particle counter
  int k;                            // Dimension counter
  int nstep;                        // Particle (integer) step size
  FLOAT dt;                         // Timestep since start of step

  debug2("[SphLeapfrogDKD::AdvanceParticles]");

  ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sph_gen);

  // Advance positions and velocities of all SPH particles
  //---------------------------------------------------------------------------
#pragma omp parallel for default(none) private(dn,dt,i,k,nstep)\
  shared(n,sphdata,timestep,npart)
  for (i=0; i<npart; i++) {
    SphParticle<ndim>& part = sphdata[i];

    // Compute time since beginning of current step
    nstep = part.nstep;
    dn = n - part.nlast;
    dt = timestep*(FLOAT) dn;

    // Advance particle positions and velocities
    for (k=0; k<ndim; k++) part.r[k] =
      part.r0[k] + part.v0[k]*dt;
    for (k=0; k<ndim; k++) part.v[k] =
      part.v0[k] + part.a[k]*dt;

    // Set particle as active at end of step
    if (dn == nstep/2) part.active = true;
    else part.active = false;
  }
  //---------------------------------------------------------------------------

  return;
}
 


//=============================================================================
//  SphLeapfrogDKD::CorrectionTerms
/// Empty function.  No correction terms for Leapfrog drift-kick-drift scheme.
//=============================================================================
template <int ndim, template <int> class ParticleType>
void SphLeapfrogDKD<ndim, ParticleType>::CorrectionTerms
(int n,                             ///< [in] Integer time in block time struct
 FLOAT timestep,                    ///< [in] Base timestep value
 int npart,                         ///< [in] Number of particles
 SphParticle<ndim>* sph_gen)        ///< [inout] Pointer to SPH particle array
{
  return;
}



//=============================================================================
//  SphLeapfrogDKD::EndTimestep
/// Record all important SPH particle quantities at the end of the step for  
/// the start of the new timestep.
//=============================================================================
template <int ndim, template <int> class ParticleType>
void SphLeapfrogDKD<ndim, ParticleType>::EndTimestep
(int n,                             ///< [in] Integer time in block time struct
 FLOAT timestep,                    ///< [in] Base timestep value
 int npart,                         ///< [in] Number of particles
 SphParticle<ndim>* sph_gen)        ///< [inout] Pointer to SPH particle array
{
  int dn;                           // Integer time since beginning of step
  int i;                            // Particle counter
  int k;                            // Dimension counter
  int nstep;                        // Particle (integer) step size

  ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sph_gen);

  debug2("[SphLeapfrogDKD::EndTimestep]");

  //---------------------------------------------------------------------------
#pragma omp parallel for default(none) private(dn,i,k,nstep) shared(n,sphdata,npart)
  for (i=0; i<npart; i++) {
    SphParticle<ndim>& part = sphdata[i];

    dn = n - part.nlast;
    nstep = part.nstep;

    if (dn == nstep) {
      for (k=0; k<ndim; k++) part.r0[k] = part.r[k];
      for (k=0; k<ndim; k++) part.v0[k] = part.v[k];
      for (k=0; k<ndim; k++) part.a0[k] = part.a[k];
      part.nlast = n;
      part.active = false;
    }
  }
  //---------------------------------------------------------------------------

  return;
}



//=============================================================================
//  SphLeapfrogDKD::CheckTimesteps
/// Record all important SPH particle quantities at the end of the step for  
/// the start of the new timestep.
//=============================================================================
template <int ndim, template <int> class ParticleType>
int SphLeapfrogDKD<ndim, ParticleType>::CheckTimesteps
(int level_diff_max,                ///< [in] Max. allowed SPH neib dt diff
 int n,                             ///< [in] Integer time in block time struct
 int npart,                         ///< [in] Number of particles
 SphParticle<ndim>* sph_gen)        ///< [inout] Pointer to SPH particle array
{
  int activecount = 0;              // ..
  int dn;                           // Integer time since beginning of step
  int i;                            // Particle counter
  int k;                            // Dimension counter
  int level_new;                    // ..
  int nnewstep;                     // ..
  int nstep;                        // Particle (integer) step size

  debug2("[SphLeapfrogDKD::CheckTimesteps]");

  ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sph_gen);

  //---------------------------------------------------------------------------
#pragma omp parallel for default(none) reduction(+:activecount)\
  private(dn,i,k,level_new,nnewstep,nstep)\
  shared(level_diff_max,n,sphdata,npart)
  for (i=0; i<npart; i++) {
    SphParticle<ndim>& part = sphdata[i];

    dn = n - part.nlast;
    nstep = part.nstep;

    // Check if neighbour timesteps are too small.  If so, then reduce 
    // timestep if possible
    if (part.levelneib - part.level > level_diff_max) {
      level_new = part.levelneib - level_diff_max;
      nnewstep = part.nstep/pow(2,level_new - part.level);

      // If new level is correctly synchronised, then change all quantities
      if (n%nnewstep == 0) {
        nstep = dn;
        part.level = level_new;
        if (dn > 0) part.nstep = dn; //nstep;
        if (dn == nnewstep/2) part.active = true;
        activecount++;
      }
    }
  }
  //---------------------------------------------------------------------------

  return activecount;
}



// Template class instances for each dimensionality value (1, 2 and 3)
template class SphLeapfrogDKD<1, GradhSphParticle>;
template class SphLeapfrogDKD<2, GradhSphParticle>;
template class SphLeapfrogDKD<3, GradhSphParticle>;
template class SphLeapfrogDKD<1, SM2012SphParticle>;
template class SphLeapfrogDKD<2, SM2012SphParticle>;
template class SphLeapfrogDKD<3, SM2012SphParticle>;
