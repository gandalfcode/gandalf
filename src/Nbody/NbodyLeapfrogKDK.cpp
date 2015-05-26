//=================================================================================================
//  NbodyLeapfrogKDK.cpp
//  Contains functions for integrating star particle positions and velocities
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
#include "Precision.h"
#include "NbodyParticle.h"
#include "StarParticle.h"
#include "Parameters.h"
#include "Nbody.h"
#include "SmoothingKernel.h"
#include "Debug.h"
#include "Exception.h"
#include "InlineFuncs.h"
#include "Sph.h"
using namespace std;



//=================================================================================================
//  NbodyLeapfrogKDK::NbodyLeapfrogKDK()
/// N-body leapfrog KDK class constructor
//=================================================================================================
template <int ndim, template<int> class kernelclass>
NbodyLeapfrogKDK<ndim, kernelclass>::NbodyLeapfrogKDK
(int nbody_softening_aux, int sub_systems_aux,
 DOUBLE nbody_mult_aux, string KernelName) :
  Nbody<ndim>(nbody_softening_aux, sub_systems_aux, nbody_mult_aux, KernelName, 1),
  kern(kernelclass<ndim>(KernelName))
{
  this->kernp = &kern;
}



//=================================================================================================
//  NbodyLeapfrogKDK::~NbodyLeapfrog()
/// N-body leapfrog KDK class destructor
//=================================================================================================
template <int ndim, template<int> class kernelclass>
NbodyLeapfrogKDK<ndim, kernelclass>::~NbodyLeapfrogKDK()
{
  //DeallocateMemory();
}



//=================================================================================================
//  NbodyLeapfrogKDK::CalculateDirectSmoothedGravForces
/// Calculate all star-star force contributions for active systems using
/// direct summation with unsoftened gravity.
//=================================================================================================
template <int ndim, template<int> class kernelclass>
void NbodyLeapfrogKDK<ndim, kernelclass>::CalculateDirectSmoothedGravForces
 (int N,                               ///< Number of stars
  NbodyParticle<ndim> **star)          ///< Array of stars/systems
{
  int i,j,k;                           // Star and dimension counters
  DOUBLE dr[ndim];                     // Relative position vector
  DOUBLE drdt;                         // Rate of change of distance
  DOUBLE drmag;                        // ..
  DOUBLE drsqd;                        // Distance squared
  DOUBLE dv[ndim];                     // Relative velocity vector
  DOUBLE invdrmag;                     // 1 / drmag
  DOUBLE invhmean;                     // ..
  DOUBLE paux;                         // ..
  DOUBLE wmean;                        // ..

  debug2("[NbodyLeapfrogKDK::CalculateDirectSmoothedGravForces]");

  // Loop over all (active) stars
  //-----------------------------------------------------------------------------------------------
  for (i=0; i<N; i++) {
    if (star[i]->active == 0) continue;

    // Sum grav. contributions for all other stars (excluding star itself)
    //---------------------------------------------------------------------------------------------
    for (j=0; j<N; j++) {
      if (i == j) continue;

      for (k=0; k<ndim; k++) dr[k] = star[j]->r[k] - star[i]->r[k];
      for (k=0; k<ndim; k++) dv[k] = star[j]->v[k] - star[i]->v[k];
      drsqd    = DotProduct(dr,dr,ndim);
      drmag    = sqrt(drsqd) + small_number;
      invdrmag = 1.0/drmag;
      invhmean = 2.0/(star[i]->h + star[j]->h);
      drdt     = DotProduct(dv,dr,ndim)*invdrmag;
      paux     = star[j]->m*invhmean*invhmean*kern.wgrav(drmag*invhmean)*invdrmag;
      wmean    = kern.w0(drmag*invhmean)*powf(invhmean,ndim);

      // Add contribution to main star array
      star[i]->gpot += star[j]->m*invhmean*kern.wpot(drmag*invhmean);
      for (k=0; k<ndim; k++) star[i]->a[k] += paux*dr[k];
      for (k=0; k<ndim; k++) star[i]->adot[k] += paux*dv[k] -
        3.0*paux*drdt*invdrmag*dr[k] + 2.0*twopi*star[j]->m*drdt*wmean*invdrmag*dr[k];

    }
    //---------------------------------------------------------------------------------------------

  }
  //-----------------------------------------------------------------------------------------------

  return;
}



//=================================================================================================
//  NbodyLeapfrogKDK::CalculateDirectSPHForces
/// Calculate all forces on star due to neighbouring and direct-sum gas ptcls.
//=================================================================================================
template <int ndim, template<int> class kernelclass>
void NbodyLeapfrogKDK<ndim, kernelclass>::CalculateDirectSPHForces
 (NbodyParticle<ndim> *star,           ///< [inout] Pointer to star
  int Nhydro,                            ///< [in] Number of SPH gas particles
  int Ndirect,                         ///< [in] No. of direct sum gas particles
  int *sphlist,                        ///< [in] List of SPH neib. particles
  int *directlist,                     ///< [in] List of direct cum gas particles
  Sph<ndim> *sph)                      ///< [in] Array of SPH particles
{
  int j,jj,k;                          // SPH particle and dimension counters
  DOUBLE dr[ndim];                     // Relative position vector
  DOUBLE drmag;                        // Distance
  DOUBLE drsqd;                        // Distance squared
  DOUBLE invhmean;                     // 1 / hmean
  DOUBLE invdrmag;                     // 1 / drmag
  DOUBLE paux;                         // Aux. force variable

  debug2("[NbodyLeapfrogKDK::CalculateDirectSPHForces]");


  // Sum grav. contributions from all neighbouring SPH particles
  //-----------------------------------------------------------------------------------------------
  for (jj=0; jj<Nhydro; jj++) {

    j = sphlist[jj];
    SphParticle<ndim>& part = sph->GetSphParticlePointer(j);

    for (k=0; k<ndim; k++) dr[k] = part.r[k] - star->r[k];
    drsqd    = DotProduct(dr,dr,ndim);
    drmag    = sqrt(drsqd) + small_number;
    invdrmag = 1.0/drmag;
    invhmean = 2.0/(star->h + part.h);
    paux     = part.m*invhmean*invhmean*kern.wgrav(drmag*invhmean)*invdrmag;

    // Add contribution to main star array
    for (k=0; k<ndim; k++) star->a[k] += paux*dr[k];
    star->gpot += part.m*invhmean*kern.wpot(drmag*invhmean);

  }
  //-----------------------------------------------------------------------------------------------


  // Now include contributions from distant, non-SPH neighbours
  // (i.e. direct summation with Newton's law of gravity)
  //-----------------------------------------------------------------------------------------------
  for (jj=0; jj<Ndirect; jj++) {

    j = directlist[jj];
    SphParticle<ndim>& part = sph->GetSphParticlePointer(j);

    for (k=0; k<ndim; k++) dr[k] = part.r[k] - star->r[k];
    drsqd    = DotProduct(dr,dr,ndim);
    drmag    = sqrt(drsqd) + small_number;
    invdrmag = 1.0/drmag;
    paux     = part.m*pow(invdrmag,3);

    // Add contribution to main star array
    for (k=0; k<ndim; k++) star->a[k] += paux*dr[k];
    star->gpot += part.m*invdrmag;

  }
  //-----------------------------------------------------------------------------------------------


  return;
}



//=================================================================================================
//  NbodyLeapfrogKDK::AdvanceParticles
/// Integrate star positions to 2nd order, and star velocities to 1st order
/// from the beginning of the step to the current simulation time, i.e.
/// $r(t+dt) = r(t) + v(t)*dt + 0.5*a(t)*dt^2$,
/// $v(t+dt) = v(t) + a(t)*dt$.
/// Also set particles at the end of step as 'active' in order to compute
/// the end-of-step force computation and velocity correction step.
//=================================================================================================
template <int ndim, template<int> class kernelclass>
void NbodyLeapfrogKDK<ndim, kernelclass>::AdvanceParticles
 (int n,                               ///< Integer time
  int N,                               ///< No. of stars/systems
  DOUBLE t,                            ///< Current time
  DOUBLE timestep,                     ///< Smallest timestep value
  NbodyParticle<ndim> **star)          ///< Main star/system array
{
  int dn;                              // Integer time since beginning of step
  int i;                               // Particle counter
  int k;                               // Dimension counter
  int nstep;                           // Particle (integer) step size
  DOUBLE dt;                           // Timestep since start of step

  debug2("[NbodyLeapfrogKDK::AdvanceParticles]");
  timing->StartTimingSection("NBODY_ADVANCE");

  // Advance positions and velocities of all star particles
  //-----------------------------------------------------------------------------------------------
  for (i=0; i<N; i++) {

    // Compute time since beginning of step
    nstep = star[i]->nstep;
    dn = n - star[i]->nlast;
    //dt = timestep*(FLOAT) dn;
    dt = t - star[i]->tlast;

    // Advance positions to second order and velocities to first order
    for (k=0; k<ndim; k++) star[i]->r[k] = star[i]->r0[k] +
      star[i]->v0[k]*dt + 0.5*star[i]->a0[k]*dt*dt;
    for (k=0; k<vdim; k++) star[i]->v[k] = star[i]->v0[k] + star[i]->a0[k]*dt;

    // If at end of step, set system particle as active
    if (dn == nstep) star[i]->active = true;
    else star[i]->active = false;
  }
  //-----------------------------------------------------------------------------------------------

  timing->EndTimingSection("NBODY_ADVANCE");

  return;
}



//=================================================================================================
//  NbodyLeapfrogKDK::CorrectionTerms
/// Compute velocity integration to second order at the end of the step by
/// adding a second order correction term,
/// $v(t+dt) -> v(t+dt) + 0.5*(a(t+dt) - a(t))*dt$.
/// The full integration therefore becomes
/// $v(t+dt) = v(t) + 0.5*(a(t) + a(t+dt))*dt$.
//=================================================================================================
template <int ndim, template<int> class kernelclass>
void NbodyLeapfrogKDK<ndim, kernelclass>::CorrectionTerms
 (int n,                               ///< Integer time
  int N,                               ///< No. of stars/systems
  DOUBLE t,                            ///< Current time
  DOUBLE timestep,                     ///< Smallest timestep value
  NbodyParticle<ndim> **star)          ///< Main star/system array
{
  int dn;                              // Integer time since beginning of step
  int i;                               // Particle counter
  int k;                               // Dimension counter
  int nstep;                           // Particle (integer) step size

  debug2("[NbodyLeapfrogKDK::CorrectionTerms]");
  timing->StartTimingSection("NBODY_CORRECTION_TERMS");

  // Loop over all star particles and calculate correction terms only for
  // those at end of step.
  //-----------------------------------------------------------------------------------------------
  for (i=0; i<N; i++) {
    dn = n - star[i]->nlast;
    nstep = star[i]->nstep;
    if (dn == nstep) {
      for (k=0; k<ndim; k++) star[i]->v[k] +=
        0.5*(star[i]->a[k] - star[i]->a0[k])*(t - star[i]->tlast);
    }
  }
  //-----------------------------------------------------------------------------------------------

  timing->EndTimingSection("NBODY_CORRECTION_TERMS");

  return;
}



//=================================================================================================
//  NbodyLeapfrogKDK::EndTimestep
/// Record all important star particle quantities at the end of the step
/// for the start of the new timestep.
//=================================================================================================
template <int ndim, template<int> class kernelclass>
void NbodyLeapfrogKDK<ndim, kernelclass>::EndTimestep
 (int n,                               ///< Integer time
  int N,                               ///< No. of stars/systems
  DOUBLE t,                            ///< Current time
  DOUBLE timestep,                     ///< Smallest timestep value
  NbodyParticle<ndim> **star)          ///< Main star/system array
{
  int dn;                              // Integer time since beginning of step
  int i;                               // Particle counter
  int k;                               // Dimension counter
  int nstep;                           // Particle (integer) step size

  debug2("[NbodyLeapfrogKDK::EndTimestep]");
  timing->StartTimingSection("NBODY_END_TIMESTEP");

  // Loop over all star particles and set values for those at end of step
  //-----------------------------------------------------------------------------------------------
  for (i=0; i<N; i++) {
    dn = n - star[i]->nlast;
    nstep = star[i]->nstep;

    if (dn == nstep) {
      for (k=0; k<ndim; k++) star[i]->r0[k] = star[i]->r[k];
      for (k=0; k<ndim; k++) star[i]->v0[k] = star[i]->v[k];
      for (k=0; k<ndim; k++) star[i]->a0[k] = star[i]->a[k];
      for (k=0; k<ndim; k++) star[i]->adot0[k] = star[i]->adot[k];
      for (k=0; k<ndim; k++) star[i]->apert[k] = 0.0;
      for (k=0; k<ndim; k++) star[i]->adotpert[k] = 0.0;
      star[i]->active = false;
      star[i]->nlast = n;
      star[i]->tlast = t;
    }
  }
  //-----------------------------------------------------------------------------------------------

  timing->EndTimingSection("NBODY_END_TIMESTEP");

  return;
}



//=================================================================================================
//  NbodyLeapfrogKDK::Timestep
/// Default timestep size for SPH particles.  Takes the minimum of :
/// (i)  const*h/(sound_speed + h*|div_v|)    (Courant condition)
/// (ii) const*sqrt(h/|a|)                    (Acceleration condition)
//=================================================================================================
template <int ndim, template<int> class kernelclass>
DOUBLE NbodyLeapfrogKDK<ndim, kernelclass>::Timestep
 (NbodyParticle<ndim> *star)           ///< Reference to SPH particle
{
  DOUBLE timestep;                     // Minimum value of particle timesteps
  DOUBLE amag;                         // Magnitude of star acceleration
  //DOUBLE adotmag;                      // Magnitude of star jerk

  // Acceleration condition
  amag = sqrt(DotProduct(star->a,star->a,ndim));
  timestep = nbody_mult*sqrt(star->h/(amag + small_number_dp));
  timestep = min(timestep,star->dt_internal);

  // Rate of change of acceleration condition
  //adotmag = sqrt(DotProduct(star->adot,star->adot,ndim));
  //if (amag > 0.0 && adotmag > 0.0)
  //  timestep = min(timestep,nbody_mult*amag/(adotmag + small_number_dp));

  return timestep;
}



//=================================================================================================
//  NbodyLeapfrogKDK::PerturberCorrectionTerms
/// Calculate correction terms on stars due to perturbers.
//=================================================================================================
template <int ndim, template<int> class kernelclass>
void NbodyLeapfrogKDK<ndim, kernelclass>::PerturberCorrectionTerms
 (int n,                               ///< Integer time
  int N,                               ///< No. of stars/systems
  DOUBLE t,                            ///< Current time
  DOUBLE timestep,                     ///< Smallest timestep value
  NbodyParticle<ndim> **star)          ///< Main star/system array
{
  int dn;                              // Integer time since beginning of step
  int i;                               // Particle counter
  int k;                               // Dimension counter
  int nstep;                           // Particle (integer) step size
  DOUBLE dt;                           // Physical time step size
  DOUBLE invdt;                        // 1 / dt

  debug2("[NbodyHermite4::PerturberCorrectionTerms]");

  // Loop over all system particles
  //-----------------------------------------------------------------------------------------------
   for (i=0; i<N; i++) {
     dn = n - star[i]->nlast;
     nstep = star[i]->nstep;

     if (dn == nstep) {
       //dt = timestep*(DOUBLE) nstep;
       dt = t - star[i]->tlast;
       invdt = 1.0 / dt;
       for (k=0; k<ndim; k++) star[i]->a[k] += star[i]->apert[k]*invdt;
     }

   }
   //----------------------------------------------------------------------------------------------

  return;
}




// Template class instances for each dimensionality value (1, 2 and 3) and
// employed kernel (M4, Quintic, Gaussian and tabulated).
template class NbodyLeapfrogKDK<1, M4Kernel>;
template class NbodyLeapfrogKDK<1, QuinticKernel>;
template class NbodyLeapfrogKDK<1, GaussianKernel>;
template class NbodyLeapfrogKDK<1, TabulatedKernel>;
template class NbodyLeapfrogKDK<2, M4Kernel>;
template class NbodyLeapfrogKDK<2, QuinticKernel>;
template class NbodyLeapfrogKDK<2, GaussianKernel>;
template class NbodyLeapfrogKDK<2, TabulatedKernel>;
template class NbodyLeapfrogKDK<3, M4Kernel>;
template class NbodyLeapfrogKDK<3, QuinticKernel>;
template class NbodyLeapfrogKDK<3, GaussianKernel>;
template class NbodyLeapfrogKDK<3, TabulatedKernel>;
