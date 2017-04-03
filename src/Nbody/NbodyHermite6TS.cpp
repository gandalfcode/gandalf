//=================================================================================================
//  NbodyHermite6TS.cpp
//  Contains functions for integrating star particle positions and velocities
//  using the 4th-order Hermite scheme (Makino & Aarseth 1992).
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
//  NbodyHermite6TS::NbodyHermite6TS()
/// N-body 4th-order Hermite class constructor
//=================================================================================================
template <int ndim, template<int> class kernelclass>
NbodyHermite6TS<ndim, kernelclass>::NbodyHermite6TS
(int nbody_softening_aux, int _perturbers, int sub_systems_aux,
 DOUBLE nbody_mult_aux, string KernelName, int Npec) :
  Nbody<ndim>(nbody_softening_aux, _perturbers, sub_systems_aux, nbody_mult_aux, KernelName, Npec),
  kern(kernelclass<ndim>(KernelName))
{
  this->kernp = &kern;
}



//=================================================================================================
//  NbodyHermite6TS::~NbodyHermite6TS()
/// N-body 4th-order Hermite class destructor
//=================================================================================================
template <int ndim, template<int> class kernelclass>
NbodyHermite6TS<ndim, kernelclass>::~NbodyHermite6TS()
{
}



//=================================================================================================
//  NbodyHermite6TS::CalculateDirectGravForces
/// Calculate all star-star force contributions for active systems using
/// direct summation with unsoftened gravity.
//=================================================================================================
template <int ndim, template<int> class kernelclass>
void NbodyHermite6TS<ndim,kernelclass>::CalculateDirectGravForces
 (int N,                               ///< Number of stars
  NbodyParticle<ndim> **star,          ///< Array of stars/systems
  DomainBox<ndim> &simbox,             ///< [in] Simulation domain box
  Ewald<ndim> *ewald)                  ///< [in] Ewald gravity object pointer
{
  FLOAT a[ndim];                       // Acceleration
  FLOAT adot[ndim];                    // 1st time derivative of accel (jerk)
  FLOAT afac,bfac;                     // Aux. summation variables
  FLOAT aperiodic[ndim];               // Ewald periodic grav. accel correction
  FLOAT da[ndim];                      // Relative acceleration
  FLOAT dr[ndim];                      // Relative position vector
  FLOAT dr_corr[ndim];                 // Periodic corrected position vector
  FLOAT drdt;                          // Rate of change of distance
  FLOAT drsqd;                         // Distance squared
  FLOAT dv[ndim];                      // Relative velocity vector
  FLOAT dvsqd;                         // Velocity squared
  FLOAT invdrmag;                      // 1 / drmag
  FLOAT invdrsqd;                      // 1 / drsqd
  FLOAT potperiodic;                   // Periodic correction for grav. potential

  debug2("[NbodyHermite6TS::CalculateDirectGravForces]");

  // Loop over all (active) stars
  //-----------------------------------------------------------------------------------------------
#pragma omp parallel for if (N > this->maxNbodyOpenMp) default(none) shared(ewald, N, simbox, star) \
private(a, adot, afac, aperiodic, bfac, da, dr, dr_corr, drdt, drsqd, dv, dvsqd, invdrmag, invdrsqd, potperiodic)
  for (int i=0; i<N; i++) {
    if (not star[i]->flags.check(active)) continue;

    // Sum grav. contributions for all other stars (excluding star itself)
    //---------------------------------------------------------------------------------------------
    for (int j=0; j<N; j++) {
      if (i == j) continue;

      for (int k=0; k<ndim; k++) dr[k] = star[j]->r[k] - star[i]->r[k];
      for (int k=0; k<ndim; k++) dv[k] = star[j]->v[k] - star[i]->v[k];
      NearestPeriodicVector(simbox, dr, dr_corr);
      drsqd = DotProduct(dr,dr,ndim) + small_number_dp;
      invdrmag = 1.0/sqrt(drsqd);
      drdt = DotProduct(dv,dr,ndim)*invdrmag;

      star[i]->gpot += star[j]->m*invdrmag;
      for (int k=0; k<ndim; k++) star[i]->a[k] += star[j]->m*dr[k]*pow(invdrmag,3);
      for (int k=0; k<ndim; k++) star[i]->adot[k] +=
        star[j]->m*pow(invdrmag,3)*(dv[k] - 3.0*drdt*invdrmag*dr[k]);

      // Add periodic gravity contribution (if activated)
      if (simbox.PeriodicGravity) {
        ewald->CalculatePeriodicCorrection(star[j]->m, dr, aperiodic, potperiodic);
        for (int k=0; k<ndim; k++) star[i]->a[k] += aperiodic[k];
        star[i]->gpot += potperiodic;
      }

    }
    //---------------------------------------------------------------------------------------------

  }
  //-----------------------------------------------------------------------------------------------


  // Loop over all stars a second time to compute 2nd time derivative
  //-----------------------------------------------------------------------------------------------
#pragma omp parallel for if (N > this->maxNbodyOpenMp) default(none) shared(ewald, N, simbox, star) \
private(a, adot, afac, aperiodic, bfac, da, dr, dr_corr, drdt, drsqd, dv, dvsqd, invdrmag, invdrsqd, potperiodic)
  for (int i=0; i<N; i++) {
    if (not star[i]->flags.check(active)) continue;
    for (int k=0; k<ndim; k++) star[i]->a2dot[k] = 0.0;

    // Sum grav. contributions for all other stars (excluding star itself)
    //---------------------------------------------------------------------------------------------
    for (int j=0; j<N; j++) {
      if (i == j) continue;

      for (int k=0; k<ndim; k++) dr[k] = star[j]->r[k] - star[i]->r[k];
      for (int k=0; k<ndim; k++) dv[k] = star[j]->v[k] - star[i]->v[k];
      for (int k=0; k<ndim; k++) da[k] = star[j]->a[k] - star[i]->a[k];
      drsqd = DotProduct(dr,dr,ndim) + small_number_dp;
      dvsqd = DotProduct(dv,dv,ndim);
      invdrsqd = 1.0/drsqd;
      invdrmag = sqrt(invdrsqd);
      drdt = DotProduct(dv,dr,ndim)*invdrmag;
      for (int k=0; k<ndim; k++) a[k] = star[j]->m*dr[k]*pow(invdrmag,3);
      for (int k=0; k<ndim; k++) adot[k] =
        star[j]->m*pow(invdrmag,3)*(dv[k] - 3.0*drdt*invdrmag*dr[k]);

      // Now compute 2nd and 3rd order derivatives
      afac = DotProduct(dv,dr,ndim)*invdrsqd;
      bfac = dvsqd*invdrsqd + afac*afac + DotProduct(da,dr,ndim)*invdrsqd;

      for (int k=0; k<ndim; k++) star[i]->a2dot[k] =
        star[j]->m*da[k]*invdrsqd*invdrmag - 6.0*afac*adot[k] - 3.0*bfac*a[k];

    }
    //---------------------------------------------------------------------------------------------

  }
  //-----------------------------------------------------------------------------------------------

  return;
}



//=================================================================================================
//  NbodyHermite6TS::CalculateDirectSmoothedGravForces
/// Calculate all star-star force contributions for active systems using kernel-softened
/// softened gravity.
//=================================================================================================
template <int ndim, template<int> class kernelclass>
void NbodyHermite6TS<ndim, kernelclass>::CalculateDirectSmoothedGravForces
 (int N,                               ///< [in] Number of stars
  NbodyParticle<ndim> **star,          ///< [inout] Array of stars/systems
  DomainBox<ndim> &simbox,             ///< [in] Simulation domain box
  Ewald<ndim> *ewald)                  ///< [in] Ewald gravity object pointer
{
  FLOAT aperiodic[ndim];               // Ewald periodic grav. accel correction
  FLOAT dr[ndim];                      // Relative position vector
  FLOAT dr_corr[ndim];                 // Periodic corrected position vector
  FLOAT drdt;                          // Rate of change of distance
  FLOAT drmag;                         // Distance
  FLOAT drsqd;                         // Distance squared
  FLOAT dv[ndim];                      // Relative velocity vector
  FLOAT invdrmag;                      // 1 / drmag
  FLOAT invhmean;                      // 1 / mean of star smoothing lengths
  FLOAT paux;                          // Common force factor
  FLOAT potperiodic;                   // Periodic correction for grav. potential
  FLOAT wmean;                         // Mean-h kernel factor

  debug2("[NbodyHermite6TS::CalculateDirectSmoothedGravForces]");

  // Loop over all (active) stars
  //-----------------------------------------------------------------------------------------------
#pragma omp parallel for if (N > this->maxNbodyOpenMp) default(none) shared(ewald, N, simbox, star) \
private(aperiodic, dr, dr_corr, drdt, drmag, drsqd, dv, invdrmag, invhmean, paux, potperiodic, wmean)
  for (int i=0; i<N; i++) {
    if (not star[i]->flags.check(active)) continue;

    // Sum grav. contributions for all other stars (excluding star itself)
    //---------------------------------------------------------------------------------------------
    for (int j=0; j<N; j++) {
      if (i == j) continue;

      for (int k=0; k<ndim; k++) dr[k] = star[j]->r[k] - star[i]->r[k];
      for (int k=0; k<ndim; k++) dv[k] = star[j]->v[k] - star[i]->v[k];
      NearestPeriodicVector(simbox, dr, dr_corr);
      drsqd = DotProduct(dr,dr,ndim);
      drmag = sqrt(drsqd);
      invdrmag = 1.0/sqrt(drsqd);
      invhmean = 2.0/(star[i]->h + star[j]->h);
      drdt = DotProduct(dv,dr,ndim)*invdrmag;
      paux = star[j]->m*invhmean*invhmean*kern.wgrav(drmag*invhmean)*invdrmag;
      wmean = kern.w0(drmag*invhmean)*powf(invhmean,ndim);

      // Add contribution to main star array
      star[i]->gpot += star[j]->m*invhmean*kern.wpot(drmag*invhmean);
      for (int k=0; k<ndim; k++) star[i]->a[k] += paux*dr[k];
      for (int k=0; k<ndim; k++) star[i]->adot[k] += paux*dv[k] -
        3.0*paux*drdt*invdrmag*dr[k] + 2.0*twopi*star[j]->m*drdt*wmean*invdrmag*dr[k];

      // Add periodic gravity contribution (if activated)
      if (simbox.PeriodicGravity) {
        ewald->CalculatePeriodicCorrection(star[j]->m, dr, aperiodic, potperiodic);
        for (int k=0; k<ndim; k++) star[i]->a[k] += aperiodic[k];
        star[i]->gpot += potperiodic;
      }

    }
    //---------------------------------------------------------------------------------------------

  }
  //-----------------------------------------------------------------------------------------------

  return;
}



//=================================================================================================
//  NbodyHermite6TS::CalculateDirectHydroForces
/// Calculate all ..
//=================================================================================================
template <int ndim, template<int> class kernelclass>
void NbodyHermite6TS<ndim, kernelclass>::CalculateDirectHydroForces
 (NbodyParticle<ndim> *star,           ///< [inout] Pointer to star
  int Nhydro,                          ///< [in] Number of gas particles
  int Ndirect,                         ///< [in] ..
  int *hydrolist,                      ///< [in] ..
  int *directlist,                     ///< [in] ..
  Hydrodynamics<ndim> *hydro,          ///< [in] Hydrodynamics object
  DomainBox<ndim> &simbox,             ///< [in] Simulation domain box
  Ewald<ndim> *ewald)                  ///< [in] Ewald gravity object pointer
{
  int j,jj,k;                          // Star and dimension counters
  FLOAT aperiodic[ndim];               // Ewald periodic grav. accel correction
  FLOAT dr[ndim];                      // Relative position vector
  FLOAT dr_corr[ndim];                 // Periodic corrected position vector
  FLOAT drmag;                         // Distance
  FLOAT drsqd;                         // Distance squared
  FLOAT drdt;                          // Rate of change of distance
  FLOAT dv[ndim];                      // Relative velocity vector
  FLOAT invhmean;                      // 1 / hmean
  FLOAT invdrmag;                      // 1 / drmag
  FLOAT paux;                          // Aux. force variable
  FLOAT potperiodic;                   // Periodic correction for grav. potential
  FLOAT wkern;                         // SPH kernel value

  debug2("[NbodyHermite6TS::CalculateDirectHydroForces]");


  // Sum grav. contributions from all neighbouring SPH particles
  //-----------------------------------------------------------------------------------------------
  for (jj=0; jj<Nhydro; jj++) {

    j = hydrolist[jj];
    Particle<ndim>& part = hydro->GetParticlePointer(j);

    for (k=0; k<ndim; k++) dr[k] = part.r[k] - star->r[k];
    for (k=0; k<ndim; k++) dv[k] = part.v[k] - star->v[k];
    NearestPeriodicVector(simbox, dr, dr_corr);
    drsqd    = DotProduct(dr,dr,ndim);
    drmag    = sqrt(drsqd);
    invdrmag = 1.0/drmag;
    invhmean = 2.0/(star->h + part.h);
    drdt     = DotProduct(dv,dr,ndim)*invdrmag;
    paux     = part.m*invhmean*invhmean*kern.wgrav(drmag*invhmean)*invdrmag;
    wkern    = kern.w0(drmag*invhmean)*powf(invhmean,ndim);

    // Add contribution to main star array
    for (k=0; k<ndim; k++) star->a[k] += paux*dr[k];
    for (k=0; k<ndim; k++) star->adot[k] += paux*dv[k] -
      3.0*paux*drdt*invdrmag*dr[k] + 2.0*twopi*part.m*drdt*wkern*invdrmag*dr[k];
    star->gpot += part.m*invhmean*kern.wpot(drmag*invhmean);

    // Add periodic gravity contribution (if activated)
    if (simbox.PeriodicGravity) {
      ewald->CalculatePeriodicCorrection(part.m, dr, aperiodic, potperiodic);
      for (k=0; k<ndim; k++) star->a[k] += aperiodic[k];
      star->gpot += potperiodic;
    }

  }
  //-----------------------------------------------------------------------------------------------


  // Now include contributions from distant, non-SPH neighbours
  // (i.e. direct summation with Newton's law of gravity)
  //-----------------------------------------------------------------------------------------------
  for (jj=0; jj<Ndirect; jj++) {

    j = directlist[jj];
    Particle<ndim>& part = hydro->GetParticlePointer(j);

    for (k=0; k<ndim; k++) dr[k] = part.r[k] - star->r[k];
    for (k=0; k<ndim; k++) dv[k] = part.v[k] - star->v[k];
    NearestPeriodicVector(simbox, dr, dr_corr);
    drsqd = DotProduct(dr,dr,ndim);
    drmag = sqrt(drsqd);
    invdrmag = 1.0/drmag;
    drdt = DotProduct(dv,dr,ndim)*invdrmag;

    // Add contribution to main star array
    for (k=0; k<ndim; k++) star->a[k] += part.m*dr[k]*pow(invdrmag,3);
    for (k=0; k<ndim; k++) star->adot[k] +=
      part.m*pow(invdrmag,3)*(dv[k] - 3.0*drdt*invdrmag*dr[k]);
    star->gpot += part.m*invdrmag;

    // Add periodic gravity contribution (if activated)
    if (simbox.PeriodicGravity) {
      ewald->CalculatePeriodicCorrection(part.m, dr, aperiodic, potperiodic);
      for (k=0; k<ndim; k++) star->a[k] += aperiodic[k];
      star->gpot += potperiodic;
    }

  }
  //-----------------------------------------------------------------------------------------------


  return;
}



//=================================================================================================
//  NbodyHermite6TS::CalculateAllStartupQuantities
/// Calculate all initial properties that are required before the first
/// timestep integration, i.e. 2nd and 3rd acceleration time derivatives.
//=================================================================================================
template <int ndim, template<int> class kernelclass>
void NbodyHermite6TS<ndim, kernelclass>::CalculateAllStartupQuantities
 (int N,                               ///< Number of stars
  NbodyParticle<ndim> **star,          ///< Array of stars/systems
  DomainBox<ndim> &simbox,             ///< [in] Simulation domain box
  Ewald<ndim> *ewald)                  ///< [in] Ewald gravity object pointer
{
  int i,j,k;                           // Star and dimension counters
  FLOAT a[ndim];                       // Acceleration
  FLOAT adot[ndim];                    // 1st time derivative of accel (jerk)
  FLOAT a2dot[ndim];                   // 2nd time deriivative of acceleration
  FLOAT afac,bfac,cfac;                // Aux. summation variables
  FLOAT da[ndim];                      // Relative acceleration
  FLOAT dadot[ndim];                   // Relative jerk
  FLOAT dr[ndim];                      // Relative position vector
  FLOAT dr_corr[ndim];                 // Periodic corrected position vector
  FLOAT drdt;                          // Rate of change of distance
  FLOAT drsqd;                         // Distance squared
  FLOAT dv[ndim];                      // Relative velocity vector
  FLOAT dvsqd;                         // Velocity squared
  FLOAT invdrmag;                      // 1 / drmag
  FLOAT invdrsqd;                      // 1 / drsqd

  debug2("[NbodyHermite4::CalculateAllStartupQuantities]");

  // Loop over all stars
  //-----------------------------------------------------------------------------------------------
  for (i=0; i<N; i++) {

    for (k=0; k<ndim; k++) star[i]->a2dot[k] = (FLOAT) 0.0;
    for (k=0; k<ndim; k++) star[i]->a3dot[k] = (FLOAT) 0.0;

    // Sum grav. contributions for all other stars (excluding star itself)
    //---------------------------------------------------------------------------------------------
    for (j=0; j<N; j++) {
      if (i == j) continue;

      for (k=0; k<ndim; k++) dr[k] = star[j]->r[k] - star[i]->r[k];
      for (k=0; k<ndim; k++) dv[k] = star[j]->v[k] - star[i]->v[k];
      for (k=0; k<ndim; k++) da[k] = star[j]->a[k] - star[i]->a[k];
      for (k=0; k<ndim; k++) dadot[k] = star[j]->adot[k] - star[i]->adot[k];
      NearestPeriodicVector(simbox, dr, dr_corr);
      drsqd = DotProduct(dr,dr,ndim) + small_number_dp;
      dvsqd = DotProduct(dv,dv,ndim);
      invdrsqd = 1.0/drsqd;
      invdrmag = sqrt(invdrsqd);
      drdt = DotProduct(dv,dr,ndim)*invdrmag;
      for (k=0; k<ndim; k++) a[k] = star[j]->m*dr[k]*pow(invdrmag,3);
      for (k=0; k<ndim; k++) adot[k] =
        star[j]->m*pow(invdrmag,3)*(dv[k] - 3.0*drdt*invdrmag*dr[k]);

      // Now compute 2nd and 3rd order derivatives
      afac = DotProduct(dv,dr,ndim)*invdrsqd;
      bfac = dvsqd*invdrsqd + afac*afac + DotProduct(da,dr,ndim)*invdrsqd;
      cfac = 3.0*DotProduct(dv,da,ndim)*invdrsqd +
        DotProduct(dr,dadot,ndim)*invdrsqd + afac*(3.0*bfac - 4.0*afac*afac);

      for (k=0; k<ndim; k++) a2dot[k] =
        star[j]->m*da[k]*invdrsqd*invdrmag - 6.0*afac*adot[k] - 3.0*bfac*a[k];
      for (k=0; k<ndim; k++) star[i]->a2dot[k] += a2dot[k];
      for (k=0; k<ndim; k++) star[i]->a3dot[k] += star[j]->m*dadot[k]*invdrsqd*invdrmag -
        9.0*afac*a2dot[k] - 9.0*bfac*adot[k] - 3.0*cfac*a[k];

    }
    //---------------------------------------------------------------------------------------------

  }
  //-----------------------------------------------------------------------------------------------


  return;
}



//=================================================================================================
//  NbodyHermite6TS::AdvanceParticles
/// Integrate particle positions to 3nd order, and particle velocities to 2nd
/// order from the beginning of the step to the current simulation time, i.e.
/// r(t+dt) = r(t) + v(t)*dt + a(t)*dt^2/2 + adot(t)*dt^3/6,
/// v(t+dt) = v(t) + a(t)*dt + adot(t)*dt^2/2.
/// Also set particles at the end of step as 'active' in order to compute
/// the end-of-step force computation.
//=================================================================================================
template <int ndim, template<int> class kernelclass>
void NbodyHermite6TS<ndim, kernelclass>::AdvanceParticles
(int n,                             ///< Integer time
 int N,                             ///< No. of stars/systems
 FLOAT t,                          ///< Current time
 FLOAT timestep,                   ///< Smallest timestep value
 NbodyParticle<ndim> **star)        ///< Main star/system array
{
  int dn;                           // Integer time since beginning of step
  int i;                            // Particle counter
  int k;                            // Dimension counter
  int nstep;                        // Particle (integer) step size
  FLOAT dt;                        // Timestep since start of step
  FLOAT dt2;                       // dt*dt

  debug2("[NbodyHermite6TS::AdvanceParticles]");

  // Advance positions and velocities of all systems
  //-----------------------------------------------------------------------------------------------
  for (i=0; i<N; i++) {

    // Compute time since beginning of step
    nstep = star[i]->nstep;
    dn    = n - star[i]->nlast;
    //dt = timestep*(FLOAT) dn;
    dt    = t - star[i]->tlast;
    dt2   = dt*dt;

    // Advance positions to third order and velocities to second order
    for (k=0; k<ndim; k++) star[i]->r[k] = star[i]->r0[k] + star[i]->v0[k]*dt +
      (FLOAT) 0.5*star[i]->a0[k]*dt2 + onesixth_dp*star[i]->adot0[k]*dt2*dt +
      (FLOAT) 0.5*onetwelfth_dp*star[i]->a2dot0[k]*dt2*dt2;
    for (k=0; k<vdim; k++) star[i]->v[k] = star[i]->v0[k] + star[i]->a0[k]*dt +
      (FLOAT) 0.5*star[i]->adot0[k]*dt2 + onesixth_dp*star[i]->a2dot0[k]*dt2*dt;

    // If at end of step, set system particle as active
    if (dn == nstep) star[i]->flags.set(active);
    else star[i]->flags.unset(active);
  }
  //-----------------------------------------------------------------------------------------------

  return;
}



//=================================================================================================
//  NbodyHermite6TS::CorrectionTerms
/// Compute 2nd and 3rd time derivatives of the acceleration using Hermite interpolation.  Finally
/// correct positions to 5th order and velocities to 4th order using higher-order derivatives.
//=================================================================================================
template <int ndim, template<int> class kernelclass>
void NbodyHermite6TS<ndim, kernelclass>::CorrectionTerms
(int n,                             ///< Integer time
 int N,                             ///< No. of stars/systems
 FLOAT t,                          ///< Current time
 FLOAT timestep,                   ///< Smallest timestep value
 NbodyParticle<ndim> **star)        ///< Main star/system array
{
  int dn;                           // Integer time since beginning of step
  int i;                            // Particle counter
  int k;                            // Dimension counter
  int nstep;                        // Particle (integer) step size
  FLOAT dt;                        // Physical time step size
  FLOAT dt3;                       // dt*dt*dt
  FLOAT invdt;                     // 1 / dt
  static const FLOAT one120 = 1.0/120.0;  // 1/120

  debug2("[NbodyHermite6TS::CorrectionTerms]");

  // Loop over all system particles
  //-----------------------------------------------------------------------------------------------
  for (i=0; i<N; i++) {
    dn = n - star[i]->nlast;
    nstep = star[i]->nstep;

    if (dn == nstep) {
      //dt = timestep*(FLOAT) nstep;
      dt = t - star[i]->tlast;
      dt3 = powf(dt,3);
      invdt = 1.0 / dt;

      for (k=0; k<ndim; k++) {
        star[i]->a3dot[k] = (12.0*(star[i]->a0[k] - star[i]->a[k]) +
                             6.0*dt*(star[i]->adot0[k] + star[i]->adot[k]))*invdt*invdt*invdt;
      }

      for (k=0; k<ndim; k++) {
        star[i]->v[k] = star[i]->v0[k] + 0.5*(star[i]->a0[k] + star[i]->a[k])*dt -
          0.1*(star[i]->adot[k] - star[i]->adot0[k])*dt*dt +
          one120*(star[i]->a2dot[k] + star[i]->a2dot0[k])*dt3;
        star[i]->r[k] = star[i]->r0[k] + 0.5*(star[i]->v0[k] + star[i]->v[k])*dt -
          0.1*(star[i]->a[k] - star[i]->a0[k])*dt*dt +
          one120*(star[i]->adot[k] + star[i]->adot0[k])*dt3;
      }

    }

  }
  //-----------------------------------------------------------------------------------------------

  return;
}



//=================================================================================================
//  NbodyHermite6TS::PerturberCorrectionTerms
/// ..
//=================================================================================================
template <int ndim, template<int> class kernelclass>
void NbodyHermite6TS<ndim, kernelclass>::PerturberCorrectionTerms
(int n,                             ///< Integer time
 int N,                             ///< No. of stars/systems
 FLOAT t,                          ///< Current time
 FLOAT timestep,                   ///< Smallest timestep value
 NbodyParticle<ndim> **star)        ///< Main star/system array
{
  return;
}



//=================================================================================================
//  NbodyHermite6TS::EndTimestep
/// Record all important star particle quantities at the end of the step
/// for the start of the new timestep.
//=================================================================================================
template <int ndim, template<int> class kernelclass>
void NbodyHermite6TS<ndim, kernelclass>::EndTimestep
(int n,                             ///< Integer time
 int N,                             ///< No. of stars/systems
 FLOAT t,                          ///< Current time
 FLOAT timestep,                   ///< Smallest timestep value
 NbodyParticle<ndim> **star)        ///< Main star/system array
{
  int i;                            // Particle counter
  int k;                            // Dimension counter

  debug2("[NbodyHermite6TS::EndTimestep]");

  // Loop over all system particles
  //-----------------------------------------------------------------------------------------------
  for (i=0; i<N; i++) {

    // If at end of the current step, set quantites for start of new step
    if (star[i]->flags.check(end_timestep)) {
      for (k=0; k<ndim; k++) star[i]->r0[k] = star[i]->r[k];
      for (k=0; k<ndim; k++) star[i]->v0[k] = star[i]->v[k];
      for (k=0; k<ndim; k++) star[i]->a0[k] = star[i]->a[k];
      for (k=0; k<ndim; k++) star[i]->adot0[k] = star[i]->adot[k];
      for (k=0; k<ndim; k++) star[i]->a2dot0[k] = star[i]->a2dot[k];
      for (k=0; k<ndim; k++) star[i]->apert[k] = 0.0;
      for (k=0; k<ndim; k++) star[i]->adotpert[k] = 0.0;
      star[i]->nlast = n;
      star[i]->tlast = t;
      star[i]->dt = star[i]->dt_next ;
      star[i]->dt_next = 0 ;
      star[i]->flags.unset(active);
      star[i]->flags.unset(end_timestep);
    }

  }
  //-----------------------------------------------------------------------------------------------

  return;
}



//=================================================================================================
//  NbodyHermite6TS::Timestep
/// Calculate the N-body timestep for a given star using the standard Aarseth timestep, i.e.
/// $dt = gamma*\sqrt{\frac{a*a2 + a1^2}{a1*a3 + a2^2}}$.
//=================================================================================================
template <int ndim, template<int> class kernelclass>
DOUBLE NbodyHermite6TS<ndim, kernelclass>::Timestep
(NbodyParticle<ndim> *star)         ///< Reference to star/system particle
{
  DOUBLE timestep;                  // Minimum value of particle timesteps
  DOUBLE asqd;                      // Magnitude of particle acceleration
  DOUBLE a1sqd;                     // Magnitude of particle acceleration derivative
  DOUBLE a2sqd;                     // Magnitude of particle acceleration 2nd derivative
  DOUBLE a3sqd;                     // Magnitude of particle acceleration 3rd derivative

  asqd  = DotProduct(star->a,star->a,ndim);
  a1sqd = DotProduct(star->adot,star->adot,ndim);
  a2sqd = DotProduct(star->a2dot,star->a2dot,ndim);
  a3sqd = DotProduct(star->a3dot,star->a3dot,ndim);

  // Normal case of all four accel quantities being defined
  if (a1sqd > small_number_dp && a2sqd > small_number_dp) {
    timestep = (sqrt(asqd*a2sqd) + a1sqd)/(sqrt(a1sqd*a3sqd) + a2sqd);
    timestep = nbody_mult*sqrt(timestep);
  }
  // Special case when 1st derivative falls to zero
  else if (asqd > small_number_dp && a2sqd > small_number_dp) {
    timestep = asqd/(a2sqd + small_number_dp);
    timestep = nbody_mult*sqrt(timestep);
  }
  // If all else fails, use simple criterion
  else if (asqd > small_number_dp) {
    timestep = sqrt(star->h/(sqrt(asqd) + small_number_dp));
  }
  else {
    timestep = big_number_dp;
  }

  timestep = min(timestep,star->dt_internal);

  return timestep;
}



// Template class instances for each dimensionality value (1, 2 and 3) and
// employed kernel (M4, Quintic, Gaussian and tabulated).
template class NbodyHermite6TS<1, M4Kernel>;
template class NbodyHermite6TS<2, M4Kernel>;
template class NbodyHermite6TS<3, M4Kernel>;
template class NbodyHermite6TS<1, QuinticKernel>;
template class NbodyHermite6TS<2, QuinticKernel>;
template class NbodyHermite6TS<3, QuinticKernel>;
template class NbodyHermite6TS<1, GaussianKernel>;
template class NbodyHermite6TS<2, GaussianKernel>;
template class NbodyHermite6TS<3, GaussianKernel>;
template class NbodyHermite6TS<1, TabulatedKernel>;
template class NbodyHermite6TS<2, TabulatedKernel>;
template class NbodyHermite6TS<3, TabulatedKernel>;
