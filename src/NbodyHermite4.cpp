//=============================================================================
//  NbodyHermite4.cpp
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
//=============================================================================


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
#include "SphKernel.h"
#include "Debug.h"
#include "Exception.h"
#include "InlineFuncs.h"
using namespace std;



//=============================================================================
//  NbodyHermite4::NbodyHermite4()
/// N-body 4th-order Hermite class constructor
//=============================================================================
template <int ndim, template<int> class kernelclass>
NbodyHermite4<ndim, kernelclass>::NbodyHermite4
(int nbody_softening_aux, int sub_systems_aux, 
 DOUBLE nbody_mult_aux, string KernelName, int Npec) :
  Nbody<ndim>(nbody_softening_aux, sub_systems_aux, 
              nbody_mult_aux, KernelName, Npec),
  kern(kernelclass<ndim>(KernelName))
{
}


//=============================================================================
//  NbodyHermite4::~NbodyHermite4()
/// N-body 4th-order Hermite class destructor
//=============================================================================
template <int ndim, template<int> class kernelclass>
NbodyHermite4<ndim, kernelclass>::~NbodyHermite4()
{
}



//=============================================================================
//  NbodyHermite4::CalculateDirectGravForces
/// Calculate all star-star force contributions for active systems using 
/// direct summation with unsoftened gravity.
//=============================================================================
template <int ndim, template<int> class kernelclass>
void NbodyHermite4<ndim, kernelclass>::CalculateDirectGravForces
(int N,                             ///< Number of stars
 NbodyParticle<ndim> **star)        ///< Array of stars/systems
{
  int i,j,k;                        // Star and dimension counters
  DOUBLE dr[ndim];                  // Relative position vector
  DOUBLE drdt;                      // Rate of change of distance
  DOUBLE drsqd;                     // Distance squared
  DOUBLE dv[ndim];                  // Relative velocity vector
  DOUBLE invdrmag;                  // 1 / drmag

  debug2("[NbodyHermite4::CalculateDirectGravForces]");

  // Loop over all (active) stars
  //---------------------------------------------------------------------------
  for (i=0; i<N; i++) {
    if (star[i]->active == 0) continue;

    // Sum grav. contributions for all other stars (excluding star itself)
    //-------------------------------------------------------------------------
    for (j=0; j<N; j++) {
      if (i == j) continue;

      for (k=0; k<ndim; k++) dr[k] = star[j]->r[k] - star[i]->r[k];
      for (k=0; k<ndim; k++) dv[k] = star[j]->v[k] - star[i]->v[k];
      drsqd = DotProduct(dr,dr,ndim);
      invdrmag = 1.0/sqrt(drsqd);
      drdt = DotProduct(dv,dr,ndim)*invdrmag;

      star[i]->gpot += star[j]->m*invdrmag;
      for (k=0; k<ndim; k++) star[i]->a[k] += star[j]->m*dr[k]*pow(invdrmag,3);
      for (k=0; k<ndim; k++) star[i]->adot[k] +=
        star[j]->m*pow(invdrmag,3)*(dv[k] - 3.0*drdt*invdrmag*dr[k]);

    }
    //-------------------------------------------------------------------------

  }
  //---------------------------------------------------------------------------

  return;
}



//=============================================================================
//  NbodyHermite4::CalculatePerturberForces
/// Calculate perturber forces on all stars in a N-body sub-system.
//=============================================================================
template <int ndim, template<int> class kernelclass>
void NbodyHermite4<ndim, kernelclass>::CalculatePerturberForces
(int N,                             ///< Number of stars
 int Npert,                         ///< Number of perturbing stars
 NbodyParticle<ndim> **star,        ///< Array of stars/systems
 NbodyParticle<ndim> *perturber,    ///< Array of perturbing stars/systems
 DOUBLE *apert,                     ///< Acceleration due to perturbers
 DOUBLE *adotpert)                  ///< Jerk due to perturbers
{
  int i,j,k;                        // Star and dimension counters
  DOUBLE dr[ndim];                  // Relative position vector
  DOUBLE drdt;                      // Rate of change of distance
  DOUBLE drsqd;                     // Distance squared
  DOUBLE dv[ndim];                  // Relative velocity vector
  DOUBLE invdrmag;                  // 1 / drmag
  DOUBLE rcom[ndim];                // Position of centre-of-mass
  DOUBLE vcom[ndim];                // Velocity of centre-of-mass
  DOUBLE msystot = 0.0;             // Total system mass
  //DOUBLE apertcom[ndim*Npertmax];   // ..
  //DOUBLE apertdotcom[ndim*Npertmax]; // ..

  debug2("[NbodyHermite4::CalculatePerturberForces]");

  // First, compute position and velocity of system COM
  for (k=0; k<ndim; k++) rcom[k] = 0.0;
  for (k=0; k<ndim; k++) vcom[k] = 0.0;
  for (i=0; i<N; i++) {
    msystot += star[i]->m;
    for (k=0; k<ndim; k++) rcom[k] += star[i]->m*star[i]->r[k];
    for (k=0; k<ndim; k++) vcom[k] += star[i]->m*star[i]->v[k];
  }
  for (k=0; k<ndim; k++) rcom[k] /= msystot;
  for (k=0; k<ndim; k++) vcom[k] /= msystot;

  for (j=0; j<Npert; j++) {
    for (k=0; k<ndim; k++) dr[k] = rcom[k] - perturber[j].r[k];
    for (k=0; k<ndim; k++) dv[k] = vcom[k] - perturber[j].v[k];
    drsqd = DotProduct(dr,dr,ndim);
    invdrmag = 1.0/sqrt(drsqd);
    drdt = DotProduct(dv,dr,ndim)*invdrmag; 
    for (k=0; k<ndim; k++) apert[ndim*j + k] = -msystot*dr[k]*pow(invdrmag,3);
    for (k=0; k<ndim; k++) adotpert[ndim*j + k] =
      -msystot*pow(invdrmag,3)*(dv[k] - 3.0*drdt*invdrmag*dr[k]);
  }

  // Loop over all (active) stars
  //---------------------------------------------------------------------------
  for (i=0; i<N; i++) {
    if (star[i]->active == 0) continue;

    // Sum grav. contributions for all perturbing stars.  
    //-------------------------------------------------------------------------
    for (j=0; j<Npert; j++) {

      for (k=0; k<ndim; k++) dr[k] = perturber[j].r[k] - star[i]->r[k];
      for (k=0; k<ndim; k++) dv[k] = perturber[j].v[k] - star[i]->v[k];
      drsqd = DotProduct(dr,dr,ndim);
      invdrmag = 1.0/sqrt(drsqd);
      drdt = DotProduct(dv,dr,ndim)*invdrmag;

      // First, add contribution of perturber to star
      star[i]->gpe_pert += perturber[j].m*invdrmag;
      for (k=0; k<ndim; k++) star[i]->a[k] += 
        perturber[j].m*dr[k]*pow(invdrmag,3);
      for (k=0; k<ndim; k++) star[i]->adot[k] += perturber[j].m*
        pow(invdrmag,3)*(dv[k] - 3.0*drdt*invdrmag*dr[k]);

      // Next, add contribution of star to perturber
      for (k=0; k<ndim; k++) apert[ndim*j + k] -=
        star[i]->m*dr[k]*pow(invdrmag,3);
      for (k=0; k<ndim; k++) adotpert[ndim*j + k] -=
        star[i]->m*pow(invdrmag,3)*(dv[k] - 3.0*drdt*invdrmag*dr[k]);

    }
    //-------------------------------------------------------------------------

  }
  //---------------------------------------------------------------------------

  return;
}



//=============================================================================
//  NbodyHermite4::CalculateDirectSPHForces
/// Calculate all ..
//=============================================================================
template <int ndim, template<int> class kernelclass>
void NbodyHermite4<ndim, kernelclass>::CalculateDirectSPHForces
(NbodyParticle<ndim> *star,         ///< [inout] Pointer to star
 int Nsph,                          ///< [in] Number of gas particles
 int Ndirect,                       ///< [in] ..
 int *sphlist,                      ///< [in] ..
 int *directlist,                   ///< [in] ..
 SphParticle<ndim> *sphdata)        ///< [in] Array of SPH particles
{
  int j,jj,k;                        // Star and dimension counters
  DOUBLE dr[ndim];                  // Relative position vector
  DOUBLE drmag;                     // Distance
  DOUBLE drsqd;                     // Distance squared
  DOUBLE drdt;                      // Rate of change of distance
  DOUBLE dv[ndim];                  // Relative velocity vector
  DOUBLE invhmean;                  // 1 / hmean
  DOUBLE invdrmag;                  // 1 / drmag
  DOUBLE paux;                      // Aux. force variable
  DOUBLE wkern;                     // SPH kernel value

  debug2("[NbodyHermite4::CalculateDirectSPHForces]");


  // Sum grav. contributions from all neighbouring SPH particles
  //---------------------------------------------------------------------------
  for (jj=0; jj<Nsph; jj++) {

    j = sphlist[jj];
    for (k=0; k<ndim; k++) dr[k] = sphdata[j].r[k] - star->r[k];
    for (k=0; k<ndim; k++) dv[k] = sphdata[j].v[k] - star->v[k];
    drsqd = DotProduct(dr,dr,ndim);
    drmag = sqrt(drsqd);
    invdrmag = 1.0/drmag;
    invhmean = 2.0/(star->h + sphdata[j].h);
    drdt = DotProduct(dv,dr,ndim)*invdrmag;
    
    paux = sphdata[j].m*invhmean*invhmean*
      kern.wgrav(drmag*invhmean)*invdrmag;
    wkern = kern.w0(drmag*invhmean)*powf(invhmean,ndim);
    
    // Add contribution to main star array
    for (k=0; k<ndim; k++) star->a[k] += paux*dr[k];
    for (k=0; k<ndim; k++) star->adot[k] += paux*dv[k] - 
      3.0*paux*drdt*invdrmag*dr[k] + 
      2.0*twopi*sphdata[j].m*drdt*wkern*invdrmag*dr[k];
    star->gpot += sphdata[j].m*invhmean*kern.wpot(drmag*invhmean);

  }
  //---------------------------------------------------------------------------


  // Now include contributions from distant, non-SPH neighbours
  // (i.e. direct summation with Newton's law of gravity)
  //---------------------------------------------------------------------------
  for (jj=0; jj<Ndirect; jj++) {

    j = directlist[jj];
    for (k=0; k<ndim; k++) dr[k] = sphdata[j].r[k] - star->r[k];
    for (k=0; k<ndim; k++) dv[k] = sphdata[j].v[k] - star->v[k];
    drsqd = DotProduct(dr,dr,ndim);
    drmag = sqrt(drsqd);
    invdrmag = 1.0/drmag;
    drdt = DotProduct(dv,dr,ndim)*invdrmag;

    // Add contribution to main star array
    for (k=0; k<ndim; k++) star->a[k] += sphdata[j].m*dr[k]*pow(invdrmag,3);
    for (k=0; k<ndim; k++) star->adot[k] +=
      sphdata[j].m*pow(invdrmag,3)*(dv[k] - 3.0*drdt*invdrmag*dr[k]);
    star->gpot += sphdata[j].m*invdrmag;
    
  }
  //---------------------------------------------------------------------------


  return;
}



//=============================================================================
//  NbodyHermite4::CalculateAllStartupQuantities
/// Calculate all initial properties that are required before the first 
/// timestep integration, i.e. 2nd and 3rd acceleration time derivatives.
//=============================================================================
template <int ndim, template<int> class kernelclass>
void NbodyHermite4<ndim, kernelclass>::CalculateAllStartupQuantities
(int N,                             ///< Number of stars
 NbodyParticle<ndim> **star)        ///< Array of stars/systems
{
  int i,j,k;                        // Star and dimension counters
  DOUBLE a[ndim];                   // Acceleration
  DOUBLE adot[ndim];                // 1st time derivative of accel (jerk)
  DOUBLE a2dot[ndim];               // 2nd time deriivative of acceleration
  DOUBLE afac,bfac,cfac;            // Aux. summation variables
  DOUBLE da[ndim];                  // Relative acceleration
  DOUBLE dadot[ndim];               // Relative jerk
  DOUBLE dr[ndim];                  // Relative position vector
  DOUBLE drdt;                      // Rate of change of distance
  DOUBLE drsqd;                     // Distance squared
  DOUBLE dv[ndim];                  // Relative velocity vector
  DOUBLE invdrmag;                  // 1 / drmag
  DOUBLE invdrsqd;                  // 1 / drsqd
  DOUBLE dvsqd;                     // Velocity squared

  debug2("[NbodyHermite4::CalculateAllStartupQuantities]");


  // Loop over all stars
  //---------------------------------------------------------------------------
  for (i=0; i<N; i++) {

    for (k=0; k<ndim; k++) star[i]->a2dot[k] = 0.0;
    for (k=0; k<ndim; k++) star[i]->a3dot[k] = 0.0;

    // Sum grav. contributions for all other stars (excluding star itself)
    //-------------------------------------------------------------------------
    for (j=0; j<N; j++) {
      if (i == j) continue;

      for (k=0; k<ndim; k++) dr[k] = star[j]->r[k] - star[i]->r[k];
      for (k=0; k<ndim; k++) dv[k] = star[j]->v[k] - star[i]->v[k];
      for (k=0; k<ndim; k++) da[k] = star[j]->a[k] - star[i]->a[k];
      for (k=0; k<ndim; k++) dadot[k] = star[j]->adot[k] - star[i]->adot[k];
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
      for (k=0; k<ndim; k++) star[i]->a3dot[k] += 
        star[j]->m*dadot[k]*invdrsqd*invdrmag -
        9.0*afac*a2dot[k] - 9.0*bfac*adot[k] - 3.0*cfac*a[k];
      
    }
    //-------------------------------------------------------------------------

  }
  //---------------------------------------------------------------------------

  return;
}



//=============================================================================
//  NbodyHermite4::AdvanceParticles
/// Integrate particle positions to 3nd order, and particle velocities to 2nd
/// order from the beginning of the step to the current simulation time, i.e. 
/// r(t+dt) = r(t) + v(t)*dt + a(t)*dt^2/2 + adot(t)*dt^3/6, 
/// v(t+dt) = v(t) + a(t)*dt + adot(t)*dt^2/2.
/// Also set particles at the end of step as 'active' in order to compute 
/// the end-of-step force computation.
//=============================================================================
template <int ndim, template<int> class kernelclass>
void NbodyHermite4<ndim, kernelclass>::AdvanceParticles
(int n,                             ///< Integer time
 int N,                             ///< No. of stars/systems
 NbodyParticle<ndim> **star,        ///< Main star/system array
 DOUBLE timestep)                   ///< Smallest timestep value
{
  int dn;                           // Integer time since beginning of step
  int i;                            // Particle counter
  int k;                            // Dimension counter
  int nstep;                        // Particle (integer) step size
  DOUBLE dt;                        // Timestep since start of step

  debug2("[NbodyHermite4::AdvanceParticles]");

  // Advance positions and velocities of all systems
  //---------------------------------------------------------------------------
  for (i=0; i<N; i++) {

    // Compute time since beginning of step
    nstep = star[i]->nstep;
    dn = n - star[i]->nlast;
    dt = timestep*(FLOAT) dn;

    // Advance positions to third order and velocities to second order
    for (k=0; k<ndim; k++) star[i]->r[k] = star[i]->r0[k] +
      star[i]->v0[k]*dt + 0.5*star[i]->a0[k]*dt*dt +
      onesixth*star[i]->adot0[k]*dt*dt*dt;
    for (k=0; k<vdim; k++) star[i]->v[k] = star[i]->v0[k] +
      star[i]->a0[k]*dt + 0.5*star[i]->adot0[k]*dt*dt;

    // If at end of step, set system particle as active
    if (dn == nstep) star[i]->active = true;
    else star[i]->active = false;
  }
  //---------------------------------------------------------------------------

  return;
}
 


//=============================================================================
//  NbodyHermite4::CorrectionTerms
/// Compute 2nd and 3rd time derivatives of the acceleration using Hermite 
/// interpolation.  Finally correct positions to 5th order and velocities to 
/// 4th order using higher-order derivatives.
//=============================================================================
template <int ndim, template<int> class kernelclass>
void NbodyHermite4<ndim, kernelclass>::CorrectionTerms
(int n,                             ///< Integer time
 int N,                             ///< No. of stars/systems
 NbodyParticle<ndim> **star,        ///< Main star/system array
 DOUBLE timestep)                   ///< Smallest timestep value
{
  int dn;                           // Integer time since beginning of step
  int i;                            // Particle counter
  int k;                            // Dimension counter
  int nstep;                        // Particle (integer) step size
  DOUBLE dt;                        // Physical time step size
  DOUBLE invdt;                     // 1 / dt

  debug2("[NbodyHermite4::CorrectionTerms]");

  // Loop over all system particles
  //---------------------------------------------------------------------------
  for (i=0; i<N; i++) {
    dn = n - star[i]->nlast;
    nstep = star[i]->nstep;

    if (dn == nstep) {
      dt = timestep*(DOUBLE) nstep;
      invdt = 1.0 / dt;
    
      for (k=0; k<ndim; k++) {
        star[i]->a2dot[k] = 
          (-6.0*(star[i]->a0[k] - star[i]->a[k]) -
           dt*(4.0*star[i]->adot0[k] + 2.0*star[i]->adot[k]))*invdt*invdt;
        star[i]->a3dot[k] = 
          (12.0*(star[i]->a0[k] - star[i]->a[k]) + 6.0*dt*
           (star[i]->adot0[k] + star[i]->adot[k]))*invdt*invdt*invdt;
      }

      for (k=0; k<ndim; k++) {
        star[i]->r[k] += star[i]->a2dot[k]*dt*dt*dt*dt/24.0 +
          star[i]->a3dot[k]*dt*dt*dt*dt*dt/120.0;
        star[i]->v[k] += star[i]->a2dot[k]*dt*dt*dt/6.0 +
          star[i]->a3dot[k]*dt*dt*dt*dt/24.0;
      }
    }
    
  }
  //---------------------------------------------------------------------------

  return;
}



//=============================================================================
//  NbodyHermite4::PerturberCorrectionTerms
/// ..
//=============================================================================
template <int ndim, template<int> class kernelclass>
void NbodyHermite4<ndim, kernelclass>::PerturberCorrectionTerms
(int n,                             ///< Integer time
 int N,                             ///< No. of stars/systems
 NbodyParticle<ndim> **star,        ///< Main star/system array
 DOUBLE timestep)                   ///< Smallest timestep value
{
  int dn;                           // Integer time since beginning of step
  int i;                            // Particle counter
  int k;                            // Dimension counter
  int nstep;                        // Particle (integer) step size
  DOUBLE dt;                        // Physical time step size
  DOUBLE invdt;                     // 1 / dt

  debug2("[NbodyHermite4::PerturberCorrectionTerms]");

  // Loop over all system particles
  //---------------------------------------------------------------------------
  for (i=0; i<N; i++) {
    dn = n - star[i]->nlast;
    nstep = star[i]->nstep;
    
    if (dn == nstep) {
      dt = timestep*(DOUBLE) nstep;
      invdt = 1.0 / dt;
      
      for (k=0; k<ndim; k++) star[i]->a[k] += star[i]->apert[k]*invdt;
      for (k=0; k<ndim; k++) star[i]->adot[k] += star[i]->adotpert[k]*invdt;
    }
    
  }
  //---------------------------------------------------------------------------

  return;
}



//=============================================================================
//  NbodyHermite4::EndTimestep
/// Record all important star particle quantities at the end of the step 
/// for the start of the new timestep.
//=============================================================================
template <int ndim, template<int> class kernelclass>
void NbodyHermite4<ndim, kernelclass>::EndTimestep
(int n,                             ///< Integer time
 int N,                             ///< No. of stars/systems
 NbodyParticle<ndim> **star)        ///< Main star/system array
{
  int dn;                           // Integer time since beginning of step
  int i;                            // Particle counter
  int k;                            // Dimension counter
  int nstep;                        // Particle (integer) step size

  debug2("[NbodyHermite4::EndTimestep]");

  // Loop over all system particles
  //---------------------------------------------------------------------------
  for (i=0; i<N; i++) {
    dn = n - star[i]->nlast;
    nstep = star[i]->nstep;

    // If at end of the current step, set quantites for start of new step
    if (dn == nstep) {
      for (k=0; k<ndim; k++) star[i]->r0[k] = star[i]->r[k];
      for (k=0; k<ndim; k++) star[i]->v0[k] = star[i]->v[k];
      for (k=0; k<ndim; k++) star[i]->a0[k] = star[i]->a[k];
      for (k=0; k<ndim; k++) star[i]->adot0[k] = star[i]->adot[k];
      for (k=0; k<ndim; k++) star[i]->apert[k] = 0.0;
      for (k=0; k<ndim; k++) star[i]->adotpert[k] = 0.0;
      star[i]->active = false;
      star[i]->nlast = n;
    }

  }
  //---------------------------------------------------------------------------

  return;
}



//=============================================================================
//  NbodyHermite4::Timestep
/// Calculate the N-body timestep for a given star using the standard 
/// Aarseth timestep, i.e. 
/// $dt = gamma*\sqrt{\frac{a*a2 + a1^2}{a1*a3 + a2^2}}$.
//=============================================================================
template <int ndim, template<int> class kernelclass>
DOUBLE NbodyHermite4<ndim, kernelclass>::Timestep
(NbodyParticle<ndim> *star)         ///< Reference to star/system particle
{
  DOUBLE timestep;                  // Minimum value of particle timesteps
  DOUBLE asqd;                      // Magnitude of particle acceleration
  DOUBLE a1sqd;                     // Magnitude of particle acceleration
  DOUBLE a2sqd;                     // Magnitude of particle acceleration
  DOUBLE a3sqd;                     // Magnitude of particle acceleration

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
  else if (asqd > small_number_dp)
    timestep = sqrt(star->h/(sqrt(asqd) + small_number_dp));
  else
    timestep = big_number_dp;

  timestep = min(timestep,star->dt_internal);

  return timestep;
}



//=============================================================================
//  NbodyHermite4::IntegrateInternalMotion
/// This function integrates the internal motion of a system. First integrates
/// the internal motion of its sub-systems by recursively calling their method,
/// then integrates the COM of the sub-systems.
//=============================================================================
template <int ndim, template<int> class kernelclass>
void NbodyHermite4<ndim, kernelclass>::IntegrateInternalMotion
(SystemParticle<ndim>* systemi,     ///< [inout] System that we wish to 
                                    ///<         integrate the internal motion
 int n,                             ///< [in]    ...
 DOUBLE timestep,                   ///< [in]    ...
 DOUBLE tlocal_end)                 ///< [in]    Time to integrate the 
                                    ///<         internal motion for.
{
  int i;                                              // Particle counter
  int k;                                              // Dimension counter
  int Nchildren = systemi->Nchildren;                 // No. of child systems
  int nsteps_local = 0;                               // No. of local steps
  int nlocal = 0;                                     // Local value of n
  DOUBLE dt;                                          // Timestep
  DOUBLE tlocal=0.0;                                  // Local time
  DOUBLE msystem=0.0;                                 // Mass of system
  DOUBLE rcom[ndim];                                  // Position of COM
  DOUBLE vcom[ndim];                                  // Velocity of COM
  DOUBLE acom[ndim];                                  // Acceleration of COM
  DOUBLE adotcom[ndim];                               // Jerk of COM
  NbodyParticle<ndim>** children = systemi->children; // Child systems

  //cout << "Integrating internal motion : " << Nchildren << "   " << tlocal_end << endl;

  // Zero all COM summation variables
  for (k=0; k<ndim; k++) rcom[k] = 0.0;
  for (k=0; k<ndim; k++) vcom[k] = 0.0;
  for (k=0; k<ndim; k++) acom[k] = 0.0;
  for (k=0; k<ndim; k++) adotcom[k] = 0.0;

  // Make local copies of children and calculate COM properties
  //---------------------------------------------------------------------------
  for (i=0; i<Nchildren; i++) {
    msystem += children[i]->m;
    for (k=0; k<ndim; k++) rcom[k] += children[i]->m*children[i]->r[k];
    for (k=0; k<ndim; k++) vcom[k] += children[i]->m*children[i]->v[k];
    for (k=0; k<ndim; k++) acom[k] += children[i]->m*children[i]->a[k];
    for (k=0; k<ndim; k++) adotcom[k] += children[i]->m*children[i]->adot[k];
  }

  // Normalise COM values
  for (k=0; k<ndim; k++) rcom[k] /= msystem;
  for (k=0; k<ndim; k++) vcom[k] /= msystem;
  for (k=0; k<ndim; k++) acom[k] /= msystem;
  for (k=0; k<ndim; k++) adotcom[k] /= msystem;


  //cout << "Initial system COM : " << rcom[0] << "    " << rcom[1] << endl;

  // Now convert to COM frame
  //---------------------------------------------------------------------------
  for (i=0; i<Nchildren; i++) {
    for (k=0; k<ndim; k++) children[i]->r[k] -= rcom[k];
    for (k=0; k<ndim; k++) children[i]->r0[k] -= rcom[k];
    for (k=0; k<ndim; k++) children[i]->v[k] -= vcom[k];
    for (k=0; k<ndim; k++) children[i]->v0[k] -= vcom[k];
    for (k=0; k<ndim; k++) children[i]->a[k] = 0.0; //acom[k];
    for (k=0; k<ndim; k++) children[i]->a0[k] = 0.0; //acom[k];
    for (k=0; k<ndim; k++) children[i]->adot[k] = 0.0; //adotcom[k];
    for (k=0; k<ndim; k++) children[i]->adot0[k] = 0.0; //adotcom[k];
    children[i]->active = true;
    children[i]->nstep = 1;
    children[i]->nlast = 0;
    children[i]->level = 0;

    /*
    cout << "initial pos : " << i << "   " << children[i]->r0[0] 
	 << "    " << children[i]->r0[1] << "     " 
	 << sqrt(DotProduct(children[i]->r,children[i]->r,ndim)) << endl;
    cout << "initial vel : " << i << "   " << children[i]->v0[0] 
	 << "    " << children[i]->v0[1] << "     "
	 << sqrt(DotProduct(children[i]->v,children[i]->v,ndim)) << endl;
    cout << "initial accel : " << i << "   " << children[i]->a0[0] 
	 << "    " << children[i]->a0[1] << "     "
	 << sqrt(DotProduct(children[i]->a,children[i]->a,ndim)) << endl;
    */

  }


  // Calculate forces, derivatives and other terms
  CalculateDirectGravForces(Nchildren, children);
  for (i=0; i<Nchildren; i++) {
    for (k=0; k<ndim; k++) children[i]->a0[k] = children[i]->a[k];
    for (k=0; k<ndim; k++) children[i]->adot0[k] = children[i]->adot[k];
  }

  // Calculate higher order derivatives
  CalculateAllStartupQuantities(Nchildren, children);


  // Main time integration loop
  //===========================================================================
  do {

    // Calculate global time-step for sub-system
    nlocal = 0;
    dt = std::min(big_number, 1.000000000001*(tlocal_end - tlocal));
    for (i=0; i<Nchildren; i++) {
      children[i]->nlast = 0;
      dt = std::min(dt, Timestep(children[i]));
    }
    tlocal += dt;
    nsteps_local +=1;
    nlocal += 1;

    // Advance position and velocities
    AdvanceParticles(nlocal, Nchildren, children, dt);

    // Zero all acceleration terms
    for (i=0; i<Nchildren; i++) {
      children[i]->active = true;
      children[i]->gpot = 0.0;
      for (k=0; k<ndim; k++) children[i]->a[k] = 0.0;
      for (k=0; k<ndim; k++) children[i]->adot[k] = 0.0;
    }
    
    // Calculate forces, derivatives and other terms
    CalculateDirectGravForces(Nchildren, children);
    
    // Apply correction terms
    CorrectionTerms(nlocal, Nchildren, children, dt);

    //for (i=0; i<Nchildren; i++)
    //cout << "(correction) pos : " << i << "   " << children[i]->r[0] 
    // << "    " << children[i]->r[1] << "     " 
    // << sqrt(DotProduct(children[i]->r,children[i]->r,ndim)) << endl;
    //cin >> i;

    // Now loop over children and, if they are systems, integrate
    // their internal motion
    //-------------------------------------------------------------------------
    for (i=0; i<Nchildren; i++) {
      if (children[i]->Ncomp > 1)
	// The cast is needed because the function is defined only in
	// SystemParticle, not in NbodyParticle.  
	// The safety of the cast relies on the correctness of the Ncomp value
	IntegrateInternalMotion(static_cast<SystemParticle<ndim>* > 
                            (children[i]), n, timestep, dt);
    }

    // Set end-of-step variables
    EndTimestep(nlocal, Nchildren, children);

  } while (tlocal < tlocal_end);
  //===========================================================================


  /*cout << "Done!!" << endl;
  for (i=0; i<Nchildren; i++) {
    cout << "Final pos : " << i << "   " << children[i]->r[0] 
	 << "    " << children[i]->r[1] << "     " 
	 << sqrt(DotProduct(children[i]->r,children[i]->r,ndim)) << endl;
    cout << "Final vel : " << i << "   " << children[i]->v[0] 
	 << "    " << children[i]->v[1] << "     "
	 << sqrt(DotProduct(children[i]->v,children[i]->v,ndim)) << endl;
    cout << "Final accel : " << i << "   " << children[i]->a[0] 
	 << "    " << children[i]->a[1] << "     "
	 << sqrt(DotProduct(children[i]->a,children[i]->a,ndim)) << endl;
	 }*/
 

  // Copy children back to main coordinate system
  //---------------------------------------------------------------------------
  for (i=0; i<Nchildren; i++) {
    for (k=0; k<ndim; k++) children[i]->r[k] += systemi->r[k];
    for (k=0; k<ndim; k++) children[i]->r0[k] += systemi->r[k];
    for (k=0; k<ndim; k++) children[i]->v[k] += systemi->v[k];
    for (k=0; k<ndim; k++) children[i]->v0[k] += systemi->v[k];
    for (k=0; k<ndim; k++) children[i]->a[k] += systemi->a[k];
    for (k=0; k<ndim; k++) children[i]->a0[k] += systemi->a[k];
    for (k=0; k<ndim; k++) children[i]->adot[k] += systemi->adot[k];
    for (k=0; k<ndim; k++) children[i]->adot0[k] += systemi->adot[k];
    children[i]->gpot = children[i]->gpot + systemi->gpot;
  }

  //cin >> i;

  return;
}



//=============================================================================
//  NbodyHermite4::UpdateChildStars
/// ..
//=============================================================================
template <int ndim, template<int> class kernelclass>
void NbodyHermite4<ndim, kernelclass>::UpdateChildStars
(SystemParticle<ndim>* systemi,     ///< [inout] System that we wish to
                                    ///<         integrate the internal motion
 int n,                             ///< [in]    ...
 DOUBLE timestep,                   ///< [in]    ...
 DOUBLE tlocal_end)                 ///< [in]    Time to integrate the
                                    ///<         internal motion for.
{
  int i;                            // ..
  int k;                            // ..
  int Nchildren;                    // ..
  DOUBLE msystot=0.0;               // ..
  DOUBLE rcom[ndim];                // ..
  DOUBLE vcom[ndim];                // ..
  DOUBLE acom[ndim];                // ..
  DOUBLE adotcom[ndim];             // ..
  NbodyParticle<ndim>** children;   // Child systems

  // Only correct positions at end of step
  if (n - systemi->nlast != systemi->nstep) return;

  debug2("[NbodyHermite4::CorrectPerturbedChildStars]");

  // Allocate memory for both stars and perturbers
  Nchildren = systemi->Nchildren;
  children = systemi->children;

  // Set time variables of child stars equal to parent stars
  for (i=0; i<Nchildren; i++) {
    children[i]->nlast = systemi->nlast;
    children[i]->nstep = systemi->nstep;
    children[i]->level = systemi->level;
  }

  // If using perturbers, then correct child star positions to confirm to new 
  // perturbed COM of parent system particle
  //---------------------------------------------------------------------------
  if (perturbers == 1) {
    
    // First calculate old COM
    for (k=0; k<ndim; k++) rcom[k] = 0.0;
    for (k=0; k<ndim; k++) vcom[k] = 0.0;
    for (k=0; k<ndim; k++) acom[k] = 0.0;
    for (k=0; k<ndim; k++) adotcom[k] = 0.0;
    for (i=0; i<Nchildren; i++) {
      msystot += children[i]->m;
      for (k=0; k<ndim; k++) rcom[k] += children[i]->m*children[i]->r[k];
      for (k=0; k<ndim; k++) vcom[k] += children[i]->m*children[i]->v[k];
      for (k=0; k<ndim; k++) acom[k] += children[i]->m*children[i]->a[k];
      for (k=0; k<ndim; k++) adotcom[k] += children[i]->m*children[i]->adot[k];
    }
    for (k=0; k<ndim; k++) rcom[k] /= msystot;
    for (k=0; k<ndim; k++) vcom[k] /= msystot;
    for (k=0; k<ndim; k++) acom[k] /= msystot;
    for (k=0; k<ndim; k++) adotcom[k] /= msystot;
    
    
    // Now translate positions to new COM
    for (i=0; i<Nchildren; i++) {
      for (k=0; k<ndim; k++) children[i]->r[k] += systemi->r[k] - rcom[k];
      for (k=0; k<ndim; k++) children[i]->v[k] += systemi->v[k] - vcom[k];
      for (k=0; k<ndim; k++) children[i]->a[k] += systemi->a[k] - acom[k];
      for (k=0; k<ndim; k++) children[i]->adot[k] += systemi->adot[k] - adotcom[k];
    }
    
    // Now update the positions of any 'grand-children' (i.e. for hierarchies)
    for (i=0; i<Nchildren; i++) {
      if (children[i]->Ncomp > 1)
	UpdateChildStars(static_cast<SystemParticle<ndim>* > (children[i]), 
			 n, timestep, timestep);
    }

  }
  //---------------------------------------------------------------------------

  return;
}



// Template class instances for each dimensionality value (1, 2 and 3) and 
// employed kernel (M4, Quintic, Gaussian and tabulated).
template class NbodyHermite4<1, M4Kernel>;
template class NbodyHermite4<1, QuinticKernel>;
template class NbodyHermite4<1, GaussianKernel>;
template class NbodyHermite4<1, TabulatedKernel>;
template class NbodyHermite4<2, M4Kernel>;
template class NbodyHermite4<2, QuinticKernel>;
template class NbodyHermite4<2, GaussianKernel>;
template class NbodyHermite4<2, TabulatedKernel>;
template class NbodyHermite4<3, M4Kernel>;
template class NbodyHermite4<3, QuinticKernel>;
template class NbodyHermite4<3, GaussianKernel>;
template class NbodyHermite4<3, TabulatedKernel>;

