//=============================================================================
//  NbodyLeapfrogDKD.cpp
//  Contains functions for integrating star particle positions and velocities 
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
#include <cassert>
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
#include "Sph.h"
using namespace std;



//=============================================================================
//  NbodyLeapfrogDKD::NbodyLeapfrogDKD()
/// N-body leapfrog KDK class constructor
//=============================================================================
template <int ndim, template<int> class kernelclass>
NbodyLeapfrogDKD<ndim, kernelclass>::NbodyLeapfrogDKD
(int nbody_softening_aux, int sub_systems_aux, 
 DOUBLE nbody_mult_aux, string KernelName) : 
  Nbody<ndim>(nbody_softening_aux, sub_systems_aux, 
              nbody_mult_aux, KernelName, 1),
  kern(kernelclass<ndim>(KernelName))
{
}



//=============================================================================
//  NbodyLeapfrogDKD::~NbodyLeapfrog()
/// N-body leapfrog KDK class destructor
//=============================================================================
template <int ndim, template<int> class kernelclass>
NbodyLeapfrogDKD<ndim, kernelclass>::~NbodyLeapfrogDKD()
{
}



//=============================================================================
//  NbodyLeapfrogDKD::CalculateDirectGravForces
/// Calculate all star-star force contributions for active systems using 
/// direct summation with unsoftened gravity.
//=============================================================================
template <int ndim, template<int> class kernelclass>
void NbodyLeapfrogDKD<ndim, kernelclass>::CalculateDirectGravForces
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
//  NbodyLeapfrogDKD::CalculatePerturberForces
/// Empty function for now.
//=============================================================================
template <int ndim, template<int> class kernelclass>
void NbodyLeapfrogDKD<ndim, kernelclass>::CalculatePerturberForces
(int N,                             ///< Number of stars
 int Npert,                         ///< Number of perturbing stars
 NbodyParticle<ndim> **star,        ///< Array of stars/systems
 NbodyParticle<ndim> *perturber,    ///< Array of perturbing stars/systems
 DOUBLE *apert,                     ///< Current sub-system timestep
 DOUBLE *adotpert)                  ///< ..
{
  return;
}



//=============================================================================
//  NbodyLeapfrogDKD::CalculateDirectSPHForces
/// Calculate all forces due to SPH neighbours and distant SPH particles 
/// due to direct-sum gravity.
//=============================================================================
template <int ndim, template<int> class kernelclass>
void NbodyLeapfrogDKD<ndim, kernelclass>::CalculateDirectSPHForces
(NbodyParticle<ndim> *star,         ///< [inout] Pointer to star
 int Nsph,                          ///< [in] No. of SPH neighbour gas ptcls
 int Ndirect,                       ///< [in] No. of distant SPH ptcls.
 int *sphlist,                      ///< [in] List of neighbour ids
 int *directlist,                   ///< [in] List of distant ptcl ids
 Sph<ndim> * sph)        ///< [in] Array of SPH particles
{
  int j,jj,k;                       // Star and dimension counters
  DOUBLE dr[ndim];                  // Relative position vector
  DOUBLE drmag;                     // Distance
  DOUBLE drsqd;                     // Distance squared
  DOUBLE invhmean;                  // 1 / hmean
  DOUBLE invdrmag;                  // 1 / drmag
  DOUBLE paux;                      // Aux. force variable

  debug2("[NbodyLeapfrogDKD::CalculateDirectSPHForces]");


  // Sum grav. contributions from all neighbouring SPH particles
  //---------------------------------------------------------------------------
  for (jj=0; jj<Nsph; jj++) {
    j = sphlist[jj];

    SphParticle<ndim>& part = sph->GetParticleIPointer(j);

    assert(part.itype != dead);

    for (k=0; k<ndim; k++) dr[k] = part.r[k] - star->r[k];
    drsqd = DotProduct(dr,dr,ndim);
    drmag = sqrt(drsqd);
    invdrmag = 1.0/drmag;
    invhmean = 2.0/(star->h + part.h);
   
    //for (k=0; k<ndim; k++) star->a[k] += part.m*dr[k]*pow(invdrmag,3);
    //star->gpot += part.m*invdrmag;
    //continue;
 
    paux = part.m*invhmean*invhmean*
      kern.wgrav(drmag*invhmean)*invdrmag;
    
    // Add contribution to main star array
    for (k=0; k<ndim; k++) star->a[k] += paux*dr[k];
    star->gpot += part.m*invhmean*kern.wpot(drmag*invhmean);

  }
  //---------------------------------------------------------------------------


  // Now include contributions from distant, non-SPH neighbours
  // (i.e. direct summation with Newton's law of gravity)
  //---------------------------------------------------------------------------
  for (jj=0; jj<Ndirect; jj++) {
    j = directlist[jj];

    SphParticle<ndim>& part = sph->GetParticleIPointer(j);
    //assert(sphdata[j].itype != dead);

    for (k=0; k<ndim; k++) dr[k] = part.r[k] - star->r[k];
    drsqd = DotProduct(dr,dr,ndim);
    drmag = sqrt(drsqd);
    invdrmag = 1.0/drmag;
    paux = part.m*pow(invdrmag,3);

    // Add contribution to main star array
    for (k=0; k<ndim; k++) star->a[k] += paux*dr[k];
    star->gpot += part.m*invdrmag;
    
  }
  //---------------------------------------------------------------------------


  return;
}



//=============================================================================
//  NbodyLeapfrogDKD::CalculateAllStartupQuantities
/// Empty function for Leapfrog-KDK (no additional quantities required).
//=============================================================================
template <int ndim, template<int> class kernelclass>
void NbodyLeapfrogDKD<ndim, kernelclass>::CalculateAllStartupQuantities
(int N,                             ///< Number of stars
 NbodyParticle<ndim> **star)        ///< Array of stars/systems
{
  return;
}



//=============================================================================
//  NbodyLeapfrogDKD::AdvanceParticles
/// Integrate star positions to 2nd order, and star velocities to 1st
/// order from the beginning of the step to the current simulation time, i.e. 
/// $r(t+dt) = r(t) + v(t)*dt + 0.5*a(t)*dt^2$, 
/// $v(t+dt) = v(t) + a(t)*dt$.
/// Also set particles at the end of step as 'active' in order to compute 
/// the end-of-step force computation and velocity correction step.
//=============================================================================
template <int ndim, template<int> class kernelclass>
void NbodyLeapfrogDKD<ndim, kernelclass>::AdvanceParticles
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

  debug2("[NbodyLeapfrogDKD::AdvanceParticles]");

  // Advance positions and velocities of all star particles
  //---------------------------------------------------------------------------
  for (i=0; i<N; i++) {

    // Compute time since beginning of step
    nstep = star[i]->nstep;
    dn = n - star[i]->nlast;
    dt = timestep*(FLOAT) dn;

    // Advance positions and velocities to first order
    if (dn < nstep) { 
      for (k=0; k<ndim; k++) 
	star[i]->r[k] = star[i]->r0[k] + star[i]->v0[k]*dt;
      for (k=0; k<vdim; k++) 
	star[i]->v[k] = star[i]->v0[k] + star[i]->a[k]*dt;
    }
    else {
      for (k=0; k<vdim; k++) 
	star[i]->v[k] = star[i]->v0[k] + star[i]->a[k]*dt;
      for (k=0; k<ndim; k++) star[i]->r[k] = 
	star[i]->r0[k] + 0.5*(star[i]->v0[k] + star[i]->v[k])*dt;
    }

    // If at half-step, set system particle as active
    if (dn == nstep/2) star[i]->active = true;
    else star[i]->active = false;
  }
  //---------------------------------------------------------------------------

  return;
}
 


//=============================================================================
//  NbodyLeapfrogDKD::CorrectionTerms
/// Empty function.  No correction terms for Leapfrog-DKD
//=============================================================================
template <int ndim, template<int> class kernelclass>
void NbodyLeapfrogDKD<ndim, kernelclass>::CorrectionTerms
(int n,                             ///< Integer time
 int N,                             ///< No. of stars/systems
 NbodyParticle<ndim> **star,        ///< Main star/system array
 DOUBLE timestep)                   ///< Smallest timestep value
{
  return;
}



//=============================================================================
//  NbodyLeapfrogDKD::PerturberCorrectionTerms
/// ..
//=============================================================================
template <int ndim, template<int> class kernelclass>
void NbodyLeapfrogDKD<ndim, kernelclass>::PerturberCorrectionTerms
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
     }

   }
   //--------------------------------------------------------------------------

  return;
}



//=============================================================================
//  NbodyLeapfrogDKD::EndTimestep
/// Record all important star particle quantities at the end of the step 
/// for the start of the new timestep.
//=============================================================================
template <int ndim, template<int> class kernelclass>
void NbodyLeapfrogDKD<ndim, kernelclass>::EndTimestep
(int n,                             ///< Integer time
 int N,                             ///< No. of stars/systems
 NbodyParticle<ndim> **star)        ///< Main star/system array
{
  int dn;                           // Integer time since beginning of step
  int i;                            // Particle counter
  int k;                            // Dimension counter
  int nstep;                        // Particle (integer) step size

  debug2("[NbodyLeapfrogDKD::EndTimestep]");

  // Loop over all star particles and set values for those at end of step
  //---------------------------------------------------------------------------
  for (i=0; i<N; i++) {
    dn = n - star[i]->nlast;
    nstep = star[i]->nstep;

    if (dn == nstep) {
      for (k=0; k<ndim; k++) star[i]->r0[k] = star[i]->r[k];
      for (k=0; k<ndim; k++) star[i]->v0[k] = star[i]->v[k];
      for (k=0; k<ndim; k++) star[i]->a0[k] = star[i]->a[k];
      star[i]->active = false;
      star[i]->nlast = n;
    }
  }
  //---------------------------------------------------------------------------

  return;
}



//=============================================================================
//  NbodyLeapfrogDKD::Timestep
/// Default timestep size for SPH particles.  Takes the minimum of : 
/// (i)  const*h/(sound_speed + h*|div_v|)    (Courant condition)
/// (ii) const*sqrt(h/|a|)                    (Acceleration condition)
//=============================================================================
template <int ndim, template<int> class kernelclass>
DOUBLE NbodyLeapfrogDKD<ndim, kernelclass>::Timestep
(NbodyParticle<ndim> *star)         ///< Reference to SPH particle
{
  DOUBLE timestep;                  // Minimum value of particle timesteps
  DOUBLE amag;                      // Magnitude of star acceleration
  DOUBLE adotmag;                   // Magnitude of star jerk

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




//=============================================================================
//  NbodyLeapfrogDKD::IntegrateInternalMotion
/// This function integrates the internal motion of a system. First integrates
/// the internal motion of its sub-systems by recursively calling their method,
/// then integrates the COM of the sub-systems.
//=============================================================================
template <int ndim, template<int> class kernelclass>
void NbodyLeapfrogDKD<ndim, kernelclass>::IntegrateInternalMotion
(SystemParticle<ndim>* systemi,     ///< [inout] System that we wish to 
                                    ///<         integrate the internal motion
 int n,                             ///< [in]    ...
 DOUBLE timestep,                   ///< [in]    ...
 DOUBLE tlocal_end)                 ///< [in]    Time to integrate the 
                                    ///<         internal motion for.
{
  return;
}



//=============================================================================
//  NbodyLeapfrogDKD::UpdateChildStars
/// ..
//=============================================================================
template <int ndim, template<int> class kernelclass>
void NbodyLeapfrogDKD<ndim, kernelclass>::UpdateChildStars
(SystemParticle<ndim>* systemi,     ///< [inout] System that we wish to
                                    ///<         integrate the internal motion
 int n,                             ///< [in]    ...
 DOUBLE timestep,                   ///< [in]    ...
 DOUBLE tlocal_end)                 ///< [in]    Time to integrate the
                                    ///<         internal motion for.
{
  return;
}



// Template class instances for each dimensionality value (1, 2 and 3) and 
// employed kernel (M4, Quintic, Gaussian and tabulated).
template class NbodyLeapfrogDKD<1, M4Kernel>;
template class NbodyLeapfrogDKD<1, QuinticKernel>;
template class NbodyLeapfrogDKD<1, GaussianKernel>;
template class NbodyLeapfrogDKD<1, TabulatedKernel>;
template class NbodyLeapfrogDKD<2, M4Kernel>;
template class NbodyLeapfrogDKD<2, QuinticKernel>;
template class NbodyLeapfrogDKD<2, GaussianKernel>;
template class NbodyLeapfrogDKD<2, TabulatedKernel>;
template class NbodyLeapfrogDKD<3, M4Kernel>;
template class NbodyLeapfrogDKD<3, QuinticKernel>;
template class NbodyLeapfrogDKD<3, GaussianKernel>;
template class NbodyLeapfrogDKD<3, TabulatedKernel>;
