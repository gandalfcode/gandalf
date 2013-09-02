//=============================================================================
//  NbodyHermite4TS.cpp
//  Contains functions for integrating star particle positions and velocities 
//  using the 4th-order Hermite scheme (Makino & Aarseth 1992) using
//  time-symmetric iterations (Hut et al. 1995??).
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
//  NbodyHermite4TS::NbodyHermite4TS()
/// N-body 4th-order Hermite class constructor
//=============================================================================
template <int ndim, template<int> class kernelclass>
NbodyHermite4TS<ndim, kernelclass>::NbodyHermite4TS
(int nbody_softening_aux, int sub_systems_aux, 
 DOUBLE nbody_mult_aux, string KernelName, int Npec) :
  NbodyHermite4<ndim, kernelclass>(nbody_softening_aux, sub_systems_aux,
                      nbody_mult_aux, KernelName, Npec)
{
}



//=============================================================================
//  NbodyHermite4TS::~NbodyHermite4TS()
/// N-body 4th-order Hermite class destructor
//=============================================================================
template <int ndim, template<int> class kernelclass>
NbodyHermite4TS<ndim, kernelclass>::~NbodyHermite4TS()
{
}



//=============================================================================
//  NbodyHermite4TS::CorrectionTerms
/// Compute 2nd and 3rd time derivatives of the acceleration using Hermite 
/// interpolation.  Finally correct positions to 5th order and velocities to 
/// 4th order using higher-order derivatives.
//=============================================================================
template <int ndim, template<int> class kernelclass>
void NbodyHermite4TS<ndim, kernelclass>::CorrectionTerms
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

  debug2("[NbodyHermite4TS::CorrectionTerms]");

  // Loop over all system particles
  // --------------------------------------------------------------------------
  for (i=0; i<N; i++) {
    dn = n - star[i]->nlast;
    nstep = star[i]->nstep;

    if (dn == nstep) {
      dt = timestep*(DOUBLE) nstep;
      invdt = 1.0 / dt;

      for (k=0; k<ndim; k++) {
        star[i]->a2dot[k] = 
	  (-6.0*(star[i]->a0[k] - star[i]->a[k]) - dt*
	   (4.0*star[i]->adot0[k] + 2.0*star[i]->adot[k]))*invdt*invdt;
        star[i]->a3dot[k] = 
	  (12.0*(star[i]->a0[k] - star[i]->a[k]) + 6.0*dt*
	   (star[i]->adot0[k] + star[i]->adot[k]))*invdt*invdt*invdt;
      }

      for (k=0; k<ndim; k++) {
        star[i]->v[k] = star[i]->v0[k] 
	  + 0.5*(star[i]->a0[k] + star[i]->a[k])*dt
	  //  + 0.1*(star[i]->adot[k] - star[i]->adot0[k])*dt*dt;
	  - onetwelfth*(star[i]->adot[k] - star[i]->adot0[k])*dt*dt;
        star[i]->r[k] = star[i]->r0[k] 
	  + 0.5*(star[i]->v0[k] + star[i]->v[k])*dt
          //+ 0.1*(star[i]->a[k] - star[i]->a0[k])*dt*dt;
          - onetwelfth*(star[i]->a[k] - star[i]->a0[k])*dt*dt;
      }
    }
    
  }
  // --------------------------------------------------------------------------

  return;
}



//=============================================================================
//  NbodyHermite4TS::IntegrateInternalMotion
/// This function integrates the internal motion of a system. First integrates
/// the internal motion of its sub-systems by recursively calling their method,
/// then integrates the COM of the sub-systems.
//=============================================================================
template <int ndim, template<int> class kernelclass>
void NbodyHermite4TS<ndim, kernelclass>::IntegrateInternalMotion
(SystemParticle<ndim>* systemi,     ///< [inout] System that we wish to 
                                    ///<         integrate the internal motion
 DOUBLE tlocal_end)                 ///< [in]    Time to integrate the 
                                    ///<         internal motion for.
{
  int i;                            // Particle counter
  int it;                           // Iteration counter
  int k;                            // Dimension counter
  int Nchildren;                    // No. of child systems
  int Npert;                        // No. of perturbing systems
  int Nstar;                        // Total no. of stars
  int nlocal=0;                     // ..
  int nsteps_local=0;               // ..
  DOUBLE dt;                        // Local timestep
  DOUBLE tlocal=0.0;                // Local time counter
  DOUBLE rcom[ndim];                // Position of COM
  DOUBLE vcom[ndim];                // Velocity of COM
  DOUBLE acom[ndim];                // Acceleration of COM
  DOUBLE adotcom[ndim];             // Jerk of COM
  NbodyParticle<ndim>** children;   // Child systems
  NbodyParticle<ndim>* perturber;   // Local array of perturber properties

  // Allocate memory for both stars and perturbers
  Nchildren = systemi->Nchildren;
  Npert = systemi->Npert;
  Nstar = Nchildren + Npert;
  children = systemi->children;

  //Npert = 0;

  if (Npert > 0) {
    perturber = new NbodyParticle<ndim>[Npert];
    for (i=0; i<Npert; i++) {
      perturber[i].m = systemi->perturber[i]->m;
      for (k=0; k<ndim; k++) perturber[i].r[k] = systemi->perturber[i]->r[k];
      for (k=0; k<ndim; k++) perturber[i].v[k] = systemi->perturber[i]->v[k];
      for (k=0; k<ndim; k++) perturber[i].a[k] = systemi->perturber[i]->a[k];
      for (k=0; k<ndim; k++) 
	perturber[i].adot[k] = systemi->perturber[i]->adot[k];
    }
  }


  // Zero all COM summation variables
  for (k=0; k<ndim; k++) rcom[k] = 0.0;
  for (k=0; k<ndim; k++) vcom[k] = 0.0;
  for (k=0; k<ndim; k++) acom[k] = 0.0;
  for (k=0; k<ndim; k++) adotcom[k] = 0.0;

  // Make local copies of children and calculate COM properties
  // --------------------------------------------------------------------------
  for (i=0; i<Nchildren; i++) {
    for (k=0; k<ndim; k++) rcom[k] += children[i]->m*children[i]->r[k];
    for (k=0; k<ndim; k++) vcom[k] += children[i]->m*children[i]->v[k];
    for (k=0; k<ndim; k++) acom[k] += children[i]->m*children[i]->a[k];
    for (k=0; k<ndim; k++) adotcom[k] += children[i]->m*children[i]->adot[k];
  }

  // Normalise COM values
  for (k=0; k<ndim; k++) rcom[k] /= systemi->m;
  for (k=0; k<ndim; k++) vcom[k] /= systemi->m;
  for (k=0; k<ndim; k++) acom[k] /= systemi->m;
  for (k=0; k<ndim; k++) adotcom[k] /= systemi->m;


  // Now convert to COM frame
  // --------------------------------------------------------------------------
  for (i=0; i<Nchildren; i++) {
    for (k=0; k<ndim; k++) children[i]->r[k] -= rcom[k];
    for (k=0; k<ndim; k++) children[i]->r0[k] -= rcom[k];
    for (k=0; k<ndim; k++) children[i]->v[k] -= vcom[k];
    for (k=0; k<ndim; k++) children[i]->v0[k] -= vcom[k];
    for (k=0; k<ndim; k++) children[i]->a[k] = 0.0; //-= acom[k];
    for (k=0; k<ndim; k++) children[i]->a0[k] = 0.0; //-= acom[k];
    for (k=0; k<ndim; k++) children[i]->adot[k] = 0.0; //-= adotcom[k];
    for (k=0; k<ndim; k++) children[i]->adot0[k] = 0.0; //-= adotcom[k];
    children[i]->gpot = 0.0;
    children[i]->gpe_pert = 0.0;
    children[i]->active = true;
    children[i]->nstep = 1;
    children[i]->nlast = 0;
    children[i]->level = 0;
  }

  // Convert perturbers to COM frame
  for (i=0; i<Npert; i++) {
    for (k=0; k<ndim; k++) perturber[i].r[k] -= rcom[k];
    for (k=0; k<ndim; k++) perturber[i].v[k] -= vcom[k];
    for (k=0; k<ndim; k++) perturber[i].a[k] -= acom[k];
    for (k=0; k<ndim; k++) perturber[i].adot[k] -= adotcom[k];
    for (k=0; k<ndim; k++) perturber[i].r0[k] = perturber[i].r[k];
    for (k=0; k<ndim; k++) perturber[i].v0[k] = perturber[i].v[k];
    for (k=0; k<ndim; k++) perturber[i].a0[k] = perturber[i].a[k];
    for (k=0; k<ndim; k++) perturber[i].adot0[k] = perturber[i].adot[k];
  }


  // Calculate forces, derivatives and other terms
  CalculateDirectGravForces(Nchildren, children);
  if (Npert > 0) 
    CalculatePerturberForces(Nchildren, Npert, children, perturber);

  for (i=0; i<Nchildren; i++) {
    for (k=0; k<ndim; k++) children[i]->a0[k] = children[i]->a[k];
    for (k=0; k<ndim; k++) children[i]->adot0[k] = children[i]->adot[k];
  }

  // Calculate higher order derivatives
  this->CalculateAllStartupQuantities(Nchildren, children);


  // Main time integration loop
  // ==========================================================================
  do {

    // Calculate time-step
    nlocal = 0;
    dt = std::min(big_number, 1.00000000001*(tlocal_end - tlocal));
    for (i=0; i<Nchildren; i++) {
      children[i]->nlast = 0;
      dt = std::min(dt, Timestep(children[i]));
    }
    tlocal += dt;
    nsteps_local +=1;
    nlocal += 1;

    //cout << "Timestep : " << dt << "    tlocal : " << tlocal << endl;

    // Advance position and velocities
    AdvanceParticles(nlocal, Nchildren, children, dt);

    // Advance positions and velocities of perturbers.
    if (Npert > 0) {
      for (i=0; i<Npert; i++) {
	for (k=0; k<ndim; k++) perturber[i].r[k] = perturber[i].r0[k] +
	  perturber[i].v0[k]*tlocal + 0.5*perturber[i].a0[k]*tlocal*tlocal +
	  onesixth*perturber[i].adot0[k]*tlocal;
	for (k=0; k<ndim; k++) perturber[i].v[k] = perturber[i].v0[k] +
	  perturber[i].a0[k]*tlocal + 0.5*perturber[i].adot0[k]*tlocal*tlocal;
	for (k=0; k<ndim; k++) perturber[i].a[k] = perturber[i].a0[k] +
	  perturber[i].adot0[k]*tlocal;
      }
    }

    // Time-symmetric iteration loop
    // ------------------------------------------------------------------------
    for (it=0; it<Npec; it++) {

      // Zero all acceleration terms
      for (i=0; i<Nchildren; i++) {
        children[i]->active = true;
        children[i]->gpot = 0.0;
        children[i]->gpe_pert = 0.0;
        for (k=0; k<ndim; k++) children[i]->a[k] = 0.0;
        for (k=0; k<ndim; k++) children[i]->adot[k] = 0.0;
      }

      // Calculate forces, derivatives and other terms
      CalculateDirectGravForces(Nchildren, children);

      // Add perturbation terms
      if (Npert > 0) 
	CalculatePerturberForces(Nchildren, Npert, children, perturber);
      
      // Apply correction terms
      CorrectionTerms(nlocal, Nchildren, children, dt);
    }
    // ------------------------------------------------------------------------


    // Now loop over children and, if they are systems, integrate
    // their internal motion
    // ------------------------------------------------------------------------
    for (i=0; i<Nchildren; i++) {

      if (children[i]->Ncomp > 1)
	// The cast is needed because the function is defined only in
	// SystemParticle, not in NbodyParticle.  
	// The safety of the cast relies on the correctness of the Ncomp value
	IntegrateInternalMotion(static_cast<SystemParticle<ndim>* > 
				(children[i]), dt);
    }

    // Set end-of-step variables
    EndTimestep(nlocal, Nchildren, children);

  } while (tlocal < tlocal_end);
  // ==========================================================================


  // Copy children back to main coordinate system
  // --------------------------------------------------------------------------
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
    children[i]->gpe *= children[i]->m;
    children[i]->gpe_internal *= children[i]->gpe;
    children[i]->gpe_pert *= children[i]->m;
  }

  if (Npert > 0) delete[] perturber;

  return;
}



// Template class instances for each dimensionality value (1, 2 and 3) and
// employed kernel (M4, Quintic, Gaussian and tabulated).
template class NbodyHermite4TS<1, M4Kernel>;
template class NbodyHermite4TS<1, QuinticKernel>;
template class NbodyHermite4TS<1, GaussianKernel>;
template class NbodyHermite4TS<1, TabulatedKernel>;
template class NbodyHermite4TS<2, M4Kernel>;
template class NbodyHermite4TS<2, QuinticKernel>;
template class NbodyHermite4TS<2, GaussianKernel>;
template class NbodyHermite4TS<2, TabulatedKernel>;
template class NbodyHermite4TS<3, M4Kernel>;
template class NbodyHermite4TS<3, QuinticKernel>;
template class NbodyHermite4TS<3, GaussianKernel>;
template class NbodyHermite4TS<3, TabulatedKernel>;
