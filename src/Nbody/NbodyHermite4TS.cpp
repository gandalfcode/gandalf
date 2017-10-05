//=================================================================================================
//  NbodyHermite4TS.cpp
//  Contains functions for integrating star particle positions and velocities
//  using the 4th-order Hermite scheme (Makino & Aarseth 1992) using
//  time-symmetric iterations (Hut et al. 1995??).
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
using namespace std;



//=================================================================================================
//  NbodyHermite4TS::NbodyHermite4TS()
/// N-body 4th-order time-symmetric Hermite class constructor
//=================================================================================================
template <int ndim, template<int> class kernelclass>
NbodyHermite4TS<ndim, kernelclass>::NbodyHermite4TS
(int nbody_softening_aux, int _perturbers, int sub_systems_aux,
 DOUBLE nbody_mult_aux, string KernelName, int Npec) :
  NbodyHermite4<ndim, kernelclass>(nbody_softening_aux, _perturbers, sub_systems_aux,
                                   nbody_mult_aux, KernelName, Npec)
{
}




//=================================================================================================
//  NbodyHermite4TS::~NbodyHermite4TS()
/// N-body 4th-order Hermite class destructor
//=================================================================================================
template <int ndim, template<int> class kernelclass>
NbodyHermite4TS<ndim, kernelclass>::~NbodyHermite4TS()
{
}



//=================================================================================================
//  NbodyHermite4TS::CorrectionTerms
/// Compute 2nd and 3rd time derivatives of the acceleration using Hermite interpolation.  Finally
/// correct positions to 5th order and velocities to 4th order using higher-order derivatives.
//=================================================================================================
template <int ndim, template<int> class kernelclass>
void NbodyHermite4TS<ndim, kernelclass>::CorrectionTerms
(int n,                             ///< [in] Integer time
 int N,                             ///< [in] No. of stars/systems
 FLOAT t,                          ///< ..
 FLOAT timestep,                   ///< ..
 NbodyParticle<ndim> **star)        ///< [inout] Main star/system array
{
  int dn;                           // Integer time since beginning of step
  int i;                            // Particle counter
  int k;                            // Dimension counter
  int nstep;                        // Particle (integer) step size
  FLOAT dt;                        // Physical time step size
  FLOAT invdt;                     // 1 / dt

  debug2("[NbodyHermite4TS::CorrectionTerms]");

  // Loop over all system particles
  //-----------------------------------------------------------------------------------------------
  for (i=0; i<N; i++) {
    dn = n - star[i]->nlast;
    nstep = star[i]->nstep;

    if (dn == nstep) {
      //dt = timestep*(FLOAT) nstep;
      dt = t - star[i]->tlast;
      invdt = 1.0 / dt;

      for (k=0; k<ndim; k++) {
        star[i]->a2dot[k] = (-6.0*(star[i]->a0[k] - star[i]->a[k]) - dt*
                             (4.0*star[i]->adot0[k] + 2.0*star[i]->adot[k]))*invdt*invdt;
        star[i]->a3dot[k] = (12.0*(star[i]->a0[k] - star[i]->a[k]) + 6.0*dt*
                             (star[i]->adot0[k] + star[i]->adot[k]))*invdt*invdt*invdt;
      }

      for (k=0; k<ndim; k++) {
        star[i]->v[k] = star[i]->v0[k] + 0.5*(star[i]->a0[k] + star[i]->a[k])*dt -
          onetwelfth*(star[i]->adot[k] - star[i]->adot0[k])*dt*dt;
        star[i]->r[k] = star[i]->r0[k] + 0.5*(star[i]->v0[k] + star[i]->v[k])*dt -
          onetwelfth*(star[i]->a[k] - star[i]->a0[k])*dt*dt;
      }
    }

  }
  //-----------------------------------------------------------------------------------------------

  return;
}



//=================================================================================================
//  NbodyHermite4TS::IntegrateInternalMotion
/// This function integrates the internal motion of a system. First integrates the
/// internal motion of its sub-systems by recursively calling their method, then integrates
/// the COM of the sub-systems.
//=================================================================================================
/*template <int ndim, template<int> class kernelclass>
void NbodyHermite4TS<ndim, kernelclass>::IntegrateInternalMotion
(SystemParticle<ndim>* systemi,     ///< [inout] System to integrate the internal motion for
 int n,                             ///< [in]    Integer time
 FLOAT timestep,                   ///< [in]    Smallest timestep
 FLOAT tlocal_end)                 ///< [in]    Time to integrate the internal motion for.
{
  int i;                            // Particle counter
  int it;                           // Iteration counter
  int k;                            // Dimension counter
  int Nchildren;                    // No. of child systems
  int Npert;                        // No. of perturbing systems
  int Nstar;                        // Total no. of stars
  int nlocal=0;                     // Local block step integer time
  int nsteps_local=0;               // Local no. of steps
  FLOAT aext[ndim];                // Acceleration due to external stars
  FLOAT adotext[ndim];             // Jerk due to external stars
  FLOAT dt;                        // Local timestep
  FLOAT tlocal=0.0;                // Local time counter
  FLOAT tpert;                     // ..
  FLOAT *apert;                    // ..
  FLOAT *adotpert;                 // ..
  NbodyParticle<ndim>** children;   // Child systems
  NbodyParticle<ndim>* perturber;   // Local array of perturber properties

  // Only integrate internal motion once COM motion has finished
  if (n - systemi->nlast != systemi->nstep) return;

  debug2("[NbodyHermite4TS::IntegrateInternalMotion]");

  // Allocate memory for both stars and perturbers
  Nchildren = systemi->Nchildren;
  Npert = systemi->Npert;
  Nstar = Nchildren + Npert;
  children = systemi->children;


  // Calculate total external acceleration and jerk terms
  for (k=0; k<ndim; k++) aext[k] = systemi->m*systemi->a0[k];
  for (k=0; k<ndim; k++) adotext[k] = systemi->m*systemi->adot0[k];


  // If using perturbers, record local copies and remove contribution to
  // external acceleration and jerk terms
  //-----------------------------------------------------------------------------------------------
  if (perturbers == 1 && Npert > 0) {
    perturber = new NbodyParticle<ndim>[Npert];
    apert     = new FLOAT[ndim*Npert];
    adotpert  = new FLOAT[ndim*Npert];

    // Create local copies of perturbers
    for (i=0; i<Npert; i++) {
      perturber[i].m = systemi->perturber[i]->m;
      perturber[i].nlast = systemi->perturber[i]->nlast;

      for (k=0; k<ndim; k++) perturber[i].apert[k] = 0.0;
      for (k=0; k<ndim; k++) perturber[i].adotpert[k] = 0.0;
      for (k=0; k<ndim; k++) perturber[i].r0[k] = systemi->perturber[i]->r0[k];
      for (k=0; k<ndim; k++) perturber[i].v0[k] = systemi->perturber[i]->v0[k];
      for (k=0; k<ndim; k++) perturber[i].a0[k] = systemi->perturber[i]->a0[k];
      for (k=0; k<ndim; k++) perturber[i].adot[k] = systemi->perturber[i]->adot0[k];

      // Set properties of perturbers at beginning of current step
      tpert = (FLOAT) (n - 1 - perturber[i].nlast)*timestep;
      for (k=0; k<ndim; k++) perturber[i].r[k] = perturber[i].r0[k] + perturber[i].v0[k]*tpert +
        0.5*perturber[i].a0[k]*tpert*tpert + onesixth*perturber[i].adot0[k]*tpert*tpert*tpert;
      for (k=0; k<ndim; k++) perturber[i].v[k] = perturber[i].v0[k] +
        perturber[i].a0[k]*tpert + 0.5*perturber[i].adot0[k]*tpert*tpert;
      for (k=0; k<ndim; k++) perturber[i].a[k] = perturber[i].a0[k] + perturber[i].adot0[k]*tpert;
      for (k=0; k<ndim; k++) perturber[i].adot[k] = perturber[i].adot0[k];
    }
  }


  // Initialise children properties
  //-----------------------------------------------------------------------------------------------
  for (i=0; i<Nchildren; i++) {
    for (k=0; k<ndim; k++) children[i]->a[k]        = 0.0;
    for (k=0; k<ndim; k++) children[i]->adot[k]     = 0.0;
    for (k=0; k<ndim; k++) children[i]->apert[k]    = 0.0;
    for (k=0; k<ndim; k++) children[i]->adotpert[k] = 0.0;
    children[i]->gpot         = 0.0;
    children[i]->gpe          = 0.0;
    children[i]->gpe_pert     = 0.0;
    children[i]->gpe_internal = 0.0;
    children[i]->flags.set(active);
    children[i]->nstep        = 1;
    children[i]->nlast        = 0;
    children[i]->level        = 0;
  }

  if (perturbers == 1 && Npert > 0) {
    this->CalculatePerturberForces(Nchildren, Npert, children, perturber, apert, adotpert);
    for (i=0; i<Nchildren; i++) {
      for (k=0; k<ndim; k++) aext[k] -= children[i]->m*children[i]->a[k];
      for (k=0; k<ndim; k++) adotext[k] -= children[i]->m*children[i]->adot[k];
    }
  }
  for (k=0; k<ndim; k++) aext[k] /= systemi->m;
  for (k=0; k<ndim; k++) adotext[k] /= systemi->m;


  // Calculate forces, derivatives and other terms
  CalculateDirectGravForces(Nchildren, children);

  for (i=0; i<Nchildren; i++) {
    for (k=0; k<ndim; k++) children[i]->a[k]    += aext[k];
    for (k=0; k<ndim; k++) children[i]->adot[k] += adotext[k];
    for (k=0; k<ndim; k++) children[i]->a0[k]    = children[i]->a[k];
    for (k=0; k<ndim; k++) children[i]->adot0[k] = children[i]->adot[k];
  }

  // Calculate higher order derivatives
  this->CalculateAllStartupQuantities(Nchildren, children);


  // Main time integration loop
  //===============================================================================================
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

    // Advance position and velocities
    AdvanceParticles(nlocal, Nchildren, children, dt);

    // Advance positions and velocities of perturbers
    if (perturbers == 1 && Npert > 0) {
      for (i=0; i<Npert; i++) {
        tpert = tlocal + (FLOAT) (n - 1 - perturber[i].nlast)*timestep;
        for (k=0; k<ndim; k++) perturber[i].r[k] = perturber[i].r0[k] +
          perturber[i].v0[k]*tpert + 0.5*perturber[i].a0[k]*tpert*tpert +
          onesixth*perturber[i].adot0[k]*tpert*tpert*tpert;
        for (k=0; k<ndim; k++) perturber[i].v[k] = perturber[i].v0[k] +
          perturber[i].a0[k]*tpert + 0.5*perturber[i].adot0[k]*tpert*tpert;
        for (k=0; k<ndim; k++) perturber[i].a[k] = perturber[i].a0[k] +
          perturber[i].adot0[k]*tpert;
      }
    }


    // Time-symmetric iteration loop
    //---------------------------------------------------------------------------------------------
    for (it=0; it<Npec; it++) {

      // Zero all acceleration terms
      for (i=0; i<Nchildren; i++) {
        children[i]->flags.set(active);
        children[i]->gpot = 0.0;
        children[i]->gpe_pert = 0.0;
        for (k=0; k<ndim; k++) children[i]->a[k] = aext[k];
        for (k=0; k<ndim; k++) children[i]->adot[k] = adotext[k];
        for (k=0; k<ndim; k++) children[i]->apert[k] = 0.0;
        for (k=0; k<ndim; k++) children[i]->adotpert[k] = 0.0;
      }

      // Calculate forces, derivatives and other terms
      CalculateDirectGravForces(Nchildren, children);

      // Add perturbation terms
      if (perturbers == 1 && Npert > 0) {
        this->CalculatePerturberForces(Nchildren, Npert, children, perturber, apert, adotpert);
      }

      // Apply correction terms
      CorrectionTerms(nlocal, Nchildren, children, dt);

    }
    //---------------------------------------------------------------------------------------------

    // Add perturber forces to local arrays
    if (perturbers == 1 && Npert > 0) {
      for (i=0; i<Npert; i++) {
        for (k=0; k<ndim; k++) perturber[i].apert[k] += apert[ndim*i + k];
        for (k=0; k<ndim; k++) perturber[i].adotpert[k] += adotpert[ndim*i + k];
      }
    }

    // Now loop over children and, if they are systems, integrate
    // their internal motion
    //---------------------------------------------------------------------------------------------
    for (i=0; i<Nchildren; i++) {
      if (children[i]->Ncomp > 1)
        // The cast is needed because the function is defined only in SystemParticle, not in
        // NbodyParticle.  The safety of the cast relies on the correctness of the Ncomp value.
        IntegrateInternalMotion(static_cast<SystemParticle<ndim>* >
                                (children[i]), n, timestep, dt);
    }

    // Calculate correction terms on perturbing stars due to sub-systems
    if (perturbers == 1 && Npert > 0) {
      this->PerturberCorrectionTerms(nlocal, Nchildren, children, dt);
      CorrectionTerms(nlocal, Nchildren, children, dt);
    }

    // Correct positions of all child stars in any hierarchical systems
    for (i=0; i<Nchildren; i++) {
      if (children[i]->Ncomp > 1) {
        this->UpdateChildStars(static_cast<SystemParticle<ndim>* > (children[i]),
                               n, timestep, dt);
      }
    }

    // Set end-of-step variables
    EndTimestep(nlocal, Nchildren, children);


  } while (tlocal < tlocal_end);
  //===============================================================================================


  // Copy children back to main coordinate system
  //-----------------------------------------------------------------------------------------------
  systemi->gpe_internal = 0.0;
  //systemi->gpe_pert = 0.0;
  for (i=0; i<Nchildren; i++) {
    systemi->gpe_internal += 0.5*children[i]->m*children[i]->gpot;
    //systemi->gpe_pert += children[i]->gpe_pert;
    children[i]->gpot = children[i]->gpot + systemi->gpot;
    children[i]->gpe = children[i]->gpot*children[i]->m;
    //children[i]->gpe_internal *= children[i]->gpe;
    //children[i]->gpe_pert *= children[i]->m;
  }

  for (k=0; k<ndim; k++) systemi->r[k] = 0.0;
  for (k=0; k<ndim; k++) systemi->v[k] = 0.0;
  for (k=0; k<ndim; k++) systemi->a[k] = 0.0;
  for (k=0; k<ndim; k++) systemi->adot[k] = 0.0;
  for (i=0; i<Nchildren; i++) {
    for (k=0; k<ndim; k++) systemi->r[k] += children[i]->m*children[i]->r[k];
    for (k=0; k<ndim; k++) systemi->v[k] += children[i]->m*children[i]->r[k];
    for (k=0; k<ndim; k++) systemi->a[k] += children[i]->m*children[i]->a[k];
    for (k=0; k<ndim; k++) systemi->adot[k] += children[i]->m*children[i]->adot[k];
  }
  for (k=0; k<ndim; k++) systemi->r[k] /= systemi->m;
  for (k=0; k<ndim; k++) systemi->v[k] /= systemi->m;
  for (k=0; k<ndim; k++) systemi->a[k] /= systemi->m;
  for (k=0; k<ndim; k++) systemi->adot[k] /= systemi->m;

  //cout << "FINAL SYSTEM M : " << systemi->m << endl;
  //cout << "FINAL AEXT     : " << aext[0] << "    " << aext[1] << endl;
  //cout << "FINAL SYSTEM R : " << systemi->r[0] << "    " << systemi->r[1] << endl;
  //cout << "FINAL SYSTEM A : " << systemi->a[0] << "    " << systemi->a[1] << endl;

  // Finally, add perturbations on perturber itself to main arrays before
  // deallocating local memory
  if (perturbers == 1 && Npert > 0) {
    for (i=0; i<Npert; i++) {
      for (k=0; k<ndim; k++) systemi->perturber[i]->apert[k] += perturber[i].apert[k];
      for (k=0; k<ndim; k++) systemi->perturber[i]->adotpert[k] += perturber[i].adotpert[k];
    }
    delete[] adotpert;
    delete[] apert;
    delete[] perturber;
  }
  //cin >> i;
  return;
}*/



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
