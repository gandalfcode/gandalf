//=============================================================================
//  Nbody.cpp
//  Contains main N-body class functions.
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


#include <algorithm>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <math.h>
#include "Precision.h"
#include "Debug.h"
#include "InlineFuncs.h"
#include "NbodyParticle.h"
#include "StarParticle.h"
#include "SystemParticle.h"
#include "Parameters.h"
#include "SphKernel.h"
#include "Nbody.h"

using namespace std;



//=============================================================================
//  Nbody::Nbody
/// Nbody class constructor
//=============================================================================
template <int ndim>
Nbody<ndim>::Nbody(int nbody_softening_aux, int sub_systems_aux, 
                   DOUBLE nbody_mult_aux, string KernelName, int Npec_aux):
  nbody_softening(nbody_softening_aux),
  sub_systems(sub_systems_aux),
  nbody_mult(nbody_mult_aux),
  kerntab(TabulatedKernel<ndim>(KernelName)),
  Nstar(0),
  Nstarmax(0),
  Nsystem(0),
  Nsystemmax(0),
  Nnbody(0),
  Nnbodymax(0),
  reset_tree(0),
  allocated(false),
  Npec(Npec_aux)
{
}



//=============================================================================
//  Nbody::AllocateMemory
/// Allocate all memory required for stars and N-body system particles.
//=============================================================================
template <int ndim>
void Nbody<ndim>::AllocateMemory(int N)
{
  debug2("[Nbody::AllocateMemory]");

  if (N > Nstarmax) {
    if (allocated) DeallocateMemory();
    Nstarmax = N;
    //Nsystem = N;
    Nsystemmax = N;
    //Nnbody = N;
    Nnbodymax = Nstarmax + Nsystemmax;
    nbodydata = new struct NbodyParticle<ndim>*[Nnbodymax];
    stardata = new struct StarParticle<ndim>[Nstarmax];
    system = new struct SystemParticle<ndim>[Nsystemmax];
    allocated = true;
  }

  return;
}
 


//=============================================================================
//  Nbody::DeallocateMemory
/// Deallocate all N-body memory.
//=============================================================================
template <int ndim>
void Nbody<ndim>::DeallocateMemory(void)
{
  debug2("[Nbody::DeallocateMemory]");

  if (allocated) {
    delete[] system;
    delete[] stardata;
    delete[] nbodydata;
  }
  allocated = false;

  return;
}



//=============================================================================
//  Nbody::IntegrateInternalMotion
/// This function integrates the internal motion of a system. First integrates
/// the internal motion of its sub-systems by recursively calling their method,
/// then integrates the COM of the sub-systems.
//=============================================================================
template <int ndim>
void Nbody<ndim>::IntegrateInternalMotion
(SystemParticle<ndim>* systemi,     ///< [inout] System that we wish to 
                                    ///<         integrate the internal motion
 int n,                             ///< [in]    Integer time
 DOUBLE timestep,                   ///< [in]    Minimum timestep value
 DOUBLE tlocal_end)                 ///< [in]    Time to integrate the 
                                    ///<         internal motion for.
{
  int i;                            // Star counter
  int it;                           // Iteration counter
  int k;                            // Dimension counter
  int Nchildren;                    // No. of child systems
  int nlocal_steps = 0;             // No. of locally integrated steps
  DOUBLE dt;                        // Timestep
  DOUBLE tlocal=0.0;                // Local integration time
  DOUBLE rcom[ndim];                // Position of centre-of-mass
  DOUBLE vcom[ndim];                // Velocity of centre-of-mass
  DOUBLE acom[ndim];                // Acceleration of centre-of-mass
  DOUBLE adotcom[ndim];             // Jerk of centre-of-mass
  NbodyParticle<ndim>** children;   // Child systems


  Nchildren = systemi->Nchildren;
  children = systemi->children;

  // Zero all COM summation variables
  for (k=0; k<ndim; k++) rcom[k] = 0.0;
  for (k=0; k<ndim; k++) vcom[k] = 0.0;
  for (k=0; k<ndim; k++) acom[k] = 0.0;
  for (k=0; k<ndim; k++) adotcom[k] = 0.0;

  // Make local copies of children and calculate COM properties
  //---------------------------------------------------------------------------
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
  //---------------------------------------------------------------------------
  for (i=0; i<Nchildren; i++) {
    for (k=0; k<ndim; k++) children[i]->r[k] -= rcom[k];
    for (k=0; k<ndim; k++) children[i]->r0[k] -= rcom[k];
    for (k=0; k<ndim; k++) children[i]->v[k] -= vcom[k];
    for (k=0; k<ndim; k++) children[i]->v0[k] -= vcom[k];
    for (k=0; k<ndim; k++) children[i]->a[k] -= acom[k];
    for (k=0; k<ndim; k++) children[i]->a0[k] -= acom[k];
    for (k=0; k<ndim; k++) children[i]->adot[k] -= adotcom[k];
    for (k=0; k<ndim; k++) children[i]->adot0[k] -= adotcom[k];
    for (k=0; k<ndim; k++) children[i]->a2dot[k] -= 0.0;
    for (k=0; k<ndim; k++) children[i]->a3dot[k] -= 0.0;
    children[i]->active = true;
    children[i]->nstep = 1;
    children[i]->level = 0;
  }


  // Main time integration loop
  //===========================================================================
  do {

    // Calculate time-step
    dt = std::min(big_number, tlocal_end - tlocal);
    for (i=0; i<Nchildren; i++) {
      dt = std::min(dt, Timestep(children[i]));
    }
    tlocal += dt;
    nlocal_steps +=1;

    // Advance position and velocities
    AdvanceParticles(nlocal_steps, Nchildren, children, dt);

    // Time-symmetric iteration loop
    //-------------------------------------------------------------------------
    for (it=0; it<Npec; it++) {

      //Zero all acceleration terms
      for (i=0; i<Nchildren; i++) {
        children[i]->gpot = 0.0;
        children[i]->gpe_internal = 0.0;
        for (k=0; k<ndim; k++) children[i]->a[k] = 0.0;
        for (k=0; k<ndim; k++) children[i]->adot[k] = 0.0;
        for (k=0; k<ndim; k++) children[i]->a2dot[k] = 0.0;
        for (k=0; k<ndim; k++) children[i]->a3dot[k] = 0.0;
      }

      // Calculate forces, derivatives and other terms
      CalculateDirectGravForces(Nchildren, children);
      
      // Apply correction terms
      CorrectionTerms(nlocal_steps, Nchildren, children, dt);
    }
    //-------------------------------------------------------------------------

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
    EndTimestep(nlocal_steps, Nchildren, children);

  } while (tlocal < tlocal_end);
  //===========================================================================


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

  return;
}



template class Nbody<1>;
template class Nbody<2>;
template class Nbody<3>;
