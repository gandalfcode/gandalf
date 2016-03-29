//=================================================================================================
//  MfvMuscl.cpp
//  Contains all functions for calculating Meshless Finite-Volume Hydrodynamics quantities.
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
#include <cassert>
#include <iostream>
#include <math.h>
#include "Precision.h"
#include "MeshlessFV.h"
#include "Particle.h"
#include "Parameters.h"
#include "SmoothingKernel.h"
#include "EOS.h"
#include "Debug.h"
#include "Exception.h"
#include "InlineFuncs.h"
using namespace std;



//=================================================================================================
//  MfvMuscl::MfvMuscl
/// MfvMuscl class constructor.  Calls main SPH class constructor and also
/// sets additional kernel-related quantities
//=================================================================================================
template <int ndim, template<int> class kernelclass>
MfvMuscl<ndim, kernelclass>::MfvMuscl
 (int _hydro_forces, int _self_gravity, FLOAT _accel_mult, FLOAT _courant_mult,
  FLOAT _h_fac, FLOAT _h_converge, FLOAT _gamma, string _gas_eos, string KernelName,
  int size_part, SimUnits &units, Parameters *params):
  MfvCommon<ndim,kernelclass>(_hydro_forces, _self_gravity, _accel_mult, _courant_mult, _h_fac,
                              _h_converge, _gamma, _gas_eos, KernelName, size_part, units, params)
{
}



//=================================================================================================
//  MfvMuscl::~MfvMuscl
/// MfvMuscl class destructor
//=================================================================================================
template <int ndim, template<int> class kernelclass>
MfvMuscl<ndim, kernelclass>::~MfvMuscl()
{
}



//=================================================================================================
//  MfvMuscl::ComputeGodunovFlux
/// Calculate the Godunov flux between particle i and all neighbours storing all the partial
/// sums of conserved variables, dQ, between neighbours.
//=================================================================================================
template <int ndim, template<int> class kernelclass>
void MfvMuscl<ndim, kernelclass>::ComputeGodunovFlux
 (const int i,                         ///< [in] id of particle
  const int Nneib,                     ///< [in] No. of neins in neibpart array
  const FLOAT timestep,                ///< [in] Minimum timestep size
  int *neiblist,                       ///< [in] id of gather neibs in neibpart
  FLOAT *drmag,                        ///< [in] Distances of gather neighbours
  FLOAT *invdrmag,                     ///< [in] Inverse distances of gather neibs
  FLOAT *dr,                           ///< [in] Position vector of gather neibs
  MeshlessFVParticle<ndim> &part,      ///< [inout] Particle i data
  MeshlessFVParticle<ndim> *neibpart)  ///< [inout] Neighbour particle data
{
  int j;                               // Neighbour list id
  int jj;                              // Aux. neighbour counter
  int k;                               // Dimension counter
  int var;                             // Particle state vector variable counter
  FLOAT Aij[ndim];                     // Pseudo 'Area' vector
  FLOAT draux[ndim];                   // Position vector of part relative to neighbour
  FLOAT dr_unit[ndim];                 // Unit vector from neighbour to part
  FLOAT drsqd;                         // Distance squared
  FLOAT invdrmagaux;                   // 1 / distance
  FLOAT psitildai[ndim];               // Normalised gradient psi value for particle i
  FLOAT psitildaj[ndim];               // Normalised gradient psi value for neighbour j
  FLOAT rface[ndim];                   // Position of working face (to compute Godunov fluxes)
  FLOAT vface[ndim];                   // Velocity of working face (to compute Godunov fluxes)
  FLOAT flux[nvar][ndim];              // Flux tensor
  FLOAT Wleft[nvar];                   // Primitive vector for LHS of Riemann problem
  FLOAT Wright[nvar];                  // Primitive vector for RHS of Riemann problem
  FLOAT Wdot[nvar];                    // Time derivative of primitive vector
  FLOAT gradW[nvar][ndim];             // Gradient of primitive vector
  FLOAT dW[nvar];                      // Change in primitive quantities
  const FLOAT dt = timestep*(FLOAT) part.nstep;    // Timestep of given particle


  // Loop over all potential neighbours in the list
  //-----------------------------------------------------------------------------------------------
  for (jj=0; jj<Nneib; jj++) {
    j = neiblist[jj];

    for (k=0; k<ndim; k++) draux[k] = part.r[k] - neibpart[j].r[k];
    drsqd = DotProduct(draux, draux, ndim);
    invdrmagaux = (FLOAT) 1.0/sqrt(drsqd + small_number);
    for (k=0; k<ndim; k++) dr_unit[k] = draux[k]*invdrmagaux;

    // Calculate psitilda values
    for (k=0; k<ndim; k++) {
      psitildai[k] = (FLOAT) 0.0;
      psitildaj[k] = (FLOAT) 0.0;
      for (int kk=0; kk<ndim; kk++) {
        psitildai[k] += neibpart[j].B[k][kk]*draux[kk]*neibpart[j].hfactor*
          kern.w0_s2(drsqd*neibpart[j].invh*neibpart[j].invh)/neibpart[j].ndens;
        psitildaj[k] -= part.B[k][kk]*draux[kk]*part.hfactor*
          kern.w0_s2(drsqd*part.invh*part.invh)/part.ndens;
      }
      Aij[k] = part.volume*psitildaj[k] - neibpart[j].volume*psitildai[k];
    }

    // Calculate position and velocity of the face
    if (staticParticles) {
      for (k=0; k<ndim; k++) rface[k] = (FLOAT) 0.5*(part.r[k] + neibpart[j].r[k]);
      for (k=0; k<ndim; k++) vface[k] = (FLOAT) 0.0;
    }
    else {
      for (k=0; k<ndim; k++) rface[k] = (FLOAT) 0.5*(part.r[k] + neibpart[j].r[k]);
      //for (k=0; k<ndim; k++) rface[k] = part.r[k] +
      //  part.h*(neibpart[j].r[k] - part.r[k])/(part.h + neibpart[j].h);
      for (k=0; k<ndim; k++) draux[k] = part.r[k] - rface[k];
      for (k=0; k<ndim; k++) vface[k] = part.v[k] +
        (neibpart[j].v[k] - part.v[k])*DotProduct(draux, dr_unit, ndim)*invdrmagaux;
    }

    // Compute slope-limited values for LHS
    for (k=0; k<ndim; k++) draux[k] = rface[k] - part.r[k];
    limiter->ComputeLimitedSlopes(part, neibpart[j], draux, gradW, dW);
    for (var=0; var<nvar; var++) Wleft[var] = part.Wprim[var] + dW[var];
    for (k=0; k<ndim; k++) Wleft[k] -= vface[k];

    // Time-integrate LHS state to half-timestep value
    this->CalculatePrimitiveTimeDerivative(Wleft, gradW, Wdot);
    for (k=0; k<ndim; k++) Wdot[k] += part.a[k];
    for (var=0; var<nvar; var++) Wleft[var] -= (FLOAT) 0.5*Wdot[var]*dt;

    // Compute slope-limited values for RHS
    for (k=0; k<ndim; k++) draux[k] = rface[k] - neibpart[j].r[k];
    limiter->ComputeLimitedSlopes(neibpart[j], part, draux, gradW, dW);
    for (var=0; var<nvar; var++) Wright[var] = neibpart[j].Wprim[var] + dW[var];
    for (k=0; k<ndim; k++) Wright[k] -= vface[k];

    // Time-integrate RHS state to half-timestep value
    this->CalculatePrimitiveTimeDerivative(Wright, gradW, Wdot);
    for (k=0; k<ndim; k++) Wdot[k] += neibpart[j].a[k];
    for (var=0; var<nvar; var++) Wright[var] -= (FLOAT) 0.5*Wdot[var]*dt;

    assert(Wleft[irho] > 0.0);
    assert(Wleft[ipress] > 0.0);
    assert(Wright[irho] > 0.0);
    assert(Wright[ipress] > 0.0);

    // Calculate Godunov flux using the selected Riemann solver
    riemann->ComputeFluxes(Wright, Wleft, dr_unit, vface, flux);

    // Finally calculate flux terms for all quantities based on Lanson & Vila gradient operators
    for (var=0; var<nvar; var++) {
      part.dQ[var] -= DotProduct(flux[var], Aij, ndim)*dt;
      neibpart[j].dQ[var] += DotProduct(flux[var], Aij, ndim)*dt;
    }

    // Compute mass-loss moments for gravitational correction terms
    for (k=0; k<ndim; k++) {
      part.rdmdt[k] -= (part.r[k] - neibpart[j].r[k])*DotProduct(flux[var], Aij, ndim);
      neibpart[j].rdmdt[k] += (part.r[k] - neibpart[j].r[k])*DotProduct(flux[var], Aij, ndim);
    }

  }
  //-----------------------------------------------------------------------------------------------


  return;
}



template class MfvMuscl<1, M4Kernel>;
template class MfvMuscl<2, M4Kernel>;
template class MfvMuscl<3, M4Kernel>;
template class MfvMuscl<1, QuinticKernel>;
template class MfvMuscl<2, QuinticKernel>;
template class MfvMuscl<3, QuinticKernel>;
template class MfvMuscl<1, TabulatedKernel>;
template class MfvMuscl<2, TabulatedKernel>;
template class MfvMuscl<3, TabulatedKernel>;
