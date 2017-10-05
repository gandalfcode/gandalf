//=================================================================================================
//  MfvRungeKutta.cpp
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
//  MfvRungeKutta::MfvRungeKutta
/// MfvRungeKutta class constructor.  Calls main SPH class constructor and also
/// sets additional kernel-related quantities
//=================================================================================================
template <int ndim, template<int> class kernelclass, class SlopeLimiter>
MfvRungeKutta<ndim, kernelclass,SlopeLimiter>::MfvRungeKutta
 (int _hydro_forces, int _self_gravity, FLOAT _accel_mult, FLOAT _courant_mult,
  FLOAT _h_fac, FLOAT _h_converge, FLOAT _gamma, string _gas_eos, string KernelName,
  int size_part, SimUnits &units, Parameters *params):
  MfvCommon<ndim,kernelclass,SlopeLimiter>(_hydro_forces, _self_gravity, _accel_mult, _courant_mult, _h_fac,
                              _h_converge, _gamma, _gas_eos, KernelName, size_part, units, params)
{
}



//=================================================================================================
//  MfvRungeKutta::ComputeGodunovFlux
/// ...
//=================================================================================================
template <int ndim, template<int> class kernelclass, class SlopeLimiter>
void MfvRungeKutta<ndim, kernelclass,SlopeLimiter>::ComputeGodunovFlux
(MeshlessFVParticle<ndim>& part,                                ///< [inout] Particle data
 NeighbourList<typename MeshlessFV<ndim>::FluxNeib>& neibpart,  ///< [inout] Neighbour data
 FLOAT timestep)                                                ///< [in]    Current timstep size
{
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

  FLOAT invh_i   = 1/part.h;
  FLOAT volume_i = 1/part.ndens;

  // Loop over all potential neighbours in the list
  //-----------------------------------------------------------------------------------------------
  int Nneib = neibpart.size() ;
  for (int j=0; j<Nneib; j++) {

    FLOAT invh_j   = 1/neibpart[j].h;
    FLOAT volume_j = 1/neibpart[j].ndens;

    for (k=0; k<ndim; k++) draux[k] = part.r[k] - neibpart[j].r[k];
    drsqd = DotProduct(draux, draux, ndim);
    invdrmagaux = 1.0/sqrt(drsqd + small_number);
    for (k=0; k<ndim; k++) dr_unit[k] = draux[k]*invdrmagaux;

    // Compute psi-tilda values using integral / sph gradients.
    if (not part.flags.check(bad_gradients)) {
      for (k=0; k<ndim; k++) {
        psitildaj[k] = 0;
        for (int kk=0; kk<ndim; kk++)
          psitildaj[k] += part.B[k][kk]*draux[kk]*part.hfactor*
              kern.w0_s2(drsqd*invh_i*invh_i)*volume_i;
      }
    }
    else {
      double dr = sqrt(drsqd) + small_number ;
      double w = part.hfactor*volume_i * kern.w1(dr*invh_i);
      for (k=0; k<ndim; k++)  psitildaj[k] = - (draux[k]/dr) * w;
    }

    if (not neibpart[j].flags.check(bad_gradients)) {
      for (k=0; k<ndim; k++) {
      psitildai[k] = 0;
      for (int kk=0; kk<ndim; kk++)
        psitildai[k] -= neibpart[j].B[k][kk]*draux[kk]*neibpart[j].hfactor*
          kern.w0_s2(drsqd*invh_j*invh_j)*volume_j;
      }
    }
    else {
      double dr = sqrt(drsqd) + small_number ;
      double w = neibpart[j].hfactor*volume_j * kern.w1(dr*invh_j);
      for (k=0; k<ndim; k++) psitildai[k] = + (draux[k]/dr) * w;
    }

    // Compute the face area
    for (k=0; k<ndim; k++)
      Aij[k] = volume_i*psitildaj[k] - volume_j*psitildai[k];

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
    limiter.ComputeLimitedSlopes(part, neibpart[j], draux, gradW, dW);
    for (var=0; var<nvar; var++) Wleft[var] = part.Wprim[var] + dW[var];
    for (k=0; k<ndim; k++) Wleft[k] -= vface[k];

    // Compute slope-limited values for RHS
    for (k=0; k<ndim; k++) draux[k] = rface[k] - neibpart[j].r[k];
    limiter.ComputeLimitedSlopes(neibpart[j], part, draux, gradW, dW);
    for (var=0; var<nvar; var++) Wright[var] = neibpart[j].Wprim[var] + dW[var];
    for (k=0; k<ndim; k++) Wright[k] -= vface[k];

    assert(Wleft[irho] > 0.0);
    assert(Wleft[ipress] > 0.0);
    assert(Wright[irho] > 0.0);
    assert(Wright[ipress] > 0.0);

    // Calculate Godunov flux using the selected Riemann solver
    if (RiemannSolverType == exact)
      riemannExact.ComputeFluxes(Wright, Wleft, dr_unit, vface, flux);
    else
      riemannHLLC.ComputeFluxes(Wright, Wleft, dr_unit, vface, flux);

    // Finally calculate flux terms for all quantities based on Lanson & Vila gradient operators
    for (var=0; var<nvar; var++) {
      part.dQ[var] -= 0.5*DotProduct(flux[var], Aij, ndim)*dt;
      neibpart[j].dQ[var] += 0.5*DotProduct(flux[var], Aij, ndim)*dt;
    }

    // Time-integrate LHS state to half-timestep value
    this->CalculatePrimitiveTimeDerivative(Wleft, gradW, Wdot);
    for (k=0; k<ndim; k++) Wdot[k] += part.a[k];
    for (var=0; var<nvar; var++) Wleft[var] += (FLOAT) Wdot[var]*dt;

    // Time-integrate RHS state to half-timestep value
    this->CalculatePrimitiveTimeDerivative(Wright, gradW, Wdot);
    for (k=0; k<ndim; k++) Wdot[k] += neibpart[j].a[k];
    for (var=0; var<nvar; var++) Wright[var] += (FLOAT) Wdot[var]*dt;

    assert(Wleft[irho] > 0.0);
    assert(Wleft[ipress] > 0.0);
    assert(Wright[irho] > 0.0);
    assert(Wright[ipress] > 0.0);

    // Calculate Godunov flux using the selected Riemann solver
    if (RiemannSolverType == exact)
      riemannExact.ComputeFluxes(Wright, Wleft, dr_unit, vface, flux);
    else
      riemannHLLC.ComputeFluxes(Wright, Wleft, dr_unit, vface, flux);

    // Finally calculate flux terms for all quantities based on Lanson & Vila gradient operators
    for (var=0; var<nvar; var++) {
      part.dQ[var] -= 0.5*DotProduct(flux[var], Aij, ndim)*dt;
      neibpart[j].dQ[var] += 0.5*DotProduct(flux[var], Aij, ndim)*dt;
    }

  }
  //-----------------------------------------------------------------------------------------------


  return;
}




template class MfvRungeKutta<1, M4Kernel, NullLimiter<1> >;
template class MfvRungeKutta<2, M4Kernel, NullLimiter<2> >;
template class MfvRungeKutta<3, M4Kernel, NullLimiter<3> >;
template class MfvRungeKutta<1, QuinticKernel, NullLimiter<1> >;
template class MfvRungeKutta<2, QuinticKernel, NullLimiter<2> >;
template class MfvRungeKutta<3, QuinticKernel, NullLimiter<3> >;
template class MfvRungeKutta<1, TabulatedKernel, NullLimiter<1> >;
template class MfvRungeKutta<2, TabulatedKernel, NullLimiter<2> >;
template class MfvRungeKutta<3, TabulatedKernel, NullLimiter<3> >;


template class MfvRungeKutta<1, M4Kernel, ZeroSlopeLimiter<1> >;
template class MfvRungeKutta<2, M4Kernel, ZeroSlopeLimiter<2> >;
template class MfvRungeKutta<3, M4Kernel, ZeroSlopeLimiter<3> >;
template class MfvRungeKutta<1, QuinticKernel, ZeroSlopeLimiter<1> >;
template class MfvRungeKutta<2, QuinticKernel, ZeroSlopeLimiter<2> >;
template class MfvRungeKutta<3, QuinticKernel, ZeroSlopeLimiter<3> >;
template class MfvRungeKutta<1, TabulatedKernel, ZeroSlopeLimiter<1> >;
template class MfvRungeKutta<2, TabulatedKernel, ZeroSlopeLimiter<2> >;
template class MfvRungeKutta<3, TabulatedKernel, ZeroSlopeLimiter<3> >;

template class MfvRungeKutta<1, M4Kernel, TVDScalarLimiter<1> >;
template class MfvRungeKutta<2, M4Kernel, TVDScalarLimiter<2> >;
template class MfvRungeKutta<3, M4Kernel, TVDScalarLimiter<3> >;
template class MfvRungeKutta<1, QuinticKernel, TVDScalarLimiter<1> >;
template class MfvRungeKutta<2, QuinticKernel, TVDScalarLimiter<2> >;
template class MfvRungeKutta<3, QuinticKernel, TVDScalarLimiter<3> >;
template class MfvRungeKutta<1, TabulatedKernel,TVDScalarLimiter<1> >;
template class MfvRungeKutta<2, TabulatedKernel, TVDScalarLimiter<2> >;
template class MfvRungeKutta<3, TabulatedKernel, TVDScalarLimiter<3> >;

template class MfvRungeKutta<1, M4Kernel, ScalarLimiter<1> >;
template class MfvRungeKutta<2, M4Kernel, ScalarLimiter<2> >;
template class MfvRungeKutta<3, M4Kernel, ScalarLimiter<3> >;
template class MfvRungeKutta<1, QuinticKernel, ScalarLimiter<1> >;
template class MfvRungeKutta<2, QuinticKernel, ScalarLimiter<2> >;
template class MfvRungeKutta<3, QuinticKernel, ScalarLimiter<3> >;
template class MfvRungeKutta<1, TabulatedKernel, ScalarLimiter<1> >;
template class MfvRungeKutta<2, TabulatedKernel, ScalarLimiter<2> >;
template class MfvRungeKutta<3, TabulatedKernel, ScalarLimiter<3> >;

template class MfvRungeKutta<1, M4Kernel, Springel2009Limiter<1> >;
template class MfvRungeKutta<2, M4Kernel, Springel2009Limiter<2> >;
template class MfvRungeKutta<3, M4Kernel, Springel2009Limiter<3> >;
template class MfvRungeKutta<1, QuinticKernel, Springel2009Limiter<1> >;
template class MfvRungeKutta<2, QuinticKernel, Springel2009Limiter<2> >;
template class MfvRungeKutta<3, QuinticKernel, Springel2009Limiter<3> >;
template class MfvRungeKutta<1, TabulatedKernel, Springel2009Limiter<1> >;
template class MfvRungeKutta<2, TabulatedKernel, Springel2009Limiter<2> >;
template class MfvRungeKutta<3, TabulatedKernel, Springel2009Limiter<3> >;

template class MfvRungeKutta<1, M4Kernel, GizmoLimiter<1> >;
template class MfvRungeKutta<2, M4Kernel, GizmoLimiter<2> >;
template class MfvRungeKutta<3, M4Kernel, GizmoLimiter<3> >;
template class MfvRungeKutta<1, QuinticKernel, GizmoLimiter<1> >;
template class MfvRungeKutta<2, QuinticKernel, GizmoLimiter<2> >;
template class MfvRungeKutta<3, QuinticKernel, GizmoLimiter<3> >;
template class MfvRungeKutta<1, TabulatedKernel, GizmoLimiter<1> >;
template class MfvRungeKutta<2, TabulatedKernel, GizmoLimiter<2> >;
template class MfvRungeKutta<3, TabulatedKernel, GizmoLimiter<3> >;

