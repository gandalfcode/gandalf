//=================================================================================================
//  FV.cpp
//  Virtual base class containing all common functionality for all Finite-Volume Hydrodynamics
//  schemes in GANDALF (e.g. Meshless-FV, MovingMeshFV).
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
#include "Debug.h"
#include "EOS.h"
#include "Exception.h"
#include "FV.h"
#include "InlineFuncs.h"
#include "Particle.h"
#include "Parameters.h"
#include "Precision.h"
#include "SmoothingKernel.h"
using namespace std;


// Declare invndim constant here (prevents warnings with some compilers)
template <int ndim>
const FLOAT FV<ndim>::invndim = 1.0/ndim;


//=================================================================================================
//  FV::FV
/// FV class constructor.  Calls main Hydrodynamics class constructor and also
/// sets additional kernel-related quantities.
//=================================================================================================
template <int ndim>
FV<ndim>::FV(int _hydro_forces, int _self_gravity, FLOAT _accel_mult, FLOAT _courant_mult,
             FLOAT _h_fac, FLOAT _h_converge, FLOAT _gamma, string _gas_eos,
             string KernelName, int size_part, SimUnits &units, Parameters *params):
  Hydrodynamics<ndim>(_hydro_forces, _self_gravity, _h_fac,
                      _gas_eos, KernelName, size_part, units, params),
  gamma_eos(_gamma),
  gammam1(_gamma - 1.0),
  eta_eos(eos->eta)
{
}



//=================================================================================================
//  FV::~FV
/// FV class destructor
//=================================================================================================
template <int ndim>
FV<ndim>::~FV()
{ }



//=================================================================================================
//  FV::ConvertConservedToPrimitive
/// Convert conserved vector, Qcons (NOT Ucons) to the primitive vector, Wprim.
//=================================================================================================
template <int ndim>
void FV<ndim>::ConvertConservedToPrimitive
 (const FLOAT ndens,                   ///< [in] Inverse Effective volume of particle
  const FLOAT Qcons[nvar],             ///< [in] Conserved vector, Qcons, of particle
  FLOAT Wprim[nvar])                   ///< [out] Primitive vector, Wprim, of particle
{
  int k;
  FLOAT ekin = 0.0;

  Wprim[irho] = Qcons[irho]*ndens;
  for (k=0; k<ndim; k++) {
    Wprim[k] = Qcons[k]/Qcons[irho];
    ekin += Wprim[k]*Wprim[k];
  }
  Wprim[ipress] = (gamma_eos - 1.0)*(Qcons[ietot] - 0.5*Qcons[irho]*ekin)*ndens;

  return;
}



//=================================================================================================
//  FV::ConvertPrimitiveToConserved
/// Convert from primitive vector, Wprim, to conserved vector, Qcons (NOT Ucons).
//=================================================================================================
template <int ndim>
void FV<ndim>::ConvertPrimitiveToConserved
 (const FLOAT ndens,                   ///< [in] Inverse Effective volume of particle
  const FLOAT Wprim[nvar],             ///< [in] Primitive vector, Wprim, of particle
  FLOAT Qcons[nvar])                   ///< [out] Conserved vector, Qcons, of particle
{
  FLOAT ekin = (FLOAT) 0.0;

  Qcons[irho] = Wprim[irho]/ndens;
  for (int k=0; k<ndim; k++) {
    Qcons[k] = Wprim[k]*Wprim[irho]/ndens;
    ekin += Wprim[k]*Wprim[k];
  }
  Qcons[ietot] = (Wprim[ipress]/(gamma_eos - 1.0) + (FLOAT) 0.5*Wprim[irho]*ekin)/ndens;

  return;
}



//=================================================================================================
//  FV::CalculateFluxVectorFromPrimitive
/// Calculate the flux vector of conserved variables, U, from the primitive variables, Wprim.
//=================================================================================================
template <int ndim>
void FV<ndim>::CalculateFluxVectorFromPrimitive
 (const FLOAT Wprim[nvar],             ///< [in] Primitive vector, Wprim, of particle
  FLOAT fluxVector[nvar][ndim])        ///< [out] Flux vector
{
  FLOAT ekin = (FLOAT) 0.0;

  for (int kv=0; kv<ndim; kv++) ekin += Wprim[kv]*Wprim[kv];

  for (int k=0; k<ndim; k++) {
    for (int kv=0; kv<ndim; kv++) fluxVector[kv][k] = Wprim[irho]*Wprim[k]*Wprim[kv];
    fluxVector[k][k]     = Wprim[irho]*Wprim[k]*Wprim[k] + Wprim[ipress];
    fluxVector[irho][k]  = Wprim[irho]*Wprim[k];
    fluxVector[ietot][k] = Wprim[k]*(Wprim[ipress]/(gamma_eos - (FLOAT) 1.0) +
      (FLOAT) 0.5*Wprim[irho]*ekin + Wprim[ipress]);
  }

  return;
}



//=================================================================================================
//  FV::CalculatePrimitiveTimeDerivative
/// Calculate the time derivative of the primitive variables, dWprim/dt.  Used for extrapolating
/// the primitive variables forward in time for, e.g. the half-step for the MUSCL scheme.
//=================================================================================================
template <>
void FV<1>::CalculatePrimitiveTimeDerivative
 (const FLOAT Wprim[nvar],             ///< [in] Primitive vector, Wprim, of particle
  const FLOAT gradW[nvar][1],          ///< [in] Gradients of primitive quantities
  FLOAT Wdot[nvar])                    ///< [out] Time derivatives of primitive quantities
{
  Wdot[irho]   = -Wprim[ivx]*gradW[irho][0] - Wprim[irho]*gradW[ivx][0];
  Wdot[ivx]    = -Wprim[ivx]*gradW[ivx][0] - gradW[ipress][0]/Wprim[irho];
  //Wdot[ipress] = -Wprim[ipress]*gradW[ivx][0] - Wprim[ivx]*gradW[ipress][0];
  Wdot[ipress] = -eta_eos*Wprim[ipress]*gradW[ivx][0] - Wprim[ivx]*gradW[ipress][0];

  return;
}

template <>
void FV<2>::CalculatePrimitiveTimeDerivative
 (const FLOAT Wprim[nvar],             ///< [in] Primitive vector, Wprim, of particle
  const FLOAT gradW[nvar][2],          ///< [in] Gradients of primitive quantities
  FLOAT Wdot[nvar])                    ///< [out] Time derivatives of primitive quantities
{
  Wdot[irho]   = -Wprim[ivx]*gradW[irho][0] - Wprim[ivy]*gradW[irho][1] -
    Wprim[irho]*(gradW[ivx][0] + gradW[ivy][1]);
  Wdot[ivx]    = -Wprim[ivx]*gradW[ivx][0] - Wprim[ivy]*gradW[ivx][1] - gradW[ipress][0]/Wprim[irho];
  Wdot[ivy]    = -Wprim[ivx]*gradW[ivy][0] - Wprim[ivy]*gradW[ivy][1] - gradW[ipress][1]/Wprim[irho];
  Wdot[ipress] = -Wprim[ivx]*gradW[ipress][0] - Wprim[ivy]*gradW[ipress][1] -
    //Wprim[ipress]*(gradW[ivx][0] + gradW[ivy][1]);
    eta_eos*Wprim[ipress]*(gradW[ivx][0] + gradW[ivy][1]);
  return;
}

template <>
void FV<3>::CalculatePrimitiveTimeDerivative
 (const FLOAT Wprim[nvar],             ///< [in] Primitive vector, Wprim, of particle
  const FLOAT gradW[nvar][3],          ///< [in] Gradients of primitive quantities
  FLOAT Wdot[nvar])                    ///< [out] Time derivatives of primitive quantities
{
  Wdot[irho]   = -Wprim[ivx]*gradW[irho][0] - Wprim[ivy]*gradW[irho][1] -
   Wprim[ivz]*gradW[irho][2] - Wprim[irho]*(gradW[ivx][0] + gradW[ivy][1] + gradW[ivz][2]);
  Wdot[ivx]    = -Wprim[ivx]*gradW[ivx][0] - Wprim[ivy]*gradW[ivx][1] -
    Wprim[ivz]*gradW[ivx][2] - gradW[ipress][0]/Wprim[irho];
  Wdot[ivy]    = -Wprim[ivx]*gradW[ivy][0] - Wprim[ivy]*gradW[ivy][1] -
    Wprim[ivz]*gradW[ivy][2] - gradW[ipress][1]/Wprim[irho];
  Wdot[ivz]    = -Wprim[ivx]*gradW[ivz][0] - Wprim[ivy]*gradW[ivz][1] -
    Wprim[ivz]*gradW[ivz][2] - gradW[ipress][2]/Wprim[irho];
  Wdot[ipress] = -Wprim[ivx]*gradW[ipress][0] - Wprim[ivy]*gradW[ipress][1] -
    Wprim[ivz]*gradW[ipress][2] -
    //Wprim[ipress]*(gradW[ivx][0] + gradW[ivy][1] + gradW[ivz][2]);
    eta_eos*Wprim[ipress]*(gradW[ivx][0] + gradW[ivy][1] + gradW[ivz][2]);
  return;
}



template class FV<1>;
template class FV<2>;
template class FV<3>;
