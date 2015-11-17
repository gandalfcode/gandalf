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



//=================================================================================================
//  FV::FV
/// FV class constructor.  Calls main Hydrodynamics class constructor and also
/// sets additional kernel-related quantities.
//=================================================================================================
template <int ndim>
FV<ndim>::FV(int _hydro_forces, int _self_gravity, FLOAT _accel_mult, FLOAT _courant_mult,
             FLOAT _h_fac, FLOAT h_converge_aux, FLOAT gamma_aux, string _gas_eos,
             string KernelName, int size_part, SimUnits &units, Parameters *params):
  Hydrodynamics<ndim>(_hydro_forces, _self_gravity, _h_fac,
                      _gas_eos, KernelName, size_part, units, params),
  gamma_eos(gamma_aux),
  gammam1(gamma_aux - 1.0)
{
}



//=================================================================================================
//  FV::~FV
/// FV class destructor
//=================================================================================================
template <int ndim>
FV<ndim>::~FV()
{
  //DeallocateMemory();
}



//=================================================================================================
//  FV::ConvertConservedToQ
/// ...
//=================================================================================================
template <int ndim>
void FV<ndim>::ConvertConservedToQ
 (const FLOAT volume,
  const FLOAT Ucons[nvar],
  FLOAT Qcons[nvar])
{
  for (int var=0; var<nvar; var++) Qcons[var] = volume*Ucons[var];

  return;
}



//=================================================================================================
//  FV::ConvertQToConserved
/// ...
//=================================================================================================
template <int ndim>
void FV<ndim>::ConvertQToConserved
 (const FLOAT volume,
  const FLOAT Qcons[nvar],
  FLOAT Ucons[nvar])
{
  for (int var=0; var<nvar; var++) Ucons[var] = Qcons[var]/volume;

  return;
}




//=================================================================================================
//  FV::ConvertConservedToPrimitive
/// ...
//=================================================================================================
template <int ndim>
void FV<ndim>::ConvertConservedToPrimitive
 (const FLOAT Ucons[nvar],
  FLOAT Wprim[nvar])
{
  int k;
  FLOAT ekin = 0.0;

  Wprim[irho] = Ucons[irho];
  for (k=0; k<ndim; k++) {
    Wprim[k] = Ucons[k]/Ucons[irho];
    ekin += Wprim[k]*Wprim[k];
  }
  Wprim[ipress] = (gamma_eos - 1.0)*(Ucons[ietot] - 0.5*Ucons[irho]*ekin);
  //Wprim[irho] = Ucons[irho]

  return;
}



//=================================================================================================
//  FV::ConvertPrimitiveToConserved
/// ...
//=================================================================================================
template <int ndim>
void FV<ndim>::ConvertPrimitiveToConserved
 (const FLOAT Wprim[nvar],
  FLOAT Ucons[nvar])
{
  int k;
  FLOAT ekin = 0.0;

  Ucons[irho] = Wprim[irho];
  for (k=0; k<ndim; k++) {
    Ucons[k] = Wprim[k]*Wprim[irho];
    ekin += Wprim[k]*Wprim[k];
  }
  Ucons[ietot] = Wprim[ipress]/(gamma_eos - 1.0) + 0.5*Wprim[irho]*ekin;

  return;
}



//=================================================================================================
//  FV::CalculateConservedFluxFromConserved
/// ...
//=================================================================================================
/*template <int ndim>
void FV<ndim>::CalculateConservedFluxFromConserved
 (int k,
  FLOAT Ucons[nvar],
  FLOAT fluxVector[nvar])
{
  int kv;
  FLOAT ekin = 0.0;
  FLOAT press;

  for (kv=0; kv<ndim; kv++) {
    ekin += Ucons[kv]*Ucons[kv];
  }
  press = (gamma_eos - 1.0)*(Ucons[ietot] - 0.5*ekin/Ucons[irho]);

  // Calculate fluxes from conserved variables (NOT Riemann solver) to compute evolved boundary values
  for (kv=0; kv<ndim; kv++) fluxVector[kv] = Ucons[k]*Ucons[kv]/Ucons[irho];
  fluxVector[k]     = Ucons[k]*Ucons[k]/Ucons[irho] + press;
  fluxVector[irho]  = Ucons[k];
  fluxVector[ietot] = Ucons[k]*(Ucons[ietot] + press)/Ucons[irho];

  return;
}*/



//=================================================================================================
//  FV::CalculateFluxVector
/// ...
//=================================================================================================
template <int ndim>
void FV<ndim>::CalculateFluxVectorFromPrimitive
 (FLOAT Wprim[nvar],
  FLOAT fluxVector[nvar][ndim])
{
  int k;
  int kv;
  FLOAT ekin = 0.0;

  for (kv=0; kv<ndim; kv++) {
    ekin += Wprim[kv]*Wprim[kv];
  }

  for (k=0; k<ndim; k++) {
    for (kv=0; kv<ndim; kv++) fluxVector[kv][k] = Wprim[irho]*Wprim[k]*Wprim[kv];
    fluxVector[k][k]     = Wprim[irho]*Wprim[k]*Wprim[k] + Wprim[ipress];
    fluxVector[irho][k]  = Wprim[irho]*Wprim[k];
    fluxVector[ietot][k] =
      Wprim[k]*(Wprim[ipress]/(gamma_eos - 1.0) + 0.5*Wprim[irho]*ekin + Wprim[ipress]);
  }

  return;
}



//=================================================================================================
//  FV::CalculatePrimitiveTimeDerivative
/// ...
//=================================================================================================
template <int ndim>
void FV<ndim>::CalculatePrimitiveTimeDerivative
 (FLOAT Wprim[nvar],
  FLOAT gradW[nvar][ndim],
  FLOAT Wdot[nvar])
{
  if (ndim == 1) {
    Wdot[irho]   = Wprim[ivx]*gradW[irho][0] + Wprim[irho]*gradW[ivx][0];
    Wdot[ivx]    = Wprim[ivx]*gradW[ivx][0] + gradW[ipress][0]/Wprim[irho];
    Wdot[ipress] = gamma_eos*Wprim[ipress]*gradW[ivx][0] + Wprim[ivx]*gradW[ipress][0];
  }
  else if (ndim == 2) {
    Wdot[irho]   = Wprim[ivx]*gradW[irho][0] + Wprim[ivy]*gradW[irho][1] +
      Wprim[irho]*(gradW[ivx][0] + gradW[ivy][1]);
    Wdot[ivx]    = Wprim[ivx]*gradW[ivx][0] + Wprim[ivy]*gradW[ivx][1] + gradW[ipress][0]/Wprim[irho];
    Wdot[ivy]    = Wprim[ivx]*gradW[ivy][0] + Wprim[ivy]*gradW[ivy][1] + gradW[ipress][1]/Wprim[irho];
    Wdot[ipress] = Wprim[ivx]*gradW[ipress][0] + Wprim[ivy]*gradW[ipress][1] +
      gamma_eos*Wprim[ipress]*(gradW[ivx][0] + gradW[ivy][1]);
  }
  else if (ndim == 3) {
    Wdot[irho]   = Wprim[ivx]*gradW[irho][0] + Wprim[ivy]*gradW[irho][1] +
      Wprim[ivz]*gradW[irho][2] + Wprim[irho]*(gradW[ivx][0] + gradW[ivy][1] + gradW[ivz][2]);
    Wdot[ivx]    = Wprim[ivx]*gradW[ivx][0] + Wprim[ivy]*gradW[ivx][1] +
      Wprim[ivz]*gradW[ivx][2] + gradW[ipress][0]/Wprim[irho];
    Wdot[ivy]    = Wprim[ivx]*gradW[ivy][0] + Wprim[ivy]*gradW[ivy][1] +
      Wprim[ivz]*gradW[ivy][2] + gradW[ipress][1]/Wprim[irho];
    Wdot[ivz]    = Wprim[ivx]*gradW[ivz][0] + Wprim[ivy]*gradW[ivz][1] +
      Wprim[ivz]*gradW[ivz][2] + gradW[ipress][2]/Wprim[irho];
    Wdot[ipress] = Wprim[ivx]*gradW[ipress][0] + Wprim[ivy]*gradW[ipress][1] +
      Wprim[ivz]*gradW[ipress][2] +
      gamma_eos*Wprim[ipress]*(gradW[ivx][0] + gradW[ivy][1] + gradW[ivz][2]);
  }

  return;
}



template class FV<1>;
template class FV<2>;
template class FV<3>;
