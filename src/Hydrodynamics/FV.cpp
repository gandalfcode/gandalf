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
                      _gas_eos, KernelName, size_part, units, params)
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
//  FV::CalculatePrimitiveTimeDerivative
/// Calculate the time derivative of the primitive variables, dWprim/dt.  Used for extrapolating
/// the primitive variables forward in time for, e.g. the half-step for the MUSCL scheme.
//=================================================================================================
template <int ndim>
void FV<ndim>::CalculatePrimitiveTimeDerivative
 (const FLOAT Wprim[nvar],             ///< [in] Primitive vector, Wprim, of particle
  const FLOAT gradW[nvar][ndim],       ///< [in] Gradients of primitive quantities
  const FLOAT c_s,                     ///< [in] Sound speed
  FLOAT Wdot[nvar])                    ///< [out] Time derivatives of primitive quantities
{
  double divV = 0;
  for (int i=0; i < ndim; i++) divV += gradW[i][i] ;

  Wdot[irho]   = - DotProduct(Wprim, gradW[irho], ndim) - Wprim[irho]*divV;
  Wdot[ipress] = - DotProduct(Wprim, gradW[ipress], ndim) - Wprim[irho]*c_s*c_s*divV;
  for (int i=0; i < ndim; i++)
    Wdot[i]    = - DotProduct(Wprim, gradW[i], ndim) - gradW[ipress][i]/Wprim[irho];

  return;
}


template class FV<1>;
template class FV<2>;
template class FV<3>;
