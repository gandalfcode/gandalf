//=================================================================================================
//  EOS.cpp
//  Contains function definitions for the base Equation of state.
//
//  This file is part of GANDALF :
//  Graphical Astrophysics code for N-body Dynamics And Lagrangian Fluids
//  https://github.com/gandalfcode/gandalf
//  Contact : gandalfcode@gmail.com
//
//  Copyright (C) 2013  D. A. Hubber, G. Rosotti, R. Booth
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

#include "EOS.h"

//=================================================================================================
//  EOS::PrimitiveToConserved()
/// Convert primitive variables to conserved quantities
//=================================================================================================
template<int ndim>
void EOS<ndim>::PrimitiveToConserved
(const PrimitiveVariables<ndim>& Wprim,       ///< [in]  Primitive variables
 ConservedVariables<ndim>& Ucons)             ///< [out] Conserved variables
{
  Ucons.density = Wprim.density;

  FLOAT Ek = 0;
  for (int i=0; i<ndim; ++i) {
    Ucons.momentum[i] = Wprim.density * Wprim.velocity[i];
    Ek += 0.5 *  Wprim.density * Wprim.velocity[i] * Wprim.velocity[i];
  }

  Ucons.energy = Ek + Wprim.density * this->InternalEnergy(Wprim.density, Wprim.pressure);
}


//=================================================================================================
//  EOS::ConservedToPrimitive()
/// Convert conserved quantities to primitive variables
//=================================================================================================
template<int ndim>
void EOS<ndim>::ConservedToPrimitive
(const ConservedVariables<ndim>& Ucons,       ///< [in]  Conserved variables
 PrimitiveVariables<ndim>& Wprim)             ///< [out] Primitive variables
{
  Wprim.density = Ucons.density;

  FLOAT Ek = 0;
  for (int i=0; i<ndim; ++i) {
    Wprim.velocity[i] = Ucons.momentum[i] / Ucons.density ;
    Ek += 0.5 *   Ucons.momentum[i] * Ucons.momentum[i] / Ucons.density;
  }

  Wprim.pressure =  this->Pressure(Ucons.density, (Ucons.energy-Ek) / Wprim.density);
}


template class EOS<1>;
template class EOS<2>;
template class EOS<3>;



