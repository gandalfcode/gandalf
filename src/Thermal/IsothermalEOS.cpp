//=================================================================================================
//  IsothermalEOS.cpp
//  Contains all function definitions for the Isothermal Equation of state.
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


#include <math.h>
#include "EOS.h"
#include "Sph.h"


//=================================================================================================
//  Isothermal::Isothermal()
/// Default constructor for isothermal EOS.  Passes and sets important
/// thermal physics variables, as well as scaling to dimensionless units.
//=================================================================================================
template <int ndim>
Isothermal<ndim>::Isothermal(Parameters* simparams, SimUnits *units):
  EOS<ndim>(1, simparams->floatparams["gamma_eos"]),
  temp0(simparams->floatparams["temp0"]/units->temp.outscale),
  mu_bar(simparams->floatparams["mu_bar"])
{
}



//=================================================================================================
//  Isothermal::Isothermal()
/// Isothermal EOS destructor
//=================================================================================================
template <int ndim>
Isothermal<ndim>::~Isothermal()
{
}



//=================================================================================================
//  Isothermal::EntropicFunction
/// Calculates and returns value of Entropic function (= P/rho^gamma) for referenced particle
//=================================================================================================
template <int ndim>
FLOAT Isothermal<ndim>::EntropicFunction(Particle<ndim> &part)
{
  return gammam1*part.u*pow(part.rho,(FLOAT) 1.0 - gamma);
}



//=================================================================================================
//  Isothermal::SoundSpeed
/// Returns isothermal sound speed of referenced SPH particle
//=================================================================================================
template <int ndim>
FLOAT Isothermal<ndim>::SoundSpeed(Particle<ndim> &part)
{
  return sqrt(gammam1*part.u);
}



//=================================================================================================
//  Isothermal::SpecificInternalEnergy
/// Returns specific internal energy of referenced SPH particle
//=================================================================================================
template <int ndim>
FLOAT Isothermal<ndim>::SpecificInternalEnergy(Particle<ndim> &part)
{
  return temp0/gammam1/mu_bar;
}



//=================================================================================================
//  Isothermal::Temperature
/// Return isothermal temperature, temp0, for referenced SPH particle
//=================================================================================================
template <int ndim>
FLOAT Isothermal<ndim>::Temperature(Particle<ndim> &part)
{
  return temp0;
}



template class Isothermal<1>;
template class Isothermal<2>;
template class Isothermal<3>;
