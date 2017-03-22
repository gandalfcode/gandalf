//=================================================================================================
//  PolytropicEOS.cpp
//  Contains all function definitions for the Polytropic Equation of state.
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


//=================================================================================================
//  Polytropic::Polytropic()
/// Default constructor for Polytropic EOS.  Passes and sets important
/// thermal physics variables, as well as scaling to dimensionless units.
//=================================================================================================
template <int ndim>
Polytropic<ndim>::Polytropic(Parameters* simparams, SimUnits *units):
  EOS<ndim>(simparams->floatparams["eta_eos"], simparams->floatparams["gamma_eos"]),
  Kpoly(simparams->floatparams["Kpoly"])
{
}



//=================================================================================================
//  Polytropic::Polytropic()
/// Polytropic EOS destructor
//=================================================================================================
template <int ndim>
Polytropic<ndim>::~Polytropic()
{
}



//=================================================================================================
//  Polytropic::Pressure
/// Calculates and returns thermal pressure of referenced particle
//=================================================================================================
template <int ndim>
FLOAT Polytropic<ndim>::Pressure(Particle<ndim> &part)
{
  return Kpoly*pow(part.rho, eta);
}



//=================================================================================================
//  Polytropic::EntropicFunction
/// Calculates and returns value of Entropic function (= P/rho^gamma) for referenced particle
//=================================================================================================
template <int ndim>
FLOAT Polytropic<ndim>::EntropicFunction(Particle<ndim> &part)
{
  return gammam1*part.u*pow(part.rho,(FLOAT) 1.0 - gamma);
}



//=================================================================================================
//  Polytropic::SoundSpeed
/// Returns Polytropic sound speed of referenced SPH particle
//=================================================================================================
template <int ndim>
FLOAT Polytropic<ndim>::SoundSpeed(Particle<ndim> &part)
{
  return sqrt(gammam1*part.u);
}



//=================================================================================================
//  Polytropic::SpecificInternalEnergy
/// Returns specific internal energy of referenced SPH particle
//=================================================================================================
template <int ndim>
FLOAT Polytropic<ndim>::SpecificInternalEnergy(Particle<ndim> &part)
{
  return Kpoly*pow(part.rho, gammam1)/gammam1;
}



//=================================================================================================
//  Polytropic::Temperature
/// Return Polytropic temperature, temp0, for referenced SPH particle
//=================================================================================================
template <int ndim>
FLOAT Polytropic<ndim>::Temperature(Particle<ndim> &part)
{
  return gammam1*part.u;
}



template class Polytropic<1>;
template class Polytropic<2>;
template class Polytropic<3>;
