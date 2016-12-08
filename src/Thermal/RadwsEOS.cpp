//=================================================================================================
//  RadwsEOS.cpp
//  Contains all function definitions for the Radws Equation of state.
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
#include "Particle.h"



//=================================================================================================
//  Radws::Radws()
/// Default constructor for perfect gas EOS.  Passes and sets important
/// thermal physics variables.
//=================================================================================================
template <int ndim>
Radws<ndim>::Radws(FLOAT temp0aux, FLOAT mu_bar_aux, FLOAT gamma_aux) : EOS<ndim> (gamma_aux)
{
  mu_bar = mu_bar_aux;
  temp0 = temp0aux;
}



//=================================================================================================
//  Radws::Radws()
//=================================================================================================
template <int ndim>
Radws<ndim>::~Radws()
{
}


//=================================================================================================
//  Radws::EntropicFunction
/// Calculates and returns value of Entropic function (= P/rho^gamma) for
/// referenced particle
//=================================================================================================
template <int ndim>
FLOAT Radws<ndim>::EntropicFunction(Particle<ndim> &part)
{
  return gammam1*part.u*pow(part.rho, (FLOAT) 1.0 - gamma);
}



//=================================================================================================
//  Radws::SoundSpeed
/// Returns adiabatic sound speed of particle
//=================================================================================================
template <int ndim>
FLOAT Radws<ndim>::SoundSpeed(Particle<ndim> &part)
{
  return sqrt(gamma*gammam1*part.u);
}



//=================================================================================================
//  Radws::SpecificInternalEnergy
/// Returns specific internal energy of particle
//=================================================================================================
template <int ndim>
FLOAT Radws<ndim>::SpecificInternalEnergy(Particle<ndim> &part)
{
  return part.u;
}



//=================================================================================================
//  Radws::Temperature
/// Returns temperature of particle
//=================================================================================================
template <int ndim>
FLOAT Radws<ndim>::Temperature(Particle<ndim> &part)
{
  return gammam1*part.u*part.mu_bar;
}



template class Radws<1>;
template class Radws<2>;
template class Radws<3>;
