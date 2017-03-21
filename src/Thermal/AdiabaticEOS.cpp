//=================================================================================================
//  AdiabaticEOS.cpp
//  Contains all function definitions for the Adiabatic Equation of state.
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
//  Adiabatic::Adiabatic()
/// Default constructor for perfect gas EOS.  Passes and sets important thermal physics variables.
//=================================================================================================
template <int ndim>
Adiabatic<ndim>::Adiabatic(Parameters* simparams, SimUnits *units):
  EOS<ndim>(simparams->floatparams["gamma_eos"], simparams->floatparams["gamma_eos"]),
  mu_bar(simparams->floatparams["mu_bar"])
{
}



//=================================================================================================
//  Adiabatic::Adiabatic()
//=================================================================================================
template <int ndim>
Adiabatic<ndim>::~Adiabatic()
{
}


//=================================================================================================
//  Adiabatic::EntropicFunction
/// Calculates and returns value of Entropic function (= P/rho^gamma) for referenced particle
//=================================================================================================
template <int ndim>
FLOAT Adiabatic<ndim>::EntropicFunction(Particle<ndim> &part)
{
  return gammam1*part.u*pow(part.rho,(FLOAT) 1.0 - gamma);
}



//=================================================================================================
//  Adiabatic::SoundSpeed
/// Returns adiabatic sound speed of particle
//=================================================================================================
template <int ndim>
FLOAT Adiabatic<ndim>::SoundSpeed(Particle<ndim> &part)
{
  return sqrt(gamma*gammam1*part.u);
}



//=================================================================================================
//  Adiabatic::SpecificInternalEnergy
/// Returns specific internal energy of particle
//=================================================================================================
template <int ndim>
FLOAT Adiabatic<ndim>::SpecificInternalEnergy(Particle<ndim> &part)
{
  return part.u;
}



//=================================================================================================
//  Adiabatic::Temperature
/// Returns temperature of particle
//=================================================================================================
template <int ndim>
FLOAT Adiabatic<ndim>::Temperature(Particle<ndim> &part)
{
  return gammam1*part.u;
}



template class Adiabatic<1>;
template class Adiabatic<2>;
template class Adiabatic<3>;
