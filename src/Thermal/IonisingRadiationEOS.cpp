//=================================================================================================
//  IonisingRadiationEOS.cpp
//  Contains all function definitions for a equation of state for ionising radiation module.
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
//  IonisingRadiation::IonisingRadiation()
/// Default constructor for ionising radiation EOS.  Passes and sets important thermal physics
/// variables, as well as setting internal EOS object for neutral gas.
//=================================================================================================
template <int ndim>
IonisingRadiation<ndim>::IonisingRadiation
 (string gas_eos, FLOAT _temp0, FLOAT _mu_bar, FLOAT _gamma, FLOAT _rho_bary, SimUnits *units) :
  EOS<ndim>(_gamma)
{
  // Set 'internal' EOS for non-ionised gas
  if (gas_eos == "energy_eqn" || gas_eos == "constant_temp") {
    eos = new Adiabatic<ndim>(_temp0, _mu_bar, _gamma);
  }
  else if (gas_eos == "isothermal") {
    eos = new Isothermal<ndim>(_temp0, _mu_bar, _gamma, units);
  }
  else if (gas_eos == "barotropic") {
    eos = new Barotropic<ndim>(_temp0, _mu_bar, _gamma, _rho_bary, units);
  }
  else {
    string message = "Unrecognised parameter : gas_eos = " + gas_eos;
    ExceptionHandler::getIstance().raise(message);
  }
  temp0 = _temp0/units->temp.outscale;
  mu_bar = _mu_bar;
}



//=================================================================================================
//  IonisingRadiation::~IonisingRadiation()
/// IonisingRadiation EOS destructor
//=================================================================================================
template <int ndim>
IonisingRadiation<ndim>::~IonisingRadiation()
{
}

//=================================================================================================
//  IonisingRadiation::EntropicFunction
/// Calculates and returns value of Entropic function (= P/rho^gamma) for referenced particle
//=================================================================================================
template <int ndim>
FLOAT IonisingRadiation<ndim>::EntropicFunction(Particle<ndim> &part)
{
  //return gammam1*part.u*pow(part.rho,(FLOAT) 1.0 - gamma);
  return eos->EntropicFunction(part);
}



//=================================================================================================
//  IonisingRadiation::SoundSpeed
/// Returns isothermal sound speed of SPH particle
//=================================================================================================
template <int ndim>
FLOAT IonisingRadiation<ndim>::SoundSpeed(Particle<ndim> &part)
{
  //return sqrt(gammam1*part.u);
  return eos->SoundSpeed(part);
}



//=================================================================================================
//  IonisingRadiation::SpecificInternalEnergy
/// Returns specific internal energy
//=================================================================================================
template <int ndim>
FLOAT IonisingRadiation<ndim>::SpecificInternalEnergy(Particle<ndim> &part)
{
  // Checks if particle's internal energy has been changed by the ionisation routine
  // If it has it compares this new internal energy to that of the EOS and chooses the largest.
  if (part.ionstate != 0) {
    FLOAT non_ionised = eos->SpecificInternalEnergy(part);
    if (part.u > non_ionised) return part.u;
    else return non_ionised;
  }
  else return eos->SpecificInternalEnergy(part);
}



//=================================================================================================
//  IonisingRadiation::Temperature
/// Returns temperature of particle.
//=================================================================================================
template <int ndim>
FLOAT IonisingRadiation<ndim>::Temperature(Particle<ndim> &part)
{
  return eos->Temperature(part);
}



template class IonisingRadiation<1>;
template class IonisingRadiation<2>;
template class IonisingRadiation<3>;
