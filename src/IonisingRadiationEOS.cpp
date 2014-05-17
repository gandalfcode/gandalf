//=============================================================================
//  IonisingRadiationEOS.cpp
//  Contains all function definitions for a barotropic Equation of state of 
//  the form T = temp0*(1 + (rho/rho_bary)^{gamma - 1}).
//  Used for star formation simulations to approximate the combined isothermal 
//  and optically-thich adiabatic regimes of the gas collapse phase.
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
//=============================================================================


#include <math.h>
#include "EOS.h"
#include "Sph.h"


//=============================================================================
//  IonisingRadiation::IonisingRadiation()
/// Default constructor for barotropic EOS.  Passes and sets important 
/// thermal physics variables.
//=============================================================================
template <int ndim>
IonisingRadiation<ndim>::IonisingRadiation(string gas_eos, FLOAT temp0aux, 
                                           FLOAT mu_bar_aux, FLOAT gamma_aux,
                                           FLOAT rho_bary_aux, SimUnits *units,
                                           SphNeighbourSearch<ndim> *sphneib):
  EOS<ndim>(gamma_aux)
{
  // Set 'internal' EOS for non-ionised gas
  if (gas_eos == "energy_eqn" || gas_eos == "constant_temp")
    eos = new Adiabatic<ndim>(temp0aux,mu_bar_aux,gamma_aux);
  else if (gas_eos == "isothermal")
    eos = new Isothermal<ndim>(temp0aux,mu_bar_aux,gamma_aux,units);
  else if (gas_eos == "barotropic")
    eos = new Barotropic<ndim>(temp0aux,mu_bar_aux,gamma_aux,
                               rho_bary_aux,units);

  else {
    string message = "Unrecognised parameter : gas_eos = " + gas_eos;
    ExceptionHandler::getIstance().raise(message);
  }
  temp0 = temp0aux/units->temp.outscale;
  mu_bar = mu_bar_aux;
}



//=============================================================================
//  IonisingRadiation::~IonisingRadiation()
/// IonisingRadiation EOS destructor
//=============================================================================
template <int ndim>
IonisingRadiation<ndim>::~IonisingRadiation()
{
}



//=============================================================================
//  IonisingRadiation::Pressure
/// Calculates and returns thermal pressure of referenced particle
//=============================================================================
template <int ndim>
FLOAT IonisingRadiation<ndim>::Pressure(SphParticle<ndim> &part)
{
  //return gammam1*part.rho*part.u;
  return eos->Pressure(part);
}



//=============================================================================
//  IonisingRadiation::EntropicFunction
/// Calculates and returns value of Entropic function (= P/rho^gamma) for 
/// referenced particle
//=============================================================================
template <int ndim>
FLOAT IonisingRadiation<ndim>::EntropicFunction(SphParticle<ndim> &part)
{
  //return gammam1*part.u*pow(part.rho,(FLOAT) 1.0 - gamma);
  return eos->EntropicFunction(part);
}



//=============================================================================
//  IonisingRadiation::SoundSpeed
/// Returns isothermal sound speed of SPH particle
//=============================================================================
template <int ndim>
FLOAT IonisingRadiation<ndim>::SoundSpeed(SphParticle<ndim> &part)
{
  //return sqrt(gammam1*part.u);
  return eos->SoundSpeed(part);
}



//=============================================================================
//  IonisingRadiation::SpecificInternalEnergy
/// Returns specific internal energy
//=============================================================================
template <int ndim>
FLOAT IonisingRadiation<ndim>::SpecificInternalEnergy(SphParticle<ndim> &part)
{
  //return temp0*(1.0 + pow(part.rho*invrho_bary,gammam1))/gammam1/mu_bar;
  return eos->SpecificInternalEnergy(part);
}



//=============================================================================
//  IonisingRadiation::Temperature
/// Returns temperature of particle.  Approximates gas in the isothermal 
/// regime (T = temp0 for rho << rho_bary) and in the optically thick 
/// adiabatic phase (T = const*rho^{gamma - 1} for rho >> rho_bary).
//=============================================================================
template <int ndim>
FLOAT IonisingRadiation<ndim>::Temperature(SphParticle<ndim> &part)
{
  //return temp0*(1.0 + pow(part.rho*invrho_bary,gammam1));
  return eos->Temperature(part);
}



template class IonisingRadiation<1>;
template class IonisingRadiation<2>;
template class IonisingRadiation<3>;

