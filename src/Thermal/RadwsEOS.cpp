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
/// Default constructor for the RadWS equation of state which employs the use of an opacity table.
//=================================================================================================
template <int ndim>
Radws<ndim>::Radws(Parameters* simparams, SimUnits *units) :
  EOS<ndim>(simparams->floatparams["gamma_eos"])
{
  opacity_table = new OpacityTable<ndim>(simparams, units);
}



//=================================================================================================
//  Radws::~Radws()
//=================================================================================================
template <int ndim>
Radws<ndim>::~Radws()
{
  delete opacity_table;
}


//=================================================================================================
//  Radws::Pressure
/// Returns pressure of a particle based on a variable gamma
//=================================================================================================
template <int ndim>
FLOAT Radws<ndim>::Pressure(const EosParticleProxy<ndim>&part)
{
  FLOAT gamma = opacity_table->GetGamma(part);
  return (gamma - 1.0) * part.rho * part.u;
}



//=================================================================================================
//  Radws::EntropicFunction
/// Calculates and returns value of Entropic function (= P/rho^gamma) for
/// referenced particle
//=================================================================================================
template <int ndim>
FLOAT Radws<ndim>::EntropicFunction(const EosParticleProxy<ndim>&part)
{
  FLOAT gamma = opacity_table->GetGamma(part);
  return (gamma - 1.0)*part.u*pow(part.rho, (FLOAT) (1.0 - gamma));
}



//=================================================================================================
//  Radws::SoundSpeed
/// Returns adiabatic sound speed of particle
//=================================================================================================
template <int ndim>
FLOAT Radws<ndim>::SoundSpeed(const EosParticleProxy<ndim>&part)
{
  FLOAT gamma = opacity_table->GetGamma(part);
  FLOAT gamma1 = opacity_table->GetGamma1(part);
  return sqrt(gamma1*(gamma - 1.0)*part.u);
}



//=================================================================================================
//  Radws::SpecificInternalEnergy
/// Returns specific internal energy of particle
//=================================================================================================
template <int ndim>
FLOAT Radws<ndim>::SpecificInternalEnergy(const EosParticleProxy<ndim>&part)
{
  return part.u;
}



//=================================================================================================
//  Radws::Temperature
/// Returns temperature of particle
//=================================================================================================
template <int ndim>
FLOAT Radws<ndim>::Temperature(const EosParticleProxy<ndim>&part)
{
  FLOAT gamma = opacity_table->GetGamma(part);
  FLOAT mu_bar = opacity_table->GetMuBar(part);
  return (gamma - 1.0)*part.u*mu_bar;
}

//=================================================================================================
//  Radws::InternalEnergyFromPressure
/// Invert pressure to get internal energy
//==============================================================================================
template <int ndim>
FLOAT Radws<ndim>::InternalEnergyFromPressure(const EosParticleProxy<ndim>&p){
  return opacity_table->GetEnergyFromPressure(p.rho, p.p);
}


template class Radws<1>;
template class Radws<2>;
template class Radws<3>;
