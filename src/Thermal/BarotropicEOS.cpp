//=================================================================================================
//  BarotropicEOS.cpp
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
//=================================================================================================


#include <math.h>
#include "EOS.h"
#include "Particle.h"


//=================================================================================================
//  Barotropic::Barotropic()
/// Default constructor for barotropic EOS.  Passes and sets important thermal physics variables.
//=================================================================================================
template <int ndim>
Barotropic<ndim>::Barotropic(Parameters* simparams, SimUnits *units):
  EOS<ndim>(simparams->floatparams["gamma_eos"], simparams->floatparams["gamma_eos"])
{
  temp0    = simparams->floatparams["temp0"]/units->temp.outscale;
  mu_bar   = simparams->floatparams["mu_bar"];
  rho_bary = simparams->floatparams["rho_bary"]/units->rho.outscale/units->rho.outcgs;
  invrho_bary = 1.0/rho_bary;
}



//=================================================================================================
//  Barotropic::~Barotropic()
/// Barotropic EOS destructor
//=================================================================================================
template <int ndim>
Barotropic<ndim>::~Barotropic()
{
}



//=================================================================================================
//  Barotropic::EntropicFunction
/// Calculates and returns value of Entropic function (= P/rho^gamma) for referenced particle.
//=================================================================================================
template <int ndim>
FLOAT Barotropic<ndim>::EntropicFunction(Particle<ndim> &part)
{
  return gammam1*part.u*pow(part.rho,(FLOAT) 1.0 - gamma);
}



//=================================================================================================
//  Barotropic::SoundSpeed
/// Returns sound speed of SPH particle
//=================================================================================================
template <int ndim>
FLOAT Barotropic<ndim>::SoundSpeed(Particle<ndim> &part)
{
  return sqrt(gammam1*part.u);
}



//=================================================================================================
//  Barotropic::SpecificInternalEnergy
/// Returns specific internal energy
//=================================================================================================
template <int ndim>
FLOAT Barotropic<ndim>::SpecificInternalEnergy(Particle<ndim> &part)
{
  return temp0*(1.0 + pow(part.rho*invrho_bary,gammam1))/gammam1/mu_bar;
}



//=================================================================================================
//  Barotropic::Temperature
/// Returns temperature of particle.  Approximates gas in the isothermal
/// regime (T = temp0 for rho << rho_bary) and in the optically thick
/// adiabatic phase (T = const*rho^{gamma - 1} for rho >> rho_bary).
//=================================================================================================
template <int ndim>
FLOAT Barotropic<ndim>::Temperature(Particle<ndim> &part)
{
  return temp0*(1.0 + pow(part.rho*invrho_bary,gammam1));
}



template class Barotropic<1>;
template class Barotropic<2>;
template class Barotropic<3>;
