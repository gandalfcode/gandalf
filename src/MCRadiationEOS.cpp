//=============================================================================
//  MCRadiationEOS.cpp
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
//  MCRadiationEOS::MCRadiationEOS()
/// Default constructor for barotropic EOS.  Passes and sets important 
/// thermal physics variables.
//=============================================================================
template <int ndim>
MCRadiationEOS<ndim>::MCRadiationEOS(string gas_eos, FLOAT temp0aux, 
				     FLOAT tempionaux, FLOAT mu_bar_aux, 
                                     FLOAT mu_ion_aux, FLOAT gamma_aux,
                                     FLOAT rho_bary_aux, SimUnits *units):
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
  mu_bar = mu_bar_aux;
  mu_ion = mu_ion_aux;
  temp0 = temp0aux/units->temp.outscale;
  temp_ion = tempionaux/units->temp.outscale;
}



//=============================================================================
//  MCRadiationEOS::~MCRadiationEOS()
/// MCRadiationEOS EOS destructor
//=============================================================================
template <int ndim>
MCRadiationEOS<ndim>::~MCRadiationEOS()
{
}



//=============================================================================
//  MCRadiationEOS::Pressure
/// Calculates and returns thermal pressure of referenced particle
//=============================================================================
template <int ndim>
FLOAT MCRadiationEOS<ndim>::Pressure(SphParticle<ndim> &part)
{
  //return gammam1*part.rho*part.u;
  //;eos->Pressure(part);
  return gammam1*part.rho*part.u;
}



//=============================================================================
//  MCRadiationEOS::EntropicFunction
/// Calculates and returns value of Entropic function (= P/rho^gamma) for 
/// referenced particle
//=============================================================================
template <int ndim>
FLOAT MCRadiationEOS<ndim>::EntropicFunction(SphParticle<ndim> &part)
{
  //return gammam1*part.u*pow(part.rho,(FLOAT) 1.0 - gamma);
  return eos->EntropicFunction(part);
}



//=============================================================================
//  MCRadiationEOS::SoundSpeed
/// Returns isothermal sound speed of SPH particle
//=============================================================================
template <int ndim>
FLOAT MCRadiationEOS<ndim>::SoundSpeed(SphParticle<ndim> &part)
{
  //return sqrt(gammam1*part.u);
  return part.ionfrac*sqrt(gammam1*part.u) + 
    (1.0 - part.ionfrac)*eos->SoundSpeed(part);
}



//=============================================================================
//  MCRadiationEOS::SpecificInternalEnergy
/// Returns specific internal energy
//=============================================================================
template <int ndim>
FLOAT MCRadiationEOS<ndim>::SpecificInternalEnergy(SphParticle<ndim> &part)
{
  //cout << "u : " << part.ionfrac << "  " << temp_ion << "   " << gammam1 << "   " << mu_ion << "  "  
  //   << (1.0 - part.ionfrac) << "   " <<  eos->SpecificInternalEnergy(part) << endl;
  return part.ionfrac*temp_ion/gammam1/mu_ion + 
    (1.0 - part.ionfrac)*eos->SpecificInternalEnergy(part);
}



//=============================================================================
//  MCRadiationEOS::Temperature
/// Returns temperature of particle.  Approximates gas in the isothermal 
/// regime (T = temp0 for rho << rho_bary) and in the optically thick 
/// adiabatic phase (T = const*rho^{gamma - 1} for rho >> rho_bary).
//=============================================================================
template <int ndim>
FLOAT MCRadiationEOS<ndim>::Temperature(SphParticle<ndim> &part)
{
  return part.ionfrac*temp_ion + (1.0 - part.ionfrac)*eos->Temperature(part);
}



template class MCRadiationEOS<1>;
template class MCRadiationEOS<2>;
template class MCRadiationEOS<3>;

