//=================================================================================================
//  Particle.cpp
//  Contains all functions relating to particle type info.
//
//  This file is part of GANDALF :
//  Graphical Astrophysics code for N-body statics And Lagrangian Fluids
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
#include "Particle.h"


//=================================================================================================
//  ParticleTypeRegister:;ParticleTypeRegister
/// \brief  Constructor setting the parameters based on the simulation params
/// \author R. A. Booth
/// \date   31/3/2016
//=================================================================================================
ParticleTypeRegister::ParticleTypeRegister(Parameters* sim_params) {

  // TODO: Add simulation parameters to control these options
  for (int k=0; k<Ntypes; k++) {
    _types[k] = ParticleTypeInfo() ;
    gravmask[k] = false;
  }

  const bool hydro_forces = sim_params->intparams["hydro_forces"];

  // Set flags for gas particle type
  //-----------------------------------------------------------------------------------------------
  _types[gas_type].hydro_forces        = hydro_forces;
  _types[gas_type].self_gravity        = true;
  _types[gas_type].drag_forces         = true;
  _types[gas_type].hmask[gas_type]     = true;
  _types[gas_type].hmask[cdm_type]     = true;
  _types[gas_type].hydromask[gas_type] = hydro_forces;
  _types[gas_type].dragmask[dust_type] = true;

  // Set flags for cdm particle type
  //-----------------------------------------------------------------------------------------------
  _types[cdm_type].self_gravity        = true;
  _types[cdm_type].hmask[gas_type]     = true;
  _types[cdm_type].hmask[cdm_type]     = true;

  // Set flags for dust particle type
  //-----------------------------------------------------------------------------------------------
  _types[dust_type].self_gravity       = true;
  _types[dust_type].drag_forces        = true;
  _types[dust_type].hmask[dust_type]   = true;
  _types[dust_type].dragmask[gas_type] = true;

  // Set flags for gravity
  //----------------------------------------------------------------------------------------------
  gravmask[gas_type]  = true;
  gravmask[cdm_type]  = true;
  gravmask[dust_type] = true;

}
