//=================================================================================================
//  PolytropeIc.cpp
//  Subroutines used to generate initial conditions with a Polytrope density field.
//  Can either solve for the general of isothermal Lane-Emden equations.
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


#include <fstream>
#include <sstream>
#include "Precision.h"
#include "Debug.h"
#include "Ic.h"
using namespace std;



//=================================================================================================
//  PolytropeIc::PolytropeIc
/// Constructor for PolytropicIc class
//=================================================================================================
template <int ndim>
PolytropeIc<ndim>::PolytropeIc(Simulation<ndim>* _sim, Hydrodynamics<ndim>* _hydro, FLOAT _invndim) :
  Ic<ndim>(_sim, _hydro, _invndim)
{
  int Ntablemax = 1000;
  xiArray    = new FLOAT[Ntablemax];
  psiArray   = new FLOAT[Ntablemax];
  phiArray   = new FLOAT[Ntablemax];
  muArray    = new FLOAT[Ntablemax];
  pressArray = new FLOAT[Ntablemax];
  rhoArray   = new FLOAT[Ntablemax];
  thetaArray = new FLOAT[Ntablemax];
  massArray  = new FLOAT[Ntablemax];
}



//=================================================================================================
//  PolytropeIc::~PolytropeIc
/// Constructor for PolytropicIc class
//=================================================================================================
template <int ndim>
PolytropeIc<ndim>::~PolytropeIc()
{
}



//=================================================================================================
//  PolytropeIc::Generate
/// Create the basic particle distribution before stretching to the correct density profile
//=================================================================================================
template <int ndim>
void PolytropeIc<ndim>::Generate(void)
{
  // Only compile for 3-dimensional case
  //-----------------------------------------------------------------------------------------------
  if (ndim == 3) {

    int i;                               // Particle counter
    int k;                               // Dimension counter
    FLOAT mp;                            // Mass of individual particle
    FLOAT z;                             // z-position of newly inserted particle

    // Create local copies of initial conditions parameters
    int Npart      = simparams->intparams["Nhydro"];
    FLOAT eta_eos  = simparams->floatparams["eta_eos"];
    FLOAT gammaone = simparams->floatparams["gamma_eos"] - (FLOAT) 1.0;
    string eos     = simparams->stringparams["gas_eos"];

    debug2("[PolytropeIc::Generate]");

    // Check parameters are valid
    if (eta_eos < (FLOAT) 1.0) {
      ExceptionHandler::getIstance().raise("Invalid value for eta_eos (< 1)");
    }

    // Allocate local and main particle memory
    hydro->Nhydro = Npart;
    sim->AllocateParticleMemory();
    //mp = m_box / (FLOAT) Npart;

    // Calculate cloud parameters depending on the chosen equation of state
    if (eos == "isothermal" || fabs(eta_eos - 1.0) < (FLOAT) 0.0001) {
      std::cout << "Setting up isothermal polytrope" << std::endl;
    }
    else {
      std::cout << "Setting up general polytrope" << std::endl;
    }


    std::cout << "Nhydro : " << Npart << std::endl;

    // Create box of particles


    // Record particle properties in main memory
    // Stretch box to occupy correct volume specified in parameters file
    for (i=0; i<hydro->Nhydro; i++) {
      Particle<ndim>& part = hydro->GetParticlePointer(i);

      //for (int k=0; k<ndim; k++)
      for (k=0; k<ndim; k++) part.v[k] = (FLOAT) 0.0;
      part.m = mp;
      part.u = (FLOAT) 1.5;
      //part.h = hydro->h_fac*powf(mp/rho,invndim);
      part.ptype = gas_type;

    }

  }
  //-----------------------------------------------------------------------------------------------

  return;
}



//=================================================================================================
//  PolytropeIc::GetValue
/// ...
//=================================================================================================
template <int ndim>
FLOAT PolytropeIc<ndim>::GetValue
 (const std::string var,
  const FLOAT r[ndim])
{
  if (var == "x") {
    return r[0];
  }
  else if (ndim > 1 && var == "y") {
    return r[1];
  }
  else if (ndim > 2 && var == "z") {
    return r[2];
  }
  else if (var == "rho") {
    return (FLOAT) 0.0;
  }
  else {
    std::cout << "Invalid string variable for PolytropeIc::GetValue" << std::endl;
    return 0.0;
  }
}



template class PolytropeIc<1>;
template class PolytropeIc<2>;
template class PolytropeIc<3>;
