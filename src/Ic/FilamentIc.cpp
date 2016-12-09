//=================================================================================================
//  FilamentIc.cpp
//  Class for generating initial conditions for simple filaments.
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
//  Filament::Filament
/// Set-up Filament-type simulation initial conditions.
//=================================================================================================
template <int ndim>
FilamentIc<ndim>::FilamentIc(Simulation<ndim>* _sim, Hydrodynamics<ndim>* _hydro, FLOAT _invndim) :
  Ic<ndim>(_sim, _hydro, _invndim)
{
  FLOAT mu_bar = simparams->floatparams["mu_bar"];
  FLOAT gammaone = simparams->floatparams["gamma_eos"] - (FLOAT) 1.0;

  // Hard-wired (for now) constants + other parameters
  n0        = simparams->floatparams["n0"];
  r0        = simparams->floatparams["r0"];
  Rfilament = simparams->floatparams["Rfilament"];
  Lfilament = simparams->floatparams["Lfilament"];
  temp0     = simparams->floatparams["temp0"];

  // Constant and derived quantities
  aconst    = (FLOAT) 10.9;
  rho0      = (FLOAT) 1000.0*m_hydrogen*mu_bar*n0;

  // Some sanity checking to ensure correct units are used for these ICs
  if (simparams->stringparams["routunit"] != "pc") {
    ExceptionHandler::getIstance().raise("r unit not set to pc");
  }
  if (simparams->stringparams["rhooutunit"] != "g_cm3") {
    ExceptionHandler::getIstance().raise("sigma unit not set to g_cm3");
  }

  // Convert any parameters to code units
  Rfilament /= simunits.r.outscale;
  Lfilament /= simunits.r.outscale;
  r0        /= simunits.r.outscale;
  rho0      /= simunits.rho.outscale;
  temp0     /= simunits.temp.outscale;

  // Compute initial internal energy of all particles
  u0 = temp0/gammaone/mu_bar;

  // Compute total mass inside simulation box
  Box<ndim> box;
  for (int k=0; k<ndim; k++) box.min[k] = simbox.min[k];
  for (int k=0; k<ndim; k++) box.max[k] = simbox.max[k];
  mtot = this->CalculateMassInBox(box);

  std::cout << "rho0 : " << rho0*simunits.rho.outscale << " " << simunits.rho.outunit << std::endl;
  std::cout << "mtot : " << mtot*simunits.m.outscale << " " << simunits.m.outunit << std::endl;

}



//=================================================================================================
//  Filament::Generate
/// Set-up Filament-type simulation initial conditions.
//=================================================================================================
template <int ndim>
void FilamentIc<ndim>::Generate(void)
{
  // Only compile for 3-dimensional case
  //-----------------------------------------------------------------------------------------------
  if (ndim == 3) {

    FLOAT Rsqd;
    int Npart      = simparams->intparams["Nhydro"];
    FLOAT gammaone = simparams->floatparams["gamma_eos"] - (FLOAT) 1.0;

    debug2("[FilamentIc::Generate]");

    // Allocate local and main particle memory
    hydro->Nhydro = Npart;
    sim->AllocateParticleMemory();
    mp = 1.0*mtot / (FLOAT) Npart;

    FLOAT *r = new FLOAT[ndim*Npart];
    Ic<ndim>::AddMonteCarloDensityField(Npart, simbox, r, sim->randnumb);

    // Copy positions to main array and initialise all other variables
    for (int i=0; i<hydro->Nhydro; i++) {
      Particle<ndim>& part = hydro->GetParticlePointer(i);
      for (int k=0; k<ndim; k++) {
        part.r[k] = r[ndim*i + k];
        part.v[k] = (FLOAT) 0.0;
        part.a[k] = (FLOAT) 0.0;
      }
      part.m = mp;
      part.u = u0;
      part.iorig = i;
      part.ptype = gas;
    }

    delete[] r;

  }
  //-----------------------------------------------------------------------------------------------

  return;
}



//=================================================================================================
//  Filament::GetValue
/// Returns the value of the requested quantity at the given position.
//=================================================================================================
template <int ndim>
FLOAT FilamentIc<ndim>::GetValue
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
    FLOAT R = sqrt(r[0]*r[0] + r[1]*r[1]);
    FLOAT radsqd = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
    FLOAT z = r[2];
    if (R < Rfilament && fabs(z) < Lfilament) {
      return rho0 / ((FLOAT) 1.0 + radsqd/r0/r0 + z*z/r0/r0/aconst/aconst);
    }
    else {
      return (FLOAT) 0.0;
    }
  }
  else {
    std::cout << "Invalid string variable for Filament::GetValue" << std::endl;
    return 0.0;
  }
}



template class FilamentIc<1>;
template class FilamentIc<2>;
template class FilamentIc<3>;
