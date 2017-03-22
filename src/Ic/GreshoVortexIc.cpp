//=================================================================================================
//  GreshoVortexIc.cpp
//  Class for generating initial conditions for Gresho-vortex simulations.
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
//  Khi::Khi
/// Set-up Khi-type simulation initial conditions.
//=================================================================================================
template <int ndim>
GreshoVortexIc<ndim>::GreshoVortexIc(Simulation<ndim>* _sim, FLOAT _invndim) :
  Ic<ndim>(_sim, _invndim)
{
  // Some sanity checking to ensure correct dimensionality and (dimensionless) units are used
  if (simparams->intparams["ndim"] != 2) {
    ExceptionHandler::getIstance().raise("Gresho vortex test only runs in 2D");
  }
  if (simparams->intparams["dimensionless"] != 1) {
    ExceptionHandler::getIstance().raise("r unit not set to pc");
  }
}



//=================================================================================================
//  Khi::Generate
/// Set-up Khi-type simulation initial conditions.
//=================================================================================================
template <int ndim>
void GreshoVortexIc<ndim>::Generate(void)
{
  // Only compile for 3-dimensional case
  //-----------------------------------------------------------------------------------------------
  if (ndim == 2) {

    int i;                                 // Particle counter
    int k;                                 // Dimension counter
    int Nbox;                              // No. of particles in fluid box 1
    int Nlattice[ndim];                    // Lattice particles in fluid box 1
    FLOAT dr_unit[ndim];                   // Unit vector
    FLOAT drmag;                           // Distance from origin
    FLOAT drsqd;                           // Distance squared from origin
    FLOAT press;                           // Local pressure
    FLOAT rotspeed;                        // Local azimuthal rotational speed
    FLOAT volume;                          // Volume of fluid box
    FLOAT *r;                              // Array of particle positions

    // Record local copies of all important parameters
    FLOAT rhofluid = (FLOAT) 1.0;
    FLOAT gammaone  = simparams->floatparams["gamma_eos"] - (FLOAT) 1.0;
    Nlattice[0] = simparams->intparams["Nlattice1[0]"];
    Nlattice[1] = simparams->intparams["Nlattice1[1]"];

    debug2("[GreshoVortexIc::Generate]");


    // Compute size and range of fluid bounding boxes
    //---------------------------------------------------------------------------------------------
    volume = (icBox.max[0] - icBox.min[0])*(icBox.max[1] - icBox.min[1]);
    Nbox = Nlattice[0]*Nlattice[1];

    // Allocate local and main particle memory and calculate positions for cubic lattice
    hydro->Nhydro = Nbox;
    sim->AllocateParticleMemory();
    r = new FLOAT[ndim*hydro->Nhydro];
    Ic<ndim>::AddCubicLattice(Nbox, Nlattice, icBox, false, r);

    for (i=0; i<Nbox; i++) {
      Particle<ndim>& part = hydro->GetParticlePointer(i);
      for (k=0; k<ndim; k++) part.r[k] = r[ndim*i + k];
      drsqd = DotProduct(part.r, part.r, ndim);
      drmag = sqrt(drsqd) + small_number;
      for (k=0; k<ndim; k++) dr_unit[k] = part.r[k]/drmag;

      // Set velocity and pressure/internal energy depending on radial position
      if (drmag < 0.2) {
        rotspeed = 5.0*drmag;
        press = 5.0 + 12.5*drsqd;
      }
      else if (drmag < 0.4) {
        rotspeed = 2.0 - 5.0*drmag;
        press = 9.0 + 12.5*drsqd - 20.0*drmag + 4.0*log(drmag/0.2);
      }
      else {
        rotspeed = 0.0;
        press = 3.0 + 4.0*log(2.0);
      }
      part.v[0] = -rotspeed*dr_unit[1];
      part.v[1] = rotspeed*dr_unit[0];

      part.m = rhofluid*volume/(FLOAT) Nbox;
      part.h = hydro->h_fac*pow(part.m/rhofluid,invndim);
      part.u = press/rhofluid/gammaone;
    }

    delete[] r;

  }
  //-----------------------------------------------------------------------------------------------

  return;
}



//=================================================================================================
//  GreshoVortexIc::GetValue
/// Returns the value of the requested quantity at the given position.
//=================================================================================================
template <int ndim>
FLOAT GreshoVortexIc<ndim>::GetValue
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
    return (FLOAT) 1.0;
  }
  else {
    std::cout << "Invalid string variable for GreshoVortexIc::GetValue" << std::endl;
    return 0.0;
  }
}



template class GreshoVortexIc<1>;
template class GreshoVortexIc<2>;
template class GreshoVortexIc<3>;
