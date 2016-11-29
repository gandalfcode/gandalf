//=================================================================================================
//  TestIc.cpp
//  Class for generating initial conditions for simple turbulent core simulations.
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
//  TestIc::TestIc
/// Set-up SILCC-type simulation initial conditions.
//=================================================================================================
template <int ndim>
TestIc<ndim>::TestIc(Simulation<ndim>* _sim, Hydrodynamics<ndim>* _hydro, FLOAT _invndim) :
  Ic<ndim>(_sim, _hydro, _invndim)
{
  if (simparams->intparams["dimensionless"] == 0) {
    ExceptionHandler::getIstance().raise("dimensionless units required");
  }

  // Compute total mass inside simulation box
  Box<ndim> box;
  for (int k=0; k<ndim; k++) box.min[k] = simbox.min[k];
  for (int k=0; k<ndim; k++) box.max[k] = simbox.max[k];
  mtot = this->CalculateMassInBox(box);

  std::cout << "mtot : " << mtot*simunits.m.outscale << " " << simunits.m.outunit << std::endl;

}



//=================================================================================================
//  Silcc::Generate
/// Set-up SILCC-type simulation initial conditions.
//=================================================================================================
template <int ndim>
void TestIc<ndim>::Generate(void)
{
  string ic = simparams->stringparams["ic"];

  debug2("[TestIc::Generate");

  //-----------------------------------------------------------------------------------------------
  int i,k;                          // Particle and dimension counters
  int Nbox;                         // No. of particles in box
  int Nlattice[3];                  // Particles per dimension for LHS lattice
  FLOAT volume;                     // Volume of box
  DomainBox<ndim>& simbox = sim->simbox;

  // Local copy of important parameters
  int Npart = simparams->intparams["Nhydro"];
  Nlattice[0] = simparams->intparams["Nlattice1[0]"];
  Nlattice[1] = simparams->intparams["Nlattice1[1]"];
  Nlattice[2] = simparams->intparams["Nlattice1[2]"];

  // Compute volume and number of particles inside box
  if (ndim == 1) {
    volume = simbox.max[0] - simbox.min[0];
    Nbox = Nlattice[0];
  }
  else if (ndim == 2) {
    volume = (simbox.max[0] - simbox.min[0])*(simbox.max[1] - simbox.min[1]);
    Nbox = Nlattice[0]*Nlattice[1];
  }
  else if (ndim == 3) {
    volume = (simbox.max[0] - simbox.min[0])*
      (simbox.max[1] - simbox.min[1])*(simbox.max[2] - simbox.min[2]);
    Nbox = Nlattice[0]*Nlattice[1]*Nlattice[2];
  }

  // Add a cube of random particles defined by the simulation bounding box and
  // depending on the chosen particle distribution
  FLOAT *r = new FLOAT[ndim*Npart];
  //Ic<ndim>::AddRandomBox(Npart, simbox, r, sim->randnumb);
  Ic<ndim>::AddMonteCarloDensityField(Npart, simbox, r, sim->randnumb);


  // Allocate global and local memory for all particles
  hydro->Nhydro = Npart;
  sim->AllocateParticleMemory();

  // Copy positions to main array and initialise all other variables
  for (i=0; i<hydro->Nhydro; i++) {
    Particle<ndim>& part = hydro->GetParticlePointer(i);
    for (k=0; k<ndim; k++) {
      part.r[k] = r[ndim*i + k];
      part.v[k] = (FLOAT) 0.0;
      part.a[k] = (FLOAT) 0.0;
    }
    part.m = 1.0*mtot / (FLOAT) hydro->Nhydro;  //volume/ (FLOAT) hydro->Nhydro;
    part.h = hydro->h_fac*pow(volume / (FLOAT) hydro->Nhydro,invndim);
    part.u = (FLOAT) 1.5;
    part.iorig = i;
  }

  sim->initial_h_provided = true;

  delete[] r;

  return;
}




//=================================================================================================
//  Silcc::GetValue
/// Returns the value of the requested quantity at the given position.
//=================================================================================================
template <int ndim>
FLOAT TestIc<ndim>::GetValue
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
    if (ndim == 1) {
      if (fabs(r[0]) < 0.25) return 1.0 + 1.0*r[0];
      else return 0.0;
      //return 1.0 + 0.5*sin(r[0]*twopi);  // + r[0]*r[0] + r[1]*r[1];
      //return 1.0 + r[0]; // + r[0]*r[0];
    }
    else if (ndim == 2 || ndim == 3) {
      //std::cout << "DENSITY : " << 1.0 + 0.1*r[0] << std::endl;
      //return 1.0 + 0.5*sin(r[0]*twopi);  // + r[0]*r[0] + r[1]*r[1];
      //if (r[0] < 0.25 || r[0] > 0.75) return 0.5;
      //else return 1.0;
      if (fabs(r[0]) < 0.25) return 1.0 + 0.0*r[0];
      else return 0.0;

      //return
      //if (r[0] < 0.5) return 1.0 + 0.5*(r[0] - 0.25);
      //else return 1.0 + 0.5*(r[0] - 0.75);
    }
  }
  else {
    std::cout << "Invalid string variable for Silcc::GetValue" << std::endl;
    return 0.0;
  }
}



template class TestIc<1>;
template class TestIc<2>;
template class TestIc<3>;
