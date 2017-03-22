//=================================================================================================
//  DustyBoxIc.cpp
//  Class for generating initial conditions for ...
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
//  DustyBoxIc::DustyBoxIc
/// Set-up SILCC-type simulation initial conditions.
//=================================================================================================
template <int ndim>
DustyBoxIc<ndim>::DustyBoxIc(Simulation<ndim>* _sim, FLOAT _invndim) :
  Ic<ndim>(_sim, _invndim)
{
}



//=================================================================================================
//  Silcc::Generate
/// Set-up SILCC-type simulation initial conditions.
//=================================================================================================
template <int ndim>
void DustyBoxIc<ndim>::Generate(void)
{
  int i;                               // Particle counter
  int k;                               // Dimension counter
  int Nbox;                            // No. of particles in box
  int Nlattice[3];                     // Lattice size
  FLOAT mbox;                          // Total mass inside simulation box
  FLOAT volume;                        // Volume of box
  FLOAT *r;                            // Positions of all particles

  // Create local copies of initial conditions parameters
  Nlattice[0]    = simparams->intparams["Nlattice1[0]"];
  Nlattice[1]    = simparams->intparams["Nlattice1[1]"];
  Nlattice[2]    = simparams->intparams["Nlattice1[2]"];
  FLOAT rhofluid = simparams->floatparams["rhofluid1"];
  FLOAT press    = simparams->floatparams["press1"];
  FLOAT gamma_m1 = simparams->floatparams["gamma_eos"] - 1;
  FLOAT v_gas    = simparams->floatparams["vfluid1[0]"] ;
  FLOAT v_dust   = simparams->floatparams["vfluid2[0]"] ;
  string particle_dist = simparams->stringparams["particle_distribution"];

  debug2("[Ic::DustyBox]");


  // Compute size and range of fluid bounding boxes
  //-----------------------------------------------------------------------------------------------
  volume = 1 ;
  Nbox = 1 ;
  for (i=0; i < ndim; ++i){
	volume *= icBox.max[i] - icBox.min[i];
	Nbox *= Nlattice[i] ;
  }
  mbox  = volume*rhofluid;


  // Allocate local and main particle memory
  hydro->Nhydro = Nbox;

  bool dusty_box = simparams->stringparams["dust_forces"] != "none" ;
  if (dusty_box) {
  	hydro->Nhydro *= 2;
  }
  else {
	  ExceptionHandler::getIstance().raise("Error: Drag forces must be enabled for dusty box test");
  }

  sim->AllocateParticleMemory();
  r = new FLOAT[ndim*Nbox];

  // Add a cube of random particles defined by the simulation bounding box and
  // depending on the chosen particle distribution
  if (particle_dist == "random") {
    Ic<ndim>::AddRandomBox(Nbox, icBox, r, sim->randnumb);
  }
  else if (particle_dist == "cubic_lattice") {
    Ic<ndim>::AddCubicLattice(Nbox, Nlattice, icBox, true, r);
  }
  else if (particle_dist == "hexagonal_lattice") {
    Ic<ndim>::AddHexagonalLattice(Nbox, Nlattice, icBox, true, r);
  }
  else {
    string message = "Invalid particle distribution option";
    ExceptionHandler::getIstance().raise(message);
  }

  // Record positions in main memory
  for (i=0; i<Nbox; i++) {
    Particle<ndim>& part = hydro->GetParticlePointer(i);
    for (k=0; k<ndim; k++) part.r[k] = r[ndim*i + k];
    for (k=0; k<ndim; k++) part.v[k] = 0.0;
    part.v[0] = v_gas ;

    part.m = mbox/(FLOAT) Nbox;
    part.h = hydro->h_fac*pow(part.m/rhofluid,invndim);
    part.u = press / (rhofluid * gamma_m1);
  }

  sim->initial_h_provided = true;


  // Add a slightly offset dust lattice
  if (dusty_box){
    FLOAT d2g = simparams->floatparams["dust_mass_factor"] ;
    for (int j = 0; j < Nbox; ++j){
      Particle<ndim>& pg = hydro->GetParticlePointer(j) ;
      Particle<ndim>& pd = hydro->GetParticlePointer(j+Nbox) ;
      pd = pg ;
      pg.ptype = gas_type ;
      pd.ptype = dust_type ;
      pd.r[0] += 0.01 * pd.h ;
      pd.v[0] = v_dust ;
      pd.m *= d2g ;
      pd.u = 0 ;
      pd.h_dust = pd.h ;
    }
  }

  sim->initial_h_provided = true;

  delete[] r;

  return;
}



template class DustyBoxIc<1>;
template class DustyBoxIc<2>;
template class DustyBoxIc<3>;
