//=================================================================================================
//  Shock2DIc.cpp
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
//  Shock2DIc::Shock2DIc
/// Set-up SILCC-type simulation initial conditions.
//=================================================================================================
template <int ndim>
Shock2DIc<ndim>::Shock2DIc(Simulation<ndim>* _sim, FLOAT _invndim) :
  Ic<ndim>(_sim, _invndim)
{
  // Some sanity checking to ensure dimensionless units are used
  if (simparams->intparams["dimensionless"] == 0) {
    ExceptionHandler::getIstance().raise("dimensionless units not permitted");
  }
  if (ndim == 1) {
    ExceptionHandler::getIstance().raise("2D ShockTube problem only works in 2 or 3D");
  }
}



//=================================================================================================
//  Silcc::Generate
/// Set-up SILCC-type simulation initial conditions.
//=================================================================================================
template <int ndim>
void Shock2DIc<ndim>::Generate(void)
{
  debug2("[Shock2DIc::Generate]");

  // Setup constants
  double rho0 = 1, P0 = 1 ;
  double rho1 = 0.125, P1 = 0.14;

  double u0 = P0 / (rho0*(simparams->floatparams["gamma_eos"]-1)) ;
  double u1 = P1 / (rho1*(simparams->floatparams["gamma_eos"]-1)) ;

  // Setup the lattices for the high / low density regions
  int Nlattice1[3] = {
      simparams->intparams["Nlattice1[0]"],
      simparams->intparams["Nlattice1[1]"],
      simparams->intparams["Nlattice1[2]"] } ;
  string particle_dist = simparams->stringparams["particle_distribution"];

  // Do the high-density region first
  double volume =1;
  int Nbox1 = 1 ;
  for (int i=0; i < ndim; i++) {
    Nbox1   *= Nlattice1[i];
    volume *= icBox.max[i]-icBox.min[i] ;
  }


  std::vector<FLOAT> _r1(ndim*Nbox1,0);
  FLOAT* r1 = &(_r1[0]) ;

  // Add a cube of particles defined by the simulation bounding box and
  // depending on the chosen particle distribution
  if (particle_dist == "random") {
    Ic<ndim>::AddRandomBox(Nbox1, icBox, r1, sim->randnumb);
  }
  else if (particle_dist == "cubic_lattice") {
    Ic<ndim>::AddCubicLattice(Nbox1, Nlattice1, icBox, true, r1);
  }
  else if (particle_dist == "hexagonal_lattice") {
    Ic<ndim>::AddHexagonalLattice(Nbox1, Nlattice1, icBox, true, r1);
  }
  else {
    string message = "Invalid particle distribution option";
    ExceptionHandler::getIstance().raise(message);
  }


  // For the low density region we use either a constant separation
  // or mass (approximately).
  std::vector<FLOAT> _r2;
  int Nbox2 ;
  FLOAT* r2;
  if (simparams->intparams["use_fixed_spacing"]) {
    r2 = r1 ;
    Nbox2 = Nbox1;
  }
  else {
    double ratio = std::pow(rho1/rho0, 1./ndim) ;

    int Nlattice2[3] ;
    Nbox2 = 1 ;
    for (int i=0 ; i<ndim; ++i) {
      Nlattice2[i] = Nlattice1[i] * ratio ;
      Nbox2 *= Nlattice2[i] ;
    }
    _r2.resize(ndim*Nbox2) ;
    r2 = &(_r2[0]) ;

    // Add a cube of particles defined by the simulation bounding box and
    // depending on the chosen particle distribution
    if (particle_dist == "random") {
      Ic<ndim>::AddRandomBox(Nbox2, icBox, r2, sim->randnumb);
    }
    else if (particle_dist == "cubic_lattice") {
      Ic<ndim>::AddCubicLattice(Nbox2, Nlattice2, icBox, true, r2);
    }
    else if (particle_dist == "hexagonal_lattice") {
      Ic<ndim>::AddHexagonalLattice(Nbox2, Nlattice2, icBox, true, r2);
    }
    else {
      string message = "Invalid particle distribution option";
      ExceptionHandler::getIstance().raise(message);
    }
  }

  double scale[2] ;
  for (int i=0; i < 2; ++i)
    scale[i] = 1 / (icBox.max[i] - icBox.min[i]) ;

  // Now work out the total number of particles needed for the initial conditions
  int Ntot = 0, Ntot1 = 0, Ntot2= 0;
  for (int i=0; i < Nbox1; ++i) {
    double x = r1[ndim*i] * scale[0] + r1[ndim*i + 1] * scale[1] ;

    if (x >= 0.5)  Ntot1++ ;
  }
  for (int i=0; i < Nbox2; ++i) {
    double x = r2[ndim*i] * scale[0] + r2[ndim*i + 1] * scale[1] ;

    if (x < 0.5)  Ntot2++ ;
  }
  Ntot = Ntot1 + Ntot2 ;

  double m0 = 0.875*volume*rho0 / Ntot1;
  double m1 = 0.125*volume*rho1 / Ntot2;


  // Allocate local and main particle memory
  hydro->Nhydro = Ntot;

  bool dusty_shock = simparams->stringparams["dust_forces"] != "none";
  if (dusty_shock) hydro->Nhydro *= 2;
  sim->AllocateParticleMemory();


  // Now set the particle properties
  int j = 0;
  for (int i=0; i<Nbox1; i++) {

    double x = r1[ndim*i] * scale[0] + r1[ndim*i + 1] * scale[1] ;

    if (x >= 0.5) {
      Particle<ndim>& part = hydro->GetParticlePointer(j++);

      for (int k=0; k<ndim; k++) part.r[k] = r1[ndim*i + k];
      for (int k=0; k<ndim; k++) part.v[k] = 0.0;
      part.m = m0 ;
      part.h = hydro->h_fac * pow(part.m/rho0,1./ndim);
      part.u = u0;
      part.ptype = gas_type ;
    }
  }
  for (int i=0; i<Nbox2; i++) {

    double x = r2[ndim*i] * scale[0] + r2[ndim*i + 1] * scale[1] ;

    if (x < 0.5) {
      Particle<ndim>& part = hydro->GetParticlePointer(j++);

      for (int k=0; k<ndim; k++) part.r[k] = r2[ndim*i + k];
      for (int k=0; k<ndim; k++) part.v[k] = 0.0;
      part.m = m1 ;
      part.h = hydro->h_fac * pow(part.m/rho0,1./ndim);
      part.u = u1;
      part.ptype = gas_type ;
    }
  }
  assert(j == Ntot);

  // Add the dust lattice
  if (dusty_shock){
    FLOAT d2g = simparams->floatparams["dust_mass_factor"] ;
    for (int j = 0; j < Ntot; ++j){
      Particle<ndim>& pg = hydro->GetParticlePointer(j) ;
      Particle<ndim>& pd = hydro->GetParticlePointer(j+Ntot) ;
      pd = pg ;
      pd.ptype = dust_type ;
      pd.m *= d2g ;
      pd.u = 0 ;
      pd.h_dust = pd.h ;
    }
  }

  // Set initial smoothing lengths and create initial ghost particles
  //-----------------------------------------------------------------------------------------------
  hydro->Nghost = 0;
  hydro->Ntot = hydro->Nhydro;
  for (int i=0; i<hydro->Nhydro; i++) hydro->GetParticlePointer(i).flags.set(active);

  sim->initial_h_provided = true;
  sim->rebuild_tree = true;

  return;
}



template class Shock2DIc<1>;
template class Shock2DIc<2>;
template class Shock2DIc<3>;
