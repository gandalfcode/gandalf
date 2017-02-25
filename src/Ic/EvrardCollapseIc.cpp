//=================================================================================================
//  EvrardCollapseIc.cpp
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
//  EvrardCollapseIc::EvrardCollapseIc
/// Set-up SILCC-type simulation initial conditions.
//=================================================================================================
template <int ndim>
EvrardCollapseIc<ndim>::EvrardCollapseIc(Simulation<ndim>* _sim, FLOAT _invndim) :
  Ic<ndim>(_sim, _invndim)
{
}



//=================================================================================================
//  EvrardCollapseIc::Generate
/// ...
//=================================================================================================
template <int ndim>
void EvrardCollapseIc<ndim>::Generate(void)
{
  int i;                                   // Particle counter
  int k;                                   // Dimension counter
  FLOAT rcentre[ndim];                     // Position of sphere centre
  FLOAT *pos;                              // Positions of all particles

  // Create local copies of initial conditions parameters
  int Npart      = simparams->intparams["Nhydro"];
  FLOAT Mtot     = simparams->floatparams["mcloud"];
  FLOAT U_fac    = simparams->floatparams["thermal_energy"];
  FLOAT radius   = simparams->floatparams["radius"];
  string particle_dist = simparams->stringparams["particle_distribution"];

  cout << "Setting ICs for Evrard collapse Problem:"
     << "\n\tNumPart (desired):\t" << Npart
     << "\n\tMass:             \t" << Mtot
     << "\n\tRadius:           \t" << radius
     << "\n\tInternal Energy:  \t" << U_fac << endl;

  for (k =0; k < ndim; k++) rcentre[k] = 0 ;

  // Allocate memory
  pos = new FLOAT[ndim*Npart];

  if (particle_dist == "random") {
    Ic<ndim>::AddRandomSphere(Npart, rcentre, radius, pos, sim->randnumb);
  }
  else if (particle_dist == "cubic_lattice" || particle_dist == "hexagonal_lattice") {
    Npart = Ic<ndim>::AddLatticeSphere(Npart, rcentre, radius, particle_dist, pos, sim->randnumb);
  }
  else {
    string message = "Invalid particle distribution option";
    ExceptionHandler::getIstance().raise(message);
  }

  cout << "\n\tNpart (actual):   \t" << Npart << endl ;

  //Allocate local and main particle memory
  hydro->Nhydro = Npart;

  bool dusty_collapse = simparams->stringparams["dust_forces"] != "none" ;
  if (dusty_collapse) hydro->Nhydro *= 2;

  sim->AllocateParticleMemory();

  // Scale to the correct density profile
  for (i = 0; i < Npart; ++i) {
    Particle<ndim>& P = hydro->GetParticlePointer(i) ;

    FLOAT* ri = pos + ndim * i ;

    FLOAT r = sqrt(DotProduct(ri, ri, ndim) + small_number) ;
    FLOAT rnew = radius * r * sqrt(r) ;
    for (k=0; k < ndim; k++){
      P.r[k] = ri[k] * rnew / r ;
      P.v[k] = 0 ;
    }
    P.m = Mtot / Npart ;
    P.u = U_fac * Mtot / radius ;
    P.rho = (Mtot/ (2*pi * pow(radius, ndim))) * (radius/rnew) ;
    P.h = pow(P.m/P.rho, 1./ndim) ;
  }

  if (dusty_collapse) {
    FLOAT d2g = simparams->floatparams["dust_mass_factor"] ;
    for (i = 0; i < Npart; ++i){
      Particle<ndim>& Pg = hydro->GetParticlePointer(i) ;
      Particle<ndim>& Pd = hydro->GetParticlePointer(i+Npart) ;
      Pd = Pg ;
      Pd.m *= d2g ;

      for (k=0; k < ndim; k++) Pd.r[k] += 0.01 * Pd.h ;

      Pd.h_dust = Pd.h ;
      Pd.u = 0 ;

      Pg.ptype = gas_type ;
      Pd.ptype = dust_type ;
    }
  }

  sim->initial_h_provided = true;


  delete[] pos ;

  return;
}



template class EvrardCollapseIc<1>;
template class EvrardCollapseIc<2>;
template class EvrardCollapseIc<3>;
