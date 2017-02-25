//=================================================================================================
//  NohIc.cpp
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
//  NohIc::NohIc
/// Set-up SILCC-type simulation initial conditions.
//=================================================================================================
template <int ndim>
NohIc<ndim>::NohIc(Simulation<ndim>* _sim, FLOAT _invndim) :
  Ic<ndim>(_sim, _invndim)
{
}



//=================================================================================================
//  NohIc::Generate
/// ...
//=================================================================================================
template <int ndim>
void NohIc<ndim>::Generate(void)
{
  int i;                            // Particle counter
  int k;                            // Dimension counter
  int Nsphere;                      // Actual number of particles in sphere
  FLOAT dr[ndim];                   // Relative position vector
  FLOAT drmag;                      // Distance
  FLOAT drsqd;                      // Distance squared
  FLOAT rcentre[ndim];              // Position of sphere centre
  FLOAT volume;                     // Volume of box
  FLOAT *r;                         // Positions of all particles

  // Create local copies of initial conditions parameters
  int Npart      = simparams->intparams["Nhydro"];
  FLOAT rhofluid = simparams->floatparams["rhofluid1"];
  FLOAT press    = simparams->floatparams["press1"];
  FLOAT radius   = simparams->floatparams["radius"];
  FLOAT gammaone = simparams->floatparams["gamma_eos"] - 1.0;
  string particle_dist = simparams->stringparams["particle_distribution"];

  debug2("[NohProblemIc::Generate]");

  r = new FLOAT[ndim*Npart];

  // Add a sphere of random particles with origin 'rcentre' and radius 'radius'
  for (k=0; k<ndim; k++) rcentre[k] = (FLOAT) 0.0;

  // Create the sphere depending on the choice of initial particle distribution
  if (particle_dist == "random") {
    Ic<ndim>::AddRandomSphere(Npart, rcentre, radius, r, sim->randnumb);
  }
  else if (particle_dist == "cubic_lattice" || particle_dist == "hexagonal_lattice") {
    Nsphere = Ic<ndim>::AddLatticeSphere(Npart, rcentre, radius, particle_dist, r, sim->randnumb);
    if (Nsphere != Npart) cout << "Warning! Unable to converge to required "
                               << "no. of ptcls due to lattice symmetry" << endl;
    Npart = Nsphere;
  }
  else {
    string message = "Invalid particle distribution option";
    ExceptionHandler::getIstance().raise(message);
  }

  // Allocate local and main particle memory
  hydro->Nhydro = Npart;
  sim->AllocateParticleMemory();

  if (ndim == 1) volume = (FLOAT) 2.0*radius;
  else if (ndim == 2) volume = pi*radius*radius;
  else if (ndim == 3) volume = (FLOAT) 4.0*onethird*pi*pow(radius,3);

  // Record particle properties in main memory
  for (i=0; i<Npart; i++) {
    Particle<ndim>& part = hydro->GetParticlePointer(i);
    for (k=0; k<ndim; k++) part.r[k] = r[ndim*i + k];
    for (k=0; k<ndim; k++) dr[k] = r[ndim*i + k];
    drsqd = DotProduct(dr,dr,ndim);
    drmag = sqrt(drsqd) + small_number;
    for (k=0; k<ndim; k++) part.v[k] = -(FLOAT) 1.0*dr[k]/drmag;
    part.m = rhofluid*volume/(FLOAT) Npart;
    part.h = hydro->h_fac*pow(part.m/rhofluid,invndim);
    part.u = press/rhofluid/gammaone;
  }

  sim->initial_h_provided = true;

  delete[] r;

  return;
}



template class NohIc<1>;
template class NohIc<2>;
template class NohIc<3>;
