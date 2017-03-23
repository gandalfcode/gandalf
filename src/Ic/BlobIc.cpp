//=================================================================================================
//  BlobIc.cpp
//  Class for generating initial conditions ...
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
//  BlobIc::BlobIc
/// Set-up Khi-type simulation initial conditions.
//=================================================================================================
template <int ndim>
BlobIc<ndim>::BlobIc(Simulation<ndim>* _sim, FLOAT _invndim) :
  Ic<ndim>(_sim, _invndim)
{
  // Some sanity checking to ensure correct dimensionality and (dimensionless) units are used
  if (simparams->intparams["ndim"] == 1) {
    ExceptionHandler::getIstance().raise("Blob test only in 2D or 3D");
  }
  if (simparams->intparams["dimensionless"] != 1) {
    ExceptionHandler::getIstance().raise("Dimensionless units required for Blob test");
  }
}



//=================================================================================================
//  BlobIc::Generate
/// Set-up Khi-type simulation initial conditions.
//=================================================================================================
template <int ndim>
void BlobIc<ndim>::Generate(void)
{
  // Create local copies of initial conditions parameters
  const FLOAT radius   = simparams->floatparams["radius"];
  const FLOAT rhofluid = simparams->floatparams["rhofluid1"];
  const FLOAT rhosphere = simparams->floatparams["rhofluid2"];
  const FLOAT gammaone = simparams->floatparams["gamma_eos"] - 1.0;
  const FLOAT gamma = simparams->floatparams["gamma_eos"];
  const FLOAT press = simparams->floatparams["press1"];
  const FLOAT mach = simparams->floatparams["mach"];
  const string particle_dist = simparams->stringparams["particle_distribution"];

  debug2("[BlobIc::Generate]");

  FLOAT volume_sphere;
  if (ndim == 1) volume_sphere = (FLOAT) 2.0*radius;
  else if (ndim == 2) volume_sphere = pi*radius*radius;
  else if (ndim == 3) volume_sphere = (FLOAT) 4.0*onethird*pi*pow(radius,3);

  // Add a sphere of random particles with origin 'rcentre' and radius 'radius'
  FLOAT rcentre[ndim];
  for (int k=0; k<ndim; k++) rcentre[k] = (FLOAT) 0.0;

  FLOAT volume_box=1;
  for (int k=0; k<ndim; k++) {
    volume_box *= (icBox.max[k]-icBox.min[k]);
  }
  // Add an uniform background
  int Nlattice[3];
  Nlattice[0]=simparams->intparams["Nlattice1[0]"];
  Nlattice[1]=simparams->intparams["Nlattice1[1]"];
  Nlattice[2]=simparams->intparams["Nlattice1[2]"];
  int Nbox=1;
  for (int k=0; k<ndim; k++) {
    Nbox *= Nlattice[k];
  }
  vector<FLOAT> r_background(ndim*Nbox);
  if (particle_dist == "random") {
    Ic<ndim>::AddRandomBox(Nbox, icBox, &r_background[0], sim->randnumb);
  }
  else if (particle_dist == "cubic_lattice") {
    Ic<ndim>::AddCubicLattice(Nbox, Nlattice, icBox, true, &r_background[0]);
  }
  else if (particle_dist == "hexagonal_lattice") {
    Ic<ndim>::AddHexagonalLattice(Nbox, Nlattice, icBox, true, &r_background[0]);
  }
  else {
    string message = "Invalid particle distribution option";
    ExceptionHandler::getIstance().raise(message);
  }
  // Count how many particles are NOT inside the sphere
  vector<FLOAT> r_back_accepted;
  r_back_accepted.reserve(ndim*Nbox);
  for (int i=0; i< Nbox; i++) {
    FLOAT distance=0;
    for (int k=0; k<ndim; k++) {
      distance += r_background[ndim*i+k]*r_background[ndim*i+k];
    }
    distance = sqrt(distance);
    if (distance>radius) {
      for (int k=0; k<ndim; k++)
        r_back_accepted.push_back(r_background[ndim*i+k]);
    }
  }
  Nbox = r_back_accepted.size()/ndim;

  const FLOAT mass_box = rhofluid*(volume_box-volume_sphere);
  const FLOAT mpart = mass_box/Nbox;


  int Nsphere = rhosphere*volume_sphere/mpart;
  vector<FLOAT> r(ndim*Nsphere);
  // Create the sphere depending on the choice of initial particle distribution
  if (particle_dist == "random") {
    Ic<ndim>::AddRandomSphere(Nsphere, rcentre, radius, &r[0], sim->randnumb);
  }
  else if (particle_dist == "cubic_lattice" || particle_dist == "hexagonal_lattice") {
    Nsphere = Ic<ndim>::AddLatticeSphere(Nsphere, rcentre, radius, particle_dist, &r[0], sim->randnumb);
//    if (Nsphere != Npart) cout << "Warning! Unable to converge to required "
//                               << "no. of ptcls due to lattice symmetry" << endl;
  }
  else {
    string message = "Invalid particle distribution option";
    ExceptionHandler::getIstance().raise(message);
  }



  // Allocate local and main particle memory
  int Npart = Nsphere+Nbox;
  hydro->Nhydro = Npart;
  sim->AllocateParticleMemory();

  cout << "Box: " << Nbox << "particles" << endl;
  cout << "Sphere: " << Nsphere << "particles" << endl;
  cout << "Total: " << Npart << "particles" << endl;


  // Record particle properties in main memory
  for (int i=0; i<Npart; i++) {
    Particle<ndim>& part = hydro->GetParticlePointer(i);

    part.m = mpart;
    if (i < Nsphere) {
      for (int k=0; k<ndim; k++) part.r[k] = r[ndim*i + k];
      part.rho = rhosphere;
    }
    else {
      for (int k=0; k<ndim; k++) part.r[k] = r_back_accepted[ndim*(i-Nsphere) + k];
      part.rho = rhofluid;
    }

    part.h = hydro->h_fac*pow(part.m/part.rho,invndim);
    // IC are in pressure equilibrium
    part.u = press/part.rho/gammaone;
    if (i >= Nsphere) {
      const FLOAT sound = sqrt(gamma*gammaone*part.u);
      const FLOAT v_back = mach*sound;
      part.v[0] = v_back;
    }
  }

  sim->initial_h_provided = true;

  return;
}



template class BlobIc<1>;
template class BlobIc<2>;
template class BlobIc<3>;
