//=================================================================================================
//  UniformIc.cpp
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
//  UniformIc::UniformIc
/// Set-up SILCC-type simulation initial conditions.
//=================================================================================================
template <int ndim>
UniformIc<ndim>::UniformIc(Simulation<ndim>* _sim, FLOAT _invndim) :
  Ic<ndim>(_sim, _invndim)
{
}



//=================================================================================================
//  Silcc::Generate
/// Set-up SILCC-type simulation initial conditions.
//=================================================================================================
template <int ndim>
void UniformIc<ndim>::Generate(void)
{
  string ic = simparams->stringparams["ic"];

  debug2("[UniformIc::Generate");

  //-----------------------------------------------------------------------------------------------
  if (ic == "box") {
    if (simparams->intparams["dimensionless"] == 0) {
      ExceptionHandler::getIstance().raise("dimensionless units required");
    }
    int i,k;                          // Particle and dimension counters
    int Nbox;                         // No. of particles in box
    int Nlattice[3];                  // Particles per dimension for LHS lattice
    FLOAT volume;                     // Volume of box
    FLOAT *r = 0;                     // Position vectors of all particles

    DomainBox<ndim>& icBox = sim->icBox;

    // Local copy of important parameters
    string particle_dist = simparams->stringparams["particle_distribution"];
    int Npart = simparams->intparams["Nhydro"];
    //FLOAT rhobox = simparams->intparams["rhofluid1"];
    Nlattice[0] = simparams->intparams["Nlattice1[0]"];
    Nlattice[1] = simparams->intparams["Nlattice1[1]"];
    Nlattice[2] = simparams->intparams["Nlattice1[2]"];

    // Compute volume and number of particles inside box
    if (ndim == 1) {
      volume = icBox.max[0] - icBox.min[0];
      Nbox = Nlattice[0];
    }
    else if (ndim == 2) {
      volume = (icBox.max[0] - icBox.min[0])*(icBox.max[1] - icBox.min[1]);
      Nbox = Nlattice[0]*Nlattice[1];
    }
    else if (ndim == 3) {
      volume = (icBox.max[0] - icBox.min[0])*
        (icBox.max[1] - icBox.min[1])*(icBox.max[2] - icBox.min[2]);
      Nbox = Nlattice[0]*Nlattice[1]*Nlattice[2];
    }

    // Add a cube of random particles defined by the simulation bounding box and
    // depending on the chosen particle distribution
    if (particle_dist == "random") {
      r = new FLOAT[ndim*Npart];
      Ic<ndim>::AddRandomBox(Npart, icBox, r, sim->randnumb);
    }
    else if (particle_dist == "cubic_lattice") {
      Npart = Nbox;
      r = new FLOAT[ndim*Npart];
      Ic<ndim>::AddCubicLattice(Npart, Nlattice, icBox, true, r);
    }
    else if (particle_dist == "hexagonal_lattice") {
      Npart = Nbox;
      r = new FLOAT[ndim*Npart];
      Ic<ndim>::AddHexagonalLattice(Npart, Nlattice, icBox, true, r);
    }
    else {
      string message = "Invalid particle distribution option";
      ExceptionHandler::getIstance().raise(message);
    }

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
      part.m = volume/ (FLOAT) hydro->Nhydro;
      part.h = hydro->h_fac*pow(volume / (FLOAT) hydro->Nhydro,invndim);
      part.u = (FLOAT) 1.5;
      part.iorig = i;
    }

    sim->initial_h_provided = true;

    delete[] r;
  }
  //-----------------------------------------------------------------------------------------------
  else if (ic == "sphere") {

    int i,k;                             // Particle and dimension counters
    int Nsphere;                       // Actual number of particles in sphere
    FLOAT rcentre[ndim];                 // Position of sphere centre
    FLOAT rhofluid;                      // Density of fluid
    FLOAT volume;                        // Volume of sphere
    FLOAT *r;                            // Particle position vectors

    // Local copies of important parameters
    int Npart      = simparams->intparams["Nhydro"];
    FLOAT mcloud   = simparams->floatparams["mcloud"];
    FLOAT radius   = simparams->floatparams["radius"];
    FLOAT press    = simparams->floatparams["press1"];
    FLOAT gammaone = simparams->floatparams["gamma_eos"] - 1.0;
    string particle_dist = simparams->stringparams["particle_distribution"];

    mcloud /= simunits.m.outscale;
    radius /= simunits.r.outscale;
    press  /= simunits.press.outscale;

    r = new FLOAT[ndim*Npart];
    for (i=0; i<ndim*Npart; i++) r[i] = (FLOAT) 0.0;

    // Add a sphere of random particles with origin 'rcentre' and radius 'radius'
    for (k=0; k<ndim; k++) rcentre[k] = (FLOAT) 0.0;

    // Create the sphere depending on the choice of initial particle distribution
    if (particle_dist == "random") {
      Ic<ndim>::AddRandomSphere(Npart, rcentre, radius, r, sim->randnumb);
    }
    else if (particle_dist == "cubic_lattice" || particle_dist == "hexagonal_lattice") {
      Nsphere = Ic<ndim>::AddLatticeSphere(Npart, rcentre, radius, particle_dist, r, sim->randnumb);
      if (Nsphere != Npart) {
        cout << "Warning! Unable to converge to required "
             << "no. of ptcls due to lattice symmetry" << endl;
      }
      Npart = Nsphere;
    }
    else {
      string message = "Invalid particle distribution option";
      ExceptionHandler::getIstance().raise(message);
    }

    hydro->Nhydro = Npart;
    sim->AllocateParticleMemory();

    if (ndim == 1) volume = (FLOAT) 2.0*radius;
    else if (ndim == 2) volume = pi*radius*radius;
    else if (ndim == 3) volume = (FLOAT) 4.0*onethird*pi*pow(radius,3);
    //if (mcloud > small_number && radius > small_number)
    //  rhofluid = mcloud / volume;
    rhofluid = mcloud / volume;


    // Record particle positions and initialise all other variables
    #pragma omp parallel for default(none)\
    shared(gammaone,mcloud,Npart,press,r,rhofluid,volume) private(i,k)
    for (i=0; i<hydro->Nhydro; i++) {
      Particle<ndim>& part = hydro->GetParticlePointer(i);
      for (k=0; k<ndim; k++) {
        part.r[k] = r[ndim*i + k];
        part.v[k] = (FLOAT) 0.0;
        part.a[k] = (FLOAT) 0.0;
      }
      //part.m = rhofluid*volume / (FLOAT) Npart;
      part.m = mcloud / (FLOAT) Npart;
      part.h = hydro->h_fac*pow(part.m/rhofluid,invndim);
      part.u = press/rhofluid/gammaone;
      part.iorig = i;
    }

    sim->initial_h_provided = true;

    delete[] r;
  }
  //-----------------------------------------------------------------------------------------------

  return;
}




template class UniformIc<1>;
template class UniformIc<2>;
template class UniformIc<3>;
