//=================================================================================================
//  SpitzerExpansionIc.cpp
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
//  SpitzerExpansionIc::SpitzerExpansionIc
/// Set-up SILCC-type simulation initial conditions.
//=================================================================================================
template <int ndim>
SpitzerExpansionIc<ndim>::SpitzerExpansionIc(Simulation<ndim>* _sim, FLOAT _invndim) :
  Ic<ndim>(_sim, _invndim)
{
  // Some sanity checking to ensure correct dimensionality is used
  if (simparams->intparams["ndim"] != 3) {
    ExceptionHandler::getIstance().raise("Spitzer expansion sim only runs in 3D");
  }
  if (simparams->intparams["dimensionless"] == 0) {
    ExceptionHandler::getIstance().raise("dimensionless units not permitted");
  }
}



//=================================================================================================
//  SpitzerExpansionIc::Generate
/// ...
//=================================================================================================
template <int ndim>
void SpitzerExpansionIc<ndim>::Generate(void)
{
  // Only compile for 3-dimensional case
  //-----------------------------------------------------------------------------------------------
  if (ndim == 3) {

    int i,k;                             // Particle and dimension counters
    int Nsphere;                         // Actual number of particles in sphere
    FLOAT rcentre[ndim];                 // Position of sphere centre
    FLOAT rhofluid;                      // ..
    FLOAT volume;                        // Volume of sphere
    FLOAT *r;                            // Particle position vectors

    // Local copies of important parameters
    int Npart            = simparams->intparams["Nhydro"];
    FLOAT mcloud         = simparams->floatparams["mcloud"];
    FLOAT radius         = simparams->floatparams["radius"];
    string particle_dist = simparams->stringparams["particle_distribution"];

    debug2("[SpitzerExpansionIc::Generate]");

    mcloud /= simunits.m.outscale;
    radius /= simunits.r.outscale;

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
    rhofluid = mcloud / volume;


    // Record particle positions and initialise all other variables
    #pragma omp parallel for default(none) shared(mcloud,Npart,r,rhofluid,volume) private(i,k)
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
      part.u = small_number; //press/rhofluid/gammaone;
      part.iorig = i;
    }

    sim->initial_h_provided = true;

    delete[] r;

  }
  //-----------------------------------------------------------------------------------------------

  return;
}




template class SpitzerExpansionIc<1>;
template class SpitzerExpansionIc<2>;
template class SpitzerExpansionIc<3>;
