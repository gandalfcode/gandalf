//=================================================================================================
//  ContactDiscontinuityIc.cpp
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
//  ContactDiscontinuityIc::ContactDiscontinuityIc
/// Set-up square pressure-equilibrium test of contact discontinuities.
//=================================================================================================
template <int ndim>
ContactDiscontinuityIc<ndim>::ContactDiscontinuityIc(Simulation<ndim>* _sim, FLOAT _invndim) :
  Ic<ndim>(_sim, _invndim)
{
  // Some sanity checking to ensure correct dimensionality is used
  if (simparams->intparams["ndim"] != 2) {
    ExceptionHandler::getIstance().raise("Contact discontinuity test only runs in 2d");
  }
  if (simparams->intparams["dimensionless"] == 0) {
    ExceptionHandler::getIstance().raise("dimensionless units required");
  }
}



//=================================================================================================
//  Silcc::Generate
/// Set-up SILCC-type simulation initial conditions.
//=================================================================================================
template <int ndim>
void ContactDiscontinuityIc<ndim>::Generate(void)
{
  int i;                               // Particle counter
  int j;                               // Aux. particle counter
  int Nbox;                            // No. of particles in LHS box
  int Nlattice[3];                     // Particles per dimension for LHS lattice
  FLOAT volume;                        // Volume of box
  FLOAT *r;                            // Position vectors
  DomainBox<ndim> domainBox;           // LHS box

  // Create local copies of all parameters required to set-up problem
  FLOAT rhofluid1 = simparams->floatparams["rhofluid1"];
  FLOAT rhofluid2 = simparams->floatparams["rhofluid2"];
  FLOAT press     = simparams->floatparams["press1"];
  FLOAT gammaone  = simparams->floatparams["gamma_eos"] - (FLOAT) 1.0;
  Nlattice[0]     = simparams->intparams["Nlattice1[0]"];
  Nlattice[1]     = simparams->intparams["Nlattice1[1]"];

  debug2("[ContactDiscontinuityIc::Generate]");

  // 1D simulation
  //===============================================================================================
  if (ndim == 2) {
    volume = icBox.size[0]*icBox.size[1];
    Nbox = Nlattice[0]*Nlattice[1];
    std::cout << "Volume : " << volume << "   Nbox : " << Nbox << std::endl;

    // Allocate local and main particle memory
    hydro->Nhydro = Nbox;
    sim->AllocateParticleMemory();
    r = new FLOAT[ndim*hydro->Nhydro];

    Ic<ndim>::AddCubicLattice(Nbox, Nlattice, icBox, false, r);
    for (i=0; i<Nbox; i++) {
      Particle<ndim>& part = hydro->GetParticlePointer(i);
      for (int k=0; k<ndim; k++) part.r[k] = r[ndim*i + k];
      for (int k=0; k<ndim; k++) part.v[k] = (FLOAT) 0.0;
      if (fabs(part.r[0]) <= 0.25 && fabs(part.r[1]) <= 0.25) {
        part.m = rhofluid2*volume/(FLOAT) Nbox;
        part.u = press/rhofluid2/gammaone;
        part.h = hydro->h_fac*pow(part.m/rhofluid2,invndim);
      }
      else {
        part.m = rhofluid1*volume/(FLOAT) Nbox;
        part.u = press/rhofluid1/gammaone;
        part.h = hydro->h_fac*pow(part.m/rhofluid1,invndim);
      }
    }

    delete[] r;

  }
  //===============================================================================================


  // Set initial smoothing lengths and create initial ghost particles
  //-----------------------------------------------------------------------------------------------
  hydro->Nghost = 0;
  hydro->Ntot = hydro->Nhydro;
  for (int i=0; i<hydro->Nhydro; i++) hydro->GetParticlePointer(i).flags.set(active);

  sim->initial_h_provided = true;


  return;
}



template class ContactDiscontinuityIc<1>;
template class ContactDiscontinuityIc<2>;
template class ContactDiscontinuityIc<3>;
