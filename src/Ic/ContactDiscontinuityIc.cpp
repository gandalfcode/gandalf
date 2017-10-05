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
/// Set-up SILCC-type simulation initial conditions.
//=================================================================================================
template <int ndim>
ContactDiscontinuityIc<ndim>::ContactDiscontinuityIc(Simulation<ndim>* _sim, FLOAT _invndim) :
  Ic<ndim>(_sim, _invndim)
{
  // Some sanity checking to ensure correct dimensionality is used
  if (simparams->intparams["ndim"] == 3) {
    ExceptionHandler::getIstance().raise("Contact discontinuity test only runs in 1d and 2d");
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
  int Nbox1;                           // No. of particles in LHS box
  int Nbox2;                           // No. of particles in RHS box
  int Nlattice1[3];                    // Particles per dimension for LHS lattice
  int Nlattice2[3];                    // Particles per dimension for RHS lattice
  FLOAT volume;                        // Volume of box
  FLOAT *r;                            // Position vectors
  DomainBox<ndim> box1;                // LHS box
  DomainBox<ndim> box2;                // RHS box

  // Create local copies of all parameters required to set-up problem
  FLOAT rhofluid1 = simparams->floatparams["rhofluid1"];
  FLOAT rhofluid2 = simparams->floatparams["rhofluid2"];
  FLOAT press1    = simparams->floatparams["press1"];
  FLOAT press2    = simparams->floatparams["press2"];
  FLOAT temp0     = simparams->floatparams["temp0"];
  FLOAT mu_bar    = simparams->floatparams["mu_bar"];
  FLOAT gammaone  = simparams->floatparams["gamma_eos"] - (FLOAT) 1.0;
  Nlattice1[0]    = simparams->intparams["Nlattice1[0]"];
  Nlattice1[1]    = simparams->intparams["Nlattice1[1]"];
  Nlattice2[0]    = simparams->intparams["Nlattice2[0]"];
  Nlattice2[1]    = simparams->intparams["Nlattice2[1]"];

  debug2("[ContactDiscontinuityIc::Generate]");

  // 1D simulation
  //===============================================================================================
  if (ndim == 1) {
    box1.min[0] = icBox.min[0];
    box1.max[0] = (FLOAT) 0.8*icBox.max[0];
    box2.min[0] = (FLOAT) 0.8*icBox.max[0];
    box2.max[0] = icBox.max[0];
    volume = box1.max[0] - box1.min[0];
    Nbox1 = Nlattice1[0];
    Nbox2 = Nlattice2[0];

    // Allocate local and main particle memory
    hydro->Nhydro = Nbox1 + Nbox2;
    sim->AllocateParticleMemory();
    r = new FLOAT[ndim*hydro->Nhydro];

    //---------------------------------------------------------------------------------------------
    if (Nbox1 > 0) {
      Ic<ndim>::AddCubicLattice(Nbox1, Nlattice1, box1, false, r);
      volume = box1.max[0] - box1.min[0];
      for (i=0; i<Nbox1; i++) {
        Particle<ndim>& part = hydro->GetParticlePointer(i);
        part.r[0] = r[i] - (FLOAT) 0.4*icBox.size[0];
        if (part.r[0] < icBox.min[0]) part.r[0] += icBox.size[0];
        part.v[0] = (FLOAT) 0.0;
        part.m = rhofluid1*volume/(FLOAT) Nbox1;
        part.h = hydro->h_fac*pow(part.m/rhofluid1,invndim);
        if (hydro->gas_eos == "isothermal") {
          part.u = temp0/gammaone/mu_bar;
        }
        else {
          part.u = press1/rhofluid1/gammaone;
        }
      }
    }

    //---------------------------------------------------------------------------------------------
    if (Nbox2 > 0) {
      Ic<ndim>::AddCubicLattice(Nbox2, Nlattice2, box2, false, r);
      volume = box2.max[0] - box2.min[0];
      for (j=0; j<Nbox2; j++) {
        i = Nbox1 + j;
        Particle<ndim>& part = hydro->GetParticlePointer(i);
        part.r[0] = r[j] - (FLOAT) 0.4*icBox.size[0];
        if (part.r[0] < icBox.min[0]) part.r[0] += icBox.size[0];
        part.v[0] = (FLOAT) 0.0;
        part.m = rhofluid2*volume/(FLOAT) Nbox2;
        part.h = hydro->h_fac*pow(part.m/rhofluid2,invndim);
        if (hydro->gas_eos == "isothermal") {
          part.u = temp0/gammaone/mu_bar;
        }
        else {
          part.u = press2/rhofluid2/gammaone;
        }
      }
    }

    delete[] r;

  }
  //===============================================================================================
  else if (ndim == 2) {



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
