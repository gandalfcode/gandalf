//=================================================================================================
//  KhiIc.cpp
//  Class for generating initial conditions for Khi-like simulations.
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
//  Khi::Khi
/// Set-up Khi-type simulation initial conditions.
//=================================================================================================
template <int ndim>
KhiIc<ndim>::KhiIc(Simulation<ndim>* _sim, FLOAT _invndim) :
  Ic<ndim>(_sim, _invndim)
{
  // Some sanity checking to ensure correct dimensionality and (dimensionless) units are used
  if (simparams->intparams["ndim"] != 2) {
    ExceptionHandler::getIstance().raise("KHI test only functions in 2D");
  }
  if (simparams->intparams["dimensionless"] != 1) {
    ExceptionHandler::getIstance().raise("r unit not set to pc");
  }
}



//=================================================================================================
//  Khi::Generate
/// Set-up Khi-type simulation initial conditions.
//=================================================================================================
template <int ndim>
void KhiIc<ndim>::Generate(void)
{
  // Only compile for 3-dimensional case
  //-----------------------------------------------------------------------------------------------
  if (ndim == 2) {

    int i;                                 // Particle counter
    int j;                                 // Aux. particle counter
    int k;                                 // Dimension counter
    int Nbox1;                             // No. of particles in fluid box 1
    int Nbox2;                             // No. of particles in fluid box 2
    int Nlattice1[ndim];                   // Lattice particles in fluid box 1
    int Nlattice2[ndim];                   // Lattice particles in fluid box 2
    FLOAT volume;                          // Volume of fluid box
    FLOAT vfluid1[ndim];                   // Velocity vector of fluid 1
    FLOAT vfluid2[ndim];                   // Velocity vector of fluid 2
    FLOAT *r;                              // Array of particle positions
    DomainBox<ndim> box1;                  // Bounding box of fluid 1
    DomainBox<ndim> box2;                  // Bounding box of fluid 2

    // Record local copies of all important parameters
    FLOAT rhofluid1 = simparams->floatparams["rhofluid1"];
    FLOAT rhofluid2 = simparams->floatparams["rhofluid2"];
    FLOAT press1    = simparams->floatparams["press1"];
    FLOAT press2    = simparams->floatparams["press2"];
    FLOAT gammaone  = simparams->floatparams["gamma_eos"] - 1.0;
    FLOAT amp       = simparams->floatparams["amp"];
    FLOAT lambda    = simparams->floatparams["lambda"];
    Nlattice1[0]    = simparams->intparams["Nlattice1[0]"];
    Nlattice1[1]    = simparams->intparams["Nlattice1[1]"];
    Nlattice2[0]    = simparams->intparams["Nlattice2[0]"];
    Nlattice2[1]    = simparams->intparams["Nlattice2[1]"];
    vfluid1[0]      = simparams->floatparams["vfluid1[0]"];
    vfluid2[0]      = simparams->floatparams["vfluid2[0]"];

    debug2("[KhiIc::Generate]");


    // Compute size and range of fluid bounding boxes
    //---------------------------------------------------------------------------------------------
    box1.min[0] = icBox.min[0];
    box1.max[0] = icBox.max[0];
    box1.min[1] = icBox.min[1];
    box1.max[1] = icBox.min[1] + icBox.half[1];
    box2.min[0] = icBox.min[0];
    box2.max[0] = icBox.max[0];
    box2.min[1] = icBox.min[1] + icBox.half[1];
    box2.max[1] = icBox.max[1];
    volume = (box1.max[0] - box1.min[0])*(box1.max[1] - box1.min[1]);
    Nbox1 = Nlattice1[0]*Nlattice1[1];
    Nbox2 = Nlattice2[0]*Nlattice2[1];

    // Allocate local and main particle memory
    hydro->Nhydro = Nbox1 + Nbox2;
    sim->AllocateParticleMemory();
    r = new FLOAT[ndim*hydro->Nhydro];


    // Add particles for LHS of the shocktube
    //---------------------------------------------------------------------------------------------
    if (Nbox1 > 0) {
      Ic<ndim>::AddCubicLattice(Nbox1, Nlattice1, box1, false, r);

      for (i=0; i<Nbox1; i++) {
        Particle<ndim>& part = hydro->GetParticlePointer(i);
        for (k=0; k<ndim; k++) part.r[k] = r[ndim*i + k];
        for (k=0; k<ndim; k++) part.v[k] = 0.0;
        part.r[1] -= (FLOAT) 0.25*icBox.size[1];
        if (part.r[1] < icBox.min[1]) part.r[1] += icBox.size[1];
        part.v[0] = vfluid1[0];
        part.m = rhofluid1*volume/(FLOAT) Nbox1;
        part.h = hydro->h_fac*pow(part.m/rhofluid1,invndim);
        part.u = press1/rhofluid1/gammaone;
      }
    }

    // Add particles for RHS of the shocktube
    //-----------------------------------------------------------------------------------------------
    if (Nbox2 > 0) {
      Ic<ndim>::AddCubicLattice(Nbox2, Nlattice2, box2, false, r);

      for (j=0; j<Nbox2; j++) {
        i = Nbox1 + j;
        Particle<ndim>& part = hydro->GetParticlePointer(i);
        for (k=0; k<ndim; k++) part.r[k] = r[ndim*j + k];
        for (k=0; k<ndim; k++) part.v[k] = 0.0;
        part.r[1] -= (FLOAT) 0.25*icBox.size[1];
        if (part.r[1] < icBox.min[1]) part.r[1] += icBox.size[1];
        part.v[0] = vfluid2[0];
        part.m = rhofluid2*volume/(FLOAT) Nbox2;
        part.h = hydro->h_fac*pow(part.m/rhofluid2,invndim);
        part.u = press2/rhofluid2/gammaone;
      }
    }

    // Add velocity perturbation here
    //---------------------------------------------------------------------------------------------
    FLOAT sigmapert = (FLOAT) 0.05/sqrt((FLOAT) 2.0);
    for (i=0; i<hydro->Nhydro; i++) {
      Particle<ndim>& part = hydro->GetParticlePointer(i);
      part.v[1] = amp*sin((FLOAT) 2.0*pi*part.r[0]/lambda)*
        (exp(-pow(part.r[1] + (FLOAT) 0.25,2)/(FLOAT) 2.0/sigmapert/sigmapert) +
         exp(-pow(part.r[1] - (FLOAT) 0.25,2)/(FLOAT) 2.0/sigmapert/sigmapert));
    }

    // Set initial smoothing lengths and create initial ghost particles
    //---------------------------------------------------------------------------------------------
    hydro->Nghost = 0;
    hydro->Ntot = hydro->Nhydro;
    for (i=0; i<hydro->Nhydro; i++) hydro->GetParticlePointer(i).flags.set(active);

    sim->initial_h_provided = true;

    delete[] r;

  }
  //-----------------------------------------------------------------------------------------------

  return;
}



template class KhiIc<1>;
template class KhiIc<2>;
template class KhiIc<3>;
