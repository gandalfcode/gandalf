//=================================================================================================
//  ShearflowIc.cpp
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
//  ShearflowIc::ShearflowIc
/// Set-up SILCC-type simulation initial conditions.
//=================================================================================================
template <int ndim>
ShearflowIc<ndim>::ShearflowIc(Simulation<ndim>* _sim, FLOAT _invndim) :
  Ic<ndim>(_sim, _invndim)
{
  // Some sanity checking to ensure correct dimensionality is used
  if (simparams->intparams["ndim"] == 1) {
    ExceptionHandler::getIstance().raise("Shear flow sim only runs in 2D and 3D");
  }
  if (simparams->intparams["dimensionless"] == 0) {
    ExceptionHandler::getIstance().raise("dimensionless units not permitted");
  }
}



//=================================================================================================
//  Silcc::Generate
/// ...
//=================================================================================================
template <int ndim>
void ShearflowIc<ndim>::Generate(void)
{

  // Only compile for 3-dimensional case
  //-----------------------------------------------------------------------------------------------
  if (ndim == 2 || ndim == 3) {

    int i;                             // Particle counter
    int k;                             // Dimension counter
    int Nbox;                          // No. of particles in box
    int Nlattice1[ndim];               // Lattice size
    FLOAT lambda;                      // Wavelength if velocity perturbation
    FLOAT kwave;                       // Wavenumber
    FLOAT volume;                      // Volume of box
    FLOAT *r;                          // Positions of particles

    // Make local copies of important parameters
    Nlattice1[0]    = simparams->intparams["Nlattice1[0]"];
    Nlattice1[1]    = simparams->intparams["Nlattice1[1]"];
    FLOAT amp       = simparams->floatparams["amp"];
    FLOAT gammaone  = simparams->floatparams["gamma_eos"] - (FLOAT) 1.0;
    FLOAT press1    = simparams->floatparams["press1"];
    FLOAT rhofluid1 = simparams->floatparams["rhofluid1"];

    debug2("[ShearflowIc::Generate]");

    // Compute size and range of fluid bounding boxes
    volume = (icBox.max[0] - icBox.min[0])*(icBox.max[1] - icBox.min[1]);
    Nbox   = Nlattice1[0]*Nlattice1[1];
    lambda = icBox.max[1] - icBox.min[1];
    kwave  = twopi/lambda;

    // Allocate local and main particle memory
    hydro->Nhydro = Nbox;
    sim->AllocateParticleMemory();
    r = new FLOAT[ndim*hydro->Nhydro];

    // Add particles from cubic lattice
    if (Nbox > 0) {
      Ic<ndim>::AddCubicLattice(Nbox, Nlattice1, icBox, false, r);
      //Ic<ndim>::AddHexagonalLattice(Nbox,Nlattice1,r,icBox,false);

      for (i=0; i<Nbox; i++) {
        Particle<ndim>& part = hydro->GetParticlePointer(i);
        for (k=0; k<ndim; k++) part.r[k] = r[ndim*i + k];
        for (k=0; k<ndim; k++) part.v[k] = (FLOAT) 0.0;
        part.v[0] = amp*sin(kwave*part.r[1]);
        part.m    = rhofluid1*volume/(FLOAT) Nbox;
        part.h    = hydro->h_fac*pow(part.m/rhofluid1,invndim);
        part.u    = press1/rhofluid1/gammaone;
      }
    }

    delete[] r;
  }
  //-----------------------------------------------------------------------------------------------

  return;
}



template class ShearflowIc<1>;
template class ShearflowIc<2>;
template class ShearflowIc<3>;
