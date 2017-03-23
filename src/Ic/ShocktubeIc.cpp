//=================================================================================================
//  ShocktubeIc.cpp
//  Class for generating initial conditions for Riemann shock-tube problems.
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
//  ShocktubeIc::ShocktubeIc
/// Set-up Shocktube tests from Riemann-problem ICs.
//=================================================================================================
template <int ndim>
ShocktubeIc<ndim>::ShocktubeIc(Simulation<ndim>* _sim, FLOAT _invndim) :
  Ic<ndim>(_sim, _invndim)
{
  // Some sanity checking to ensure correct dimensionality is used
  if (simparams->intparams["ndim"] != 1) {
    ExceptionHandler::getIstance().raise("Currently can only set-up 1D shocktubes");
  }
  if (simparams->intparams["dimensionless"] != 1) {
    ExceptionHandler::getIstance().raise("dimensionless units not used");
  }
}



//=================================================================================================
//  Khi::Generate
/// Set-up Khi-type simulation initial conditions.
//=================================================================================================
template <int ndim>
void ShocktubeIc<ndim>::Generate(void)
{
  // Only compile for 3-dimensional case
  //-----------------------------------------------------------------------------------------------
  if (ndim == 1) {

    int i;                                 // Particle counter
    int j;                                 // Aux. particle counter
    int k;                                 // Dimension counter
    int Nbox1;                             // No. of particles in LHS box
    int Nbox2;                             // No. of particles in RHS box
    int Nlattice1[ndim];                   // Particles per dimension for LHS lattice
    int Nlattice2[ndim];                   // Particles per dimension for RHS lattice
    FLOAT volume;                          // Volume of box
    FLOAT vfluid1[ndim];                   // Velocity vector of LHS fluid
    FLOAT vfluid2[ndim];                   // Velocity vector of RHS fluid
    FLOAT *r;                              // Position vectors
    DomainBox<ndim> box1;                  // LHS box
    DomainBox<ndim> box2;                  // RHS box
    Parameters* simparams = sim->simparams;

    // Set local copies of various input parameters for setting-up test
    FLOAT rhofluid1 = simparams->floatparams["rhofluid1"];
    FLOAT rhofluid2 = simparams->floatparams["rhofluid2"];
    FLOAT press1    = simparams->floatparams["press1"];
    FLOAT press2    = simparams->floatparams["press2"];
    FLOAT temp0     = simparams->floatparams["temp0"];
    FLOAT mu_bar    = simparams->floatparams["mu_bar"];
    FLOAT gammaone  = simparams->floatparams["gamma_eos"] - (FLOAT) 1.0;
    Nlattice1[0]    = simparams->intparams["Nlattice1[0]"];
    Nlattice2[0]    = simparams->intparams["Nlattice2[0]"];
    vfluid1[0]      = simparams->floatparams["vfluid1[0]"];
    vfluid2[0]      = simparams->floatparams["vfluid2[0]"];

    debug2("[ShocktubeIc::Generate]");

    // Compute size and range of fluid bounding boxes
    box1.min[0]   = icBox.min[0];
    box1.max[0]   = (FLOAT) 0.0;
    box2.min[0]   = (FLOAT) 0.0;
    box2.max[0]   = icBox.max[0];
    Nbox1         = Nlattice1[0];
    Nbox2         = Nlattice2[0];
    hydro->Nhydro = Nbox1 + Nbox2;

    bool dusty_shock = simparams->stringparams["dust_forces"] != "none" ;
    if (dusty_shock) hydro->Nhydro *= 2;

    // Allocate local and main particle memory
    sim->AllocateParticleMemory();
    r = new FLOAT[ndim*hydro->Nhydro];


    // Add particles for LHS of the shocktube
    //---------------------------------------------------------------------------------------------
    if (Nbox1 > 0) {
      Ic<ndim>::AddCubicLattice(Nbox1, Nlattice1, box1, false, r);
      volume = box1.max[0] - box1.min[0];

      for (i=0; i<Nbox1; i++) {
        Particle<ndim>& part = hydro->GetParticlePointer(i);
        for (k=0; k<ndim; k++) part.r[k] = r[ndim*i + k];
        for (k=0; k<ndim; k++) part.v[k] = (FLOAT) 0.0;
        part.v[0] = vfluid1[0];
        part.rho = rhofluid1;
        part.m = rhofluid1*volume/(FLOAT) Nbox1;
        part.h = hydro->h_fac*pow(part.m/rhofluid1,invndim);
        if (hydro->gas_eos == "isothermal") part.u = temp0/gammaone/mu_bar;
        else part.u = press1/rhofluid1/gammaone;
      }
    }

    // Add particles for RHS of the shocktube
    //---------------------------------------------------------------------------------------------
    if (Nbox2 > 0) {
      Ic<ndim>::AddCubicLattice(Nbox2, Nlattice2, box2, false, r);
      volume = box2.max[0] - box2.min[0];

      for (j=0; j<Nbox2; j++) {
        i = Nbox1 + j;
        Particle<ndim>& part = hydro->GetParticlePointer(i);
        for (k=0; k<ndim; k++) part.r[k] = r[ndim*j + k];
        for (k=0; k<ndim; k++) part.v[k] = (FLOAT) 0.0;
        part.rho = rhofluid2;
        part.v[0] = vfluid2[0];
        part.m = rhofluid2*volume/(FLOAT) Nbox2;
        part.h = hydro->h_fac*pow(part.m/rhofluid2,invndim);
        if (hydro->gas_eos == "isothermal") part.u = temp0/gammaone/mu_bar;
        else part.u = press2/rhofluid2/gammaone;
      }
    }


    // Add a slightly offset dust lattice
    if (dusty_shock){
      int Ngas = Nbox1 + Nbox2 ;
      FLOAT d2g = simparams->floatparams["dust_mass_factor"] ;
      for (j = 0; j < Ngas; ++j){
        Particle<ndim>& pg = hydro->GetParticlePointer(j) ;
        Particle<ndim>& pd = hydro->GetParticlePointer(j+Ngas) ;
        pd = pg ;
        pg.ptype = gas_type ;
        pd.ptype = dust_type ;
        pd.r[0] += 0.01 * pd.h ;
        pd.m *= d2g ;
        pd.u = 0 ;
        pd.h_dust = pd.h ;
      }
    }

    sim->initial_h_provided = true;

    delete[] r;

  }
  //-----------------------------------------------------------------------------------------------

  return;
}



template class ShocktubeIc<1>;
template class ShocktubeIc<2>;
template class ShocktubeIc<3>;
