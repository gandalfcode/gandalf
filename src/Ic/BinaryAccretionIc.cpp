//=================================================================================================
//  BinaryAccretionIc.cpp
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
//  BinaryAccretionIc::BinaryAccretionIc
/// Set-up Khi-type simulation initial conditions.
//=================================================================================================
template <int ndim>
BinaryAccretionIc<ndim>::BinaryAccretionIc(Simulation<ndim>* _sim, FLOAT _invndim) :
  Ic<ndim>(_sim, _invndim)
{
  // Some sanity checking to ensure correct dimensionality is used
  if (simparams->intparams["ndim"] == 1) {
    ExceptionHandler::getIstance().raise("Binary accretion sim only runs in 2D and 3D");
  }
}



//=================================================================================================
//  Khi::Generate
/// Set-up Khi-type simulation initial conditions.
//=================================================================================================
template <int ndim>
void BinaryAccretionIc<ndim>::Generate(void)
{
  // Only compile for 3-dimensional case
  //-----------------------------------------------------------------------------------------------
  if (ndim == 2 || ndim == 3) {

    int i;                               // Particle counter
    int j;                               // Aux. particle counter
    int k;                               // Dimension counter
    int Nbox1 = 0;                       // No. of particles in fluid 1
    int Nbox2 = 0;                       // No. of particles in fluid 2
    int Nlattice1[3];                    // Lattice dimensions for fluid 1
    int Nlattice2[3];                    // Lattice dimensions for fluid 2
    int Nneib;                           // Average no. of SPH neighbours
    FLOAT hfluid1 = (FLOAT) 0.0;         // Smoothing length of fluid 1
    FLOAT hsink;                         // Smoothing length of sink
    FLOAT rbinary[ndim];                 // Initial position of binary COM
    FLOAT vbinary[ndim];                 // Initial velocity of binary COM
    FLOAT rsonic;                        // Sonic radius
    FLOAT rsink;                         // Sink radius
    FLOAT volume1;                       // Volume of box1
    FLOAT volume2;                       // Volume of box2
    FLOAT *r1;                           // Positions for particles in fluid 1
    FLOAT *r2;                           // Positions for particles in fluid 2
    DomainBox<ndim> box1;                // Bounding box for fluid 1
    DomainBox<ndim> box2;                // Bounding box for fluid 2
    Nbody<ndim>* nbody = sim->nbody;     // Pointer to Nbody object
    Sinks<ndim>* sinks = sim->sinks;     // Point to Sinks object

    // Create local copies of initial conditions parameters
    int Nstar        = simparams->intparams["Nstar"];
    FLOAT abin       = simparams->floatparams["abin"];
    FLOAT ebin       = simparams->floatparams["ebin"];
    FLOAT phirot     = simparams->floatparams["phirot"];
    FLOAT thetarot   = simparams->floatparams["thetarot"];
    FLOAT psirot     = simparams->floatparams["psirot"];
    FLOAT vmachbin   = simparams->floatparams["vmachbin"];
    FLOAT m1         = simparams->floatparams["m1"];
    FLOAT m2         = simparams->floatparams["m2"];
    FLOAT gammaone   = simparams->floatparams["gamma_eos"] - (FLOAT) 1.0;
    FLOAT rhofluid1  = simparams->floatparams["rhofluid1"];
    FLOAT rhofluid2  = simparams->floatparams["rhofluid2"];
    FLOAT press1     = simparams->floatparams["press1"];
    string part_dist = simparams->stringparams["particle_distribution"];
    Nlattice1[0]     = simparams->intparams["Nlattice1[0]"];
    Nlattice1[1]     = simparams->intparams["Nlattice1[1]"];
    Nlattice1[2]     = simparams->intparams["Nlattice1[2]"];
    Nlattice2[0]     = simparams->intparams["Nlattice2[0]"];
    Nlattice2[1]     = simparams->intparams["Nlattice2[1]"];
    Nlattice2[2]     = simparams->intparams["Nlattice2[2]"];

    debug2("[Ic::BinaryAccretion]");

    // Convert parameters to dimensionless units
    rhofluid1 /= simunits.rho.inscale;
    rhofluid2 /= simunits.rho.inscale;
    press1    /= simunits.press.inscale;

    // Compute number of particles in each fluid box
    if (ndim == 2) {
      Nbox1 = Nlattice1[0]*Nlattice1[1];
      Nbox2 = Nlattice2[0]*Nlattice2[1];
    }
    else if (ndim == 3) {
      Nbox1 = Nlattice1[0]*Nlattice1[1]*Nlattice1[2];
      Nbox2 = Nlattice2[0]*Nlattice2[1]*Nlattice2[2];
    }


    // Set box limits for both fluids, depending on the selected number
    if (Nbox1 > 0 && Nbox2 == 0) {
      box1 = icBox;
    }
    else if (Nbox1 > 0 && Nbox2 > 0) {
      box1 = icBox;
      box2 = icBox;
      box1.min[0] = icBox.min[0];
      box1.max[0] = icBox.min[0] + icBox.half[0];
      box2.min[0] = icBox.min[0] + icBox.half[0];
      box2.max[0] = icBox.max[0];
    }
    else {
      string message = "Invalid number of particles chosen";
      ExceptionHandler::getIstance().raise(message);
    }


    // Compute size and range of fluid bounding boxes
    //---------------------------------------------------------------------------------------------
    if (ndim == 2) {
      volume1 = (box1.max[0] - box1.min[0])*(box1.max[1] - box1.min[1]);
      volume2 = (box2.max[0] - box2.min[0])*(box2.max[1] - box2.min[1]);
      Nneib   = (int) (pi*pow(hydro->kernp->kernrange*hydro->h_fac,2));
      hfluid1 = sqrtf((volume1*(FLOAT) Nneib)/((FLOAT) 4.0*(FLOAT) Nbox1));
    }
    else if (ndim == 3) {
      volume1 = (box1.max[0] - box1.min[0])*
        (box1.max[1] - box1.min[1])*(box1.max[2] - box1.min[2]);
      volume2 = (box2.max[0] - box2.min[0])*
        (box2.max[1] - box2.min[1])*(box2.max[2] - box2.min[2]);
      Nneib = (int) (pi*pow(hydro->kernp->kernrange*hydro->h_fac,2));
      hfluid1 = powf(((FLOAT) 3.0*volume1*(FLOAT) Nneib)/((FLOAT) 32.0*pi*(FLOAT) Nbox1),onethird);
    }


    // Allocate main particle memory
    hydro->Nhydro = Nbox1 + Nbox2;
    sim->nbody->Nstar = Nstar;
    sim->AllocateParticleMemory();


    // Add a cube of random particles defined by the simulation bounding box and
    // depending on the chosen particle distribution
    //---------------------------------------------------------------------------------------------
    if (Nbox1 > 0) {
      r1 = new FLOAT[ndim*Nbox1];
      if (part_dist == "random") {
        Ic<ndim>::AddRandomBox(Nbox1, box1, r1, sim->randnumb);
      }
      else if (part_dist == "cubic_lattice") {
        Ic<ndim>::AddCubicLattice(Nbox1, Nlattice1, box1, true, r1);
      }
      else if (part_dist == "hexagonal_lattice") {
        Ic<ndim>::AddHexagonalLattice(Nbox1, Nlattice1, box1, true, r1);
      }
      else {
        string message = "Invalid particle distribution option";
        ExceptionHandler::getIstance().raise(message);
      }

      // Record positions in main memory
      for (i=0; i<Nbox1; i++) {
        Particle<ndim>& part = hydro->GetParticlePointer(i);
        for (k=0; k<ndim; k++) part.r[k] = r1[ndim*i + k];
        part.r[0] += (FLOAT) 0.25*icBox.size[0];
        if (part.r[0] > icBox.max[0]) part.r[0] -= icBox.size[0];
        for (k=0; k<ndim; k++) part.v[k] = (FLOAT) 0.0;
        part.m = rhofluid1*volume1/(FLOAT) Nbox1;
        part.h = hydro->h_fac*pow(part.m/rhofluid1,invndim);
        part.u = press1/rhofluid1/gammaone;
      }
      delete[] r1;
    }


    // Add a cube of random particles defined by the simulation bounding box and
    // depending on the chosen particle distribution
    //---------------------------------------------------------------------------------------------
    if (Nbox2 > 0) {
      r2 = new FLOAT[ndim*Nbox2];
      if (part_dist == "random") {
        Ic<ndim>::AddRandomBox(Nbox2, box2, r2, sim->randnumb);
      }
      else if (part_dist == "cubic_lattice") {
        Ic<ndim>::AddCubicLattice(Nbox2, Nlattice2, box2, true, r2);
      }
      else if (part_dist == "hexagonal_lattice") {
        Ic<ndim>::AddHexagonalLattice(Nbox2, Nlattice2, box2, true, r2);
      }
      else {
        string message = "Invalid particle distribution option";
        ExceptionHandler::getIstance().raise(message);
      }

      // Record positions in main memory
      for (j=0; j<Nbox2; j++) {
        i = Nbox1 + j;
        Particle<ndim>& part = hydro->GetParticlePointer(i);
        for (k=0; k<ndim; k++) part.r[k] = r2[ndim*j + k];
        part.r[0] += (FLOAT) 0.25*icBox.size[0];
        if (part.r[0] > icBox.max[0]) part.r[0] -= icBox.size[0];
        for (k=0; k<ndim; k++) part.v[k] = (FLOAT) 0.0;
        part.h = hydro->h_fac*pow(part.m/rhofluid2,invndim);
        part.m = rhofluid2*volume2/(FLOAT) Nbox1;
        part.u = press1/rhofluid2/gammaone;
      }
      delete[] r2;
    }

    sim->initial_h_provided = true;

    rsonic = (FLOAT) 0.5*(m1 + m2)/(press1/rhofluid1);
    //hsink = hydro->kernp->invkernrange*rsonic;
    //hsink = min(hsink,hfluid1);
    hsink = hfluid1/pow((FLOAT) 4.4817,invndim);
    rsink = hydro->kernp->kernrange*hsink/pow((FLOAT) 4.4817,invndim);
    FLOAT mmax = (FLOAT) 2.0*rhofluid1*pow(hydro->kernp->kernrange*hfluid1,ndim)/(FLOAT) 3.0;

    cout << "Sound speed : " << sqrt(press1/rhofluid1) << endl;
    cout << "rsonic      : " << rsonic << endl;
    cout << "rbondi      : " << (FLOAT) 4.0*rsonic << endl;
    cout << "rsink       : " << rsink << endl;
    cout << "mmax        : " << mmax << endl;
    cout << "hfluid      : " << hfluid1 << endl;
    cout << "vbin        : " << vmachbin*sqrt(press1/rhofluid1) << endl;
    cout << "Bondi accretion, dmdt : "
         << (FLOAT) 4.0*pi*rhofluid1*(m1 + m2)*(m1 + m2)/pow(press1/rhofluid1,(FLOAT) 1.5) << endl;
    cout << "No. of particles      : " << Nbox1 << "   " << Nbox2 << "   " << hydro->Nhydro << endl;

    // Set hmin_sink here, since no other sinks will be formed
    hydro->hmin_sink = hsink;


    // Add star particles to simulation
    //---------------------------------------------------------------------------------------------
    if (Nstar == 1) {
      for (k=0; k<ndim; k++) nbody->stardata[0].r[k] = (FLOAT) 0.0;
      for (k=0; k<ndim; k++) nbody->stardata[0].v[k] = (FLOAT) 0.0;
      if (vmachbin < small_number) {
        nbody->stardata[0].r[0] = icBox.min[0] + (FLOAT) 0.5*icBox.size[0];
      }
      else {
        nbody->stardata[0].r[0] = icBox.min[0] + (FLOAT) 0.0625*icBox.size[0];
      }
      nbody->stardata[0].v[0] = vmachbin*hydro->eos->SoundSpeed(hydro->GetParticlePointer(0));
      nbody->stardata[0].m = m1 + m2;
      nbody->stardata[0].h = hsink;
      nbody->stardata[0].radius = rsink;
      sinks->sink[0].star = &(nbody->stardata[0]);
      sinks->sink[0].radius = rsink;
      sinks->sink[0].mmax = mmax;
      sinks->Nsink = Nstar;
    }
    else if (Nstar == 2) {
      for (k=0; k<ndim; k++) rbinary[k] = (FLOAT) 0.0;
      for (k=0; k<ndim; k++) vbinary[k] = (FLOAT) 0.0;
      if (vmachbin < small_number) {
        rbinary[0] = icBox.min[0] + (FLOAT) 0.5*icBox.size[0];
      }
      else {
        rbinary[0] = icBox.min[0] + (FLOAT) 0.0625*icBox.size[0];
      }
      vbinary[0] = vmachbin*hydro->eos->SoundSpeed(hydro->GetParticlePointer(0));
      Ic<ndim>::AddBinaryStar(abin,ebin,m1,m2,hsink,hsink,phirot,thetarot,psirot,0.0,
                              rbinary,vbinary,nbody->stardata[0],nbody->stardata[1]);
      sinks->sink[0].star = &(nbody->stardata[0]);
      sinks->sink[1].star = &(nbody->stardata[1]);
      sinks->sink[0].radius = rsink;
      sinks->sink[1].radius = rsink;
      sinks->sink[0].mmax = mmax;
      sinks->sink[1].mmax = mmax;
      sinks->Nsink = Nstar;
    }
    else {
      string message = "Invalid number of star particles";
      ExceptionHandler::getIstance().raise(message);
    }

  }
  //-----------------------------------------------------------------------------------------------

  return;
}



template class BinaryAccretionIc<1>;
template class BinaryAccretionIc<2>;
template class BinaryAccretionIc<3>;
