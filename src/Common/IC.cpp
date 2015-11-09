//=================================================================================================
//  IC.cpp
//  Implementation of the initial conditions generation
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
#include "IC.h"
using namespace std;


//=================================================================================================
//  Ic::CheckInitialConditions
/// Performs some simple sanity checks on all initial conditions
//=================================================================================================
template <int ndim>
void Ic<ndim>::CheckInitialConditions(void)
{
  bool okflag;                               // Flag problem with current particle
  bool valid_ic = true;                      // Valid initial conditions flag
  int i,k;                                   // Particle and dimension counter
  DomainBox<ndim>& simbox = sim->simbox;


  // Loop through all particles performing various checks
  //-----------------------------------------------------------------------------------------------
  for (i=0; i<hydro->Nhydro; i++) {
    Particle<ndim>& part = hydro->GetParticlePointer(i);
    okflag = true;

    // Check that main particle properties are valid floating point numbers
    if (part.m <= (FLOAT) 0.0 || isnan(part.m) || isinf(part.m)) okflag = false;
    if (part.u < (FLOAT) 0.0 || isnan(part.u) || isinf(part.u)) okflag = false;
    for (k=0; k<ndim; k++) {
      if (isnan(part.r[k]) || isinf(part.r[k])) okflag = false;
      if (isnan(part.v[k]) || isinf(part.v[k])) okflag = false;
    }

    // Check that all particles reside inside any defined boundaries
    for (k=0; k<ndim; k++) {
      if (part.r[k] < simbox.boxmin[k]) {
        if (simbox.boundary_lhs[k] == periodicBoundary) okflag = false;
      }
      if (part.r[k] > simbox.boxmax[k]) {
        if (simbox.boundary_rhs[k] == periodicBoundary) okflag = false;
      }
    }

    // If flag indicates a problem, print error and quit
    if (!okflag) {
      cout << "Particle " << i << " not inside periodic box" << endl;
      for (k=0; k<ndim; k++)
        cout << "r[" << k << "] : " << part.r[k] << "    "
             << simbox.boxmin[k] << "    " << simbox.boxmax[k] << endl;
    }

    valid_ic = okflag;

  }
  //-----------------------------------------------------------------------------------------------

  if (!valid_ic) {
    string message = "Invalid initial conditions for SPH particles";
    ExceptionHandler::getIstance().raise(message);
  }

  return;
}



//=================================================================================================
//  Ic::BinaryAccretion
/// Create initial conditions for binary accretion simulation.
//=================================================================================================
template <int ndim>
void Ic<ndim>::BinaryAccretion(void)
{
  // Only compile for 2 or 3 dimensions
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
      box1 = simbox;
    }
    else if (Nbox1 > 0 && Nbox2 > 0) {
      box1 = simbox;
      box2 = simbox;
      box1.boxmin[0] = simbox.boxmin[0];
      box1.boxmax[0] = simbox.boxmin[0] + simbox.boxhalf[0];
      box2.boxmin[0] = simbox.boxmin[0] + simbox.boxhalf[0];
      box2.boxmax[0] = simbox.boxmax[0];
    }
    else {
      string message = "Invalid number of particles chosen";
      ExceptionHandler::getIstance().raise(message);
    }


    // Compute size and range of fluid bounding boxes
    //---------------------------------------------------------------------------------------------
    if (ndim == 2) {
      volume1 = (box1.boxmax[0] - box1.boxmin[0])*(box1.boxmax[1] - box1.boxmin[1]);
      volume2 = (box2.boxmax[0] - box2.boxmin[0])*(box2.boxmax[1] - box2.boxmin[1]);
      Nneib   = (int) (pi*pow(hydro->kernp->kernrange*hydro->h_fac,2));
      hfluid1 = sqrtf((volume1*(FLOAT) Nneib)/((FLOAT) 4.0*(FLOAT) Nbox1));
    }
    else if (ndim == 3) {
      volume1 = (box1.boxmax[0] - box1.boxmin[0])*
        (box1.boxmax[1] - box1.boxmin[1])*(box1.boxmax[2] - box1.boxmin[2]);
      volume2 = (box2.boxmax[0] - box2.boxmin[0])*
        (box2.boxmax[1] - box2.boxmin[1])*(box2.boxmax[2] - box2.boxmin[2]);
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
        AddRandomBox(Nbox1, box1, r1);
      }
      else if (part_dist == "cubic_lattice") {
        AddCubicLattice(Nbox1, Nlattice1, box1, true, r1);
      }
      else if (part_dist == "hexagonal_lattice") {
        AddHexagonalLattice(Nbox1, Nlattice1, box1, true, r1);
      }
      else {
        string message = "Invalid particle distribution option";
        ExceptionHandler::getIstance().raise(message);
      }

      // Record positions in main memory
      for (i=0; i<Nbox1; i++) {
        Particle<ndim>& part = hydro->GetParticlePointer(i);
        for (k=0; k<ndim; k++) part.r[k] = r1[ndim*i + k];
        part.r[0] += (FLOAT) 0.25*simbox.boxsize[0];
        if (part.r[0] > simbox.boxmax[0]) part.r[0] -= simbox.boxsize[0];
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
        AddRandomBox(Nbox2, box2, r2);
      }
      else if (part_dist == "cubic_lattice") {
        AddCubicLattice(Nbox2, Nlattice2, box2, true, r2);
      }
      else if (part_dist == "hexagonal_lattice") {
        AddHexagonalLattice(Nbox2, Nlattice2, box2, true, r2);
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
        part.r[0] += (FLOAT) 0.25*simbox.boxsize[0];
        if (part.r[0] > simbox.boxmax[0]) part.r[0] -= simbox.boxsize[0];
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
        nbody->stardata[0].r[0] = simbox.boxmin[0] + (FLOAT) 0.5*simbox.boxsize[0];
      }
      else {
        nbody->stardata[0].r[0] = simbox.boxmin[0] + (FLOAT) 0.0625*simbox.boxsize[0];
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
        rbinary[0] = simbox.boxmin[0] + (FLOAT) 0.5*simbox.boxsize[0];
      }
      else {
        rbinary[0] = simbox.boxmin[0] + (FLOAT) 0.0625*simbox.boxsize[0];
      }
      vbinary[0] = vmachbin*hydro->eos->SoundSpeed(hydro->GetParticlePointer(0));
      AddBinaryStar(abin,ebin,m1,m2,hsink,hsink,phirot,thetarot,psirot,0.0,
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

  return;
}



//=================================================================================================
//  Ic::ShockTube
/// Generate 1D shock-tube test problem.
//=================================================================================================
template <int ndim>
void Ic<ndim>::ShockTube(void)
{
  // Only compile for 1 dimension
  //-----------------------------------------------------------------------------------------------
  if (ndim == 1) {

    int i;                               // Particle counter
    int j;                               // Aux. particle counter
    int k;                               // Dimension counter
    int Nbox1;                           // No. of particles in LHS box
    int Nbox2;                           // No. of particles in RHS box
    int Nlattice1[ndim];                 // Particles per dimension for LHS lattice
    int Nlattice2[ndim];                 // Particles per dimension for RHS lattice
    //FLOAT dr[ndim];                      // Relative position vector
    //FLOAT drmag;                         // Distance
    //FLOAT drsqd;                         // Distance squared
    FLOAT volume;                        // Volume of box
    FLOAT vfluid1[ndim];                 // Velocity vector of LHS fluid
    FLOAT vfluid2[ndim];                 // Velocity vector of RHS fluid
    //FLOAT wnorm;                         // Kernel normalisation
    FLOAT *r;                            // Position vectors
    //FLOAT *uaux;                         // Temp. array for internal energy
    //FLOAT *vaux;                         // Temp. array for x-velocities
    DomainBox<ndim> box1;                // LHS box
    DomainBox<ndim> box2;                // RHS box

    Parameters* simparams = sim->simparams;
    DomainBox<ndim>& simbox = sim->simbox;

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

    debug2("[Ic::ShockTube]");

    // Compute size and range of fluid bounding boxes
    box1.boxmin[0] = simbox.boxmin[0];
    box1.boxmax[0] = (FLOAT) 0.0;
    box2.boxmin[0] = (FLOAT) 0.0;
    box2.boxmax[0] = simbox.boxmax[0];
    volume = box1.boxmax[0] - box1.boxmin[0];
    Nbox1 = Nlattice1[0];
    Nbox2 = Nlattice2[0];
    hydro->Nhydro = Nbox1 + Nbox2;

    // Allocate local and main particle memory
    sim->AllocateParticleMemory();
    r = new FLOAT[ndim*hydro->Nhydro];


    // Add particles for LHS of the shocktube
    //---------------------------------------------------------------------------------------------
    if (Nbox1 > 0) {
      AddCubicLattice(Nbox1, Nlattice1, box1, false, r);

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
      AddCubicLattice(Nbox2, Nlattice2, box2, false, r);

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

    sim->initial_h_provided = true;
    //bool smooth_ic = true;

    // Smooth the initial conditions
    //---------------------------------------------------------------------------------------------
    /*if (smooth_ic) {

      // Set initial smoothing lengths and create initial ghost particles
      //-------------------------------------------------------------------------------------------
      hydro->Nghost = 0;
      hydro->Nghostmax = hydro->Nhydromax - hydro->Nhydro;
      hydro->Ntot = hydro->Nhydro;
      for (i=0; i<hydro->Nhydro; i++) hydro->GetParticlePointer(i).active = true;

      //hydro->InitialSmoothingLengthGuess();
      sphneib->BuildTree(rebuild_tree,n,ntreebuildstep,ntreestockstep,
                         hydro->Ntot,hydro->Nhydromax,partdata,sph,timestep);

      // Search ghost particles
      sphneib->SearchBoundaryGhostParticles((FLOAT) 0.0,simbox,sph);

      // Update neighbour tree
      //hydro->InitialSmoothingLengthGuess();
      sphneib->BuildGhostTree(rebuild_tree,n,ntreebuildstep,ntreestockstep,
                              hydro->Ntot,hydro->Nhydromax,partdata,sph,timestep);

      // Calculate all SPH properties
      sphneib->UpdateAllSphProperties(hydro->Nhydro,hydro->Ntot,partdata,sph,nbody);

      LocalGhosts->CopyHydroDataToGhosts(simbox,sph);
      sphneib->BuildTree(rebuild_tree,n,ntreebuildstep,ntreestockstep,
                         hydro->Ntot,hydro->Nhydromax,partdata,sph,timestep);
      sphneib->BuildGhostTree(rebuild_tree,n,ntreebuildstep,ntreestockstep,
                              hydro->Ntot,hydro->Nhydromax,partdata,sph,timestep);

      // Calculate all SPH properties
      sphneib->UpdateAllSphProperties(hydro->Nhydro,hydro->Ntot,partdata,sph,nbody);


      uaux = new FLOAT[hydro->Nhydro];
      vaux = new FLOAT[hydro->Nhydro*ndim];

      // Now compute smoothed quantities
#pragma omp parallel for default(none) shared(uaux,vaux) private(dr,drmag,drsqd,i,j,k,wnorm)
      for (i=0; i<hydro->Nhydro; i++) {
        Particle<ndim>& part = hydro->GetParticlePointer(i);
        uaux[i] = (FLOAT) 0.0;
        for (k=0; k<ndim; k++) vaux[ndim*i + k] = (FLOAT) 0.0;
        wnorm = (FLOAT) 0.0;
        for (j=0; j<hydro->Ntot; j++) {
          Particle<ndim>& partj = hydro->GetParticlePointer(j);
          for (k=0; k<ndim; k++) dr[k] = partj.r[k] - part.r[k];
          drsqd = DotProduct(dr,dr,ndim);
          if (drsqd > pow(hydro->kernp->kernrange*part.h,2)) continue;
          drmag = sqrt(drsqd);
          uaux[i] += partj.m*partj.u*hydro->kernp->w0(drmag*part.invh)*pow(part.invh,ndim)*part.invrho;
          for (k=0; k<ndim; k++) vaux[ndim*i + k] +=
            partj.m*partj.v[k]*hydro->kernp->w0(drmag*part.invh)*pow(part.invh,ndim)*part.invrho;
          wnorm += partj.m*hydro->kernp->w0(drmag*part.invh)*pow(part.invh,ndim)*part.invrho;
        }
        uaux[i] /= wnorm;
        for (k=0; k<ndim; k++) vaux[ndim*i + k] /= wnorm;
      }

      for (i=0; i<hydro->Nhydro; i++) {
        Particle<ndim>& part = hydro->GetParticlePointer(i);
        part.u = uaux[i];
        for (k=0; k<ndim; k++) part.v[k] = vaux[ndim*i + k];
      }

      delete[] vaux;
      delete[] uaux;

    }*/

    delete[] r;

  }
  //-----------------------------------------------------------------------------------------------
  else {
    std::ostringstream message;
    message << "Invalid dimensionality chosen for Shocktube ICs : ndim = " << ndim;
    ExceptionHandler::getIstance().raise(message.str());
  }
  //-----------------------------------------------------------------------------------------------


  return;
}



//=================================================================================================
//  Ic::UniformBox
/// Populate the simulation bounding box with random particles.
//=================================================================================================
template <int ndim>
void Ic<ndim>::UniformBox(void)
{
  int i,k;                          // Particle and dimension counters
  int Nbox;                         // No. of particles in box
  int Nlattice[3];                  // Particles per dimension for LHS lattice
  FLOAT volume;                     // Volume of box
  FLOAT *r = 0;                     // Position vectors of all particles

  DomainBox<ndim>& simbox = sim->simbox;

  // Local copy of important parameters
  string particle_dist = simparams->stringparams["particle_distribution"];
  int Npart = simparams->intparams["Nhydro"];
  //FLOAT rhobox = simparams->intparams["rhofluid1"];
  Nlattice[0] = simparams->intparams["Nlattice1[0]"];
  Nlattice[1] = simparams->intparams["Nlattice1[1]"];
  Nlattice[2] = simparams->intparams["Nlattice1[2]"];

  debug2("[Ic::UniformBox]");

  // Compute volume and number of particles inside box
  if (ndim == 1) {
    volume = simbox.boxmax[0] - simbox.boxmin[0];
    Nbox = Nlattice[0];
  }
  else if (ndim == 2) {
    volume = (simbox.boxmax[0] - simbox.boxmin[0])*(simbox.boxmax[1] - simbox.boxmin[1]);
    Nbox = Nlattice[0]*Nlattice[1];
  }
  else if (ndim == 3) {
    volume = (simbox.boxmax[0] - simbox.boxmin[0])*
      (simbox.boxmax[1] - simbox.boxmin[1])*(simbox.boxmax[2] - simbox.boxmin[2]);
    Nbox = Nlattice[0]*Nlattice[1]*Nlattice[2];
  }

  // Add a cube of random particles defined by the simulation bounding box and
  // depending on the chosen particle distribution
  if (particle_dist == "random") {
    r = new FLOAT[ndim*Npart];
    AddRandomBox(Npart, simbox, r);
  }
  else if (particle_dist == "cubic_lattice") {
    Npart = Nbox;
    r = new FLOAT[ndim*Npart];
    AddCubicLattice(Npart, Nlattice, simbox, true, r);
  }
  else if (particle_dist == "hexagonal_lattice") {
    Npart = Nbox;
    r = new FLOAT[ndim*Npart];
    AddHexagonalLattice(Npart, Nlattice, simbox, true, r);
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

  return;
}



//=================================================================================================
//  Ic::UniformSphere
/// Create a uniform-density sphere of particles of given origin and radius.
//=================================================================================================
template <int ndim>
void Ic<ndim>::UniformSphere(void)
{
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

  debug2("[Ic::UniformSphere]");

  mcloud /= simunits.m.outscale;
  radius /= simunits.r.outscale;
  press  /= simunits.press.outscale;


  r = new FLOAT[ndim*Npart];
  for (i=0; i<ndim*Npart; i++) r[i] = (FLOAT) 0.0;

  // Add a sphere of random particles with origin 'rcentre' and radius 'radius'
  for (k=0; k<ndim; k++) rcentre[k] = (FLOAT) 0.0;

  // Create the sphere depending on the choice of initial particle distribution
  if (particle_dist == "random") {
    AddRandomSphere(Npart, rcentre, radius, r);
  }
  else if (particle_dist == "cubic_lattice" || particle_dist == "hexagonal_lattice") {
    Nsphere = AddLatticeSphere(Npart, rcentre, radius, particle_dist, r);
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

  return;
}



//=================================================================================================
//  Ic::ContactDiscontinuity
/// Set-up contact discontinuity problem.
//=================================================================================================
template <int ndim>
void Ic<ndim>::ContactDiscontinuity(void)
{
  int i;                               // Particle counter
  int j;                               // Aux. particle counter
  int Nbox1;                           // No. of particles in LHS box
  int Nbox2;                           // No. of particles in RHS box
  int Nlattice1[3];                    // Particles per dimension for LHS lattice
  int Nlattice2[3];                    // Particles per dimension for RHS lattice
  FLOAT volume;                        // Volume of box
  //FLOAT vfluid1[3];                    // Velocity vector of LHS fluid
  //FLOAT vfluid2[3];                    // Velocity vector of RHS fluid
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
  //vfluid1[0] = simparams->floatparams["vfluid1[0]"];
  //vfluid2[0] = simparams->floatparams["vfluid2[0]"];

  debug2("[Ic::ContactDiscontinuity]");


  // 1D simulation
  //===============================================================================================
  if (ndim == 1) {
    box1.boxmin[0] = simbox.boxmin[0];
    box1.boxmax[0] = (FLOAT) 0.8*simbox.boxmax[0];
    box2.boxmin[0] = (FLOAT) 0.8*simbox.boxmax[0];
    box2.boxmax[0] = simbox.boxmax[0];
    volume = box1.boxmax[0] - box1.boxmin[0];
    Nbox1 = Nlattice1[0];
    Nbox2 = Nlattice2[0];

    // Allocate local and main particle memory
    hydro->Nhydro = Nbox1 + Nbox2;
    sim->AllocateParticleMemory();
    r = new FLOAT[ndim*hydro->Nhydro];

    //---------------------------------------------------------------------------------------------
    if (Nbox1 > 0) {
      AddCubicLattice(Nbox1, Nlattice1, box1, false, r);
      volume = box1.boxmax[0] - box1.boxmin[0];
      for (i=0; i<Nbox1; i++) {
        Particle<ndim>& part = hydro->GetParticlePointer(i);
        part.r[0] = r[i] - (FLOAT) 0.4*simbox.boxsize[0];
        if (part.r[0] < simbox.boxmin[0]) part.r[0] += simbox.boxsize[0];
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
      AddCubicLattice(Nbox2, Nlattice2, box2, false, r);
      volume = box2.boxmax[0] - box2.boxmin[0];
      for (j=0; j<Nbox2; j++) {
        i = Nbox1 + j;
        Particle<ndim>& part = hydro->GetParticlePointer(i);
        part.r[0] = r[j] - (FLOAT) 0.4*simbox.boxsize[0];
        if (part.r[0] < simbox.boxmin[0]) part.r[0] += simbox.boxsize[0];
        part.v[0] = (FLOAT) 0.0;
        part.m = rhofluid2*volume/(FLOAT) Nbox2;
        part.h = hydro->h_fac*pow(part.m/rhofluid2,invndim);
        if (hydro->gas_eos == "isothermal")
          part.u = temp0/gammaone/mu_bar;
        else
          part.u = press2/rhofluid2/gammaone;
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
  hydro->Nghostmax = hydro->Nhydromax - hydro->Nhydro;
  hydro->Ntot = hydro->Nhydro;
  for (int i=0; i<hydro->Nhydro; i++) hydro->GetParticlePointer(i).active = true;

  sim->initial_h_provided = true;
  /*sphneib->BuildTree(rebuild_tree,n,ntreebuildstep,ntreestockstep,
                     hydro->Ntot,hydro->Nhydromax,partdata,sph,timestep);

  // Search ghost particles
  sphneib->SearchBoundaryGhostParticles(0.0,simbox,sph);

  // Update neighbour tree
  sphneib->BuildTree(rebuild_tree,n,ntreebuildstep,ntreestockstep,
                     hydro->Ntot,hydro->Nhydromax,partdata,sph,timestep);

  // Calculate all SPH properties
  sphneib->UpdateAllSphProperties(hydro->Nhydro,hydro->Ntot,partdata,sph,nbody);

  sphneib->BuildTree(rebuild_tree,n,ntreebuildstep,ntreestockstep,
                     hydro->Ntot,hydro->Nhydromax,partdata,sph,timestep);

  sphneib->UpdateAllSphProperties(hydro->Nhydro,hydro->Ntot,partdata,sph,nbody);

  LocalGhosts->CopyHydroDataToGhosts(simbox,sph);

  // Calculate all SPH properties
  sphneib->UpdateAllSphProperties(hydro->Nhydro,hydro->Ntot,partdata,sph,nbody);
  */

  return;
}



//=================================================================================================
//  Ic::GreshoVortex
/// Set-up Gresho & Chan (1990) vortex test problem.
//=================================================================================================
template <int ndim>
void Ic<ndim>::GreshoVortex(void)
{
  // Only compile for 2 dimensional case
  //-----------------------------------------------------------------------------------------------
  if (ndim == 2) {

    int i;                             // Particle counter
    int k;                             // Dimension counter
    int Nbox;                          // No. of particles in fluid box 1
    int Nlattice[ndim];                // Lattice particles in fluid box 1
    FLOAT dr_unit[ndim];               // ..
    FLOAT drmag;                       // ..
    FLOAT drsqd;                       // ..
    FLOAT press;
    FLOAT rhofluid;                    // Density of fluid
    FLOAT rotspeed;
    FLOAT volume;                      // Volume of fluid box
    FLOAT *r;                          // Array of particle positions

    // Record local copies of all important parameters
    FLOAT gammaone  = simparams->floatparams["gamma_eos"] - (FLOAT) 1.0;
    Nlattice[0] = simparams->intparams["Nlattice1[0]"];
    Nlattice[1] = simparams->intparams["Nlattice1[1]"];

    debug2("[Ic::GreshoVortex]");


    // Compute size and range of fluid bounding boxes
    //---------------------------------------------------------------------------------------------
    rhofluid = 1.0;
    volume = (simbox.boxmax[0] - simbox.boxmin[0])*(simbox.boxmax[1] - simbox.boxmin[1]);
    Nbox = Nlattice[0]*Nlattice[1];

    // Allocate local and main particle memory and calculate positions for cubic lattice
    hydro->Nhydro = Nbox;
    sim->AllocateParticleMemory();
    r = new FLOAT[ndim*hydro->Nhydro];
    AddCubicLattice(Nbox, Nlattice, simbox, false, r);

    for (i=0; i<Nbox; i++) {
      Particle<ndim>& part = hydro->GetParticlePointer(i);
      for (k=0; k<ndim; k++) part.r[k] = r[ndim*i + k];
      drsqd = DotProduct(part.r, part.r, ndim);
      drmag = sqrt(drsqd) + small_number;
      //invdrmag = 1.0/drmag;
      for (k=0; k<ndim; k++) dr_unit[k] = part.r[k]/drmag;

      // Set velocity and pressure/internal energy depending on radial position
      if (drmag < 0.2) {
        rotspeed = 5.0*drmag;
        press = 5.0 + 12.5*drsqd;
      }
      else if (drmag < 0.4) {
        rotspeed = 2.0 - 5.0*drmag;
        press = 9.0 + 12.5*drsqd - 20.0*drmag + 4.0*log(drmag/0.2);
      }
      else {
        rotspeed = 0.0;
        press = 3.0 + 4.0*log(2.0);
      }
      part.v[0] = -rotspeed*dr_unit[1];
      part.v[1] = rotspeed*dr_unit[0];

      part.m = rhofluid*volume/(FLOAT) Nbox;
      part.h = hydro->h_fac*pow(part.m/rhofluid,invndim);
      part.u = press/rhofluid/gammaone;
    }

  }
  //-----------------------------------------------------------------------------------------------
  else {
    string message = "Gresho vortex only available in 2D";
    ExceptionHandler::getIstance().raise(message);
  }
  //-----------------------------------------------------------------------------------------------


  return;
}



//=================================================================================================
//  Ic::KHI
/// Set-up 2D Kelvin-Helmholtz instability test.
//=================================================================================================
template <int ndim>
void Ic<ndim>::KHI(void)
{
  // Only compile for 2 dimensional case
  //-----------------------------------------------------------------------------------------------
  if (ndim == 2) {

    int i;                            // Particle counter
    int j;                            // Aux. particle counter
    int k;                            // Dimension counter
    int Nbox1;                        // No. of particles in fluid box 1
    int Nbox2;                        // No. of particles in fluid box 2
    int Nlattice1[ndim];              // Lattice particles in fluid box 1
    int Nlattice2[ndim];              // Lattice particles in fluid box 2
    FLOAT volume;                     // Volume of fluid box
    FLOAT vfluid1[ndim];              // Velocity vector of fluid 1
    FLOAT vfluid2[ndim];              // Velocity vector of fluid 2
    FLOAT *r;                         // Array of particle positions
    DomainBox<ndim> box1;             // Bounding box of fluid 1
    DomainBox<ndim> box2;             // Bounding box of fluid 2

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

    debug2("[Ic::KHI]");


    // Compute size and range of fluid bounding boxes
    //---------------------------------------------------------------------------------------------
    box1.boxmin[0] = simbox.boxmin[0];
    box1.boxmax[0] = simbox.boxmax[0];
    box1.boxmin[1] = simbox.boxmin[1];
    box1.boxmax[1] = simbox.boxmin[1] + simbox.boxhalf[1];
    box2.boxmin[0] = simbox.boxmin[0];
    box2.boxmax[0] = simbox.boxmax[0];
    box2.boxmin[1] = simbox.boxmin[1] + simbox.boxhalf[1];
    box2.boxmax[1] = simbox.boxmax[1];
    volume = (box1.boxmax[0] - box1.boxmin[0])*(box1.boxmax[1] - box1.boxmin[1]);
    Nbox1 = Nlattice1[0]*Nlattice1[1];
    Nbox2 = Nlattice2[0]*Nlattice2[1];


    // Allocate local and main particle memory
    hydro->Nhydro = Nbox1 + Nbox2;
    sim->AllocateParticleMemory();
    r = new FLOAT[ndim*hydro->Nhydro];


    // Add particles for LHS of the shocktube
    //---------------------------------------------------------------------------------------------
    if (Nbox1 > 0) {
      AddCubicLattice(Nbox1, Nlattice1, box1, false, r);

      for (i=0; i<Nbox1; i++) {
        Particle<ndim>& part = hydro->GetParticlePointer(i);
        for (k=0; k<ndim; k++) part.r[k] = r[ndim*i + k];
        for (k=0; k<ndim; k++) part.v[k] = 0.0;
        part.r[1] -= (FLOAT) 0.25*simbox.boxsize[1];
        if (part.r[1] < simbox.boxmin[1]) part.r[1] += simbox.boxsize[1];
        part.v[0] = vfluid1[0];
        part.m = rhofluid1*volume/(FLOAT) Nbox1;
        part.h = hydro->h_fac*pow(part.m/rhofluid1,invndim);
        part.u = press1/rhofluid1/gammaone;
      }
    }

    // Add particles for RHS of the shocktube
    //-----------------------------------------------------------------------------------------------
    if (Nbox2 > 0) {
      AddCubicLattice(Nbox2, Nlattice2, box2, false, r);

      for (j=0; j<Nbox2; j++) {
        i = Nbox1 + j;
        Particle<ndim>& part = hydro->GetParticlePointer(i);
        for (k=0; k<ndim; k++) part.r[k] = r[ndim*j + k];
        for (k=0; k<ndim; k++) part.v[k] = 0.0;
        part.r[1] -= (FLOAT) 0.25*simbox.boxsize[1];
        if (part.r[1] < simbox.boxmin[1]) part.r[1] += simbox.boxsize[1];
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
    hydro->Nghostmax = hydro->Nhydromax - hydro->Nhydro;
    hydro->Ntot = hydro->Nhydro;
    for (i=0; i<hydro->Nhydro; i++) hydro->GetParticlePointer(i).active = true;

    sim->initial_h_provided = true;

    // Update neighbour tree
    sim->rebuild_tree = true;
    /*sphneib->BuildTree(rebuild_tree,n,ntreebuildstep,ntreestockstep,
                       hydro->Ntot,hydro->Nhydromax,partdata,sph,timestep);

    // Search ghost particles
    sphneib->SearchBoundaryGhostParticles((FLOAT) 0.0,simbox,sph);
    sphneib->BuildGhostTree(rebuild_tree,n,ntreebuildstep,ntreestockstep,
                            hydro->Ntot,hydro->Nhydromax,partdata,sph,timestep);

    // Calculate all SPH properties
    sphneib->UpdateAllSphProperties(hydro->Nhydro,hydro->Ntot,partdata,sph,nbody);
    */
    //for (i=0; i<hydro->Nhydro; i++) {
    //  Particle<ndim>& part = hydro->GetParticlePointer(i);
    //  //part.u = press1/part.rho/gammaone;
    //}


    delete[] r;

  }
  //-----------------------------------------------------------------------------------------------
  else {
    string message = "Kelvin-Helmholtz instability only in 2D";
    ExceptionHandler::getIstance().raise(message);
  }
  //-----------------------------------------------------------------------------------------------

  return;
}



//=================================================================================================
//  Ic::NohProblem
/// Set-up Noh Problem initial conditions
//=================================================================================================
template <int ndim>
void Ic<ndim>::NohProblem(void)
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

  debug2("[Ic::NohProblem]");

  r = new FLOAT[ndim*Npart];

  // Add a sphere of random particles with origin 'rcentre' and radius 'radius'
  for (k=0; k<ndim; k++) rcentre[k] = (FLOAT) 0.0;

  // Create the sphere depending on the choice of initial particle distribution
  if (particle_dist == "random") {
    AddRandomSphere(Npart, rcentre, radius, r);
  }
  else if (particle_dist == "cubic_lattice" || particle_dist == "hexagonal_lattice") {
    Nsphere = AddLatticeSphere(Npart, rcentre, radius, particle_dist, r);
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



//=================================================================================================
//  Ic::RTI
/// Set-up 2D Rayleigh-Taylor instability test.
//=================================================================================================
template <int ndim>
void Ic<ndim>::RTI(void)
{
  // Only compile for 2 dimensional case
  //-----------------------------------------------------------------------------------------------
  if (ndim == 2) {

    int i;                             // Particle counter
    int j;                             // Aux. particle counter
    int k;                             // Dimension counter
    int Nbox1;                         // No. of particles in fluid box 1
    int Nbox2;                         // No. of particles in fluid box 2
    int Nlattice1[ndim];               // Lattice particles in fluid box 1
    int Nlattice2[ndim];               // Lattice particles in fluid box 2
    FLOAT rho;                         // Density
    FLOAT volume;                      // Volume of fluid box
    FLOAT *r;                          // Array of particle positions
    DomainBox<ndim> box1;              // Bounding box of fluid 1
    DomainBox<ndim> box2;              // Bounding box of fluid 2

    // Record local copies of all important parameters
    FLOAT rhofluid1 = simparams->floatparams["rhofluid1"];
    FLOAT rhofluid2 = simparams->floatparams["rhofluid2"];
    FLOAT press1    = simparams->floatparams["press1"];
    FLOAT press2    = simparams->floatparams["press2"];
    FLOAT gammaone  = simparams->floatparams["gamma_eos"] - 1.0;
    FLOAT amp       = simparams->floatparams["amp"];
    Nlattice1[0]    = simparams->intparams["Nlattice1[0]"];
    Nlattice1[1]    = simparams->intparams["Nlattice1[1]"];
    Nlattice2[0]    = simparams->intparams["Nlattice2[0]"];
    Nlattice2[1]    = simparams->intparams["Nlattice2[1]"];

    debug2("[Ic::RTI]");


    // Compute size and range of fluid bounding boxes
    //---------------------------------------------------------------------------------------------
    box1.boxmin[0] = simbox.boxmin[0];
    box1.boxmax[0] = simbox.boxmax[0];
    box1.boxmin[1] = simbox.boxmin[1];
    box1.boxmax[1] = simbox.boxmin[1] + simbox.boxhalf[1];
    box2.boxmin[0] = simbox.boxmin[0];
    box2.boxmax[0] = simbox.boxmax[0];
    box2.boxmin[1] = simbox.boxmin[1] + simbox.boxhalf[1];
    box2.boxmax[1] = simbox.boxmax[1];
    volume = (box1.boxmax[0] - box1.boxmin[0])*(box1.boxmax[1] - box1.boxmin[1]);
    Nbox1 = Nlattice1[0]*Nlattice1[1];
    Nbox2 = Nlattice2[0]*Nlattice2[1];


    // Allocate local and main particle memory
    hydro->Nhydro = Nbox1 + Nbox2;
    sim->AllocateParticleMemory();
    r = new FLOAT[ndim*hydro->Nhydro];


    // Add particles for LHS of the shocktube
    //---------------------------------------------------------------------------------------------
    if (Nbox1 > 0) {
      AddCubicLattice(Nbox1, Nlattice1, box1, false, r);

      for (i=0; i<Nbox1; i++) {
        Particle<ndim>& part = hydro->GetParticlePointer(i);
        for (k=0; k<ndim; k++) part.r[k] = r[ndim*i + k];
        for (k=0; k<ndim; k++) part.v[k] = 0.0;

        FLOAT delta = 0.025;
        rho = rhofluid1 + (rhofluid2 - rhofluid1)/
          ((FLOAT) 1.0 + exp(-(part.r[1] - (FLOAT) 0.5)/delta));

        part.m = rho*volume/(FLOAT) Nbox1;
        part.h = hydro->h_fac*pow(part.m/rhofluid1,invndim);
        part.u = press1/rhofluid1/gammaone;
      }
    }

    // Add particles for RHS of the shocktube
    //-----------------------------------------------------------------------------------------------
    if (Nbox2 > 0) {
      AddCubicLattice(Nbox2, Nlattice2, box2, false, r);

      for (j=0; j<Nbox2; j++) {
        i = Nbox1 + j;
        Particle<ndim>& part = hydro->GetParticlePointer(i);
        for (k=0; k<ndim; k++) part.r[k] = r[ndim*j + k];
        for (k=0; k<ndim; k++) part.v[k] = 0.0;

        FLOAT delta = 0.025;
        rho = rhofluid1 + (rhofluid2 - rhofluid1)/
          ((FLOAT) 1.0 + exp(-(part.r[1] - (FLOAT) 0.5)/delta));

        part.m = rho*volume/(FLOAT) Nbox2;
        part.h = hydro->h_fac*pow(part.m/rhofluid2,invndim);
        part.u = press2/rhofluid2/gammaone;
      }
    }

    // Add velocity perturbation here
    //---------------------------------------------------------------------------------------------
    for (i=0; i<hydro->Nhydro; i++) {
      Particle<ndim>& part = hydro->GetParticlePointer(i);
      if (part.r[1] >= 0.3 && part.r[1] <= 0.7) {
        part.v[1] = amp*((FLOAT) 1.0 + cos((FLOAT) 8.0*pi*(part.r[0] + (FLOAT) 0.25)))*
          ((FLOAT) 1.0 + cos((FLOAT) 5.0*pi*(part.r[1] - (FLOAT) 0.5)));
      }
    }

    // Set initial smoothing lengths and create initial ghost particles
    //---------------------------------------------------------------------------------------------
    hydro->Nghost = 0;
    hydro->Nghostmax = hydro->Nhydromax - hydro->Nhydro;
    hydro->Ntot = hydro->Nhydro;
    for (i=0; i<hydro->Nhydro; i++) hydro->GetParticlePointer(i).active = true;

    sim->initial_h_provided = true;

    // Update neighbour tree
    sim->rebuild_tree = true;
    /*sphneib->BuildTree(rebuild_tree,n,ntreebuildstep,ntreestockstep,
                       hydro->Ntot,hydro->Nhydromax,partdata,sph,timestep);

    // Search ghost particles
    sphneib->SearchBoundaryGhostParticles((FLOAT) 0.0,simbox,sph);
    sphneib->BuildGhostTree(rebuild_tree,n,ntreebuildstep,ntreestockstep,
                            hydro->Ntot,hydro->Nhydromax,partdata,sph,timestep);

    // Calculate all SPH properties
    sphneib->UpdateAllSphProperties(hydro->Nhydro,hydro->Ntot,partdata,sph,nbody);
    */
    /*for (i=0; i<hydro->Nhydro; i++) {
      Particle<ndim>& part = hydro->GetParticlePointer(i);
      //cout << "r[" << i << "] : " << part.r[0] << "   " << part.r[1] << endl;
      //part.u = press1/part.rho/gammaone;
    }*/


    delete[] r;

  }
  //-----------------------------------------------------------------------------------------------
  else {
    string message = "Kelvin-Helmholtz instability only in 2D";
    ExceptionHandler::getIstance().raise(message);
  }
  //-----------------------------------------------------------------------------------------------

  return;
}



//=================================================================================================
//  Ic::BossBodenheimer
/// Set-up Boss-Bodenheimer (1979) initial conditions for collapse of a
/// rotating uniform sphere with an imposed m=2 azimuthal density perturbation.
//=================================================================================================
template <int ndim>
void Ic<ndim>::BossBodenheimer(void)
{
  int i;                               // Particle counter
  int k;                               // Dimension counter
  int Nsphere;                         // Actual number of particles in sphere
  FLOAT mp;                            // Mass of one particle
  FLOAT rcentre[ndim];                 // Position of sphere centre
  FLOAT rho;                           // Fluid density
  FLOAT *r;                            // Positions of all particles
  FLOAT *v;                            // Velocities of all particles

  // Create local copies of initial conditions parameters
  int Npart      = simparams->intparams["Nhydro"];
  FLOAT amp      = simparams->floatparams["amp"];
  FLOAT angvel   = simparams->floatparams["angvel"];
  FLOAT mcloud   = simparams->floatparams["mcloud"];
  FLOAT mu_bar   = simparams->floatparams["mu_bar"];
  FLOAT press    = simparams->floatparams["press1"];
  FLOAT radius   = simparams->floatparams["radius"];
  FLOAT temp0    = simparams->floatparams["temp0"];
  FLOAT gammaone = simparams->floatparams["gamma_eos"] - 1.0;
  string particle_dist = simparams->stringparams["particle_distribution"];

  debug2("[Ic::BossBodenheimer]");

  // Convert any parameters to code units
  angvel /= simunits.angvel.outscale;
  mcloud /= simunits.m.outscale;
  press  /= simunits.press.outscale;
  radius /= simunits.r.outscale;
  temp0  /= simunits.temp.outscale;

  r = new FLOAT[ndim*Npart];
  v = new FLOAT[ndim*Npart];

  // Add a sphere of random particles with origin 'rcentre' and radius 'radius'
  for (k=0; k<ndim; k++) rcentre[k] = (FLOAT) 0.0;

  // Create the sphere depending on the choice of initial particle distribution
  if (particle_dist == "random") {
    AddRandomSphere(Npart, rcentre, radius, r);
  }
  else if (particle_dist == "cubic_lattice" || particle_dist == "hexagonal_lattice") {
    Nsphere = AddLatticeSphere(Npart, rcentre, radius, particle_dist, r);
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

  // Allocate local and main particle memory
  hydro->Nhydro = Npart;
  sim->AllocateParticleMemory();
  mp = mcloud / (FLOAT) Npart;
  rho = (FLOAT) 3.0*mcloud / ((FLOAT)4.0*pi*pow(radius,3));

  // Perturb positions of particles in cloud
  AddAzimuthalDensityPerturbation(Npart, 2, amp, rcentre, r);

  // Add solid-body rotational velocity field
  AddRotationalVelocityField(Npart, angvel, rcentre, r, v);

  // Record particle properties in main memory
  //#pragma omp parallel for default(none) shared(gammaone,Npart,mp,mu_bar,r,rho,temp0,v) private(i,k)
  for (i=0; i<Npart; i++) {
    Particle<ndim>& part = hydro->GetParticlePointer(i);
    for (k=0; k<ndim; k++) part.r[k] = r[ndim*i + k];
    for (k=0; k<ndim; k++) part.v[k] = v[ndim*i + k];
    part.m = mp;
    part.h = hydro->h_fac*pow(part.m/rho,invndim);
    part.u = temp0/gammaone/mu_bar;
    //if (hydro->gas_eos == "isothermal" || hydro->gas_eos == "barotropic")
    //  part.u = temp0/gammaone/mu_bar;
    //else
    //  part.u = press/rho/gammaone;
  }

  sim->initial_h_provided = true;

  delete[] v;
  delete[] r;


  return;
}



//=================================================================================================
//  Ic::TurbulentCore
/// Set-up Boss-Bodenheimer (1979) initial conditions for collapse of a
/// rotating uniform sphere with an imposed m=2 azimuthal density perturbation.
//=================================================================================================
template <int ndim>
void Ic<ndim>::TurbulentCore(void)
{
  int i;                               // Particle counter
  int k;                               // Dimension counter
  int Nsphere;                         // Actual number of particles in sphere
  FLOAT gpecloud;                      // Total grav. potential energy of entire cloud
  FLOAT keturb;                        // Total turbulent kinetic energy of entire cloud
  FLOAT mp;                            // Mass of one particle
  FLOAT rcentre[ndim];                 // Position of sphere centre
  FLOAT rho;                           // Fluid density
  FLOAT xmin;                          // Minimum coordinate value
  FLOAT vfactor;                       // Velocity scaling factor (to scale correct alpha_turb)
  FLOAT *r;                            // Positions of all particles
  FLOAT *v;                            // Velocities of all particles
  FLOAT dxgrid;                        // Grid spacing
  FLOAT rmax[ndim];                    // Maximum size of bounding box
  FLOAT rmin[ndim];                    // Minimum size of bounding box
  DOUBLE *vfield;                      // Table with turbulent velocity field from

  // Create local copies of initial conditions parameters
  int field_type   = simparams->intparams["field_type"];
  int gridsize     = simparams->intparams["gridsize"];
  int Npart        = simparams->intparams["Nhydro"];
  FLOAT alpha_turb = simparams->floatparams["alpha_turb"];
  FLOAT gammaone   = simparams->floatparams["gamma_eos"] - 1.0;
  FLOAT mcloud     = simparams->floatparams["mcloud"];
  FLOAT mu_bar     = simparams->floatparams["mu_bar"];
  FLOAT power_turb = simparams->floatparams["power_turb"];
  FLOAT radius     = simparams->floatparams["radius"];
  FLOAT temp0      = simparams->floatparams["temp0"];
  string particle_dist = simparams->stringparams["particle_distribution"];

#if !defined(FFTW_TURBULENCE)
  string message = "FFTW turbulence flag not set";
  ExceptionHandler::getIstance().raise(message);
#endif

  debug2("[Ic::TurbulentCore]");

  // Convert any parameters to code units
  mcloud /= simunits.m.outscale;
  radius /= simunits.r.outscale;
  temp0  /= simunits.temp.outscale;

  // Calculate gravitational potential energy of uniform density spherical cloud
  gpecloud = (FLOAT) 0.6*mcloud*mcloud/radius;

  r = new FLOAT[ndim*Npart];
  v = new FLOAT[ndim*Npart];

  // Add a sphere of random particles with origin 'rcentre' and radius 'radius'
  for (k=0; k<ndim; k++) rcentre[k] = (FLOAT) 0.0;

  // Create the sphere depending on the choice of initial particle distribution
  if (particle_dist == "random") {
    AddRandomSphere(Npart, rcentre, radius, r);
  }
  else if (particle_dist == "cubic_lattice" || particle_dist == "hexagonal_lattice") {
    Nsphere = AddLatticeSphere(Npart, rcentre, radius, particle_dist, r);
    assert(Nsphere <= Npart);
    if (Nsphere != Npart)
      cout << "Warning! Unable to converge to required "
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
  mp = mcloud / (FLOAT) Npart;
  rho = (FLOAT) 3.0*mcloud / ((FLOAT) 4.0*pi*pow(radius,3));


  // Record particle properties in main memory
  for (i=0; i<hydro->Nhydro; i++) {
    Particle<ndim>& part = hydro->GetParticlePointer(i);
    for (k=0; k<ndim; k++) part.r[k] = r[ndim*i + k];
    for (k=0; k<ndim; k++) part.v[k] = 0.0;
    part.m = mp;
    part.h = hydro->h_fac*powf(mp/rho,invndim);
    part.u = temp0/gammaone/mu_bar;
    part.itype = gas;
  }

  sim->initial_h_provided = true;


  // Generate turbulent velocity field for given power spectrum slope
  vfield = new DOUBLE[ndim*gridsize*gridsize*gridsize];

  // Calculate bounding box of SPH smoothing kernels
  for (k=0; k<ndim; k++) rmin[k] = big_number;
  for (k=0; k<ndim; k++) rmax[k] = -big_number;
  for (i=0; i<Npart; i++) {
    Particle<ndim>& part = hydro->GetParticlePointer(i);
    for (k=0; k<ndim; k++) rmin[k] = min(rmin[k], part.r[k] - hydro->kernrange*part.h);
    for (k=0; k<ndim; k++) rmax[k] = max(rmax[k], part.r[k] + hydro->kernrange*part.h);
  }

  xmin = (FLOAT) 9.9e20;
  dxgrid = (FLOAT) 0.0;
  for (k=0; k<ndim; k++) {
    xmin = min(xmin, rmin[k]);
    dxgrid = max(dxgrid, (rmax[k] - rmin[k])/(FLOAT) (gridsize - 1));
    //xmin = min(xmin,rmin[k]);
  }
  dxgrid = max(dxgrid, (FLOAT) 2.0*fabs(xmin)/(FLOAT) (gridsize - 1));

  // Generate gridded velocity field
  GenerateTurbulentVelocityField(field_type, gridsize, power_turb, vfield);

  // Now interpolate generated field onto particle positions
  InterpolateVelocityField(hydro->Nhydro, gridsize, xmin, dxgrid, r, vfield, v);

  // Finally, copy velocities to main SPH particle array
  for (i=0; i<hydro->Nhydro; i++) {
    Particle<ndim>& part = hydro->GetParticlePointer(i);
    for (k=0; k<ndim; k++) part.v[k] = v[ndim*i + k];
  }

  // Change to COM frame of reference
  sim->SetComFrame();

  // Calculate total kinetic energy of turbulent velocity field
  keturb = (FLOAT) 0.0;
  for (i=0; i<hydro->Nhydro; i++) {
    Particle<ndim>& part = hydro->GetParticlePointer(i);
    keturb += part.m*DotProduct(part.v, part.v, ndim);
  }
  keturb *= (FLOAT) 0.5;

  vfactor = sqrt(alpha_turb*gpecloud/keturb);
  cout << "Scaling factor : " << vfactor << endl;

  // Now rescale velocities to give required turbulent energy in cloud
  for (i=0; i<hydro->Nhydro; i++) {
    Particle<ndim>& part = hydro->GetParticlePointer(i);
    for (k=0; k<ndim; k++) part.v[k] *= vfactor;
  }


  delete[] v;
  delete[] r;

  return;
}



//=================================================================================================
//  Ic::BondiAccretion
/// Set-up spherically symmetric Bondi accretion test problem.  First, numerically solve Bondi's
/// differential equations about the sonic point to obtain the numerical solution.  Next stretch
/// a uniform particle distribution to the correct radial profile.  Finally add a sink of the
/// correct size and mass and simulate accretion onto the sink.
//=================================================================================================
template <int ndim>
void Ic<ndim>::BondiAccretion(void)
{
  int i;                            // Particle counter
  int k;                            // Dimension counter
  int Nsphere;                      // Actual number of particles in sphere
  FLOAT dr[ndim];                   // Relative position vector
  FLOAT drsqd;                      // Distance squared
  FLOAT mp;                         // Mass of one particle
  FLOAT racc;                       // Bondi accretion radius
  FLOAT rcentre[ndim];              // Position of sphere centre
  FLOAT rsonic;                     // Sonic radius
  FLOAT *r;                         // Positions of all particles
  FLOAT *v;                         // Velocities of all particles
  FLOAT *w,*x,*y,*z;                // Arrays for numerical Bondi solution

  Nbody<ndim>* nbody = sim->nbody;
  Sinks<ndim>* sinks = sim->sinks;

  // Create local copies of initial conditions parameters
  int Npart     = simparams->intparams["Nhydro"];
  FLOAT temp0   = simparams->floatparams["temp0"];
  FLOAT mu_bar  = simparams->floatparams["mu_bar"];
  FLOAT mcloud  = simparams->floatparams["mcloud"];
  FLOAT msink   = simparams->floatparams["m1"];
  FLOAT rsink   = simparams->floatparams["sink_radius"];
  FLOAT rhogas  = simparams->floatparams["rhofluid1"];
  string particle_dist = simparams->stringparams["particle_distribution"];

  debug2("[Ic::BondiAccretion]");

  // Convert any parameters to code units
  temp0  /= simunits.temp.outscale;
  mcloud /= simunits.m.outscale;
  msink  /= simunits.m.outscale;
  FLOAT radius = 1.0;

  FLOAT asound = sqrtf(temp0/mu_bar);
  racc = (FLOAT) 2.0*msink/asound/asound;
  rsonic = (FLOAT) 0.5*msink/asound/asound;
  FLOAT dmdt = exp((FLOAT) 1.5)*pi*msink*msink*rhogas/powf(asound,3);

  r = new FLOAT[ndim*Npart];
  v = new FLOAT[ndim*Npart];

  // Add a sphere of random particles with origin 'rcentre' and radius 'radius'
  for (k=0; k<ndim; k++) rcentre[k] = (FLOAT) 0.0;

  // Create the sphere depending on the choice of initial particle distribution
  if (particle_dist == "random") {
    AddRandomSphere(Npart, rcentre, radius, r);
  }
  else if (particle_dist == "cubic_lattice" || particle_dist == "hexagonal_lattice") {
    Nsphere = AddLatticeSphere(Npart, rcentre, radius, particle_dist, r);
    if (Nsphere != Npart)
      cout << "Warning! Unable to converge to required "
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
  //mp = mcloud / (FLOAT) Npart;
  //mp = 4.0*pi*powf(rsonic,3)*rhogas*mcloud / (FLOAT) Npart;

  // Find Bondi solution
  int Ntablemax = 4000;
  w = new FLOAT[Ntablemax];
  x = new FLOAT[Ntablemax];
  y = new FLOAT[Ntablemax];
  z = new FLOAT[Ntablemax];
  ComputeBondiSolution(Ntablemax, w, x, y, z);

  // First check that mass of cloud is not greater than that provided by the table
  if (mcloud > z[Ntablemax-1]) {
    cout << "Cloud mass too big" << endl;
    exit(0);
  }

  // Find radius of cloud
  int iradius;
  for (i=0; i<Ntablemax; i++) {
    iradius = i;
    if (z[i] >= mcloud) break;
  }
  FLOAT remainder = (mcloud - z[iradius-1])/(z[iradius] - z[iradius - 1]);
  radius = x[iradius - 1] + remainder*(x[iradius] - x[iradius - 1]);
  cout << "Cloud mass : " << mcloud        << "    iradius   : " << iradius << endl;
  cout << "Radius     : " << radius*rsonic << "    remainder : " << remainder << endl;
  cout << "rsonic     : " << rsonic        << "    racc      : " << racc << endl;
  cout << "temp0      : " << temp0         << "    asound    : " << asound << endl;
  cout << "dmdt       : " << dmdt << endl;


  // Now stretch particles in radial direction and rescale masses to reproduce the required
  // density distribution provided by the numiercal solution.
  //-----------------------------------------------------------------------------------------------
  //mp = 4.0_PR*pi*pow(rsonic,3)*rhogas*mcloud / (FLOAT) Npart;
  //mp = 0.5_PR*PI*rhogas*mcloud*msink**3 / (agas**6) / real(pgas,PR)
  //mp = mcloud/(FLOAT) Npart;
  mp = 4.0*pi*powf(rsonic,3)*rhogas*mcloud / (FLOAT) Npart;
  cout << "mp : " << mp << "    "
       << 0.5*pi*rhogas*mcloud*powf(msink,3)/powf(asound,6)/(FLOAT) Npart << endl;


  for (i=0; i<Npart; i++) {
    Particle<ndim>& part = hydro->GetParticlePointer(i);
    for (k=0; k<ndim; k++) dr[k] = r[ndim*i + k] - rcentre[k];
    FLOAT drmag = sqrtf(DotProduct(dr, dr, ndim)) + small_number;
    FLOAT mint = mcloud*powf(drmag,3);
    for (int j=0; j<Ntablemax; j++) {
      iradius = j;
      if (z[j] >= mint) break;
    }
    remainder = (mint - z[iradius-1])/(z[iradius] - z[iradius - 1]);
    FLOAT radp = x[iradius - 1] + remainder*(x[iradius] - x[iradius - 1]);
    FLOAT vradp = w[iradius - 1] + remainder*(w[iradius] - w[iradius - 1]);
    for (k=0; k<ndim; k++) r[ndim*i + k] = rsonic*r[ndim*i + k]*radp/drmag;
    for (k=0; k<ndim; k++) part.r[k] = r[ndim*i + k];
    for (k=0; k<ndim; k++) part.v[k] = -asound*vradp*dr[k]/drmag;
    part.m = mp;
  }


  // Now add star/sink particle
  rsink *= rsonic;
  nbody->Nstar = 1;
  sinks->Nsink = 1;
  for (k=0; k<ndim; k++) nbody->stardata[0].r[k] = 0.0;
  for (k=0; k<ndim; k++) nbody->stardata[0].v[k] = 0.0;
  nbody->stardata[0].m      = msink;
  nbody->stardata[0].radius = rsink;
  nbody->stardata[0].h      = nbody->kernp->invkernrange*nbody->stardata[0].radius;
  nbody->stardata[0].invh   = 1.0/nbody->stardata[0].h;
  sinks->sink[0].star        = &(nbody->stardata[0]);
  sinks->sink[0].radius      = rsink;
  sinks->sink[0].racc        = rsink;
  sinks->sink[0].mmax        = 0.0;
  sinks->sink[0].menc        = 0.0;


  // Find total mass inside sink and set to mmax
  for (i=0; i<Npart; i++) {
    Particle<ndim>& part = hydro->GetParticlePointer(i);
    for (k=0; k<ndim; k++) dr[k] = part.r[k] - rcentre[k];
    drsqd = DotProduct(dr, dr, ndim);
    if (drsqd > rsink*rsink) continue;
    sinks->sink[0].mmax += part.m;
  }

  cout << "mmax : " << sinks->sink[0].mmax << "    " << sinks->sink[0].mmax/mp << endl;
  cout << "rsink : " << rsink << "     rsink/rsonic : " << rsink/rsonic << endl;
  //exit(0);

  sim->initial_h_provided = false;

  delete[] v;
  delete[] r;

  return;
}



//=================================================================================================
//  Ic::EwaldDensity
/// Set-up simple density distributions to test Ewald periodic gravity field.
//=================================================================================================
template <int ndim>
void Ic<ndim>::EwaldDensity(void)
{
  // Only compile for 3 dimensions
  //-----------------------------------------------------------------------------------------------
  if (ndim == 3) {

    int i,k;                             // Particle and dimension counters
    int Nlattice1[3];                    // Lattice size
    int Npart;                           // No. of particles in lattice
    int periodicity;                     // orientation of the axis of symmetry
    FLOAT csound;                        // (Isothermal) sound speed
    FLOAT h0;                            // Slab scale height
    FLOAT a2inv;                         // Squared inverse scale height for cylinder
    FLOAT lambda;                        // Wavelength of perturbation
    FLOAT kwave;                         // Wave number of perturbing sound wave
    //FLOAT omegawave;                     // Angular frequency of sound wave
    FLOAT ugas;                          // Internal energy of gas
    FLOAT volume;                        // Simulation box volume
    FLOAT *r;                            // Particle positions

    // Make local copies of parameters for setting up problem
    Nlattice1[0]    = simparams->intparams["Nlattice1[0]"];
    Nlattice1[1]    = simparams->intparams["Nlattice1[1]"];
    Nlattice1[2]    = simparams->intparams["Nlattice1[2]"];
    string x_boundary_lhs = simparams->stringparams["boundary_lhs[0]"];
    string x_boundary_rhs = simparams->stringparams["boundary_rhs[0]"];
    string y_boundary_lhs = simparams->stringparams["boundary_lhs[1]"];
    string y_boundary_rhs = simparams->stringparams["boundary_rhs[1]"];
    string z_boundary_lhs = simparams->stringparams["boundary_lhs[2]"];
    string z_boundary_rhs = simparams->stringparams["boundary_rhs[2]"];
    string ic       = simparams->stringparams["ic"];
    FLOAT rhofluid1 = simparams->floatparams["rhofluid1"];
    FLOAT press1    = simparams->floatparams["press1"];
    FLOAT gamma     = simparams->floatparams["gamma_eos"];
    FLOAT gammaone  = gamma - (FLOAT) 1.0;
    FLOAT amp       = simparams->floatparams["amp"];
    FLOAT temp0     = simparams->floatparams["temp0"];
    FLOAT mu_bar    = simparams->floatparams["mu_bar"];
    //FLOAT zmax      = simparams->floatparams["zmax"];

    debug2("[Ic::EwaldDensity]");

    if (hydro->gas_eos == "isothermal") {
      ugas   = temp0/gammaone/mu_bar;
      press1 = gammaone*rhofluid1*ugas;
      csound = sqrt(press1/rhofluid1);
    }
    else {
      ugas   = press1/rhofluid1/gammaone;
      csound = sqrt(gamma*press1/rhofluid1);
    }

    lambda = simbox.boxsize[0];
    volume = simbox.boxsize[0]*simbox.boxsize[1]*simbox.boxsize[2];
    kwave  = twopi/lambda;
    //omegawave = twopi*csound/lambda;

    // Allocate local and main particle memory
    if (ndim == 3) {
      Npart = Nlattice1[0]*Nlattice1[1]*Nlattice1[2];
    }
    hydro->Nhydro = Npart;
    //Nlattice1[0] = Npart;
    sim->AllocateParticleMemory();
    r = new FLOAT[ndim*hydro->Nhydro];

    // determine orientation of the problem (variable periodicity)
    // controls orientation of density structure in the case of
    // mixed boundary conditions (i.e. normal to the plane or axis
    // of the cylinder)
    periodicity = 0;
    if (x_boundary_lhs == "periodic" && x_boundary_rhs == "periodic") {
      periodicity += 1;
    }
    if (y_boundary_lhs == "periodic" && y_boundary_rhs == "periodic") {
      periodicity += 2;
    }
    if (z_boundary_lhs == "periodic" && z_boundary_rhs == "periodic") {
      periodicity += 4;
    }


    // 1D sinusoidal density perturbation
    //=============================================================================================
    if (ic == "ewaldsine") {

      lambda    = simbox.boxmax[0] - simbox.boxmin[0];
      kwave     = twopi/lambda;
      //omegawave = twopi*csound/lambda;

      // Add regular distribution of SPH particles
      AddCubicLattice(Npart, Nlattice1, simbox, false, r);
      //AddRandomBox(Npart,r,simbox);

      // Add sinusoidal density perturbation to particle distribution
      AddSinusoidalDensityPerturbation(Npart, amp, lambda, r);

      // Set all other particle quantities
      //-------------------------------------------------------------------------------------------
      for (i=0; i<Npart; i++) {
        Particle<ndim>& part = hydro->GetParticlePointer(i);

        // Set positions in main array with corresponind velocity perturbation
        for (k=0; k<ndim; k++) part.r[k] = r[ndim*i + k];
        for (k=0; k<ndim; k++) part.v[k] = (FLOAT) 0.0;
        part.m = rhofluid1*volume/(FLOAT) Npart;
        part.h = hydro->h_fac*pow(part.m/rhofluid1,invndim);

        if (hydro->gas_eos == "isothermal") part.u = temp0/gammaone/mu_bar;
        else part.u = press1/rhofluid1/gammaone;

      }
      //-------------------------------------------------------------------------------------------


    }
    //=============================================================================================
    else if (ic == "ewaldsine2") {

      lambda = simbox.boxmax[0] - simbox.boxmin[0];
      kwave = twopi/lambda;

      // Add regular distribution of SPH particles
      AddCubicLattice(Npart, Nlattice1, simbox, false, r);

      volume = simbox.boxsize[0]*simbox.boxsize[1]*simbox.boxsize[2];
      FLOAT volp = volume/(FLOAT) Npart;

      // Set all other particle quantities
      //-------------------------------------------------------------------------------------------
      for (i=0; i<Npart; i++) {
        Particle<ndim>& part = hydro->GetParticlePointer(i);

        // Set positions in main array with corresponind velocity perturbation
        for (k=0; k<ndim; k++) part.r[k] = r[ndim*i + k];
        for (k=0; k<ndim; k++) part.v[k] = 0.0;
        part.rho = rhofluid1*(1.0+amp*sin(kwave*part.r[0]));
        part.m   = part.rho*volp;
        part.h   = hydro->h_fac*pow(part.m/part.rho,invndim);

        if (hydro->gas_eos == "isothermal") part.u = temp0/gammaone/mu_bar;
        else part.u = press1/rhofluid1/gammaone;

      }
      //-------------------------------------------------------------------------------------------

    }
    //=============================================================================================
    else if (ic == "ewaldslab") {

      // it is assumed that the slab is selfgravitating, so boundary conditions for gravity
      // must be set properly
      if ((periodicity != 3) && (periodicity != 5) && (periodicity != 6)) {
        printf("For this test boundary conditions must be periodic in two directions");
        exit(0);
      }

      h0 = csound/sqrtf(twopi*rhofluid1);

      // Add regular distribution of SPH particles
      AddCubicLattice(Npart, Nlattice1, simbox, false, r);

      volume = simbox.boxsize[0]*simbox.boxsize[1]*simbox.boxsize[2];
      FLOAT volp = volume/(FLOAT) Npart;

      // Set all other particle quantities
      //-------------------------------------------------------------------------------------------
      for (i=0; i<Npart; i++) {
        Particle<ndim>& part = hydro->GetParticlePointer(i);

        // Set positions in main array with corresponind velocity perturbation
        for (k=0; k<ndim; k++) part.r[k] = r[ndim*i + k];
        for (k=0; k<ndim; k++) part.v[k] = 0.0;

        // normal to the slab depends on chosen BCs
        switch (periodicity) {
          case 3:
            part.rho = rhofluid1/powf(cosh(part.r[2]/h0),2);
            break;
          case 5:
            part.rho = rhofluid1/powf(cosh(part.r[1]/h0),2);
            break;
          case 6:
            part.rho = rhofluid1/powf(cosh(part.r[0]/h0),2);
            break;
        }
        part.m   = part.rho*volp;
        part.h   = hydro->h_fac*pow(part.m/part.rho,invndim);

        if (hydro->gas_eos == "isothermal") part.u = temp0/gammaone/mu_bar;
        else part.u = press1/rhofluid1/gammaone;

      }

      //-------------------------------------------------------------------------------------------


    }
    //=============================================================================================
    else if (ic == "ewaldcylinder") {

      // it is assumed that cylinder is selfgravitating, so boundary conditions for gravity
      // must be set properly
      if ((periodicity != 1) && (periodicity != 2) && (periodicity != 4)) {
        printf("For this test boundary conditions must be periodic in one direction");
        exit(0);
      }

      a2inv = pi*rhofluid1*0.5/pow(csound,2);

      // Add regular distribution of SPH particles
      AddCubicLattice(Npart, Nlattice1, simbox, false, r);

      volume = simbox.boxsize[0]*simbox.boxsize[1]*simbox.boxsize[2];
      FLOAT volp = volume/(FLOAT) Npart;

      // Set all other particle quantities
      //-------------------------------------------------------------------------------------------
      for (i=0; i<Npart; i++) {
        Particle<ndim>& part = hydro->GetParticlePointer(i);

        // Set positions in main array with corresponind velocity perturbation
        for (k=0; k<ndim; k++) part.r[k] = r[ndim*i + k];
        for (k=0; k<ndim; k++) part.v[k] = 0.0;
        switch (periodicity) {
          case 1:
            part.rho = rhofluid1/pow((1.0+a2inv*(pow(part.r[1],2) + pow(part.r[2],2))),2);
            break;
          case 2:
            part.rho = rhofluid1/pow((1.0+a2inv*(pow(part.r[0],2) + pow(part.r[2],2))),2);
            break;
          case 4:
            part.rho = rhofluid1/pow((1.0+a2inv*(pow(part.r[0],2) + pow(part.r[1],2))),2);
            break;
        }
        part.m   = part.rho*volp;
        part.h   = hydro->h_fac*pow(part.m/part.rho,invndim);

        if (hydro->gas_eos == "isothermal") part.u = temp0/gammaone/mu_bar;
        else part.u = press1/rhofluid1/gammaone;

      }

      //-------------------------------------------------------------------------------------------

    }
    //=============================================================================================

    sim->initial_h_provided = true;

  }
  //-----------------------------------------------------------------------------------------------
  else {
    string message = "Error : Ewald tests only run in 3D";
    ExceptionHandler::getIstance().raise(message);
  }
  //-----------------------------------------------------------------------------------------------

  return;
}



//=================================================================================================
//  Ic::PlummerSphere
/// Generate a Plummer sphere containing either stars, gas or a mixture of both.
/// Uses the algorithm described by Aarseth et al. (197?).  Only valid for 3 dimensions.
//=================================================================================================
template <int ndim>
void Ic<ndim>::PlummerSphere(void)
{

  Nbody<ndim>* nbody = sim->nbody;

  // Only compile for 3 dimensions
  //-----------------------------------------------------------------------------------------------
  if (ndim == 3) {

    bool flag;                           // Aux. flag
    int i,j,k;                           // Particle and dimension counter
    FLOAT raux;                          // Aux. float variable
    //FLOAT rcentre[ndim];                 // Position of centre of Plummer sphere
    FLOAT vplummer;                      // Velocity of Plummer components
    FLOAT x1,x2,x3,x4,x5,x6,x7;          // Aux. random number variables
    FLOAT rad,vm,ve,t1,t2,w,z;           // Other variables

    // Local copies of important parameters
    int Nhydro        = simparams->intparams["Nhydro"];
    int Nstar       = simparams->intparams["Nstar"];
    FLOAT gamma_eos = simparams->floatparams["gamma_eos"];
    FLOAT gasfrac   = simparams->floatparams["gasfrac"];
    FLOAT starfrac  = simparams->floatparams["starfrac"];
    FLOAT mplummer  = simparams->floatparams["mplummer"];
    FLOAT rplummer  = simparams->floatparams["rplummer"];
    FLOAT radius    = simparams->floatparams["radius"];
    FLOAT rstar     = simparams->floatparams["rstar"];

    debug1("[Ic::PlummerSphere]");

    hydro->Nhydro = Nhydro;
    hydro->Ntot = Nhydro;
    nbody->Nstar = Nstar;
    sim->AllocateParticleMemory();

    //for (k=0; k<ndim; k++) rcentre[k] = 0.0;
    raux = gasfrac + starfrac;
    gasfrac /= raux;
    starfrac /= raux;


    // Loop over all particles (gas and stars)
    //=============================================================================================
    for (j=0; j<Nhydro+Nstar; j++) {

      do {
        flag = false;
        x1 = sim->randnumb->floatrand();
        x2 = sim->randnumb->floatrand();
        x3 = sim->randnumb->floatrand();

        if (x1 == (FLOAT) 0.0 && x2 == (FLOAT) 0.0 && x3 == (FLOAT) 0.0) flag = true;
        rad = (FLOAT) 1.0 / sqrt(pow(x1,-(FLOAT) 2.0/ (FLOAT) 3.0) - (FLOAT) 1.0);
        if (rad > radius/rplummer) flag = true;
      } while (flag);

      z = ((FLOAT) 1.0 - (FLOAT) 2.0*x2)*rad;


      // Set position depending on particle type
      //-------------------------------------------------------------------------------------------
      if (j >= Nstar && j < Nstar + Nhydro) {
        i = j - Nstar;
        Particle<ndim>& part = hydro->GetParticlePointer(i);
        part.r[2] = z;
        part.r[0] = sqrt(rad*rad - z*z)*cos(twopi*x3);
        part.r[1] = sqrt(rad*rad - z*z)*sin(twopi*x3);
        part.m = gasfrac / (FLOAT) Nhydro;
      }
      else {
        i = j;
        nbody->stardata[i].r[0] = sqrt(rad*rad - z*z)*cos(twopi*x3);
        nbody->stardata[i].r[1] = sqrt(rad*rad - z*z)*sin(twopi*x3);
        nbody->stardata[i].r[2] = z;
        nbody->stardata[i].m    = starfrac / (FLOAT) Nstar;
      }

      // Maximum velocity for this distance
      ve = sqrt((FLOAT) 2.0 / sqrt((FLOAT) 1.0 + rad*rad));


      // Velocity of particle
      //-------------------------------------------------------------------------------------------
      do {
        x4 = sim->randnumb->floatrand();
        x5 = sim->randnumb->floatrand();
        t1 = (FLOAT) 0.1*x5;
        t2 = x4*x4*pow((FLOAT) 1.0 - x4*x4,(FLOAT) 3.5);
      } while (t1 > t2);

      vm = ve*x4;
      x6 = sim->randnumb->floatrand();
      x7 = sim->randnumb->floatrand();
      w = ((FLOAT) 1.0 - (FLOAT) 2.0*x6)*vm;


      // Set velocity depending on particle type
      //-------------------------------------------------------------------------------------------
      if (j >= Nstar && j < Nstar + Nhydro) {
        i = j - Nstar;
        Particle<ndim>& part = hydro->GetParticlePointer(i);
        for (k=0; k<ndim; k++) part.v[k] = (FLOAT) 0.0;
        FLOAT sound = sqrt((FLOAT) 0.16666666666666666 / sqrt((FLOAT) 1.0 + rad*rad));
        part.rho   = (FLOAT) 1.0;
        part.u     = sound*sound/(gamma_eos - (FLOAT) 1.0);
      }
      else {
        i = j;
        nbody->stardata[i].v[0] = sqrt(vm*vm - w*w)*cos(twopi*x7);
        nbody->stardata[i].v[1] = sqrt(vm*vm - w*w)*sin(twopi*x7);
        nbody->stardata[i].v[2] = w;
      }

    }
    //=============================================================================================


    // Instanly move to COM
    //ConvertToComFrame();
    vplummer = sqrt(mplummer/rplummer);

    // Now scale variables to required physical size
    for (i=0; i<Nhydro; i++) {
      Particle<ndim>& part = hydro->GetParticlePointer(i);
      for (k=0; k<ndim; k++) {
        part.r[k] = part.r[k]*rplummer;
        part.v[k] = part.v[k]*vplummer;
      }
      part.m = part.m*mplummer;
      if (i < Nhydro) part.u = part.u*(mplummer/rplummer);
    }

    for (i=0; i<Nstar; i++) {
      for (k=0; k<ndim; k++) {
        nbody->stardata[i].r[k] = nbody->stardata[i].r[k]*rplummer;
        nbody->stardata[i].v[k] = nbody->stardata[i].v[k]*vplummer;
      }
      nbody->stardata[i].m      = nbody->stardata[i].m*mplummer;
      nbody->stardata[i].radius = rstar;
      nbody->stardata[i].h      = nbody->kernp->invkernrange*rstar;
      nbody->stardata[i].invh   = 1.0 / (nbody->stardata[i].h + small_number_dp);
    }

  }
  //-----------------------------------------------------------------------------------------------
  else {
    ExceptionHandler::getIstance().raise("PlummerSphere can only be generated in 3d");
  }
  //-----------------------------------------------------------------------------------------------

  return;
}



//=================================================================================================
//  Ic::BlastWave
/// Set-up spherical blast wave test.
//=================================================================================================
template <int ndim>
void Ic<ndim>::BlastWave(void)
{
  int i;                               // Particle counter
  int k;                               // Dimension counter
  int Nbox;                            // No. of particles in box
  int Nlattice[3];                     // Lattice size
  FLOAT drsqd;                         // Distance squared
  FLOAT mbox;                          // Total mass inside simulation box
  FLOAT volume;                        // Volume of box
  FLOAT *r;                            // Positions of all particles

  // Create local copies of initial conditions parameters
  Nlattice[0]          = simparams->intparams["Nlattice1[0]"];
  Nlattice[1]          = simparams->intparams["Nlattice1[1]"];
  Nlattice[2]          = simparams->intparams["Nlattice1[2]"];
  FLOAT gammaone       = simparams->floatparams["gamma_eos"] - (FLOAT) 1.0;
  FLOAT rhofluid       = simparams->floatparams["rhofluid1"];
  FLOAT press1         = simparams->floatparams["press1"];
  FLOAT press2         = simparams->floatparams["press2"];
  FLOAT radius         = simparams->floatparams["radius"];
  string particle_dist = simparams->stringparams["particle_distribution"];

  debug2("[Ic::BlastWave]");


  // Compute size and range of fluid bounding boxes
  //-----------------------------------------------------------------------------------------------
  if (ndim == 1) {
    volume = simbox.boxmax[0] - simbox.boxmin[0];
    Nbox = Nlattice[0];
  }
  else if (ndim == 2) {
    volume = (simbox.boxmax[0] - simbox.boxmin[0])*(simbox.boxmax[1] - simbox.boxmin[1]);
    Nbox = Nlattice[0]*Nlattice[1];
  }
  else if (ndim == 3) {
    volume = (simbox.boxmax[0] - simbox.boxmin[0])*
      (simbox.boxmax[1] - simbox.boxmin[1])*(simbox.boxmax[2] - simbox.boxmin[2]);
    Nbox = Nlattice[0]*Nlattice[1]*Nlattice[2];
  }
  mbox  = volume*rhofluid;

  // Allocate local and main particle memory
  hydro->Nhydro = Nbox;
  sim->AllocateParticleMemory();
  r = new FLOAT[ndim*hydro->Nhydro];

  // Set pointer to SPH particle data


  // Add a cube of random particles defined by the simulation bounding box and
  // depending on the chosen particle distribution
  if (particle_dist == "random") {
    AddRandomBox(Nbox, simbox, r);
  }
  else if (particle_dist == "cubic_lattice") {
    AddCubicLattice(Nbox, Nlattice, simbox, true, r);
  }
  else if (particle_dist == "hexagonal_lattice") {
    AddHexagonalLattice(Nbox, Nlattice, simbox, true, r);
  }
  else {
    string message = "Invalid particle distribution option";
    ExceptionHandler::getIstance().raise(message);
  }

  // Record positions in main memory
  for (i=0; i<Nbox; i++) {
    Particle<ndim>& part = hydro->GetParticlePointer(i);
    for (k=0; k<ndim; k++) part.r[k] = r[ndim*i + k];
    drsqd = DotProduct(part.r, part.r, ndim);
    for (k=0; k<ndim; k++) part.v[k] = 0.0;
    part.m = mbox/(FLOAT) Nbox;
    part.h = hydro->h_fac*pow(part.m/rhofluid,invndim);
    if (drsqd < radius*radius) {
      part.u = press1/rhofluid/gammaone;
    }
    else {
      part.u = press2/rhofluid/gammaone;
    }
  }

  // Set initial smoothing lengths and create initial ghost particles
  //-----------------------------------------------------------------------------------------------
  hydro->Nghost = 0;
  hydro->Nghostmax = hydro->Nhydromax - hydro->Nhydro;
  hydro->Ntot = hydro->Nhydro;
  for (i=0; i<hydro->Nhydro; i++) hydro->GetParticlePointer(i).active = true;

  // Search ghost particles
  //sim->sphneib->SearchBoundaryGhostParticles(0.0,simbox,sph);

  sim->initial_h_provided = true;
  sim->rebuild_tree = true;
  //sphneib->BuildTree(rebuild_tree,n,ntreebuildstep,ntreestockstep,
                     //hydro->Ntot,hydro->Nhydromax,partdata,sph,timestep);

  //sphneib->UpdateAllSphProperties(hydro->Nhydro,hydro->Ntot,partdata,sph,nbody);

  // Update neighbour tre
  sim->rebuild_tree = true;
  //sphneib->BuildTree(rebuild_tree,n,ntreebuildstep,ntreestockstep,
                     //hydro->Ntot,hydro->Nhydromax,partdata,sph,timestep);

  // Calculate all SPH properties
  //sphneib->UpdateAllSphProperties(hydro->Nhydro,hydro->Ntot,partdata,sph,nbody);


  delete[] r;

  return;
}


//=================================================================================================
//  Ic::SedovBlastWave
/// Set-up Sedov blast wave test
//=================================================================================================
template <int ndim>
void Ic<ndim>::SedovBlastWave(void)
{
  int i;                               // Particle counter
  int k;                               // Dimension counter
  int Nbox;                            // No. of particles in box
  int Ncold;                           // No. of cold particles
  int Nhot;                            // No. of hot particles
  int Nlattice[3];                     // Lattice size
  int *hotlist;                        // List of 'hot' particles
  FLOAT drmag;                         // Distance
  FLOAT drsqd;                         // Distance squared
  FLOAT mbox;                          // Total mass inside simulation box
  FLOAT r_hot;                         // Size of 'hot' region
  FLOAT ufrac;                         // Internal energy fraction
  FLOAT umax;                          // Maximum u of all particles
  FLOAT utot;                          // Total internal energy
  FLOAT volume;                        // Volume of box
  FLOAT *r;                            // Positions of all particles

  // Create local copies of initial conditions parameters
  Nlattice[0]    = simparams->intparams["Nlattice1[0]"];
  Nlattice[1]    = simparams->intparams["Nlattice1[1]"];
  Nlattice[2]    = simparams->intparams["Nlattice1[2]"];
  int smooth_ic  = simparams->intparams["smooth_ic"];
  FLOAT rhofluid = simparams->floatparams["rhofluid1"];
  FLOAT kefrac   = simparams->floatparams["kefrac"];
  string particle_dist = simparams->stringparams["particle_distribution"];

  debug2("[Ic::SedovBlastWave]");


  // Compute size and range of fluid bounding boxes
  //-----------------------------------------------------------------------------------------------
  if (ndim == 1) {
    volume = simbox.boxmax[0] - simbox.boxmin[0];
    Nbox = Nlattice[0];
  }
  else if (ndim == 2) {
    volume = (simbox.boxmax[0] - simbox.boxmin[0])*(simbox.boxmax[1] - simbox.boxmin[1]);
    Nbox = Nlattice[0]*Nlattice[1];
  }
  else if (ndim == 3) {
    volume = (simbox.boxmax[0] - simbox.boxmin[0])*
      (simbox.boxmax[1] - simbox.boxmin[1])*(simbox.boxmax[2] - simbox.boxmin[2]);
    Nbox = Nlattice[0]*Nlattice[1]*Nlattice[2];
  }
  mbox  = volume*rhofluid;
  ufrac = max((FLOAT) 0.0,(FLOAT) 1.0 - kefrac);
  Ncold = 0;
  Nhot  = 0;
  r_hot = hydro->h_fac*hydro->kernrange*simbox.boxsize[0]/Nlattice[0];

  // Allocate local and main particle memory
  hydro->Nhydro = Nbox;
  sim->AllocateParticleMemory();
  r = new FLOAT[ndim*hydro->Nhydro];
  hotlist = new int[hydro->Nhydro];

  // Add a cube of random particles defined by the simulation bounding box and
  // depending on the chosen particle distribution
  if (particle_dist == "random") {
    AddRandomBox(Nbox, simbox, r);
  }
  else if (particle_dist == "cubic_lattice") {
    AddCubicLattice(Nbox, Nlattice, simbox, true, r);
  }
  else if (particle_dist == "hexagonal_lattice") {
    AddHexagonalLattice(Nbox, Nlattice, simbox, true, r);
  }
  else {
    string message = "Invalid particle distribution option";
    ExceptionHandler::getIstance().raise(message);
  }

  // Record positions in main memory
  for (i=0; i<Nbox; i++) {
    Particle<ndim>& part = hydro->GetParticlePointer(i);
    for (k=0; k<ndim; k++) part.r[k] = r[ndim*i + k];
    for (k=0; k<ndim; k++) part.v[k] = 0.0;
    part.m = mbox/(FLOAT) Nbox;
    part.h = hydro->h_fac*pow(part.m/rhofluid,invndim);
    part.u = small_number;
    //part.v[0] = 0.1*part.r[0] + 0.2*part.r[1] - 0.7*part.r[2];
  }

  // Set initial smoothing lengths and create initial ghost particles
  //-----------------------------------------------------------------------------------------------
  hydro->Nghost = 0;
  hydro->Nghostmax = hydro->Nhydromax - hydro->Nhydro;
  hydro->Ntot = hydro->Nhydro;
  for (i=0; i<hydro->Nhydro; i++) hydro->GetParticlePointer(i).active = true;

  // Search ghost particles
  //sim->sphneib->SearchBoundaryGhostParticles(0.0,simbox,sph);

  sim->initial_h_provided = true;
  sim->rebuild_tree = true;
  //sphneib->BuildTree(rebuild_tree,n,ntreebuildstep,ntreestockstep,
                     //hydro->Ntot,hydro->Nhydromax,partdata,sph,timestep);

  //sphneib->UpdateAllSphProperties(hydro->Nhydro,hydro->Ntot,partdata,sph,nbody);

  // Update neighbour tre
  sim->rebuild_tree = true;
  //sphneib->BuildTree(rebuild_tree,n,ntreebuildstep,ntreestockstep,
                     //hydro->Ntot,hydro->Nhydromax,partdata,sph,timestep);

  // Calculate all SPH properties
  //sphneib->UpdateAllSphProperties(hydro->Nhydro,hydro->Ntot,partdata,sph,nbody);

  // Now calculate which particles are hot
  //-----------------------------------------------------------------------------------------------
  umax = (FLOAT) 0.0;
  utot = (FLOAT) 0.0;
  for (i=0; i<hydro->Nhydro; i++) {
    Particle<ndim>& part = hydro->GetParticlePointer(i);
    drsqd = DotProduct(part.r,part.r,ndim);
    if (drsqd < r_hot*r_hot) {
      hotlist[i] = 1;
      if (smooth_ic == 1) part.u = part.m*hydro->kernp->w0(hydro->kernp->kernrange*sqrt(drsqd)/r_hot);
      else part.u = part.m;
      utot += part.u;
      umax = max(umax,part.u);
      Nhot++;
    }
    else {
      hotlist[i] = 0;
      Ncold++;
    }
  }

  // Normalise the energies
  //-----------------------------------------------------------------------------------------------
  for (i=0; i<hydro->Nhydro; i++) {
    Particle<ndim>& part = hydro->GetParticlePointer(i);
    if (hotlist[i] == 1) {
      drmag = sqrt(DotProduct(part.r,part.r,ndim));
      part.u = part.u/utot/part.m;
      //,1.0e-6*umax/part.m);
      for (k=0; k<ndim; k++) part.v[k] = sqrt(2.0*kefrac*part.u)*part.r[k]/(drmag + small_number);
      part.u = ufrac*part.u;
      //,1.0e-6*umax/part.m);
    }
    else {
      part.u = 1.0e-6/part.m;
      //for (k=0; k<ndim; k++) part.v[k] = (FLOAT) 0.0;
    }
  }

  delete[] hotlist;
  delete[] r;

  return;
}



//=================================================================================================
//  Ic::ShearFlow
/// Create shear-flow to test effective shear viscosity.
//=================================================================================================
template <int ndim>
void Ic<ndim>::ShearFlow(void)
{
  // Only compile for 2 dimensions
  //-----------------------------------------------------------------------------------------------
  if (ndim == 2) {

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

    debug2("[Ic::ShearFlow]");


    // Compute size and range of fluid bounding boxes
    volume = (simbox.boxmax[0] - simbox.boxmin[0])*(simbox.boxmax[1] - simbox.boxmin[1]);
    Nbox   = Nlattice1[0]*Nlattice1[1];
    lambda = simbox.boxmax[1] - simbox.boxmin[1];
    kwave  = twopi/lambda;

    // Allocate local and main particle memory
    hydro->Nhydro = Nbox;
    sim->AllocateParticleMemory();
    r = new FLOAT[ndim*hydro->Nhydro];


    // Add particles from cubic lattice
    if (Nbox > 0) {
      AddCubicLattice(Nbox, Nlattice1, simbox, false, r);
      //AddHexagonalLattice(Nbox,Nlattice1,r,simbox,false);

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
  else {
    ExceptionHandler::getIstance().raise("Shear-flow test only in 2D");
  }
  //-----------------------------------------------------------------------------------------------

  return;
}



//=================================================================================================
//  Ic::SoundWave
/// Set-up isothermal sound-wave test.
//=================================================================================================
template <int ndim>
void Ic<ndim>::SoundWave(void)
{
  int i,k;                          // Particle and dimension counters
  int Nlattice1[ndim];              // Lattice size
  FLOAT csound;                     // (Isothermal) sound speed
  FLOAT lambda;                     // Wavelength of perturbation
  FLOAT kwave;                      // Wave number of perturbing sound wave
  //FLOAT omegawave;                  // Angular frequency of sound wave
  FLOAT ugas;                       // Internal energy of gas
  FLOAT *r;                         // Particle positions

  // Make local copies of parameters for setting up problem
  int Npart       = simparams->intparams["Nhydro"];
  FLOAT rhofluid1 = simparams->floatparams["rhofluid1"];
  FLOAT press1    = simparams->floatparams["press1"];
  FLOAT gamma     = simparams->floatparams["gamma_eos"];
  FLOAT gammaone  = gamma - (FLOAT) 1.0;
  FLOAT amp       = simparams->floatparams["amp"];
  FLOAT temp0     = simparams->floatparams["temp0"];
  FLOAT mu_bar    = simparams->floatparams["mu_bar"];
  //Nlattice1[0]    = simparams->intparams["Nlattice1[0]"];

  debug2("[Ic::SoundWave]");

  if (ndim != 1) {
    string message = "Sound wave only available in 1D";
    ExceptionHandler::getIstance().raise(message);
  }

  if (hydro->gas_eos == "isothermal") {
    ugas   = temp0/gammaone/mu_bar;
    press1 = gammaone*rhofluid1*ugas;
    csound = sqrt(press1/rhofluid1);
  }
  else {
    ugas   = press1/rhofluid1/gammaone;
    csound = sqrt(gamma*press1/rhofluid1);
  }

  lambda = simbox.boxmax[0] - simbox.boxmin[0];
  kwave = twopi/lambda;
  //omegawave = twopi*csound/lambda;

  // Allocate local and main particle memory
  for (k=0; k<ndim; k++) Nlattice1[k] = 1;
  Nlattice1[0] = Npart;
  hydro->Nhydro = Npart;
  sim->AllocateParticleMemory();
  r = new FLOAT[ndim*hydro->Nhydro];

  // Add regular distribution of SPH particles
  AddCubicLattice(Npart, Nlattice1, simbox, false, r);

  // Add sinusoidal density perturbation to particle distribution
  AddSinusoidalDensityPerturbation(Npart, amp, lambda, r);

  // Set all other particle quantities
  //----------------------------------------------------------------------------------------------
  for (i=0; i<Npart; i++) {
    Particle<ndim>& part = hydro->GetParticlePointer(i);

    // Set positions in main array with corresponind velocity perturbation
    for (k=0; k<ndim; k++) part.r[k] = r[ndim*i];
    for (k=0; k<ndim; k++) part.v[k] = csound*amp*sin(kwave*r[ndim*i]);
    part.m = rhofluid1*lambda/(FLOAT) Npart;
    part.h = hydro->h_fac*pow(part.m/rhofluid1,invndim);

    if (hydro->gas_eos == "isothermal") {
      part.u = temp0/gammaone/mu_bar;
    }
    else {
      part.u = press1/rhofluid1/gammaone;
    }
  }
  //-----------------------------------------------------------------------------------------------

  sim->initial_h_provided = true;
  delete[] r;

  return;
}



//=================================================================================================
//  Ic::SpitzerExpansion
/// Set-up Spitzer expansion simulation for single ionising source
//=================================================================================================
template <int ndim>
void Ic<ndim>::SpitzerExpansion(void)
{
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

  debug2("[Ic::SpitzerExpansion]");

  mcloud /= simunits.m.outscale;
  radius /= simunits.r.outscale;


  r = new FLOAT[ndim*Npart];
  for (i=0; i<ndim*Npart; i++) r[i] = (FLOAT) 0.0;

  // Add a sphere of random particles with origin 'rcentre' and radius 'radius'
  for (k=0; k<ndim; k++) rcentre[k] = (FLOAT) 0.0;

  // Create the sphere depending on the choice of initial particle distribution
  if (particle_dist == "random") {
    AddRandomSphere(Npart, rcentre, radius, r);
  }
  else if (particle_dist == "cubic_lattice" || particle_dist == "hexagonal_lattice") {
    Nsphere = AddLatticeSphere(Npart, rcentre, radius, particle_dist, r);
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

  return;
}


//=================================================================================================
//  Ic::IsothermSphere
/// Create a r^-2 density sphere of particles of given origin and radius.
//=================================================================================================
template <int ndim>
void Ic<ndim>::IsothermSphere(void)
{
  int i,k;                          // Particle and dimension counters
  FLOAT rcentre[ndim];              // Position of sphere centre
  FLOAT volume;                     // Volume of sphere
  FLOAT *r;                         // Particle position vectors

  // Local copies of important parameters
  int Npart            = simparams->intparams["Nsph"];
  FLOAT mcloud         = simparams->floatparams["mcloud"];
  FLOAT radius         = simparams->floatparams["radius"];
  FLOAT rhofluid       = simparams->floatparams["rhofluid1"];
  FLOAT press          = simparams->floatparams["press1"];
  FLOAT temp0          = simparams->floatparams["temp0"];
  FLOAT mu_bar         = simparams->floatparams["mu_bar"];
  FLOAT gammaone       = simparams->floatparams["gamma_eos"] - 1.0;
  string particle_dist = simparams->stringparams["particle_distribution"];

  debug2("[Ic::IsothermSphere]");

  mcloud   /= simunits.m.outscale;
  radius   /= simunits.r.outscale;
  rhofluid /= simunits.rho.outscale;
  press    /= simunits.press.outscale;
  r = new FLOAT[ndim*Npart];

  // Add a sphere of random particles with origin 'rcentre' and radius 'radius'
  for (k=0; k<ndim; k++) rcentre[k] = (FLOAT) 0.0;
  Addr2Sphere(Npart, r, rcentre, radius);

  hydro->Nhydro = Npart;
  sim->AllocateParticleMemory();

  if (ndim == 1) volume = 2.0*radius;
  else if (ndim == 2) volume = pi*radius*radius;
  else if (ndim == 3) volume = 4.0*onethird*pi*pow(radius,3);
  //if (mcloud > small_number && radius > small_number)
  //  rhofluid = mcloud / volume;
  rhofluid = mcloud / volume;


  // Record particle positions and initialise all other variables
#pragma omp parallel for default(none)\
  shared(gammaone,mcloud,Npart,press,r,rhofluid,volume,temp0,mu_bar) private(i,k)
  for (i=0; i<hydro->Nhydro; i++) {
    Particle<ndim>& part = hydro->GetParticlePointer(i);
    for (k=0; k<ndim; k++) {
      part.r[k] = r[ndim*i + k];
      part.v[k] = (FLOAT) 0.0;
      part.a[k] = (FLOAT) 0.0;
    }
    if (hydro->gas_eos == "isothermal") {
      part.u = temp0/gammaone/mu_bar;
    }
    else {
      part.u = press/rhofluid/gammaone;
    }
    part.m = mcloud / (FLOAT) Npart;
    part.h = hydro->h_fac*pow(part.m/rhofluid,invndim);

    part.iorig = i;
  }

  sim->initial_h_provided = true;

  delete[] r;

  return;
}



//=================================================================================================
//  Ic::RotIsothermSphere
/// Create a r^-2 density sphere of particles of given origin and radius with solid body rotation.
//=================================================================================================
template <int ndim>
void Ic<ndim>::RotIsothermSphere(void)
{
  int i,k;                             // Particle and dimension counters
  FLOAT rcentre[ndim];                 // Position of sphere centre
  FLOAT volume;                        // Volume of sphere
  FLOAT *r;                            // Particle position vectors
  FLOAT r_perp[3];                     // Perpendicular vector
  FLOAT r_center;                      // Radius to centre
  FLOAT norm;                          // ..
  FLOAT pc = 3.08567758*pow(10,16);    // ..

  // Local copies of important parameters
  int Npart            = simparams->intparams["Nsph"];
  FLOAT mcloud         = simparams->floatparams["mcloud"];
  FLOAT radius         = simparams->floatparams["radius"];
  FLOAT rhofluid       = simparams->floatparams["rhofluid1"];
  FLOAT press          = simparams->floatparams["press1"];
  FLOAT temp0          = simparams->floatparams["temp0"];
  FLOAT mu_bar         = simparams->floatparams["mu_bar"];
  FLOAT gammaone       = simparams->floatparams["gamma_eos"] - 1.0;
  FLOAT omega          = simparams->floatparams["omega"];
  FLOAT omega_inv      = 1.0/omega;
  string particle_dist = simparams->stringparams["particle_distribution"];

  debug2("[Ic::RotIsothermSphere]");

  mcloud    /= simunits.m.outscale;
  radius    /= simunits.r.outscale;
  rhofluid  /= simunits.rho.outscale;
  press     /= simunits.press.outscale;
  omega_inv /= simunits.t.outscale;
  std::cout << omega << endl;
  omega = 1 / omega_inv;
  std::cout << omega << endl;
  r = new FLOAT[ndim*Npart];

  // Add a sphere of random particles with origin 'rcentre' and radius 'radius'
  for (k=0; k<ndim; k++) rcentre[k] = (FLOAT) 0.0;
  Addr2Sphere(Npart, r, rcentre, radius);

  hydro->Nhydro = Npart;
  sim->AllocateParticleMemory();

  if (ndim == 1) volume = 2.0*radius;
  else if (ndim == 2) volume = pi*radius*radius;
  else if (ndim == 3) volume = 4.0*onethird*pi*pow(radius,3);
  //if (mcloud > small_number && radius > small_number)
  //  rhofluid = mcloud / volume;
  rhofluid = mcloud / volume;


  // Record particle positions and initialise all other variables
#pragma omp parallel for default(none)\
  shared(gammaone,mcloud,Npart,press,r,rhofluid,volume,temp0,mu_bar) private(i,k)
  for (i=0; i<hydro->Nhydro; i++) {
    Particle<ndim>& part = hydro->GetParticlePointer(i);
    for (k=0; k<ndim; k++) {
      part.r[k] = r[ndim*i + k];
      part.a[k] = (FLOAT) 0.0;
    }

    if (hydro->gas_eos == "isothermal") {
      part.u = temp0/gammaone/mu_bar;
    }
    else {
      part.u = press/rhofluid/gammaone;
    }
    part.m = mcloud / (FLOAT) Npart;
    part.h = hydro->h_fac*pow(part.m/rhofluid, invndim);
  }

  if (ndim == 3) {
    for (i=0; i<hydro->Nhydro; i++) {
      Particle<ndim>& part = hydro->GetParticlePointer(i);
      part.iorig = i;

      r_perp[0] = -part.r[1];
      r_perp[1] = part.r[0];
      r_perp[2] = (FLOAT) 0.0;
      norm      = pow(pow(r_perp[0],2) + pow(r_perp[1],2),0.5);
      r_perp[0] = r_perp[0]/norm;
      r_perp[1] = r_perp[1]/norm;
      r_perp[2] = r_perp[2]/norm;
      r_center  = pow(pow(part.r[0],2) + pow(part.r[1],2) + pow(part.r[2],2),0.5);
      for (k=0; k<ndim; k++) part.v[k] = r_perp[k]*r_center*omega*pc;
    }
  }

  sim->initial_h_provided = true;

  delete[] r;

  return;
}



//=================================================================================================
//  Ic::TurbIsothermSphere
/// Create a r^-2 density sphere of particles of given origin and radius with turbulent velocity.
//=================================================================================================
template <int ndim>
void Ic<ndim>::TurbIsothermSphere(void)
{
  int i,k;                          // Particle and dimension counters
  FLOAT rcentre[ndim];              // Position of sphere centre
  FLOAT gpecloud;                   // ..
  FLOAT keturb;                     // ..
  FLOAT mp;                         // Mass of one particle
  FLOAT rho;                        // Fluid density
  FLOAT xmin;                       // ..
  FLOAT vfactor;                    // ..
  FLOAT *r;                         // Positions of all particles
  FLOAT *v;                         // Velocities of all particles
  FLOAT dxgrid;                     // ..
  FLOAT rmax[ndim];                 // ..
  FLOAT rmin[ndim];                 // ..
  DOUBLE *vfield;                   // ..

  // Local copies of important parameters
  int Npart            = simparams->intparams["Nsph"];
  int field_type       = simparams->intparams["field_type"];
  int gridsize         = simparams->intparams["gridsize"];
  FLOAT mcloud         = simparams->floatparams["mcloud"];
  FLOAT radius         = simparams->floatparams["radius"];
  FLOAT temp0          = simparams->floatparams["temp0"];
  FLOAT mu_bar         = simparams->floatparams["mu_bar"];
  FLOAT gammaone       = simparams->floatparams["gamma_eos"] - 1.0;
  FLOAT alpha_turb     = simparams->floatparams["alpha_turb"];
  FLOAT power_turb     = simparams->floatparams["power_turb"];
  string particle_dist = simparams->stringparams["particle_distribution"];

  debug2("[Ic::TurbIsothermSphere]");

#if !defined(FFTW_TURBULENCE)
  string message = "FFTW turbulence flag not set";
  ExceptionHandler::getIstance().raise(message);
#endif

  // Convert any parameters to code units
  mcloud /= simunits.m.outscale;
  radius /= simunits.r.outscale;
  temp0  /= simunits.temp.outscale;

  // Calculate gravitational potential energy of uniform cloud
  gpecloud = 0.6*mcloud*mcloud/radius;

  r = new FLOAT[ndim*Npart];
  v = new FLOAT[ndim*Npart];

  // Add a sphere of random particles with origin 'rcentre' and radius 'radius'
  for (k=0; k<ndim; k++) rcentre[k] = (FLOAT) 0.0;
  Addr2Sphere(Npart, r, rcentre, radius);

  // Allocate local and main particle memory
  hydro->Nhydro = Npart;
  sim->AllocateParticleMemory();
  mp = mcloud / (FLOAT) Npart;
  rho = 3.0*mcloud / (4.0*pi*pow(radius,3));

  // Record particle properties in main memory
  for (i=0; i<Npart; i++) {
    Particle<ndim>& part = hydro->GetParticlePointer(i);
    for (k=0; k<ndim; k++) part.r[k] = r[ndim*i + k];
    for (k=0; k<ndim; k++) part.v[k] = v[ndim*i + k];
    part.m = mp;
    part.h = hydro->h_fac*pow(part.m/rho,invndim);
    part.u = temp0/gammaone/mu_bar;
  }

  sim->initial_h_provided = true;


  // Generate turbulent velocity field for given power spectrum slope
  vfield = new DOUBLE[ndim*gridsize*gridsize*gridsize];
  hydro->ComputeBoundingBox(rmax, rmin, hydro->Nhydro);
  xmin = 9.9e20;
  dxgrid = 0.0;
  for (k=0; k<ndim; k++) {
    dxgrid = max(dxgrid,(rmax[k] - rmin[k])/(FLOAT) (gridsize - 1));
    xmin = min(xmin, rmin[k]);
  }

  // Generate gridded velocity field
  GenerateTurbulentVelocityField(field_type, gridsize, power_turb, vfield);

  // Now interpolate generated field onto particle positions
  InterpolateVelocityField(hydro->Nhydro, gridsize, xmin, dxgrid, r, vfield, v);

  // Finally, copy velocities to main SPH particle array
  for (i=0; i<hydro->Nhydro; i++) {
    Particle<ndim>& part = hydro->GetParticlePointer(i);
    for (k=0; k<ndim; k++) part.v[k] = v[ndim*i + k];
  }

  // Change to COM frame of reference
  sim->SetComFrame();

  // Calculate total kinetic energy of turbulent velocity field
  keturb = 0.0;
  for (i=0; i<hydro->Nhydro; i++) {
    Particle<ndim>& part = hydro->GetParticlePointer(i);
    keturb += part.m*DotProduct(part.v, part.v, ndim);
  }
  keturb *= 0.5;

  vfactor = sqrt(alpha_turb*gpecloud/keturb);
  cout << "Scaling factor : " << vfactor << endl;

  // Now rescale velocities to give required turbulent energy in cloud
  for (i=0; i<hydro->Nhydro; i++) {
    Particle<ndim>& part = hydro->GetParticlePointer(i);
    for (k=0; k<ndim; k++) part.v[k] *= vfactor;
  }

  delete[] v;
  delete[] r;

  return;
}



//=================================================================================================
//  Ic::BinaryStar
/// Create a simple binary star problem
//=================================================================================================
template <int ndim>
void Ic<ndim>::BinaryStar(void)
{
  int k;                               // Dimension counter
  FLOAT rbinary[ndim];                // Position of binary COM
  FLOAT vbinary[ndim];                // Velocity of binary COM

  // Binary star parameters
  FLOAT abin     = simparams->floatparams["abin"];
  FLOAT ebin     = simparams->floatparams["ebin"];
  FLOAT m1       = simparams->floatparams["m1"];
  FLOAT m2       = simparams->floatparams["m2"];
  FLOAT phirot   = simparams->floatparams["phirot"];
  FLOAT thetarot = simparams->floatparams["thetarot"];
  FLOAT psirot   = simparams->floatparams["psirot"];

  Nbody<ndim>* nbody = sim->nbody;

  debug2("[Ic::BinaryStar]");

  if (ndim == 1) {
    string message = "Binary test not available in 1D";
    ExceptionHandler::getIstance().raise(message);
  }

  // Allocate local and main particle memory
  //hydro->Nhydro = 0;
  //hydro->Ntot = 0;
  nbody->Nstar = 2;
  sim->AllocateParticleMemory();

  // Add binary star
  for (k=0; k<ndim; k++) rbinary[k] = 0.0;
  for (k=0; k<ndim; k++) vbinary[k] = 0.0;
  AddBinaryStar(abin,ebin,m1,m2,0.01,0.01,phirot,thetarot,psirot,0.0,
                rbinary,vbinary,nbody->stardata[0],nbody->stardata[1]);

  return;
}



//=================================================================================================
//  Ic::TripleStar
/// Create a simple triple star problem
//=================================================================================================
template <int ndim>
void Ic<ndim>::TripleStar(void)
{
  int k;                               // Dimension counter
  FLOAT rbinary[ndim];                // Position of binary COM
  FLOAT vbinary[ndim];                // Velocity of binary COM
  NbodyParticle<ndim> b1;              // Inner binary COM particle

  // Triple star parameters
  FLOAT abin1    = simparams->floatparams["abin"];
  FLOAT abin2    = simparams->floatparams["abin2"];
  FLOAT ebin1    = simparams->floatparams["ebin"];
  FLOAT ebin2    = simparams->floatparams["ebin2"];
  FLOAT m1       = simparams->floatparams["m1"];
  FLOAT m2       = simparams->floatparams["m2"];
  FLOAT m3       = simparams->floatparams["m3"];
  FLOAT phirot   = simparams->floatparams["phirot"];
  FLOAT psirot   = simparams->floatparams["psirot"];
  FLOAT thetarot = simparams->floatparams["thetarot"];

  Nbody<ndim>* nbody = sim->nbody;

  debug2("[SphIc::TripleStar]");

  if (ndim == 1) {
    string message = "Quadruple test not available in 1D";
    ExceptionHandler::getIstance().raise(message);
  }

  // Allocate local and main particle memory
  //hydro->Nhydro = 0;
  //hydro->Ntot = 0;
  nbody->Nstar = 3;
  sim->AllocateParticleMemory();

  // Compute main binary orbit
  for (k=0; k<ndim; k++) rbinary[k] = 0.0;
  for (k=0; k<ndim; k++) vbinary[k] = 0.0;
  AddBinaryStar(abin1,ebin1,m1+m2,m3,0.0001,0.0001,phirot,thetarot,psirot,0.0,
                rbinary,vbinary,b1,nbody->stardata[2]);

  // Now compute both components
  AddBinaryStar(abin2,ebin2,m1,m2,0.0001,0.0001,phirot,thetarot,psirot,0.0,
                b1.r,b1.v,nbody->stardata[0],nbody->stardata[1]);

  return;
}



//=================================================================================================
//  Ic::QuadrupleStar
/// Create a simple quadruple star problem
//=================================================================================================
template <int ndim>
void Ic<ndim>::QuadrupleStar(void)
{
  int k;                               // Dimension counter
  FLOAT rbinary[ndim];                // Position of binary COM
  FLOAT vbinary[ndim];                // Velocity of binary COM
  NbodyParticle<ndim> b1;              // Star/binary 1
  NbodyParticle<ndim> b2;              // Star/binary 2

  // Quadruple star parameters
  FLOAT abin1    = simparams->floatparams["abin"];
  FLOAT abin2    = simparams->floatparams["abin2"];
  FLOAT ebin1    = simparams->floatparams["ebin"];
  FLOAT ebin2    = simparams->floatparams["ebin2"];
  FLOAT m1       = simparams->floatparams["m1"];
  FLOAT m2       = simparams->floatparams["m2"];
  FLOAT m3       = simparams->floatparams["m3"];
  FLOAT m4       = simparams->floatparams["m4"];
  FLOAT phirot   = simparams->floatparams["phirot"];
  FLOAT psirot   = simparams->floatparams["psirot"];
  FLOAT thetarot = simparams->floatparams["thetarot"];

  Nbody<ndim>* nbody = sim->nbody;

  debug2("[SphIc::QuadrupleStar]");

  if (ndim == 1) {
    string message = "Quadruple test not available in 1D";
    ExceptionHandler::getIstance().raise(message);
  }

  // Allocate local and main particle memory
  //hydro->Nhydro = 0;
  //hydro->Ntot = 0;
  nbody->Nstar = 4;
  sim->AllocateParticleMemory();

  // Compute main binary orbit
  for (k=0; k<ndim; k++) rbinary[k] = 0.0;
  for (k=0; k<ndim; k++) vbinary[k] = 0.0;
  AddBinaryStar(abin1,ebin1,m1+m2,m3+m4,0.01,0.01,phirot,thetarot,psirot,
                0.0,rbinary,vbinary,b1,b2);

  // Now compute components of both inner binaries
  AddBinaryStar(abin2,ebin2,m1,m2,0.0001,0.0001,phirot,thetarot,psirot,
                0.0,b1.r,b1.v,nbody->stardata[0],nbody->stardata[1]);
  AddBinaryStar(abin2,ebin2,m3,m4,0.0001,0.0001,phirot,thetarot,psirot,
                0.0,b2.r,b2.v,nbody->stardata[2],nbody->stardata[3]);

  return;
}



//=================================================================================================
//  Ic::AddBinaryStar
/// Add a binary star of given mass, eccentricity and separation.
/// (Code provided courtesy of S. P. Goodwin; 29/09/2013)
//=================================================================================================
template <int ndim>
void Ic<ndim>::AddBinaryStar
 (FLOAT sma,                          ///< Semi-major axis
  FLOAT eccent,                       ///< Orbital eccentricity
  FLOAT m1,                           ///< Mass of star 1
  FLOAT m2,                           ///< Mass of star 2
  FLOAT h1,                           ///< Smoothing length of star 1
  FLOAT h2,                           ///< Smoothing length of star 2
  FLOAT phirot,                       ///< 'phi' Euler rotation angle
  FLOAT thetarot,                     ///< 'theta' Euler rotation angle
  FLOAT phase,                        ///< Phase angle
  FLOAT psirot,                       ///< 'tpsi' rotation angle
  FLOAT *rbinary,                     ///< Position of COM of binary
  FLOAT *vbinary,                     ///< Velocity of COM of binary
  NbodyParticle<ndim> &s1,             ///< Star 1
  NbodyParticle<ndim> &s2)             ///< Star 2
{

  //-----------------------------------------------------------------------------------------------
  if (ndim > 1) {
    debug2("[Ic::AddBinaryStar]");

    int k;                               // Dimension counter
    FLOAT mbin = m1 + m2;                // Total binary mass

    // randomly sample M to get theta
    FLOAT M = 2.0*pi*sim->randnumb->floatrand();

    // from this solve to get eccentric anomoly E - e sin(theta) = M
    // N-R method x_1 = x_0 - f(x_0)/f'(x_0)
    FLOAT Ee = M;
    FLOAT old = M;
    do {
      old = Ee;
      Ee = Ee - (Ee - eccent*sin(Ee) - M)/(1.0 - eccent*cos(Ee));
    } while (fabs(old - Ee) > (FLOAT) 1.0e-6);

    // Next get theta tan(t/2) = sqrt((1+e)/(1-e))*tan(E/2)
    FLOAT theta = sqrt(((FLOAT) 1.0 + eccent)/((FLOAT) 1.0 - eccent))*tan((FLOAT) 0.5*Ee);
    theta = (FLOAT) 2.0*atan(theta);

    // Total separation
    FLOAT sep = sma*((FLOAT) 1.0 - eccent*eccent)/((FLOAT) 1.0 + eccent*cos(theta));

    // Get velocity
    FLOAT vel = (m1 + m2)*((FLOAT) 2.0/sep - (FLOAT) 1.0/sma);
    vel = sqrt(vel);

    // ..
    FLOAT hc = sqrt(((FLOAT) 1.0 + eccent*cos(theta))/((FLOAT) 2.0 - sep/sma));
    FLOAT phi = acos(hc);

    // Set properties of star 1
    // put on x-y plane to start
    for (k=0; k<ndim; k++) s1.r[k] = (FLOAT) 0.0;
    for (k=0; k<ndim; k++) s1.v[k] = (FLOAT) 0.0;
    s1.m = m1;
    s1.h = h1;
    s1.invh = (FLOAT) 1.0 / s1.h;
    s1.r[0] += sep*cos(theta)*m2/mbin;
    s1.r[1] += sep*sin(theta)*m2/mbin;
    s1.v[0] += -vel*cos((FLOAT) 0.5*pi - theta + phi)*m2/mbin;
    s1.v[1] += vel*sin((FLOAT) 0.5*pi - theta + phi)*m2/mbin;

    // Set properties of star 2
    for (k=0; k<ndim; k++) s2.r[k] = (FLOAT) 0.0;
    for (k=0; k<ndim; k++) s2.v[k] = (FLOAT) 0.0;
    s2.m = m2;
    s2.h = h2;
    s2.invh = (FLOAT) 1.0 / s2.h;
    s2.r[0] -= sep*cos(theta)*m1/mbin;
    s2.r[1] -= sep*sin(theta)*m1/mbin;
    s2.v[0] -= -vel*cos((FLOAT) 0.5*pi - theta + phi)*m1/mbin;
    s2.v[1] -= vel*sin((FLOAT) 0.5*pi - theta + phi)*m1/mbin;

    // Rotate binary to given orientation using Euler angles
    EulerAngleRotation(phirot,thetarot,psirot,s1.r);
    EulerAngleRotation(phirot,thetarot,psirot,s1.v);
    EulerAngleRotation(phirot,thetarot,psirot,s2.r);
    EulerAngleRotation(phirot,thetarot,psirot,s2.v);

    // Now move binary to given centre of position and velocity
    for (k=0; k<ndim; k++) s1.r[k] += rbinary[k];
    for (k=0; k<ndim; k++) s1.v[k] += vbinary[k];
    for (k=0; k<ndim; k++) s2.r[k] += rbinary[k];
    for (k=0; k<ndim; k++) s2.v[k] += vbinary[k];

  }
  //-----------------------------------------------------------------------------------------------
  else {
    string message = "Binary test not available in 1D";
    ExceptionHandler::getIstance().raise(message);
  }
  //-----------------------------------------------------------------------------------------------

  return;
}



//=================================================================================================
//  Ic::AddRandomBox
/// Populate given bounding box with random particles.
//=================================================================================================
template <int ndim>
void Ic<ndim>::AddRandomBox
 (const int Npart,                     ///< [in] No. of particles
  const DomainBox<ndim> box,           ///< [in] Bounding box containing particles
  FLOAT *r)                            ///< [out] Positions of particles
{
  debug2("[Ic::AddRandomBox]");
  assert(r);

  for (int i=0; i<Npart; i++) {
    for (int k=0; k<ndim; k++) {
      r[ndim*i + k] = box.boxmin[k] + (box.boxmax[k] - box.boxmin[k])*sim->randnumb->floatrand();
    }
  }

  return;
}



//=================================================================================================
//  Ic::AddRandomSphere
/// Add random sphere of particles
//=================================================================================================
template <int ndim>
void Ic<ndim>::AddRandomSphere
 (const int Npart,                     ///< [in] No. of particles in sphere
  const FLOAT rcentre[ndim],           ///< [in] Position of sphere centre
  const FLOAT radius,                  ///< [in] Radius of sphere
  FLOAT *r)                            ///< [out] Positions of particles in sphere
{
  int i,k;                             // Particle and dimension counters
  FLOAT rad;                           // Radius of particle
  FLOAT rpos[ndim];                    // Random position of new particle

  debug2("[Ic::AddRandomSphere]");
  assert(r);

  // Loop over all required particles
  //-----------------------------------------------------------------------------------------------
  for (i=0; i<Npart; i++) {

    // Continously loop until random particle lies inside sphere
    do {
      for (k=0; k<ndim; k++)
      rpos[k] = (FLOAT) 1.0 - (FLOAT) 2.0*sim->randnumb->floatrand();
      rad = DotProduct(rpos,rpos,ndim);
    } while (rad > 1.0);

    for (k=0; k<ndim; k++) r[ndim*i + k] = rcentre[k] + radius*rpos[k];
  }
  //-----------------------------------------------------------------------------------------------

  return;
}



//=================================================================================================
//  Ic::AddLatticeSphere
/// Add sphere of particles cut-out of regular lattice
//=================================================================================================
template <int ndim>
int Ic<ndim>::AddLatticeSphere
 (const int Npart,                     ///< [in] Desired/maximum no. of particles in sphere
  const FLOAT rcentre[ndim],           ///< [in] Position of sphere centre
  const FLOAT radius,                  ///< [in] Radius of sphere
  const string particle_dist,          ///< [in] String of lattice type
  FLOAT *r)                            ///< [out] Positions of particles in sphere
{
  int i,k;                             // Particle and dimension counters
  int Naux;                            // Aux. particle number
  int Nlattice[3];                     // Lattice size
  FLOAT theta;                         // Euler angle for random rotation of lattice
  FLOAT phi;                           //        "                      "
  FLOAT psi;                           //        "                      "
  FLOAT *raux;                         // Temp. array to hold particle positions
  DomainBox<ndim> box1;                // Bounding box

  debug2("[Ic::AddLatticeSphere]");
  assert(r);


  // Set parameters for box and lattice to ensure it contains enough particles
  for (k=0; k<3; k++) Nlattice[k] = 1;
  for (k=0; k<ndim; k++) Nlattice[k] = (int) (1.5*powf((FLOAT) Npart, invndim));
  for (k=0; k<ndim; k++) box1.boxmin[k] = -2.0;
  for (k=0; k<ndim; k++) box1.boxmax[k] = 2.0;
  Naux = Nlattice[0]*Nlattice[1]*Nlattice[2];
  raux = new FLOAT[ndim*Naux];

  // Create a bounding box to hold lattice sphere
  if (particle_dist == "cubic_lattice") {
    AddCubicLattice(Naux, Nlattice, box1, true, raux);
  }
  else if (particle_dist == "hexagonal_lattice") {
    AddHexagonalLattice(Naux, Nlattice, box1, true, raux);
  }
  else {
    string message = "Invalid particle distribution option";
    ExceptionHandler::getIstance().raise(message);
  }

  // Now cut-out sphere from lattice containing exact number of particles
  // (unless lattice structure prevents this).
  Naux = CutSphere(Npart, Naux, box1, false, raux);
  assert(Naux <= Npart);

  // Rotate sphere through random Euler angles (to prevent alignment problems
  // during tree construction)
  if (ndim == 2) {
    FLOAT rtemp[ndim];
    theta = twopi*sim->randnumb->floatrand();
    for (i=0; i<Naux; i++) {
      for (k=0; k<ndim; k++) rtemp[k] = raux[ndim*i + k];
      raux[ndim*i] = rtemp[0]*cos(theta) - rtemp[1]*sin(theta);
      raux[ndim*i + 1] = rtemp[0]*sin(theta) + rtemp[1]*cos(theta);
    }
  }
  else if (ndim == 3) {
    theta = acos(sqrtf(sim->randnumb->floatrand()));
    phi   = twopi*sim->randnumb->floatrand();
    psi   = twopi*sim->randnumb->floatrand();
    EulerAngleArrayRotation(Naux, phi, theta, psi, raux);
  }

  // Copy particle positions to main position array to be returned
  for (i=0; i<Naux; i++) {
    for (k=0; k<ndim; k++) r[ndim*i + k] = rcentre[k] + radius*raux[ndim*i + k];
  }

  // Free allocated memory
  delete[] raux;

  return Naux;
}



//=================================================================================================
//  Ic::Addr2Sphere
/// Add r^-2 sphere of particles
//=================================================================================================
template <int ndim>
void Ic<ndim>::Addr2Sphere
 (int Npart,                           ///< [in] No. of particles in sphere
  FLOAT *r,                            ///< [out] Positions of particles in sphere
  FLOAT rcentre[ndim],                 ///< [in] Position of sphere centre
  FLOAT radius)                        ///< [in] Radius of sphere
{
  int i;                               // Particle counter
  FLOAT phi;                           // ..
  FLOAT theta;                         // ..
  FLOAT sintheta;                      // ..
  FLOAT costheta;                      // ..
  FLOAT rpart;                         // ..

  debug2("[Ic::Addr2Sphere]");
  assert(r);

  // Loop over all required particles
  //-----------------------------------------------------------------------------------------------
  for (i=0; i<Npart; i++) {

    // Continously loop until random particle lies inside sphere
    phi      = (FLOAT) 2.0*pi*randnumb->floatrand();
    costheta = (FLOAT) 2.0*randnumb->floatrand() - (FLOAT) 1.0;
    theta    = acos(costheta);
    sintheta = sin(theta);
    rpart    = radius*randnumb->floatrand();
    r[ndim*i + 0] = rpart*sintheta*cos(phi);
    r[ndim*i + 1] = rpart*sintheta*sin(phi);
    r[ndim*i + 2] = rpart*costheta;

  }
  //-----------------------------------------------------------------------------------------------

  return;
}



//=================================================================================================
//  Ic::AddCubicLattice
/// Add regular (cubic) lattice of particles
//=================================================================================================
template <int ndim>
void Ic<ndim>::AddCubicLattice
 (const int Npart,                     ///< [in] No. of particles in lattice
  const int Nlattice[ndim],            ///< [in] Ptcls per dimension in lattice
  const DomainBox<ndim> box,           ///< [in] Bounding box of particles
  const bool normalise,                ///< [in] Normalise lattice shape and size
  FLOAT *r)                            ///< [out] Positions of particles
{
  int i,k;                             // Particle and dimension counters
  int ii,jj,kk;                        // Aux. lattice counters
  FLOAT spacing[ndim];                 // Lattice spacing in each direction

  debug2("[Ic::AddCubicLattice]");
  assert(r);

  // If normalised, ensure equal spacing between all lattice layers.
  // Otherwise set spacing to fit bounding box
  if (normalise) {
    for (k=0; k<ndim; k++) {
      spacing[k] = (box.boxmax[0] - box.boxmin[0])/(FLOAT) Nlattice[0];
    }
  }
  else {
    for (k=0; k<ndim; k++) {
      spacing[k] = (box.boxmax[k] - box.boxmin[k])/(FLOAT) Nlattice[k];
    }
  }


  // Create lattice depending on dimensionality
  //-----------------------------------------------------------------------------------------------
  if (ndim == 1) {
    for (ii=0; ii<Nlattice[0]; ii++) {
      i = ii;
      r[i] = box.boxmin[0] + ((FLOAT)ii + (FLOAT) 0.5)*spacing[0];
    }
  }
  //-----------------------------------------------------------------------------------------------
  else if (ndim == 2) {
    for (jj=0; jj<Nlattice[1]; jj++) {
      for (ii=0; ii<Nlattice[0]; ii++) {
        i = jj*Nlattice[0] + ii;
        r[ndim*i] = box.boxmin[0] + ((FLOAT)ii + (FLOAT) 0.5)*spacing[0];
        r[ndim*i + 1] = box.boxmin[1] + ((FLOAT)jj + (FLOAT) 0.5)*spacing[1];
      }
    }
  }
  //-----------------------------------------------------------------------------------------------
  else if (ndim == 3) {
#pragma omp parallel for default(none) shared(Nlattice,r,spacing) private(i,ii,jj,kk)
    for (kk=0; kk<Nlattice[2]; kk++) {
      for (jj=0; jj<Nlattice[1]; jj++) {
        for (ii=0; ii<Nlattice[0]; ii++) {
          i = kk*Nlattice[0]*Nlattice[1] + jj*Nlattice[0] + ii;
          r[ndim*i] = box.boxmin[0] + ((FLOAT)ii + (FLOAT) 0.5)*spacing[0];
          r[ndim*i + 1] = box.boxmin[1] + ((FLOAT)jj + (FLOAT) 0.5)*spacing[1];
          r[ndim*i + 2] = box.boxmin[2] + ((FLOAT)kk + (FLOAT) 0.5)*spacing[2];
        }
      }
    }
  }

  return;
}



//=================================================================================================
//  Ic::AddHexagonalLattice
/// Create simple hexagonal-packed lattice using A-B-A-B pattern in z-direction
/// N.B. the box is scaled to fit to the x-boxsize.
//=================================================================================================
template <int ndim>
void Ic<ndim>::AddHexagonalLattice
 (const int Npart,                     ///< [in] No. of particles in lattice
  const int Nlattice[ndim],            ///< [in] Ptcls per dimension in lattice
  const DomainBox<ndim> box,           ///< [in] Bounding box of particles
  const bool normalise,                ///< [in] Normalise lattice shape and size
  FLOAT *r)                            ///< [out] Positions of particles
{
  int i,k;                             // Particle and dimension counters
  int ii,jj,kk;                        // Aux. lattice counters
  FLOAT rad[ndim];                     // 'Radius' of particle in lattice

  debug2("[Ic::AddHexagonalLattice]");
  assert(r);

  // If normalised, ensure equal spacing between all particles.
  // Otherwise set spacing to fit bounding box.
  if (normalise) {
    for (k=0; k<ndim; k++) rad[k] = (FLOAT) 0.5*(box.boxmax[0] - box.boxmin[0])/(FLOAT) Nlattice[0];
  }
  else {
    for (k=0; k<ndim; k++) rad[k] = (FLOAT) 0.5*(box.boxmax[k] - box.boxmin[k])/(FLOAT) Nlattice[k];
  }


  // Create lattice depending on dimensionality
  //-----------------------------------------------------------------------------------------------
  if (ndim == 1) {
    for (ii=0; ii<Nlattice[0]; ii++) {
      i = ii;
      r[i] = box.boxmin[0] + (FLOAT) 0.5*rad[0] + (FLOAT) 2.0*(FLOAT)ii*rad[0];
    }
  }

  //-----------------------------------------------------------------------------------------------
  else if (ndim == 2) {
    for (jj=0; jj<Nlattice[1]; jj++) {
      for (ii=0; ii<Nlattice[0]; ii++) {
        i = jj*Nlattice[0] + ii;
        r[ndim*i] = box.boxmin[0] +
          (FLOAT) 0.5*rad[0] + ((FLOAT) 2.0*(FLOAT)ii + (FLOAT)(jj%2))*rad[0];
        r[ndim*i + 1] = box.boxmin[1] +
          (FLOAT) 0.5*sqrt((FLOAT) 3.0)*rad[1] + (FLOAT)jj*sqrt(3.0)*rad[1];
      }
    }
  }

  //-----------------------------------------------------------------------------------------------
  else if (ndim == 3) {
#pragma omp parallel for default(none) shared(Nlattice,r,rad) private(i,ii,jj,kk)
    for (kk=0; kk<Nlattice[2]; kk++) {
      for (jj=0; jj<Nlattice[1]; jj++) {
        for (ii=0; ii<Nlattice[0]; ii++) {
          i = kk*Nlattice[0]*Nlattice[1] + jj*Nlattice[0] + ii;
          r[ndim*i] = box.boxmin[0] + (FLOAT) 0.5*rad[0] +
            ((FLOAT) 2.0*(FLOAT) ii + (FLOAT) (jj%2) + (FLOAT) ((kk+1)%2))*rad[0];
          r[ndim*i + 1] = box.boxmin[1] + (FLOAT) 0.5*sqrt((FLOAT) 3.0)*rad[1] +
            (FLOAT) jj*sqrt((FLOAT) 3.0)*rad[1] + (FLOAT) (kk%2)*rad[1]/sqrt((FLOAT) 3.0);
          r[ndim*i + 2] = box.boxmin[2] + sqrt((FLOAT) 6.0)*rad[2]/(FLOAT) 3.0 +
            (FLOAT) kk*(FLOAT) 2.0*sqrt((FLOAT) 6.0)*rad[2]/(FLOAT) 3.0;
        }
      }
    }
  }

  return;
}



//=================================================================================================
//  Ic::CutSphere
/// Cut-out a sphere containing exactly 'Nsphere' particles from a uniform box of particles.
//=================================================================================================
template <int ndim>
int Ic<ndim>::CutSphere
 (const int Nsphere,                   ///< [in] Desired no. of particles in sphere
  const int Naux,                      ///< [in] No. of input particles
  const DomainBox<ndim> box,           ///< [in] Bounding box of particles
  const bool exact,                    ///< [in] Flag if we require exact numbers (not used)
  FLOAT *r)                            ///< [inout] Positions of particles
{
  int i,k;                             // Particle and dimension counters
  int Ninterior = 0;                   // No. of particle
  FLOAT dr[ndim];                      // Relative position vector
  FLOAT drsqd;                         // Distance squared
  FLOAT r_low = 0.0;                   // Lower-bound for bisection iteration
  FLOAT r_high;                        // Upper-bound for bisection iteration
  FLOAT radius;                        // Current radius containing Nsphere ptcls
  FLOAT rcentre[ndim];                 // Centre of sphere

  debug2("[Ic::CutSphere]");

  assert(r);

  // Find centre and shortest edge-length of bounding box
  r_high = (FLOAT) big_number;
  for (k=0; k<ndim; k++) {
    rcentre[k] = (FLOAT) 0.5*(box.boxmin[k] + box.boxmax[k]);
    r_high = min(r_high, (FLOAT) 0.5*(box.boxmax[k] - box.boxmin[k]));
  }

  // Bisection iteration to determine the radius containing the desired
  // number of particles
  //-----------------------------------------------------------------------------------------------
  do {
    radius = (FLOAT) 0.5*(r_low + r_high);
    Ninterior = 0;

    // Count how many particles lie inside current radius
    for (i=0; i<Naux; i++) {
      for (k=0; k<ndim; k++) dr[k] = r[ndim*i + k] - rcentre[k];
      drsqd = DotProduct(dr,dr,ndim);
      if (drsqd <= radius*radius) Ninterior++;
    }

    // If it's impossible to converge on the desired number of particles, due
    // to lattice configurations, then exit iteration with approximate number
    // of particles (must be less than Nsphere due to memory).
    if (Ninterior < Nsphere && fabs(r_high - r_low)/radius < 1.e-8) break;

    // Otherwise, continue bisection iteration to find radius
    if (Ninterior > Nsphere) r_high = radius;
    if (Ninterior < Nsphere) r_low = radius;

  } while (Ninterior != Nsphere);


  // Now that the radius containing require number has been identified,
  // record only the particles inside the sphere.
  //-----------------------------------------------------------------------------------------------
  Ninterior = 0;
  for (i=0; i<Naux; i++) {
    for (k=0; k<ndim; k++) dr[k] = r[ndim*i + k] - rcentre[k];
    drsqd = DotProduct(dr,dr,ndim);
    if (drsqd <= radius*radius) {
      for (k=0; k<ndim; k++) r[ndim*Ninterior + k] = r[ndim*i + k]/radius;
      Ninterior++;
    }
  }

  return Ninterior;
}



//=================================================================================================
//  Ic::AddAzimuthalDensityPerturbation
/// Add an azimuthal density perturbation for implementing Boss-Bodenheimer-type initial conditions
//=================================================================================================
template <int ndim>
void Ic<ndim>::AddAzimuthalDensityPerturbation
 (const int Npart,                       ///< [in] No. of particles in sphere
  const int mpert,                       ///< [in] Perturbation mode
  const FLOAT amp,                       ///< [in] Amplitude of perturbation
  const FLOAT rcentre[ndim],             ///< [in] Position of sphere centre
  FLOAT *r)                              ///< [inout] Positions of particles
{
  if (ndim > 1) {
    int i,k;                             // Particle and dimension counters
    int j;                               // Aux. counter
    FLOAT phi,phi1,phi2,phiprime;        // Aux. azimuthal angle variables
    FLOAT Rsqd;                          // Radial distance (from z-axis) squared
    FLOAT Rmag;                          // Radial distance (from z-axis)
    FLOAT rpos[2];                       // Random position of new particle
    const int tabtot = 2048;             // No of elements in tables
    const FLOAT invmpert = (FLOAT) 1.0/(FLOAT) mpert;  // 1 / mpert
    const FLOAT spacing = twopi/(FLOAT)(tabtot - 1);   // Table angular spacing

    debug2("[Ic::AddAzimuthalDensityPerturbation]");
    assert(r);


    // Loop over all required particles
    //---------------------------------------------------------------------------------------------
#pragma omp parallel for default(none) shared(r,rcentre)\
  private(i,j,k,phi,phiprime,phi1,phi2,rpos,Rmag,Rsqd)
    for (i=0; i<Npart; i++) {
      for (k=0; k<2; k++) rpos[k] = r[ndim*i + k] - rcentre[k];

      // Calculate distance from z-axis and
      Rsqd = rpos[0]*rpos[0] + rpos[1]*rpos[1];
      Rmag = sqrt(Rsqd);

      // Find azimuthal angle around z-axis correcting for which quadrant
      if (Rmag > small_number) phi = asin(fabs(rpos[1])/Rmag);
      else phi = (FLOAT) 0.0;
      phiprime = (FLOAT) 0.0;

      if (rpos[0] < (FLOAT) 0.0 && rpos[1] > (FLOAT) 0.0) phi = pi - phi;
      else if (rpos[0] < (FLOAT) 0.0 && rpos[1] < (FLOAT) 0.0) phi = pi + phi;
      else if (rpos[0] > (FLOAT) 0.0 && rpos[1] < (FLOAT) 0.0) phi = twopi - phi;

      // Wrap angle to fit between 0 and two*pi
      if (phi < amp*invmpert) phi = phi + twopi;

      // Numerically find new phi angle for perturbation.  Search through grid of values,
      // find upper and lower bounds, then use linear interpolation to find new value of phi.
      for (j=1; j<tabtot; j++) {
        phi1 = spacing*(FLOAT) (j - 1);
        phi2 = spacing*(FLOAT) j;
        phi1 = phi1 + amp*cos((FLOAT) mpert*phi1)*invmpert;
        phi2 = phi2 + amp*cos((FLOAT) mpert*phi2)*invmpert;

        if (phi2 >= phi && phi1 < phi) {
          phiprime = spacing*(FLOAT)(j - 1) + spacing*(phi - phi1) / (phi2 - phi1);
          r[ndim*i] = rcentre[0] + Rmag*cos(phiprime);
          r[ndim*i + 1] = rcentre[1] + Rmag*sin(phiprime);
          break;
        }
      }

      // Reposition particle with new angle
      //r[ndim*i] = rcentre[0] + Rmag*cos(phiprime);
      //r[ndim*i + 1] = rcentre[1] + Rmag*sin(phiprime);

    }
    //---------------------------------------------------------------------------------------------

  }

  return;
}



//=================================================================================================
//  Ic::AddSinusoidalDensityPerturbation
/// Add a 1D sinusoidal density perturbation (in x-direction) to given uniform density field.
//=================================================================================================
template <int ndim>
void Ic<ndim>::AddSinusoidalDensityPerturbation
 (const int Npart,                     ///< [in] No. of particles in sphere
  const FLOAT amp,                     ///< [in] Amplitude of perturbation
  const FLOAT lambda,                  ///< [in] Wave number of perturbation
  FLOAT *r)                            ///< [inout] Positions of particles
{
  int i;                               // Particle counter
  FLOAT diff;                          // Convergence error/difference
  FLOAT xold;                          // Old particle x-position in iteration
  FLOAT xnew;                          // New particle position
  FLOAT kwave = twopi/lambda;          // Sine wave-number

  debug2("[Ic::AddSinusoidalDensityPerturbation]");
  assert(r);


  // Loop over all required particles
  //-----------------------------------------------------------------------------------------------
#pragma omp parallel for default(none) shared(kwave,r) private(diff,i,xnew,xold)
  for (i=0; i<Npart; i++) {
    xnew = r[ndim*i];

    // Solve iterative procedure for particle positions in sinusoid
    do {
      xold = xnew;
      xnew = r[ndim*i] - amp*((FLOAT) 1.0 - cos(kwave*xnew))/kwave;
      diff = fabs((xnew - xold)/lambda);
    } while (diff > (FLOAT) 1.0e-6);

    if (xnew > simbox.boxmax[0]) xnew -= simbox.boxsize[0];
    if (xnew < simbox.boxmin[0]) xnew += simbox.boxsize[0];

    r[ndim*i] = xnew;
  }
  //-----------------------------------------------------------------------------------------------


  return;
}



//=================================================================================================
//  Ic::AddRotationalVelocityField
/// Add a solid-body rotational velocity field
//=================================================================================================
template <int ndim>
void Ic<ndim>::AddRotationalVelocityField
 (const int Npart,                     ///< [in] No. of particles in sphere
  const FLOAT angvelaux,               ///< [in] Angular velocity of cloud
  const FLOAT rcentre[ndim],           ///< [in] Position of sphere centre
  const FLOAT *r,                      ///< [in] Positions of particles
  FLOAT *v)                            ///< [out] Velocities of particles
{
  // Only compile for 2 or 3 dimensions
  //-----------------------------------------------------------------------------------------------
  if (ndim > 1) {

    int i,k;                           // Particle and dimension counters
    FLOAT Rmag;                        // Distance from z-axis
    FLOAT Rsqd;                        // Distance squared from z-axis
    FLOAT dr[3];                       // Relative position vector

    debug2("[Ic::AddRotationalVelocityField]");
    assert(r);
    assert(v);


    // Loop over all required particles
    //---------------------------------------------------------------------------------------------
#pragma omp parallel for default(none) shared(r,rcentre,v) private(dr,i,k,Rmag,Rsqd)
    for (i=0; i<Npart; i++) {
      for (k=0; k<ndim; k++) dr[k] = r[ndim*i + k] - rcentre[k];
      for (k=0; k<ndim; k++) v[ndim*i + k] = (FLOAT) 0.0;

      // Calculate distance from z-axis and
      Rsqd = dr[0]*dr[0] + dr[1]*dr[1] + small_number;
      Rmag = sqrt(Rsqd);

      // Find azimuthal angle around z-axis correcting for which quadrant
      if (Rmag > small_number) {
        dr[0] = dr[0]/Rmag;
        dr[1] = dr[1]/Rmag;

        v[ndim*i] = -angvelaux*Rmag*dr[1];
        v[ndim*i + 1] = angvelaux*Rmag*dr[0];
      }
    }
    //---------------------------------------------------------------------------------------------

  }
  //-----------------------------------------------------------------------------------------------

  return;
}



//=================================================================================================
//  Ic::ComputeBondiSolution
/// Compute the numerical solution to the Bondi accretion problem.
/// Translated from F90 subroutine written by A. P. Whitworth.
//=================================================================================================
template <int ndim>
void Ic<ndim>::ComputeBondiSolution
 (int Ntable,                          ///< ..
  FLOAT *w,                            ///< ..
  FLOAT *x,                            ///< ..
  FLOAT *y,                            ///< ..
  FLOAT *z)                            ///< ..
{
  const int isonic = 2*Ntable/3;       // Position in table for sonic point
  const FLOAT asympt = exp(1.5);       // ..
  FLOAT epsilon = 0.001;               // ..

  int i;                               // ..
  FLOAT disc,dx,hdx;                   // ..
  FLOAT f1,f2,f3,f4,f5;                // ..
  FLOAT g1,g2,g3,g4,g5;                // ..
  //FLOAT w0,x0,y0,z0;                   // ..
  FLOAT x1,x2,x4;                      // ..
  // FLOAT x3;
  FLOAT w1,w2,w3,w4;                   // ..
  FLOAT z1;                            // ..
  FLOAT *dlnx,*lw,*lx,*ly,*lz;         // ..

  debug2("[Ic::ComputeBondiSolution]");

  dlnx = new FLOAT[Ntable];
  lw   = new FLOAT[Ntable];
  lx   = new FLOAT[Ntable];
  ly   = new FLOAT[Ntable];
  lz   = new FLOAT[Ntable];

  // Set values for the sonic point directly
  x[isonic]  = (FLOAT) 1.0;
  lx[isonic] = (FLOAT) 0.0;
  w[isonic]  = (FLOAT) 1.0;
  lw[isonic] = (FLOAT) 0.0;
  y[isonic]  = asympt;
  ly[isonic] = log10(y[isonic]);
  z[isonic]  = (FLOAT) 2.4102434440;
  lz[isonic] = log10(z[isonic]);


  f5 = (FLOAT) 0.0;
  g5 = (FLOAT) 0.0;
  x1 = (FLOAT) 1.0 + epsilon;
  w1 = (FLOAT) 1.0 - epsilon;
  z1 = (FLOAT) 2.4102434440 + asympt*epsilon;
  i = isonic + 1;

  //-----------------------------------------------------------------------------------------------
  do {
    disc = (FLOAT) 100.0*log10(x1) + (FLOAT) isonic - (FLOAT) i;
    if (disc > (FLOAT) 0.0) {
      lx[i]   = (FLOAT) 0.01*(FLOAT)(i - isonic);
      x[i]    = powf((FLOAT) 10.0,lx[i]);
      dlnx[i] = (x[i] - x1)/x1;
      w[i]    = w1 + f5*(x[i] - x1);
      y[i]    = asympt/(x[i]*x[i]*w[i]);
      z[i]    = z1 + g5*(x[i] - x1);
      lw[i]   = log10(w[i]);
      ly[i]   = log10(y[i]);
      lz[i]   = log10(z[i]);
      i++;
    }
    dx = x1*epsilon;
    hdx = (FLOAT) 0.5*dx;
    f1 = (FLOAT) 2.0*(((FLOAT) 1.0/x1) - ((FLOAT) 1.0/(x1*x1)))/(w1 - ((FLOAT) 1.0/w1));
    g1 = asympt/w1;
    x2 = x1 + hdx;
    w2 = w1 + f1*hdx;
    f2 = (FLOAT) 2.0*(((FLOAT) 1.0/x2) - ((FLOAT) 1.0/(x2*x2)))/(w2 - ((FLOAT) 1.0/w2));
    g2 = asympt/w2;
    w3 = w1 + f2*hdx;
    f3 = (FLOAT) 2.0*(((FLOAT) 1.0/x2) - ((FLOAT) 1.0/(x2*x2)))/(w3 - ((FLOAT) 1.0/w3));
    g3 = asympt/w3;
    x4 = x1 + dx;
    w4 = w1 + f3*dx;
    f4 = (FLOAT) 2.0*(((FLOAT) 1.0/x4) - ((FLOAT) 1.0/(x4*x4)))/(w4 - ((FLOAT) 1.0/w4));
    g4 = asympt/w4;
    x1 = x4;
    f5 = (f1 + (FLOAT) 2.0*f2 + (FLOAT) 2.0*f3 + f4)/(FLOAT) 6.0;
    w1 = w1 + f5*dx;
    g5 = (g1 + (FLOAT) 2.0*g2 + (FLOAT) 2.0*g3 + g4)/(FLOAT) 6.0;
    z1 = z1 + g5*dx;
  } while (i < Ntable);


  epsilon = -epsilon;
  x1 = (FLOAT) 1.0 + epsilon;
  w1 = (FLOAT) 1.0 - epsilon;
  z1 = (FLOAT) 2.4102434440 + asympt*epsilon;
  i = isonic - 1;

  do {
    disc = (FLOAT) i - (FLOAT) 100.0*log10(x1) - (FLOAT) isonic;
    if (disc > (FLOAT) 0.0) {
      lx[i]   = (FLOAT) 0.01*(FLOAT) (i - isonic);
      x[i]    = powf((FLOAT) 10.0,lx[i]);
      dlnx[i] = (x[i] - x1)/x1;
      w[i]    = w1 + f5*(x[i] - x1);
      lw[i]   = log10(w[i]);
      y[i]    = asympt/(x[i]*x[i]*w[i]);
      ly[i]   = log10(y[i]);
      z[i]    = z1 + g5*(x[i] - x1);
      lz[i]   = -(FLOAT) 8.0;
      if (z[i] > (FLOAT) 0.1E-07) lz[i] = log10(z[i]);
      i--;
    }
    dx = x1*epsilon;
    hdx = (FLOAT) 0.5*dx;
    f1 = (FLOAT) 2.0*(((FLOAT) 1.0/x1) - ((FLOAT) 1.0/(x1*x1)))/(w1 - ((FLOAT) 1.0/w1));
    g1 = asympt/w1;
    x2 = x1 + hdx;
    w2 = w1 + f1*hdx;
    f2 = (FLOAT) 2.0*(((FLOAT) 1.0/x2) - ((FLOAT) 1.0/(x2*x2)))/(w2 - ((FLOAT) 1.0/w2));
    g2 = asympt/w2;
    w3 = w1 + f2*hdx;
    f3 = (FLOAT) 2.0*(((FLOAT) 1.0/x2) - ((FLOAT) 1.0/(x2*x2)))/(w3 - ((FLOAT) 1.0/w3));
    g3 = asympt/w3;
    x4 = x1 + dx;
    w4 = w1 + f3*dx;
    f4 = (FLOAT) 2.0*(((FLOAT) 1.0/x4) - ((FLOAT) 1.0/(x4*x4)))/(w4 - ((FLOAT) 1.0/w4));
    g4 = asympt/w4;
    x1 = x4;
    f5 = (f1 + (FLOAT) 2.0*f2 + (FLOAT) 2.0*f3 + f4)/(FLOAT) 6.0;
    w1 = w1 + f5*dx;
    g5 = (g1 + (FLOAT) 2.0*g2 + (FLOAT) 2.0*g3 + g4)/(FLOAT) 6.0;
    z1 = z1 + g5*dx;
  } while (i >= 0);


  // Write to file (for checking)
  ofstream outfile;
  string filename = "BONDI_SOLUTION.dat";
  outfile.open(filename.c_str());
  for (i=0; i<Ntable; i++) {
    outfile << x[i] << "    " << z[i] << "    " << y[i] << "    " << w[i] << endl;
  }
  outfile.close();


  // Free up all locally allocated arrays
  delete[] lz;
  delete[] ly;
  delete[] lx;
  delete[] lw;
  delete[] dlnx;

  return;
}



//=================================================================================================
//  Ic::GenerateTurbulentVelocityField
/// Generates turbulent velocity field using FFTW library.
/// Based on original code by A. McLeod.
//=================================================================================================
template <int ndim>
void Ic<ndim>::GenerateTurbulentVelocityField
 (const int field_type,                ///< [in] Type of turbulent velocity field
  const int gridsize,                  ///< [in] Size of velocity grid
  const DOUBLE power_turb,             ///< [in] Power spectrum index
  DOUBLE *vfield)                      ///< [out] Array containing velocity field grid
{
#if defined(FFTW_TURBULENCE)

  // Only valid for 3 dimensions
  //-----------------------------------------------------------------------------------------------
  if (ndim == 3) {

    bool divfree;                      // Select div-free turbulence
    bool curlfree;                     // Select curl-free turbulence
    int kmax;                          // Max. extent of k (in 3D)
    int kmin;                          // Min. extent of k (in 3D)
    int krange;                        // Range of k values (kmax - kmin + 1)
    int i,j,k;                         // Grid counters
    int ii,jj,kk;                      // Aux. grid counters
    int kmod;                          // ..
    int d;                             // Dimension counter
    DOUBLE F[ndim];                    // Fourier vector component
    DOUBLE unitk[3];                   // Unit k-vector
    DOUBLE **power,**phase;            // ..
    DOUBLE Rnd[3];                     // Random numbers, variable in Gaussian calculation
    //DOUBLE k_rot[ndim];                // bulk rotation modes
    //DOUBLE k_com[ndim];                // bulk compression modes
    fftw_plan plan;                    // ??
    fftw_complex *incomplexfield;      // ..
    fftw_complex *outcomplexfield;     // ..

    debug2("[Ic::GenerateTurbulentVelocityField]");


    divfree = false;
    curlfree = false;
    if (field_type == 1) curlfree = true;
    if (field_type == 2) divfree = true;

    kmin = -(gridsize/2 - 1);
    kmax = gridsize/2;
    krange = kmax - kmin + 1;
//    cout << "kmin : " << kmin << "    kmax : " << kmax << "    krange : "
//         << krange << "    gridsize : " << gridsize << endl;
    if (krange != gridsize) exit(0);

    power = new DOUBLE*[3];  phase = new DOUBLE*[3];
    for (d=0; d<3; d++) power[d] = new DOUBLE[krange*krange*krange];
    for (d=0; d<3; d++) phase[d] = new DOUBLE[krange*krange*krange];
    incomplexfield = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*gridsize*gridsize*gridsize);
    outcomplexfield = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*gridsize*gridsize*gridsize);

    for (d=0; d<3; d++) {
      for (i=0; i<krange*krange*krange; i++) power[d][i] = 0.0;
      for (i=0; i<krange*krange*krange; i++) phase[d][i] = 0.0;
    }
    for (i=0; i<3*gridsize*gridsize*gridsize; i++) vfield[i] = 0.0;


    // Define wave vectors in Fourier space.  Each wave vector has coordinates in Fourier space,
    // random phases in three dimensions, and power in three dimensions, giving a power and a
    // phase vector field in Fourier space.  With a 'natural' field type there is no coordination
    // between x,y,z components.  With a div-free or curl-free field, there is a relation between
    // coordinates in Fourier space and power vector. For a div-free field, power vectors are
    // perpendicular to the Fourier coordinate vector, and for a curl-free field power vectors
    // are (anti-)parallel to the Fourier coordinate vector
    // (i,j,k) is the Fourier coordinate vector
    // power(1:3,i,j,k) is the power vector at Fourier coordinates i,j,k
    // phase(1:3,i,j,k) is the phase vector at Fourier coordinates i,j,k
    //---------------------------------------------------------------------------------------------
    for (kmod=0; kmod<=kmax; kmod++) {

      for (i=kmin; i<=kmax; i++) {
        for (j=kmin; j<=kmax; j++) {
          for (k=kmin; k<=kmax; k++) {

            // Skip any k-vectors that have already been calculated
            if (abs(i) != kmod && abs(j) != kmod && abs(k) != kmod) continue;
            if (abs(i) > kmod || abs(j) > kmod || abs(k) > kmod) continue;

            //continue;
            ii = (i + krange)%krange;
            jj = (j + krange)%krange;
            kk = (k + krange)%krange;

            // cycle antiparallel k-vectors
            //if (k < 0) continue;
            //if (k == 0) {
            //  if (j < 0) continue;
            //  if (j == 0 && i < 0) continue;
            //}

            // Central power = 0
            if (i == 0 && j == 0 && k == 0) continue;
            if (i*i + j*j + k*k >= kmax*kmax) continue;

            // Power value, to be multipled by random power chosen from a Gaussian
            // This is what gives the slope of the eventual power spectrum
            for (d=0; d<3; d++) F[d] = sqrt(pow(sqrt((DOUBLE)(i*i + j*j + k*k)),power_turb));

            for (d=0; d<3; d++) {
              Rnd[0] = randnumb->floatrand();

              // Random phase between 0 and 2*pi (actually -pi and +pi).
              phase[d][ii + krange*jj + krange*krange*kk] = (2.0*Rnd[0] - 1.0)*pi;

              // Create Gaussian distributed random numbers
              Rnd[1] = randnumb->gaussrand(0.0,1.0);
              Rnd[2] = randnumb->gaussrand(0.0,1.0);
              F[d] = Rnd[1]*F[d];
            }

            // Helmholtz decomposition!
            unitk[0] = (DOUBLE) i;
            unitk[1] = (DOUBLE) j;
            unitk[2] = (DOUBLE) k;
            DOUBLE ksqd = DotProduct(unitk,unitk,3);
            for (d=0; d<3; d++) unitk[d] /= sqrt(ksqd);

            // For curl free turbulence, vector F should be
            // parallel/anti-parallel to vector k
            if (curlfree) {
              for (d=0; d<3; d++) {
                power[d][ii + krange*jj + krange*krange*kk] = unitk[d]*DotProduct(F,unitk,3);
              }
            }

            // For divergence free turbulence, vector F should be perpendicular to vector k
            else if (divfree) {
              for (d=0; d<3; d++) {
                power[d][ii + krange*jj + krange*krange*kk] = F[d] - unitk[d]*DotProduct(F,unitk,3);
              }
            }
            else {
              for (d=0; d<3; d++) power[d][ii + krange*jj + krange*krange*kk] = F[d];
            }

          }
        }
      }

    }
    //---------------------------------------------------------------------------------------------


    DOUBLE power_spectrum[kmax+1][3];
    for (i=0; i<kmax+1; i++) {
      for (d=0; d<3; d++) power_spectrum[i][d] = 0.0;
    }

    for (i=kmin; i<=kmax; i++) {
      for (j=kmin; j<=kmax; j++) {
        for (k=kmin; k<=kmax; k++) {
          ii = (i + krange)%krange;
          jj = (j + krange)%krange;
          kk = (k + krange)%krange;
          DOUBLE kmag = sqrt((DOUBLE)(i*i + j*j + k*k));
          if (kmag >= kmax) kmag = kmax;
          int ibin = (int) (1.0001*kmag);
          power_spectrum[ibin][0] += pow(power[0][ii + krange*jj + krange*krange*kk],2);
          power_spectrum[ibin][1] += pow(power[1][ii + krange*jj + krange*krange*kk],2);
          power_spectrum[ibin][2] += pow(power[2][ii + krange*jj + krange*krange*kk],2);
        }
      }
    }

    plan = fftw_plan_dft_3d(gridsize, gridsize, gridsize, incomplexfield,
                            outcomplexfield, FFTW_BACKWARD, FFTW_ESTIMATE);


    // reorder array: positive wavenumbers are placed in ascending order along
    // first half of dimension, i.e. 0 to k_max, negative wavenumbers are placed
    // along second half of dimension, i.e. -k_min to 1.
    for (d=0; d<3; d++) {

      for (i=kmin; i<=kmax; i++) {
        for (j=kmin; j<=kmax; j++) {
          for (k=kmin; k<=kmax; k++) {
            ii = (i + krange)%krange;
            jj = (j + krange)%krange;
            kk = (k + krange)%krange;
            //ii = i - kmin;
            //jj = j - kmin;
            //kk = k - kmin;
            incomplexfield[ii +krange*jj + krange*krange*kk][0] =
              power[d][ii + krange*jj + krange*krange*kk]*
              cos(phase[d][ii + krange*jj + krange*krange*kk]);
            incomplexfield[ii +krange*jj + krange*krange*kk][1] =
              power[d][ii + krange*jj + krange*krange*kk]*
              sin(phase[d][ii + krange*jj + krange*krange*kk]);

          }
        }
      }

      //fftw_execute_dft(plan, complexfield, complexfield);
      fftw_execute_dft(plan, incomplexfield, outcomplexfield);


      for (i=0; i<krange; i++) {
        for (j=0; j<krange; j++) {
          for (k=0; k<krange; k++) {
            vfield[d + 3*i + 3*krange*j + 3*krange*krange*k] =
              outcomplexfield[i + krange*j + krange*krange*k][0];

          }
        }
      }

    }


    fftw_destroy_plan(plan);

    for (d=0; d<3; d++) delete[] phase[d];
    for (d=0; d<3; d++) delete[] power[d];
    delete[] phase;
    delete[] power;
    fftw_free(outcomplexfield);
    fftw_free(incomplexfield);


  }
  //-----------------------------------------------------------------------------------------------

#endif

  return;
}



//=================================================================================================
//  Ic::InterpolateVelocityField
/// Calculate Interpolated velocity from uniform grid onto particle positions.
//=================================================================================================
template <int ndim>
void Ic<ndim>::InterpolateVelocityField
 (const int Npart,                       ///< [in] No of particles
  const int Ngrid,                       ///< [in] Size (per dim) of velocity grid
  const FLOAT xmin,                      ///< [in] Minimum position
  const FLOAT dxgrid,                    ///< [in] Grid size
  const FLOAT *r,                        ///< [in] Positions of particles
  const DOUBLE *vfield,                   ///< [in] Tabulated velocity field
  FLOAT *v)                              ///< [out] Interpolated particle velocity
{

  // Only compile routine for 3 dimensions
  //-----------------------------------------------------------------------------------------------
  if (ndim == 3) {

    int i,j,k;                           // Grid coordinates
    int kk;                              // Dimension counter
    int p;                               // Particle counter
    FLOAT dx[3];                         // Position relative to grid point
    FLOAT vint[8];                       // Interpolated velocity


    // Now interpolate velocity field onto particle positions
    //---------------------------------------------------------------------------------------------
#pragma omp parallel for default(none) private(dx,i,j,k,kk,p,vint) \
    shared(cout,r,v,vfield)
    for (p=0; p<Npart; p++) {
      for (kk=0; kk<ndim; kk++) dx[kk] = (r[ndim*p + kk] - xmin)/dxgrid;

      i = (int) dx[0];
      j = (int) dx[1];
      k = (int) dx[2];

      if (i > Ngrid || j > Ngrid || k > Ngrid || i < 0 || j < 0 || k < 0)  {
        cout << "Problem with velocity interpolation grid!! : "
             << i << "    " << j << "    " << k << "   " << Ngrid << endl;
        exit(0);
      }

      for (kk=0; kk<ndim; kk++) dx[kk] -= (int) dx[kk];

      // Interpolate to get more accurate velocities
      vint[0] = ((FLOAT) 1.0 - dx[0])*((FLOAT) 1.0 - dx[1])*((FLOAT) 1.0 - dx[2]);
      vint[1] = ((FLOAT) 1.0 - dx[0])*((FLOAT) 1.0 - dx[1])*dx[2];
      vint[2] = ((FLOAT) 1.0 - dx[0])*dx[1]*((FLOAT) 1.0 - dx[2]);
      vint[3] = ((FLOAT) 1.0 - dx[0])*dx[1]*dx[2];
      vint[4] = dx[0]*((FLOAT) 1.0 - dx[1])*((FLOAT) 1.0 - dx[2]);
      vint[5] = dx[0]*((FLOAT) 1.0 - dx[1])*dx[2];
      vint[6] = dx[0]*dx[1]*((FLOAT) 1.0 - dx[2]);
      vint[7] = dx[0]*dx[1]*dx[2];

      v[ndim*p]   = vint[0]*vfield[3*i + 3*Ngrid*j + 3*Ngrid*Ngrid*k] +
                    vint[1]*vfield[3*i + 3*Ngrid*j + 3*Ngrid*Ngrid*(k+1)] +
                    vint[2]*vfield[3*i + 3*Ngrid*(j+1) + 3*Ngrid*Ngrid*k] +
                    vint[3]*vfield[3*i + 3*Ngrid*(j+1) + 3*Ngrid*Ngrid*(k+1)] +
                    vint[4]*vfield[3*(i+1) + 3*Ngrid*j + 3*Ngrid*Ngrid*k] +
                    vint[5]*vfield[3*(i+1) + 3*Ngrid*j + 3*Ngrid*Ngrid*(k+1)] +
                    vint[6]*vfield[3*(i+1) + 3*Ngrid*(j+1) + 3*Ngrid*Ngrid*k] +
                    vint[7]*vfield[3*(i+1) + 3*Ngrid*(j+1) + 3*Ngrid*Ngrid*(k+1)];

      v[ndim*p+1] = vint[0]*vfield[1 + 3*i + 3*Ngrid*j + 3*Ngrid*Ngrid*k] +
                    vint[1]*vfield[1 + 3*i + 3*Ngrid*j + 3*Ngrid*Ngrid*(k+1)] +
                    vint[2]*vfield[1 + 3*i + 3*Ngrid*(j+1) + 3*Ngrid*Ngrid*k] +
                    vint[3]*vfield[1 + 3*i + 3*Ngrid*(j+1) + 3*Ngrid*Ngrid*(k+1)] +
                    vint[4]*vfield[1 + 3*(i+1) + 3*Ngrid*j + 3*Ngrid*Ngrid*k] +
                    vint[5]*vfield[1 + 3*(i+1) + 3*Ngrid*j + 3*Ngrid*Ngrid*(k+1)] +
                    vint[6]*vfield[1 + 3*(i+1) + 3*Ngrid*(j+1) + 3*Ngrid*Ngrid*k] +
                    vint[7]*vfield[1 + 3*(i+1) + 3*Ngrid*(j+1) + 3*Ngrid*Ngrid*(k+1)];

      v[ndim*p+2] = vint[0]*vfield[2 + 3*i + 3*Ngrid*j + 3*Ngrid*Ngrid*k] +
                    vint[1]*vfield[2 + 3*i + 3*Ngrid*j + 3*Ngrid*Ngrid*(k+1)] +
                    vint[2]*vfield[2 + 3*i + 3*Ngrid*(j+1) + 3*Ngrid*Ngrid*k] +
                    vint[3]*vfield[2 + 3*i + 3*Ngrid*(j+1) + 3*Ngrid*Ngrid*(k+1)] +
                    vint[4]*vfield[2 + 3*(i+1) + 3*Ngrid*j + 3*Ngrid*Ngrid*k] +
                    vint[5]*vfield[2 + 3*(i+1) + 3*Ngrid*j + 3*Ngrid*Ngrid*(k+1)] +
                    vint[6]*vfield[2 + 3*(i+1) + 3*Ngrid*(j+1) + 3*Ngrid*Ngrid*k] +
                    vint[7]*vfield[2 + 3*(i+1) + 3*Ngrid*(j+1) + 3*Ngrid*Ngrid*(k+1)];

    }
    //---------------------------------------------------------------------------------------------

  }
  //-----------------------------------------------------------------------------------------------

  return;
}


template class Ic<1>;
template class Ic<2>;
template class Ic<3>;
