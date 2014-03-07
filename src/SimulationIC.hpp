//=============================================================================
//  SimulationIC.hpp
//  Contains all routines for generating initial conditions on the fly.
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
//=============================================================================


#include <iostream>
#include <string>
#include <cstdio>
#include <cstring>
#include <math.h>
#include "Precision.h"
#include "Exception.h"
#include "Simulation.h"
#include "Sph.h"
#include "Parameters.h"
#include "InlineFuncs.h"
#include "Debug.h"
#include "Ghosts.h"
#if defined(FFTW_TURBULENCE)
#include "fftw3.h"
#endif
using namespace std;



//=============================================================================
//  Simulation::GenerateIC
/// Generate initial conditions for SPH simulation chosen in parameters file.
//=============================================================================
template <int ndim>
void Simulation<ndim>::GenerateIC(void)
{
  string in_file;                   // Restart snapshot filename
  string in_file_form;              // Restart snapshot file format
  string filename;                  // Simulation '.restart' filename
  ifstream f;                       // Stream of input file

  debug2("[Simulation::GenerateIC]");


  // First, check special case of restarting a simulation, in which case 
  // determine the name of the last snapshot file to be re-read
  //---------------------------------------------------------------------------
  if (restart) {
    filename = run_id + ".restart";
    f.open(filename.c_str());

    // If file opens successfully, read snapshot and return
    if (!f.fail()) {
      f >> in_file_form;
      f >> in_file;
      f.close();
      ReadSnapshotFile(in_file,in_file_form);
      ConvertToCodeUnits();
      return;
    }
  }


  // If not a restart, generate initial conditions either from external 
  // file or created on the fly.
  //---------------------------------------------------------------------------
  if (simparams->stringparams["ic"] == "file") {
    ReadSnapshotFile(simparams->stringparams["in_file"],
		     simparams->stringparams["in_file_form"]);
    rescale_particle_data = true;
  }
  //---------------------------------------------------------------------------
  else if (simparams->stringparams["ic"] == "binaryacc")
    BinaryAccretion();
  else if (simparams->stringparams["ic"] == "binary")
    BinaryStar();
  else if (simparams->stringparams["ic"] == "bb")
    BossBodenheimer();
  else if (simparams->stringparams["ic"] == "box")
    UniformBox();
  else if (simparams->stringparams["ic"] == "cdiscontinuity")
    ContactDiscontinuity();
  else if (simparams->stringparams["ic"] == "khi")
    KHI();
  else if (simparams->stringparams["ic"] == "noh")
    NohProblem();
  else if (simparams->stringparams["ic"] == "plummer")
    PlummerSphere();
  else if (simparams->stringparams["ic"] == "quadruple")
    QuadrupleStar();
  else if (simparams->stringparams["ic"] == "sedov")
    SedovBlastWave();
  else if (simparams->stringparams["ic"] == "shearflow")
    ShearFlow();
  else if (simparams->stringparams["ic"] == "shocktube")
    ShockTube();
  else if (simparams->stringparams["ic"] == "soundwave")
    SoundWave();
  else if (simparams->stringparams["ic"] == "sphere")
    UniformSphere();
  else if (simparams->stringparams["ic"] == "turbcore")
    TurbulentCore();
  else if (simparams->stringparams["ic"] == "triple")
    TripleStar();
  else if (simparams->stringparams["ic"] == "python")
    return;
  //---------------------------------------------------------------------------
  else {
    string message = "Unrecognised parameter : ic = " 
      + simparams->stringparams["ic"];
    ExceptionHandler::getIstance().raise(message);
  }
  //---------------------------------------------------------------------------

  // Scale particle data to dimensionless code units if required
  if (rescale_particle_data) ConvertToCodeUnits();

  // Check that the initial conditions are valid
  CheckInitialConditions();

  return;
}



//=============================================================================
//  Simulation::CheckInitialConditions
/// Performs some simple sanity checks on all initial conditions
//=============================================================================
template <int ndim>
void Simulation<ndim>::CheckInitialConditions(void)
{
  bool okflag;                      // Flag problem with current particle
  bool valid_ic = true;             // Valid initial conditions flag
  int i,k;                          // Particle and dimension counter
  SphParticle<ndim> *part;          // Pointer to SPH particle data


  // Check that all particles reside inside any defined boundaries
  //---------------------------------------------------------------------------
  for (i=0; i<sph->Nsph; i++) {
    part = sph->GetParticleIPointer(i);
    okflag = true;

    if (part->r[0] < simbox.boxmin[0])
      if (simbox.x_boundary_lhs == "periodic") okflag = false;
    if (part->r[0] > simbox.boxmax[0])
      if (simbox.x_boundary_rhs == "periodic") okflag = false;

    if (ndim >= 2 && part->r[1] < simbox.boxmin[1])
      if (simbox.y_boundary_lhs == "periodic") okflag = false;
    if (ndim >= 2 && part->r[1] > simbox.boxmax[1])
      if (simbox.y_boundary_rhs == "periodic") okflag = false;

    if (ndim == 3 && part->r[2] < simbox.boxmin[2])
      if (simbox.z_boundary_lhs == "periodic") okflag = false;
    if (ndim == 3 && part->r[2] > simbox.boxmax[2])
      if (simbox.z_boundary_rhs == "periodic") okflag = false;

    // If flag indicates a problem, print error and quit
    if (!okflag) {
      cout << "Particle " << i << " not inside periodic box" << endl;
      for (k=0; k<ndim; k++)
	cout << "r[" << k << "] : " << part->r[k] << "    " 
	     << simbox.boxmin[k] << "    " << simbox.boxmax[k] << endl;
    }

    valid_ic = okflag;

  }
  //---------------------------------------------------------------------------

  if (!valid_ic) {
    string message = "Invalid initial conditions for SPH particles";
    ExceptionHandler::getIstance().raise(message);
  }

  return;
}



//=============================================================================
//  Simulation::BinaryAccretion
/// Create initial conditions for binary accretion simulation.
//=============================================================================
template <int ndim>
void Simulation<ndim>::BinaryAccretion(void)
{
  int i;                            // Particle counter
  int j;                            // Aux. particle counter
  int k;                            // Dimension counter
  int Nbox1;                        // No. of particles in fluid 1
  int Nbox2;                        // No. of particles in fluid 2
  int Nlattice1[3];                 // Lattice dimensions for fluid 1
  int Nlattice2[3];                 // Lattice dimensions for fluid 2
  int Nneib;                        // Average no. of SPH neighbours
  FLOAT hfluid1;                    // Smoothing length of fluid 1
  FLOAT hsink;                      // Smoothing length of sink
  FLOAT rbinary[ndim];              // Initial position of binary COM
  FLOAT vbinary[ndim];              // Initial velocity of binary COM
  FLOAT rsonic;                     // Sonic radius
  FLOAT rsink;                      // Sink radius
  FLOAT volume1;                    // Volume of box1
  FLOAT volume2;                    // Volume of box1
  FLOAT *r1;                        // Positions for particles in fluid 1
  FLOAT *r2;                        // Positions for particles in fluid 2
  DomainBox<ndim> box1;             // Bounding box for fluid 1
  DomainBox<ndim> box2;             // Bounding box for fluid 2

  // Create local copies of initial conditions parameters
  int Nstar = simparams->intparams["Nstar"];
  FLOAT abin = simparams->floatparams["abin"];
  FLOAT ebin = simparams->floatparams["ebin"];
  FLOAT phirot = simparams->floatparams["phirot"];
  FLOAT thetarot = simparams->floatparams["thetarot"];
  FLOAT psirot = simparams->floatparams["psirot"];
  FLOAT vmachbin = simparams->floatparams["vmachbin"];
  FLOAT m1 = simparams->floatparams["m1"];
  FLOAT m2 = simparams->floatparams["m2"];
  FLOAT gammaone = simparams->floatparams["gamma_eos"] - 1.0;
  FLOAT rhofluid1 = simparams->floatparams["rhofluid1"];
  FLOAT rhofluid2 = simparams->floatparams["rhofluid2"];
  FLOAT press1 = simparams->floatparams["press1"];
  string particle_dist = simparams->stringparams["particle_distribution"];
  Nlattice1[0] = simparams->intparams["Nlattice1[0]"];
  Nlattice1[1] = simparams->intparams["Nlattice1[1]"];
  Nlattice1[2] = simparams->intparams["Nlattice1[2]"];
  Nlattice2[0] = simparams->intparams["Nlattice2[0]"];
  Nlattice2[1] = simparams->intparams["Nlattice2[1]"];
  Nlattice2[2] = simparams->intparams["Nlattice2[2]"];

  debug2("[Simulation::BinaryAccretion]");

  // Convert parameters to dimensionless units
  rhofluid1 /= simunits.rho.inscale;
  rhofluid2 /= simunits.rho.inscale;
  press1 /= simunits.press.inscale;

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
  //---------------------------------------------------------------------------
  if (ndim == 2) {
    volume1 = (box1.boxmax[0] - box1.boxmin[0])*
      (box1.boxmax[1] - box1.boxmin[1]);
    volume2 = (box2.boxmax[0] - box2.boxmin[0])*
      (box2.boxmax[1] - box2.boxmin[1]);
    Nneib = (int) (pi*pow(sph->kernp->kernrange*sph->h_fac,2));
    hfluid1 = sqrtf((volume1*(FLOAT) Nneib)/(4.0*(FLOAT) Nbox1));
  }
  else if (ndim == 3) {
    volume1 = (box1.boxmax[0] - box1.boxmin[0])*
      (box1.boxmax[1] - box1.boxmin[1])*
      (box1.boxmax[2] - box1.boxmin[2]);
    volume2 = (box2.boxmax[0] - box2.boxmin[0])*
      (box2.boxmax[1] - box2.boxmin[1])*
      (box2.boxmax[2] - box2.boxmin[2]);
    Nneib = (int) (pi*pow(sph->kernp->kernrange*sph->h_fac,2));
    hfluid1 = powf((3.0*volume1*(FLOAT) Nneib)/
 		   (32.0*pi*(FLOAT) Nbox1),onethird);
  }


  // Allocate main particle memory
  sph->Nsph = Nbox1 + Nbox2;
  nbody->Nstar = Nstar;
  AllocateParticleMemory();


  // Add a cube of random particles defined by the simulation bounding box and 
  // depending on the chosen particle distribution
  //---------------------------------------------------------------------------
  if (Nbox1 > 0) {
    r1 = new FLOAT[ndim*Nbox1];
    if (particle_dist == "random")
      AddRandomBox(Nbox1,r1,box1);
    else if (particle_dist == "cubic_lattice")
      AddCubicLattice(Nbox1,Nlattice1,r1,box1,true);
    else if (particle_dist == "hexagonal_lattice")
      AddHexagonalLattice(Nbox1,Nlattice1,r1,box1,true);
    else {
      string message = "Invalid particle distribution option";
      ExceptionHandler::getIstance().raise(message);
    }
    
    // Record positions in main memory
    for (i=0; i<Nbox1; i++) {
      for (k=0; k<ndim; k++) sph->sphdata[i].r[k] = r1[ndim*i + k];
      sph->sphdata[i].r[0] += 0.25*simbox.boxsize[0];
      if (sph->sphdata[i].r[0] > simbox.boxmax[0])
        sph->sphdata[i].r[0] -= simbox.boxsize[0];
      for (k=0; k<ndim; k++) sph->sphdata[i].v[k] = 0.0;
      sph->sphdata[i].m = rhofluid1*volume1/(FLOAT) Nbox1;
      sph->sphdata[i].h = sph->h_fac*pow(sph->sphdata[i].m/rhofluid1,invndim);
      sph->sphdata[i].u = press1/rhofluid1/gammaone;
    }
    delete[] r1;
  }


  // Add a cube of random particles defined by the simulation bounding box and 
  // depending on the chosen particle distribution
  //---------------------------------------------------------------------------
  if (Nbox2 > 0) {
    r2 = new FLOAT[ndim*Nbox2];
    if (particle_dist == "random")
      AddRandomBox(Nbox2,r2,box2);
    else if (particle_dist == "cubic_lattice")
      AddCubicLattice(Nbox2,Nlattice2,r2,box2,true);
    else if (particle_dist == "hexagonal_lattice")
      AddHexagonalLattice(Nbox2,Nlattice2,r2,box2,true);
    else {
      string message = "Invalid particle distribution option";
      ExceptionHandler::getIstance().raise(message);
    }
    
    // Record positions in main memory
    for (j=0; j<Nbox2; j++) {
      i = Nbox1 + j;
      for (k=0; k<ndim; k++) sph->sphdata[i].r[k] = r2[ndim*j + k];
      sph->sphdata[i].r[0] += 0.25*simbox.boxsize[0];
      if (sph->sphdata[i].r[0] > simbox.boxmax[0])
        sph->sphdata[i].r[0] -= simbox.boxsize[0];
      for (k=0; k<ndim; k++) sph->sphdata[i].v[k] = 0.0;
      sph->sphdata[i].h = sph->h_fac*pow(sph->sphdata[i].m/rhofluid2,invndim);
      sph->sphdata[i].m = rhofluid2*volume2/(FLOAT) Nbox1;
      sph->sphdata[i].u = press1/rhofluid2/gammaone;
    }
    delete[] r2;
  }

  initial_h_provided = true;

  rsonic = 0.5*m1/(press1/rhofluid1);
  //hsink = sph->kernp->invkernrange*rsonic;
  //hsink = min(hsink,hfluid1);
  hsink = hfluid1/pow(4.4817,invndim);
  rsink = sph->kernp->kernrange*hsink/pow(4.4817,invndim);
  FLOAT mmax = 0.5*4.0*rhofluid1*pow(sph->kernp->kernrange*hfluid1,ndim)/3.0;


  cout << "Sound speed : " << sqrt(press1/rhofluid1) << endl;
  cout << "rsonic      : " << rsonic << endl;
  cout << "rbondi      : " << 4.0*rsonic << endl;
  cout << "rsink       : " << rsink << endl;
  cout << "mmax        : " << mmax << endl;
  cout << "hfluid      : " << hfluid1 << endl;
  cout << "vbin        : " << vmachbin*sqrt(press1/rhofluid1) << endl;
  cout << "Bondi accretion, dmdt : " << 4.0*pi*rhofluid1*
    (m1 + m2)*(m1 + m2)/pow(press1/rhofluid1,1.5) << endl;

  // Set hmin_sink here, since no other sinks will be formed
  sph->hmin_sink = hsink;


  // Add star particles to simulation
  //---------------------------------------------------------------------------
  if (Nstar == 1) {
    for (k=0; k<ndim; k++) nbody->stardata[0].r[k] = 0.0;
    for (k=0; k<ndim; k++) nbody->stardata[0].v[k] = 0.0;
    if (vmachbin < small_number)
      nbody->stardata[0].r[0] = simbox.boxmin[0] + 0.5*simbox.boxsize[0];
    else
      nbody->stardata[0].r[0] = simbox.boxmin[0] + 0.25*simbox.boxsize[0];
    nbody->stardata[0].v[0] = vmachbin*sph->eos->SoundSpeed(sph->sphdata[0]);
    nbody->stardata[0].m = m1 + m2;
    nbody->stardata[0].h = hsink;
    nbody->stardata[0].radius = rsink;
    sinks.sink[0].star = &(nbody->stardata[0]);
    sinks.sink[0].radius = rsink;
    sinks.sink[0].mmax = mmax;
    sinks.Nsink = Nstar;
  }
  else if (Nstar == 2) {
    for (k=0; k<ndim; k++) rbinary[k] = 0.0;
    for (k=0; k<ndim; k++) vbinary[k] = 0.0;
    if (vmachbin < small_number)
      rbinary[0] = simbox.boxmin[0] + 0.5*simbox.boxsize[0];
    else
      rbinary[0] = simbox.boxmin[0] + 0.25*simbox.boxsize[0];
    vbinary[0] = vmachbin*sph->eos->SoundSpeed(sph->sphdata[0]);
    AddBinaryStar(abin,ebin,m1,m2,hsink,hsink,phirot,thetarot,psirot,0.0,
                  rbinary,vbinary,nbody->stardata[0],nbody->stardata[1]);
    sinks.sink[0].star = &(nbody->stardata[0]);
    sinks.sink[1].star = &(nbody->stardata[1]);
    sinks.sink[0].radius = rsink;
    sinks.sink[1].radius = rsink;
    sinks.sink[0].mmax = mmax;
    sinks.sink[1].mmax = mmax;
    sinks.Nsink = Nstar;
  }
  else {
    string message = "Invalid number of star particles";
    ExceptionHandler::getIstance().raise(message);
  }


  return;
}



//=============================================================================
//  Simulation::ShockTube
/// Generate 1D shock-tube test problem.
//=============================================================================
template <int ndim>
void Simulation<ndim>::ShockTube(void)
{
  int i;                            // Particle counter
  int j;                            // Aux. particle counter
  int k;                            // Dimension counter
  int Nbox1;                        // No. of particles in LHS box
  int Nbox2;                        // No. of particles in RHS box
  int Nlattice1[ndim];              // Particles per dimension for LHS lattice
  int Nlattice2[ndim];              // Particles per dimension for RHS lattice
  FLOAT dr[ndim];                   // Relative position vector
  FLOAT drmag;                      // Distance
  FLOAT drsqd;                      // Distance squared
  FLOAT volume;                     // Volume of box
  FLOAT vfluid1[ndim];              // Velocity vector of LHS fluid
  FLOAT vfluid2[ndim];              // Velocity vector of RHS fluid
  FLOAT wnorm;                      // Kernel normalisation
  FLOAT *r;                         // Position vectors
  FLOAT *uaux;                      // Temp. array for internal energy
  FLOAT *vaux;                      // Temp. array for x-velocities
  DomainBox<ndim> box1;             // LHS box
  DomainBox<ndim> box2;             // RHS box

  // Set local copies of various input parameters for setting-up test
  FLOAT rhofluid1 = simparams->floatparams["rhofluid1"];
  FLOAT rhofluid2 = simparams->floatparams["rhofluid2"];
  FLOAT press1 = simparams->floatparams["press1"];
  FLOAT press2 = simparams->floatparams["press2"];
  FLOAT temp0 = simparams->floatparams["temp0"];
  FLOAT mu_bar = simparams->floatparams["mu_bar"];
  FLOAT gammaone = simparams->floatparams["gamma_eos"] - 1.0;
  Nlattice1[0] = simparams->intparams["Nlattice1[0]"];
  Nlattice2[0] = simparams->intparams["Nlattice2[0]"];
  vfluid1[0] = simparams->floatparams["vfluid1[0]"];
  vfluid2[0] = simparams->floatparams["vfluid2[0]"];

  debug2("[Simulation::ShockTube]");

  if (ndim != 1) {
    string message = "Wrong dimensionality : " + ndim;
    ExceptionHandler::getIstance().raise(message);
  }

  // Compute size and range of fluid bounding boxes
  //---------------------------------------------------------------------------
  if (ndim == 1) {
    box1.boxmin[0] = simbox.boxmin[0];
    box1.boxmax[0] = 0.0;
    box2.boxmin[0] = 0.0;
    box2.boxmax[0] = simbox.boxmax[0];
    volume = box1.boxmax[0] - box1.boxmin[0];
    Nbox1 = Nlattice1[0];
    Nbox2 = Nlattice2[0];
  }

  // Allocate local and main particle memory
  sph->Nsph = Nbox1 + Nbox2;
  AllocateParticleMemory();
  r = new FLOAT[ndim*sph->Nsph];


  // Add particles for LHS of the shocktube
  //---------------------------------------------------------------------------
  if (Nbox1 > 0) {
    AddCubicLattice(Nbox1,Nlattice1,r,box1,false);

    for (i=0; i<Nbox1; i++) {
      for (k=0; k<ndim; k++) sph->sphdata[i].r[k] = r[ndim*i + k];
      for (k=0; k<ndim; k++) sph->sphdata[i].v[k] = 0.0;
      sph->sphdata[i].v[0] = vfluid1[0];
      sph->sphdata[i].m = rhofluid1*volume/(FLOAT) Nbox1;
      sph->sphdata[i].h = sph->h_fac*pow(sph->sphdata[i].m/rhofluid1,invndim);
      if (sph->gas_eos == "isothermal")
    	sph->sphdata[i].u = temp0/gammaone/mu_bar;
      else
        sph->sphdata[i].u = press1/rhofluid1/gammaone;
    }
  }

  // Add particles for RHS of the shocktube
  //---------------------------------------------------------------------------
  if (Nbox2 > 0) {
    AddCubicLattice(Nbox2,Nlattice2,r,box2,false);

    for (j=0; j<Nbox2; j++) {
      i = Nbox1 + j;
      for (k=0; k<ndim; k++) sph->sphdata[i].r[k] = r[ndim*j + k];
      for (k=0; k<ndim; k++) sph->sphdata[i].v[k] = 0.0;
      sph->sphdata[i].v[0] = vfluid2[0];
      sph->sphdata[i].m = rhofluid2*volume/(FLOAT) Nbox2;
      sph->sphdata[i].h = sph->h_fac*pow(sph->sphdata[i].m/rhofluid2,invndim);
      if (sph->gas_eos == "isothermal")
        sph->sphdata[i].u = temp0/gammaone/mu_bar;
      else
        sph->sphdata[i].u = press2/rhofluid2/gammaone;
    }
  }

  initial_h_provided = true;
  bool smooth_ic = true;

  // Smooth the initial conditions
  //---------------------------------------------------------------------------
  if (smooth_ic) {

    // Set initial smoothing lengths and create initial ghost particles
    //-------------------------------------------------------------------------
    sph->Nghost = 0;
    sph->Nghostmax = sph->Nsphmax - sph->Nsph;
    sph->Ntot = sph->Nsph;
    for (i=0; i<sph->Nsph; i++) sph->sphdata[i].active = true;
    
    //sph->InitialSmoothingLengthGuess();
    sphneib->BuildTree(rebuild_tree,n,ntreebuildstep,ntreestockstep,timestep,sph);
    
    // Search ghost particles
    LocalGhosts->SearchGhostParticles(0.0,simbox,sph);
    
    // Update neighbour tree
    sphneib->BuildTree(rebuild_tree,n,ntreebuildstep,ntreestockstep,timestep,sph);
    
    // Calculate all SPH properties
    sphneib->UpdateAllSphProperties(sph,nbody);
    
    sphneib->BuildTree(rebuild_tree,n,ntreebuildstep,ntreestockstep,timestep,sph);
    sphneib->UpdateAllSphProperties(sph,nbody);
    
    LocalGhosts->CopySphDataToGhosts(simbox,sph);
    
    // Calculate all SPH properties
    sphneib->UpdateAllSphProperties(sph,nbody);
    

    uaux = new FLOAT[sph->Nsph];
    vaux = new FLOAT[sph->Nsph*ndim];
    
    // Now compute smoothed quantities
    for (i=0; i<sph->Nsph; i++) {
      uaux[i] = 0.0;
      for (k=0; k<ndim; k++) vaux[ndim*i + k] = 0.0;
      wnorm = 0.0;
      for (j=0; j<sph->Ntot; j++) {
        for (k=0; k<ndim; k++)
          dr[k] = sph->sphdata[j].r[k] - sph->sphdata[i].r[k];
        drsqd = DotProduct(dr,dr,ndim);
        if (drsqd > pow(sph->kernp->kernrange*sph->sphdata[i].h,2)) continue;
        drmag = sqrt(drsqd);
        uaux[i] += sph->sphdata[j].m*sph->sphdata[j].u*
          sph->kernp->w0(drmag*sph->sphdata[i].invh)*
          pow(sph->sphdata[i].invh,ndim)*sph->sphdata[i].invrho;
        for (k=0; k<ndim; k++) vaux[ndim*i + k] +=
          sph->sphdata[j].m*sph->sphdata[j].v[k]*
          sph->kernp->w0(drmag*sph->sphdata[i].invh)*
          pow(sph->sphdata[i].invh,ndim)*sph->sphdata[i].invrho;
        wnorm += sph->sphdata[j].m*sph->kernp->w0(drmag*sph->sphdata[i].invh)*
          pow(sph->sphdata[i].invh,ndim)*sph->sphdata[i].invrho;
      }
      uaux[i] /= wnorm;
      for (k=0; k<ndim; k++) vaux[ndim*i + k] /= wnorm;
    }

    for (i=0; i<sph->Nsph; i++) {
      sph->sphdata[i].u = uaux[i];
      for (k=0; k<ndim; k++) sph->sphdata[i].v[k] = vaux[ndim*i + k];
    }

    delete[] vaux;
    delete[] uaux;

  }

  delete[] r;

  return;
}



//=============================================================================
//  Simulation::UniformBox
/// Populate the simulation bounding box with random particles.
//=============================================================================
template <int ndim>
void Simulation<ndim>::UniformBox(void)
{
  int i,k;                          // Particle and dimension counters
  int Nbox;                         // No. of particles in box
  int Nlattice[3];                  // Particles per dimension for LHS lattice
  FLOAT volume;                     // Volume of box
  FLOAT *r;                         // Position vectors of all particles

  // Local copy of important parameters
  string particle_dist = simparams->stringparams["particle_distribution"];
  int Npart = simparams->intparams["Nsph"];
  //FLOAT rhobox = simparams->intparams["rhofluid1"];
  Nlattice[0] = simparams->intparams["Nlattice1[0]"];
  Nlattice[1] = simparams->intparams["Nlattice1[1]"];
  Nlattice[2] = simparams->intparams["Nlattice1[2]"];

  debug2("[Simulation::UniformBox]");

  // Compute volume and number of particles inside box
  if (ndim == 1) {
    volume = simbox.boxmax[0] - simbox.boxmin[0];
    Nbox = Nlattice[0];
  }
  else if (ndim == 2) {
    volume = (simbox.boxmax[0] - simbox.boxmin[0])*
      (simbox.boxmax[1] - simbox.boxmin[1]);
    Nbox = Nlattice[0]*Nlattice[1];
  }
  else if (ndim == 3) {
    volume = (simbox.boxmax[0] - simbox.boxmin[0])*
      (simbox.boxmax[1] - simbox.boxmin[1])*
      (simbox.boxmax[2] - simbox.boxmin[2]);
    Nbox = Nlattice[0]*Nlattice[1]*Nlattice[2];
  }

  // Add a cube of random particles defined by the simulation bounding box and 
  // depending on the chosen particle distribution
  if (particle_dist == "random") {
    r = new FLOAT[ndim*Npart]; 
    AddRandomBox(Npart,r,simbox);
  }
  else if (particle_dist == "cubic_lattice") {
    Npart = Nbox;
    r = new FLOAT[ndim*Npart];
    AddCubicLattice(Npart,Nlattice,r,simbox,true);
  }
  else if (particle_dist == "hexagonal_lattice") {
    Npart = Nbox;
    r = new FLOAT[ndim*Npart];
    AddHexagonalLattice(Npart,Nlattice,r,simbox,true);
  }
  else {
    string message = "Invalid particle distribution option";
    ExceptionHandler::getIstance().raise(message);
  }

  // Allocate global and local memory for all particles
  sph->Nsph = Npart;
  AllocateParticleMemory();

  // Copy positions to main array and initialise all other variables
  for (i=0; i<sph->Nsph; i++) {
    for (k=0; k<ndim; k++) {
      sph->sphdata[i].r[k] = r[ndim*i + k];
      sph->sphdata[i].v[k] = (FLOAT) 0.0;
      sph->sphdata[i].a[k] = (FLOAT) 0.0;
    }
    sph->sphdata[i].m = volume/ (FLOAT) sph->Nsph;
    sph->sphdata[i].h = sph->h_fac*pow(volume / (FLOAT) sph->Nsph,invndim);
    sph->sphdata[i].invomega = (FLOAT) 1.0;
    sph->sphdata[i].iorig = i;
    sph->sphdata[i].u = (FLOAT) 1.5;
  }

  initial_h_provided = true;

  delete[] r;

  return;
}



//=============================================================================
//  Simulation::UniformSphere
/// Create a uniform-density sphere of particles of given origin and radius.
//=============================================================================
template <int ndim>
void Simulation<ndim>::UniformSphere(void)
{
  int i,k;                          // Particle and dimension counters
  int Nsphere;                      // Actual number of particles in sphere
  FLOAT rcentre[ndim];              // Position of sphere centre
  FLOAT volume;                     // Volume of sphere
  FLOAT *r;                         // Particle position vectors

  // Local copies of important parameters
  int Npart = simparams->intparams["Nsph"];
  FLOAT radius = simparams->floatparams["radius"];
  FLOAT rhofluid = simparams->floatparams["rhofluid1"];
  FLOAT press = simparams->floatparams["press1"];
  FLOAT gammaone = simparams->floatparams["gamma_eos"] - 1.0;
  string particle_dist = simparams->stringparams["particle_distribution"];

  debug2("[Simulation::UniformSphere]");

  r = new FLOAT[ndim*Npart];

  // Add a sphere of random particles with origin 'rcentre' and radius 'radius'
  for (k=0; k<ndim; k++) rcentre[k] = (FLOAT) 0.0;

  // Create the sphere depending on the choice of initial particle distribution
  if (particle_dist == "random")
    AddRandomSphere(Npart,r,rcentre,radius);
  else if (particle_dist == "cubic_lattice" || 
	   particle_dist == "hexagonal_lattice") {
    Nsphere = AddLatticeSphere(Npart,r,rcentre,radius,particle_dist);
    if (Nsphere != Npart) 
      cout << "Warning! Unable to converge to required " 
	   << "no. of ptcls due to lattice symmetry" << endl;
    Npart = Nsphere;
  }
  else {
    string message = "Invalid particle distribution option";
    ExceptionHandler::getIstance().raise(message);
  }

  sph->Nsph = Npart;
  AllocateParticleMemory();

  if (ndim == 1) volume = 2.0*radius;
  else if (ndim == 2) volume = pi*radius*radius;
  else if (ndim == 3) volume = 4.0*onethird*pi*pow(radius,3);

  // Record particle positions and initialise all other variables
#pragma omp parallel for default(none)\
  shared(gammaone,Npart,press,r,rhofluid,volume) private(i,k)
  for (i=0; i<sph->Nsph; i++) {
    for (k=0; k<ndim; k++) {
      sph->sphdata[i].r[k] = r[ndim*i + k];
      sph->sphdata[i].v[k] = (FLOAT) 0.0;
      sph->sphdata[i].a[k] = (FLOAT) 0.0;
    }
    sph->sphdata[i].m = rhofluid*volume / (FLOAT) Npart;
    sph->sphdata[i].h = sph->h_fac*pow(sph->sphdata[i].m/rhofluid,invndim);
    sph->sphdata[i].u = press/rhofluid/gammaone;
    sph->sphdata[i].invomega = (FLOAT) 1.0;
    sph->sphdata[i].zeta = (FLOAT) 0.0;
    sph->sphdata[i].iorig = i;
  }

  initial_h_provided = true;

  delete[] r;

  return;
}



//=============================================================================
//  Simulation::ContactDiscontinuity
/// Set-up contact discontinuity problem.
//=============================================================================
template <int ndim>
void Simulation<ndim>::ContactDiscontinuity(void)
{
  int i;                            // Particle counter
  int j;                            // Aux. particle counter
  int Nbox1;                        // No. of particles in LHS box
  int Nbox2;                        // No. of particles in RHS box
  int Nlattice1[ndim];              // Particles per dimension for LHS lattice
  int Nlattice2[ndim];              // Particles per dimension for RHS lattice
  FLOAT volume;                     // Volume of box
  FLOAT vfluid1[ndim];              // Velocity vector of LHS fluid
  FLOAT vfluid2[ndim];              // Velocity vector of RHS fluid
  FLOAT *r;                         // Position vectors
  DomainBox<ndim> box1;             // LHS box
  DomainBox<ndim> box2;             // RHS box

  // Create local copies of all parameters required to set-up problem
  FLOAT rhofluid1 = simparams->floatparams["rhofluid1"];
  FLOAT rhofluid2 = simparams->floatparams["rhofluid2"];
  FLOAT press1 = simparams->floatparams["press1"];
  FLOAT press2 = simparams->floatparams["press2"];
  FLOAT temp0 = simparams->floatparams["temp0"];
  FLOAT mu_bar = simparams->floatparams["mu_bar"];
  FLOAT gammaone = simparams->floatparams["gamma_eos"] - 1.0;
  Nlattice1[0] = simparams->intparams["Nlattice1[0]"];
  Nlattice1[1] = simparams->intparams["Nlattice1[1]"];
  Nlattice2[0] = simparams->intparams["Nlattice2[0]"];
  Nlattice2[1] = simparams->intparams["Nlattice2[1]"];
  vfluid1[0] = simparams->floatparams["vfluid1[0]"];
  vfluid2[0] = simparams->floatparams["vfluid2[0]"];

  debug2("[Simulation::ContactDiscontinuity]");


  // 1D simulation
  //===========================================================================
  if (ndim == 1) {
    box1.boxmin[0] = simbox.boxmin[0];
    box1.boxmax[0] = 0.8*simbox.boxmax[0];
    box2.boxmin[0] = 0.8*simbox.boxmax[0];
    box2.boxmax[0] = simbox.boxmax[0];
    volume = box1.boxmax[0] - box1.boxmin[0];
    Nbox1 = Nlattice1[0];
    Nbox2 = Nlattice2[0];

    // Allocate local and main particle memory
    sph->Nsph = Nbox1 + Nbox2;
    AllocateParticleMemory();
    r = new FLOAT[ndim*sph->Nsph];


    //-------------------------------------------------------------------------
    if (Nbox1 > 0) {
      AddCubicLattice(Nbox1,Nlattice1,r,box1,false);
      volume = box1.boxmax[0] - box1.boxmin[0];
      for (i=0; i<Nbox1; i++) {
        sph->sphdata[i].r[0] = r[i] - 0.4*simbox.boxsize[0];
        if (sph->sphdata[i].r[0] < simbox.boxmin[0])
          sph->sphdata[i].r[0] += simbox.boxsize[0];
        sph->sphdata[i].v[0] = 0.0;
        sph->sphdata[i].m = rhofluid1*volume/(FLOAT) Nbox1;
	sph->sphdata[i].h = sph->h_fac*pow(sph->sphdata[i].m/rhofluid1,invndim);
        if (sph->gas_eos == "isothermal")
          sph->sphdata[i].u = temp0/gammaone/mu_bar;
        else
          sph->sphdata[i].u = press1/rhofluid1/gammaone;
      }
    }

    //-------------------------------------------------------------------------
    if (Nbox2 > 0) {
      AddCubicLattice(Nbox2,Nlattice2,r,box2,false);
      volume = box2.boxmax[0] - box2.boxmin[0];
      for (j=0; j<Nbox2; j++) {
        i = Nbox1 + j;
        sph->sphdata[i].r[0] = r[j] - 0.4*simbox.boxsize[0];
        if (sph->sphdata[i].r[0] < simbox.boxmin[0])
          sph->sphdata[i].r[0] += simbox.boxsize[0];
        sph->sphdata[i].v[0] = 0.0;
        sph->sphdata[i].m = rhofluid2*volume/(FLOAT) Nbox2;
	sph->sphdata[i].h = sph->h_fac*pow(sph->sphdata[i].m/rhofluid2,invndim);
        if (sph->gas_eos == "isothermal")
          sph->sphdata[i].u = temp0/gammaone/mu_bar;
        else
          sph->sphdata[i].u = press2/rhofluid2/gammaone;
      }
    }

  }
  //===========================================================================
  else if (ndim == 2) {



  }
  //===========================================================================


  // Set initial smoothing lengths and create initial ghost particles
  //---------------------------------------------------------------------------
  sph->Nghost = 0;
  sph->Nghostmax = sph->Nsphmax - sph->Nsph;
  sph->Ntot = sph->Nsph;
  for (int i=0; i<sph->Nsph; i++) sph->sphdata[i].active = true;

  initial_h_provided = true;
  sphneib->BuildTree(rebuild_tree,n,ntreebuildstep,ntreestockstep,timestep,sph);

  // Search ghost particles
  LocalGhosts->SearchGhostParticles(0.0,simbox,sph);


  // Update neighbour tre
  sphneib->BuildTree(rebuild_tree,n,ntreebuildstep,ntreestockstep,timestep,sph);

  // Calculate all SPH properties
  sphneib->UpdateAllSphProperties(sph,nbody);

  sphneib->BuildTree(rebuild_tree,n,ntreebuildstep,ntreestockstep,timestep,sph);
  sphneib->UpdateAllSphProperties(sph,nbody);

  LocalGhosts->CopySphDataToGhosts(simbox,sph);

  // Calculate all SPH properties
  sphneib->UpdateAllSphProperties(sph,nbody);

  delete[] r;

  return;
}



//=============================================================================
//  Simulation::KHI
/// Set-up 2D Kelvin-Helmholtz instability test.
//=============================================================================
template <int ndim>
void Simulation<ndim>::KHI(void)
{
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
  FLOAT press1 = simparams->floatparams["press1"];
  FLOAT press2 = simparams->floatparams["press2"];
  FLOAT gammaone = simparams->floatparams["gamma_eos"] - 1.0;
  FLOAT amp = simparams->floatparams["amp"];
  FLOAT lambda = simparams->floatparams["lambda"];
  Nlattice1[0] = simparams->intparams["Nlattice1[0]"];
  Nlattice1[1] = simparams->intparams["Nlattice1[1]"];
  Nlattice2[0] = simparams->intparams["Nlattice2[0]"];
  Nlattice2[1] = simparams->intparams["Nlattice2[1]"];
  vfluid1[0] = simparams->floatparams["vfluid1[0]"];
  vfluid2[0] = simparams->floatparams["vfluid2[0]"];

  debug2("[Simulation::ShockTube]");

  if (ndim != 2) {
    string message = "Kelvin-Helmholtz instability only in 2D";
    ExceptionHandler::getIstance().raise(message);
  }

  // Compute size and range of fluid bounding boxes
  //---------------------------------------------------------------------------
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
  sph->Nsph = Nbox1 + Nbox2;
  AllocateParticleMemory();
  r = new FLOAT[ndim*sph->Nsph];


  // Add particles for LHS of the shocktube
  //---------------------------------------------------------------------------
  if (Nbox1 > 0) {
    AddCubicLattice(Nbox1,Nlattice1,r,box1,false);

    for (i=0; i<Nbox1; i++) {
      for (k=0; k<ndim; k++) sph->sphdata[i].r[k] = r[ndim*i + k];
      for (k=0; k<ndim; k++) sph->sphdata[i].v[k] = 0.0;
      sph->sphdata[i].r[1] -= 0.25*simbox.boxsize[1];
      if (sph->sphdata[i].r[1] < simbox.boxmin[1]) 
        sph->sphdata[i].r[1] += simbox.boxsize[1];
      sph->sphdata[i].v[0] = vfluid1[0];
      sph->sphdata[i].m = rhofluid1*volume/(FLOAT) Nbox1;
      sph->sphdata[i].h = sph->h_fac*pow(sph->sphdata[i].m/rhofluid1,invndim);
      sph->sphdata[i].u = press1/rhofluid1/gammaone;
    }
  }

  // Add particles for RHS of the shocktube
  //---------------------------------------------------------------------------
  if (Nbox2 > 0) {
    AddCubicLattice(Nbox2,Nlattice2,r,box2,false);

    for (j=0; j<Nbox2; j++) {
      i = Nbox1 + j;
      for (k=0; k<ndim; k++) sph->sphdata[i].r[k] = r[ndim*j + k];
      for (k=0; k<ndim; k++) sph->sphdata[i].v[k] = 0.0;
      sph->sphdata[i].r[1] -= 0.25*simbox.boxsize[1];
      if (sph->sphdata[i].r[1] < simbox.boxmin[1]) 
        sph->sphdata[i].r[1] += simbox.boxsize[1];
      sph->sphdata[i].v[0] = vfluid2[0];
      sph->sphdata[i].m = rhofluid2*volume/(FLOAT) Nbox2;
      sph->sphdata[i].h = sph->h_fac*pow(sph->sphdata[i].m/rhofluid2,invndim);
      sph->sphdata[i].u = press2/rhofluid2/gammaone;
    }
  }

  // Add velocity perturbation here
  //---------------------------------------------------------------------------
  FLOAT sigmapert = 0.05/sqrt(2.0);
  for (i=0; i<sph->Nsph; i++) {
    sph->sphdata[i].v[1] = amp*sin(2.0*pi*sph->sphdata[i].r[0]/lambda)*
      (exp(-pow(sph->sphdata[i].r[1] + 0.25,2)/2.0/sigmapert/sigmapert) +  
       exp(-pow(sph->sphdata[i].r[1] - 0.25,2)/2.0/sigmapert/sigmapert));
  }

  // Set initial smoothing lengths and create initial ghost particles
  //---------------------------------------------------------------------------
  sph->Nghost = 0;
  sph->Nghostmax = sph->Nsphmax - sph->Nsph;
  sph->Ntot = sph->Nsph;
  for (i=0; i<sph->Nsph; i++) sph->sphdata[i].active = true;
  
  initial_h_provided = true;
  
  // Search ghost particles
  LocalGhosts->SearchGhostParticles(0.0,simbox,sph);

  // Update neighbour tree
  rebuild_tree = true;
  sphneib->BuildTree(rebuild_tree,n,ntreebuildstep,ntreestockstep,timestep,sph);

  // Calculate all SPH properties
  sphneib->UpdateAllSphProperties(sph,nbody);
  
  for (i=0; i<sph->Nsph; i++) 
    sph->sphdata[i].u = press1/sph->sphdata[i].rho/gammaone;

  delete[] r;

  return;
}



//=============================================================================
//  Simulation::NohProblem
/// Set-up Noh Problem initial conditions
//=============================================================================
template <int ndim>
void Simulation<ndim>::NohProblem(void)
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
  int Npart = simparams->intparams["Nsph"];
  FLOAT rhofluid = simparams->floatparams["rhofluid1"];
  FLOAT press = simparams->floatparams["press1"];
  FLOAT radius = simparams->floatparams["radius"];
  FLOAT gammaone = simparams->floatparams["gamma_eos"] - 1.0;
  string particle_dist = simparams->stringparams["particle_distribution"];

  debug2("[Simulation::NohProblem]");

  r = new FLOAT[ndim*Npart];

  // Add a sphere of random particles with origin 'rcentre' and radius 'radius'
  for (k=0; k<ndim; k++) rcentre[k] = (FLOAT) 0.0;

  // Create the sphere depending on the choice of initial particle distribution
  if (particle_dist == "random")
    AddRandomSphere(Npart,r,rcentre,radius);
  else if (particle_dist == "cubic_lattice" || 
	   particle_dist == "hexagonal_lattice") {
    Nsphere = AddLatticeSphere(Npart,r,rcentre,radius,particle_dist);
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
  sph->Nsph = Npart;
  AllocateParticleMemory();

  if (ndim == 1) volume = 2.0*radius;
  else if (ndim == 2) volume = pi*radius*radius;
  else if (ndim == 3) volume = 4.0*onethird*pi*pow(radius,3);

  // Record particle properties in main memory
  for (i=0; i<Npart; i++) {
    for (k=0; k<ndim; k++) sph->sphdata[i].r[k] = r[ndim*i + k];
    for (k=0; k<ndim; k++) dr[k] = r[ndim*i + k];
    drsqd = DotProduct(dr,dr,ndim);
    drmag = sqrt(drsqd) + small_number;
    for (k=0; k<ndim; k++) 
      sph->sphdata[i].v[k] = -1.0*dr[k]/drmag;
    sph->sphdata[i].m = rhofluid*volume/(FLOAT) Npart;
    sph->sphdata[i].h = sph->h_fac*pow(sph->sphdata[i].m/rhofluid,invndim);
    sph->sphdata[i].u = press/rhofluid/gammaone;
  }

  initial_h_provided = true;

  delete[] r;

  return;
}



//=============================================================================
//  Simulation::BossBodenheimer
/// Set-up Boss-Bodenheimer (1979) initial conditions for collapse of a 
/// rotating uniform sphere with an imposed m=2 azimuthal density perturbation.
//=============================================================================
template <int ndim>
void Simulation<ndim>::BossBodenheimer(void)
{
  int i;                            // Particle counter
  int k;                            // Dimension counter
  int Nsphere;                      // Actual number of particles in sphere
  FLOAT mp;                         // Mass of one particle
  FLOAT rcentre[ndim];              // Position of sphere centre
  FLOAT rho;                        // Fluid density
  FLOAT *r;                         // Positions of all particles
  FLOAT *v;                         // Velocities of all particles

  // Create local copies of initial conditions parameters
  int Npart = simparams->intparams["Nsph"];
  FLOAT amp = simparams->floatparams["amp"];
  FLOAT angvel = simparams->floatparams["angvel"];
  FLOAT mcloud = simparams->floatparams["mcloud"];
  FLOAT mu_bar = simparams->floatparams["mu_bar"];
  FLOAT press = simparams->floatparams["press1"];
  FLOAT radius = simparams->floatparams["radius"];
  FLOAT temp0 = simparams->floatparams["temp0"];
  FLOAT gammaone = simparams->floatparams["gamma_eos"] - 1.0;
  string particle_dist = simparams->stringparams["particle_distribution"];

  debug2("[Simulation::BossBodenheimer]");

  // Convert any parameters to code units
  angvel /= simunits.angvel.outscale;
  mcloud /= simunits.m.outscale;
  press /= simunits.press.outscale;
  radius /= simunits.r.outscale;
  temp0 /= simunits.temp.outscale;

  r = new FLOAT[ndim*Npart];
  v = new FLOAT[ndim*Npart];

  // Add a sphere of random particles with origin 'rcentre' and radius 'radius'
  for (k=0; k<ndim; k++) rcentre[k] = (FLOAT) 0.0;

  // Create the sphere depending on the choice of initial particle distribution
  if (particle_dist == "random")
    AddRandomSphere(Npart,r,rcentre,radius);
  else if (particle_dist == "cubic_lattice" || 
	   particle_dist == "hexagonal_lattice") {
    Nsphere = AddLatticeSphere(Npart,r,rcentre,radius,particle_dist);
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
  sph->Nsph = Npart;
  AllocateParticleMemory();
  mp = mcloud / (FLOAT) Npart;
  rho = 3.0*mcloud / (4.0*pi*pow(radius,3));

  // Perturb positions of particles in cloud
  AddAzimuthalDensityPerturbation(Npart,2,amp,rcentre,r);

  // Add solid-body rotational velocity field
  AddRotationalVelocityField(Npart,angvel,rcentre,r,v);

  // Record particle properties in main memory
#pragma omp parallel for default(none)\
  shared(gammaone,Npart,mp,mu_bar,r,rho,temp0,v) private(i,k)
  for (i=0; i<Npart; i++) {
    for (k=0; k<ndim; k++) sph->sphdata[i].r[k] = r[ndim*i + k];
    for (k=0; k<ndim; k++) sph->sphdata[i].v[k] = v[ndim*i + k];
    sph->sphdata[i].m = mp;
    sph->sphdata[i].h = sph->h_fac*pow(sph->sphdata[i].m/rho,invndim);
    sph->sphdata[i].u = temp0/gammaone/mu_bar;
    //if (sph->gas_eos == "isothermal" || sph->gas_eos == "barotropic")
    //  sph->sphdata[i].u = temp0/gammaone/mu_bar;
    //else
    //  sph->sphdata[i].u = press/rho/gammaone;
  }

  initial_h_provided = true;

  delete[] v;
  delete[] r;

  return;
}



//=============================================================================
//  Simulation::TurbulentCore
/// Set-up Boss-Bodenheimer (1979) initial conditions for collapse of a 
/// rotating uniform sphere with an imposed m=2 azimuthal density perturbation.
//=============================================================================
template <int ndim>
void Simulation<ndim>::TurbulentCore(void)
{
  int i;                            // Particle counter
  int j;                            // ..
  int k;                            // Dimension counter
  int kk;                           // ..
  int p;                            // ..
  int Nsphere;                      // Actual number of particles in sphere
  FLOAT dx[3];                      // ..
  FLOAT dxgrid;                     // ..
  FLOAT gpecloud;                   // ..
  FLOAT keturb;                     // ..
  FLOAT mp;                         // Mass of one particle
  FLOAT rcentre[ndim];              // Position of sphere centre
  FLOAT rmax[ndim];                 // ..
  FLOAT rmin[ndim];                 // ..
  FLOAT rho;                        // Fluid density
  FLOAT vint[8];                    // ..
  FLOAT xmin;                       // ..
  FLOAT vfactor;                    // ..
  FLOAT *r;                         // Positions of all particles
  FLOAT *v;                         // Velocities of all particles
  DOUBLE *vfield;                   // ..

  // Create local copies of initial conditions parameters
  int field_type = simparams->intparams["field_type"];
  int gridsize = simparams->intparams["gridsize"];
  int Npart = simparams->intparams["Nsph"];
  FLOAT alpha_turb = simparams->floatparams["alpha_turb"];
  FLOAT gammaone = simparams->floatparams["gamma_eos"] - 1.0;
  FLOAT mcloud = simparams->floatparams["mcloud"];
  FLOAT mu_bar = simparams->floatparams["mu_bar"];
  FLOAT power_turb = simparams->floatparams["power_turb"];
  FLOAT radius = simparams->floatparams["radius"];
  FLOAT temp0 = simparams->floatparams["temp0"];
  string particle_dist = simparams->stringparams["particle_distribution"];

  debug2("[Simulation::TurbulentCore]");

  // Convert any parameters to code units
  mcloud /= simunits.m.outscale;
  radius /= simunits.r.outscale;
  temp0 /= simunits.temp.outscale;

  // Calculate gravitational potential energy of uniform cloud
  gpecloud = 0.6*mcloud*mcloud/radius;

  r = new FLOAT[ndim*Npart];
  v = new FLOAT[ndim*Npart];

  // Add a sphere of random particles with origin 'rcentre' and radius 'radius'
  for (k=0; k<ndim; k++) rcentre[k] = (FLOAT) 0.0;
  
  // Create the sphere depending on the choice of initial particle distribution
  if (particle_dist == "random")
    AddRandomSphere(Npart,r,rcentre,radius);
  else if (particle_dist == "cubic_lattice" || 
	   particle_dist == "hexagonal_lattice") {
    Nsphere = AddLatticeSphere(Npart,r,rcentre,radius,particle_dist);
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
  sph->Nsph = Npart;
  AllocateParticleMemory();
  mp = mcloud / (FLOAT) Npart;
  rho = 3.0*mcloud / (4.0*pi*pow(radius,3));


  // Record particle properties in main memory
  for (i=0; i<Npart; i++) {
    for (k=0; k<ndim; k++) sph->sphdata[i].r[k] = r[ndim*i + k];
    for (k=0; k<ndim; k++) sph->sphdata[i].v[k] = v[ndim*i + k];
    sph->sphdata[i].m = mp;
    sph->sphdata[i].h = sph->h_fac*pow(sph->sphdata[i].m/rho,invndim);
    sph->sphdata[i].u = temp0/gammaone/mu_bar;
    //if (sph->gas_eos == "isothermal" || sph->gas_eos == "barotropic")
    //  sph->sphdata[i].u = temp0/gammaone/mu_bar;
    //else
    //  sph->sphdata[i].u = press/rho/gammaone;
  }

  initial_h_provided = true;


  // Generate turbulent velocity field for given power spectrum slope
#if defined(FFTW_TURBULENCE)
  vfield = new DOUBLE[ndim*gridsize*gridsize*gridsize];
  sph->SphBoundingBox(rmax,rmin,sph->Nsph);
  xmin = 9.9e20;
  dxgrid = 0.0;
  for (k=0; k<ndim; k++) {
    dxgrid = max(dxgrid,(rmax[k] - rmin[k]/(FLOAT) (gridsize - 1)));
    xmin = min(xmin,rmin[k]);
  }

  // Generate gridded velocity field
  GenerateTurbulentVelocityField(field_type,gridsize,power_turb,vfield);

  // Now interpolate velocity field onto particle positions
  //---------------------------------------------------------------------------
  for (p=0; p<sph->Nsph; p++) {
    for (kk=0; kk<ndim; kk++) dx[kk] = (sph->sphdata[p].r[kk] - xmin)/dxgrid;

    i = (int) dx[0];
    j = (int) dx[1];
    k = (int) dx[2];
    
    if (i > gridsize || j > gridsize || k > gridsize) {
      cout << "Grid too big!! : " << i << "    " << j << "    " << k 
	   << "   " << gridsize << endl;
      exit(0);
    }
    
    for (kk=0; kk<ndim; kk++) dx[kk] -= (int) dx[kk];
    
    // Interpolate to get more accurate velocities
    if (ndim == 3) {
      vint[0] = (1.0 - dx[0])*(1.0 - dx[1])*(1.0 - dx[2]);
      vint[1] = (1.0 - dx[0])*(1.0 - dx[1])*dx[2];
      vint[2] = (1.0 - dx[0])*dx[1]*(1.0 - dx[2]);
      vint[3] = (1.0 - dx[0])*dx[1]*dx[2];
      vint[4] = dx[0]*(1.0 - dx[1])*(1.0 - dx[2]);
      vint[5] = dx[0]*(1.0 - dx[1])*dx[2];
      vint[6] = dx[0]*dx[1]*(1.0 - dx[2]);
      vint[7] = dx[0]*dx[1]*dx[2];

      v[ndim*p] = 
	vint[0]*vfield[3*i + 3*gridsize*j + 3*gridsize*gridsize*k] + 
	vint[1]*vfield[3*i + 3*gridsize*j + 3*gridsize*gridsize*(k+1)] + 
	vint[2]*vfield[3*i + 3*gridsize*(j+1) + 3*gridsize*gridsize*k] + 
	vint[3]*vfield[3*i + 3*gridsize*(j+1) + 3*gridsize*gridsize*(k+1)] + 
	vint[4]*vfield[3*(i+1) + 3*gridsize*j + 3*gridsize*gridsize*k] + 
	vint[5]*vfield[3*(i+1) + 3*gridsize*j + 3*gridsize*gridsize*(k+1)] + 
	vint[6]*vfield[3*(i+1) + 3*gridsize*(j+1) + 3*gridsize*gridsize*k] + 
	vint[7]*vfield[3*(i+1) + 3*gridsize*(j+1) + 3*gridsize*gridsize*(k+1)];

      v[ndim*p+1] = 
	vint[0]*vfield[1 + 3*i + 3*gridsize*j + 3*gridsize*gridsize*k] + 
	vint[1]*vfield[1 + 3*i + 3*gridsize*j + 3*gridsize*gridsize*(k+1)] + 
	vint[2]*vfield[1 + 3*i + 3*gridsize*(j+1) + 3*gridsize*gridsize*k] + 
	vint[3]*vfield[1 + 3*i + 3*gridsize*(j+1) + 3*gridsize*gridsize*(k+1)] + 
	vint[4]*vfield[1 + 3*(i+1) + 3*gridsize*j + 3*gridsize*gridsize*k] + 
	vint[5]*vfield[1 + 3*(i+1) + 3*gridsize*j + 3*gridsize*gridsize*(k+1)] + 
	vint[6]*vfield[1 + 3*(i+1) + 3*gridsize*(j+1) + 3*gridsize*gridsize*k] + 
	vint[7]*vfield[1 + 3*(i+1) + 3*gridsize*(j+1) + 3*gridsize*gridsize*(k+1)];

      v[ndim*p+2] = 
	vint[0]*vfield[2 + 3*i + 3*gridsize*j + 3*gridsize*gridsize*k] + 
	vint[1]*vfield[2 + 3*i + 3*gridsize*j + 3*gridsize*gridsize*(k+1)] + 
	vint[2]*vfield[2 + 3*i + 3*gridsize*(j+1) + 3*gridsize*gridsize*k] + 
	vint[3]*vfield[2 + 3*i + 3*gridsize*(j+1) + 3*gridsize*gridsize*(k+1)] + 
	vint[4]*vfield[2 + 3*(i+1) + 3*gridsize*j + 3*gridsize*gridsize*k] + 
	vint[5]*vfield[2 + 3*(i+1) + 3*gridsize*j + 3*gridsize*gridsize*(k+1)] + 
	vint[6]*vfield[2 + 3*(i+1) + 3*gridsize*(j+1) + 3*gridsize*gridsize*k] + 
	vint[7]*vfield[2 + 3*(i+1) + 3*gridsize*(j+1) + 3*gridsize*gridsize*(k+1)];

      for (kk=0; kk<ndim; kk++) sph->sphdata[p].v[kk] = v[ndim*p + kk];
      //cout << "ijk : " << i << "   " << j << "    " << k << endl;
      //cout << "dx2  : " << dx[0] << "     " << dx[1] << "    " << dx[2] << endl;
      //cout << "Part position : " << p << "    " << sph->sphdata[p].r[0] << "    " 
      //     << sph->sphdata[p].r[1] << "   " << sph->sphdata[p].r[2] << endl;
      //cout << "Part velocity : " << p << "    " << sph->sphdata[p].v[0] << "    " 
      //     << sph->sphdata[p].v[1] << "   " << sph->sphdata[p].v[2] << endl;
	
    }
  }
  //---------------------------------------------------------------------------


  cout << "xmin : " << xmin << "     " << dxgrid << "    gridsize : " << gridsize << endl;
  cout << "rmin : " << rmin[0] << "    " << rmin[1] << "    " << rmin[2] << endl;
  cout << "rmax : " << rmax[0] << "    " << rmax[1] << "    " << rmax[2] << endl;


#else
  string message = "FFTW turbulence flag not set";
  ExceptionHandler::getIstance().raise(message);
#endif


  // Change to COM frame of reference
  SetComFrame();

  // Calculate total kinetic energy of turbulent velocity field
  keturb = 0.0;
  for (i=0; i<sph->Nsph; i++) {
    keturb += 
      sph->sphdata[i].m*DotProduct(sph->sphdata[i].v,sph->sphdata[i].v,ndim);
  }
  keturb *= 0.5;

  vfactor = sqrt(alpha_turb*gpecloud/keturb);
  cout << "Scaling factor : " << vfactor << endl;

  // Now rescale velocities to give required turbulent energy in cloud
  for (i=0; i<sph->Nsph; i++) {
    for (k=0; k<ndim; k++) sph->sphdata[i].v[k] *= vfactor;
  }


  delete[] v;
  delete[] r;

  return;
}



//=============================================================================
//  Simulation::PlummerSphere
/// Generate a Plummer sphere containing either stars, gas or a mixture of 
/// both.  Uses the algorithm described by Aarseth et al. (197?).
//=============================================================================
template <int ndim>
void Simulation<ndim>::PlummerSphere(void)
{
  bool flag;                        // Aux. flag
  int i,j,k;                        // Particle and dimension counter
  FLOAT raux;                       // Aux. float variable
  FLOAT rcentre[ndim];              // Position of centre of Plummer sphere
  FLOAT vplummer;                   // ..
  FLOAT x1,x2,x3,x4,x5,x6,x7;       // ..
  FLOAT rad,vm,ve,t1,t2,w,z;        // ..

  // Local copies of important parameters
  int Nsph = simparams->intparams["Nsph"];
  int Nstar = simparams->intparams["Nstar"];
  FLOAT gamma_eos = simparams->floatparams["gamma_eos"];
  FLOAT gasfrac = simparams->floatparams["gasfrac"];
  FLOAT starfrac = simparams->floatparams["starfrac"];
  FLOAT mplummer = simparams->floatparams["mplummer"];
  FLOAT rplummer = simparams->floatparams["rplummer"];
  FLOAT radius = simparams->floatparams["radius"];
  FLOAT rstar = simparams->floatparams["rstar"];

  debug1("[Simulation::PlummerSphere]");
    
  sph->Nsph = Nsph;
  sph->Ntot = Nsph;
  nbody->Nstar = Nstar;
  AllocateParticleMemory();

  for (k=0; k<ndim; k++) rcentre[k] = 0.0;
  raux = gasfrac + starfrac;
  gasfrac /= raux;
  starfrac /= raux;
    
  // Loop over all particles (gas and stars)
  //===========================================================================
  for (j=0; j<Nsph+Nstar; j++) {

    do {
      flag = false;
      x1 = (FLOAT)(rand()%RAND_MAX)/(FLOAT)RAND_MAX;
      x2 = (FLOAT)(rand()%RAND_MAX)/(FLOAT)RAND_MAX;
      x3 = (FLOAT)(rand()%RAND_MAX)/(FLOAT)RAND_MAX;

      if (x1 == 0.0 && x2 == 0.0 && x3 == 0.0) flag = true;
      rad = 1.0 / sqrt(pow(x1,-2.0/3.0) - 1.0);
      if (rad > radius/rplummer) flag = true;
    } while (flag);

    z = (1.0 - 2.0*x2)*rad;


    // Set position depending on particle type
    //-------------------------------------------------------------------------
    if (j >= Nstar && j < Nstar + Nsph) {
      i = j - Nstar;
      sph->sphdata[i].r[2] = z;
      sph->sphdata[i].r[0] = sqrt(rad*rad - z*z)*cos(twopi*x3);
      sph->sphdata[i].r[1] = sqrt(rad*rad - z*z)*sin(twopi*x3);
      sph->sphdata[i].m = gasfrac / (FLOAT) Nsph;
    }
    else {
      i = j;
      nbody->stardata[i].r[2] = z;
      nbody->stardata[i].r[0] = sqrt(rad*rad - z*z)*cos(twopi*x3);
      nbody->stardata[i].r[1] = sqrt(rad*rad - z*z)*sin(twopi*x3);
      nbody->stardata[i].m = starfrac / (FLOAT) Nstar;
    }
       
    // Maximum velocity for this distance 
    ve = sqrt(2.0 / sqrt(1.0 + rad*rad));


    // Velocity of particle
    //-------------------------------------------------------------------------
    do {
      x4 = (FLOAT)(rand()%RAND_MAX)/(FLOAT)RAND_MAX;
      x5 = (FLOAT)(rand()%RAND_MAX)/(FLOAT)RAND_MAX;
      t1 = 0.1*x5;
      t2 = x4*x4*pow(1.0 - x4*x4,3.5);
    } while (t1 > t2);

    vm = ve*x4;
    x6 = (FLOAT)(rand()%RAND_MAX)/(FLOAT)RAND_MAX;
    x7 = (FLOAT)(rand()%RAND_MAX)/(FLOAT)RAND_MAX;
    w = (1.0 - 2.0*x6)*vm;
       

    // Set velocity depending on particle type
    //-------------------------------------------------------------------------
    if (j >= Nstar && j < Nstar + Nsph) {
      i = j - Nstar;
      for(k=0; k<ndim; k++) sph->sphdata[i].v[k] = 0.0;
      sph->sphdata[i].sound = sqrt(0.1666666 / sqrt(1.0 + rad*rad));
      sph->sphdata[i].rho = 1.0;
      sph->sphdata[i].u = sph->sphdata[i].sound*
        sph->sphdata[i].sound/(gamma_eos - 1.0);
    }
    else {
      i = j;
      nbody->stardata[i].v[2] = w;
      nbody->stardata[i].v[0] = sqrt(vm*vm - w*w)*cos(twopi*x7);
      nbody->stardata[i].v[1] = sqrt(vm*vm - w*w)*sin(twopi*x7);
    }
      
  }
  //===========================================================================

  // Instanly move to COM
  //ConvertToComFrame();
  vplummer = sqrt(mplummer/rplummer);

  // Now scale variables to required physical size
  for (i=0; i<Nsph; i++) {
    for (k=0; k<ndim; k++) {
      sph->sphdata[i].r[k] = sph->sphdata[i].r[k]*rplummer;
      sph->sphdata[i].v[k] = sph->sphdata[i].v[k]*vplummer;
    }
    sph->sphdata[i].m    = sph->sphdata[i].m*mplummer;
    if (i < Nsph) sph->sphdata[i].u = sph->sphdata[i].u*(mplummer/rplummer);
  }
  for (i=0; i<Nstar; i++) {
    for (k=0; k<ndim; k++) {
      nbody->stardata[i].r[k] = nbody->stardata[i].r[k]*rplummer;
      nbody->stardata[i].v[k] = nbody->stardata[i].v[k]*vplummer;
    }
    nbody->stardata[i].m      = nbody->stardata[i].m*mplummer;
    //nbody->stardata[i].radius = rstar;
    nbody->stardata[i].h      = sph->kernp->invkernrange*rstar;
    nbody->stardata[i].invh   = 1.0 / nbody->stardata[i].h;
  }

  return;
}



//=============================================================================
//  Simulation::SedovBlastWave
/// Set-up Sedov blast wave test
//=============================================================================
template <int ndim>
void Simulation<ndim>::SedovBlastWave(void)
{
  int i;                            // Particle counter
  int k;                            // Dimension counter
  int Nbox;                         // No. of particles in box
  int Ncold;                        // No. of cold particles
  int Nhot;                         // No. of hot particles
  int Nlattice[3];                  // Lattice size
  int *hotlist;                     // List of 'hot' particles
  FLOAT drmag;                      // Distance
  FLOAT drsqd;                      // Distance squared
  FLOAT mbox;                       // Total mass inside simulation box
  FLOAT r_hot;                      // Size of 'hot' region
  FLOAT ufrac;                      // Internal energy fraction
  FLOAT umax;                       // Maximum u of all particles
  FLOAT utot;                       // Total internal energy
  FLOAT volume;                     // Volume of box
  FLOAT *r;                         // Positions of all particles

  // Create local copies of initial conditions parameters
  int smooth_ic = simparams->intparams["smooth_ic"];
  FLOAT rhofluid = simparams->floatparams["rhofluid1"];
  FLOAT kefrac = simparams->floatparams["kefrac"];
  string particle_dist = simparams->stringparams["particle_distribution"];
  Nlattice[0] = simparams->intparams["Nlattice1[0]"];
  Nlattice[1] = simparams->intparams["Nlattice1[1]"];
  Nlattice[2] = simparams->intparams["Nlattice1[2]"];

  debug2("[Simulation::SedovBlastWave]");


  // Compute size and range of fluid bounding boxes
  //---------------------------------------------------------------------------
  if (ndim == 1) {
    volume = simbox.boxmax[0] - simbox.boxmin[0];
    Nbox = Nlattice[0];
  }
  else if (ndim == 2) {
    volume = (simbox.boxmax[0] - simbox.boxmin[0])*
      (simbox.boxmax[1] - simbox.boxmin[1]);
    Nbox = Nlattice[0]*Nlattice[1];
  }
  else if (ndim == 3) {
    volume = (simbox.boxmax[0] - simbox.boxmin[0])*
      (simbox.boxmax[1] - simbox.boxmin[1])*
      (simbox.boxmax[2] - simbox.boxmin[2]);
    Nbox = Nlattice[0]*Nlattice[1]*Nlattice[2];
  }
  mbox = volume*rhofluid;
  ufrac = max((FLOAT) 0.0,(FLOAT) 1.0 - kefrac);
  Ncold = 0;
  Nhot = 0;
  r_hot = powf(powf(4.0,ndim)/(FLOAT) Nbox,sph->invndim);

  // Allocate local and main particle memory
  sph->Nsph = Nbox;
  AllocateParticleMemory();
  r = new FLOAT[ndim*sph->Nsph];
  hotlist = new int[sph->Nsph];

  // Add a cube of random particles defined by the simulation bounding box and 
  // depending on the chosen particle distribution
  if (particle_dist == "random")
    AddRandomBox(Nbox,r,simbox);
  else if (particle_dist == "cubic_lattice")
    AddCubicLattice(Nbox,Nlattice,r,simbox,true);
  else if (particle_dist == "hexagonal_lattice")
    AddHexagonalLattice(Nbox,Nlattice,r,simbox,true);
  else {
    string message = "Invalid particle distribution option";
    ExceptionHandler::getIstance().raise(message);
  }

  // Record positions in main memory
  for (i=0; i<Nbox; i++) {
    for (k=0; k<ndim; k++) sph->sphdata[i].r[k] = r[ndim*i + k];
    for (k=0; k<ndim; k++) sph->sphdata[i].v[k] = 0.0;
    sph->sphdata[i].m = mbox/(FLOAT) Nbox;
    sph->sphdata[i].h = sph->h_fac*pow(sph->sphdata[i].m/rhofluid,invndim);
    sph->sphdata[i].u = small_number;
  }

  // Set initial smoothing lengths and create initial ghost particles
  //---------------------------------------------------------------------------
  sph->Nghost = 0;
  sph->Nghostmax = sph->Nsphmax - sph->Nsph;
  sph->Ntot = sph->Nsph;
  for (i=0; i<sph->Nsph; i++) sph->sphdata[i].active = true;

  // Search ghost particles
  LocalGhosts->SearchGhostParticles(0.0,simbox,sph);
  
  initial_h_provided = true;
  rebuild_tree = true;
  sphneib->BuildTree(rebuild_tree,n,ntreebuildstep,ntreestockstep,timestep,sph);
  
  sphneib->UpdateAllSphProperties(sph,nbody);

  // Update neighbour tre
  rebuild_tree = true;
  sphneib->BuildTree(rebuild_tree,n,ntreebuildstep,ntreestockstep,timestep,sph);

  // Calculate all SPH properties
  sphneib->UpdateAllSphProperties(sph,nbody);

  // Now calculate which particles are hot
  //---------------------------------------------------------------------------
  umax = (FLOAT) 0.0;
  utot = (FLOAT) 0.0;
  for (i=0; i<sph->Nsph; i++) {
    drsqd = DotProduct(sph->sphdata[i].r,sph->sphdata[i].r,ndim);
    if (drsqd < r_hot*r_hot) {
      hotlist[i] = 1;
      if (smooth_ic == 1)
	sph->sphdata[i].u = sph->sphdata[i].m*
	  sph->kernp->w0(sph->kernp->kernrange*sqrt(drsqd)/r_hot);
      else
	sph->sphdata[i].u = sph->sphdata[i].m;
      utot += sph->sphdata[i].u;
      umax = max(umax,sph->sphdata[i].u);
      Nhot++;
    }
    else {
      hotlist[i] = 0;
      Ncold++;
    }
  }

  // Normalise the energies
  //---------------------------------------------------------------------------
  for (i=0; i<sph->Nsph; i++) {
    if (hotlist[i] == 1) {
      drmag = sqrt(DotProduct(sph->sphdata[i].r,sph->sphdata[i].r,ndim));
      sph->sphdata[i].u = sph->sphdata[i].u/utot/sph->sphdata[i].m;
      //,1.0e-6*umax/sph->sphdata[i].m);
      for (k=0; k<ndim; k++) sph->sphdata[i].v[k] = 
	sqrt(2.0*kefrac*sph->sphdata[i].u)*
	sph->sphdata[i].r[k]/(drmag + small_number);
      sph->sphdata[i].u = ufrac*sph->sphdata[i].u;
      //,1.0e-6*umax/sph->sphdata[i].m);
    }
    else {
      sph->sphdata[i].u = 1.0e-6/sph->sphdata[i].m;
      for (k=0; k<ndim; k++) sph->sphdata[i].v[k] = (FLOAT) 0.0;
    }
  }

  delete[] hotlist;
  delete[] r;

  return;
}



//=============================================================================
//  Simulation::ShearFlow
/// Create shear-flow to test effective shear viscosity.
//=============================================================================
template <int ndim>
void Simulation<ndim>::ShearFlow(void)
{
  int i;                            // Particle counter
  int k;                            // Dimension counter
  int Nbox;                         // No. of particles in box
  int Nlattice1[ndim];              // Lattice size
  FLOAT lambda;                     // Wavelength if velocity perturbation
  FLOAT kwave;                      // Wavenumber
  FLOAT volume;                     // Volume of box
  FLOAT *r;                         // Positions of particles

  // Make local copies of important parameters
  FLOAT rhofluid1 = simparams->floatparams["rhofluid1"];
  FLOAT press1 = simparams->floatparams["press1"];
  FLOAT gammaone = simparams->floatparams["gamma_eos"] - 1.0;
  FLOAT amp = simparams->floatparams["amp"];
  Nlattice1[0] = simparams->intparams["Nlattice1[0]"];
  Nlattice1[1] = simparams->intparams["Nlattice1[1]"];

  debug2("[Simulation::ShockTube]");

  if (ndim != 2) {
    string message = "Shear-flow test only in 2D";
    ExceptionHandler::getIstance().raise(message);
  }

  // Compute size and range of fluid bounding boxes
  //---------------------------------------------------------------------------
  volume = (simbox.boxmax[0] - simbox.boxmin[0])*
    (simbox.boxmax[1] - simbox.boxmin[1]);
  Nbox = Nlattice1[0]*Nlattice1[1];

  lambda = simbox.boxmax[1] - simbox.boxmin[1];
  kwave = twopi/lambda;

  // Allocate local and main particle memory
  sph->Nsph = Nbox;
  AllocateParticleMemory();
  r = new FLOAT[ndim*sph->Nsph];


  // Add particles for LHS of the shocktube
  //---------------------------------------------------------------------------
  if (Nbox > 0) {
    AddCubicLattice(Nbox,Nlattice1,r,simbox,false);
    //AddHexagonalLattice(Nbox,Nlattice1,r,simbox,false);

    for (i=0; i<Nbox; i++) {
      for (k=0; k<ndim; k++) sph->sphdata[i].r[k] = r[ndim*i + k];
      for (k=0; k<ndim; k++) sph->sphdata[i].v[k] = 0.0;
      sph->sphdata[i].v[0] = amp*sin(kwave*sph->sphdata[i].r[1]);
      sph->sphdata[i].m = rhofluid1*volume/(FLOAT) Nbox;
      sph->sphdata[i].h = sph->h_fac*pow(sph->sphdata[i].m/rhofluid1,invndim);
      sph->sphdata[i].u = press1/rhofluid1/gammaone;
    }
  }

  delete[] r;

  return;
}



//=============================================================================
//  Simulation::SoundWave
/// Set-up isothermal sound-wave test.
//=============================================================================
template <int ndim>
void Simulation<ndim>::SoundWave(void)
{
  int i,k;                          // Particle and dimension counters
  int Nlattice1[ndim];              // Lattice size
  FLOAT csound;                     // (Isothermal) sound speed
  FLOAT diff;                       // ??
  FLOAT lambda;                     // Wavelength of perturbation
  FLOAT kwave;                      // Wave number of perturbing sound wave
  FLOAT omegawave;                  // Angular frequency of sound wave
  FLOAT ugas;                       // Internal energy of gas
  FLOAT xold;                       // Old x-position (for iteration)
  FLOAT xnew;                       // New x-position (for iteration)
  FLOAT *r;                         // Particle positions

  // Make local copies of parameters for setting up problem
  int Npart = simparams->intparams["Nsph"];
  FLOAT rhofluid1 = simparams->floatparams["rhofluid1"];
  FLOAT press1 = simparams->floatparams["press1"];
  FLOAT gamma = simparams->floatparams["gamma_eos"];
  FLOAT gammaone = gamma - 1.0;
  FLOAT amp = simparams->floatparams["amp"];
  FLOAT temp0 = simparams->floatparams["temp0"];
  FLOAT mu_bar = simparams->floatparams["mu_bar"];
  Nlattice1[0] = simparams->intparams["Nlattice1[0]"];

  debug2("[Simulation::SoundWave]");

  if (ndim != 1) {
    string message = "Sound wave only available in 1D";
    ExceptionHandler::getIstance().raise(message);
  }

  if (sph->gas_eos == "isothermal") {
    ugas = temp0/gammaone/mu_bar;
    press1 = gammaone*rhofluid1*ugas;
    csound = sqrt(press1/rhofluid1);
  }
  else {
    ugas = press1/rhofluid1/gammaone;
    csound = sqrt(gamma*press1/rhofluid1);
  }

  lambda = simbox.boxmax[0] - simbox.boxmin[0];
  kwave = twopi/lambda;
  omegawave = twopi*csound/lambda;

  // Allocate local and main particle memory
  sph->Nsph = Npart;
  Nlattice1[0] = Npart;
  AllocateParticleMemory();
  r = new FLOAT[ndim*sph->Nsph];

  // Add regular distribution of SPH particles
  AddCubicLattice(Npart,Nlattice1,r,simbox,false);


  // Set positions of all particles to produce density perturbation
  //---------------------------------------------------------------------------
  for (i=0; i<Npart; i++) {
    xnew = r[ndim*i];

    // Solve iterative procedure for particle positions in sound wave
    do {
      xold = xnew;
      xnew = r[ndim*i] - amp*(1.0 - cos(kwave*xnew))/kwave;
      diff = fabs((xnew - xold)/lambda);
    } while (diff > 1.0e-6);

    if (xnew > simbox.boxmax[0]) xnew -= simbox.boxsize[0];
    if (xnew < simbox.boxmin[0]) xnew += simbox.boxsize[0];

    // Set positions in main array with corresponind velocity perturbation
    for (k=0; k<ndim; k++) sph->sphdata[i].r[k] = xnew;
    for (k=0; k<ndim; k++) 
      sph->sphdata[i].v[k] = csound*amp*sin(kwave*xnew);
    sph->sphdata[i].m = rhofluid1*lambda/(FLOAT) Npart;
    sph->sphdata[i].h = sph->h_fac*pow(sph->sphdata[i].m/rhofluid1,invndim);

    if (sph->gas_eos == "isothermal")
      sph->sphdata[i].u = temp0/gammaone/mu_bar;
    else
      sph->sphdata[i].u = press1/rhofluid1/gammaone;
  }
  //---------------------------------------------------------------------------

  initial_h_provided = true;

  return;
}



//=============================================================================
//  Simulation::BinaryStar
/// Create a simple binary star problem
//=============================================================================
template <int ndim>
void Simulation<ndim>::BinaryStar(void)
{
  int k;                           // Dimension counter
  DOUBLE rbinary[ndim];            // Position of binary COM
  DOUBLE vbinary[ndim];            // Velocity of binary COM

  // Binary star parameters
  FLOAT abin = simparams->floatparams["abin"];
  FLOAT ebin = simparams->floatparams["ebin"];
  FLOAT m1 = simparams->floatparams["m1"];
  FLOAT m2 = simparams->floatparams["m2"];
  FLOAT phirot = simparams->floatparams["phirot"];
  FLOAT thetarot = simparams->floatparams["thetarot"];
  FLOAT psirot = simparams->floatparams["psirot"];

  debug2("[Simulation::BinaryStar]");

  if (ndim == 1) {
    string message = "Binary test not available in 1D";
    ExceptionHandler::getIstance().raise(message);
  }

  // Allocate local and main particle memory
  sph->Nsph = 0;
  sph->Ntot = 0;
  nbody->Nstar = 2;
  AllocateParticleMemory();

  // Add binary star
  for (k=0; k<ndim; k++) rbinary[k] = 0.0;
  for (k=0; k<ndim; k++) vbinary[k] = 0.0;
  AddBinaryStar(abin,ebin,m1,m2,0.01,0.01,phirot,thetarot,psirot,0.0,
                rbinary,vbinary,nbody->stardata[0],nbody->stardata[1]);

  return;
}



//=============================================================================
//  SphSimulation::TripleStar
/// Create a simple quadruple star problem
//=============================================================================
template <int ndim>
void Simulation<ndim>::TripleStar(void)
{
  int k;                           // Dimension counter
  DOUBLE rbinary[ndim];            // Position of binary COM
  DOUBLE vbinary[ndim];            // Velocity of binary COM
  NbodyParticle<ndim> b1;          // Star/binary 1
  NbodyParticle<ndim> b2;          // Star/binary 2

  // Triple star parameters
  FLOAT m1 = simparams->floatparams["m1"];
  FLOAT m2 = simparams->floatparams["m2"];
  FLOAT m3 = simparams->floatparams["m3"];
  FLOAT abin1 = simparams->floatparams["abin"];
  FLOAT abin2 = simparams->floatparams["abin2"];
  FLOAT ebin1 = simparams->floatparams["ebin"];
  FLOAT ebin2 = simparams->floatparams["ebin2"];
  FLOAT phirot = simparams->floatparams["phirot"];
  FLOAT thetarot = simparams->floatparams["thetarot"];
  FLOAT psirot = simparams->floatparams["psirot"];

  debug2("[SphSimulation::TripleStar]");

  if (ndim == 1) {
    string message = "Quadruple test not available in 1D";
    ExceptionHandler::getIstance().raise(message);
  }

  // Allocate local and main particle memory
  sph->Nsph = 0;
  sph->Ntot = 0;
  nbody->Nstar = 3;
  AllocateParticleMemory();

  // Compute main binary orbit
  for (k=0; k<ndim; k++) rbinary[k] = 0.0;
  for (k=0; k<ndim; k++) vbinary[k] = 0.0;
  AddBinaryStar(abin1,ebin1,m1,m2+m3,0.0001,0.0001,phirot,thetarot,psirot,
                0.0,rbinary,vbinary,b1,nbody->stardata[2]);

  // Now compute both components
  AddBinaryStar(abin2,ebin2,m2,m3,0.0001,0.0001,phirot,thetarot,psirot,0.0,
                b1.r,b1.v,nbody->stardata[0],nbody->stardata[1]);

  return;
}



//=============================================================================
//  SphSimulation::QuadrupleStar
/// Create a simple quadruple star problem
//=============================================================================
template <int ndim>
void Simulation<ndim>::QuadrupleStar(void)
{
  int k;                           // Dimension counter
  DOUBLE rbinary[ndim];            // Position of binary COM
  DOUBLE vbinary[ndim];            // Velocity of binary COM
  NbodyParticle<ndim> b1;          // Star/binary 1
  NbodyParticle<ndim> b2;          // Star/binary 2

  // Quadruple star parameters
  FLOAT m1 = simparams->floatparams["m1"];
  FLOAT m2 = simparams->floatparams["m2"];
  FLOAT m3 = simparams->floatparams["m3"];
  FLOAT m4 = simparams->floatparams["m3"];
  FLOAT abin1 = simparams->floatparams["abin"];
  FLOAT abin2 = simparams->floatparams["abin2"];
  FLOAT ebin1 = simparams->floatparams["ebin"];
  FLOAT ebin2 = simparams->floatparams["ebin2"];
  FLOAT phirot = simparams->floatparams["phirot"];
  FLOAT thetarot = simparams->floatparams["thetarot"];
  FLOAT psirot = simparams->floatparams["psirot"];

  debug2("[SphSimulation::QuadrupleStar]");

  if (ndim == 1) {
    string message = "Quadruple test not available in 1D";
    ExceptionHandler::getIstance().raise(message);
  }

  // Allocate local and main particle memory
  sph->Nsph = 0;
  sph->Ntot = 0;
  nbody->Nstar = 4;
  AllocateParticleMemory();

  // Compute main binary orbit
  for (k=0; k<ndim; k++) rbinary[k] = 0.0;
  for (k=0; k<ndim; k++) vbinary[k] = 0.0;
  AddBinaryStar(abin1,ebin1,m1+m2,m3+m4,0.01,0.01,phirot,thetarot,psirot,0.0,
                rbinary,vbinary,b1,b2);

  // Now compute components of both inner binaries
  AddBinaryStar(abin2,ebin2,m1,m2,0.0001,0.0001,phirot,thetarot,psirot,0.0,
                b1.r,b1.v,nbody->stardata[0],nbody->stardata[1]);
  AddBinaryStar(abin2,ebin2,m3,m4,0.0001,0.0001,phirot,thetarot,psirot,0.0,
                b2.r,b2.v,nbody->stardata[2],nbody->stardata[3]);

  return;
}



//=============================================================================
//  Simulation::AddBinaryStar
/// Add a binary star of given mass, eccentricity and separation.
/// (Code provided courtesy of S. P. Goodwin; 29/09/2013)
//=============================================================================
template <int ndim>
void Simulation<ndim>::AddBinaryStar
(DOUBLE sma,                       ///< Semi-major axis
 DOUBLE eccent,                    ///< Orbital eccentricity
 DOUBLE m1,                        ///< Mass of star 1
 DOUBLE m2,                        ///< Mass of star 2
 DOUBLE h1,                        ///< Smoothing length of star 1
 DOUBLE h2,                        ///< Smoothing length of star 2
 DOUBLE phirot,                    ///< 'phi' Euler rotation angle
 DOUBLE thetarot,                  ///< 'theta' Euler rotation angle
 DOUBLE phase,                     ///< Phase angle
 DOUBLE psirot,                    ///< 'tpsi' rotation angle
 DOUBLE *rbinary,                  ///< Position of COM of binary
 DOUBLE *vbinary,                  ///< Velocity of COM of binary
 NbodyParticle<ndim> &s1,          ///< Star 1
 NbodyParticle<ndim> &s2)          ///< Star 2
{
  int k;                           // Dimension counter
  FLOAT mbin = m1 + m2;            // Total binary mass

  debug2("[Simulation::AddBinaryStar]");

  if (ndim == 1) {
    string message = "Binary test not available in 1D";
    ExceptionHandler::getIstance().raise(message);
  }

  // randomly sample M to get theta
  FLOAT M = 2.0*pi*(FLOAT)(rand()%RAND_MAX)/(FLOAT)RAND_MAX;

  // from this solve to get eccentric anomoly E - e sin(theta) = M
  // N-R method x_1 = x_0 - f(x_0)/f'(x_0)
  FLOAT Ee = M;
  FLOAT old = M;
  do {
    old = Ee;
    Ee = Ee - (Ee - eccent*sin(Ee) - M)/(1.0 - eccent*cos(Ee));
  } while (fabs(old - Ee) > 1.0e-6);

  // Next get theta tan(t/2) = sqrt((1+e)/(1-e))*tan(E/2)
  FLOAT theta = sqrt((1.0 + eccent)/(1.0 - eccent))*tan(0.5*Ee);
  theta = 2.0*atan(theta);

  // Total separation
  FLOAT sep = sma*(1.0 - eccent*eccent)/(1.0 + eccent*cos(theta));

  // Get velocity
  FLOAT vel = (m1 + m2)*(2.0/sep - 1.0/sma);
  vel = sqrt(vel);

  // ..
  FLOAT hc = sqrt((1.0 + eccent*cos(theta))/(2.0 - sep/sma));
  FLOAT phi = acos(hc);

  // Set properties of star 1
  // put on x-y plane to start
  for (k=0; k<ndim; k++) s1.r[k] = 0.0;
  for (k=0; k<ndim; k++) s1.v[k] = 0.0;
  s1.m = m1;
  s1.h = h1;
  s1.invh = 1.0 / s1.h;
  s1.r[0] += sep*cos(theta)*m2/mbin;
  s1.r[1] += sep*sin(theta)*m2/mbin;
  s1.v[0] += -vel*cos(0.5*pi - theta + phi)*m2/mbin;
  s1.v[1] += vel*sin(0.5*pi - theta + phi)*m2/mbin;

  // Set properties of star 2
  for (k=0; k<ndim; k++) s2.r[k] = 0.0;
  for (k=0; k<ndim; k++) s2.v[k] = 0.0;
  s2.m = m2;
  s2.h = h2;
  s2.invh = 1.0 / s2.h;
  s2.r[0] -= sep*cos(theta)*m2/mbin;
  s2.r[1] -= sep*sin(theta)*m2/mbin;
  s2.v[0] -= -vel*cos(0.5*pi - theta + phi)*m1/mbin;
  s2.v[1] -= vel*sin(0.5*pi - theta + phi)*m1/mbin;

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



  return;
}




//=============================================================================
//  Simulation::AddRandomBox
/// Populate given bounding box with random particles.
//=============================================================================
template <int ndim>
void Simulation<ndim>::AddRandomBox
(int Npart,                         ///< [in] No. of particles
 FLOAT *r,                          ///< [out] Positions of particles
 DomainBox<ndim> box)               ///< [in] Bounding box containing particles
{
  debug2("[Simulation::AddRandomBox]");

  for (int i=0; i<Npart; i++) {
    for (int k=0; k<ndim; k++) {
      r[ndim*i + k] = box.boxmin[k] + (box.boxmax[k] - box.boxmin[k])*
	(FLOAT)(rand()%RAND_MAX)/(FLOAT)RAND_MAX;
    }
  }

  return;
}



//=============================================================================
//  Simulation::AddRandomsphere
/// Add random sphere of particles
//=============================================================================
template <int ndim>
void Simulation<ndim>::AddRandomSphere
(int Npart,                         ///< [in] No. of particles in sphere
 FLOAT *r,                          ///< [out] Positions of particles in sphere
 FLOAT *rcentre,                    ///< [in] Position of sphere centre
 FLOAT radius)                      ///< [in] Radius of sphere
{
  int i,k;                          // Particle and dimension counters
  FLOAT rad;                        // Radius of particle
  FLOAT rpos[ndim];                 // Random position of new particle

  debug2("[Simulation::AddRandomSphere]");

  // Loop over all required particles
  //---------------------------------------------------------------------------
  for (i=0; i<Npart; i++) {

    // Continously loop until random particle lies inside sphere
    do {
      for (k=0; k<ndim; k++) 
	rpos[k] = 1.0 - 2.0*(FLOAT)(rand()%RAND_MAX)/(FLOAT)RAND_MAX;
      rad = DotProduct(rpos,rpos,ndim);
    } while (rad > radius);

    for (k=0; k<ndim; k++) r[ndim*i + k] = rcentre[k] + rpos[k];
  }
  //---------------------------------------------------------------------------

  return;
}



//=============================================================================
//  Simulation::AddLatticeSphere
/// Add sphere of particles cut-out of regular lattice
//=============================================================================
template <int ndim>
int Simulation<ndim>::AddLatticeSphere
(int Npart,                         ///< [in] No. of particles in sphere
 FLOAT *r,                          ///< [out] Positions of particles in sphere
 FLOAT *rcentre,                    ///< [in] Position of sphere centre
 FLOAT radius,                      ///< [in] Radius of sphere
 string particle_dist)              ///< [in] String of lattice type
{
  int i,k;                          // Particle and dimension counters
  int Naux;                         // Aux. particle number
  int Nlattice[3];                  // Lattice size
  FLOAT *raux;                      // Temp. array to hold particle positions
  DomainBox<ndim> box1;             // Bounding box

  debug2("[Simulation::AddLatticeSphere]");

  // Set parameters for box and lattice to ensure it contains enough particles
  for (k=0; k<3; k++) Nlattice[k] = 1;
  for (k=0; k<ndim; k++) Nlattice[k] = 3*(int) powf((FLOAT) Npart, invndim);
  for (k=0; k<ndim; k++) box1.boxmin[k] = -2.0;
  for (k=0; k<ndim; k++) box1.boxmax[k] = 2.0;
  Naux = Nlattice[0]*Nlattice[1]*Nlattice[2];
  raux = new FLOAT[ndim*Naux];

  // Create a bounding box to hold lattice sphere
  if (particle_dist == "cubic_lattice")
    AddCubicLattice(Naux,Nlattice,raux,box1,true);
  else if (particle_dist == "hexagonal_lattice")
    AddHexagonalLattice(Naux,Nlattice,raux,box1,true);
  else {
    string message = "Invalid particle distribution option";
    ExceptionHandler::getIstance().raise(message);
  }

  // Now cut-out sphere from lattice containing exact number of particles 
  // (unless lattice structure prevents this).
  Naux = CutSphere(Npart,Naux,radius,raux,box1,false);

  // Copy particle positions to main position array to be returned
  for (i=0; i<Naux; i++)
    for (k=0; k<ndim; k++) r[ndim*i + k] = radius*raux[ndim*i + k];

  // Free allocated memory
  delete[] raux;

  return Naux;
}




//=============================================================================
//  Simulation::AddCubicLattice
/// Add regular (cubic) lattice of particles
//=============================================================================
template <int ndim>
void Simulation<ndim>::AddCubicLattice
(int Npart,                         ///< [in] No. of particles in lattice
 int Nlattice[ndim],                ///< [in] Ptcls per dimension in lattice
 FLOAT *r,                          ///< [out] Positions of particles
 DomainBox<ndim> box,               ///< [in] Bounding box of particles
 bool normalise)                    ///< [in] Normalise lattice shape and size
{
  int i,k;                          // Particle and dimension counters
  int ii,jj,kk;                     // Aux. lattice counters
  FLOAT spacing[ndim];              // Lattice spacing in each direction

  debug2("[Simulation::AddCubicLattice]");

  // If normalised, ensure equal spacing between all lattice layers.  
  // Otherwise set spacing to fit bounding box
  if (normalise) {
    for (k=0; k<ndim; k++) 
      spacing[k] = (box.boxmax[0] - box.boxmin[0])/(FLOAT) Nlattice[0];
  }
  else {
    for (k=0; k<ndim; k++) 
      spacing[k] = (box.boxmax[k] - box.boxmin[k])/(FLOAT) Nlattice[k];
  }

  
  // Create lattice depending on dimensionality
  //---------------------------------------------------------------------------
  if (ndim == 1) {
    for (ii=0; ii<Nlattice[0]; ii++) {
      i = ii;
      r[i] = box.boxmin[0] + ((FLOAT)ii + 0.5)*spacing[0];
    }
  }
  //---------------------------------------------------------------------------
  else if (ndim == 2) {
    for (jj=0; jj<Nlattice[1]; jj++) {
      for (ii=0; ii<Nlattice[0]; ii++) {
        i = jj*Nlattice[0] + ii;
        r[ndim*i] = box.boxmin[0] + ((FLOAT)ii + 0.5)*spacing[0];
        r[ndim*i + 1] = box.boxmin[1] + ((FLOAT)jj + 0.5)*spacing[1];
      }
    }
  }
  //---------------------------------------------------------------------------
  else if (ndim == 3) {
#pragma omp parallel for default(none)\
  shared(box,Nlattice,r,spacing) private(i,ii,jj,kk)
    for (kk=0; kk<Nlattice[2]; kk++) {
      for (jj=0; jj<Nlattice[1]; jj++) {
	for (ii=0; ii<Nlattice[0]; ii++) {
	  i = kk*Nlattice[0]*Nlattice[1] + jj*Nlattice[0] + ii;
	  r[ndim*i] = box.boxmin[0] + ((FLOAT)ii + 0.5)*spacing[0];
	  r[ndim*i + 1] = box.boxmin[1] + ((FLOAT)jj + 0.5)*spacing[1];
	  r[ndim*i + 2] = box.boxmin[2] + ((FLOAT)kk + 0.5)*spacing[2];
	}
      }
    }
  }

  return;
}



//=============================================================================
//  Simulation::AddHexagonalLattice
/// Create simple hexagonal-packed lattice using A-B-A-B pattern in z-direction
/// N.B. the box is scaled to fit to the x-boxsize.
//=============================================================================
template <int ndim>
void Simulation<ndim>::AddHexagonalLattice
(int Npart,                         ///< [in] No. of particles in lattice
 int Nlattice[ndim],                ///< [in] Ptcls per dimension in lattice
 FLOAT *r,                          ///< [out] Positions of particles
 DomainBox<ndim> box,               ///< [in] Bounding box of particles
 bool normalise)                    ///< [in] Normalise lattice shape and size
{
  int i,k;                          // Particle and dimension counters
  int ii,jj,kk;                     // Aux. lattice counters
  FLOAT rad[ndim];                  // 'Radius' of particle in lattice

  debug2("[Simulation::AddHexagonalLattice]");

  // If normalised, ensure equal spacing between all particles.  
  // Otherwise set spacing to fit bounding box.
  if (normalise) {
    for (k=0; k<ndim; k++) 
      rad[k] = 0.5*(box.boxmax[0] - box.boxmin[0])/(FLOAT) Nlattice[0];
  }
  else {
    for (k=0; k<ndim; k++) 
      rad[k] = 0.5*(box.boxmax[k] - box.boxmin[k])/(FLOAT) Nlattice[k];
  }
  

  // Create lattice depending on dimensionality
  //---------------------------------------------------------------------------
  if (ndim == 1) {
    for (ii=0; ii<Nlattice[0]; ii++) {
      i = ii;
      r[i] = box.boxmin[0] + 0.5*rad[0] + 2.0*(FLOAT)ii*rad[0];
    }
  }

  //---------------------------------------------------------------------------
  else if (ndim == 2) {
    for (jj=0; jj<Nlattice[1]; jj++) {
      for (ii=0; ii<Nlattice[0]; ii++) {
	i = jj*Nlattice[0] + ii;
	r[ndim*i] = box.boxmin[0] + 0.5*rad[0] + 
	  (2.0*(FLOAT)ii + (FLOAT)(jj%2))*rad[0];
	r[ndim*i + 1] = box.boxmin[1] + 0.5*sqrt(3.0)*rad[1] + 
	  (FLOAT)jj*sqrt(3.0)*rad[1];
      }
    }
  }

  //---------------------------------------------------------------------------
  else if (ndim == 3) {
#pragma omp parallel for default(none)\
  shared(box,Nlattice,r,rad) private(i,ii,jj,kk)
    for (kk=0; kk<Nlattice[2]; kk++) {
      for (jj=0; jj<Nlattice[1]; jj++) {
	for (ii=0; ii<Nlattice[0]; ii++) {
	  i = kk*Nlattice[0]*Nlattice[1] + jj*Nlattice[0] + ii;
	  r[ndim*i] = box.boxmin[0] + 0.5*rad[0] + 
	    (2.0*(FLOAT)ii + (FLOAT)(jj%2) + (FLOAT)((kk+1)%2))*rad[0];
	  r[ndim*i + 1] = box.boxmin[1] + 0.5*sqrt(3.0)*rad[1] + 
	    (FLOAT)jj*sqrt(3.0)*rad[1] + (FLOAT)(kk%2)*rad[1]/sqrt(3.0);
	  r[ndim*i + 2] = box.boxmin[2] + sqrt(6.0)*rad[2]/3.0 + 
	    (FLOAT)kk*2.0*sqrt(6.0)*rad[2]/3.0;
	}
      }
    }
  }

  return;
}



//=============================================================================
//  Simulation::CutSphere
/// Cut-out a sphere containing exactly 'Nsphere' particles from a uniform 
/// box of particles.
//=============================================================================
template <int ndim>
int Simulation<ndim>::CutSphere
(int Nsphere,                       ///< [in] Desired np of particles in sphere
 int Npart,                         ///< [in] No. of particles in cube
 FLOAT radsphere,                   ///< [in] ??
 FLOAT *r,                          ///< [inout] Positions of particles
 DomainBox<ndim> box,               ///< [in] Bounding box of particles
 bool exact)                        ///< [in] ??
{
  int i,k;                          // Particle and dimension counters
  int Ninterior = 0;                // No. of particle
  FLOAT dr[ndim];                   // Relative position vector
  FLOAT drsqd;                      // Distance squared
  FLOAT r_low = 0.0;                // Lower-bound for bisection iteration
  FLOAT r_high;                     // Upper-bound for bisection iteration
  FLOAT radius;                     // Current radius containing Nsphere ptcls
  FLOAT rcentre[ndim];              // Centre of sphere

  debug2("[Simulation::CutSphere]");

  // Find centre and shortest edge-length of bounding box
  r_high = (FLOAT) big_number;
  for (k=0; k<ndim; k++) {
    rcentre[k] = 0.5*(box.boxmin[k] + box.boxmax[k]);
    r_high = min(r_high,(FLOAT)0.5*(box.boxmax[k] - box.boxmin[k]));
  }

  // Bisection iteration to determine the radius containing the desired 
  // number of particles
  //---------------------------------------------------------------------------
  do {
    radius = 0.5*(r_low + r_high);
    Ninterior = 0;

    // Count how many particles lie inside current radius
    for (i=0; i<Npart; i++) {
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
  //---------------------------------------------------------------------------
  Ninterior = 0;
  for (i=0; i<Npart; i++) {
    for (k=0; k<ndim; k++) dr[k] = r[ndim*i + k] - rcentre[k];
    drsqd = DotProduct(dr,dr,ndim);
    if (drsqd <= radius*radius) {
      for (k=0; k<ndim; k++) r[ndim*Ninterior + k] = r[ndim*i + k]/radius;
      Ninterior++;
    }
  }

  return Ninterior;
}



//=============================================================================
//  Simulation::AddAzimuthalDensityPerturbation
/// Add an azimuthal density perturbation for implementing Boss-Bodenheimer
/// type initial conditions
//=============================================================================
template <int ndim>
void Simulation<ndim>::AddAzimuthalDensityPerturbation
(int Npart,                         ///< [in] No. of particles in sphere
 int mpert,                         ///< [in] Perturbation mode
 FLOAT amp,                         ///< [in] Amplitude of perturbation
 FLOAT *rcentre,                    ///< [in] Position of sphere centre
 FLOAT *r)                          ///< [inout] Positions of particles
{
  int i,k;                          // Particle and dimension counters
  int j;
  int tabtot;                       // ..
  FLOAT phi,phi1,phi2,phiprime;     // ..
  FLOAT Rsqd;                       // ..
  FLOAT Rmag;                       // ..
  FLOAT rpos[ndim];                 // Random position of new particle
  FLOAT spacing;                    // ..

  debug2("[Simulation::AddAzimuthalDensityPerturbation]");

  tabtot = 10000;
  spacing = twopi/(FLOAT)(tabtot - 1);

  // Loop over all required particles
  //---------------------------------------------------------------------------
#pragma omp parallel for default(none)\
  shared(amp,mpert,Npart,r,rcentre,spacing,tabtot)\
  private(i,j,k,phi,phiprime,phi1,phi2,rpos,Rmag,Rsqd)
  for (i=0; i<Npart; i++) {
    for (k=0; k<ndim; k++) rpos[k] = r[ndim*i + k] - rcentre[k];

    // Calculate distance from z-axis and 
    Rsqd = rpos[0]*rpos[0] + rpos[1]*rpos[1];
    Rmag = sqrt(Rsqd);

    // Find azimuthal angle around z-axis correcting for which quadrant
    if (Rmag > small_number) phi = asin(fabs(rpos[1])/Rmag);
    else phi = 0.0;

    if (rpos[0] < 0.0 && rpos[1] > 0.0) phi = pi - phi;
    else if (rpos[0] < 0.0 && rpos[1] < 0.0) phi = pi + phi;
    else if (rpos[0] > 0.0 && rpos[1] < 0.0) phi = 2.0*pi - phi;

    // Wrap angle to fit between 0 and two*pi
    if (phi < amp/(FLOAT) mpert) phi = phi + twopi;

    // Numerically find new phi angle for perturbation.  Search through
    // grid of values, find upper and lower bounds, then use linear 
    // interpolation to find new value of phi.
    for (j=1; j<tabtot; j++) {
      phi1 = spacing*(FLOAT) (j - 1);
      phi2 = spacing*(FLOAT) j;
      phi1 = phi1 + amp*cos((FLOAT) mpert*phi1)/(FLOAT) mpert;
      phi2 = phi2 + amp*cos((FLOAT) mpert*phi2)/(FLOAT) mpert;

      if (phi2 >= phi && phi1 < phi) {
        phiprime = spacing*(FLOAT)(j - 1) +
          spacing*(phi - phi1) / (phi2 - phi1);
        break;
      }
    }

    // Reposition particle with new angle
    r[ndim*i] = rcentre[0] + Rmag*cos(phiprime);
    r[ndim*i + 1] = rcentre[1] + Rmag*sin(phiprime);

  }
  //---------------------------------------------------------------------------

  return;
}



//=============================================================================
//  Simulation::AddRotationalVelocityField
/// Add a solid-body rotational velocity field
//=============================================================================
template <int ndim>
void Simulation<ndim>::AddRotationalVelocityField
(int Npart,                         ///< [in] No. of particles in sphere
 FLOAT angvelaux,                   ///< [in] Angular velocity of cloud
 FLOAT *rcentre,                    ///< [in] Position of sphere centre
 FLOAT *r,                          ///< [in] Positions of particles
 FLOAT *v)                          ///< [out] Velocities of particles
{
  int i,k;                          // Particle and dimension counters
  FLOAT Rmag;                       // ..
  FLOAT Rsqd;                       // ..
  FLOAT dr[ndim];                   // ..

  debug2("[Simulation::AddAzimuthalDensityPerturbation]");


  // Loop over all required particles
  //---------------------------------------------------------------------------
#pragma omp parallel for default(none)\
  shared(angvelaux,Npart,r,rcentre,v) private(dr,i,k,Rmag,Rsqd)
  for (i=0; i<Npart; i++) {
    for (k=0; k<ndim; k++) dr[k] = r[ndim*i + k] - rcentre[k];
    for (k=0; k<ndim; k++) v[ndim*i + k] = 0.0;

    // Calculate distance from z-axis and 
    Rsqd = dr[0]*dr[0] + dr[1]*dr[1];
    Rmag = sqrt(Rsqd);

    // Find azimuthal angle around z-axis correcting for which quadrant
    if (Rmag > small_number) {
      dr[0] = dr[0]/Rmag;
      dr[1] = dr[1]/Rmag;

      v[ndim*i] = -angvelaux*Rmag*dr[1];
      v[ndim*i + 1] = angvelaux*Rmag*dr[0];
    }
  }
  //---------------------------------------------------------------------------

  return;
}



#if defined(FFTW_TURBULENCE)
//=============================================================================
//  Simulation::GenerateTurbulentVelocityField
/// ..
/// ..
//=============================================================================
template <int ndim>
void Simulation<ndim>::GenerateTurbulentVelocityField
(int field_type,                    ///< Type of turbulent velocity field
 int gridsize,                      ///< Size of velocity grid
 DOUBLE power_turb,                 ///< Power spectrum index
 DOUBLE *vfield)                    ///< Array containing velocity field grid
{
  bool divfree;                     // Select div-free turbulence
  bool curlfree;                    // Select curl-free turbulence
  int kmax;                         // Max. extent of k (in 3D)
  int kmin;                         // Min. extent of k (in 3D)
  int krange;                       // Range of k values (kmax - kmin + 1)
  int shift;                        // Power grid shift
  int i,j,k;                        // Grid counters
  int ii,jj,kk;                     // Aux. grid counters
  int k1,k2,k3;                     // ??
  int d;                            // Dimension counter
  DOUBLE F[3];                      // Fourier vector component
  DOUBLE unitk[3];                  // Unit k-vector
  DOUBLE *power, *phase;            // Fourier components
  DOUBLE *dummy1, *dummy2;          // Dummy arrays for array resequencing
  DOUBLE Rnd[3],w;                  // Random numbers, variable in Gaussian calculation
  DOUBLE k_rot[3];                  // bulk rotation modes
  DOUBLE k_com[3];                  // bulk compression modes
  fftw_plan plan;                   // ??
  fftw_complex *complexfield;       // ..


  debug2("[Simulation::GenerateTurbulentVelocityField]");

  // Initalise random number seed
  for (i=0; i<simparams->intparams["randseed"]; i++) j = rand()%RAND_MAX;

  divfree = false;
  curlfree = false;
  if (field_type == 1) curlfree = true;
  if (field_type == 2) divfree = true;

  kmin = -(gridsize/2 - 1);
  kmax = gridsize/2;
  krange = kmax - kmin + 1;
  cout << "kmin : " << kmin << "    kmax : " << kmax << "    krange : " 
       << krange << "    gridsize : " << gridsize << endl;
  if (krange != gridsize) exit(0);

  dummy1 = new DOUBLE[3*krange*krange*krange];
  dummy2 = new DOUBLE[3*krange*krange*krange];
  power = new DOUBLE[3*krange*krange*krange];
  phase = new DOUBLE[3*krange*krange*krange];
  complexfield = new fftw_complex[gridsize*gridsize*gridsize];
  //complexfield = new fftw_complex[krangep1*krangep1*krangep1];

  for (i=0; i<3*krange*krange*krange; i++) power[i] = 0.0;
  for (i=0; i<3*krange*krange*krange; i++) phase[i] = 0.0;
  for (i=0; i<3*gridsize*gridsize*gridsize; i++) vfield[i] = 0.0;


  // Define wave vectors in Fourier space
  // Each wave vector has coordinates in Fourier space, random phases in
  // three dimensions, and power in three dimensions, giving a power
  // and a phase vector field in Fourier space
  // With a 'natural' field type there is no coordination between x,y,z components
  // With a div-free or curl-free field, there is a relation between
  // coordinates in Fourier space and power vector. For a div-free field,
  // power vectors are perpendicular to the Fourier coordinate vector, and
  // for a curl-free field power vectors are (anti-)parallel to the Fourier
  // coordinate vector
  // (i,j,k) is the Fourier coordinate vector
  // power(1:3,i,j,k) is the power vector at Fourier coordinates i,j,k
  // phase(1:3,i,j,k) is the phase vector at Fourier coordinates i,j,k
  //----------------------------------------------------------------------------
  for (i=kmin; i<=kmax; i++) {
    for (j=kmin; j<=kmax; j++) {
      for (k=kmin; k<=kmax; k++) {
	ii = i + kmin;
	jj = j + kmin;
	kk = k + kmin;
	  
	// cycle antiparallel k-vectors
	//if (k < 0) continue;            
	//if (k == 0) {
	//  if (j < 0) continue;
	//  if (j == 0 && i < 0) continue;
	//}

	// Central power = 0
	if (i == 0 && j == 0 && k == 0) continue;
        if (i*i + j*j + k*k > kmax*kmax) continue;

	// Power value, to be multipled by random power chosen from a Gaussian
	// This is what gives the slope of the eventual power spectrum
	for (d=0; d<3; d++)
	  F[d] = sqrt(pow(sqrt((DOUBLE)(i*i + j*j + k*k)),power_turb));

	for (d=0; d<3; d++) {
	  Rnd[0] = (FLOAT)(rand()%RAND_MAX)/(FLOAT)RAND_MAX;

	  // Random phase between 0 and 2*pi
          phase[d + 3*ii + 3*krange*jj + 3*krange*krange*kk] =
	    (2.0*Rnd[0] - 1.0)*pi;	  

	  // Create Gaussian distributed random numbers
	  Rnd[1] = GaussRand(0.0,1.0);
	  Rnd[2] = GaussRand(0.0,1.0);
	  //cout << "Rnd : " << Rnd[0] << "    " << Rnd[1] << "   " << Rnd[2] << endl;
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
	  for (d=0; d<3; d++)
	    power[d + 3*ii + 3*krange*jj + 3*krange*krange*kk] 
	      = unitk[d]*DotProduct(F,unitk,3);
	}
	// For divergence free turbulence, vector F should be perpendicular to vector k
	else if (divfree) {
	  for (d=0; d<3; d++)
	    power[d + 3*ii + 3*krange*jj + 3*krange*krange*kk] 
	      = F[d] - unitk[d]*DotProduct(F,unitk,3);
	}
	else {
	  for (d=0; d<3; d++)
	    power[d + 3*ii + 3*krange*jj + 3*krange*krange*kk] 
	      = F[d];
	}
	cout << "POWER : " << d << "  " << sqrt(i*i + j*j + k*k) << "   " << F[d] << endl;
	/*cout << "POWER[" << i << "," << j << "," << k << "] : " 
	     << power[3*ii + 3*krange*jj + 3*krange*krange*kk] << "   "
	     << power[1 + 3*ii + 3*krange*jj + 3*krange*krange*kk] << "   "
	     << power[2 + 3*ii + 3*krange*jj + 3*krange*krange*kk] << endl;
	cout << "PHASE[" << i << "," << j << "," << k << "] : " 
	     << phase[3*ii + 3*krange*jj + 3*krange*krange*kk] << "   "
	     << phase[1 + 3*ii + 3*krange*jj + 3*krange*krange*kk] << "   "
	     << phase[2 + 3*ii + 3*krange*jj + 3*krange*krange*kk] << endl;*/
      }
    }
  }
  //----------------------------------------------------------------------------


  DOUBLE power_spectrum[gridsize+1];
  for (i=0; i<kmax+1; i++) power_spectrum[i] = 0.0;
  for (i=kmin; i<=kmax; i++) {
    for (j=kmin; j<=kmax; j++) {
      for (k=kmin; k<=kmax; k++) {
	ii = i + kmin;
	jj = j + kmin;
	kk = k + kmin;
	DOUBLE kmag = sqrt((DOUBLE)(i*i + j*j + k*k));
	int ibin = (int) kmag;
	if (kmag > kmax) kmag = kmax;
	power_spectrum[ibin] += 
	  sqrt(pow(power[3*ii + 3*krange*jj + 3*krange*krange*kk],2) + 
	       pow(power[1 + 3*ii + 3*krange*jj + 3*krange*krange*kk],2) +
	       pow(power[2 + 3*ii + 3*krange*jj + 3*krange*krange*kk],2));
      }
    }
  }
  for (i=0; i<kmax+1; i++) 
    cout << "POWER : " << i << "   " << power_spectrum[i] << endl;





  plan = fftw_plan_dft_3d(gridsize, gridsize, gridsize, complexfield,
			  complexfield, FFTW_BACKWARD, FFTW_ESTIMATE);

  for (i=kmin; i<=kmax; i++) {
    for (j=kmin; j<=kmax; j++) {
      for (k=kmin; k<=kmax; k++) {
	ii = i + kmin;
	jj = j + kmin;
	kk = k + kmin;
	for (d=0; d<3; d++) {
	  dummy1[d + 3*ii + 3*krange*jj + 3*krange*krange*kk] =
	    phase[d + 3*ii + 3*krange*jj + 3*krange*krange*kk];
	  dummy2[d + 3*ii + 3*krange*jj + 3*krange*krange*kk] =
	    power[d + 3*ii + 3*krange*jj + 3*krange*krange*kk];
	}
	/*cout << "DUMMY1[" << i << "," << j << "," << k << "] : " 
	     << dummy1[3*ii + 3*krange*jj + 3*krange*krange*kk] << "   "
	     << dummy1[1 + 3*ii + 3*krange*jj + 3*krange*krange*kk] << "   "
	     << dummy1[2 + 3*ii + 3*krange*jj + 3*krange*krange*kk] << endl;
	cout << "DUMMY2[" << i << "," << j << "," << k << "] : " 
	     << dummy2[3*ii + 3*krange*jj + 3*krange*krange*kk] << "   "
	     << dummy2[1 + 3*ii + 3*krange*jj + 3*krange*krange*kk] << "   "
	     << dummy2[2 + 3*ii + 3*krange*jj + 3*krange*krange*kk] << endl;*/

      }
    }
  }


  shift = -kmin;


  // reorder array: positive wavenumbers are placed in ascending order along
  // first half of dimension, i.e. 0 to k_max, negative wavenumbers are placed
  // along second half of dimension, i.e. -k_min to 1.
  for (d=0; d<3; d++) {

    for (i=kmin; i<=kmax; i++) {
      for (j=kmin; j<=kmax; j++) {
	for (k=kmin; k<=kmax; k++) {
	  ii = i + kmin;
	  jj = j + kmin;
	  kk = k + kmin;
	  phase[d + 3*ii + 3*krange*jj + 3*krange*krange*kk] = 
	    dummy1[d + 3*((ii + shift)%krange) + 
		   3*krange*((jj + shift)%krange) + 
		   3*krange*krange*((kk + shift)%krange)];
	  power[d + 3*ii + 3*krange*jj + 3*krange*krange*kk] = 
	    dummy2[d + 3*((ii + shift)%krange) + 
		   3*krange*((jj + shift)%krange) + 
		   3*krange*krange*((kk + shift)%krange)];

	  complexfield[0][ii +krange*jj +krange*krange*kk]
	    = power[d + 3*ii + 3*krange*jj + 3*krange*krange*kk]*
	    cos(phase[d + 3*ii + 3*krange*jj + 3*krange*krange*kk]);
	  complexfield[1][ii + krange*jj + krange*krange*kk]
	    = power[d + 3*ii + 3*krange*jj + 3*krange*krange*kk]*
	    sin(phase[d + 3*ii + 3*krange*jj + 3*krange*krange*kk]);
	  /*cout << "COMPLEX[" << i << "," << j << "," << k << "] : " 
	       << complexfield[0][ii +krange*jj +krange*krange*kk] << "   "
	       << complexfield[1][ii +krange*jj +krange*krange*kk] << endl;
	  cout << "POWER[" << i << "," << j << "," << k << "] : " 
	       << power[d + 3*ii + 3*krange*jj + 3*krange*krange*kk] << endl;
	  cout << "PHASE[" << i << "," << j << "," << k << "] : " 
	  << phase[d + 3*ii + 3*krange*jj + 3*krange*krange*kk] << endl;*/

	}
      }
    }
    
    fftw_execute_dft(plan, complexfield, complexfield);
    

    for (i=0; i<=krange; i++) {
      for (j=0; j<=krange; j++) {
	for (k=0; k<=krange; k++) {
	  vfield[d + 3*i + 3*krange*j + 3*krange*krange*k] = 
	    complexfield[0][i + krange*j + krange*krange*k];
	}
      }
    }

  }


  /*for (i=0; i<krange; i++) {
    for (j=0; j<krange; j++) {
      for (k=0; k<krange; k++) {
  	cout << "v[" << i << "," << j << "," << k << "] : " 
	     << vfield[3*i + 3*krange*j + 3*krange*krange*k] << "    " 
	     << vfield[1 + 3*i + 3*krange*j + 3*krange*krange*k] << "    " 
	     << vfield[2 + 3*i + 3*krange*j + 3*krange*krange*k] << "    " 
	     << endl;
      }
    }
    }*/
  

  fftw_destroy_plan(plan);


  delete[] complexfield;
  delete[] power;
  delete[] phase;
  delete[] dummy2;
  delete[] dummy1;

  return;
}
#endif

