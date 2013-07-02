//=============================================================================
//  SimulationIC.hpp
//  Contains all routines for generating initial conditions on the fly.
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
using namespace std;



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
  // --------------------------------------------------------------------------
  for (i=0; i<sph->Nsph; i++) {
    part = &sph->sphdata[i];
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
  // --------------------------------------------------------------------------

  if (!valid_ic) {
    string message = "Invalid initial conditions for SPH particles";
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
  FLOAT volume;                     // Volume of box
  FLOAT vfluid1[ndim];              // Velocity vector of LHS fluid
  FLOAT vfluid2[ndim];              // Velocity vector of RHS fluid
  FLOAT *r;                         // Position vectors
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
  // --------------------------------------------------------------------------
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
  cout << "Allocating memory : " << sph->Nsph << endl;


  // Add particles for LHS of the shocktube
  // --------------------------------------------------------------------------
  if (Nbox1 > 0) {
    AddCubicLattice(Nbox1,Nlattice1,r,box1,false);

    for (i=0; i<Nbox1; i++) {
      for (k=0; k<ndim; k++) sph->sphdata[i].r[k] = r[ndim*i + k];
      for (k=0; k<ndim; k++) sph->sphdata[i].v[k] = 0.0;
      sph->sphdata[i].v[0] = vfluid1[0];
      sph->sphdata[i].m = rhofluid1*volume/(FLOAT) Nbox1;
      if (sph->gas_eos == "isothermal")
    	sph->sphdata[i].u = temp0/gammaone/mu_bar;
      else
        sph->sphdata[i].u = press1/rhofluid1/gammaone;
    }
  }

  // Add particles for RHS of the shocktube
  // --------------------------------------------------------------------------
  if (Nbox2 > 0) {
    AddCubicLattice(Nbox2,Nlattice2,r,box2,false);

    for (j=0; j<Nbox2; j++) {
      i = Nbox1 + j;
      for (k=0; k<ndim; k++) sph->sphdata[i].r[k] = r[ndim*j + k];
      for (k=0; k<ndim; k++) sph->sphdata[i].v[k] = 0.0;
      sph->sphdata[i].v[0] = vfluid2[0];
      sph->sphdata[i].m = rhofluid2*volume/(FLOAT) Nbox2;
      if (sph->gas_eos == "isothermal")
	    sph->sphdata[i].u = temp0/gammaone/mu_bar;
      else
	    sph->sphdata[i].u = press2/rhofluid2/gammaone;
    }
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
  int Npart = simparams->intparams["Npart"];
  FLOAT rhobox = simparams->intparams["rhofluid1"];
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
    sph->sphdata[i].invomega = (FLOAT) 1.0;
    sph->sphdata[i].iorig = i;
    sph->sphdata[i].u = (FLOAT) 1.5;
  }

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
  int Npart = simparams->intparams["Npart"];
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
  for (i=0; i<sph->Nsph; i++) {
    for (k=0; k<ndim; k++) {
      sph->sphdata[i].r[k] = r[ndim*i + k];
      sph->sphdata[i].v[k] = (FLOAT) 0.0;
      sph->sphdata[i].a[k] = (FLOAT) 0.0;
    }
    sph->sphdata[i].m = rhofluid*volume / (FLOAT) Npart;
    sph->sphdata[i].u = press/rhofluid/gammaone;
    sph->sphdata[i].invomega = (FLOAT) 1.0;
    sph->sphdata[i].zeta = (FLOAT) 0.0;
    sph->sphdata[i].iorig = i;
  }

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
  int k;                            // Dimension counter
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
  FLOAT amp = simparams->floatparams["amp"];
  FLOAT lambda = simparams->floatparams["lambda"];
  Nlattice1[0] = simparams->intparams["Nlattice1[0]"];
  Nlattice1[1] = simparams->intparams["Nlattice1[1]"];
  Nlattice2[0] = simparams->intparams["Nlattice2[0]"];
  Nlattice2[1] = simparams->intparams["Nlattice2[1]"];
  vfluid1[0] = simparams->floatparams["vfluid1[0]"];
  vfluid2[0] = simparams->floatparams["vfluid2[0]"];

  debug2("[Simulation::ContactDiscontinuity]");


  // 1D simulation
  // ==========================================================================
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
    cout << "Allocating memory : " << sph->Nsph << endl;


    // ------------------------------------------------------------------------
    if (Nbox1 > 0) {
      AddCubicLattice(Nbox1,Nlattice1,r,box1,false);
      volume = box1.boxmax[0] - box1.boxmin[0];
      cout << "Vol1 : " << volume << "   m1 : " 
	   << rhofluid1*volume/(FLOAT) Nbox1 << "   rho : "
	   << rhofluid1*volume << "   Nbox : " << Nbox1 << endl;
      for (i=0; i<Nbox1; i++) {
        sph->sphdata[i].r[0] = r[i] - 0.4*simbox.boxsize[0];
        if (sph->sphdata[i].r[0] < simbox.boxmin[0])
          sph->sphdata[i].r[0] += simbox.boxsize[0];
        sph->sphdata[i].v[0] = 0.0;
        sph->sphdata[i].m = rhofluid1*volume/(FLOAT) Nbox1;
        if (sph->gas_eos == "isothermal")
          sph->sphdata[i].u = temp0/gammaone/mu_bar;
        else
          sph->sphdata[i].u = press1/rhofluid1/gammaone;
      }
    }

    // ------------------------------------------------------------------------
    if (Nbox2 > 0) {
      AddCubicLattice(Nbox2,Nlattice2,r,box2,false);
      volume = box2.boxmax[0] - box2.boxmin[0];
      cout << "Vol2 : " << volume << "   m2 : " 
	   << rhofluid2*volume/(FLOAT) Nbox2 
	   << "   rho : " << rhofluid2*volume << "   Nbox : " << Nbox2 << endl;
      for (j=0; j<Nbox2; j++) {
        i = Nbox1 + j;
        sph->sphdata[i].r[0] = r[j] - 0.4*simbox.boxsize[0];
        if (sph->sphdata[i].r[0] < simbox.boxmin[0])
          sph->sphdata[i].r[0] += simbox.boxsize[0];
        sph->sphdata[i].v[0] = 0.0;
        sph->sphdata[i].m = rhofluid2*volume/(FLOAT) Nbox2;
        if (sph->gas_eos == "isothermal")
          sph->sphdata[i].u = temp0/gammaone/mu_bar;
        else
          sph->sphdata[i].u = press2/rhofluid2/gammaone;
      }
    }

  }
  // ==========================================================================
  else if (ndim == 2) {



  }
  // ==========================================================================


  // Set initial smoothing lengths and create initial ghost particles
  // --------------------------------------------------------------------------
  sph->Nghost = 0;
  sph->Nghostmax = sph->Nsphmax - sph->Nsph;
  sph->Ntot = sph->Nsph;
  for (int i=0; i<sph->Nsph; i++) sph->sphdata[i].active = true;

  sph->InitialSmoothingLengthGuess();
  sphneib->UpdateTree(sph,*simparams);

  // Search ghost particles
  ghosts.SearchGhostParticles(simbox,sph);


  // Update neighbour tre
  sphneib->UpdateTree(sph,*simparams);

  // Calculate all SPH properties
  sphneib->UpdateAllSphProperties(sph,nbody);

  sphneib->UpdateTree(sph,*simparams);
  sphneib->UpdateAllSphProperties(sph,nbody);

  ghosts.CopySphDataToGhosts(sph);

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
  int i;
  int j;
  int k;
  int Nbox1;
  int Nbox2;
  FLOAT volume;
  FLOAT *r;
  DomainBox<ndim> box1;
  DomainBox<ndim> box2;
  int Nlattice1[ndim];
  int Nlattice2[ndim];
  FLOAT vfluid1[ndim];
  FLOAT vfluid2[ndim];

  // Record local copies of all important parameters
  FLOAT rhofluid1 = simparams->floatparams["rhofluid1"];
  FLOAT rhofluid2 = simparams->floatparams["rhofluid2"];
  FLOAT press1 = simparams->floatparams["press1"];
  FLOAT press2 = simparams->floatparams["press2"];
  FLOAT temp0 = simparams->floatparams["temp0"];
  FLOAT mu_bar = simparams->floatparams["mu_bar"];
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
  // --------------------------------------------------------------------------
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
  cout << "Nbox1 : " << Nbox1 << "    Nbox2 : " << Nbox2 << endl;
  cout << "Allocating memory : " << sph->Nsph << endl;


  // Add particles for LHS of the shocktube
  // --------------------------------------------------------------------------
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
      sph->sphdata[i].u = press1/rhofluid1/gammaone;
    }
  }

  // Add particles for RHS of the shocktube
  // --------------------------------------------------------------------------
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
      sph->sphdata[i].u = press2/rhofluid2/gammaone;
    }
  }

  // Add velocity perturbation here
  // --------------------------------------------------------------------------
  FLOAT sigmapert = 0.05/sqrt(2.0);
  for (i=0; i<sph->Nsph; i++) {
    sph->sphdata[i].v[1] = amp*sin(2.0*pi*sph->sphdata[i].r[0]/lambda)*
      (exp(-pow(sph->sphdata[i].r[1] + 0.25,2)/2.0/sigmapert/sigmapert) +  
       exp(-pow(sph->sphdata[i].r[1] - 0.25,2)/2.0/sigmapert/sigmapert));
  }

  // Set initial smoothing lengths and create initial ghost particles
  // --------------------------------------------------------------------------
  sph->Nghost = 0;
  sph->Nghostmax = sph->Nsphmax - sph->Nsph;
  sph->Ntot = sph->Nsph;
  for (int i=0; i<sph->Nsph; i++) sph->sphdata[i].active = true;
  
  sph->InitialSmoothingLengthGuess();
  sphneib->UpdateTree(sph,*simparams);
  
  // Search ghost particles
  ghosts.SearchGhostParticles(simbox,sph);

  sphneib->UpdateAllSphProperties(sph,nbody);

  // Update neighbour tre
  sphneib->UpdateTree(sph,*simparams);

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
  int j;                            // Aux. particle counter
  int k;                            // Dimension counter
  int Nsphere;                      // Actual number of particles in sphere
  int Nlattice[3];                  // Lattice size
  FLOAT dr[ndim];                   // Relative position vector
  FLOAT drmag;                      // Distance
  FLOAT drsqd;                      // Distance squared
  FLOAT rcentre[ndim];              // Position of sphere centre
  FLOAT volume;                     // Volume of box
  FLOAT *r;                         // Positions of all particles

  // Create local copies of initial conditions parameters
  int Npart = simparams->intparams["Npart"];
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
    sph->sphdata[i].u = press/rhofluid/gammaone;
  }

  delete[] r;

  return;
}



//=============================================================================
//  Simulation::BossBodenheimer
/// Set-up Noh Problem initial conditions
//=============================================================================
template <int ndim>
void Simulation<ndim>::BossBodenheimer(void)
{
  int i;                            // Particle counter
  int j;                            // Aux. particle counter
  int k;                            // Dimension counter
  int Nsphere;                      // Actual number of particles in sphere
  int Nlattice[3];                  // Lattice size
  FLOAT dr[ndim];                   // Relative position vector
  FLOAT drmag;                      // Distance
  FLOAT drsqd;                      // Distance squared
  FLOAT mp;                         // Mass of one particle
  FLOAT rcentre[ndim];              // Position of sphere centre
  FLOAT rho;                        // Fluid density
  FLOAT volume;                     // Volume of box
  FLOAT *r;                         // Positions of all particles
  FLOAT *v;                         // Velocities of all particles

  // Create local copies of initial conditions parameters
  int Npart = simparams->intparams["Npart"];
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

  cout << "ANGVEL : " << angvel*radius*simunits.v.outscale << simunits.v.outunit << endl;

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
  for (i=0; i<Npart; i++) {
    for (k=0; k<ndim; k++) sph->sphdata[i].r[k] = r[ndim*i + k];
    for (k=0; k<ndim; k++) sph->sphdata[i].v[k] = v[ndim*i + k];
    sph->sphdata[i].m = mp;
    sph->sphdata[i].u = temp0/gammaone/mu_bar;
    //if (sph->gas_eos == "isothermal" || sph->gas_eos == "barotropic")
    //  sph->sphdata[i].u = temp0/gammaone/mu_bar;
    //else
    //  sph->sphdata[i].u = press/rho/gammaone;
  }

  delete[] v;
  delete[] r;

  return;
}



//=============================================================================
//  Simulation::PlummerSphere
/// ..
//=============================================================================
template <int ndim>
void Simulation<ndim>::PlummerSphere(void)
{
  bool flag;                        // Aux. flag
  bool istar;                       // Are particles 'stars'?
  int i,j,k;                        // Particle and dimension counter
  int N;                            // ??
  int s;                            // Star counter
  int *porder;                      // ..
  FLOAT dr[ndim];                   // ..
  FLOAT drmag;                      // Distance
  FLOAT drsqd;                      // Distance squared
  FLOAT mp;                         // Mass of particle
  FLOAT raux;                       // Aux. float variable
  FLOAT rcentre[ndim];              // Position of centre of Plummer sphere
  FLOAT vplummer;                   // 
  FLOAT *radsqd;                    // ..

  int idum;
  int mcount = 0;
  FLOAT rpl,rlim,tcr,g,gp,rpold;
  FLOAT mrpl,mrlim;
  FLOAT x1,x2,x3,x4,x5,x6,x7;
  FLOAT rad,vm,ve,t1,t2,w,z;

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
  // ==========================================================================
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
    // ------------------------------------------------------------------------
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
    // ------------------------------------------------------------------------
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
    // ------------------------------------------------------------------------
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
  // ==========================================================================

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

  cout << "Finished generating Plummer sphere" << endl;

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
  int j;                            // Aux. particle counter
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
  FLOAT press = simparams->floatparams["press1"];
  FLOAT mu_bar = simparams->floatparams["mu_bar"];
  FLOAT gammaone = simparams->floatparams["gamma_eos"] - 1.0;
  FLOAT kefrac = simparams->floatparams["kefrac"];
  string particle_dist = simparams->stringparams["particle_distribution"];
  Nlattice[0] = simparams->intparams["Nlattice1[0]"];
  Nlattice[1] = simparams->intparams["Nlattice1[1]"];
  Nlattice[2] = simparams->intparams["Nlattice1[2]"];

  debug2("[Simulation::SedovBlastWave]");


  // Compute size and range of fluid bounding boxes
  // --------------------------------------------------------------------------
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
  cout << "Here??? : " << endl;
  cout << "mbox : " << mbox << "   r_hot : " << r_hot 
       << "   Nbox : " << Nbox << endl;

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
    sph->sphdata[i].u = small_number;
  }

  // Set initial smoothing lengths and create initial ghost particles
  // --------------------------------------------------------------------------
  sph->Nghost = 0;
  sph->Nghostmax = sph->Nsphmax - sph->Nsph;
  sph->Ntot = sph->Nsph;
  for (i=0; i<sph->Nsph; i++) sph->sphdata[i].active = true;
  
  sph->InitialSmoothingLengthGuess();
  sphneib->UpdateTree(sph,*simparams);
  
  // Search ghost particles
  ghosts.SearchGhostParticles(simbox,sph);

  sphneib->UpdateAllSphProperties(sph,nbody);

  // Update neighbour tre
  sphneib->UpdateTree(sph,*simparams);

  // Calculate all SPH properties
  sphneib->UpdateAllSphProperties(sph,nbody);

  // Now calculate which particles are hot
  // --------------------------------------------------------------------------
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

  cout << "utot : " << utot << endl;

  // Normalise the energies
  // --------------------------------------------------------------------------
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
  int j;                            // Aux. particle counter
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
  FLOAT temp0 = simparams->floatparams["temp0"];
  FLOAT mu_bar = simparams->floatparams["mu_bar"];
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
  // --------------------------------------------------------------------------
  volume = (simbox.boxmax[0] - simbox.boxmin[0])*
    (simbox.boxmax[1] - simbox.boxmin[1]);
  Nbox = Nlattice1[0]*Nlattice1[1];

  lambda = simbox.boxmax[1] - simbox.boxmin[1];
  kwave = twopi/lambda;

  // Allocate local and main particle memory
  sph->Nsph = Nbox;
  AllocateParticleMemory();
  r = new FLOAT[ndim*sph->Nsph];
  cout << "Nbox1 : " << Nbox << endl;
  cout << "Allocating memory : " << sph->Nsph << endl;


  // Add particles for LHS of the shocktube
  // --------------------------------------------------------------------------
  if (Nbox > 0) {
    AddCubicLattice(Nbox,Nlattice1,r,simbox,false);
    //AddHexagonalLattice(Nbox,Nlattice1,r,simbox,false);

    for (i=0; i<Nbox; i++) {
      for (k=0; k<ndim; k++) sph->sphdata[i].r[k] = r[ndim*i + k];
      for (k=0; k<ndim; k++) sph->sphdata[i].v[k] = 0.0;
      sph->sphdata[i].v[0] = amp*sin(kwave*sph->sphdata[i].r[1]);
      sph->sphdata[i].m = rhofluid1*volume/(FLOAT) Nbox;
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
  FLOAT volume;                     // Total gas volume
  FLOAT xold;                       // Old x-position (for iteration)
  FLOAT xnew;                       // New x-position (for iteration)
  FLOAT *r;                         // Particle positions

  // Make local copies of parameters for setting up problem
  int Npart = simparams->intparams["Npart"];
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

  cout << "gas_eos : " << sph->gas_eos << endl;
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

  cout << "amp : " << amp << "    rhofluid1" << rhofluid1 << endl;
  cout << "kwave : " << kwave << "    csound : " << csound << endl;

  // Allocate local and main particle memory
  sph->Nsph = Npart;
  Nlattice1[0] = Npart;
  AllocateParticleMemory();
  r = new FLOAT[ndim*sph->Nsph];
  cout << "Allocating memory : " << sph->Nsph << endl;

  AddCubicLattice(Npart,Nlattice1,r,simbox,false);

  // Set positions of all particles to produce density perturbation
  // --------------------------------------------------------------------------
  for (i=0; i<Npart; i++) {
    xnew = r[ndim*i];

    // Solve iterative procedure for particle positions in sound wave
    do {
      xold = xnew;
      xnew = r[ndim*i] - amp*(1.0 - cos(kwave*xnew))/kwave;
      diff = fabs((xnew - xold)/lambda);
    } while(diff > 1.0e-6);

    if (xnew > simbox.boxmax[0]) xnew -= simbox.boxsize[0];
    if (xnew < simbox.boxmin[0]) xnew += simbox.boxsize[0];

    // Set positions in main array with corresponind velocity perturbation
    for (k=0; k<ndim; k++) sph->sphdata[i].r[k] = xnew;
    for (k=0; k<ndim; k++) 
      sph->sphdata[i].v[k] = csound*amp*sin(kwave*xnew);
    sph->sphdata[i].m = rhofluid1*lambda/(FLOAT) Npart;

    if (sph->gas_eos == "isothermal")
      sph->sphdata[i].u = temp0/gammaone/mu_bar;
    else
      sph->sphdata[i].u = press1/rhofluid1/gammaone;
  }
  // --------------------------------------------------------------------------

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
  FLOAT sma = 1.0;                 // ..
  FLOAT eccent = 0.0;              // ..
  FLOAT m1 = 0.5;                  // ..
  FLOAT m2 = 0.5;                  // ..
  DOUBLE rbinary[ndim];            // ..
  DOUBLE vbinary[ndim];            // ..

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
  AddBinaryStar(sma,eccent,m1,m2,0.01,0.01,rbinary,vbinary,
                nbody->stardata[0],nbody->stardata[1]);

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
  FLOAT sma1 = 1.0;                // ..
  FLOAT sma2 = 0.001;              // ..
  FLOAT eccent = 0.0;              // ..
  FLOAT m1 = 0.5;                  // ..
  FLOAT m2 = 0.5;                  // ..
  DOUBLE rbinary[ndim];            // ..
  DOUBLE vbinary[ndim];            // ..
  NbodyParticle<ndim> b1;          // ..
  NbodyParticle<ndim> b2;          // ..

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
  AddBinaryStar(sma1,eccent,m1,m2,0.01,0.01,rbinary,vbinary,b1,b2);

  cout << "b1, r : " << b1.r[0] << "    " << b1.r[1] << endl;
  cout << "b1, v : " << b1.v[0] << "    " << b1.v[1] << endl;
  cout << "b2, r : " << b2.r[0] << "    " << b2.r[1] << endl;
  cout << "b2, v : " << b2.v[0] << "    " << b2.v[1] << endl;

  // Now compute both components
  AddBinaryStar(sma2,eccent,0.5*m1,0.5*m1,0.0001,0.0001,b1.r,b1.v,
		        nbody->stardata[0],nbody->stardata[1]);
  AddBinaryStar(sma2,eccent,0.5*m2,0.5*m2,0.0001,0.0001,b2.r,b2.v,
		        nbody->stardata[2],nbody->stardata[3]);

  for (int i=0; i<nbody->Nstar; i++) {
    cout << "Star : " << i << "   r : " << nbody->stardata[i].r[0] << "  " 
	 << nbody->stardata[i].r[1] << "     v : " << nbody->stardata[i].v[0]
	 << "    " << nbody->stardata[i].v[1] << endl;
  }

  return;
}



//=============================================================================
//  Simulation::AddBinaryStar
/// Create a simple binary star problem
//=============================================================================
template <int ndim>
void Simulation<ndim>::AddBinaryStar
(DOUBLE sma,                       ///< ..
 DOUBLE eccent,                    ///< ..
 DOUBLE m1,                        ///< ..
 DOUBLE m2,                        ///< ..
 DOUBLE h1,                        ///< ..
 DOUBLE h2,                        ///< ..
 DOUBLE *rbinary,                  ///< ..
 DOUBLE *vbinary,                  ///< ..
 NbodyParticle<ndim> &s1,          ///< ..
 NbodyParticle<ndim> &s2)          ///< ..
{
  int i;                           // Particle counter
  int j;                           // Aux. particle counter
  int k;                           // Dimension counter

  FLOAT mbin = m1 + m2;
  FLOAT period = twopi*sqrt(sma*sma*sma/mbin);
  FLOAT vbin = twopi*sma/period;

  debug2("[Simulation::AddBinaryStar]");

  if (ndim == 1) {
    string message = "Binary test not available in 1D";
    ExceptionHandler::getIstance().raise(message);
  }

  cout << "Adding binary with : " << m1 << "   " << m2 << "    " << sma << endl;

  // Set properties of star 1
  for (k=0; k<ndim; k++) s1.r[k] = rbinary[k];
  for (k=0; k<ndim; k++) s1.v[k] = vbinary[k];
  s1.m = m1;
  s1.h = h1;
  s1.invh = 1.0 / s1.h;
  s1.r[0] += sma*m2/mbin;
  s1.v[1] += vbin*m2/mbin;

  // Set properties of star 2
  for (k=0; k<ndim; k++) s2.r[k] = rbinary[k];
  for (k=0; k<ndim; k++) s2.v[k] = vbinary[k];
  s2.m = m2;
  s2.h = h2;
  s2.invh = 1.0 / s2.h;
  s2.r[0] -= sma*m1/mbin;
  s2.v[1] -= vbin*m1/mbin;

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
  // --------------------------------------------------------------------------
  for (i=0; i<Npart; i++) {

    // Continously loop until random particle lies inside sphere
    do {
      for (k=0; k<ndim; k++) 
	rpos[k] = 1.0 - 2.0*(FLOAT)(rand()%RAND_MAX)/(FLOAT)RAND_MAX;
      rad = DotProduct(rpos,rpos,ndim);
    } while (rad > radius);

    for (k=0; k<ndim; k++) r[ndim*i + k] = rcentre[k] + rpos[k];
  }
  // --------------------------------------------------------------------------

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

  debug2("[Simulation::AddCubicLattice]");
  
  // Create lattice depending on dimensionality
  // --------------------------------------------------------------------------
  if (ndim == 1) {
    for (ii=0; ii<Nlattice[0]; ii++) {
      i = ii;
      r[i] = box.boxmin[0] + ((FLOAT)ii + 0.5)*spacing[0];
    }
  }
  // --------------------------------------------------------------------------
  else if (ndim == 2) {
    for (jj=0; jj<Nlattice[1]; jj++) {
      for (ii=0; ii<Nlattice[0]; ii++) {
	i = jj*Nlattice[0] + ii;
	r[ndim*i] = box.boxmin[0] + ((FLOAT)ii + 0.5)*spacing[0];
	r[ndim*i + 1] = box.boxmin[1] + ((FLOAT)jj + 0.5)*spacing[1];
      }
    }
  }
  // --------------------------------------------------------------------------
  else if (ndim == 3) {
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
  // --------------------------------------------------------------------------
  if (ndim == 1) {
    for (ii=0; ii<Nlattice[0]; ii++) {
      i = ii;
      r[i] = box.boxmin[0] + 0.5*rad[0] + 2.0*(FLOAT)ii*rad[0];
    }
  }

  // --------------------------------------------------------------------------
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

  // --------------------------------------------------------------------------
  else if (ndim == 3) {
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
  // --------------------------------------------------------------------------
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
  // --------------------------------------------------------------------------
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
  // --------------------------------------------------------------------------
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
  // --------------------------------------------------------------------------

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
  FLOAT spacing;                    // ..

  debug2("[Simulation::AddAzimuthalDensityPerturbation]");


  // Loop over all required particles
  // --------------------------------------------------------------------------
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
  // --------------------------------------------------------------------------

  return;
}
