//=================================================================================================
//  MeshlessFVSimulation.cpp
//  Contains all main functions controlling Meshless Finite-Volume simulation work-flow.
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


#include <iostream>
#include <iomanip>
#include <ostream>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>
#include <time.h>
#include <cstdio>
#include <cstring>
#include "Precision.h"
#include "CodeTiming.h"
#include "Exception.h"
#include "Debug.h"
#include "InlineFuncs.h"
#include "Simulation.h"
#include "Parameters.h"
#include "Nbody.h"
#include "RandomNumber.h"
#include "RiemannSolver.h"
#include "Ghosts.h"
#include "Sinks.h"
using namespace std;



// Create template class instances of the main MeshlessFVSimulation object for
// each dimension used (1, 2 and 3)
template class MeshlessFVSimulation<1>;
template class MeshlessFVSimulation<2>;
template class MeshlessFVSimulation<3>;



// Construct the meshless integration object.
template<int ndim, class MeshlessType>
MeshlessFV<ndim>* MeshlessFactoryFull
(Parameters* simparams,
 SimUnits& simunits)
 {
  // Local references to parameter variables for brevity
  map<string, int> &intparams = simparams->intparams;
  map<string, double> &floatparams = simparams->floatparams;
  map<string, string> &stringparams = simparams->stringparams;

  string KernelName = stringparams["kernel"];

  return new MeshlessType
      (intparams["hydro_forces"], intparams["self_gravity"], floatparams["accel_mult"],
          floatparams["courant_mult"], floatparams["h_fac"], floatparams["h_converge"],
          floatparams["gamma_eos"], stringparams["gas_eos"], KernelName,
          sizeof(MeshlessFVParticle<ndim>), simunits, simparams);
 }


// Select the slope limiter
template<int ndim, template<int> class Kernel, template<int, template<int> class, class> class MeshlessType>
MeshlessFV<ndim>*  _MeshlessFactorySlopes
(Parameters* simparams,
 SimUnits& simunits)
 {
  string limiter = simparams->stringparams["slope_limiter"];

  if (limiter == "null") {
    typedef NullLimiter<ndim> Limiter;
    return MeshlessFactoryFull<ndim, MeshlessType<ndim,Kernel,Limiter> >(simparams,simunits) ;
  }
  else if (limiter == "zeroslope") {
    typedef ZeroSlopeLimiter<ndim> Limiter;
    return MeshlessFactoryFull<ndim, MeshlessType<ndim,Kernel,Limiter> >(simparams,simunits) ;
  }
  else if (limiter == "tvdscalar" || limiter == "tess2011") {
    typedef TVDScalarLimiter<ndim> Limiter;
    return MeshlessFactoryFull<ndim, MeshlessType<ndim,Kernel,Limiter> >(simparams,simunits) ;
  }
  else if (limiter == "scalar" || limiter == "balsara2004") {
    typedef ScalarLimiter<ndim> Limiter;
    return MeshlessFactoryFull<ndim, MeshlessType<ndim,Kernel,Limiter> >(simparams,simunits) ;
  }
  else if (limiter == "springel2009") {
    typedef Springel2009Limiter<ndim> Limiter;
    return MeshlessFactoryFull<ndim, MeshlessType<ndim,Kernel,Limiter> >(simparams,simunits) ;
  }
  else if (limiter == "gizmo") {
    typedef GizmoLimiter<ndim> Limiter;
    return MeshlessFactoryFull<ndim, MeshlessType<ndim,Kernel,Limiter> >(simparams,simunits) ;
  }

  string message = "Unrecognised parameter : slope_limiter = " + limiter;
  ExceptionHandler::getIstance().raise(message);

  return NULL ;
 }

// Select the Meshless time integration type.
template<int ndim, template<int> class Kernel>
MeshlessFV<ndim>* _MeshlessTimeIntegFactory
(Parameters* simparams,
 SimUnits& simunits)
{
  string sim = simparams->stringparams["sim"];

  if (sim == "meshlessfv" || sim == "mfvmuscl") {
    return _MeshlessFactorySlopes<ndim, Kernel, MfvMuscl>(simparams, simunits) ;
  }
  else if (sim == "mfvrk") {
    return _MeshlessFactorySlopes<ndim, Kernel, MfvRungeKutta>(simparams, simunits) ;
  }

  string message = "Invalid option for the simulation type parameter: " + sim;
  ExceptionHandler::getIstance().raise(message);

  return NULL ;
}


// Create the full Meshless forces object in a type by type way.
template <int ndim>
MeshlessFV<ndim>* MeshlessFactory
(Parameters* simparams,
    SimUnits& simunits)
    {
  // Local references to parameter variables for brevity
  map<string, int> &intparams = simparams->intparams;
  map<string, string> &stringparams = simparams->stringparams;
  string KernelName = stringparams["kernel"];

  if (intparams["tabulated_kernel"] == 1) {
    return _MeshlessTimeIntegFactory<ndim, TabulatedKernel>(simparams, simunits) ;
  }
  else {
    if (KernelName == "m4") {
      return _MeshlessTimeIntegFactory<ndim, M4Kernel>(simparams, simunits) ;
    }
    else if (KernelName == "quintic") {
      return _MeshlessTimeIntegFactory<ndim, QuinticKernel>(simparams, simunits) ;
    }
  }

  string message = "Unrecognised parameter : kernel = " + simparams->stringparams["kernel"];
  ExceptionHandler::getIstance().raise(message);

  return NULL ;

    }


//=================================================================================================
//  MeshlessFVSimulation::ProcessParameters
/// Process all the options chosen in the parameters file, setting various
/// simulation variables and creating important simulation objects.
/// Meshless specific version
//=================================================================================================
template <int ndim>
void MeshlessFVSimulation<ndim>::ProcessParameters(void)
{
  // Local references to parameter variables for brevity
  map<string, int> &intparams = simparams->intparams;
  map<string, double> &floatparams = simparams->floatparams;
  map<string, string> &stringparams = simparams->stringparams;
  string sim = stringparams["sim"];
  string gas_eos = stringparams["gas_eos"];
  string gas_radiation = stringparams["radiation"];
  string KernelName = stringparams["kernel"];

  debug2("[MeshlessFVSimulation::ProcessParameters]");

  // Common set-up
  Simulation<ndim>::ProcessParameters();

  if (intparams["self_gravity"] && !intparams["zero_mass_flux"]) {
    string message = "The meshless scheme only works with self-gravity when MFM is activated";
    ExceptionHandler::getIstance().raise(message);
  }

  if (intparams["sink_particles"] && !intparams["zero_mass_flux"]) {
    string message = "The meshless scheme only works with sink particles when MFM is activated";
    ExceptionHandler::getIstance().raise(message);
  }


#ifdef MPI_PARALLEL
  if (stringparams["mpi_decomposition"] == "kdtree") {
    mpicontrol = new MpiKDTreeDecomposition<ndim,MeshlessFVParticle>(simbox);
  }
  else {
    string message = "Unrecognised parameter : mpi_decomposition = "
        + simparams->stringparams["mpi_decomposition"];
    ExceptionHandler::getIstance().raise(message);
  }
  mpicontrol->timing = timing;
  rank = mpicontrol->rank;
  Nmpi = mpicontrol->Nmpi;
#endif

  // Create Meshless Finite-Volume object depending on choice of kernel
  //===============================================================================================
  hydro = mfv = MeshlessFactory<ndim>(simparams, simunits) ;


  // Meshless-time integration object
  hydroint = new MfvIntegration<ndim, MeshlessFVParticle>(simparams);

  // Set-up the cooling / energy integration object
  //===============================================================================================
  radfb = NULL ;
  if (stringparams["energy_integration"] == "radws") {
    // Radiative feedback object
    if (intparams["rad_fb"])
      radfb = new RadiativeFB<ndim>(&simunits, simparams);

    uint = new EnergyRadws<ndim, MeshlessFVParticle>
        (simparams, &simunits, (Radws<ndim> *)mfv->eos, radfb);
  }
  else if (stringparams["energy_integration"] == "null" ||
      stringparams["energy_integration"] == "none") {
    uint = new NullEnergy<ndim>(floatparams["energy_mult"]);
  }
  else {
    string message = "Unrecognised parameter : energy_integration = "
        + simparams->stringparams["energy_integration"];
    ExceptionHandler::getIstance().raise(message);
  }

  // Create neighbour searching object based on chosen method in params file
  //-----------------------------------------------------------------------------------------------
  string tree_type = stringparams["neib_search"] ;

  mfvneib = new MeshlessFVTree<ndim,MeshlessFVParticle>
  (tree_type, intparams["Nleafmax"], Nmpi, intparams["pruning_level_min"], intparams["pruning_level_max"],
      floatparams["thetamaxsqd"], hydro->kernp->kernrange, floatparams["macerror"],
      stringparams["gravity_mac"], stringparams["multipole"], &simbox, mfv->kernp, timing, mfv->types);

  neib = mfvneib ;
  // Here I do a horrible hack to get at the underlying tree, needed for the dust.
  TreeBase<ndim>
    * t  = mfvneib->GetTree(),
    * gt = mfvneib->GetGhostTree(),
    *mpit = NULL ;
#ifdef MPI_PARALLEL
  mpit = mfvneib->GetMPIGhostTree();
#endif

  // Setup the dust
  mfvdust = DustFactory<ndim, MeshlessFVParticle>::ProcessParameters(simparams, timing, simunits,
                                                                     mfv->types, simbox,
                                                                     t, gt, mpit) ;

#if defined MPI_PARALLEL
 if (mfvdust != NULL) mfvdust->SetMpiControl(mpicontrol);
#endif

  // Depending on the dimensionality, calculate expected neighbour number
  //-----------------------------------------------------------------------------------------------
  if (ndim == 1) {
    mfv->Ngather = (int) (2.0*mfv->kernp->kernrange*mfv->h_fac);
  }
  else if (ndim == 2) {
    mfv->Ngather = (int) (pi*pow(mfv->kernp->kernrange*mfv->h_fac,2));
  }
  else if (ndim == 3) {
    mfv->Ngather = (int) (4.0*pi*pow(mfv->kernp->kernrange*mfv->h_fac,3)/3.0);
  }


  // Process all N-body parameters and set-up main N-body objects
  this->ProcessNbodyParameters();


  // Set pointers to external potential field object
  mfv->extpot = extpot;
  nbody->extpot = extpot;

  // Set eos pointer to nbody
  mfv->eos->set_nbody_data(nbody);



  // Set all other hydro parameter variables
  mfv->Nhydromax       = intparams["Nhydromax"];
  mfv->create_sinks    = intparams["create_sinks"];
  mfv->fixed_sink_mass = intparams["fixed_sink_mass"];
  mfv->msink_fixed     = floatparams["m1"];


  // Set important variables for N-body objects
  //nbody->Nstar          = intparams["Nstar"];
  nbody->Nstarmax       = intparams["Nstarmax"];
  nbody_single_timestep = intparams["nbody_single_timestep"];
  nbodytree.gpehard     = floatparams["gpehard"];
  nbodytree.gpesoft     = floatparams["gpesoft"];
  //nbody->perturbers     = intparams["perturbers"];
  //if (intparams["sub_systems"] == 1) subsystem->perturbers = intparams["perturbers"];

  // Sink particles
  //-----------------------------------------------------------------------------------------------
  sinks = new Sinks<ndim>(mfvneib);
  sink_particles             = intparams["sink_particles"];
  sinks->sink_particles      = intparams["sink_particles"];
  sinks->create_sinks        = intparams["create_sinks"];
  sinks->smooth_accretion    = intparams["smooth_accretion"];
  sinks->Nsinkfixed          = intparams["Nsinkfixed"];
  sinks->alpha_ss            = floatparams["alpha_ss"];
  sinks->smooth_accrete_frac = floatparams["smooth_accrete_frac"];
  sinks->smooth_accrete_dt   = floatparams["smooth_accrete_dt"];
  sinks->sink_radius_mode    = stringparams["sink_radius_mode"];
  sinks->rho_sink            = floatparams["rho_sink"];
  sinks->rho_sink            /= simunits.rho.outscale/simunits.rho.outcgs;

  if (sinks->sink_radius_mode == "fixed") {
    sinks->sink_radius = floatparams["sink_radius"]/simunits.r.outscale;
  }
  else {
    sinks->sink_radius = floatparams["sink_radius"];
  }

  // Sanity-check for various sink particle values
  if (intparams["sink_particles"] == 1 &&
      (stringparams["nbody"] != "lfkdk" && stringparams["nbody"] != "lfdkd")) {
    string message = "Invalid parameter : nbody must use lfkdk or lfdkd when "
        "using accreting sink particles";
    ExceptionHandler::getIstance().raise(message);
  }

  // Supernova feedback
  //-----------------------------------------------------------------------------------------------
  if (stringparams["supernova_feedback"] == "none") {
    snDriver = new NullSupernovaDriver<ndim>(this);
  }
  else if (stringparams["supernova_feedback"] == "single") {
    snDriver = new SedovTestDriver<ndim>(this,simparams, simunits);
  }
  else if (stringparams["supernova_feedback"] == "random") {
    snDriver = new RandomSedovTestDriver<ndim>(this,simparams, simunits, simbox);
  }
  else {
    string message = "Unrecognised parameter : external_potential = "
      + simparams->stringparams["supernova_feedback"];
    ExceptionHandler::getIstance().raise(message);
  }

  if (radfb) radfb->SetSinks(sinks);


  // Set pointers to timing object
  nbody->timing   = timing;
  //if (sim == "sph" || sim == "gradhsph" || sim == "sm2012sph" || sim == "godunov_hydro") {
  sinks->timing    = timing;
  mfvneib->SetTimingObject(timing);
  mfv->timing = timing;
  hydroint->timing = timing;
  uint->timing = timing;
  //}*/

  // Create ghost particle object
  //-----------------------------------------------------------------------------------------------
  if (IsAnyBoundarySpecial(simbox)) {
    LocalGhosts = new PeriodicGhostsSpecific<ndim,MeshlessFVParticle >();
  }
  else {
    LocalGhosts = new NullGhosts<ndim>();
  }


#if defined MPI_PARALLEL
  mpicontrol->SetNeibSearch(mfvneib);
  sinks->SetMpiControl(mpicontrol);
  if (stringparams["out_file_form"]=="sf") {
    string message = "The sf format is not supported with MPI! Use the column "
        "or (better) the su format";
    ExceptionHandler::getIstance().raise(message);
  }
  MpiGhosts = new MpiGhostsSpecific<ndim, MeshlessFVParticle>(mpicontrol);
#endif

  time_step_limiter_type = stringparams["time_step_limiter"] ;

  // Flag that we've processed all parameters already
  ParametersProcessed = true;


  return;
}



//=================================================================================================
//  MeshlessFVSimulation::PostInitialConditionsSetup
/// Call routines for calculating all initial hydro and N-body quantities
/// once initial conditions have been set-up.
//=================================================================================================
template <int ndim>
void MeshlessFVSimulation<ndim>::PostInitialConditionsSetup(void)
{
  int i;                               // Particle counter
  int k;                               // Dimension counter

  debug2("[MeshlessFVSimulation::PostInitialConditionsSetup]");

  if (nbody->Nstar > 0 && !simparams->intparams["zero_mass_flux"]) {
    string message = "The meshless scheme only works with star particles when MFM is activated";
    ExceptionHandler::getIstance().raise(message);
  }

  // Set iorig
  if (rank == 0) {
    MeshlessFVParticle<ndim> *partdata = mfv->GetMeshlessFVParticleArray();
    for (i=0; i<mfv->Nhydro; i++) {
      if (not restart) {
        partdata[i].iorig = i;
      }
      partdata[i].flags.set(active);
      partdata[i].flags.set(update_density) ;
      partdata[i].gpot = big_number;
    }
  }

  // Perform initial MPI decomposition
  //-----------------------------------------------------------------------------------------------
#ifdef MPI_PARALLEL
  mpicontrol->CreateInitialDomainDecomposition(mfv, nbody, simparams, this->initial_h_provided);
  this->AllocateParticleMemory();
#endif

  // Copy information from the stars to the sinks
  if (simparams->intparams["sink_particles"]==1) {
    sinks->Nsink = nbody->Nstar;
    sinks->AllocateMemory(sinks->Nsink);
    for (int i=0; i<sinks->Nsink; i++) {
      sinks->sink[i].star   = &(nbody->stardata[i]);
      sinks->sink[i].istar  = i;
      sinks->sink[i].radius = hydro->kernp->kernrange*nbody->stardata[i].h;
      //sinks->sink[i].radius = simparams->floatparams["sink_radius"];
    }
  }


  // Set time variables here (for now)
  nresync = 0;   // DAVID : Need to adapt this for block timesteps
  n = 0;
  integration_step = 1;
  level_step = 1;

  // Set all relevant particle counters
  mfv->Nghost = 0;
  mfv->Ntot = mfv->Nhydro;

  // Initialise conserved variables
  for (i=0; i<mfv->Nhydro; i++) {
    MeshlessFVParticle<ndim>& part = mfv->GetMeshlessFVParticlePointer(i);
    part.Qcons0[MeshlessFV<ndim>::irho] = part.m;
    for (int k=0; k<ndim; k++) part.Qcons0[k] = part.m*part.v[k];
    FLOAT ekin = (FLOAT) 0.0;
    for (int k=0; k<ndim; k++) ekin += part.v[k]*part.v[k];
    part.Qcons0[MeshlessFV<ndim>::ietot] = part.u*part.m + (FLOAT) 0.5*part.m*ekin;
  }



  // If the smoothing lengths have not been provided beforehand, then
  // calculate the initial values here

  mfvneib->ToggleNeighbourCheck(false);
  if (!this->initial_h_provided) {
    mfv->InitialSmoothingLengthGuess();
    mfvneib->BuildTree(rebuild_tree, 0, ntreebuildstep, ntreestockstep, timestep, mfv);
    mfvneib->UpdateAllProperties(mfv, nbody);

    for (i=0; i<mfv->Nhydro; i++) {
      MeshlessFVParticle<ndim>& part = mfv->GetMeshlessFVParticlePointer(i);
      part.flags.set(update_density) ;
    }
  }
  else {
    mfvneib->BuildTree(rebuild_tree, 0, ntreebuildstep, ntreestockstep, timestep, mfv);
  }

#ifdef MPI_PARALLEL
  mpicontrol->UpdateAllBoundingBoxes(mfv->Nhydro, mfv, mfv->kernp);
#endif

  // Search ghost particles
  mfvneib->SearchBoundaryGhostParticles((FLOAT) 0.0, simbox, mfv);
  mfvneib->BuildGhostTree(true, 0, ntreebuildstep, ntreestockstep, timestep, mfv);
#ifdef MPI_PARALLEL
  mpicontrol->UpdateAllBoundingBoxes(mfv->Nhydro+mfv->NPeriodicGhost, mfv, mfv->kernp);
  for (int i=0; i<mfv->Nhydro+mfv->NPeriodicGhost; i++) {
    MeshlessFVParticle<ndim>& parti =  mfv->GetMeshlessFVParticlePointer(i);
    parti.hrangesqd = mfv->kernp->kernrangesqd*parti.h*parti.h;
  }
  MpiGhosts->SearchGhostParticles((FLOAT) 0.0, simbox, mfv);
  mfvneib->BuildMpiGhostTree(true, 0, ntreebuildstep, ntreestockstep, timestep, mfv);
#endif

  // Zero accelerations
  for (i=0; i<mfv->Nhydro; i++) {
    MeshlessFVParticle<ndim>& part = mfv->GetMeshlessFVParticlePointer(i);
    part.flags.set(active);
    part.flags.set(update_density) ;
    part.gpot = big_number;
  }

  // Update neighbour tree
  rebuild_tree = true;
  mfvneib->BuildTree(rebuild_tree, 0, ntreebuildstep, ntreestockstep, timestep, mfv);
#ifdef MPI_PARALLEL
  mfvneib->InitialiseCellWorkCounters();
#endif


  // Update neighbour tree
  rebuild_tree = true;
  mfvneib->BuildTree(rebuild_tree, 0, ntreebuildstep, ntreestockstep, timestep, mfv);

  // Calculate all hydro properties
  mfvneib->UpdateAllProperties(mfv, nbody);

#ifdef MPI_PARALLEL
  mpicontrol->UpdateAllBoundingBoxes(mfv->Nhydro, mfv, mfv->kernp);
#endif

  // Search ghost particles
  mfvneib->SearchBoundaryGhostParticles((FLOAT) 0.0, simbox, mfv);
  mfvneib->BuildGhostTree(true, 0, ntreebuildstep, ntreestockstep, timestep, mfv);
#ifdef MPI_PARALLEL
  mpicontrol->UpdateAllBoundingBoxes(mfv->Nhydro + mfv->NPeriodicGhost, mfv, mfv->kernp);
  MpiGhosts->SearchGhostParticles((FLOAT) 0.0, simbox, mfv);
  mfvneib->BuildMpiGhostTree(true, 0, ntreebuildstep, ntreestockstep, timestep, mfv);
#endif


  // Update neighbour tree
  rebuild_tree = true;
  mfvneib->BuildTree(rebuild_tree, 0, ntreebuildstep, ntreestockstep, timestep, mfv);
#ifdef MPI_PARALLEL
  mfvneib->InitialiseCellWorkCounters();
#endif
  mfvneib->SearchBoundaryGhostParticles((FLOAT) 0.0, simbox, mfv);
  mfvneib->BuildGhostTree(true, 0, ntreebuildstep, ntreestockstep, timestep, mfv);
  // Communicate pruned trees for MPI
#ifdef MPI_PARALLEL
  mfvneib->BuildPrunedTree(rank, simbox, mpicontrol->mpinode, mfv);
#endif
  mfvneib->ToggleNeighbourCheck(true);

  // Compute mean mass (used for smooth sink accretion)
  if (!restart) {
    MeshlessFVParticle<ndim> *partdata = mfv->GetMeshlessFVParticleArray();
    mfv->mmean = (FLOAT) 0.0;
    for (i=0; i<mfv->Nhydro; i++) mfv->mmean += partdata[i].m;
    mfv->mmean /= (FLOAT) mfv->Nhydro;
  }

  // Compute all initial N-body terms
  //-----------------------------------------------------------------------------------------------
  for (i=0; i<nbody->Nstar; i++) {
    for (k=0; k<ndim; k++) nbody->stardata[i].a[k]     = (FLOAT) 0.0;
    for (k=0; k<ndim; k++) nbody->stardata[i].adot[k]  = (FLOAT) 0.0;
    for (k=0; k<ndim; k++) nbody->stardata[i].a2dot[k] = (FLOAT) 0.0;
    for (k=0; k<ndim; k++) nbody->stardata[i].a3dot[k] = (FLOAT) 0.0;
    nbody->stardata[i].gpot   = (FLOAT) 0.0;
    nbody->stardata[i].gpe    = (FLOAT) 0.0;
    nbody->stardata[i].tlast  = t;
    nbody->stardata[i].flags.set(active);
    nbody->stardata[i].level  = 0;
    nbody->stardata[i].nstep  = 0;
    nbody->stardata[i].nlast  = 0;
    nbody->nbodydata[i]       = &(nbody->stardata[i]);
  }
  nbody->Nnbody = nbody->Nstar;

  // Read-in N-body table here
  nbody->LoadStellarPropertiesTable(&simunits);
  nbody->UpdateStellarProperties();

  for (i=0; i<mfv->Ntot; i++) {
    MeshlessFVParticle<ndim>& part = mfv->GetMeshlessFVParticlePointer(i);
    for (k=0; k<ndim; k++) {
      part.a[k] = (FLOAT) 0.0;
      part.atree[k] = (FLOAT) 0.0;
      part.rdmdt[k] = 0.0;
      part.rdmdt0[k] = 0.0;
    }
    for (k=0; k<ndim+2; k++) part.dQ[k] = (FLOAT) 0.0;
    for (k=0; k<ndim+2; k++) part.dQdt[k] = (FLOAT) 0.0;
    part.level  = 0;
    part.nstep  = 0;
    part.nlast  = 0;
    part.tlast  = t;
    part.dt     = 0;
    part.flags.set(active);
  }

  // Compute all initial hydro terms
  // We will need to iterate if we are going to use a relative opening criterion
  //-----------------------------------------------------------------------------------------------
  MAC_Type mac_type = mfvneib->GetOpeningCriterion() ;
  int n_iter = 1 + (mfv->self_gravity == 1 && mac_type != geometric) ;
  for (int iter=0; iter < n_iter; iter++) {

    if (iter==0 && mac_type != geometric)
      mfvneib->SetOpeningCriterion(geometric);
    else
      mfvneib->SetOpeningCriterion(mac_type) ;

    // Copy all other data from real hydro particles to ghosts
    LocalGhosts->CopyHydroDataToGhosts(simbox,mfv);

    mfvneib->BuildTree(true, 0, ntreebuildstep, ntreestockstep, timestep, mfv);
#ifdef MPI_PARALLEL
    mfvneib->InitialiseCellWorkCounters();
#endif
    mfvneib->SearchBoundaryGhostParticles((FLOAT) 0.0, simbox, mfv);
    mfvneib->BuildGhostTree(true, 0, ntreebuildstep, ntreestockstep, timestep, mfv);

#ifdef MPI_PARALLEL
    mpicontrol->UpdateAllBoundingBoxes(mfv->Nhydro + mfv->NPeriodicGhost, mfv, mfv->kernp);
    MpiGhosts->SearchGhostParticles((FLOAT) 0.0, simbox, mfv);
    mfvneib->BuildMpiGhostTree(true, 0, ntreebuildstep, ntreestockstep, timestep, mfv);
#endif

    if (iter == 0) {
      mfvneib->UpdateAllProperties(mfv, nbody);

      LocalGhosts->CopyHydroDataToGhosts(simbox,mfv);
#ifdef MPI_PARALLEL
      MpiGhosts->CopyHydroDataToGhosts(simbox,mfv);
#endif

      mfvneib->UpdateGradientMatrices(mfv, nbody, simbox);
      LocalGhosts->CopyHydroDataToGhosts(simbox,mfv);
#ifdef MPI_PARALLEL
      MpiGhosts->CopyHydroDataToGhosts(simbox,mfv);
#endif
    }

    if (mfv->self_gravity == 1 || nbody->Nnbody > 0) {
      mfv->ZeroAccelerations();

#ifdef MPI_PARALLEL
      if (mfv->self_gravity ==1 ) {
        mfvneib->UpdateGravityExportList(rank, mfv, nbody, simbox, ewald);
        mpicontrol->ExportParticlesBeforeForceLoop(mfv);
      }
#endif
      mfvneib->UpdateAllGravForces(mfv, nbody, simbox, ewald);
#ifdef MPI_PARALLEL
      if (mfv->self_gravity ==1 ) {
        mpicontrol->GetExportedParticlesAccelerations(mfv);
      }
#endif
    }

    /*
    mfvneib->UpdateGradientMatrices(mfv, nbody, simbox);
    LocalGhosts->CopyHydroDataToGhosts(simbox,mfv);
#ifdef MPI_PARALLEL
    MpiGhosts->CopyHydroDataToGhosts(simbox,mfv);
#endif
     */
  } // End of force iteration.

  // ..
  for (i=0; i<mfv->Nhydro; i++) {
    MeshlessFVParticle<ndim>& part = mfv->GetMeshlessFVParticlePointer(i);
    for (k=0; k<ndim; k++) part.r0[k] = part.r[k];
    for (k=0; k<ndim; k++) part.v0[k] = part.v[k];
    for (k=0; k<ndim; k++) part.a0[k] = part.a[k];
  }

  // Compute initial N-body forces
  //-----------------------------------------------------------------------------------------------
  if (mfv->self_gravity == 1 && mfv->Nhydro > 0) {
    mfvneib->UpdateAllStarGasForces(mfv, nbody, simbox, ewald);  //, simbox, ewald);

    // We need to sum up the contributions from the different domains
#if defined MPI_PARALLEL
    mpicontrol->ComputeTotalStarGasForces(nbody);
#endif
  }

  if (nbody->nbody_softening == 1) {
    nbody->CalculateDirectSmoothedGravForces(nbody->Nnbody, nbody->nbodydata, simbox, ewald);
  }
  else {
    nbody->CalculateDirectGravForces(nbody->Nnbody, nbody->nbodydata, simbox, ewald);
  }
  nbody->CalculateAllStartupQuantities(nbody->Nnbody, nbody->nbodydata, simbox, ewald);


  for (i=0; i<nbody->Nnbody; i++) {
    if (nbody->nbodydata[i]->flags.check(active)) {
      nbody->extpot->AddExternalPotential(nbody->nbodydata[i]->r, nbody->nbodydata[i]->v,
          nbody->nbodydata[i]->a, nbody->nbodydata[i]->adot,
          nbody->nbodydata[i]->gpot);
    }
  }

  // Compute the dust forces if present.
  if (mfvdust != NULL) {

    // Copy properties from original particles to ghost particles
    LocalGhosts->CopyHydroDataToGhosts(simbox,mfv);
#ifdef MPI_PARALLEL
    MpiGhosts->CopyHydroDataToGhosts(simbox, mfv);
#endif
    mfvdust->UpdateAllDragForces(mfv) ;

    for (i=0; i<mfv->Nhydro; i++) {
      MeshlessFVParticle<ndim>& part = mfv->GetMeshlessFVParticlePointer(i);
      for (k=0; k<ndim+2; k++) part.dQdt[k] = 0 ;
    }
  }


  // Compute timesteps for all particles
  if (Nlevels == 1) {
    this->ComputeGlobalTimestep();
  }
  else {
    if (time_step_limiter_type == "conservative") {
      mfvneib->UpdateTimestepsLimitsFromDistantParticles(mfv,false);
  #ifdef MPI_PARALLEL
      mpicontrol->ExportParticlesBeforeForceLoop(mfv);
      mfvneib->UpdateTimestepsLimitsFromDistantParticles(mfv,true);
      mpicontrol->GetExportedParticlesAccelerations(mfv);
  #endif
    }
    this->ComputeBlockTimesteps();
  }


  // Call EndTimestep to set all 'beginning-of-step' variables
  uint->EndTimestep(n, t, timestep, mfv);
  hydroint->EndTimestep(n, t, timestep, mfv);
  nbody->EndTimestep(n, nbody->Nstar, t, timestep, nbody->nbodydata);

  this->CalculateDiagnostics();
  this->diag0 = this->diag;
  this->setup = true;

  rebuild_tree=false;

  return;
}



//=================================================================================================
//  MeshlessFVSimulation::WriteExtraSinkOutput
/// For any simulations loaded into memory via a snapshot file, all particle
/// variables are converted into dimensionless code units here.
//=================================================================================================
template <int ndim>
void MeshlessFVSimulation<ndim>::WriteExtraSinkOutput(void)
{
  return;
}



//=================================================================================================
//  MeshlessFVSimulation::WriteExtraSinkOutput
/// For any simulations loaded into memory via a snapshot file, all particle
/// variables are converted into dimensionless code units here.
//=================================================================================================
template <int ndim>
void MeshlessFVSimulation<ndim>::FinaliseSimulation(void)
{
  /*
  MeshlessFVParticle<ndim> *partdata = mfv->GetMeshlessFVParticleArray();

  for (int i=0; i<mfv->Nhydro; i++) {
    MeshlessFVParticle<ndim> &part = partdata[i];
    if (part.flags.is_dead()) continue;
    // TODO: Check this.
    //   Qcons is now interpretted as the predicted value of Q at the current time, surely this
    //   can be used instead without updating?
    for (int var=0; var<ndim+2; var++) part.Qcons[var] = part.Qcons0[var] + part.dQ[var];
    mfv->ConvertConservedToPrimitive(part.ndens, part.Qcons, part.Wprim);
    mfv->UpdateArrayVariables(part);
  }
   */
  return;
}
