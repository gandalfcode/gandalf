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

  // Create neighbour searching object based on chosen method in params file
  //-----------------------------------------------------------------------------------------------
  string tree_type = stringparams["neib_search"] ;

  mfvneib = new MeshlessFVTree<ndim,MeshlessFVParticle>
  (tree_type, intparams["Nleafmax"], Nmpi, intparams["pruning_level_min"], intparams["pruning_level_max"],
      floatparams["thetamaxsqd"], hydro->kernp->kernrange, floatparams["macerror"],
      stringparams["gravity_mac"], stringparams["multipole"], &simbox, mfv->kernp, timing, mfv->types);


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


  // Set pointers to timing object
  nbody->timing   = timing;
  //if (sim == "sph" || sim == "gradhsph" || sim == "sm2012sph" || sim == "godunov_hydro") {
  sinks->timing    = timing;
  mfvneib->SetTimingObject(timing);
  mfv->timing = timing;
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

  // Set iorig
  if (rank == 0) {
    MeshlessFVParticle<ndim> *partdata = mfv->GetMeshlessFVParticleArray();
    for (i=0; i<mfv->Nhydro; i++) {
      partdata[i].iorig = i;
      partdata[i].flags.set_flag(active);
      partdata[i].flags.set_flag(update_density) ;
      partdata[i].gpot = big_number;
    }
  }

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

  // Perform initial MPI decomposition
  //-----------------------------------------------------------------------------------------------
#ifdef MPI_PARALLEL
  mpicontrol->CreateInitialDomainDecomposition(mfv, nbody, simparams, simbox,
      this->initial_h_provided);
  this->AllocateParticleMemory();
#endif

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
    mfvneib->UpdateAllProperties(mfv, nbody, simbox);

    for (i=0; i<mfv->Nhydro; i++) {
      MeshlessFVParticle<ndim>& part = mfv->GetMeshlessFVParticlePointer(i);
      part.flags.set_flag(update_density) ;
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
    part.flags.set_flag(active);
    part.flags.set_flag(update_density) ;
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
  mfvneib->UpdateAllProperties(mfv, nbody, simbox);

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
    nbody->stardata[i].active = true;
    nbody->stardata[i].level  = 0;
    nbody->stardata[i].nstep  = 0;
    nbody->stardata[i].nlast  = 0;
    nbody->nbodydata[i]       = &(nbody->stardata[i]);
  }
  nbody->Nnbody = nbody->Nstar;

  // Read-in N-body table here
  nbody->LoadStellarPropertiesTable(&simunits);
  nbody->UpdateStellarProperties();


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

    for (i=0; i<mfv->Ntot; i++) {
      MeshlessFVParticle<ndim>& part = mfv->GetMeshlessFVParticlePointer(i);
      for (k=0; k<ndim; k++) {
        part.a[k] = (FLOAT) 0.0;
        part.rdmdt[k] = 0.0;
        part.rdmdt0[k] = 0.0;
      }
      for (k=0; k<ndim+2; k++) part.dQ[k] = (FLOAT) 0.0;
      part.level  = 0;
      part.nstep  = 0;
      part.nlast  = 0;
      part.tlast  = t;
      part.flags.unset_flag(active);
    }
    for (i=0; i<mfv->Nhydro; i++) mfv->GetMeshlessFVParticlePointer(i).flags.set_flag(active);

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

    // ..
    for (i=0; i<mfv->Nhydro; i++) {
      MeshlessFVParticle<ndim>& part = mfv->GetMeshlessFVParticlePointer(i);
      for (k=0; k<ndim; k++) part.r0[k] = part.r[k];
      for (k=0; k<ndim; k++) part.v0[k] = part.v[k];
      for (k=0; k<ndim; k++) part.a0[k] = part.a[k];
      part.flags.set_flag(active);
    }

    mfvneib->UpdateAllProperties(mfv, nbody, simbox);


    LocalGhosts->CopyHydroDataToGhosts(simbox,mfv);
#ifdef MPI_PARALLEL
    MpiGhosts->CopyHydroDataToGhosts(simbox,mfv);
#endif


#ifdef MPI_PARALLEL
    //mfvneib->UpdateHydroExportList(rank, mfv->Nhydro, mfv->Ntot, partdata, mfv, nbody, simbox);

    //mpicontrol->ExportParticlesBeforeForceLoop(mfv);
    // Update pointer in case there has been a reallocation
    //partdata = mfv->GetMeshlessFVParticleArray();
#endif


    mfvneib->UpdateGradientMatrices(mfv, nbody, simbox);
    LocalGhosts->CopyHydroDataToGhosts(simbox,mfv);
#ifdef MPI_PARALLEL
    MpiGhosts->CopyHydroDataToGhosts(simbox,mfv);
#endif

    if (mfv->self_gravity == 1 || nbody->Nnbody > 0) {

      mfv->ZeroAccelerations();

#ifdef MPI_PARALLEL
      if (mfv->self_gravity ==1 ) {
        MeshlessFVParticle<ndim> *partdata = mfv->GetMeshlessFVParticleArray();
        for (int i=0; i< mfv->Nhydro; i++) {
          for (int k=0; k<ndim; k++) partdata[i].a[k] = 0.0;
          partdata[i].gpot=0.0;
        }
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
  } // End of force iteration.

  // Compute initial N-body forces
  //-----------------------------------------------------------------------------------------------
  if (mfv->self_gravity == 1 && mfv->Nhydro > 0) {
    mfvneib->UpdateAllStarGasForces(mfv, nbody);  //, simbox, ewald);

    // We need to sum up the contributions from the different domains
#if defined MPI_PARALLEL
    mpicontrol->ComputeTotalStarGasForces(nbody);
#endif
  }

  if (nbody->nbody_softening == 1) {
    nbody->CalculateDirectSmoothedGravForces(nbody->Nnbody, nbody->nbodydata);
  }
  else {
    nbody->CalculateDirectGravForces(nbody->Nnbody, nbody->nbodydata);
  }
  nbody->CalculateAllStartupQuantities(nbody->Nnbody, nbody->nbodydata);

  for (i=0; i<nbody->Nnbody; i++) {
    if (nbody->nbodydata[i]->active) {
      nbody->extpot->AddExternalPotential(nbody->nbodydata[i]->r, nbody->nbodydata[i]->v,
          nbody->nbodydata[i]->a, nbody->nbodydata[i]->adot,
          nbody->nbodydata[i]->gpot);
    }
  }


  // Call EndTimestep to set all 'beginning-of-step' variables
  mfv->EndTimestep(n, mfv->Nhydro, t, timestep, mfv->GetMeshlessFVParticleArray());
  nbody->EndTimestep(n, nbody->Nstar, t, timestep, nbody->nbodydata);

  this->CalculateDiagnostics();
  this->diag0 = this->diag;
  this->setup = true;

  rebuild_tree=false;

  return;
}



//=================================================================================================
//  MeshlessFVSimulation::ComputeGlobalTimestep
/// Computes global timestep for MFV simulation.  Calculates the minimum
/// timestep for all hydro and N-body particles in the simulation.
//=================================================================================================
template <int ndim>
void MeshlessFVSimulation<ndim>::ComputeGlobalTimestep(void)
{
  int i;                               // Particle counter
  DOUBLE dt_min = big_number_dp;       // Local copy of minimum timestep
  DOUBLE dt;


  debug2("[MeshlessFVSimulation::ComputeGlobalTimestep]");
  CodeTiming::BlockTimer timer = timing->StartNewTimer("GLOBAL_TIMESTEPS");


  // Only update timestep when all particles are synced at end of last step.
  //-----------------------------------------------------------------------------------------------
  if (n == nresync) {

    n            = 0;
    level_max    = 0;
    level_step   = level_max + integration_step - 1;
    nresync      = integration_step;
    dt_min_hydro = big_number_dp;
    dt_min_nbody = big_number_dp;


    // Find minimum timestep from all hydro particles
    //---------------------------------------------------------------------------------------------
#pragma omp parallel default(none) private(i,dt) shared(dt_min)
    {
      dt       = big_number_dp;           // Aux. minimum timestep
      DOUBLE dt_nbody = big_number_dp;           // Aux. minimum N-body timestep
      DOUBLE dt_hydro = big_number_dp;           // Aux. minimum hydro timestep

#pragma omp for
      for (i=0; i<mfv->Nhydro; i++) {
        MeshlessFVParticle<ndim>& part = mfv->GetMeshlessFVParticlePointer(i);
        part.level     = 0;
        part.levelneib = 0;
        part.nstep     = pow(2,level_step - part.level);
        part.nlast     = n;
        part.tlast     = t;
        part.dt        = mfv->Timestep(part);
        dt             = min(dt, part.dt);
        dt_hydro       = min(dt_hydro, part.dt);
      }

      // Now compute minimum timestep due to stars/systems
#pragma omp for
      for (i=0; i<nbody->Nnbody; i++) {
        nbody->nbodydata[i]->level = 0;
        nbody->nbodydata[i]->nstep = pow(2, level_step - nbody->nbodydata[i]->level);
        nbody->nbodydata[i]->nlast = n;
        nbody->nbodydata[i]->tlast = t;
        nbody->nbodydata[i]->dt    = nbody->Timestep(nbody->nbodydata[i]);
        dt       = min(dt,nbody->nbodydata[i]->dt);
        dt_nbody = min(dt_nbody, nbody->nbodydata[i]->dt);
      }

#pragma omp critical
      {
        if (dt < dt_min) dt_min = dt;
        if (dt_hydro < dt_min_hydro) dt_min_hydro = dt_hydro;
        if (dt_nbody < dt_min_nbody) dt_min_nbody = dt_nbody;
      }

    }
    //---------------------------------------------------------------------------------------------

#ifdef MPI_PARALLEL
    dt = dt_min;
    MPI_Allreduce(&dt, &dt_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
#endif
    timestep = dt_min;

    // Set minimum timestep for all hydro and N-body particles
    for (i=0; i<mfv->Nhydro; i++) mfv->GetMeshlessFVParticlePointer(i).dt = timestep;
    for (i=0; i<nbody->Nnbody; i++) nbody->nbodydata[i]->dt = timestep;

  }
  //-----------------------------------------------------------------------------------------------

  return;
}



//=================================================================================================
//  MeshlessFVSimulation::ComputeBlockTimesteps
/// Compute timesteps for all particles using hierarchical block timesteps.
//=================================================================================================
template <int ndim>
void MeshlessFVSimulation<ndim>::ComputeBlockTimesteps(void)
{
  int istep;                                 // Aux. variable for changing steps
  int last_level;                            // Previous timestep level
  int level;                                 // Particle timestep level
  int level_max_aux;                         // Aux. maximum level variable
  int level_max_nbody = 0;                   // level_max for star particles only
  int level_max_old = level_max;             // Old level_max
  int level_max_hydro = 0;                   // level_max for hydro particles only
  int level_min_hydro = 9999999;             // level_min for hydro particles
  int level_nbody;                           // local thread var. for N-body level
  int level_hydro;                           // local thread var. for hydro level
  int nfactor;                               // Increase/decrease factor of n
  //int nstep;                                 // Particle integer step-size
  DOUBLE dt;                                 // Aux. timestep variable
  DOUBLE dt_min = big_number_dp;             // Minimum timestep
  //DOUBLE dt_min_aux;                         // Aux. minimum timestep variable
  //DOUBLE dt_nbody;                           // Aux. minimum N-body timestep
  //DOUBLE dt_hydro;                           // Aux. minimum hydro timestep

  debug2("[MeshlessFVSimulation::ComputeBlockTimesteps]");
  CodeTiming::BlockTimer timer = timing->StartNewTimer("BLOCK_TIMESTEPS");


  dt_min_nbody = big_number_dp;
  dt_min_hydro = big_number_dp;


  // Synchronise all timesteps and reconstruct block timestep structure.
  //===============================================================================================
  if (n == nresync) {

    n = 0;
    timestep = big_number_dp;

#pragma omp parallel default(none) shared(dt_min) private(dt, level) //private(dt,dt_min_aux,dt_nbody,dt_hydro,i)
    {
      // Initialise all timestep and min/max variables
      DOUBLE dt_min_aux = big_number_dp;
      DOUBLE dt_hydro   = big_number_dp;
      DOUBLE dt_nbody   = big_number_dp;

      // Find minimum timestep from all hydro particles
#pragma omp for
      for (int i=0; i<mfv->Nhydro; i++) {
        MeshlessFVParticle<ndim>& part = mfv->GetMeshlessFVParticlePointer(i);
        if (part.flags.is_dead()) continue;
        dt         = mfv->Timestep(part);
        dt_min_aux = min(dt_min_aux, dt);
        dt_hydro   = min(dt_hydro, dt);
        part.dt    = dt;
      }

      // Now compute minimum timestep due to stars/systems
#pragma omp for
      for (int i=0; i<nbody->Nnbody; i++) {
        dt         = nbody->Timestep(nbody->nbodydata[i]);
        dt_min_aux = min(dt_min_aux, dt);
        dt_nbody   = min(dt_nbody, dt);
        nbody->nbodydata[i]->dt = dt;
      }

#pragma omp critical
      {
        timestep     = min(timestep, dt_min_aux);
        dt_min_hydro = min(dt_min_hydro, dt_hydro);
        dt_min_nbody = min(dt_min_nbody, dt_nbody);
      }
#pragma omp barrier
    }


    // For MPI, determine the global minimum timestep over all processors
#ifdef MPI_PARALLEL
    dt = timestep;
    MPI_Allreduce(&dt, &timestep, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    dt = dt_min_hydro;
    MPI_Allreduce(&dt, &dt_min_hydro, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    dt = dt_min_nbody;
    MPI_Allreduce(&dt, &dt_min_nbody, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
#endif


    // Calculate new block timestep levels
    level_max  = Nlevels - 1;
    level_step = level_max + integration_step - 1;
    dt_max     = timestep*powf(2.0, level_max);

    // Calculate the maximum level occupied by all hydro particles
    level_max_hydro = min(ComputeTimestepLevel(dt_min_hydro, dt_max), level_max);
    level_max_nbody = min(ComputeTimestepLevel(dt_min_nbody, dt_max), level_max);

    // Populate timestep levels with N-body particles.
    // Ensures that N-body particles occupy levels lower than all hydro particles
    for (int i=0; i<nbody->Nnbody; i++) {
      dt = nbody->nbodydata[i]->dt;
      level = min(ComputeTimestepLevel(dt, dt_max), level_max);
      nbody->nbodydata[i]->level = max(level, level_max_hydro);
      nbody->nbodydata[i]->nlast = n;
      nbody->nbodydata[i]->nstep = pow(2, level_step - nbody->nbodydata[i]->level);
      nbody->nbodydata[i]->tlast = t;
    }

    // If particles are sink neighbours, set to same timesteps as sinks
    for (int i=0; i<mfv->Nhydro; i++) {
      MeshlessFVParticle<ndim>& part = mfv->GetMeshlessFVParticlePointer(i);
      if (part.sinkid != -1) {
        if (sinks->sink[part.sinkid].star->level - part.level > level_diff_max) {
          part.level      = sinks->sink[part.sinkid].star->level - level_diff_max;
          part.levelneib  = sinks->sink[part.sinkid].star->level;
          level_max_hydro = max(level_max_hydro, part.level);
        }
      }
    }

    // If enforcing a single hydro timestep, set it here.
    // Otherwise, populate the timestep levels with hydro particles.
    if (sph_single_timestep == 1) {
      for (int i=0; i<mfv->Nhydro; i++) {
        MeshlessFVParticle<ndim>& part = mfv->GetMeshlessFVParticlePointer(i);
        if (part.flags.is_dead()) continue;
        part.flags.set_flag(active);
        part.level     = level_max_hydro;
        part.levelneib = level_max_hydro;
        part.nlast     = n;
        part.tlast     = t;
        part.nstep     = pow(2, level_step - part.level);
      }
      level_min_hydro = level_max_hydro;
    }
    else {
      for (int i=0; i<mfv->Nhydro; i++) {
        MeshlessFVParticle<ndim>& part = mfv->GetMeshlessFVParticlePointer(i);
        if (part.flags.is_dead()) continue;
        dt              = part.dt;
        level           = min((int) (invlogetwo*log(dt_max/dt)) + 1, level_max);
        level           = max(level, 0);
        part.flags.set_flag(active);
        part.level      = level;
        part.levelneib  = level;
        part.nlast      = n;
        part.tlast      = t;
        part.nstep      = pow(2, level_step - part.level);
        level_min_hydro = min(level_min_hydro, part.level);
      }
    }

    nresync = pow(2,level_step);
    assert(nresync > 0);
    timestep = dt_max / (DOUBLE) nresync;

  }
  // If not resynchronising, check if any hydro/N-body particles need to move
  // up or down timestep levels.
  //===============================================================================================
  else {

    level_max_old   = level_max;
    level_max       = 0;
    level_max_nbody = 0;
    level_max_hydro = 0;


#pragma omp parallel default(shared) private(dt, level) \
    shared(dt_min, istep, level_max_hydro, level_max_nbody, level_max_old, nfactor)
    {
      int last_level;
      int nstep;
      int level_max_aux = 0;
      int level_nbody   = 0;
      int level_hydro   = 0;
      DOUBLE dt;
      DOUBLE dt_hydro   = big_number_dp;
      DOUBLE dt_nbody   = big_number_dp;


      // Find all hydro particles at the beginning of a new timestep
      //-------------------------------------------------------------------------------------------
#pragma omp for
      for (int i=0; i<mfv->Nhydro; i++) {
        MeshlessFVParticle<ndim>& part = mfv->GetMeshlessFVParticlePointer(i);
        if (part.flags.is_dead()) continue;
        part.flags.unset_flag(active);

        if (part.nlast == n) {
          nstep      = part.nstep;
          last_level = part.level;

          // Compute new timestep value and level number
          dt    = mfv->Timestep(part);
          level = max(ComputeTimestepLevel(dt, dt_max), part.levelneib - level_diff_max);

          // Move up one level (if levels are correctly synchronised) or
          // down several levels if required
          if (level < last_level && last_level > 1 && n%(2*nstep) == 0) {
            part.level = last_level - 1;
          }
          else if (level > last_level) {
            part.level = level;
          }
          else {
            part.level = last_level;
          }

          part.flags.set_flag(active);
          part.levelneib = level;
          part.dt        = dt;
          part.nlast     = n;
          part.tlast     = t;
          part.nstep     = pow(2, level_step - part.level);
        }

        // Find maximum level of all hydro particles
        level_hydro   = max(level_hydro, part.level);
        level_max_aux = max(level_max_aux, part.level);

        dt_hydro = min(dt_hydro, part.dt);
      }
      //-------------------------------------------------------------------------------------------


#pragma omp critical
      {
        dt_min          = min(dt_min, dt_hydro);
        dt_min_hydro    = min(dt_min_hydro, dt_hydro);
        level_max       = max(level_max, level_max_aux);
        level_max_hydro = max(level_max_hydro, level_hydro);
      }
#pragma omp barrier

#if defined MPI_PARALLEL
#pragma omp master
      {
        level = level_max_hydro;
        MPI_Allreduce(&level, &level_max_hydro, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
      }
#pragma omp barrier
#endif

      // Now find all N-body particles at the beginning of a new timestep
      //-------------------------------------------------------------------------------------------
#pragma omp for
      for (int i=0; i<nbody->Nnbody; i++) {

        // Skip particles that are not at end of step
        if (nbody->nbodydata[i]->nlast == n) {
          nstep = nbody->nbodydata[i]->nstep;
          last_level = nbody->nbodydata[i]->level;

          // Compute new timestep value and level number
          dt    = nbody->Timestep(nbody->nbodydata[i]);
          level = max(ComputeTimestepLevel(dt, dt_max), level_max_hydro);

          // Move up one level (if levels are correctly synchronised) or
          // down several levels if required
          if (level < last_level && level > level_max_hydro && last_level > 1 && n%(2*nstep) == 0) {
            nbody->nbodydata[i]->level = last_level - 1;
          }
          else if (level > last_level) {
            nbody->nbodydata[i]->level = level;
          }
          else {
            nbody->nbodydata[i]->level = last_level;
          }

          nbody->nbodydata[i]->dt    = dt;
          nbody->nbodydata[i]->nlast = n;
          nbody->nbodydata[i]->nstep = pow(2, level_step - nbody->nbodydata[i]->level);
          nbody->nbodydata[i]->tlast = t;
        }

        // Find maximum level of all N-body particles
        level_nbody   = max(level_nbody, nbody->nbodydata[i]->level);
        level_max_aux = max(level_max_aux, nbody->nbodydata[i]->level);
        dt_nbody      = min(dt_nbody, nbody->nbodydata[i]->dt);
      }
      //-------------------------------------------------------------------------------------------

#pragma omp critical
      {
        dt_min          = min(dt_min, dt_nbody);
        dt_min_nbody    = min(dt_min_nbody, dt_nbody);
        level_max       = max(level_max, level_max_aux);
        level_max_nbody = max(level_max_nbody, level_nbody);
      }
#pragma omp barrier


      // Correct timestep levels for any particles that have entered a sink
      //-------------------------------------------------------------------------------------------
      /*#pragma omp for
      for (i=0; i<mfv->Nhydro; i++) {
        SphParticle<ndim>& part = mfv->GetSphParticlePointer(i);
        if (part.flags.is_dead() || part.nlast != n) continue;
        if (part.sinkid != -1) {
          if (sinks.sink[part.sinkid].star->level - part.level > level_diff_max) {
            part.level = sinks.sink[part.sinkid].star->level - level_diff_max;
          }
        }
      }*/

    }

    // For MPI, find the global maximum timestep levels for each processor
#ifdef MPI_PARALLEL
    level = level_max;
    MPI_Allreduce(&level, &level_max, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    level = level_max_nbody;
    MPI_Allreduce(&level, &level_max_nbody, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    assert(level_max_hydro >= 0);
#endif

    assert(!(isnan(dt_min)) && !(isinf(dt_min)));
    assert(!(isnan(dt_max)) && !(isinf(dt_max)));

    // Set fixed hydro timestep level here in case maximum has changed
    if (sph_single_timestep == 1) {
      for (int i=0; i<mfv->Nhydro; i++) {
        MeshlessFVParticle<ndim>& part = mfv->GetMeshlessFVParticlePointer(i);
        if (part.flags.is_dead()) continue;
        if (part.nlast == n) part.level = level_max_hydro;
      }
    }


    istep = pow(2, level_step - level_max_old + 1);

    // Adjust integer time if levels are added or removed
    //---------------------------------------------------------------------------------------------
    if (level_max > level_max_old) {
      nfactor = pow(2, level_max - level_max_old);
      n *= nfactor;
      for (int i=0; i<mfv->Nhydro; i++) {
        MeshlessFVParticle<ndim>& part = mfv->GetMeshlessFVParticlePointer(i);
        if (part.flags.is_dead()) continue;
        part.nstep *= nfactor;
        part.nlast *= nfactor;
      }
      for (int i=0; i<nbody->Nnbody; i++) nbody->nbodydata[i]->nstep *= nfactor;
      for (int i=0; i<nbody->Nnbody; i++) nbody->nbodydata[i]->nlast *= nfactor;
    }
    else if (level_max <= level_max_old - 1 && level_max_old > 1 && n%istep == 0) {
      level_max = level_max_old - 1;

      nfactor = pow(2, level_max_old - level_max);
      assert(n%nfactor == 0);
      n /= nfactor;
      for (int i=0; i<mfv->Nhydro; i++) {
        MeshlessFVParticle<ndim>& part = mfv->GetMeshlessFVParticlePointer(i);
        if (part.flags.is_dead()) continue;
        part.nlast /= nfactor;
        part.nstep /= nfactor;
      }
      for (int i=0; i<nbody->Nnbody; i++) nbody->nbodydata[i]->nlast /= nfactor;
      for (int i=0; i<nbody->Nnbody; i++) nbody->nbodydata[i]->nstep /= nfactor;
    }
    else {
      level_max = level_max_old;
    }
    //---------------------------------------------------------------------------------------------


    level_step = level_max + integration_step - 1;
    nresync    = pow(2,level_step);
    timestep   = dt_max / (DOUBLE) nresync;

    // Update values of nstep for both hydro and star particles
    for (int i=0; i<mfv->Nhydro; i++) {
      MeshlessFVParticle<ndim>& part = mfv->GetMeshlessFVParticlePointer(i);
      if (part.flags.is_dead()) continue;
      if (part.nlast == n) part.nstep = pow(2,level_step - part.level);
    }
    for (int i=0; i<nbody->Nnbody; i++) {
      if (nbody->nbodydata[i]->nlast == n) {
        nbody->nbodydata[i]->nstep = pow(2, level_step - nbody->nbodydata[i]->level);
      }
    }

    assert(level_max >= level_max_old - 1);

  }
  //===============================================================================================


  // Various asserts for debugging
  assert(timestep >= 0.0);
  assert(level_step == level_max + integration_step - 1);
  assert(level_max_hydro <= level_max);
  assert(level_max_nbody <= level_max);
  assert(n <= nresync);
  for (int i=0; i<mfv->Nhydro; i++) {
    MeshlessFVParticle<ndim>& part = mfv->GetMeshlessFVParticlePointer(i);
    if (part.flags.is_dead()) continue;
    assert(part.level <= level_max);
    assert(part.nlast <= n);
    assert(part.tlast <= t);
    assert(part.nstep == pow(2,level_step - part.level));
    assert(pow(2,level_step - part.level) >= n - part.nlast);
    assert(part.nlast != n || n%part.nstep == 0);
  }
  for (int i=0; i<nbody->Nnbody; i++) {
    assert(nbody->nbodydata[i]->level <= level_max);
    assert(nbody->nbodydata[i]->nlast <= n);
    assert(nbody->nbodydata[i]->nstep == pow(2,level_step - nbody->nbodydata[i]->level));
    assert(pow(2,level_step - nbody->nbodydata[i]->level) >= n - nbody->nbodydata[i]->nlast);
    assert(nbody->nbodydata[i]->nlast != n || n%nbody->nbodydata[i]->nstep == 0);
    assert(nbody->nbodydata[i]->level >= level_max_hydro);
    assert(nbody->nbodydata[i]->tlast <= t);
  }
  if (timestep <= 0.0) {
    cout << "Timestep fallen to zero : " << timestep << endl;
    ExceptionHandler::getIstance().raise("Timestep fallen to zero");
  }

  return;


  // Some validations
  //-----------------------------------------------------------------------------------------------
  /*int *ninlevel;
  int Nactive=0;
  ninlevel = new int[level_max+1];
  MeshlessFVParticle<ndim>& part = mfv->GetMeshlessFVParticlePointer(imin);
  cout << "-----------------------------------------------------" << endl;
  cout << "Checking timesteps : " << level_max << "   " << level_max_hydro << "    "
       << level_max_nbody << "    " << level_step << "   " << level_max_old << endl;
  cout << "n : " << n << endl;
  cout << "dt_min_hydro : " << dt_min_hydro << "    dt_min_nbody : " << dt_min_nbody
       << "    timestep : " << timestep << endl;
  cout << "imin : " << imin << "    " << part.dt << "     " << part.h << "    " << "    "
       << part.sound << "     " << part.div_v << "     "
       << part.h/(part.sound + part.h*fabs(part.div_v)) << endl;
  for (int l=0; l<=level_max; l++) ninlevel[l] = 0;
  for (i=0; i<mfv->Nhydro; i++) {
    MeshlessFVParticle<ndim>& part = mfv->GetMeshlessFVParticlePointer(i);
    if (part.flags.check_flag(active)) Nactive++;
    ninlevel[part.level]++;
  }
  cout << "No. of active Hydro particles : " << Nactive << endl;
  cout << "Hydro level occupancy" << endl;
  for (int l=0; l<=level_max; l++) cout << "level : " << l << "     N : " << ninlevel[l] << endl;
  for (int l=0; l<=level_max; l++) ninlevel[l] = 0;
  for (i=0; i<nbody->Nstar; i++) ninlevel[nbody->nbodydata[i]->level]++;
  cout << "N-body level occupancy" << endl;
  for (int l=0; l<=level_max; l++) cout << "level : " << l << "     N : " << ninlevel[l] << endl;

  delete[] ninlevel;

  if (timestep <= 0.0) {
    cout << "Timestep fallen to zero : " << timestep << endl;
    ExceptionHandler::getIstance().raise("Timestep fallen to zero");
  }

  return;*/
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
