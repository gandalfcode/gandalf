//=================================================================================================
//  SphSimulation.cpp
//  Contains all main functions controlling SPH simulation work-flow.
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
#include "Sph.h"
#include "RiemannSolver.h"
#include "Ghosts.h"
#include "Sinks.h"
using namespace std;



// Create template class instances of the main SphSimulation object for
// each dimension used (1, 2 and 3)
template class SphSimulation<1>;
template class SphSimulation<2>;
template class SphSimulation<3>;



//=================================================================================================
//  SphSimulation::ProcessParameters
/// Process all the options chosen in the parameters file, setting various
/// simulation variables and creating important simulation objects.
//=================================================================================================
template <int ndim>
void SphSimulation<ndim>::ProcessParameters(void)
{
  // Local references to parameter variables for brevity
  map<string, int> &intparams = simparams->intparams;
  map<string, double> &floatparams = simparams->floatparams;
  map<string, string> &stringparams = simparams->stringparams;
  string sim = stringparams["sim"];
  string gas_eos = stringparams["gas_eos"];
  string gas_radiation = stringparams["radiation"];

  debug2("[SphSimulation::ProcessParameters]");


  // Sanity check for valid dimensionality
  if (ndim < 1 || ndim > 3) {
    string message = "Invalid dimensionality chosen : ndim = " + ndim;
    ExceptionHandler::getIstance().raise(message);
  }

  // Set-up random number generator object
  //-----------------------------------------------------------------------------------------------
  if (stringparams["rand_algorithm"] == "xorshift")
    randnumb = new XorshiftRand(intparams["randseed"]);
  else if (stringparams["rand_algorithm"] == "none")
    randnumb = new DefaultSystemRand(intparams["randseed"]);
  else {
    string message = "Unrecognised parameter : rand_algorithm= " +
      stringparams["rand_algorithm"];
    ExceptionHandler::getIstance().raise(message);
  }


  // Set-up all output units for scaling parameters
  simunits.SetupUnits(simparams);

  // Boundary condition variables
  //-----------------------------------------------------------------------------------------------
  simbox.x_boundary_lhs = setBoundaryType(stringparams["x_boundary_lhs"]);
  simbox.x_boundary_rhs = setBoundaryType(stringparams["x_boundary_rhs"]);
  simbox.boxmin[0] = floatparams["boxmin[0]"]/simunits.r.outscale;
  simbox.boxmax[0] = floatparams["boxmax[0]"]/simunits.r.outscale;

  if (ndim > 1) {
    simbox.y_boundary_lhs = setBoundaryType(stringparams["y_boundary_lhs"]);
    simbox.y_boundary_rhs = setBoundaryType(stringparams["y_boundary_rhs"]);
    simbox.boxmin[1] = floatparams["boxmin[1]"]/simunits.r.outscale;
    simbox.boxmax[1] = floatparams["boxmax[1]"]/simunits.r.outscale;
  }

  if (ndim == 3) {
    simbox.z_boundary_lhs = setBoundaryType(stringparams["z_boundary_lhs"]);
    simbox.z_boundary_rhs = setBoundaryType(stringparams["z_boundary_rhs"]);
    simbox.boxmin[2] = floatparams["boxmin[2]"]/simunits.r.outscale;
    simbox.boxmax[2] = floatparams["boxmax[2]"]/simunits.r.outscale;
  }

  for (int k=0; k<ndim; k++) {
    simbox.boxsize[k] = simbox.boxmax[k] - simbox.boxmin[k];
    simbox.boxhalf[k] = 0.5*simbox.boxsize[k];
  }


  // Set-up main SPH objects depending on which SPH algorithm we are using
  ProcessSphParameters();
  hydro = sph;

  // Process all N-body parameters and set-up main N-body objects
  this->ProcessNbodyParameters();


  // Thermal physics object.  If energy equation is chosen, also initiate
  // the energy integration object.
  //-----------------------------------------------------------------------------------------------
  if ((gas_eos == "energy_eqn" || gas_eos == "constant_temp" ||
       gas_eos == "isothermal" || gas_eos == "barotropic" ||
       gas_eos == "barotropic2") && gas_radiation == "ionisation") {
    sph->eos = new IonisingRadiation<ndim>
      (gas_eos, floatparams["temp0"], floatparams["mu_bar"],
      floatparams["gamma_eos"], floatparams["rho_bary"], &simunits, sphneib);
  }
  else if ((gas_eos == "energy_eqn" || gas_eos == "constant_temp" ||
             gas_eos == "isothermal" || gas_eos == "barotropic" ||
             gas_eos == "barotropic2") && gas_radiation == "monoionisation") {
    sph->eos = new MCRadiationEOS<ndim>
      (gas_eos, floatparams["temp0"], floatparams["temp_ion"], floatparams["mu_bar"],
       floatparams["mu_ion"], floatparams["gamma_eos"], floatparams["rho_bary"], &simunits);
  }
  else if (gas_eos == "energy_eqn" || gas_eos == "constant_temp") {
    sph->eos = new Adiabatic<ndim>
      (floatparams["temp0"], floatparams["mu_bar"], floatparams["gamma_eos"]);
  }
  else if (gas_eos == "isothermal") {
    sph->eos = new Isothermal<ndim>
      (floatparams["temp0"], floatparams["mu_bar"], floatparams["gamma_eos"], &simunits);
  }
  else if (gas_eos == "barotropic") {
    sph->eos = new Barotropic<ndim>(floatparams["temp0"], floatparams["mu_bar"],
                                    floatparams["gamma_eos"], floatparams["rho_bary"], &simunits);
  }
  else if (gas_eos == "barotropic2") {
    sph->eos = new Barotropic2<ndim>(floatparams["temp0"], floatparams["mu_bar"],
                                     floatparams["gamma_eos"], floatparams["rho_bary"], &simunits);
  }
  else {
    string message = "Unrecognised parameter : gas_eos = " + gas_eos;
    ExceptionHandler::getIstance().raise(message);
  }


  // Set external potential field object and set pointers to object
  if (stringparams["external_potential"] == "none") {
    extpot = new NullPotential<ndim>();
  }
  else if (stringparams["external_potential"] == "plummer") {
    extpot = new PlummerPotential<ndim>(floatparams["mplummer"], floatparams["rplummer"]);
  }
  else {
    string message = "Unrecognised parameter : external_potential = "
      + simparams->stringparams["external_potential"];
    ExceptionHandler::getIstance().raise(message);
  }
  sph->extpot = extpot;
  nbody->extpot = extpot;


  // Create Ewald periodic gravity object
  periodicBoundaries = IsAnyBoundarySpecial(simbox);
  if (periodicBoundaries && intparams["self_gravity"] == 1) {
    ewaldGravity = true;
    ewald = new Ewald<ndim>(simbox,intparams["gr_bhewaldseriesn"],intparams["in"],
                            intparams["nEwaldGrid"],floatparams["ewald_mult"],
                            floatparams["ixmin"],floatparams["ixmax"],timing);
  }


  // Set all other SPH parameter variables
  //sph->Nhydro            = intparams["Nhydro"];
  sph->Nhydromax         = intparams["Nhydromax"];
  sph->create_sinks    = intparams["create_sinks"];
  sph->fixed_sink_mass = intparams["fixed_sink_mass"];
  sph->msink_fixed     = floatparams["m1"];
  sph->alpha_visc_min  = floatparams["alpha_visc_min"];


  // Set important variables for N-body objects
  //nbody->Nstar          = intparams["Nstar"];
  nbody->Nstarmax       = intparams["Nstarmax"];
  nbody_single_timestep = intparams["nbody_single_timestep"];
  nbodytree.gpehard     = floatparams["gpehard"];
  nbodytree.gpesoft     = floatparams["gpesoft"];
  nbody->perturbers     = intparams["perturbers"];
  if (intparams["sub_systems"] == 1) subsystem->perturbers = intparams["perturbers"];


  // Sink particles
  //-----------------------------------------------------------------------------------------------
  sink_particles            = intparams["sink_particles"];
  sinks.sink_particles      = intparams["sink_particles"];
  sinks.create_sinks        = intparams["create_sinks"];
  sinks.smooth_accretion    = intparams["smooth_accretion"];
  sinks.alpha_ss            = floatparams["alpha_ss"];
  sinks.smooth_accrete_frac = floatparams["smooth_accrete_frac"];
  sinks.smooth_accrete_dt   = floatparams["smooth_accrete_dt"];
  sinks.sink_radius_mode    = stringparams["sink_radius_mode"];
  sinks.rho_sink            = floatparams["rho_sink"];
  sinks.rho_sink            /= simunits.rho.outscale/simunits.rho.outcgs;

  if (sinks.sink_radius_mode == "fixed") {
    sinks.sink_radius = floatparams["sink_radius"]/simunits.r.outscale;
  }
  else {
    sinks.sink_radius = floatparams["sink_radius"];
  }

  // Sanity-check for various sink particle values
  if (intparams["sink_particles"] == 1 &&
      (stringparams["nbody"] != "lfkdk" && stringparams["nbody"] != "lfdkd")) {
    string message = "Invalid parameter : nbody must use lfkdk or lfdkd when "
      "using accreting sink particles";
    ExceptionHandler::getIstance().raise(message);
  }


  // Set other important simulation variables
  dt_litesnap         = floatparams["dt_litesnap"]/simunits.t.outscale;
  dt_python           = floatparams["dt_python"];
  dt_snap             = floatparams["dt_snap"]/simunits.t.outscale;
  extra_sink_output   = intparams["extra_sink_output"];
  level_diff_max      = intparams["level_diff_max"];
  litesnap            = intparams["litesnap"];
  Nlevels             = intparams["Nlevels"];
  ndiagstep           = intparams["ndiagstep"];
  noutputstep         = intparams["noutputstep"];
  nrestartstep        = intparams["nrestartstep"];
  ntreebuildstep      = intparams["ntreebuildstep"];
  ntreestockstep      = intparams["ntreestockstep"];
  Nstepsmax           = intparams["Nstepsmax"];
  out_file_form       = stringparams["out_file_form"];
  pruning_level       = intparams["pruning_level"];
  run_id              = stringparams["run_id"];
  sph_single_timestep = intparams["sph_single_timestep"];
  tmax_wallclock      = floatparams["tmax_wallclock"];
  tend                = floatparams["tend"]/simunits.t.outscale;
  tlitesnapnext       = floatparams["tlitesnapfirst"]/simunits.t.outscale;
  tsnapnext           = floatparams["tsnapfirst"]/simunits.t.outscale;


  // Set pointers to timing object
  nbody->timing   = timing;
  if (sim == "sph" || sim == "gradhsph" || sim == "sm2012sph" || sim == "godunov_sph") {
    sinks.timing    = timing;
    sphint->timing  = timing;
    sphneib->timing = timing;
    uint->timing    = timing;
    radiation->timing = timing;
  }


  // Flag that we've processed all parameters already
  ParametersProcessed = true;


  return;
}



//=================================================================================================
//  SphSimulation::PostInitialConditionsSetup
/// Call routines for calculating all initial SPH and N-body quantities
/// once initial conditions have been set-up.
//=================================================================================================
template <int ndim>
void SphSimulation<ndim>::PostInitialConditionsSetup(void)
{
  int i;                               // Particle counter
  int k;                               // Dimension counter
  SphParticle<ndim> *partdata;         // Pointer to main SPH data array

  debug2("[SphSimulation::PostInitialConditionsSetup]");

  // Set iorig
  if (rank == 0) {
    for (i=0; i<sph->Nhydro; i++) sph->GetSphParticlePointer(i).iorig = i;
  }


  // Perform initial MPI decomposition
  //-----------------------------------------------------------------------------------------------
#ifdef MPI_PARALLEL
  mpicontrol->CreateInitialDomainDecomposition(sph,nbody,simparams,simbox,this->initial_h_provided);
#endif

  // Set pointer to SPH particle data
  partdata = sph->GetSphParticleArray();

  // Set time variables here (for now)
  nresync = 0;
  n = 0;

  // Set initial smoothing lengths and create initial ghost particles
  //-----------------------------------------------------------------------------------------------
  if (sph->Nhydro > 0) {

    // Set all relevant particle counters
    sph->Nghost = 0;
    sph->Nghostmax = sph->Nhydromax - sph->Nhydro;
    sph->Ntot = sph->Nhydro;
    for (i=0; i<sph->Nhydro; i++) sph->GetSphParticlePointer(i).active = true;

    // Set initial artificial viscosity alpha values
    if (simparams->stringparams["time_dependent_avisc"] == "none") {
      for (i=0; i<sph->Nhydro; i++) sph->GetSphParticlePointer(i).alpha = sph->alpha_visc;
    }
    else {
      for (i=0; i<sph->Nhydro; i++) sph->GetSphParticlePointer(i).alpha = sph->alpha_visc_min;
    }

    // Compute mean mass (used for smooth sink accretion)
    if (!restart) {
      sph->mmean = 0.0;
      for (i=0; i<sph->Nhydro; i++) sph->mmean += sph->GetSphParticlePointer(i).m;
      sph->mmean /= (FLOAT) sph->Nhydro;
    }

    // Compute minimum smoothing length of sinks
    sph->hmin_sink = big_number;
    for (i=0; i<sinks.Nsink; i++) sph->hmin_sink = min(sph->hmin_sink,sinks.sink[i].star->h);

    // If the smoothing lengths have not been provided beforehand, then
    // calculate the initial values here
    sphneib->neibcheck = false;
    if (!this->initial_h_provided) {
      sph->InitialSmoothingLengthGuess();
      sphneib->BuildTree(rebuild_tree,0,ntreebuildstep,ntreestockstep,sph->Ntot,
                         sph->Nhydromax,sph->GetSphParticleArray(),sph,timestep);
      sphneib->UpdateAllSphProperties(sph->Nhydro,sph->Ntot,sph->GetSphParticleArray(),sph,nbody);
    }
    else {
      sphneib->BuildTree(rebuild_tree,0,ntreebuildstep,ntreestockstep,sph->Ntot,
                         sph->Nhydromax,sph->GetSphParticleArray(),sph,timestep);
    }

#ifdef MPI_PARALLEL
    mpicontrol->UpdateAllBoundingBoxes(sph->Nhydro, sph, sph->kernp);
#endif

    // Search ghost particles
    sphneib->SearchBoundaryGhostParticles(0.0,simbox,sph);
    sphneib->BuildGhostTree(true,0,ntreebuildstep,ntreestockstep,sph->Ntot,
                            sph->Nhydromax,sph->GetSphParticleArray(),sph,timestep);
#ifdef MPI_PARALLEL
    mpicontrol->UpdateAllBoundingBoxes(sph->Nhydro+sph->NPeriodicGhost, sph, sph->kernp);
    for (int i=0; i<sph->Nhydro+sph->NPeriodicGhost; i++) {
      SphParticle<ndim>& parti = sph->GetSphParticlePointer(i);
      parti.hrangesqd = sph->kernfacsqd*sph->kernp->kernrangesqd*parti.h*parti.h;
    }
    MpiGhosts->SearchGhostParticles(0.0,simbox,sph);
    sphneib->BuildMpiGhostTree(true,0,ntreebuildstep,ntreestockstep,sph->Ntot,
                               sph->Nhydromax,sph->GetSphParticleArray(),sph,timestep);
#endif


    // Zero accelerations
    for (i=0; i<sph->Nhydro; i++) sph->GetSphParticlePointer(i).active = true;

    // Update neighbour tree
    rebuild_tree = true;
    sphneib->BuildTree(rebuild_tree,0,ntreebuildstep,ntreestockstep,sph->Ntot,
                       sph->Nhydromax,sph->GetSphParticleArray(),sph,timestep);
    level_step = 1;


    // For Eigenvalue MAC, need non-zero values
    for (i=0; i<sph->Nhydro; i++) sph->GetSphParticlePointer(i).gpot = big_number;

    // Calculate all SPH properties
    sphneib->UpdateAllSphProperties(sph->Nhydro,sph->Ntot,sph->GetSphParticleArray(),sph,nbody);

#ifdef MPI_PARALLEL
    mpicontrol->UpdateAllBoundingBoxes(sph->Nhydro, sph, sph->kernp);
#endif

    // Search ghost particles
    sphneib->SearchBoundaryGhostParticles(0.0,simbox,sph);
    sphneib->BuildGhostTree(true,0,ntreebuildstep,ntreestockstep,sph->Ntot,
                            sph->Nhydromax,sph->GetSphParticleArray(),sph,timestep);
#ifdef MPI_PARALLEL
    mpicontrol->UpdateAllBoundingBoxes(sph->Nhydro+sph->NPeriodicGhost, sph, sph->kernp);
    MpiGhosts->SearchGhostParticles(0.0,simbox,sph);
    sphneib->BuildMpiGhostTree(true,0,ntreebuildstep,ntreestockstep,sph->Ntot,
                               sph->Nhydromax,sph->GetSphParticleArray(),sph,timestep);
#endif

    // Update neighbour tree
    rebuild_tree = true;
    sphneib->BuildTree(rebuild_tree,0,ntreebuildstep,ntreestockstep,
                       sph->Ntot,sph->Nhydromax,sph->GetSphParticleArray(),sph,timestep);
    sphneib->SearchBoundaryGhostParticles(0.0,simbox,sph);
    sphneib->BuildGhostTree(true,0,ntreebuildstep,ntreestockstep,sph->Ntot,
                            sph->Nhydromax,sph->GetSphParticleArray(),sph,timestep);
    sphneib->neibcheck = true;

    // Communicate pruned trees for MPI
#ifdef MPI_PARALLEL
    sphneib->BuildPrunedTree(pruning_level,rank);
    mpicontrol->CommunicatePrunedTrees();
//    exit(0);
#endif

  }


  // Compute all initial N-body terms
  //-----------------------------------------------------------------------------------------------
  if (nbody->Nstar > 0) {

    // Zero all acceleration terms
    for (i=0; i<nbody->Nstar; i++) {
      for (k=0; k<ndim; k++) nbody->stardata[i].a[k] = 0.0;
      for (k=0; k<ndim; k++) nbody->stardata[i].adot[k] = 0.0;
      for (k=0; k<ndim; k++) nbody->stardata[i].a2dot[k] = 0.0;
      for (k=0; k<ndim; k++) nbody->stardata[i].a3dot[k] = 0.0;
      nbody->stardata[i].gpot   = 0.0;
      nbody->stardata[i].gpe    = 0.0;
      nbody->stardata[i].tlast  = t;
      nbody->stardata[i].active = true;
      nbody->stardata[i].level  = 0;
      nbody->stardata[i].nstep  = 0;
      nbody->stardata[i].nlast  = 0;
      nbody->nbodydata[i]       = &(nbody->stardata[i]);
    }
    nbody->Nnbody = nbody->Nstar;

  }


  // Read-in N-body table here
//  nbody->LoadStellarPropertiesTable(&simunits);
//  nbody->UpdateStellarProperties();


  // Compute all initial SPH force terms
  //-----------------------------------------------------------------------------------------------
  if (sph->Nhydro > 0) {

    // Zero accelerations (here for now)
    for (i=0; i<sph->Ntot; i++) {
      SphParticle<ndim>& part = sph->GetSphParticlePointer(i);
      part.level  = 0;
      part.nstep  = 0;
      part.nlast  = 0;
      part.tlast  = t;
      part.active = false;
    }
    for (i=0; i<sph->Nhydro; i++) sph->GetSphParticlePointer(i).active = true;

    // Copy all other data from real SPH particles to ghosts
    LocalGhosts->CopySphDataToGhosts(simbox,sph);

    sphneib->BuildTree(true,0,ntreebuildstep,ntreestockstep,sph->Ntot,
                       sph->Nhydromax,sph->GetSphParticleArray(),sph,timestep);
    sphneib->SearchBoundaryGhostParticles(0.0,simbox,sph);
    sphneib->BuildGhostTree(true,0,ntreebuildstep,ntreestockstep,sph->Ntot,
                            sph->Nhydromax,sph->GetSphParticleArray(),sph,timestep);

#ifdef MPI_PARALLEL
    if (sph->self_gravity == 1)
      sphneib->UpdateGravityExportList(rank,sph->Nhydro,sph->Ntot,
                                       sph->GetSphParticleArray(),sph,nbody);
    else
      sphneib->UpdateHydroExportList(rank,sph->Nhydro,sph->Ntot,
                                     sph->GetSphParticleArray(),sph,nbody);

    mpicontrol->ExportParticlesBeforeForceLoop(sph);
#endif


    // Update the radiation field
    for (int jj=0; jj<10; jj++) {
      radiation->UpdateRadiationField(sph->Nhydro, nbody->Nnbody, sinks.Nsink,
                                      partdata, nbody->nbodydata, sinks.sink);
    }


    // Update thermal properties (if radiation field has altered them)
    for (i=0; i<sph->Nhydro;i++) {
      SphParticle<ndim>& part = sph->GetSphParticlePointer(i);
      sph->ComputeThermalProperties(part);
    }


    // Calculate SPH gravity and hydro forces, depending on which are activated
    if (ewaldGravity && sph->hydro_forces == 1 && sph->self_gravity == 1)
      sphneib->UpdateAllSphPeriodicForces(sph->Nhydro,sph->Ntot,partdata,sph,nbody,simbox,ewald);
    else if (ewaldGravity && sph->self_gravity == 1)
      sphneib->UpdateAllSphPeriodicGravForces(sph->Nhydro,sph->Ntot,partdata,sph,nbody,simbox,ewald);
    else if (sph->hydro_forces == 1 && sph->self_gravity == 1)
      sphneib->UpdateAllSphForces(sph->Nhydro,sph->Ntot,partdata,sph,nbody);
    else if (sph->hydro_forces == 1)
      sphneib->UpdateAllSphHydroForces(sph->Nhydro,sph->Ntot,sph->GetSphParticleArray(),sph,nbody);
    else if (sph->self_gravity == 1)
      sphneib->UpdateAllSphGravForces(sph->Nhydro,sph->Ntot,sph->GetSphParticleArray(),sph,nbody);

#if defined MPI_PARALLEL
    mpicontrol->GetExportedParticlesAccelerations(sph);
#endif

    // Set initial accelerations
    for (i=0; i<sph->Nhydro; i++) {
      SphParticle<ndim>& part = sph->GetSphParticlePointer(i);
      for (k=0; k<ndim; k++) part.r0[k] = part.r[k];
      for (k=0; k<ndim; k++) part.v0[k] = part.v[k];
      for (k=0; k<ndim; k++) part.a0[k] = part.a[k];
      part.active = false;
    }

    LocalGhosts->CopySphDataToGhosts(simbox,sph);
#ifdef MPI_PARALLEL
//    MpiGhosts->CopySphDataToGhosts(simbox,sph);
#endif

  }


  // Compute initial N-body forces
  //-----------------------------------------------------------------------------------------------
  if (nbody->Nstar > 0) {

    if (nbody->nbody_softening == 1) {
      nbody->CalculateDirectSmoothedGravForces(nbody->Nnbody,nbody->nbodydata);
    }
    else {
      nbody->CalculateDirectGravForces(nbody->Nnbody,nbody->nbodydata);
    }
    if (sph->self_gravity == 1 && sph->Nhydro > 0)
      sphneib->UpdateAllStarGasForces(sph->Nhydro,sph->Ntot,sph->GetSphParticleArray(),sph,nbody);
    nbody->CalculateAllStartupQuantities(nbody->Nnbody,nbody->nbodydata);

  }


  // Set particle values for initial step (e.g. r0, v0, a0, u0)
  uint->EndTimestep(n,sph->Nhydro,t,timestep,sph->GetSphParticleArray());
  sphint->EndTimestep(n,sph->Nhydro,t,timestep,sph->GetSphParticleArray());
  nbody->EndTimestep(n,nbody->Nstar,t,timestep,nbody->nbodydata);

  this->CalculateDiagnostics();
  this->diag0 = this->diag;
  this->setup = true;


  return;
}



//=================================================================================================
//  SphSimulation::MainLoop
/// Main SPH simulation integration loop.
//=================================================================================================
template <int ndim>
void SphSimulation<ndim>::MainLoop(void)
{
  int activecount = 0;                 // Flag if we need to recompute particles
  int i;                               // Particle loop counter
  int it;                              // Time-symmetric iteration counter
  int k;                               // Dimension counter
  FLOAT tghost;                        // Approx. ghost particle lifetime
  SphParticle<ndim> *partdata;         // Pointer to main SPH data array

  debug2("[SphSimulation::MainLoop]");


  // Set pointer for SPH data array
  partdata = sph->GetSphParticleArray();

  // Compute timesteps for all particles
  if (Nlevels == 1) this->ComputeGlobalTimestep();
  else this->ComputeBlockTimesteps();

  // Advance time variables
  n = n + 1;
  Nsteps = Nsteps + 1;
  t = t + timestep;
  if (n == nresync) Nblocksteps = Nblocksteps + 1;
  if (n%integration_step == 0) Nfullsteps = Nfullsteps + 1;

  // Advance SPH and N-body particles' positions and velocities
  uint->EnergyIntegration(n,sph->Nhydro,(FLOAT) t,(FLOAT) timestep,sph->GetSphParticleArray());
  sphint->AdvanceParticles(n,sph->Nhydro,(FLOAT) t,(FLOAT) timestep,sph->GetSphParticleArray());
  nbody->AdvanceParticles(n,nbody->Nnbody,t,timestep,nbody->nbodydata);

  // Check all boundary conditions
  // (DAVID : Move this function to sphint and create an analagous one
  //  for N-body.  Also, only check this on tree-build steps)
  if (Nsteps%ntreebuildstep == 0 || rebuild_tree) sphint->CheckBoundaries(simbox,sph);


  //-----------------------------------------------------------------------------------------------
  // MPI : On tree re-build step, determine load balancing for all MPI nodes.
  //       (How is this done?  All computed on root node??)
  //       Send/receive particles to their new nodes.
  //       Compute and transmit all bounding boxes (e.g. all particles, active
  //       particles, h-extent, ghosts, etc..) to all other MPI nodes
  //-----------------------------------------------------------------------------------------------
#ifdef MPI_PARALLEL
  if (Nsteps%ntreebuildstep == 0 || rebuild_tree) {
    sphneib->BuildPrunedTree(pruning_level,rank);
    mpicontrol->UpdateAllBoundingBoxes(sph->Nhydro, sph, sph->kernp);
    mpicontrol->CommunicatePrunedTrees();
    mpicontrol->LoadBalancing(sph,nbody);
  }
#endif


  // Compute all SPH quantities
  //-----------------------------------------------------------------------------------------------
  if (sph->Nhydro > 0) {

    // Rebuild or update local neighbour and gravity tree
    sphneib->BuildTree(rebuild_tree,Nsteps,ntreebuildstep,ntreestockstep,
                       sph->Ntot,sph->Nhydromax,partdata,sph,timestep);


    // Search for new ghost particles and create on local processor
    if (Nsteps%ntreebuildstep == 0 || rebuild_tree) {
      tghost = timestep*(FLOAT)(ntreebuildstep - 1);
      sphneib->SearchBoundaryGhostParticles(tghost,simbox,sph);
      sphneib->BuildGhostTree(rebuild_tree,Nsteps,ntreebuildstep,ntreestockstep,
                              sph->Ntot,sph->Nhydromax,partdata,sph,timestep);
#ifdef MPI_PARALLEL
      sphneib->BuildPrunedTree(pruning_level,rank);
      mpicontrol->CommunicatePrunedTrees();
      mpicontrol->UpdateAllBoundingBoxes(sph->Nhydro+sph->NPeriodicGhost,sph,sph->kernp);
      MpiGhosts->SearchGhostParticles(tghost,simbox,sph);
      sphneib->BuildMpiGhostTree(rebuild_tree,Nsteps,ntreebuildstep,ntreestockstep,
                                 sph->Ntot,sph->Nhydromax,partdata,sph,timestep);
#endif
    }
    // Otherwise copy properties from original particles to ghost particles
    else {
      LocalGhosts->CopySphDataToGhosts(simbox,sph);
#ifdef MPI_PARALLEL
      MpiGhosts->CopySphDataToGhosts(simbox,sph);
#endif
    }


    // Reorder particles to tree-walk order (not implemented yet)


    // Iterate if we need to immediately change SPH particle timesteps
    // (e.g. due to feedback, or sudden change in neighbour timesteps)
    //---------------------------------------------------------------------------------------------
    do {

      // Update cells containing active particles
      if (activecount > 0) sphneib->UpdateActiveParticleCounters(partdata,sph);

      // Calculate all SPH properties
      sphneib->UpdateAllSphProperties(sph->Nhydro,sph->Ntot,partdata,sph,nbody);

      // Update the radiation field
      /*if (Nsteps%1 == 0) {
        radiation->UpdateRadiationField(sph->Nhydro, nbody->Nnbody, sinks.Nsink,
                                        partdata, nbody->nbodydata, sinks.sink);
        for (i=0; i<sph->Nhydro; i++) {
          SphParticle<ndim>& part = sph->GetSphParticlePointer(i);
          sph->ComputeThermalProperties(part);
        }
      }*/

      // Update the radiation field
      for (int jj=0; jj<2; jj++) {
        radiation->UpdateRadiationField(sph->Nhydro, nbody->Nnbody, sinks.Nsink,
                                        partdata, nbody->nbodydata, sinks.sink);
        for (i=0; i<sph->Nhydro; i++) {
          SphParticle<ndim>& part = sph->GetSphParticlePointer(i);
          sph->ComputeThermalProperties(part);
        }
      }


      // Copy properties from original particles to ghost particles
      LocalGhosts->CopySphDataToGhosts(simbox,sph);

      // Calculate gravitational forces from other distant MPI nodes.
      // Also determines particles that must be exported to other nodes
      // if too close to the domain boundaries
#ifdef MPI_PARALLEL
      if (sph->self_gravity == 1)
        sphneib->UpdateGravityExportList(rank,sph->Nhydro,sph->Ntot,
                                         sph->GetSphParticleArray(),sph,nbody);
      else
        sphneib->UpdateHydroExportList(rank,sph->Nhydro,sph->Ntot,
                                       sph->GetSphParticleArray(),sph,nbody);

      // If active particles need forces from other domains, export particles
      mpicontrol->ExportParticlesBeforeForceLoop(sph);
#endif


      // Calculate SPH gravity and hydro forces, depending on which are activated
      if (ewaldGravity && sph->hydro_forces == 1 && sph->self_gravity == 1)
        sphneib->UpdateAllSphPeriodicForces(sph->Nhydro,sph->Ntot,partdata,sph,nbody,simbox,ewald);
      else if (ewaldGravity && sph->self_gravity == 1)
        sphneib->UpdateAllSphPeriodicGravForces(sph->Nhydro,sph->Ntot,partdata,
                                                sph,nbody,simbox,ewald);
      else if (sph->hydro_forces == 1 && sph->self_gravity == 1)
        sphneib->UpdateAllSphForces(sph->Nhydro,sph->Ntot,partdata,sph,nbody);
      else if (sph->hydro_forces == 1)
        sphneib->UpdateAllSphHydroForces(sph->Nhydro,sph->Ntot,partdata,sph,nbody);
      else if (sph->self_gravity == 1)
        sphneib->UpdateAllSphGravForces(sph->Nhydro,sph->Ntot,partdata,sph,nbody);


      // Check if all neighbouring timesteps are acceptable
      if (Nlevels > 1)
        activecount = sphint->CheckTimesteps(level_diff_max,level_step,n,
                                             sph->Nhydro,sph->GetSphParticleArray());
      else activecount = 0;
      //activecount = 0;

#if defined MPI_PARALLEL
      mpicontrol->GetExportedParticlesAccelerations(sph);
#endif


    } while (activecount > 0);
    //---------------------------------------------------------------------------------------------

  }
  //-----------------------------------------------------------------------------------------------


  // Compute N-body forces
  //-----------------------------------------------------------------------------------------------
  if (nbody->Nnbody > 0) {

    // Iterate for P(EC)^n schemes
    //---------------------------------------------------------------------------------------------
    for (it=0; it<nbody->Npec; it++) {

      // Zero all acceleration terms
      for (i=0; i<nbody->Nnbody; i++) {
        if (nbody->nbodydata[i]->active) {
          for (k=0; k<ndim; k++) nbody->nbodydata[i]->a[k] = 0.0;
          for (k=0; k<ndim; k++) nbody->nbodydata[i]->adot[k] = 0.0;
          for (k=0; k<ndim; k++) nbody->nbodydata[i]->a2dot[k] = 0.0;
          for (k=0; k<ndim; k++) nbody->nbodydata[i]->a3dot[k] = 0.0;
          nbody->nbodydata[i]->gpot = 0.0;
          nbody->nbodydata[i]->gpe = 0.0;
        }
      }
      if (sink_particles == 1) {
        for (i=0; i<sinks.Nsink; i++) {
          if (sinks.sink[i].star->active) {
            for (k=0; k<ndim; k++) sinks.sink[i].fhydro[k] = 0.0;
          }
        }
      }

      if (sph->self_gravity == 1 && sph->Nhydro > 0) {
        sphneib->UpdateAllStarGasForces(sph->Nhydro,sph->Ntot,partdata,sph,nbody);
#if defined MPI_PARALLEL
        // We need to sum up the contributions from the different domains
        mpicontrol->ComputeTotalStarGasForces(nbody);
#endif
      }

      // Calculate forces, force derivatives etc.., for active stars/systems
      if (nbody->nbody_softening == 1) {
        nbody->CalculateDirectSmoothedGravForces(nbody->Nnbody,nbody->nbodydata);
      }
      else {
        nbody->CalculateDirectGravForces(nbody->Nnbody,nbody->nbodydata);
      }

      // Calculate correction step for all stars at end of step, except the
      // final iteration (since correction is computed in EndStep also).
      //if (it < nbody->Npec - 1)
      nbody->CorrectionTerms(n,nbody->Nnbody,t,timestep,nbody->nbodydata);

    }
    //---------------------------------------------------------------------------------------------

  }
  //-----------------------------------------------------------------------------------------------


  rebuild_tree = false;


  // End-step terms for all SPH particles
  if (sph->Nhydro > 0) {
    uint->EndTimestep(n,sph->Nhydro,(FLOAT) t,(FLOAT) timestep,sph->GetSphParticleArray());
    sphint->EndTimestep(n,sph->Nhydro,(FLOAT) t,(FLOAT) timestep,sph->GetSphParticleArray());
  }

  // End-step terms for all star particles
  if (nbody->Nstar > 0) nbody->EndTimestep(n,nbody->Nnbody,t,timestep,nbody->nbodydata);

  // Search for new sink particles (if activated) and accrete to existing sinks
  if (sink_particles == 1) {
    if (sinks.create_sinks == 1 && (rebuild_tree || Nfullsteps%ntreebuildstep == 0)) {
      sinks.SearchForNewSinkParticles(n,t,sph,nbody);
    }
    if (sinks.Nsink > 0) {
      sinks.AccreteMassToSinks(sph,nbody,n,timestep);
      nbody->UpdateStellarProperties();
      if (extra_sink_output) WriteExtraSinkOutput();
    }
    // If we will output a snapshot (regular or for restarts), then delete all accreted particles
    if ((t >= tsnapnext && sinks.Nsink > 0) || n == nresync || kill_simulation ||
         timing->WallClockTime() - timing->tstart_wall > 0.99*tmax_wallclock) {
      sph->DeleteDeadParticles();
      rebuild_tree = true;
    }
  }


  return;
}



//=================================================================================================
//  SphSimulation::ComputeGlobalTimestep
/// Computes global timestep for SPH simulation.  Calculates the minimum
/// timestep for all SPH and N-body particles in the simulation.
//=================================================================================================
template <int ndim>
void SphSimulation<ndim>::ComputeGlobalTimestep(void)
{
  int i;                               // Particle counter
  DOUBLE dt;                           // Particle timestep
  DOUBLE dt_min = big_number_dp;       // Local copy of minimum timestep
  DOUBLE dt_nbody;                     // Aux. minimum N-body timestep
  DOUBLE dt_sph;                       // Aux. minimum SPH timestep

  debug2("[SphSimulation::ComputeGlobalTimestep]");
  timing->StartTimingSection("GLOBAL_TIMESTEPS",2);


  // Only update timestep when all particles are synced at end of last step.
  //-----------------------------------------------------------------------------------------------
  if (n == nresync) {

    n            = 0;
    level_max    = 0;
    level_step   = level_max + integration_step - 1;
    nresync      = integration_step;
    dt_min_nbody = big_number_dp;
    dt_min_sph   = big_number_dp;

    // Find minimum timestep from all SPH particles
    //---------------------------------------------------------------------------------------------
#pragma omp parallel default(none) private(i,dt,dt_nbody,dt_sph) \
  shared(dt_min) //,dt_min_nbody,dt_min_sph)
    {
      dt       = big_number_dp;
      dt_nbody = big_number_dp;
      dt_sph   = big_number_dp;

#pragma omp for
      for (i=0; i<sph->Nhydro; i++) {
        SphParticle<ndim>& part = sph->GetSphParticlePointer(i);
        part.level     = 0;
        part.levelneib = 0;
        part.nstep     = pow(2,level_step - part.level);
        part.nlast     = n;
        part.tlast     = t;
        part.dt        = sphint->Timestep(part,sph);
        dt             = min(dt,part.dt);
        dt_sph         = min(dt_sph,part.dt);
      }

      // Now compute minimum timestep due to stars/systems
#pragma omp for
      for (i=0; i<nbody->Nnbody; i++) {
        nbody->nbodydata[i]->level = 0;
        nbody->nbodydata[i]->nstep = pow(2,level_step - nbody->nbodydata[i]->level);
        nbody->nbodydata[i]->nlast = n;
        nbody->nbodydata[i]->tlast = t;
        nbody->nbodydata[i]->dt    = nbody->Timestep(nbody->nbodydata[i]);
        dt       = min(dt,nbody->nbodydata[i]->dt);
        dt_nbody = min(dt_nbody,nbody->nbodydata[i]->dt);
      }

#pragma omp critical
      {
        if (dt < dt_min) dt_min = dt;
        if (dt_sph < dt_min_sph) dt_min_sph = dt_sph;
        if (dt_nbody < dt_min_nbody) dt_min_nbody = dt_nbody;
      }

    }
    //---------------------------------------------------------------------------------------------

#ifdef MPI_PARALLEL
    dt = dt_min;
    MPI_Allreduce(&dt,&dt_min,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
#endif
    timestep = dt_min;

    // Set minimum timestep for all SPH and N-body particles
    for (i=0; i<sph->Nhydro; i++) sph->GetSphParticlePointer(i).dt = timestep;
    for (i=0; i<nbody->Nnbody; i++) nbody->nbodydata[i]->dt = timestep;

  }
  //-----------------------------------------------------------------------------------------------

  timing->EndTimingSection("GLOBAL_TIMESTEPS");

  return;
}



//=================================================================================================
//  SphSimulation::ComputeBlockTimesteps
/// Compute timesteps for all particles using hierarchical block timesteps.
//=================================================================================================
template <int ndim>
void SphSimulation<ndim>::ComputeBlockTimesteps(void)
{
  int i;                                // Particle counter
  int imin;                             // i.d. of ptcl with minimum timestep
  int imin_aux;                         // ..
  int istep;                            // Aux. variable for changing steps
  int last_level;                       // Previous timestep level
  int level;                            // Particle timestep level
  int level_max_aux;                    // Aux. maximum level variable
  int level_max_nbody = 0;              // level_max for star particles only
  int level_max_old;                    // Old level_max
  int level_max_sph = 0;                // level_max for SPH particles only
  int level_min_sph = 9999999;          // level_min for SPH particles
  int level_nbody;                      // local thread var. for N-body level
  int level_sph;                        // local thread var. for SPH level
  int nfactor;                          // Increase/decrease factor of n
  int nstep;                            // Particle integer step-size
  DOUBLE dt;                            // Aux. timestep variable
  DOUBLE dt_min = big_number_dp;        // Minimum timestep
  DOUBLE dt_min_aux;                    // Aux. minimum timestep variable
  DOUBLE dt_nbody;                      // Aux. minimum N-body timestep
  DOUBLE dt_sph;                        // Aux. minimum SPH timestep

  debug2("[SphSimulation::ComputeBlockTimesteps]");
  timing->StartTimingSection("BLOCK_TIMESTEPS",2);


  dt_min_nbody = big_number_dp;
  dt_min_sph = big_number_dp;


  // Synchronise all timesteps and reconstruct block timestep structure.
  //===============================================================================================
  if (n == nresync) {

    n = 0;
    timestep = big_number_dp;

#pragma omp parallel default(none) shared(imin) private(dt,dt_min_aux,dt_nbody,dt_sph,i,imin_aux)
    {
      // Initialise all timestep and min/max variables
      dt_min_aux = big_number_dp;
      dt_sph     = big_number_dp;
      dt_nbody   = big_number_dp;

      // Find minimum timestep from all SPH particles
#pragma omp for
      for (i=0; i<sph->Nhydro; i++) {
        SphParticle<ndim>& part = sph->GetSphParticlePointer(i);
        if (part.itype == dead) continue;
        dt = sphint->Timestep(part,sph);
        if (dt < dt_sph) imin_aux = i;
        dt_min_aux = min(dt_min_aux,dt);
        dt_sph     = min(dt_sph,dt);
        part.dt    = dt;
      }

      // Now compute minimum timestep due to stars/systems
#pragma omp for
      for (i=0; i<nbody->Nnbody; i++) {
        dt         = nbody->Timestep(nbody->nbodydata[i]);
        dt_min_aux = min(dt_min_aux,dt);
        dt_nbody   = min(dt_nbody,dt);
        nbody->nbodydata[i]->dt = dt;
      }

#pragma omp critical
      {
        if (dt_min_aux < timestep) imin = imin_aux;
        timestep     = min(timestep,dt_min_aux);
        dt_min_sph   = min(dt_min_sph,dt_sph);
        dt_min_nbody = min(dt_min_nbody,dt_nbody);
      }
#pragma omp barrier
    }


    // For MPI, determine the global minimum timestep over all processors
#ifdef MPI_PARALLEL
    dt = timestep;
    MPI_Allreduce(&dt,&timestep,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
    dt = dt_min_sph;
    MPI_Allreduce(&dt,&dt_min_sph,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
    dt = dt_min_nbody;
	MPI_Allreduce(&dt,&dt_min_nbody,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
#endif


    // Calculate new block timestep levels
    level_max  = Nlevels - 1;
    level_step = level_max + integration_step - 1;
    dt_max     = timestep*powf(2.0,level_max);

    // Calculate the maximum level occupied by all SPH particles
    level_max_sph   = min((int) (invlogetwo*log(dt_max/dt_min_sph)) + 1, level_max);
    level_max_nbody = min((int) (invlogetwo*log(dt_max/dt_min_nbody)) + 1, level_max);

    // Populate timestep levels with N-body particles.
    // Ensures that N-body particles occupy levels lower than all SPH particles
    for (i=0; i<nbody->Nnbody; i++) {
      dt = nbody->nbodydata[i]->dt;
      level = min((int) (invlogetwo*log(dt_max/dt)) + 1, level_max);
      level = max(level,0);
      nbody->nbodydata[i]->level = max(level,level_max_sph);
      nbody->nbodydata[i]->nlast = n;
      nbody->nbodydata[i]->nstep = pow(2,level_step - nbody->nbodydata[i]->level);
      nbody->nbodydata[i]->tlast = t;
    }

    // If particles are sink neighbours, set to same timesteps as sinks
    for (i=0; i<sph->Nhydro; i++) {
      SphParticle<ndim>& part = sph->GetSphParticlePointer(i);
      if (part.sinkid != -1) {
        if (sinks.sink[part.sinkid].star->level - part.level > level_diff_max) {
          part.level     = sinks.sink[part.sinkid].star->level - level_diff_max;
          part.levelneib = sinks.sink[part.sinkid].star->level;
          level_max_sph  = max(level_max_sph,part.level);
        }
      }
    }

    // If enforcing a single SPH timestep, set it here.
    // Otherwise, populate the timestep levels with SPH particles.
    if (sph_single_timestep == 1) {
      for (i=0; i<sph->Nhydro; i++) {
        SphParticle<ndim>& part = sph->GetSphParticlePointer(i);
        if (part.itype == dead) continue;
        part.level     = level_max_sph;
        part.levelneib = level_max_sph;
        part.nlast     = n;
        part.tlast     = t;
        part.nstep     = pow(2,level_step - part.level);
      }
      level_min_sph = level_max_sph;
    }
    else {
      for (i=0; i<sph->Nhydro; i++) {
        SphParticle<ndim>& part = sph->GetSphParticlePointer(i);
        if (part.itype == dead) continue;
        dt             = part.dt;
        level          = min((int) (invlogetwo*log(dt_max/dt)) + 1, level_max);
        level          = max(level,0);
        part.level     = level;
        part.levelneib = level;
        part.nlast     = n;
        part.tlast     = t;
        part.nstep     = pow(2,level_step - part.level);
        level_min_sph  = min(level_min_sph,part.level);
      }
    }


    nresync = pow(2,level_step);
    assert(nresync>0);
    timestep = dt_max / (DOUBLE) nresync;

  }
  // If not resynchronising, check if any SPH/N-body particles need to move
  // up or down timestep levels.
  //===============================================================================================
  else {

    level_max_old   = level_max;
    level_max       = 0;
    level_max_nbody = 0;
    level_max_sph   = 0;


#pragma omp parallel default(none) private(dt,dt_nbody,dt_sph,i,imin_aux) \
  private(istep,last_level,level,level_max_aux,level_nbody,level_sph,nstep,nfactor) \
  shared(dt_min,imin,level_max_nbody,level_max_sph,level_min_sph)
    {
      dt_sph        = big_number_dp;
      dt_nbody      = big_number_dp;
      level_max_aux = 0;
      level_nbody   = 0;
      level_sph     = 0;


      // Find all SPH particles at the beginning of a new timestep
      //-------------------------------------------------------------------------------------------
#pragma omp for
      for (i=0; i<sph->Nhydro; i++) {
        SphParticle<ndim>& part = sph->GetSphParticlePointer(i);
        if (part.itype == dead) continue;

        // SPH particles whose timestep has been artificially reduced by Saitoh & Makino scheme.
        if (part.nlast == n && part.nstep != pow(2,level_step - part.level)) {
          dt    = sphint->Timestep(part,sph);
          level = max((int) (invlogetwo*log(dt_max/dt)) + 1, 0);
          level = max(level,part.levelneib - level_diff_max);

          part.level     = max(part.level,level);
          part.levelneib = part.level;
          part.dt        = dt;
          part.nlast     = n;
          part.tlast     = t;
          part.nstep     = pow(2,level_step - part.level);
        }
        // SPH particles that have naturally reached the end of their step
        else if (part.nlast == n) {
          nstep      = part.nstep;
          last_level = part.level;

          // Compute new timestep value and level number
          dt    = sphint->Timestep(part,sph);
          level = max((int) (invlogetwo*log(dt_max/dt)) + 1, 0);
          level = max(level,part.levelneib - level_diff_max);

          // Move up one level (if levels are correctly synchronised) or
          // down several levels if required
          if (level < last_level && last_level > 1 && n%(2*nstep) == 0)
            part.level = last_level - 1;
          else if (level > last_level)
            part.level = level;
          else
            part.level = last_level;

          part.levelneib = level;
          part.dt        = dt;
          part.nlast     = n;
          part.tlast     = t;
          part.nstep     = pow(2,level_step - part.level);
        }

        // Find maximum level of all SPH particles
        level_sph = max(level_sph,part.level);
        if (part.dt < dt_sph) imin_aux = i;
        level_max_aux = max(level_max_aux,part.level);

        dt_sph = min(dt_sph,part.dt);
      }
      //-------------------------------------------------------------------------------------------


#pragma omp critical
      {
        if (dt_sph < dt_min_sph) imin = imin_aux;
        dt_min        = min(dt_min,dt_sph);
        dt_min_sph    = min(dt_min_sph,dt_sph);
        level_max     = max(level_max,level_max_aux);
        level_max_sph = max(level_max_sph,level_sph);
      }
#pragma omp barrier


      // Now find all N-body particles at the beginning of a new timestep
      //-------------------------------------------------------------------------------------------
#pragma omp for
      for (i=0; i<nbody->Nnbody; i++) {

        // Skip particles that are not at end of step
        if (nbody->nbodydata[i]->nlast == n) {
          nstep = nbody->nbodydata[i]->nstep;
          last_level = nbody->nbodydata[i]->level;

          // Compute new timestep value and level number
          dt    = nbody->Timestep(nbody->nbodydata[i]);
          level = max((int) (invlogetwo*log(dt_max/dt)) + 1, 0);
          level = max(level,level_max_sph);

          // Move up one level (if levels are correctly synchronised) or
          // down several levels if required
          if (level < last_level && level > level_max_sph && last_level > 1 && n%(2*nstep) == 0)
            nbody->nbodydata[i]->level = last_level - 1;
          else if (level > last_level)
            nbody->nbodydata[i]->level = level;
          else
            nbody->nbodydata[i]->level = last_level;

          nbody->nbodydata[i]->dt    = dt;
          nbody->nbodydata[i]->nlast = n;
          nbody->nbodydata[i]->nstep = pow(2,level_step - nbody->nbodydata[i]->level);
          nbody->nbodydata[i]->tlast = t;
        }

        // Find maximum level of all N-body particles
        level_nbody   = max(level_nbody,nbody->nbodydata[i]->level);
        level_max_aux = max(level_max_aux,nbody->nbodydata[i]->level);
        dt_nbody      = min(dt_nbody,nbody->nbodydata[i]->dt);
      }
      //-------------------------------------------------------------------------------------------

      #pragma omp critical
      {
        dt_min          = min(dt_min,dt_nbody);
        dt_min_nbody    = min(dt_min_nbody,dt_nbody);
        level_max       = max(level_max,level_max_aux);
        level_max_nbody = max(level_max_nbody,level_nbody);
      }
      #pragma omp barrier


      // Correct timestep levels for any particles that have entered a sink
      //-------------------------------------------------------------------------------------------
  /*#pragma omp for
      for (i=0; i<sph->Nhydro; i++) {
        SphParticle<ndim>& part = sph->GetSphParticlePointer(i);
        if (part.itype == dead || part.nlast != n) continue;
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
    MPI_Allreduce(&level,&level_max,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
    level = level_max_sph;
    MPI_Allreduce(&level,&level_max_sph,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
    level = level_max_nbody;
    MPI_Allreduce(&level,&level_max_nbody,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
    assert(level_max_sph>=0);
#endif


    // Set fixed SPH timestep level here in case maximum has changed
    if (sph_single_timestep == 1) {
      for (i=0; i<sph->Nhydro; i++) {
        SphParticle<ndim>& part = sph->GetSphParticlePointer(i);
        if (part.itype == dead) continue;
        if (part.nlast == n) part.level = level_max_sph;
      }
    }

    // Update all timestep variables if we have removed or added any levels
    //---------------------------------------------------------------------------------------------
    if (level_max != level_max_old) {

      // Increase maximum timestep level if correctly synchronised
      istep = pow(2,level_step - level_max_old + 1);
      if (level_max <= level_max_old - 1 && level_max_old > 1 && n%istep == 0) {
        level_max = level_max_old - 1;
      }
      else if (level_max < level_max_old) {
        level_max = level_max_old;
      }

      level_step = level_max + integration_step - 1;

      // Adjust integer time if levels added or removed
      if (level_max > level_max_old) {
        nfactor = pow(2,level_max - level_max_old);
        n *= nfactor;
        for (i=0; i<sph->Nhydro; i++) {
          SphParticle<ndim>& part = sph->GetSphParticlePointer(i);
          if (part.itype == dead) continue;
          part.nstep *= nfactor;
          part.nlast *= nfactor;
        }
        for (i=0; i<nbody->Nnbody; i++) nbody->nbodydata[i]->nstep *= nfactor;
        for (i=0; i<nbody->Nnbody; i++) nbody->nbodydata[i]->nlast *= nfactor;
      }
      else if (level_max < level_max_old) {
        nfactor = pow(2,level_max_old - level_max);
        n /= nfactor;
        for (i=0; i<sph->Nhydro; i++) {
          SphParticle<ndim>& part = sph->GetSphParticlePointer(i);
          if (part.itype == dead) continue;
          part.nlast /= nfactor;
          part.nstep /= nfactor;
        }
        for (i=0; i<nbody->Nnbody; i++) nbody->nbodydata[i]->nlast /= nfactor;
        for (i=0; i<nbody->Nnbody; i++) nbody->nbodydata[i]->nstep /= nfactor;
      }

    }
    //---------------------------------------------------------------------------------------------

    level_step = level_max + integration_step - 1;
    nresync    = pow(2,level_step);
    timestep   = dt_max / (DOUBLE) nresync;

    // Update values of nstep for both SPH and star particles
    for (i=0; i<sph->Nhydro; i++) {
      SphParticle<ndim>& part = sph->GetSphParticlePointer(i);
      if (part.itype == dead) continue;
      if (part.nlast == n) part.nstep = pow(2,level_step - part.level);
    }
    for (i=0; i<nbody->Nnbody; i++) {
      if (nbody->nbodydata[i]->nlast == n) nbody->nbodydata[i]->nstep =
        pow(2,level_step - nbody->nbodydata[i]->level);
    }

    assert(level_max >= level_max_old - 1);

  }
  //===============================================================================================


  // Various asserts for debugging
  assert(timestep >= 0.0);
  assert(level_step == level_max + integration_step - 1);
  assert(level_max_sph <= level_max);
  assert(level_max_nbody <= level_max);
  assert(n <= nresync);
  for (i=0; i<sph->Nhydro; i++) {
    SphParticle<ndim>& part = sph->GetSphParticlePointer(i);
    if (part.itype == dead) continue;
    assert(part.level <= level_max);
    assert(part.nlast <= n);
    assert(part.tlast <= t);
    assert(part.nstep == pow(2,level_step - part.level));
    assert(part.nlast != n || n%part.nstep == 0);
  }
  for (i=0; i<nbody->Nnbody; i++) {
    assert(nbody->nbodydata[i]->level <= level_max);
    assert(nbody->nbodydata[i]->nlast <= n);
    assert(nbody->nbodydata[i]->nstep == pow(2,level_step - nbody->nbodydata[i]->level));
    assert(nbody->nbodydata[i]->nlast != n || n%nbody->nbodydata[i]->nstep == 0);
    assert(nbody->nbodydata[i]->level >= level_max_sph);
    assert(nbody->nbodydata[i]->tlast <= t);
  }
  if (timestep <= 0.0) {
    cout << "Timestep fallen to zero : " << timestep << endl;
    exit(0);
  }


  timing->EndTimingSection("BLOCK_TIMESTEPS");

  return;



  // Some validations
  //-----------------------------------------------------------------------------------------------
  int *ninlevel;
  int Nactive=0;
  ninlevel = new int[level_max+1];
  SphParticle<ndim>& part = sph->GetSphParticlePointer(imin);
  cout << "-----------------------------------------------------" << endl;
  cout << "Checking timesteps : " << level_max << "   " << level_max_sph << "    "
       << level_max_nbody << "    " << level_step << "   " << level_max_old << endl;
  cout << "n : " << n << endl;
  cout << "dt_min_sph : " << dt_min_sph << "    dt_min_nbody : " << dt_min_nbody
       << "    timestep : " << timestep << endl;
  cout << "imin : " << imin << "    " << part.dt << "     " << part.m/sph->mmean << "    "
       << part.h << "    " << "    " << part.sound << "     " << part.div_v << "     "
       << part.h/(part.sound + part.h*fabs(part.div_v)) << "     "
       << sqrt(part.h/sqrt(DotProduct(part.a,part.a,ndim)))
       << endl;
  for (int l=0; l<=level_max; l++) ninlevel[l] = 0;
  for (i=0; i<sph->Nhydro; i++) {
    SphParticle<ndim>& part = sph->GetSphParticlePointer(i);
    if (part.active) Nactive++;
    ninlevel[part.level]++;
  }
  cout << "No. of active SPH particles : " << Nactive << endl;
  cout << "SPH level occupancy" << endl;
  for (int l=0; l<=level_max; l++) cout << "level : " << l << "     N : " << ninlevel[l] << endl;

  for (int l=0; l<=level_max; l++) ninlevel[l] = 0;
  for (i=0; i<nbody->Nstar; i++) ninlevel[nbody->nbodydata[i]->level]++;
  cout << "N-body level occupancy" << endl;
  for (int l=0; l<=level_max; l++) cout << "level : " << l << "     N : " << ninlevel[l] << endl;

  delete[] ninlevel;

  if (timestep <= 0.0) {
    cout << "Timestep fallen to zero : " << timestep << endl;
    exit(0);
  }


  return;

}



//=================================================================================================
//  SphSimulation::WriteExtraSinkOutput
/// For any simulations loaded into memory via a snapshot file, all particle
/// variables are converted into dimensionless code units here.
//=================================================================================================
template <int ndim>
void SphSimulation<ndim>::WriteExtraSinkOutput(void)
{
  int k;                               // ..
  int s;                               // ..
  string filename;                     // Output snapshot filename
  string nostring;                     // String of number of snapshots
  stringstream ss;                     // Stream object for preparing filename
  ofstream outfile;                    // Stream of restart file


  //-----------------------------------------------------------------------------------------------
  for (s=0; s<sinks.Nsink; s++) {

    SinkParticle<ndim> &sink = sinks.sink[s];
    nostring = "";
    ss << setfill('0') << setw(5) << s;
    nostring = ss.str();
    filename = run_id + ".sink." + nostring;
    ss.str(std::string());

    outfile.open(filename.c_str(), std::ofstream::app);
    outfile << t << "    ";
    outfile << Nsteps << "    ";
    for (k=0; k<ndim; k++) outfile << sink.star->r[k] << "    ";
    for (k=0; k<ndim; k++) outfile << sink.star->v[k] << "    ";
    for (k=0; k<ndim; k++) outfile << sink.star->a[k] << "    ";
    for (k=0; k<3; k++) outfile << sink.angmom[k] << "    ";
    outfile << sink.star->m  << "    ";
    outfile << sink.menc     << "    ";
    outfile << sink.mmax     << "    ";
    outfile << sink.macctot  << "    ";
    outfile << sink.dmdt     << "    ";
    outfile << sink.ketot    << "    ";
    outfile << sink.gpetot   << "    ";
    outfile << sink.rotketot << "    ";
    outfile << sink.utot     << "    ";
    outfile << sink.taccrete << "    ";
    outfile << sink.trad     << "    ";
    outfile << sink.trot     << "    ";
    outfile << sink.tvisc    << "    ";
    outfile << endl;
    outfile.close();

  }
  //-----------------------------------------------------------------------------------------------

  return;
}
