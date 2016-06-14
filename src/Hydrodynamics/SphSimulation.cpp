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
#if defined(NDIM_1)
template class SphSimulation<1>;
#endif
#if defined(NDIM_2)
template class SphSimulation<2>;
#endif
#if defined(NDIM_3)
template class SphSimulation<3>;
#endif



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
    std::ostringstream message;
    message << "Invalid dimensionality chosen : ndim = " << ndim;
    ExceptionHandler::getIstance().raise(message.str());
  }

  // Set-up random number generator object
  //-----------------------------------------------------------------------------------------------
  if (stringparams["rand_algorithm"] == "xorshift") {
    randnumb = new XorshiftRand(intparams["randseed"]);
  }
  else if (stringparams["rand_algorithm"] == "none") {
    randnumb = new DefaultSystemRand(intparams["randseed"]);
  }
  else {
    string message = "Unrecognised parameter : rand_algorithm= " +
      stringparams["rand_algorithm"];
    ExceptionHandler::getIstance().raise(message);
  }


  // Set-up all output units for scaling parameters
  simunits.SetupUnits(simparams);

  // Boundary condition variables
  //-----------------------------------------------------------------------------------------------
  simbox.boundary_lhs[0] = setBoundaryType(stringparams["boundary_lhs[0]"]);
  simbox.boundary_rhs[0] = setBoundaryType(stringparams["boundary_rhs[0]"]);
  simbox.boxmin[0] = floatparams["boxmin[0]"]/simunits.r.outscale;
  simbox.boxmax[0] = floatparams["boxmax[0]"]/simunits.r.outscale;

  if (ndim > 1) {
    simbox.boundary_lhs[1] = setBoundaryType(stringparams["boundary_lhs[1]"]);
    simbox.boundary_rhs[1] = setBoundaryType(stringparams["boundary_rhs[1]"]);
    simbox.boxmin[1] = floatparams["boxmin[1]"]/simunits.r.outscale;
    simbox.boxmax[1] = floatparams["boxmax[1]"]/simunits.r.outscale;
  }

  if (ndim == 3) {
    simbox.boundary_lhs[2] = setBoundaryType(stringparams["boundary_lhs[2]"]);
    simbox.boundary_rhs[2] = setBoundaryType(stringparams["boundary_rhs[2]"]);
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


  // Set external potential field object and set pointers to object
  if (stringparams["external_potential"] == "none") {
    extpot = new NullPotential<ndim>();
  }
  else if (stringparams["external_potential"] == "vertical") {
    extpot = new VerticalPotential<ndim>
      (intparams["kgrav"], floatparams["avert"], simbox.boxmin[intparams["kgrav"]]);
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
  periodicBoundaries = IsAnyBoundaryPeriodic(simbox);
  if (periodicBoundaries && intparams["self_gravity"] == 1) {
    ewaldGravity = true;
    ewald = new Ewald<ndim>
      (simbox, intparams["gr_bhewaldseriesn"], intparams["in"], intparams["nEwaldGrid"],
       floatparams["ewald_mult"], floatparams["ixmin"], floatparams["ixmax"],
       floatparams["EFratio"], timing);
    simbox.PeriodicGravity = true ;
  }
  else{
    simbox.PeriodicGravity = false ;
    if (IsAnyBoundaryReflecting(simbox) && intparams["self_gravity"]){
      ExceptionHandler::getIstance().raise("Error: Reflecting boundaries and self-gravity is not "
    		                               "supported") ;
    }
  }


  // Set all other SPH parameter variables
  sph->Nhydromax       = intparams["Nhydromax"];
  sph->create_sinks    = intparams["create_sinks"];
  sph->fixed_sink_mass = intparams["fixed_sink_mass"];
  sph->msink_fixed     = floatparams["m1"];
  sph->alpha_visc_min  = floatparams["alpha_visc_min"];


  // Set important variables for N-body objects
  nbody->Nstarmax       = intparams["Nstarmax"];
  nbody_single_timestep = intparams["nbody_single_timestep"];
  nbodytree.gpehard     = floatparams["gpehard"];
  nbodytree.gpesoft     = floatparams["gpesoft"];
  //nbody->perturbers     = intparams["perturbers"];
  //if (intparams["sub_systems"] == 1) subsystem->perturbers = intparams["perturbers"];


  // Sink particles
  //-----------------------------------------------------------------------------------------------
  sinks = new Sinks<ndim>(sphneib);
  sink_particles             = intparams["sink_particles"];
  sinks->sink_particles      = intparams["sink_particles"];
  sinks->create_sinks        = intparams["create_sinks"];
  sinks->smooth_accretion    = intparams["smooth_accretion"];
  sinks->alpha_ss            = floatparams["alpha_ss"];
  sinks->smooth_accrete_frac = floatparams["smooth_accrete_frac"];
  sinks->smooth_accrete_dt   = floatparams["smooth_accrete_dt"];
  sinks->sink_radius_mode    = stringparams["sink_radius_mode"];
  sinks->rho_sink            = floatparams["rho_sink"];
  sinks->rho_sink            /= (simunits.rho.outscale*simunits.rho.outcgs);

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

#if defined MPI_PARALLEL
  sinks->SetMpiControl(mpicontrol);
#endif

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
  nradstep            = intparams["nradstep"];
  nrestartstep        = intparams["nrestartstep"];
  ntreebuildstep      = intparams["ntreebuildstep"];
  ntreestockstep      = intparams["ntreestockstep"];
  Nstepsmax           = intparams["Nstepsmax"];
  out_file_form       = stringparams["out_file_form"];
  pruning_level_min   = intparams["pruning_level_min"];
  pruning_level_max   = intparams["pruning_level_max"];
  run_id              = stringparams["run_id"];
  sph_single_timestep = intparams["sph_single_timestep"];
  tmax_wallclock      = floatparams["tmax_wallclock"];
  tend                = floatparams["tend"]/simunits.t.outscale;
  tlitesnapnext       = floatparams["tlitesnapfirst"]/simunits.t.outscale;
  tsnapnext           = floatparams["tsnapfirst"]/simunits.t.outscale;


  // Set pointers to timing object
  nbody->timing   = timing;
  if (sim == "sph" || sim == "gradhsph" || sim == "sm2012sph") {
    sinks->timing    = timing;
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
  FLOAT adot[ndim];                    // Dummy adot array
  SphParticle<ndim> *partdata;         // Pointer to main SPH data array

  debug2("[SphSimulation::PostInitialConditionsSetup]");

  // Set iorig
  if (rank == 0) {
    for (i=0; i<sph->Nhydro; i++) sph->GetSphParticlePointer(i).iorig = i;
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
  mpicontrol->CreateInitialDomainDecomposition(sph, nbody, simparams, simbox,
                                               this->initial_h_provided);
  this->AllocateParticleMemory();
#endif

  // Set pointer to SPH particle data
  partdata = sph->GetSphParticleArray();

  // Set time variables here (for now)
  nresync = 0;
  n = 0;

  // Set initial smoothing lengths and create initial ghost particles
  //-----------------------------------------------------------------------------------------------
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

  // Compute instantaneous mean mass (used for smooth sink accretion)
  sph->mmean = (FLOAT) 0.0;
  for (i=0; i<sph->Nhydro; i++) sph->mmean += sph->GetSphParticlePointer(i).m;
  sph->mmean /= (FLOAT) sph->Nhydro;

  // Compute minimum smoothing length of sinks
  sph->hmin_sink = big_number;
  for (i=0; i<sinks->Nsink; i++) {
    sph->hmin_sink = min(sph->hmin_sink, (FLOAT) sinks->sink[i].star->h);
  }

  // If the smoothing lengths have not been provided beforehand, then
  // calculate the initial values here
  //sphneib->neibcheck = false;
  if (!this->initial_h_provided) {
    sph->InitialSmoothingLengthGuess();
    sphneib->BuildTree(rebuild_tree, 0, ntreebuildstep, ntreestockstep,
                       sph->Ntot, sph->Nhydromax, timestep, partdata, sph);
    sphneib->UpdateAllSphProperties(sph->Nhydro, sph->Ntot, partdata, sph, nbody);
  }
  else {
    sphneib->BuildTree(rebuild_tree, 0, ntreebuildstep, ntreestockstep,
                       sph->Ntot, sph->Nhydromax, timestep, partdata, sph);
  }


#ifdef MPI_PARALLEL
  mpicontrol->UpdateAllBoundingBoxes(sph->Nhydro, sph, sph->kernp);
#endif

  // Search ghost particles
  sphneib->SearchBoundaryGhostParticles((FLOAT) 0.0, simbox, sph);
  sphneib->BuildGhostTree(true, 0, ntreebuildstep, ntreestockstep, sph->Ntot,
                          sph->Nhydromax, timestep, partdata, sph);
#ifdef MPI_PARALLEL
  mpicontrol->UpdateAllBoundingBoxes(sph->Nhydro+sph->NPeriodicGhost, sph, sph->kernp);
  for (int i=0; i<sph->Nhydro+sph->NPeriodicGhost; i++) {
    SphParticle<ndim>& parti = sph->GetSphParticlePointer(i);
    parti.hrangesqd = sph->kernfacsqd*sph->kernp->kernrangesqd*parti.h*parti.h;
  }
  MpiGhosts->SearchGhostParticles((FLOAT) 0.0, simbox, sph);
  sphneib->BuildMpiGhostTree(true, 0, ntreebuildstep, ntreestockstep, sph->Ntot,
                             sph->Nhydromax, timestep, partdata, sph);
#endif

  // Zero accelerations
  for (i=0; i<sph->Nhydro; i++) sph->GetSphParticlePointer(i).active = true;

  // Update neighbour tree
  rebuild_tree = true;
  sphneib->BuildTree(rebuild_tree, 0, ntreebuildstep, ntreestockstep,
                     sph->Ntot, sph->Nhydromax, timestep, partdata, sph);
  level_step = 1;


  // For Eigenvalue MAC, need non-zero values
  for (i=0; i<sph->Nhydro; i++) sph->GetSphParticlePointer(i).gpot = big_number;

  // Calculate all SPH properties
  sphneib->UpdateAllSphProperties(sph->Nhydro, sph->Ntot, partdata, sph, nbody);

#ifdef MPI_PARALLEL
  mpicontrol->UpdateAllBoundingBoxes(sph->Nhydro, sph, sph->kernp);
#endif

  // Regularise particle positions (if selected in parameters file)
  if (simparams->intparams["regularise_particle_ics"] == 1) {
    RegulariseParticleDistribution(simparams->intparams["Nreg"]);
  }

  // Search ghost particles
  sphneib->SearchBoundaryGhostParticles((FLOAT) 0.0, simbox, sph);
  sphneib->BuildGhostTree(true, 0, ntreebuildstep, ntreestockstep, sph->Ntot,
                          sph->Nhydromax, timestep, partdata, sph);
#ifdef MPI_PARALLEL
  mpicontrol->UpdateAllBoundingBoxes(sph->Nhydro + sph->NPeriodicGhost, sph, sph->kernp);
  MpiGhosts->SearchGhostParticles((FLOAT) 0.0, simbox, sph);
  sphneib->BuildMpiGhostTree(true, 0, ntreebuildstep, ntreestockstep, sph->Ntot,
                             sph->Nhydromax, timestep, partdata, sph);
#endif

  // Update neighbour tree
  rebuild_tree = true;
  sphneib->BuildTree(rebuild_tree, 0, ntreebuildstep, ntreestockstep,
                     sph->Ntot, sph->Nhydromax, timestep, partdata, sph);
  sphneib->SearchBoundaryGhostParticles((FLOAT) 0.0, simbox, sph);
  sphneib->BuildGhostTree(true, 0, ntreebuildstep, ntreestockstep, sph->Ntot,
                          sph->Nhydromax, timestep, partdata, sph);
  //sphneib->neibcheck = true;

    // Communicate pruned trees for MPI
#ifdef MPI_PARALLEL
  sphneib->BuildPrunedTree(rank, sph->Nhydromax, simbox, mpicontrol->mpinode, partdata);
#endif


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


  // Read-in stellar properties table here
  nbody->LoadStellarPropertiesTable(&simunits);
  nbody->UpdateStellarProperties();


  // Compute all initial SPH force terms
  //-----------------------------------------------------------------------------------------------

  // Zero accelerations (here for now)
  for (i=0; i<sph->Ntot; i++) {
    SphParticle<ndim>& part = sph->GetSphParticlePointer(i);
    part.tlast     = t;
    part.active    = false;
    part.level     = 0;
    part.levelneib = 0;
    part.nstep     = 0;
    part.nlast     = 0;
    part.dalphadt  = (FLOAT) 0.0;
    part.div_v     = (FLOAT) 0.0;
    part.dudt      = (FLOAT) 0.0;
    part.gpot      = (FLOAT) 0.0;
    part.mu_bar    = (FLOAT) simparams->floatparams["mu_bar"];
    for (k=0; k<ndim; k++) part.a[k] = (FLOAT) 0.0;
  }
  for (i=0; i<sph->Nhydro; i++) sph->GetSphParticlePointer(i).active = true;

  // Copy all other data from real SPH particles to ghosts
  LocalGhosts->CopyHydroDataToGhosts(simbox, sph);

  sphneib->BuildTree(true, 0, ntreebuildstep, ntreestockstep, sph->Ntot,
                     sph->Nhydromax, timestep, partdata, sph);
  sphneib->SearchBoundaryGhostParticles((FLOAT) 0.0, simbox, sph);
  sphneib->BuildGhostTree(true, 0, ntreebuildstep, ntreestockstep, sph->Ntot,
                          sph->Nhydromax, timestep, partdata, sph);
#ifdef MPI_PARALLEL
  mpicontrol->UpdateAllBoundingBoxes(sph->Nhydro + sph->NPeriodicGhost, sph, sph->kernp);
  MpiGhosts->SearchGhostParticles((FLOAT) 0.0, simbox, sph);
  sphneib->BuildMpiGhostTree(true, 0, ntreebuildstep, ntreestockstep, sph->Ntot,
                             sph->Nhydromax, timestep, partdata, sph);
#endif

  // Calculate all SPH properties
  sphneib->UpdateAllSphProperties(sph->Nhydro, sph->Ntot, partdata, sph, nbody);


#ifdef MPI_PARALLEL
  if (sph->self_gravity == 1) {
    sphneib->UpdateGravityExportList(rank, sph->Nhydro, sph->Ntot, partdata, sph, nbody, simbox);
  }
  else {
    sphneib->UpdateHydroExportList(rank, sph->Nhydro, sph->Ntot, partdata, sph, nbody, simbox);
  }

  mpicontrol->ExportParticlesBeforeForceLoop(sph);
#endif


  for (i=0; i<sph->Nhydro; i++) {
    SphParticle<ndim>& part = sph->GetSphParticlePointer(i);
    part.ionfrac = (FLOAT) 0.9999999;
  }
  // Update the radiation field
  for (int jj=0; jj<10; jj++) {
    radiation->UpdateRadiationField(sph->Nhydro, nbody->Nnbody, sinks->Nsink,
                                    partdata, nbody->nbodydata, sinks->sink);
  }


  // Update thermal properties (if radiation field has altered them)
  for (i=0; i<sph->Nhydro; i++) {
    SphParticle<ndim>& part = sph->GetSphParticlePointer(i);
    sph->ComputeThermalProperties(part);
  }


  // Calculate SPH gravity and hydro forces, depending on which are activated
  if (sph->hydro_forces == 1 && sph->self_gravity == 1) {
    sphneib->UpdateAllSphForces(sph->Nhydro, sph->Ntot, partdata, sph, nbody, simbox, ewald);
  }
  else if (sph->self_gravity == 1) {
    sphneib->UpdateAllSphGravForces(sph->Nhydro, sph->Ntot, partdata, sph, nbody, simbox, ewald);
  }
  else if (sph->hydro_forces == 1) {
    sphneib->UpdateAllSphHydroForces(sph->Nhydro, sph->Ntot, partdata, sph, nbody, simbox);
  }
  else{
    ExceptionHandler::getIstance().raise("Error: No forces included in simulation");
  }

#if defined MPI_PARALLEL
  mpicontrol->GetExportedParticlesAccelerations(sph);
#endif

  // Add external potential for all active SPH particles
  for (i=0; i<sph->Nhydro; i++) {
    SphParticle<ndim>& part = sph->GetSphParticlePointer(i);
    sph->extpot->AddExternalPotential(part.r, part.v, part.a, adot, part.gpot);
  }

  // Compute the dust forces if present.
  if (sphdust != NULL){
    // Copy properties from original particles to ghost particles
    LocalGhosts->CopyHydroDataToGhosts(simbox, sph);
#ifdef MPI_PARALLEL
    MpiGhosts->CopyHydroDataToGhosts(simbox, sph);
#endif
    sphdust->UpdateAllDragForces(sph->Nhydro, sph->Ntot, partdata) ;
  }

  // Set initial accelerations
  for (i=0; i<sph->Nhydro; i++) {
      SphParticle<ndim>& part = sph->GetSphParticlePointer(i);
      for (k=0; k<ndim; k++) part.r0[k] = part.r[k];
      for (k=0; k<ndim; k++) part.v0[k] = part.v[k];
      for (k=0; k<ndim; k++) part.a0[k] = part.a[k];
      part.active = false;
    }

  LocalGhosts->CopyHydroDataToGhosts(simbox,sph);
#ifdef MPI_PARALLEL
  MpiGhosts->CopyHydroDataToGhosts(simbox,sph);
#endif


  // Compute initial N-body forces
  //-----------------------------------------------------------------------------------------------
  if (sph->self_gravity == 1 && sph->Nhydro > 0) {
    sphneib->UpdateAllStarGasForces(sph->Nhydro, sph->Ntot, partdata, sph, nbody);
#if defined MPI_PARALLEL
    // We need to sum up the contributions from the different domains
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


  // Set particle values for initial step (e.g. r0, v0, a0, u0)
  uint->EndTimestep(n, sph->Nhydro, t, timestep, partdata);
  sphint->EndTimestep(n, sph->Nhydro, t, timestep, partdata);
  nbody->EndTimestep(n, nbody->Nstar, t, timestep, nbody->nbodydata);

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
  FLOAT adot[ndim];                    // Dummy adot variable
  FLOAT tghost;                        // Approx. ghost particle lifetime
  SphParticle<ndim> *partdata = sph->GetSphParticleArray();

  debug2("[SphSimulation::MainLoop]");

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
  uint->EnergyIntegration(n, sph->Nhydro, (FLOAT) t, (FLOAT) timestep, partdata);
  sphint->AdvanceParticles(n, sph->Nhydro, (FLOAT) t, (FLOAT) timestep, partdata);
  nbody->AdvanceParticles(n, nbody->Nnbody, t, timestep, nbody->nbodydata);

  // Check all boundary conditions
  // (DAVID : Move this function to sphint and create an analagous one
  //  for N-body.  Also, only check this on tree-build steps)
  if (Nsteps%ntreebuildstep == 0 || rebuild_tree) sphint->CheckBoundaries(simbox,sph);


  // Perform the load-balancing step for MPI simulations.  First update the pruned trees on all
  // processors, then compute the new load-balanced MPI domains and finally transfer the
  // particles to the new domains.
  //-----------------------------------------------------------------------------------------------
#ifdef MPI_PARALLEL
  if (Nsteps%ntreebuildstep == 0 || rebuild_tree) {
    sphneib->BuildPrunedTree(rank, sph->Nhydromax, simbox, mpicontrol->mpinode, partdata);
    mpicontrol->UpdateAllBoundingBoxes(sph->Nhydro, sph, sph->kernp);
    mpicontrol->LoadBalancing(sph, nbody);
    sphneib->InitialiseCellWorkCounters();
  }
#endif


  // Rebuild or update local neighbour and gravity tree
  sphneib->BuildTree(rebuild_tree, Nsteps, ntreebuildstep, ntreestockstep,
                     sph->Ntot, sph->Nhydromax, timestep, partdata, sph);

  // Search for new ghost particles and create on local processor
  //-----------------------------------------------------------------------------------------------
  if (Nsteps%ntreebuildstep == 0 || rebuild_tree) {
    tghost = timestep*(FLOAT)(ntreebuildstep - 1);
    sphneib->SearchBoundaryGhostParticles(tghost, simbox, sph);
    sphneib->BuildGhostTree(rebuild_tree, Nsteps, ntreebuildstep, ntreestockstep,
                            sph->Ntot, sph->Nhydromax, timestep, partdata, sph);

  // Re-build and communicate the new pruned trees (since the trees will necessarily change
  // once there has been communication of particles to new domains)
#ifdef MPI_PARALLEL
    sphneib->BuildPrunedTree(rank, sph->Nhydromax, simbox, mpicontrol->mpinode, partdata);
    mpicontrol->UpdateAllBoundingBoxes(sph->Nhydro + sph->NPeriodicGhost, sph, sph->kernp);
    MpiGhosts->SearchGhostParticles(tghost, simbox, sph);
    sphneib->BuildMpiGhostTree(rebuild_tree, Nsteps, ntreebuildstep, ntreestockstep,
                               sph->Ntot, sph->Nhydromax, timestep, partdata, sph);
#endif
  }
  // Otherwise copy properties from original particles to ghost particles
  else {
    LocalGhosts->CopyHydroDataToGhosts(simbox, sph);
#ifdef MPI_PARALLEL
    MpiGhosts->CopyHydroDataToGhosts(simbox, sph);
#endif
    sphneib->BuildGhostTree(rebuild_tree, Nsteps, ntreebuildstep, ntreestockstep,
                            sph->Ntot, sph->Nhydromax, timestep, partdata, sph);
  }


  // Iterate if we need to immediately change SPH particle timesteps
  // (e.g. due to feedback, or sudden change in neighbour timesteps)
  //-----------------------------------------------------------------------------------------------
  do {

    // Update cells containing active particles
    if (activecount > 0) sphneib->UpdateActiveParticleCounters(partdata, sph);

      // Zero accelerations (here for now)
      for (i=0; i<sph->Nhydro; i++) {
        SphParticle<ndim>& part = sph->GetSphParticlePointer(i);
        if (part.active) {
          part.levelneib = 0;
          part.dalphadt  = (FLOAT) 0.0;
          part.div_v     = (FLOAT) 0.0;
          part.dudt      = (FLOAT) 0.0;
          part.gpot      = (FLOAT) 0.0;
          for (k=0; k<ndim; k++) part.a[k] = (FLOAT) 0.0;
          //for (k=0; k<ndim; k++) part.a_dust[k] = (FLOAT) 0.0;
        }
      }

      // Calculate all SPH properties
      sphneib->UpdateAllSphProperties(sph->Nhydro, sph->Ntot, partdata, sph, nbody);

      // Update the radiation field
      if (Nsteps%nradstep == 0 || recomputeRadiation) {
        radiation->UpdateRadiationField(sph->Nhydro, nbody->Nnbody, sinks->Nsink,
                                        partdata, nbody->nbodydata, sinks->sink);
        for (i=0; i<sph->Nhydro; i++) {
          SphParticle<ndim>& part = sph->GetSphParticlePointer(i);
          sph->ComputeThermalProperties(part);
        }
      }

      // Copy properties from original particles to ghost particles
      LocalGhosts->CopyHydroDataToGhosts(simbox, sph);

      // Calculate gravitational forces from other distant MPI nodes.
      // Also determines particles that must be exported to other nodes
      // if too close to the domain boundaries
#ifdef MPI_PARALLEL
      if (sph->self_gravity == 1) {
        sphneib->UpdateGravityExportList(rank, sph->Nhydro, sph->Ntot, partdata, sph, nbody, simbox);
      }
      else {
        sphneib->UpdateHydroExportList(rank, sph->Nhydro, sph->Ntot, partdata, sph, nbody, simbox);
      }

      // If active particles need forces from other domains, export particles
      mpicontrol->ExportParticlesBeforeForceLoop(sph);
#endif


      // Calculate SPH gravity and hydro forces, depending on which are activated
      if (sph->hydro_forces == 1 && sph->self_gravity == 1) {
        sphneib->UpdateAllSphForces(sph->Nhydro, sph->Ntot, partdata, sph, nbody, simbox, ewald);
      }
      else if (sph->self_gravity == 1) {
        sphneib->UpdateAllSphGravForces(sph->Nhydro, sph->Ntot, partdata, sph, nbody, simbox, ewald);
      }
      else if (sph->hydro_forces == 1) {
        sphneib->UpdateAllSphHydroForces(sph->Nhydro, sph->Ntot, partdata, sph, nbody, simbox);
      }

      // Add external potential for all active SPH particles
      for (i=0; i<sph->Nhydro; i++) {
        SphParticle<ndim>& part = sph->GetSphParticlePointer(i);
        if (part.active) {
          sph->extpot->AddExternalPotential(part.r, part.v, part.a, adot, part.gpot);
        }
      }

      // Checking if acceleration or other values are invalid
      /*for (i=0; i<sph->Nhydro; i++) {
        SphParticle<ndim>& part = sph->GetSphParticlePointer(i);
        if (part.active) {
          for (k=0; k<ndim; k++) assert(part.r[k] == part.r[k]);
          for (k=0; k<ndim; k++) assert(part.v[k] == part.v[k]);
          for (k=0; k<ndim; k++) assert(part.a[k] == part.a[k]);
          assert(part.gpot == part.gpot);
        }
	}*/

#if defined MPI_PARALLEL
      mpicontrol->GetExportedParticlesAccelerations(sph);
#endif

      // Compute the dust forces if present.
      if (sphdust != NULL){
        // Copy properties from original particles to ghost particles
        LocalGhosts->CopyHydroDataToGhosts(simbox, sph);
#ifdef MPI_PARALLEL
        MpiGhosts->CopyHydroDataToGhosts(simbox, sph);
#endif
        sphdust->UpdateAllDragForces(sph->Nhydro, sph->Ntot, partdata) ;
      }


    // Zero all active flags once accelerations have been computed
    for (i=0; i<sph->Nhydro; i++) {
      SphParticle<ndim>& part = sph->GetSphParticlePointer(i);
      part.active = false;
    }

    // Check if all neighbouring timesteps are acceptable.  If not, then set any
    // invalid particles to active to recompute forces immediately.
    if (Nlevels > 1) {
      activecount = sphint->CheckTimesteps(level_diff_max, level_step, n, sph->Nhydro, partdata);
    }
    else {
      activecount = 0;
    }

#if defined MPI_PARALLEL
    MPI_Allreduce(MPI_IN_PLACE, &activecount, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif


  } while (activecount > 0);
  //-----------------------------------------------------------------------------------------------


  /* Check that we have sensible smoothing lengths */
  if (periodicBoundaries) {
    double hmax = sphneib->GetMaximumSmoothingLength() ;
    hmax *= sph->kernp->kernrange ;
    for (i=0; i < ndim; i++)
      if (simbox.boxhalf[i] < 2*hmax){
        string message = "Error: Smoothing length too large, self-interaction will occur" ;
    	ExceptionHandler::getIstance().raise(message);
      }
  }

  // Iterate for P(EC)^n schemes for N-body particles
  //-----------------------------------------------------------------------------------------------
  for (it=0; it<nbody->Npec; it++) {

    // Zero all acceleration terms
    for (i=0; i<nbody->Nnbody; i++) {
      if (nbody->nbodydata[i]->active) {
        for (k=0; k<ndim; k++) nbody->nbodydata[i]->a[k]     = (FLOAT) 0.0;
        for (k=0; k<ndim; k++) nbody->nbodydata[i]->adot[k]  = (FLOAT) 0.0;
        for (k=0; k<ndim; k++) nbody->nbodydata[i]->a2dot[k] = (FLOAT) 0.0;
        for (k=0; k<ndim; k++) nbody->nbodydata[i]->a3dot[k] = (FLOAT) 0.0;
        nbody->nbodydata[i]->gpot = (FLOAT) 0.0;
        nbody->nbodydata[i]->gpe = (FLOAT) 0.0;
      }
    }
    if (sink_particles == 1) {
      for (i=0; i<sinks->Nsink; i++) {
        if (sinks->sink[i].star->active) {
          for (k=0; k<ndim; k++) sinks->sink[i].fhydro[k] = (FLOAT) 0.0;
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

    for (i=0; i<nbody->Nnbody; i++) {
      if (nbody->nbodydata[i]->active) {
        nbody->extpot->AddExternalPotential(nbody->nbodydata[i]->r, nbody->nbodydata[i]->v,
                                            nbody->nbodydata[i]->a, nbody->nbodydata[i]->adot,
                                            nbody->nbodydata[i]->gpot);
      }
    }

    // Calculate correction step for all stars at end of step, except the
    // final iteration (since correction is computed in EndStep also).
    //if (it < nbody->Npec - 1)
    nbody->CorrectionTerms(n,nbody->Nnbody,t,timestep,nbody->nbodydata);

  }
  //-----------------------------------------------------------------------------------------------


  rebuild_tree = false;
  recomputeRadiation = false;


  // End-step terms for all SPH particles
  if (sph->Nhydro > 0) {
    uint->EndTimestep(n, sph->Nhydro, (FLOAT) t, (FLOAT) timestep, partdata);
    sphint->EndTimestep(n, sph->Nhydro, (FLOAT) t, (FLOAT) timestep, partdata);
  }

  // End-step terms for all star particles
  if (nbody->Nstar > 0) nbody->EndTimestep(n, nbody->Nnbody, t, timestep, nbody->nbodydata);

  // Search for new sink particles (if activated) and accrete to existing sinks
  if (sink_particles == 1) {
    if (sinks->create_sinks == 1 && (rebuild_tree || Nfullsteps%ntreebuildstep == 0)) {
      sinks->SearchForNewSinkParticles(n, t ,sph, nbody);
    }
    if (sinks->Nsink > 0) {
      sph->mmean = (FLOAT) 0.0;
      for (i=0; i<sph->Nhydro; i++) sph->mmean += sph->GetSphParticlePointer(i).m;
      sph->mmean /= (FLOAT) sph->Nhydro;
      sph->hmin_sink = big_number;
      for (i=0; i<sinks->Nsink; i++) {
        sph->hmin_sink = min(sph->hmin_sink, (FLOAT) sinks->sink[i].star->h);
      }

      sinks->AccreteMassToSinks(n, timestep, partdata, sph, nbody);
      nbody->UpdateStellarProperties();
      if (extra_sink_output) WriteExtraSinkOutput();
    }
    // If we will output a snapshot (regular or for restarts), then delete all accreted particles
    if ((t >= tsnapnext && sinks->Nsink > 0) || n == nresync || kill_simulation ||
         timing->WallClockTime() - timing->tstart_wall > (FLOAT) 0.99*tmax_wallclock) {
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
  DOUBLE dt_min = big_number_dp;           // Local copy of minimum timestep

  debug2("[SphSimulation::ComputeGlobalTimestep]");
  timing->StartTimingSection("GLOBAL_TIMESTEPS");


  // Only update timestep when all particles are synced at end of last step.
  //-----------------------------------------------------------------------------------------------
  if (n == nresync) {

    n            = 0;
    level_max    = 0;
    level_step   = level_max + integration_step - 1;
    nresync      = integration_step;
    dt_min_nbody = big_number_dp;
    dt_min_hydro   = big_number_dp;

    // Find minimum timestep from all SPH particles
    //---------------------------------------------------------------------------------------------
#pragma omp parallel default(none) shared(dt_min)
    {
      DOUBLE dt_part;                           // Timestep of current particle
      DOUBLE dt = big_number_dp;                // Particle timestep
      DOUBLE dt_nbody = big_number_dp;          // Aux. minimum N-body timestep
      DOUBLE dt_sph = big_number_dp;            // Aux. minimum SPH timestep

      dt       = big_number_dp;
      dt_nbody = big_number_dp;
      dt_sph   = big_number_dp;

#pragma omp for
      for (int i=0; i<sph->Nhydro; i++) {
        SphParticle<ndim>& part = sph->GetSphParticlePointer(i);
        part.level     = 0;
        part.levelneib = 0;
        part.nstep     = pow(2,level_step - part.level);
        part.nlast     = n;
        part.tlast     = t;
        dt_part        = sphint->Timestep(part, sph);
        dt             = min(dt, dt_part);
        dt_sph         = min(dt_sph, dt_part);
      }

      // Now compute minimum timestep due to stars/systems
#pragma omp for
      for (int i=0; i<nbody->Nnbody; i++) {
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
        if (dt_sph < dt_min_hydro) dt_min_hydro = dt_sph;
        if (dt_nbody < dt_min_nbody) dt_min_nbody = dt_nbody;
      }

    }
    //---------------------------------------------------------------------------------------------

#ifdef MPI_PARALLEL
    dt = dt_min;
    MPI_Allreduce(&dt, &dt_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
#endif
    timestep = dt_min;

    // Set minimum timestep for all SPH and N-body particles
    for (int i=0; i<nbody->Nnbody; i++) nbody->nbodydata[i]->dt = timestep;

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
  int i;                                       // Particle counter
  int istep;                                   // Aux. variable for changing steps
  int last_level;                              // Previous timestep level
  int level;                                   // Particle timestep level
  int level_max_aux;                           // Aux. maximum level variable
  int level_max_nbody = 0;                     // level_max for star particles only
  int level_max_old;                           // Old level_max
  int level_max_sph = 0;                       // level_max for SPH particles only
  int level_min_sph = 9999999;                 // level_min for SPH particles
  int level_nbody;                             // local thread var. for N-body level
  int level_sph;                               // local thread var. for SPH level
  int nfactor;                                 // Increase/decrease factor of n
  int nstep;                                   // Particle integer step-size
  DOUBLE dt;                                   // Aux. timestep variable
  DOUBLE dt_part;                              // Particle timestep
  DOUBLE dt_min = big_number_dp;               // Minimum timestep
  DOUBLE dt_min_aux;                           // Aux. minimum timestep variable
  DOUBLE dt_nbody;                             // Aux. minimum N-body timestep
  DOUBLE dt_sph;                               // Aux. minimum SPH timestep

  debug2("[SphSimulation::ComputeBlockTimesteps]");
  timing->StartTimingSection("BLOCK_TIMESTEPS");


  dt_min_nbody = big_number_dp;
  dt_min_hydro = big_number_dp;


  // Synchronise all timesteps and reconstruct block timestep structure.
  //===============================================================================================
  if (n == nresync) {

    n = 0;
    timestep = big_number_dp;

#pragma omp parallel default(none) private(dt,dt_min_aux,dt_nbody,dt_sph,i)
    {
      // Initialise all timestep and min/max variables
      dt_min_aux = big_number_dp;
      dt_sph     = big_number_dp;
      dt_nbody   = big_number_dp;

      // Find minimum timestep from all SPH particles
#pragma omp for
      for (i=0; i<sph->Nhydro; i++) {
        SphParticle<ndim>& part = sph->GetSphParticlePointer(i);
        if (part.flags.is_dead()) continue;
        dt          = sphint->Timestep(part, sph);
        dt_min_aux  = min(dt_min_aux, dt);
        dt_sph      = min(dt_sph, dt);
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
        timestep     = min(timestep,dt_min_aux);
        dt_min_hydro = min(dt_min_hydro,dt_sph);
        dt_min_nbody = min(dt_min_nbody,dt_nbody);
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

    // Calculate the maximum level occupied by all SPH particles
    level_max_sph   = min(ComputeTimestepLevel(dt_min_hydro, dt_max), level_max);
    level_max_nbody = min(ComputeTimestepLevel(dt_min_nbody, dt_max), level_max);

    // Populate timestep levels with N-body particles.
    // Ensures that N-body particles occupy levels lower than all SPH particles
    for (i=0; i<nbody->Nnbody; i++) {
      dt = nbody->nbodydata[i]->dt;
      level = min(ComputeTimestepLevel(dt, dt_max), level_max);
      nbody->nbodydata[i]->level = max(level, level_max_sph);
      nbody->nbodydata[i]->nlast = n;
      nbody->nbodydata[i]->nstep = pow(2, level_step - nbody->nbodydata[i]->level);
      nbody->nbodydata[i]->tlast = t;
    }

    // If particles are sink neighbours, set to same timesteps as sinks
    for (i=0; i<sph->Nhydro; i++) {
      SphParticle<ndim>& part = sph->GetSphParticlePointer(i);
      if (part.sinkid != -1) {
        if (sinks->sink[part.sinkid].star->level - part.level > level_diff_max) {
          part.level     = sinks->sink[part.sinkid].star->level - level_diff_max;
          part.levelneib = sinks->sink[part.sinkid].star->level;
          level_max_sph  = max(level_max_sph, part.level);
        }
      }
    }

    // If enforcing a single SPH timestep, set it here.
    // Otherwise, populate the timestep levels with SPH particles.
    if (sph_single_timestep == 1) {
      for (i=0; i<sph->Nhydro; i++) {
        SphParticle<ndim>& part = sph->GetSphParticlePointer(i);
        if (part.flags.is_dead()) continue;
        part.level     = level_max_sph;
        part.levelneib = level_max_sph;
        part.nlast     = n;
        part.tlast     = t;
        part.nstep     = pow(2, level_step - part.level);
      }
      level_min_sph = level_max_sph;
    }
    else {
      for (i=0; i<sph->Nhydro; i++) {
        SphParticle<ndim>& part = sph->GetSphParticlePointer(i);
        if (part.flags.is_dead()) continue;
        dt             = sphint->Timestep(part, sph);
        level          = min(ComputeTimestepLevel(dt, dt_max), level_max);
        part.level     = level;
        part.levelneib = level;
        part.nlast     = n;
        part.tlast     = t;
        part.nstep     = pow(2, level_step - part.level);
        level_min_sph  = min(level_min_sph, part.level);
      }
    }


    nresync = pow(2, level_step);
    assert(nresync > 0);
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


#pragma omp parallel default(shared) private(dt,dt_nbody,dt_sph,i)\
     private(istep,last_level,level,level_max_aux,level_nbody,level_sph,nstep,nfactor)
    //shared(dt_min,imin,level_max_nbody,level_max_sph,level_min_sph)
    //shared(cout)
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
        if (part.flags.is_dead()) continue;

        // SPH particles whose timestep has been artificially reduced by Saitoh & Makino scheme.
        if (part.nlast == n && part.nstep != pow(2,level_step - part.level)) {
          dt             = sphint->Timestep(part, sph);
          level          = max(ComputeTimestepLevel(dt, dt_max), part.levelneib - level_diff_max);
          part.level     = max(part.level, level);
          part.levelneib = part.level;
          part.nlast     = n;
          part.tlast     = t;
          part.nstep     = pow(2, level_step - part.level);
        }
        // SPH particles that have naturally reached the end of their step
        else if (part.nlast == n) {
          nstep      = part.nstep;
          last_level = part.level;

          // Compute new timestep value and level number
          dt    = sphint->Timestep(part, sph);
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

          part.levelneib = level;
          part.nlast     = n;
          part.tlast     = t;
          part.nstep     = pow(2, level_step - part.level);
        }

        // Find maximum level of all SPH particles
        level_sph     = max(level_sph, part.level);
        level_max_aux = max(level_max_aux, part.level);

        dt_sph = min(dt_sph, dt);
      }
      //-------------------------------------------------------------------------------------------


#pragma omp critical
      {
        dt_min        = min(dt_min, dt_sph);
        dt_min_hydro  = min(dt_min_hydro, dt_sph);
        level_max     = max(level_max, level_max_aux);
        level_max_sph = max(level_max_sph, level_sph);
      }
#pragma omp barrier

#if defined MPI_PARALLEL
#pragma omp master
      {
        level = level_max_sph;
        MPI_Allreduce(&level, &level_max_sph, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
      }
#pragma omp barrier
#endif

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
          level = max(ComputeTimestepLevel(dt, dt_max), level_max_sph);

          // Move up one level (if levels are correctly synchronised) or
          // down several levels if required
          if (level < last_level && level > level_max_sph && last_level > 1 && n%(2*nstep) == 0) {
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
      for (i=0; i<sph->Nhydro; i++) {
        SphParticle<ndim>& part = sph->GetSphParticlePointer(i);
        if (part.flags.is_dead() || part.nlast != n) continue;
        if (part.sinkid != -1) {
          if (sinks->sink[part.sinkid].star->level - part.level > level_diff_max) {
            part.level = sinks->sink[part.sinkid].star->level - level_diff_max;
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
    assert(level_max_sph >= 0);
#endif

    assert(!(isnan(dt_min)) && !(isinf(dt_min)));
    assert(!(isnan(dt_max)) && !(isinf(dt_max)));

    // Set fixed SPH timestep level here in case maximum has changed
    if (sph_single_timestep == 1) {
      for (i=0; i<sph->Nhydro; i++) {
        SphParticle<ndim>& part = sph->GetSphParticlePointer(i);
        if (part.flags.is_dead()) continue;
        if (part.nlast == n) part.level = level_max_sph;
      }
    }

    // Update all timestep variables if we have removed or added any levels
    /*if (level_max != level_max_old) {

      // Increase maximum timestep level if correctly synchronised
      istep = pow(2, level_step - level_max_old + 1);
      if (level_max <= level_max_old - 1 && level_max_old > 1 && n%istep == 0) {
        level_max = level_max_old - 1;
      }
      else if (level_max < level_max_old) {
        level_max = level_max_old;
      }
    }*/

    istep = pow(2, level_step - level_max_old + 1);

    // Adjust integer time if levels are added or removed
    if (level_max > level_max_old) {
      nfactor = pow(2, level_max - level_max_old);
      n *= nfactor;
      for (i=0; i<sph->Nhydro; i++) {
        SphParticle<ndim>& part = sph->GetSphParticlePointer(i);
        if (part.flags.is_dead()) continue;
        part.nstep *= nfactor;
        part.nlast *= nfactor;
      }
      for (i=0; i<nbody->Nnbody; i++) nbody->nbodydata[i]->nstep *= nfactor;
      for (i=0; i<nbody->Nnbody; i++) nbody->nbodydata[i]->nlast *= nfactor;
    }
    else if (level_max <= level_max_old - 1 && level_max_old > 1 && n%istep == 0) {
      level_max = level_max_old - 1;

      nfactor = pow(2, level_max_old - level_max);
      assert(n%nfactor == 0);
      n /= nfactor;
      for (i=0; i<sph->Nhydro; i++) {
        SphParticle<ndim>& part = sph->GetSphParticlePointer(i);
        if (part.flags.is_dead()) continue;
        part.nlast /= nfactor;
        part.nstep /= nfactor;
      }
      for (i=0; i<nbody->Nnbody; i++) nbody->nbodydata[i]->nlast /= nfactor;
      for (i=0; i<nbody->Nnbody; i++) nbody->nbodydata[i]->nstep /= nfactor;
    }
    else {
      level_max = level_max_old;
    }

    level_step = level_max + integration_step - 1;
    nresync    = pow(2, level_step);
    timestep   = dt_max / (DOUBLE) nresync;

    // Update values of nstep for both SPH and star particles
    for (i=0; i<sph->Nhydro; i++) {
      SphParticle<ndim>& part = sph->GetSphParticlePointer(i);
      if (part.flags.is_dead()) continue;
      if (part.nlast == n) part.nstep = pow(2, level_step - part.level);
    }
    for (i=0; i<nbody->Nnbody; i++) {
      if (nbody->nbodydata[i]->nlast == n) {
        nbody->nbodydata[i]->nstep = pow(2, level_step - nbody->nbodydata[i]->level);
      }
    }

    assert(level_max >= level_max_old - 1);

  }
  //===============================================================================================


  // Various asserts for debugging
  assert(timestep >= 0.0 && !(isinf(timestep)) && !(isnan(timestep)));
  assert(dt_max > 0.0 && !(isinf(dt_max)) && !(isnan(dt_max)));
  assert(level_step == level_max + integration_step - 1);
  assert(level_max_sph <= level_max);
  assert(level_max_nbody <= level_max);
  assert(n <= nresync);
  for (i=0; i<sph->Nhydro; i++) {
    SphParticle<ndim>& part = sph->GetSphParticlePointer(i);
    if (part.flags.is_dead()) continue;
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
    cout << "Timestep fallen to zero : " << timestep << "    dtmax: " << dt_max
         << "    nresync " << nresync << endl;
    ExceptionHandler::getIstance().raise("Error : timestep fallen to zero");
  }

  timing->EndTimingSection("BLOCK_TIMESTEPS");

  return;


  // Some validations
  //-----------------------------------------------------------------------------------------------
  /*int *ninlevel;
  int Nactive=0;
  ninlevel = new int[level_max+1];
  //SphParticle<ndim>& part = sph->GetSphParticlePointer(imin);
  cout << "-----------------------------------------------------" << endl;
  cout << "Checking timesteps : " << level_max << "   " << level_max_sph << "    "
       << level_max_nbody << "    " << level_step << "   " << level_max_old << endl;
  cout << "n : " << n << endl;
  cout << "dt_min_hydro : " << dt_min_hydro << "    dt_min_nbody : " << dt_min_nbody
       << "    timestep : " << timestep << endl;
  cout << "hmin_sink : " << sph->hmin_sink << endl;
  //cout << "imin : " << imin << "    " << part.dt << "     " << part.m/sph->mmean << "    "
  //     << part.h << "    " << "    " << part.sound << "     " << part.div_v << "     "
  //     << part.h/(part.sound + part.h*fabs(part.div_v)) << "     "
  //     << sqrt(part.h/sqrt(DotProduct(part.a,part.a,ndim)))
  //     << endl;
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
    ExceptionHandler::getIstance().raise("Error : Timestep fallen to zero");
  }
  cin >> i;*/


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
  for (s=0; s<sinks->Nsink; s++) {

    SinkParticle<ndim> &sink = sinks->sink[s];
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



//=================================================================================================
//  SphSimulation::RegulariseParticleDistribution
/// ...
//=================================================================================================
template <int ndim>
void SphSimulation<ndim>::RegulariseParticleDistribution
 (const int Nreg)                                  ///< [in] No. of regularisation steps
{
  const FLOAT alphaReg = 0.2;                      // Particle displacement magnitude
  FLOAT *rreg = new FLOAT[ndim*sph->Nhydromax];    // Arrya of particle positions
  SphParticle<ndim> *partdata = sph->GetSphParticleArray();


  //===============================================================================================
  for (int ireg=0; ireg<Nreg; ireg++) {

    // Buid/re-build tree, create ghosts and update particle properties
    rebuild_tree = true;
    sphneib->BuildTree(rebuild_tree, 0, ntreebuildstep, ntreestockstep,
                       sph->Ntot, sph->Nhydromax, timestep, partdata, sph);
    sphneib->SearchBoundaryGhostParticles((FLOAT) 0.0, simbox, sph);
    sphneib->BuildGhostTree(true, 0, ntreebuildstep, ntreestockstep, sph->Ntot,
                            sph->Nhydromax, timestep, partdata, sph);
    sphneib->UpdateAllSphProperties(sph->Nhydro, sph->Ntot, partdata, sph, nbody);


    //=============================================================================================
#pragma omp parallel default(none) shared(partdata, rreg)
    {
      int k;
      FLOAT dr[ndim];
      FLOAT drsqd;
      int *neiblist = new int[sph->Nhydromax];


      //-------------------------------------------------------------------------------------------
#pragma omp for
      for (int i=0; i<sph->Nhydro; i++) {
        SphParticle<ndim> &part = sph->GetSphParticlePointer(i);
        FLOAT invhsqd = part.invh*part.invh;
        for (k=0; k<ndim; k++) rreg[ndim*i + k] = (FLOAT) 0.0;

        // Find list of gather neighbours
        int Nneib = sphneib->GetGatherNeighbourList(part.r, sph->kernrange*part.h, partdata,
                                                    sph->Ntot, sph->Nhydromax, neiblist);

        // Loop over all neighbours and calculate position correction for regularisation
        //-----------------------------------------------------------------------------------------
        for (int jj=0; jj<Nneib; jj++) {
          int j = neiblist[jj];
          SphParticle<ndim> &neibpart = sph->GetSphParticlePointer(j);

          for (k=0; k<ndim; k++) dr[k] = neibpart.r[k] - part.r[k];
          drsqd = DotProduct(dr, dr, ndim);
          if (drsqd >= part.hrangesqd) continue;
          for (k=0; k<ndim; k++) rreg[ndim*i + k] += dr[k]*sph->kernp->w0_s2(drsqd*invhsqd);

        }
        //-----------------------------------------------------------------------------------------

      }
      //-------------------------------------------------------------------------------------------


      // Apply all regularisation corrections to particle positions
      //-------------------------------------------------------------------------------------------
#pragma omp for
      for (int i=0; i<sph->Nhydro; i++) {
        SphParticle<ndim> &part = sph->GetSphParticlePointer(i);
        for (k=0; k<ndim; k++) part.r[k] -= alphaReg*rreg[ndim*i + k];
      }

      delete[] neiblist;

    }
    //=============================================================================================

    // Check that new positions don't fall outside the domain box
    sphint->CheckBoundaries(simbox, sph);

  }
  //================================================================================================

  delete[] rreg;

  return;
}



//=================================================================================================
//  SphSimulation::SmoothParticleQuantity
/// ...
//=================================================================================================
template <int ndim>
void SphSimulation<ndim>::SmoothParticleQuantity
 (const int Npart,
  FLOAT *values)
{

  return;
}
