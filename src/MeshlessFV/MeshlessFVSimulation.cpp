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
#include "Sph.h"
#include "RiemannSolver.h"
#include "Ghosts.h"
#include "Sinks.h"
using namespace std;



// Create template class instances of the main MeshlessFVSimulation object for
// each dimension used (1, 2 and 3)
template class MeshlessFVSimulation<1>;
template class MeshlessFVSimulation<2>;
template class MeshlessFVSimulation<3>;



//=================================================================================================
//  MeshlessFVSimulation::ProcessParameters
/// Process all the options chosen in the parameters file, setting various
/// simulation variables and creating important simulation objects.
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


  // Sanity check for valid dimensionality
  if (ndim < 1 || ndim > 3) {
    std::ostringstream message;
    message << "Invalid dimensionality chosen : ndim = " << ndim;
    ExceptionHandler::getIstance().raise(message.str());
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


  // Create Meshless Finite-Volume object depending on choice of kernel
  //===============================================================================================
  if (sim == "meshlessfv" || sim == "mfvmuscl") {
    if (intparams["tabulated_kernel"] == 1) {
      mfv = new MfvMuscl<ndim, TabulatedKernel>
        (intparams["hydro_forces"], intparams["self_gravity"], floatparams["accel_mult"],
         floatparams["courant_mult"], floatparams["h_fac"], floatparams["h_converge"],
         floatparams["gamma_eos"], stringparams["gas_eos"],
         KernelName, sizeof(MeshlessFVParticle<ndim>));
    }
    else if (intparams["tabulated_kernel"] == 0) {
      if (KernelName == "m4") {
        mfv = new MfvMuscl<ndim, M4Kernel>
          (intparams["hydro_forces"], intparams["self_gravity"], floatparams["accel_mult"],
           floatparams["courant_mult"], floatparams["h_fac"], floatparams["h_converge"],
           floatparams["gamma_eos"], stringparams["gas_eos"],
           KernelName, sizeof(MeshlessFVParticle<ndim>));
      }
      else if (KernelName == "quintic") {
        mfv = new MfvMuscl<ndim, QuinticKernel>
          (intparams["hydro_forces"], intparams["self_gravity"], floatparams["accel_mult"],
           floatparams["courant_mult"], floatparams["h_fac"], floatparams["h_converge"],
           floatparams["gamma_eos"], stringparams["gas_eos"],
           KernelName, sizeof(MeshlessFVParticle<ndim>));
      }
      else if (KernelName == "gaussian") {
        mfv = new MfvMuscl<ndim, GaussianKernel>
          (intparams["hydro_forces"], intparams["self_gravity"], floatparams["accel_mult"],
           floatparams["courant_mult"], floatparams["h_fac"], floatparams["h_converge"],
           floatparams["gamma_eos"], stringparams["gas_eos"],
           KernelName, sizeof(MeshlessFVParticle<ndim>));
      }
      else {
        string message = "Unrecognised parameter : kernel = " + simparams->stringparams["kernel"];
        ExceptionHandler::getIstance().raise(message);
      }
    }
    else {
      string message = "Invalid option for the tabulated_kernel parameter: " +
        stringparams["tabulated_kernel"];
      ExceptionHandler::getIstance().raise(message);
    }
  }
  //===============================================================================================
  else if (sim == "mfvrk") {
    if (intparams["tabulated_kernel"] == 1) {
      mfv = new MfvRungeKutta<ndim, TabulatedKernel>
        (intparams["hydro_forces"], intparams["self_gravity"], floatparams["accel_mult"],
         floatparams["courant_mult"], floatparams["h_fac"], floatparams["h_converge"],
         floatparams["gamma_eos"], stringparams["gas_eos"],
         KernelName, sizeof(MeshlessFVParticle<ndim>));
    }
    else if (intparams["tabulated_kernel"] == 0) {
      if (KernelName == "m4") {
        mfv = new MfvRungeKutta<ndim, M4Kernel>
          (intparams["hydro_forces"], intparams["self_gravity"], floatparams["accel_mult"],
           floatparams["courant_mult"], floatparams["h_fac"], floatparams["h_converge"],
           floatparams["gamma_eos"], stringparams["gas_eos"],
           KernelName, sizeof(MeshlessFVParticle<ndim>));
      }
      else if (KernelName == "quintic") {
        mfv = new MfvRungeKutta<ndim, QuinticKernel>
          (intparams["hydro_forces"], intparams["self_gravity"], floatparams["accel_mult"],
           floatparams["courant_mult"], floatparams["h_fac"], floatparams["h_converge"],
           floatparams["gamma_eos"], stringparams["gas_eos"],
           KernelName, sizeof(MeshlessFVParticle<ndim>));
      }
      else if (KernelName == "gaussian") {
        mfv = new MfvRungeKutta<ndim, GaussianKernel>
          (intparams["hydro_forces"], intparams["self_gravity"], floatparams["accel_mult"],
           floatparams["courant_mult"], floatparams["h_fac"], floatparams["h_converge"],
           floatparams["gamma_eos"], stringparams["gas_eos"],
           KernelName, sizeof(MeshlessFVParticle<ndim>));
      }
      else {
        string message = "Unrecognised parameter : kernel = " + simparams->stringparams["kernel"];
        ExceptionHandler::getIstance().raise(message);
      }
    }
    else {
      string message = "Invalid option for the tabulated_kernel parameter: " +
        stringparams["tabulated_kernel"];
      ExceptionHandler::getIstance().raise(message);
    }
  }
  //===============================================================================================
  else {
    string message = "Invalid option for the simulation type parameter: " + sim;
    ExceptionHandler::getIstance().raise(message);
  }
  //===============================================================================================
  hydro = mfv;


  // Riemann solver object
  //-----------------------------------------------------------------------------------------------
  string riemann = stringparams["riemann_solver"];
  if (riemann == "exact") {
    mfv->riemann = new ExactRiemannSolver<ndim>(floatparams["gamma_eos"], intparams["zero_mass_flux"]);
  }
  else if (riemann == "hllc") {
    mfv->riemann = new HllcRiemannSolver<ndim>(floatparams["gamma_eos"], intparams["zero_mass_flux"]);
  }
  else {
    string message = "Unrecognised parameter : riemann_solver = " + riemann;
    ExceptionHandler::getIstance().raise(message);
  }


  // Slope limiter
  //-----------------------------------------------------------------------------------------------
  string limiter = stringparams["slope_limiter"];
  if (limiter == "null") {
    mfv->limiter = new NullLimiter<ndim,MeshlessFVParticle>();
  }
  else if (limiter == "zeroslope") {
    mfv->limiter = new ZeroSlopeLimiter<ndim,MeshlessFVParticle>();
  }
  else if (limiter == "balsara2004") {
    mfv->limiter = new Balsara2004Limiter<ndim,MeshlessFVParticle>();
  }
  else if (limiter == "springel2009") {
    mfv->limiter = new Springel2009Limiter<ndim,MeshlessFVParticle>();
  }
  else if (limiter == "tess2011") {
    mfv->limiter = new TESS2011Limiter<ndim,MeshlessFVParticle>();
  }
  else if (limiter == "gizmo") {
    mfv->limiter = new GizmoLimiter<ndim,MeshlessFVParticle>();
  }
  else if (limiter == "minmod") {
    mfv->limiter = new MinModLimiter<ndim,MeshlessFVParticle>();
  }
  else {
    string message = "Unrecognised parameter : slope_limiter = " + limiter;
    ExceptionHandler::getIstance().raise(message);
  }


  // Create neighbour searching object based on chosen method in params file
  //-----------------------------------------------------------------------------------------------
  if (stringparams["neib_search"] == "bruteforce") {
    mfvneib = new MeshlessFVBruteForce<ndim,MeshlessFVParticle>
      (mfv->kernp->kernrange, &simbox, mfv->kernp, timing);
  }
  else if (stringparams["neib_search"] == "kdtree") {
    mfvneib = new MeshlessFVKDTree<ndim,MeshlessFVParticle,KDTreeCell>
     (intparams["Nleafmax"], Nmpi, intparams["pruning_level_min"], intparams["pruning_level_max"],
      floatparams["thetamaxsqd"], hydro->kernp->kernrange, floatparams["macerror"],
      stringparams["gravity_mac"], stringparams["multipole"], &simbox, mfv->kernp, timing);
  }
  else if (stringparams["neib_search"] == "octtree") {
    mfvneib = new MeshlessFVOctTree<ndim,MeshlessFVParticle,OctTreeCell>
     (intparams["Nleafmax"], Nmpi, intparams["pruning_level_min"], intparams["pruning_level_max"],
      floatparams["thetamaxsqd"], hydro->kernp->kernrange, floatparams["macerror"],
      stringparams["gravity_mac"], stringparams["multipole"], &simbox, mfv->kernp, timing);
;
  }
  else {
    string message = "Unrecognised parameter : neib_search = "
      + simparams->stringparams["neib_search"];
    ExceptionHandler::getIstance().raise(message);
  }
  ////mfvneib->kernp = mfv->kernp;
  mfvneib->kernfac = mfv->kernfac;
  ////mfvneib->kernrange = mfv->kernp->kernrange;


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

  // Create ghost particle object
  //-----------------------------------------------------------------------------------------------
  /*if (IsAnyBoundarySpecial(simbox)) {
    LocalGhosts = new PeriodicGhostsSpecific<ndim,MeshlessFVParticle >();
  }
  else {
    LocalGhosts = new NullGhosts<ndim>();
  }*/

  // Process all N-body parameters and set-up main N-body objects
  this->ProcessNbodyParameters();


  // Thermal physics object.  If energy equation is chosen, also initiate
  // the energy integration object.
  //-----------------------------------------------------------------------------------------------
  if (gas_eos == "energy_eqn") {
    mfv->eos = new Adiabatic<ndim>
      (floatparams["temp0"], floatparams["mu_bar"], floatparams["gamma_eos"]);
  }
  else {
    string message = "Unrecognised or invalid parameter : gas_eos = " + gas_eos;
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
  mfv->extpot = extpot;
  nbody->extpot = extpot;


  // Create Ewald periodic gravity object
  /*periodicBoundaries = IsAnyBoundarySpecial(simbox);
  if (periodicBoundaries && intparams["self_gravity"] == 1) {
    ewaldGravity = true;
    ewald = new Ewald<ndim>(simbox,intparams["gr_bhewaldseriesn"],intparams["in"],
                            intparams["nEwaldGrid"],floatparams["ewald_mult"],
                            floatparams["ixmin"],floatparams["ixmax"],timing);
  }*/


  // Set all other SPH parameter variables
  //mfv->Nhydro            = intparams["Nhydro"];
  mfv->Nhydromax       = intparams["Nhydromax"];
  //mfv->create_sinks    = intparams["create_sinks"];
  //mfv->fixed_sink_mass = intparams["fixed_sink_mass"];
  //mfv->msink_fixed     = floatparams["m1"];


  // Set important variables for N-body objects
  //nbody->Nstar          = intparams["Nstar"];
  nbody->Nstarmax       = intparams["Nstarmax"];
  nbody_single_timestep = intparams["nbody_single_timestep"];
  nbodytree.gpehard     = floatparams["gpehard"];
  nbodytree.gpesoft     = floatparams["gpesoft"];
  //nbody->perturbers     = intparams["perturbers"];
  //if (intparams["sub_systems"] == 1) subsystem->perturbers = intparams["perturbers"];


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
  /*if (sim == "sph" || sim == "gradhsph" || sim == "sm2012sph" || sim == "godunov_hydro") {
    sinks.timing    = timing;
    sphint->timing  = timing;
    mfvneib->timing = timing;
    uint->timing    = timing;
    radiation->timing = timing;
  }*/

  sinks = new Sinks<ndim>(mfvneib);


  // Flag that we've processed all parameters already
  ParametersProcessed = true;


  return;
}



//=================================================================================================
//  MeshlessFVSimulation::PostInitialConditionsSetup
/// Call routines for calculating all initial SPH and N-body quantities
/// once initial conditions have been set-up.
//=================================================================================================
template <int ndim>
void MeshlessFVSimulation<ndim>::PostInitialConditionsSetup(void)
{
  int i;                               // Particle counter
  int k;                               // Dimension counter
  MeshlessFVParticle<ndim> *partdata;  // Pointer to main SPH data array

  debug2("[MeshlessFVSimulation::PostInitialConditionsSetup]");

  // Set iorig
  if (rank == 0) {
    for (i=0; i<mfv->Nhydro; i++) mfv->GetMeshlessFVParticlePointer(i).iorig = i;
  }

  // Set pointer to SPH particle data
  partdata = mfv->GetMeshlessFVParticleArray();

  // Set time variables here (for now)
  nresync = 0;   // DAVID : Need to adapt this for block timesteps
  n = 0;
  integration_step = 1;

  // Set initial smoothing lengths and create initial ghost particles
  //-----------------------------------------------------------------------------------------------
  if (mfv->Nhydro > 0) {

    // Set all relevant particle counters
    mfv->Nghost = 0;
    mfv->Nghostmax = mfv->Nhydromax - mfv->Nhydro;
    mfv->Ntot = mfv->Nhydro;
    for (i=0; i<mfv->Nhydro; i++) mfv->GetMeshlessFVParticlePointer(i).active = true;

    // Compute mean mass (used for smooth sink accretion)
    if (!restart) {
      mfv->mmean = 0.0;
      for (i=0; i<mfv->Nhydro; i++) mfv->mmean += mfv->GetMeshlessFVParticlePointer(i).m;
      mfv->mmean /= (FLOAT) mfv->Nhydro;
    }

    // If the smoothing lengths have not been provided beforehand, then
    // calculate the initial values here
    mfvneib->neibcheck = false;
    if (!this->initial_h_provided) {
      mfv->InitialSmoothingLengthGuess();
      mfvneib->BuildTree(rebuild_tree, 0, ntreebuildstep, ntreestockstep, mfv->Ntot,
                         mfv->Nhydromax, timestep, mfv->GetMeshlessFVParticleArray(), mfv);
      mfvneib->UpdateAllProperties(mfv->Nhydro, mfv->Ntot,
                                   mfv->GetMeshlessFVParticleArray(), mfv, nbody);
    }
    else {
      mfvneib->BuildTree(rebuild_tree, 0, ntreebuildstep, ntreestockstep, mfv->Ntot,
                         mfv->Nhydromax, timestep, mfv->GetMeshlessFVParticleArray(), mfv);
    }

    // Search ghost particles
    mfvneib->SearchBoundaryGhostParticles(0.0,simbox,mfv);
    mfvneib->BuildGhostTree(true, 0, ntreebuildstep, ntreestockstep ,mfv->Ntot,
                            mfv->Nhydromax, timestep, mfv->GetMeshlessFVParticleArray(), mfv);


    // Zero accelerations
    for (i=0; i<mfv->Nhydro; i++) mfv->GetMeshlessFVParticlePointer(i).active = true;

    // Update neighbour tree
    rebuild_tree = true;
    mfvneib->BuildTree(rebuild_tree, 0, ntreebuildstep, ntreestockstep, mfv->Ntot,
                       mfv->Nhydromax, timestep, mfv->GetMeshlessFVParticleArray(), mfv);
    level_step = 1;


    // For Eigenvalue MAC, need non-zero values
    for (i=0; i<mfv->Nhydro; i++) mfv->GetMeshlessFVParticlePointer(i).gpot = big_number;

    // Calculate all SPH properties
    mfvneib->UpdateAllProperties(mfv->Nhydro ,mfv->Ntot,
                                 mfv->GetMeshlessFVParticleArray(), mfv, nbody);

    // Search ghost particles
    mfvneib->SearchBoundaryGhostParticles(0.0, simbox, mfv);
    mfvneib->BuildGhostTree(true, 0, ntreebuildstep, ntreestockstep, mfv->Ntot,
                            mfv->Nhydromax, timestep, mfv->GetMeshlessFVParticleArray(), mfv);

    // Update neighbour tree
    rebuild_tree = true;
    mfvneib->BuildTree(rebuild_tree, 0, ntreebuildstep, ntreestockstep, mfv->Ntot,
                       mfv->Nhydromax, timestep, mfv->GetMeshlessFVParticleArray(), mfv);
    mfvneib->SearchBoundaryGhostParticles(0.0,simbox,mfv);
    mfvneib->BuildGhostTree(true, 0, ntreebuildstep, ntreestockstep, mfv->Ntot,
                            mfv->Nhydromax, timestep, mfv->GetMeshlessFVParticleArray(), mfv);
    mfvneib->neibcheck = true;

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
  nbody->LoadStellarPropertiesTable(&simunits);
  nbody->UpdateStellarProperties();


  // Compute all initial SPH force terms
  //-----------------------------------------------------------------------------------------------
  if (mfv->Nhydro > 0) {

    // Zero accelerations (here for now)
    for (i=0; i<mfv->Ntot; i++) {
      MeshlessFVParticle<ndim>& part = mfv->GetMeshlessFVParticlePointer(i);
      part.level  = 0;
      part.nstep  = 0;
      part.nlast  = 0;
      part.tlast  = t;
      part.active = false;
    }
    for (i=0; i<mfv->Nhydro; i++) mfv->GetMeshlessFVParticlePointer(i).active = true;

    // Copy all other data from real SPH particles to ghosts
    mfv->CopyDataToGhosts(simbox, partdata);
    //LocalGhosts->CopyHydroDataToGhosts(simbox,sph);

    mfvneib->BuildTree(true, 0, ntreebuildstep, ntreestockstep, mfv->Ntot,
                       mfv->Nhydromax, timestep, mfv->GetMeshlessFVParticleArray(), mfv);
    mfvneib->SearchBoundaryGhostParticles(0.0, simbox, mfv);
    mfvneib->BuildGhostTree(true, 0, ntreebuildstep, ntreestockstep, mfv->Ntot,
                            mfv->Nhydromax, timestep, mfv->GetMeshlessFVParticleArray(), mfv);


    // Set initial accelerations
    for (i=0; i<mfv->Nhydro; i++) {
      MeshlessFVParticle<ndim>& part = mfv->GetMeshlessFVParticlePointer(i);
      for (k=0; k<ndim; k++) part.r0[k] = part.r[k];
      for (k=0; k<ndim; k++) part.v0[k] = part.v[k];
      for (k=0; k<ndim; k++) part.a0[k] = part.a[k];
      part.active = true;
    }

    mfv->CopyDataToGhosts(simbox, partdata);
#ifdef MPI_PARALLEL
//    MpiGhosts->CopyHydroDataToGhosts(simbox,sph);
#endif

    // Update the primitive vectors for all particles
    for (i=0; i<mfv->Ntot; i++) {
      mfv->ComputeThermalProperties(partdata[i]);
      mfv->UpdatePrimitiveVector(partdata[i]);
      mfv->ConvertPrimitiveToConserved(partdata[i].Wprim, partdata[i].Ucons);
      mfv->ConvertConservedToQ(partdata[i].volume, partdata[i].Ucons, partdata[i].Qcons);
      partdata[i].Utot = partdata[i].u*partdata[i].m;
      for (k=0; k<ndim+2; k++) partdata[i].dQ[k] = (FLOAT) 0.0;
      //for (k=0; k<ndim; k++) partdata[i].v0[k] = partdata[i].v[k];
    }

    mfvneib->UpdateGradientMatrices(mfv->Nhydro, mfv->Ntot, partdata, mfv, nbody);
    mfv->CopyDataToGhosts(simbox, partdata);
    //mfvneib->UpdateGodunovFluxes(mfv->Nhydro, mfv->Ntot, timestep, partdata, mfv, nbody);

    if (mfv->self_gravity == 1) {
      mfvneib->UpdateAllGravForces(mfv->Nhydro, mfv->Ntot, partdata, mfv, nbody);
    }

  }
  //-----------------------------------------------------------------------------------------------


  // Compute initial N-body forces
  //-----------------------------------------------------------------------------------------------
  if (nbody->Nstar > 0) {
    if (mfv->self_gravity == 1 && mfv->Nhydro > 0) {
      mfvneib->UpdateAllStarGasForces(mfv->Nhydro, mfv->Ntot,
                                      mfv->GetMeshlessFVParticleArray(), mfv, nbody);
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

  }
  //-----------------------------------------------------------------------------------------------


  // Call EndTimestep to set all 'beginning-of-step' variables
  mfv->EndTimestep(n, mfv->Nhydro, t, timestep, mfv->GetMeshlessFVParticleArray());
  nbody->EndTimestep(n, nbody->Nstar , t, timestep, nbody->nbodydata);

  this->CalculateDiagnostics();
  this->diag0 = this->diag;
  this->setup = true;


  return;
}



//=================================================================================================
//  MeshlessFVSimulation::ComputeGlobalTimestep
/// Computes global timestep for SPH simulation.  Calculates the minimum
/// timestep for all SPH and N-body particles in the simulation.
//=================================================================================================
template <int ndim>
void MeshlessFVSimulation<ndim>::ComputeGlobalTimestep(void)
{
  int i;                               // Particle counter
  //DOUBLE dt;                           // Particle timestep
  DOUBLE dt_min = big_number_dp;       // Local copy of minimum timestep
  //DOUBLE dt_nbody;                     // Aux. minimum N-body timestep
  //DOUBLE dt_hydro;                     // Aux. minimum SPH timestep

  debug2("[MeshlessFVSimulation::ComputeGlobalTimestep]");
  timing->StartTimingSection("GLOBAL_TIMESTEPS");


  // Only update timestep when all particles are synced at end of last step.
  //-----------------------------------------------------------------------------------------------
  if (n == nresync) {

    n            = 0;
    level_max    = 0;
    level_step   = level_max + integration_step - 1;
    nresync      = integration_step;
    dt_min_hydro = big_number_dp;


    // Find minimum timestep from all hydro particles
    //---------------------------------------------------------------------------------------------
#pragma omp parallel default(none) private(i) shared(dt_min) //,dt_min_nbody,dt_min_hydro)
    {

#pragma omp for
      for (i=0; i<mfv->Nhydro; i++) {
        MeshlessFVParticle<ndim>& part = mfv->GetMeshlessFVParticlePointer(i);
        part.level     = 0;
        part.levelneib = 0;
        part.nstep     = pow(2,level_step - part.level);
        part.nlast     = n;
        part.tlast     = t;
        part.dt        = mfv->Timestep(part);
        dt_min         = min(dt_min, part.dt);
      }

    }
    //---------------------------------------------------------------------------------------------

    timestep = dt_min;

    // Set minimum timestep for all SPH and N-body particles
    for (i=0; i<mfv->Nhydro; i++) mfv->GetMeshlessFVParticlePointer(i).dt = timestep;

  }
  //-----------------------------------------------------------------------------------------------

  timing->EndTimingSection("GLOBAL_TIMESTEPS");

  return;
}



//=================================================================================================
//  MeshlessFVSimulation::ComputeBlockTimesteps
/// Compute timesteps for all particles using hierarchical block timesteps.
//=================================================================================================
template <int ndim>
void MeshlessFVSimulation<ndim>::ComputeBlockTimesteps(void)
{
  int i;                                     // Particle counter
  int istep;                                 // Aux. variable for changing steps
  int last_level;                            // Previous timestep level
  int level;                                 // Particle timestep level
  int level_max_aux;                         // Aux. maximum level variable
  int level_max_nbody = 0;                   // level_max for star particles only
  int level_max_old = level_max;             // Old level_max
  int level_max_hydro = 0;                   // level_max for SPH particles only
  int level_min_hydro = 9999999;             // level_min for SPH particles
  //int level_nbody;                         // local thread var. for N-body level
  int level_hydro;                           // local thread var. for SPH level
  int nfactor;                               // Increase/decrease factor of n
  int nstep;                                 // Particle integer step-size
  DOUBLE dt;                                 // Aux. timestep variable
  DOUBLE dt_min = big_number_dp;             // Minimum timestep
  DOUBLE dt_min_aux;                         // Aux. minimum timestep variable
  DOUBLE dt_nbody;                           // Aux. minimum N-body timestep
  DOUBLE dt_hydro;                           // Aux. minimum SPH timestep

  debug2("[MeshlessFVSimulation::ComputeBlockTimesteps]");
  timing->StartTimingSection("BLOCK_TIMESTEPS");


  dt_min_nbody = big_number_dp;
  dt_min_hydro = big_number_dp;


  // Synchronise all timesteps and reconstruct block timestep structure.
  //===============================================================================================
  if (n == nresync) {

    n = 0;
    timestep = big_number_dp;

#pragma omp parallel default(none) private(dt,dt_min_aux,dt_nbody,dt_hydro,i)
    {
      // Initialise all timestep and min/max variables
      dt_min_aux = big_number_dp;
      dt_hydro   = big_number_dp;
      dt_nbody   = big_number_dp;

      // Find minimum timestep from all hydro particles
#pragma omp for
      for (i=0; i<mfv->Nhydro; i++) {
        MeshlessFVParticle<ndim>& part = mfv->GetMeshlessFVParticlePointer(i);
        if (part.itype == dead) continue;
        dt         = mfv->Timestep(part);
        dt_min_aux = min(dt_min_aux, dt);
        dt_hydro   = min(dt_hydro, dt);
        part.dt    = dt;
      }

      // Now compute minimum timestep due to stars/systems
/*#pragma omp for
      for (i=0; i<nbody->Nnbody; i++) {
        dt         = nbody->Timestep(nbody->nbodydata[i]);
        dt_min_aux = min(dt_min_aux, dt);
        dt_nbody   = min(dt_nbody, dt);
        nbody->nbodydata[i]->dt = dt;
      }*/

#pragma omp critical
      {
        timestep     = min(timestep, dt_min_aux);
        dt_min_hydro = min(dt_min_hydro, dt_hydro);
        dt_min_nbody = min(dt_min_nbody, dt_nbody);
      }
#pragma omp barrier
    }


    // For MPI, determine the global minimum timestep over all processors
/*#ifdef MPI_PARALLEL
    dt = timestep;
    MPI_Allreduce(&dt,&timestep,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
    dt = dt_min_hydro;
    MPI_Allreduce(&dt,&dt_min_hydro,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
    dt = dt_min_nbody;
    MPI_Allreduce(&dt,&dt_min_nbody,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
#endif*/


    // Calculate new block timestep levels
    level_max  = Nlevels - 1;
    level_step = level_max + integration_step - 1;
    dt_max     = timestep*powf(2.0,level_max);

    // Calculate the maximum level occupied by all SPH particles
    level_max_hydro   = min((int) (invlogetwo*log(dt_max/dt_min_hydro)) + 1, level_max);
    /*level_max_nbody = min((int) (invlogetwo*log(dt_max/dt_min_nbody)) + 1, level_max);

    // Populate timestep levels with N-body particles.
    // Ensures that N-body particles occupy levels lower than all SPH particles
    for (i=0; i<nbody->Nnbody; i++) {
      dt = nbody->nbodydata[i]->dt;
      level = min((int) (invlogetwo*log(dt_max/dt)) + 1, level_max);
      level = max(level,0);
      nbody->nbodydata[i]->level = max(level,level_max_hydro);
      nbody->nbodydata[i]->nlast = n;
      nbody->nbodydata[i]->nstep = pow(2,level_step - nbody->nbodydata[i]->level);
      nbody->nbodydata[i]->tlast = t;
    }

    // If particles are sink neighbours, set to same timesteps as sinks
    for (i=0; i<mfv->Nhydro; i++) {
      SphParticle<ndim>& part = mfv->GetSphParticlePointer(i);
      if (part.sinkid != -1) {
        if (sinks.sink[part.sinkid].star->level - part.level > level_diff_max) {
          part.level     = sinks.sink[part.sinkid].star->level - level_diff_max;
          part.levelneib = sinks.sink[part.sinkid].star->level;
          level_max_hydro  = max(level_max_hydro,part.level);
        }
      }
    }*/

    // If enforcing a single SPH timestep, set it here.
    // Otherwise, populate the timestep levels with SPH particles.
    if (sph_single_timestep == 1) {
      for (i=0; i<mfv->Nhydro; i++) {
        MeshlessFVParticle<ndim>& part = mfv->GetMeshlessFVParticlePointer(i);
        if (part.itype == dead) continue;
        part.active    = true;
        part.level     = level_max_hydro;
        part.levelneib = level_max_hydro;
        part.nlast     = n;
        part.tlast     = t;
        part.nstep     = pow(2,level_step - part.level);
      }
      level_min_hydro = level_max_hydro;
    }
    else {
      for (i=0; i<mfv->Nhydro; i++) {
        MeshlessFVParticle<ndim>& part = mfv->GetMeshlessFVParticlePointer(i);
        if (part.itype == dead) continue;
        dt              = part.dt;
        level           = min((int) (invlogetwo*log(dt_max/dt)) + 1, level_max);
        level           = max(level, 0);
        part.active     = true;
        part.level      = level;
        part.levelneib  = level;
        part.nlast      = n;
        part.tlast      = t;
        part.nstep      = pow(2,level_step - part.level);
        level_min_hydro = min(level_min_hydro, part.level);
      }
    }

    nresync = pow(2,level_step);
    assert(nresync > 0);
    timestep = dt_max / (DOUBLE) nresync;

  }
  // If not resynchronising, check if any SPH/N-body particles need to move
  // up or down timestep levels.
  //===============================================================================================
  else {

    level_max       = 0;
    level_max_nbody = 0;
    level_max_hydro = 0;


#pragma omp parallel default(none) private(dt,dt_nbody,dt_hydro,i) \
  private(istep,last_level,level,level_max_aux,level_hydro,nstep,nfactor) \
  shared(dt_min,level_max_nbody,level_max_hydro,level_min_hydro)
    {
      dt_hydro      = big_number_dp;
      dt_nbody      = big_number_dp;
      level_max_aux = 0;
      //level_nbody   = 0;
      level_hydro   = 0;


      // Find all hydro particles at the beginning of a new timestep
      //-------------------------------------------------------------------------------------------
#pragma omp for
      for (i=0; i<mfv->Nhydro; i++) {
        MeshlessFVParticle<ndim>& part = mfv->GetMeshlessFVParticlePointer(i);
        if (part.itype == dead) continue;
        part.active = false;

        // SPH particles whose timestep has been artificially reduced by Saitoh & Makino scheme.
        if (part.nlast == n && part.nstep != pow(2, level_step - part.level)) {
          dt    = mfv->Timestep(part);
          level = max((int) (invlogetwo*log(dt_max/dt)) + 1, 0);
          level = max(level, part.levelneib - level_diff_max);

          part.active    = true;
          part.level     = max(part.level, level);
          part.levelneib = part.level;
          part.dt        = dt;
          part.nlast     = n;
          part.tlast     = t;
          part.nstep     = pow(2, level_step - part.level);
        }
        // Hydro particles that have naturally reached the end of their step
        else if (part.nlast == n) {
          nstep      = part.nstep;
          last_level = part.level;

          // Compute new timestep value and level number
          dt    = mfv->Timestep(part);
          level = max((int) (invlogetwo*log(dt_max/dt)) + 1, 0);
          level = max(level, part.levelneib - level_diff_max);

          // Move up one level (if levels are correctly synchronised) or
          // down several levels if required
          if (level < last_level && last_level > 1 && n%(2*nstep) == 0)
            part.level = last_level - 1;
          else if (level > last_level)
            part.level = level;
          else
            part.level = last_level;

          part.active    = true;
          part.levelneib = level;
          part.dt        = dt;
          part.nlast     = n;
          part.tlast     = t;
          part.nstep     = pow(2u,level_step - part.level);
        }

        // Find maximum level of all SPH particles
        level_hydro   = max(level_hydro,part.level);
        level_max_aux = max(level_max_aux,part.level);

        dt_hydro = min(dt_hydro,part.dt);
      }
      //-------------------------------------------------------------------------------------------


#pragma omp critical
      {
        dt_min          = min(dt_min,dt_hydro);
        dt_min_hydro    = min(dt_min_hydro,dt_hydro);
        level_max       = max(level_max,level_max_aux);
        level_max_hydro = max(level_max_hydro,level_hydro);
      }
#pragma omp barrier


      // Now find all N-body particles at the beginning of a new timestep
      //-------------------------------------------------------------------------------------------
/*#pragma omp for
      for (i=0; i<nbody->Nnbody; i++) {

        // Skip particles that are not at end of step
        if (nbody->nbodydata[i]->nlast == n) {
          nstep = nbody->nbodydata[i]->nstep;
          last_level = nbody->nbodydata[i]->level;

          // Compute new timestep value and level number
          dt    = nbody->Timestep(nbody->nbodydata[i]);
          level = max((int) (invlogetwo*log(dt_max/dt)) + 1, 0);
          level = max(level,level_max_hydro);

          // Move up one level (if levels are correctly synchronised) or
          // down several levels if required
          if (level < last_level && level > level_max_hydro && last_level > 1 && n%(2*nstep) == 0)
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
      }*/
      #pragma omp barrier


      // Correct timestep levels for any particles that have entered a sink
      //-------------------------------------------------------------------------------------------
  /*#pragma omp for
      for (i=0; i<mfv->Nhydro; i++) {
        SphParticle<ndim>& part = mfv->GetSphParticlePointer(i);
        if (part.itype == dead || part.nlast != n) continue;
        if (part.sinkid != -1) {
          if (sinks.sink[part.sinkid].star->level - part.level > level_diff_max) {
            part.level = sinks.sink[part.sinkid].star->level - level_diff_max;
          }
        }
      }*/

    }

    // For MPI, find the global maximum timestep levels for each processor
/*#ifdef MPI_PARALLEL
    level = level_max;
    MPI_Allreduce(&level,&level_max,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
    level = level_max_hydro;
    MPI_Allreduce(&level,&level_max_hydro,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
    level = level_max_nbody;
    MPI_Allreduce(&level,&level_max_nbody,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
    assert(level_max_hydro>=0);
#endif*/


    // Set fixed SPH timestep level here in case maximum has changed
    if (sph_single_timestep == 1) {
      for (i=0; i<mfv->Nhydro; i++) {
        MeshlessFVParticle<ndim>& part = mfv->GetMeshlessFVParticlePointer(i);
        if (part.itype == dead) continue;
        if (part.nlast == n) part.level = level_max_hydro;
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
        for (i=0; i<mfv->Nhydro; i++) {
          MeshlessFVParticle<ndim>& part = mfv->GetMeshlessFVParticlePointer(i);
          if (part.itype == dead) continue;
          part.nstep *= nfactor;
          part.nlast *= nfactor;
        }
        //for (i=0; i<nbody->Nnbody; i++) nbody->nbodydata[i]->nstep *= nfactor;
        //for (i=0; i<nbody->Nnbody; i++) nbody->nbodydata[i]->nlast *= nfactor;
      }
      else if (level_max < level_max_old) {
        nfactor = pow(2,level_max_old - level_max);
        n /= nfactor;
        for (i=0; i<mfv->Nhydro; i++) {
          MeshlessFVParticle<ndim>& part = mfv->GetMeshlessFVParticlePointer(i);
          if (part.itype == dead) continue;
          part.nlast /= nfactor;
          part.nstep /= nfactor;
        }
        //for (i=0; i<nbody->Nnbody; i++) nbody->nbodydata[i]->nlast /= nfactor;
        //for (i=0; i<nbody->Nnbody; i++) nbody->nbodydata[i]->nstep /= nfactor;
      }

    }
    //---------------------------------------------------------------------------------------------

    level_step = level_max + integration_step - 1;
    nresync    = pow(2,level_step);
    timestep   = dt_max / (DOUBLE) nresync;

    // Update values of nstep for both SPH and star particles
    for (i=0; i<mfv->Nhydro; i++) {
      MeshlessFVParticle<ndim>& part = mfv->GetMeshlessFVParticlePointer(i);
      if (part.itype == dead) continue;
      if (part.nlast == n) part.nstep = pow(2,level_step - part.level);
    }
    /*for (i=0; i<nbody->Nnbody; i++) {
      if (nbody->nbodydata[i]->nlast == n) nbody->nbodydata[i]->nstep =
        pow(2,level_step - nbody->nbodydata[i]->level);
    }*/

    assert(level_max >= level_max_old - 1);

  }
  //===============================================================================================


  // Various asserts for debugging
  assert(timestep >= 0.0);
  assert(level_step == level_max + integration_step - 1);
  assert(level_max_hydro <= level_max);
  assert(level_max_nbody <= level_max);
  assert(n <= nresync);
  for (i=0; i<mfv->Nhydro; i++) {
    MeshlessFVParticle<ndim>& part = mfv->GetMeshlessFVParticlePointer(i);
    if (part.itype == dead) continue;
    assert(part.level <= level_max);
    assert(part.nlast <= n);
    assert(part.tlast <= t);
    assert(part.nstep == pow(2,level_step - part.level));
    assert(part.nlast != n || n%part.nstep == 0);
  }
  /*for (i=0; i<nbody->Nnbody; i++) {
    assert(nbody->nbodydata[i]->level <= level_max);
    assert(nbody->nbodydata[i]->nlast <= n);
    assert(nbody->nbodydata[i]->nstep == pow(2,level_step - nbody->nbodydata[i]->level));
    assert(nbody->nbodydata[i]->nlast != n || n%nbody->nbodydata[i]->nstep == 0);
    assert(nbody->nbodydata[i]->level >= level_max_hydro);
    assert(nbody->nbodydata[i]->tlast <= t);
  }*/
  if (timestep <= 0.0) {
    cout << "Timestep fallen to zero : " << timestep << endl;
    exit(0);
  }

  timing->EndTimingSection("BLOCK_TIMESTEPS");

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
    if (part.active) Nactive++;
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
    exit(0);
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
  MeshlessFVParticle<ndim> *partdata = mfv->GetMeshlessFVParticleArray();

  for (int i=0; i<mfv->Nhydro; i++) {
    MeshlessFVParticle<ndim> &part = partdata[i];
    if (part.itype == dead) continue;
    for (int var=0; var<ndim+2; var++) part.Qcons[var] += part.dQ[var];
    mfv->ConvertQToConserved(part.volume, part.Qcons, part.Ucons);
    mfv->ConvertConservedToPrimitive(part.Ucons, part.Wprim);
    mfv->UpdateArrayVariables(part);
  }

  return;
}
