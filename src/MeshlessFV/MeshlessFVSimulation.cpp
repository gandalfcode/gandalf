//=================================================================================================
//  MeshlessFVSimulation.cpp
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
  //ProcessSphParameters();
  // Create 'grad-h' SPH object depending on choice of kernel
  //===============================================================================================
  if (intparams["tabulated_kernel"] == 1) {
    mfv = new LV2008MFV<ndim, TabulatedKernel>
      (intparams["hydro_forces"], intparams["self_gravity"], floatparams["h_fac"],
       floatparams["h_converge"], floatparams["gamma_eos"], stringparams["gas_eos"],
       KernelName, sizeof(MeshlessFVParticle<ndim>));
  }
  else if (intparams["tabulated_kernel"] == 0) {
    if (KernelName == "m4") {
      mfv = new LV2008MFV<ndim, M4Kernel>
        (intparams["hydro_forces"], intparams["self_gravity"], floatparams["h_fac"],
         floatparams["h_converge"], floatparams["gamma_eos"], stringparams["gas_eos"],
         KernelName, sizeof(MeshlessFVParticle<ndim>));
    }
    else if (KernelName == "quintic") {
      mfv = new LV2008MFV<ndim, QuinticKernel>
        (intparams["hydro_forces"], intparams["self_gravity"], floatparams["h_fac"],
         floatparams["h_converge"], floatparams["gamma_eos"], stringparams["gas_eos"],
         KernelName, sizeof(MeshlessFVParticle<ndim>));
    }
    else if (KernelName == "gaussian") {
      mfv = new LV2008MFV<ndim, GaussianKernel>
        (intparams["hydro_forces"], intparams["self_gravity"], floatparams["h_fac"],
         floatparams["h_converge"], floatparams["gamma_eos"], stringparams["gas_eos"],
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
  hydro = mfv;


  // Riemann solver object
  //-----------------------------------------------------------------------------------------------
  string riemann = stringparams["riemann_solver"];
  if (riemann == "exact")
    mfv->riemann = new ExactRiemannSolver<ndim>(floatparams["gamma_eos"]);
  else if (riemann == "hllc")
    mfv->riemann = new HllcRiemannSolver<ndim>(floatparams["gamma_eos"]);
  else {
    string message = "Unrecognised parameter : riemann_solver = " + riemann;
    ExceptionHandler::getIstance().raise(message);
  }



  // Create neighbour searching object based on chosen method in params file
  //-----------------------------------------------------------------------------------------------
  if (stringparams["neib_search"] == "bruteforce")
    mfvneib = new MeshlessFVBruteForce<ndim,MeshlessFVParticle>
      (mfv->kernp->kernrange, &simbox, mfv->kernp, timing);
  /*else if (stringparams["neib_search"] == "kdtree") {
    mfvneib = new GradhSphKDTree<ndim,GradhMeshlessFVParticle,KDTreeCell>
     (intparams["Nleafmax"],Nmpi,floatparams["thetamaxsqd"],mfv->kernp->kernrange,
      floatparams["macerror"],stringparams["gravity_mac"],stringparams["multipole"],
      &simbox,mfv->kernp,timing);
  }
  else if (stringparams["neib_search"] == "octtree") {
    mfvneib = new GradhSphOctTree<ndim,GradhMeshlessFVParticle,OctTreeCell>
     (intparams["Nleafmax"],Nmpi,floatparams["thetamaxsqd"],mfv->kernp->kernrange,
      floatparams["macerror"],stringparams["gravity_mac"],stringparams["multipole"],
      &simbox,mfv->kernp,timing);
  }*/
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
  nbody->perturbers     = intparams["perturbers"];
  if (intparams["sub_systems"] == 1) subsystem->perturbers = intparams["perturbers"];


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
  /*if (sim == "sph" || sim == "gradhsph" || sim == "sm2012sph" || sim == "godunov_sph") {
    sinks.timing    = timing;
    sphint->timing  = timing;
    mfvneib->timing = timing;
    uint->timing    = timing;
    radiation->timing = timing;
  }*/


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
  MeshlessFVParticle<ndim> *partdata;         // Pointer to main SPH data array

  debug2("[MeshlessFVSimulation::PostInitialConditionsSetup]");

  // Set iorig
  if (rank == 0) {
    for (i=0; i<mfv->Nhydro; i++) mfv->GetMeshlessFVParticlePointer(i).iorig = i;
  }

  // Set pointer to SPH particle data
  partdata = mfv->GetMeshlessFVParticleArray();

  // Set time variables here (for now)
  nresync = 0;
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
      mfvneib->BuildTree(rebuild_tree,0,ntreebuildstep,ntreestockstep,mfv->Ntot,
                         mfv->Nhydromax,mfv->GetMeshlessFVParticleArray(),mfv,timestep);
      mfvneib->UpdateAllProperties(mfv->Nhydro,mfv->Ntot,mfv->GetMeshlessFVParticleArray(),mfv,nbody);
    }
    else {
      mfvneib->BuildTree(rebuild_tree,0,ntreebuildstep,ntreestockstep,mfv->Ntot,
                         mfv->Nhydromax,mfv->GetMeshlessFVParticleArray(),mfv,timestep);
    }

    // Search ghost particles
    mfvneib->SearchBoundaryGhostParticles(0.0,simbox,mfv);
    //mfvneib->BuildGhostTree(true,0,ntreebuildstep,ntreestockstep,mfv->Ntot,
    //                        mfv->Nhydromax,mfv->GetMeshlessFVParticleArray(),sph,timestep);


    // Zero accelerations
    for (i=0; i<mfv->Nhydro; i++) mfv->GetMeshlessFVParticlePointer(i).active = true;

    // Update neighbour tree
    rebuild_tree = true;
    mfvneib->BuildTree(rebuild_tree, 0, ntreebuildstep, ntreestockstep, mfv->Ntot,
                       mfv->Nhydromax, mfv->GetMeshlessFVParticleArray(), mfv, timestep);
    level_step = 1;


    // For Eigenvalue MAC, need non-zero values
    for (i=0; i<mfv->Nhydro; i++) mfv->GetMeshlessFVParticlePointer(i).gpot = big_number;

    // Calculate all SPH properties
    mfvneib->UpdateAllProperties(mfv->Nhydro,mfv->Ntot,mfv->GetMeshlessFVParticleArray(),mfv,nbody);

    // Search ghost particles
    mfvneib->SearchBoundaryGhostParticles(0.0,simbox,mfv);
    //mfvneib->BuildGhostTree(true,0,ntreebuildstep,ntreestockstep,mfv->Ntot,
    //                        mfv->Nhydromax,mfv->GetMeshlessFVParticleArray(),sph,timestep);

    // Update neighbour tree
    rebuild_tree = true;
    mfvneib->BuildTree(rebuild_tree,0,ntreebuildstep,ntreestockstep,
                       mfv->Ntot,mfv->Nhydromax,mfv->GetMeshlessFVParticleArray(),mfv,timestep);
    //mfvneib->SearchBoundaryGhostParticles(0.0,simbox,sph);
    //mfvneib->BuildGhostTree(true,0,ntreebuildstep,ntreestockstep,mfv->Ntot,
    //                        mfv->Nhydromax,mfv->GetMeshlessFVParticleArray(),sph,timestep);
    mfvneib->neibcheck = true;


  }


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
    //LocalGhosts->CopySphDataToGhosts(simbox,sph);

    mfvneib->BuildTree(true,0,ntreebuildstep,ntreestockstep,mfv->Ntot,
                       mfv->Nhydromax,mfv->GetMeshlessFVParticleArray(),mfv,timestep);
    mfvneib->SearchBoundaryGhostParticles(0.0,simbox,mfv);
    //mfvneib->BuildGhostTree(true,0,ntreebuildstep,ntreestockstep,mfv->Ntot,
    //                        mfv->Nhydromax,mfv->GetMeshlessFVParticleArray(),sph,timestep);


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
//    MpiGhosts->CopySphDataToGhosts(simbox,sph);
#endif

    // Update the primitive vectors for all particles
    for (i=0; i<mfv->Ntot; i++) {
      mfv->ComputeThermalProperties(partdata[i]);
      mfv->UpdatePrimitiveVector(partdata[i]);
      mfv->ConvertPrimitiveToConserved(partdata[i].Wprim, partdata[i].Ucons);
      mfv->ConvertConservedToQ(partdata[i].volume, partdata[i].Ucons, partdata[i].Qcons);
      //for (k=0; k<ndim; k++) partdata[i].v0[k] = partdata[i].v[k];
    }

  }


  // Set particle values for initial step (e.g. r0, v0, a0, u0)
  /*uint->EndTimestep(n,mfv->Nhydro,t,timestep,mfv->GetMeshlessFVParticleArray());
  sphint->EndTimestep(n,mfv->Nhydro,t,timestep,mfv->GetMeshlessFVParticleArray());
  nbody->EndTimestep(n,nbody->Nstar,t,timestep,nbody->nbodydata);*/

  this->CalculateDiagnostics();
  this->diag0 = this->diag;
  this->setup = true;


  return;
}



//=================================================================================================
//  MeshlessFVSimulation::MainLoop
/// Main SPH simulation integration loop.
//=================================================================================================
template <int ndim>
void MeshlessFVSimulation<ndim>::MainLoop(void)
{
  int activecount = 0;                 // Flag if we need to recompute particles
  int i;                               // Particle loop counter
  int it;                              // Time-symmetric iteration counter
  int k;                               // Dimension counter
  FLOAT tghost;                        // Approx. ghost particle lifetime
  MeshlessFVParticle<ndim> *partdata;         // Pointer to main SPH data array

  debug2("[MeshlessFVSimulation::MainLoop]");

//cin >> i;
  // Set pointer for SPH data array
  partdata = mfv->GetMeshlessFVParticleArray();

  // Compute timesteps for all particles
  this->ComputeGlobalTimestep();
  //if (Nlevels == 1) this->ComputeGlobalTimestep();
  //else this->ComputeBlockTimesteps();

  // Advance time variables
  n = n + 1;
  Nsteps = Nsteps + 1;
  t = t + timestep;
  if (n == nresync) Nblocksteps = Nblocksteps + 1;
  if (n%integration_step == 0) Nfullsteps = Nfullsteps + 1;

  // Advance SPH and N-body particles' positions and velocities
  /*uint->EnergyIntegration(n,mfv->Nhydro,(FLOAT) t,(FLOAT) timestep,mfv->GetMeshlessFVParticleArray());
  sphint->AdvanceParticles(n,mfv->Nhydro,(FLOAT) t,(FLOAT) timestep,mfv->GetMeshlessFVParticleArray());
  nbody->AdvanceParticles(n,nbody->Nnbody,t,timestep,nbody->nbodydata);*/

  // Check all boundary conditions
  // (DAVID : Move this function to sphint and create an analagous one
  //  for N-body.  Also, only check this on tree-build steps)
  //if (Nsteps%ntreebuildstep == 0 || rebuild_tree) sphint->CheckBoundaries(simbox,sph);


  // Compute all SPH quantities
  //-----------------------------------------------------------------------------------------------
  if (mfv->Nhydro > 0) {

    // Rebuild or update local neighbour and gravity tree
    mfvneib->BuildTree(rebuild_tree,Nsteps,ntreebuildstep,ntreestockstep,
                       mfv->Ntot,mfv->Nhydromax,partdata,mfv,timestep);


    // Search for new ghost particles and create on local processor
    //if (Nsteps%ntreebuildstep == 0 || rebuild_tree) {
      tghost = timestep*(FLOAT)(ntreebuildstep - 1);
      mfvneib->SearchBoundaryGhostParticles(tghost,simbox,mfv);
      //mfvneib->BuildGhostTree(rebuild_tree,Nsteps,ntreebuildstep,ntreestockstep,
        //                      mfv->Ntot,mfv->Nhydromax,partdata,mfv,timestep);
    //}
    // Otherwise copy properties from original particles to ghost particles
    /*else {
      LocalGhosts->CopySphDataToGhosts(simbox, sph);
#ifdef MPI_PARALLEL
      MpiGhosts->CopySphDataToGhosts(simbox, sph);
#endif
    }*/

    for (i=0; i<mfv->Nhydro; i++) {
      MeshlessFVParticle<ndim>& part = mfv->GetMeshlessFVParticlePointer(i);
      part.active = true;
    }

    // Calculate all properties
    mfvneib->UpdateAllProperties(mfv->Nhydro, mfv->Ntot, partdata, mfv, nbody);

    mfv->CopyDataToGhosts(simbox, partdata);

    mfvneib->UpdateGradientMatrices(mfv->Nhydro, mfv->Ntot, partdata, mfv, nbody);

    mfv->CopyDataToGhosts(simbox, partdata);

    mfvneib->UpdateGodunovFluxes(mfv->Nhydro, mfv->Ntot, timestep, partdata, mfv, nbody);

    // Integrate all conserved variables to end of timestep
    for (i=0; i<mfv->Nhydro; i++) {
      mfv->IntegrateConservedVariables(partdata[i], timestep);
      mfv->ConvertQToConserved(partdata[i].volume, partdata[i].Qcons, partdata[i].Ucons);
      mfv->ConvertConservedToPrimitive(partdata[i].Ucons, partdata[i].Wprim);
      mfv->UpdateArrayVariables(partdata[i]);

      //for (int k=0; k<ndim; k++) partdata[i].r[k] += partdata[i].v[k]*timestep;
      for (k=0; k<ndim; k++) {
        partdata[i].r[k] += (FLOAT) 0.5*(partdata[i].v0[k] + partdata[i].v[k])*timestep;
        partdata[i].v0[k] = partdata[i].v[k];
        //partdata[i].r[k] += partdata[i].v[k]*timestep;
      }
      if (partdata[i].r[0] < simbox.boxmin[0])
        if (simbox.boundary_lhs[0] == periodicBoundary) {
          partdata[i].r[0] += simbox.boxsize[0];
        }
      if (partdata[i].r[0] > simbox.boxmax[0])
        if (simbox.boundary_rhs[0] == periodicBoundary) {
          partdata[i].r[0] -= simbox.boxsize[0];
        }

    }

    this->CalculateDiagnostics();
    this->OutputDiagnostics();
    this->UpdateDiagnostics();

  }
  //-----------------------------------------------------------------------------------------------



  rebuild_tree = false;


  // End-step terms for all SPH particles
  /*if (mfv->Nhydro > 0) {
    uint->EndTimestep(n,mfv->Nhydro,(FLOAT) t,(FLOAT) timestep,mfv->GetMeshlessFVParticleArray());
    sphint->EndTimestep(n,mfv->Nhydro,(FLOAT) t,(FLOAT) timestep,mfv->GetMeshlessFVParticleArray());
  }*/

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
  DOUBLE dt;                           // Particle timestep
  DOUBLE dt_min = big_number_dp;       // Local copy of minimum timestep
  DOUBLE dt_nbody;                     // Aux. minimum N-body timestep
  DOUBLE dt_sph;                       // Aux. minimum SPH timestep

  debug2("[MeshlessFVSimulation::ComputeGlobalTimestep]");
  timing->StartTimingSection("GLOBAL_TIMESTEPS",2);


  // Only update timestep when all particles are synced at end of last step.
  //-----------------------------------------------------------------------------------------------
  if (n == nresync) {

    n            = 0;
    level_max    = 0;
    level_step   = level_max + integration_step - 1;
    nresync      = integration_step;
    dt_min_sph   = big_number_dp;

    // Find minimum timestep from all SPH particles
    //---------------------------------------------------------------------------------------------
#pragma omp parallel default(none) private(i,dt,dt_nbody,dt_sph) \
  shared(dt_min) //,dt_min_nbody,dt_min_sph)
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
