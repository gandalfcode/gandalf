//=============================================================================
//  GodunovSphSimulation.cpp
//  Contains all main functions controlling Godunov SPH simulation work-flow.
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
#include <iomanip>
#include <sstream>
#include <string>
#include <math.h>
#include <time.h>
#include <cstdio>
#include <cstring>
#include "Precision.h"
#include "Exception.h"
#include "Simulation.h"
#include "Parameters.h"
#include "InlineFuncs.h"
#include "Debug.h"
#include "Nbody.h"
#include "Sph.h"
#include "RiemannSolver.h"
#include "Sinks.h"
using namespace std;



// Create template class instances of the main GodunovSphSimulation object for
// each dimension used (1, 2 and 3)
template class GodunovSphSimulation<1>;
template class GodunovSphSimulation<2>;
template class GodunovSphSimulation<3>;



//=============================================================================
//  GodunovSphSimulation::ProcessSphParameters
/// Process parameter particular to setting up a Godunov SPH simulation object.
//=============================================================================
template <int ndim>
void GodunovSphSimulation<ndim>::ProcessSphParameters(void)
{
  aviscenum avisc;                  // Artificial viscosity enum
  acondenum acond;                  // Artificial conductivity enum
  eosenum gas_eos;                  // Gas EOS enum
  tdaviscenum tdavisc;              // Time-dependent viscosity enum

  map<string, int> &intparams = simparams->intparams;
  map<string, double> &floatparams = simparams->floatparams;
  map<string, string> &stringparams = simparams->stringparams;
  string KernelName = stringparams["kernel"];

  // Create SPH object based on chosen method in params file
  //---------------------------------------------------------------------------
  if (intparams["tabulated_kernel"] == 1) {
    sph = new GodunovSph<ndim, TabulatedKernel>
      (intparams["hydro_forces"], intparams["self_gravity"],
       floatparams["alpha_visc"], floatparams["beta_visc"],
       floatparams["h_fac"], floatparams["h_converge"],
       avisc, acond, tdavisc, stringparams["gas_eos"], KernelName);
  }
  else if (intparams["tabulated_kernel"] == 0) {
    if (KernelName == "gaussian") {
      sph = new GodunovSph<ndim, GaussianKernel>
        (intparams["hydro_forces"], intparams["self_gravity"],
         floatparams["alpha_visc"], floatparams["beta_visc"],
         floatparams["h_fac"], floatparams["h_converge"],
         avisc, acond, tdavisc, stringparams["gas_eos"], KernelName);
    }
    else {
      string message = "Unrecognised parameter : kernel = " +
	simparams->stringparams["kernel"];
      ExceptionHandler::getIstance().raise(message);
    }
  }
  else {
    string message = "Invalid option for the tabulated_kernel parameter: " +
      stringparams["tabulated_kernel"];
    ExceptionHandler::getIstance().raise(message);
  }


  // Riemann solver object
  //---------------------------------------------------------------------------
  string riemann = stringparams["riemann_solver"];
  if (riemann == "exact")
    sph->riemann = new ExactRiemannSolver(floatparams["gamma_eos"]);
  else if (riemann == "hllc")
    sph->riemann = new HllcRiemannSolver(floatparams["gamma_eos"]);
  else {
    string message = "Unrecognised parameter : riemann_solver = "
      + riemann;
    ExceptionHandler::getIstance().raise(message);
  }


  // Create SPH particle integration object
  //---------------------------------------------------------------------------
  sphint = new SphGodunovIntegration<ndim, GodunovSphParticle>(floatparams["accel_mult"],
					   floatparams["courant_mult"],
			               floatparams["energy_mult"],
					   gas_eos, tdavisc);


  // Energy integration object
  //---------------------------------------------------------------------------
  uint = new EnergyGodunovIntegration<ndim>(floatparams["energy_mult"]);

#if defined MPI_PARALLEL

  mpicontrol = new MpiControlType<ndim, GodunovSphParticle>;

#endif

  // Create neighbour searching object based on chosen method in params file
  //-------------------------------------------------------------------------
  if (stringparams["neib_search"] == "bruteforce")
    sphneib = new GodunovSphBruteForce<ndim,GodunovSphParticle>
     (sph->kernp->kernrange,&simbox,sph->kernp,this->timing);
  else if (stringparams["neib_search"] == "kdtree") {
    sphneib = new GodunovSphKDTree<ndim,GodunovSphParticle,KDTreeCell>
     (intparams["Nleafmax"],Nmpi,floatparams["thetamaxsqd"],sph->kernp->kernrange,
      floatparams["macerror"],stringparams["gravity_mac"],stringparams["multipole"],
      &simbox,sph->kernp,timing);
  }
  else if (stringparams["neib_search"] == "octtree") {
    sphneib = new GodunovSphOctTree<ndim,GodunovSphParticle,OctTreeCell>
     (intparams["Nleafmax"],Nmpi,floatparams["thetamaxsqd"],sph->kernp->kernrange,
      floatparams["macerror"],stringparams["gravity_mac"],stringparams["multipole"],
      &simbox,sph->kernp,timing);
  }
  else {
    string message = "Unrecognised parameter : neib_search = "
      + simparams->stringparams["neib_search"];
    ExceptionHandler::getIstance().raise(message);
  }

#if defined MPI_PARALLEL
  mpicontrol->SetNeibSearch(sphneib);
#endif


  // Set other important parameters
  sph->riemann_solver = stringparams["riemann_solver"];
  sph->slope_limiter  = stringparams["slope_limiter"];
  sph->riemann_order  = intparams["riemann_order"];

  // Create ghosts object
  if (IsAnyBoundarySpecial(simbox))
    LocalGhosts = new PeriodicGhostsSpecific<ndim,GodunovSphParticle >();
  else
    LocalGhosts = new NullGhosts<ndim>();
#ifdef MPI_PARALLEL
  MpiGhosts = new MPIGhostsSpecific<ndim, GodunovSphParticle>(mpicontrol);
#endif


  return;
}



//TODO: make this mess more modular (note: initial h computation
//should be done inside the neighbour search)
//=============================================================================
//  SphSimulation::PostInitialConditionsSetup
/// ..
//=============================================================================
template <int ndim>
void GodunovSphSimulation<ndim>::PostInitialConditionsSetup(void)
{
  int i;                            // Particle counter
  int k;                            // Dimension counter
  SphParticle<ndim> *partdata = sph->GetParticlesArray();

  debug2("[SphSimulation::PostInitialConditionsSetup]");

  // Set time variables here (for now)
  Noutsnap = 0;
  nresync = 0;
  //tsnapnext = dt_snap;

  // Set initial smoothing lengths and create initial ghost particles
  //---------------------------------------------------------------------------
  if (sph->Nsph > 0) {

    // Set all relevant particle counters
    sph->Nghost = 0;
    sph->Nghostmax = sph->Nsphmax - sph->Nsph;
    sph->Ntot = sph->Nsph;
    for (i=0; i<sph->Nsph; i++) {
      SphParticle<ndim>& part = sph->GetParticleIPointer(i);
      part.active = true;
    }

    // Compute mean mass
    sph->mmean = 0.0;
    for (i=0; i<sph->Nsph; i++) {
      SphParticle<ndim>& part = sph->GetParticleIPointer(i);
      sph->mmean += part.m;
    }
    sph->mmean /= (FLOAT) sph->Nsph;

    sph->InitialSmoothingLengthGuess();
    //sphneib->BuildTree(rebuild_tree,n,ntreebuildstep,ntreestockstep,timestep,sph);

    sphneib->neibcheck = false;
    sphneib->UpdateAllSphProperties(sph->Nsph,sph->Ntot,partdata,sph,nbody);

    // Search ghost particles
    sphneib->SearchBoundaryGhostParticles(0.0,simbox,sph);
#ifdef MPI_PARALLEL
    MpiGhosts->SearchGhostParticles(0.0,simbox,sph);
#endif

    level_step = 1;

    // Zero accelerations (perhaps here)
    for (i=0; i<sph->Ntot; i++) {
      SphParticle<ndim>& part = sph->GetParticleIPointer(i);
      part.active = true;
    }

    // Calculate all SPH properties
    //sphneib->BuildTree(rebuild_tree,n,ntreebuildstep,ntreestockstep,timestep,sph);
    sphneib->UpdateAllSphProperties(sph->Nsph,sph->Ntot,partdata,sph,nbody);

    // Search ghost particles
    sphneib->SearchBoundaryGhostParticles(0.0,simbox,sph);
#ifdef MPI_PARALLEL
    MpiGhosts->SearchGhostParticles(0.0,simbox,sph);
#endif

    // Update neighbour tre
    //sphneib->BuildTree(rebuild_tree,n,ntreebuildstep,ntreestockstep,timestep,sph);
    sphneib->neibcheck = true;
    sphneib->UpdateAllSphProperties(sph->Nsph,sph->Ntot,partdata,sph,nbody);
    sphneib->UpdateAllSphDerivatives(sph->Nsph,sph->Ntot,partdata,sph);

  }


  // Compute all initial N-body terms
  //---------------------------------------------------------------------------
  if (nbody->Nstar > 0) {

    // Zero all acceleration terms
    for (i=0; i<nbody->Nstar; i++) {
      for (k=0; k<ndim; k++) nbody->stardata[i].a[k] = 0.0;
      for (k=0; k<ndim; k++) nbody->stardata[i].adot[k] = 0.0;
      for (k=0; k<ndim; k++) nbody->stardata[i].a2dot[k] = 0.0;
      for (k=0; k<ndim; k++) nbody->stardata[i].a3dot[k] = 0.0;
      nbody->stardata[i].gpot = 0.0;
      nbody->stardata[i].active = true;
      nbody->stardata[i].level = level_step;
      nbody->stardata[i].nstep = 0;
      nbody->stardata[i].nlast = 0;
      nbody->nbodydata[i] = &(nbody->stardata[i]);
    }

    nbody->Nnbody = nbody->Nstar;
    nbody->CalculateDirectGravForces(nbody->Nnbody,nbody->nbodydata);
    //nbody->CalculateDirectSPHForces(nbody->Nnbody,sph->Nsph,
    //                              sph->sphdata,nbody->nbodydata);
    nbody->CalculateAllStartupQuantities(nbody->Nnbody,nbody->nbodydata);

  }


  // Compute all initial SPH force terms
  //---------------------------------------------------------------------------
  if (sph->Nsph > 0) {

    // Zero accelerations (here for now)
//    for (i=0; i<sph->Ntot; i++) {
//      for (k=0; k<ndim; k++) sph->sphdata[i].a[k] = (FLOAT) 0.0;
//      for (k=0; k<ndim; k++) sph->sphdata[i].agrav[k] = (FLOAT) 0.0;
//      sph->sphdata[i].gpot = (FLOAT) 0.0;
//      sph->sphdata[i].dudt = (FLOAT) 0.0;
//      sph->sphdata[i].active = true;
//      sph->sphdata[i].level = 0;
//      sph->sphdata[i].nstep = 0;
//      sph->sphdata[i].nlast = 0;
//    }

    LocalGhosts->CopySphDataToGhosts(simbox,sph);
#ifdef MPI_PARALLEL
    MpiGhosts->CopySphDataToGhosts(simbox,sph);
#endif
    //sphneib->BuildTree(rebuild_tree,n,ntreebuildstep,ntreestockstep,timestep,sph);

    // Compute timesteps for all particles
    if (Nlevels == 1)
      this->ComputeGlobalTimestep();
    else
      this->ComputeBlockTimesteps();

    // Calculate SPH gravity and hydro forces, depending on which are activated
    if (sph->hydro_forces == 1 && sph->self_gravity == 1)
      sphneib->UpdateAllSphForces(sph->Nsph,sph->Ntot,partdata,sph,nbody);
    else if (sph->hydro_forces == 1)
      sphneib->UpdateAllSphHydroForces(sph->Nsph,sph->Ntot,partdata,sph,nbody);
    else if (sph->self_gravity == 1)
      sphneib->UpdateAllSphGravForces(sph->Nsph,sph->Ntot,partdata,sph,nbody);

    // Compute contribution to grav. accel from stars
    for (i=0; i<sph->Nsph; i++) {
      SphParticle<ndim>& part = sph->GetParticleIPointer(i);
      if (part.active)
        sph->ComputeStarGravForces(nbody->Nnbody,nbody->nbodydata,part);
    }

    sphneib->UpdateAllSphDudt(sph->Nsph,sph->Ntot,partdata,sph);

    // Add accelerations
    for (i=0; i<sph->Nsph; i++) {
      SphParticle<ndim>& part = sph->GetParticleIPointer(i);
      part.active = false;
      for (k=0; k<ndim; k++)
        part.a[k] += part.agrav[k];
    }

    LocalGhosts->CopySphDataToGhosts(simbox,sph);
#ifdef MPI_PARALLEL
    MpiGhosts->CopySphDataToGhosts(simbox,sph);
#endif

  }

  // Set particle values for initial step (e.g. r0, v0, a0)
  if (simparams->stringparams["gas_eos"] == "energy_eqn")
    uint->EndTimestep(n,sph->Nsph,timestep,sph->GetParticlesArray());
  sphint->EndTimestep(n,timestep,sph->Nsph,sph->GetParticlesArray());
  nbody->EndTimestep(n,nbody->Nstar,t,timestep,nbody->nbodydata);

  // Compute timesteps for all particles
  nresync = n;
  if (Nlevels == 1)
    this->ComputeGlobalTimestep();
  else
    this->ComputeBlockTimesteps();

  this->CalculateDiagnostics();
  this->OutputDiagnostics();
  this->diag0 = this->diag;
  this->setup = true;

  return;
}



//=============================================================================
//  GodunovSphSimulation::MainLoop
/// Main Godunov SPH simulation integration loop.
//=============================================================================
template <int ndim>
void GodunovSphSimulation<ndim>::MainLoop(void)
{
  //int activecount;                  // ..
  int i;                            // Particle loop counter
  int it;                           // Time-symmetric iteration counter
  int k;                            // Dimension counter
  SphParticle<ndim> *partdata = sph->GetParticlesArray();

  debug2("[GodunovSphSimulation::MainLoop]");

  //cin >> i;

  // Advance time variables
  n = n + 1;
  Nsteps = Nsteps + 1;
  t = t + timestep;

  // Advance SPH particles positions and velocities
  sphint->AdvanceParticles(n,(FLOAT) timestep,sph->Nsph,sph->GetParticlesArray());
  if (simparams->stringparams["gas_eos"] == "energy_eqn")
    uint->EnergyIntegration(n,sph->Nsph,sph->GetParticlesArray(),(FLOAT) timestep);
  nbody->AdvanceParticles(n,nbody->Nstar,t,timestep,nbody->nbodydata);

  // Check all boundary conditions
  sphint->CheckBoundaries(simbox,sph);

  //---------------------------------------------------------------------------
  if (sph->Nsph > 0) {

    // Reorder particles

    // Search ghost particles
    sphneib->SearchBoundaryGhostParticles(0.0,simbox,sph);
#ifdef MPI_PARALLEL
    MpiGhosts->SearchGhostParticles(0.0,simbox,sph);
#endif

    // Update neighbour tree
    //sphneib->BuildTree(rebuild_tree,n,ntreebuildstep,ntreestockstep,timestep,sph);

    // Update cells containing active particles
    sphneib->UpdateActiveParticleCounters(partdata,sph);

    // Calculate all SPH properties
    sphneib->UpdateAllSphProperties(sph->Nsph,sph->Ntot,partdata,sph,nbody);

    // Compute timesteps for all particles
    if (Nlevels == 1)
      this->ComputeGlobalTimestep();
    else
      this->ComputeBlockTimesteps();

    sphneib->UpdateAllSphDerivatives(sph->Nsph,sph->Ntot,partdata,sph);

    // Copy properties from original particles to ghost particles
    LocalGhosts->CopySphDataToGhosts(simbox,sph);
#ifdef MPI_PARALLEL
    MpiGhosts->CopySphDataToGhosts(simbox,sph);
#endif

//    // Zero accelerations (perhaps)
//    for (i=0; i<sph->Ntot; i++) {
//      if (sph->sphdata[i].active) {
//        for (k=0; k<ndim; k++) sph->sphdata[i].a[k] = (FLOAT) 0.0;
//        for (k=0; k<ndim; k++) sph->sphdata[i].agrav[k] = (FLOAT) 0.0;
//        sph->sphdata[i].gpot = (FLOAT) 0.0;
//        sph->sphdata[i].gpe = (FLOAT) 0.0;
//        sph->sphdata[i].dudt = (FLOAT) 0.0;
//        sph->sphdata[i].levelneib = 0;
//      }
//    }

    // Calculate SPH gravity and hydro forces, depending on which are activated
    if (sph->hydro_forces == 1 && sph->self_gravity == 1)
      sphneib->UpdateAllSphForces(sph->Nsph,sph->Ntot,partdata,sph,nbody);
    else if (sph->hydro_forces == 1)
      sphneib->UpdateAllSphHydroForces(sph->Nsph,sph->Ntot,partdata,sph,nbody);
    else if (sph->self_gravity == 1)
      sphneib->UpdateAllSphGravForces(sph->Nsph,sph->Ntot,partdata,sph,nbody);

    // Compute contribution to grav. accel from stars
    for (i=0; i<sph->Nsph; i++) {
      SphParticle<ndim>& part = sph->GetParticleIPointer(i);
      if (part.active) {
        sph->ComputeStarGravForces(nbody->Nnbody,nbody->nbodydata,part);
      }
    }

    // Add accelerations
    if (sph->self_gravity == 1) {
      for (i=0; i<sph->Nsph; i++) {
        SphParticle<ndim>& part = sph->GetParticleIPointer(i);
        if (part.active) {
          for (k=0; k<ndim; k++)
            part.a[k] += part.agrav[k];
        }
      }
    }

    // Now update compressional heating rate
    sphneib->UpdateAllSphDudt(sph->Nsph,sph->Ntot,partdata,sph);

    // Set all end-of-step variables
    if (simparams->stringparams["gas_eos"] == "energy_eqn")
      uint->EndTimestep(n,sph->Nsph,timestep,sph->GetParticlesArray());
    sphint->EndTimestep(n,timestep,sph->Nsph,sph->GetParticlesArray());

  }
  //---------------------------------------------------------------------------


  // Compute N-body forces
  //---------------------------------------------------------------------------
  if (nbody->Nnbody > 0) {

    // Iterate end-of-step
    //-------------------------------------------------------------------------
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

      // Calculate forces, force derivatives etc.., for active stars/systems
      nbody->CalculateDirectGravForces(nbody->Nnbody,nbody->nbodydata);
      //nbody->CalculateDirectSPHForces(nbody->Nnbody,sph->Nsph,
      //                              sph->sphdata,nbody->nbodydata);

      // Calculate correction step for all stars at end of step
      nbody->CorrectionTerms(n,nbody->Nnbody,t,timestep,nbody->nbodydata);

    }
    //-------------------------------------------------------------------------

    nbody->EndTimestep(n,nbody->Nnbody,t,timestep,nbody->nbodydata);

  }
  //---------------------------------------------------------------------------


  return;
}



//=============================================================================
//  GodunovSphSimulation::ComputeGlobalTimestep
/// Computes global timestep for Godunov SPH simulation.  Calculates the
/// minimum timestep for all SPH and N-body particles in the simulation.
//=============================================================================
template <int ndim>
void GodunovSphSimulation<ndim>::ComputeGlobalTimestep(void)
{
  int i;                            // Particle counter
  DOUBLE dt;                        // Particle timestep
  DOUBLE dt_min = big_number_dp;    // Local copy of minimum timestep

  debug2("[SphSimulation::ComputeGlobalTimestep]");

  //---------------------------------------------------------------------------
  if (n == nresync) {

    n = 0;
    level_max = 0;
    level_step = level_max + integration_step - 1;
    nresync = integration_step;

    // Find minimum timestep from all SPH particles
    //-------------------------------------------------------------------------
#pragma omp parallel default(none) private(i,dt) shared(dt_min)
    {
      dt = big_number_dp;
#pragma omp for
      for (i=0; i<sph->Nsph; i++) {
        SphParticle<ndim>& part = sph->GetParticleIPointer(i);
        part.dt = sphint->Timestep(part,sph);
        dt = min(dt,part.dt);
      }

      // If integrating energy equation, include energy timestep
      if (simparams->stringparams["gas_eos"] == "energy_eqn") {
#pragma omp for
        for (i=0; i<sph->Nsph; i++) {
          SphParticle<ndim>& part = sph->GetParticleIPointer(i);
          part.dt = min(part.dt,
                                   uint->Timestep(part));
          dt = min(dt,part.dt);
        }

      }

#pragma omp critical
      {
        if (dt < dt_min) dt_min = dt;
      }
    }
    //-------------------------------------------------------------------------


    // Now compute minimum timestep due to stars/systems
    for (i=0; i<nbody->Nnbody; i++)
      dt_min = min(dt_min,nbody->Timestep(nbody->nbodydata[i]));


    // Set all particles to same timestep
    timestep = dt_min;
    for (i=0; i<sph->Nsph; i++) {
      SphParticle<ndim>& part = sph->GetParticleIPointer(i);
      part.level = 0;
      part.levelneib = 0;
      part.dt = timestep;
      part.nstep = pow(2,level_step - part.level);
      part.nlast = n;

    }
    for (i=0; i<nbody->Nnbody; i++) {
      nbody->nbodydata[i]->level = 0;
      nbody->nbodydata[i]->nstep =
        pow(2,level_step - nbody->nbodydata[i]->level);
      nbody->nbodydata[i]->nlast = n;
      nbody->nbodydata[i]->dt = timestep;
    }

  }
  //---------------------------------------------------------------------------

  assert(timestep > 0.0);

  return;
}



//=============================================================================
//  GodunovSphSimulation::ComputeBlockTimesteps
/// ..
//=============================================================================
template <int ndim>
void GodunovSphSimulation<ndim>::ComputeBlockTimesteps(void)
{
  int i;                                // Particle counter
  int istep;                            // ??
  int level;                            // Particle timestep level
  int last_level;                       // Previous timestep level
  int level_max_old;                    // Old level_max
  int level_max_sph = 0;                // level_max for SPH particles only
  int level_max_nbody = 0;              // level_max for SPH particles only
  int nstep;                            // ??
  int nfactor;                          // ??
  DOUBLE dt;                            // Aux. timestep variable
  DOUBLE dt_min_sph = big_number_dp;    // Maximum SPH particle timestep
  DOUBLE dt_min_nbody = big_number_dp;  // Maximum N-body particle timestep

  debug2("[SphSimulation::ComputeBlockTimesteps]");

  // Synchronise all timesteps and reconstruct block timestep structure.
  //===========================================================================
  if (n == nresync) {

    n = 0;
    timestep = big_number_dp;
    for (i=0; i<sph->Nsph; i++) {
      SphParticle<ndim>& part = sph->GetParticleIPointer(i);
      part.dt = big_number_dp;
    }
    for (i=0; i<nbody->Nnbody; i++) nbody->nbodydata[i]->dt = big_number_dp;

    // If integrating energy equation, calculate the explicit energy timestep
    if (sph->gas_eos == "energy_eqn") {
      for (i=0; i<sph->Nsph; i++) {
        SphParticle<ndim>& part = sph->GetParticleIPointer(i);
        part.dt = uint->Timestep(part);
      }
    }

    // Find minimum timestep from all SPH particles
    for (i=0; i<sph->Nsph; i++) {
      SphParticle<ndim>& part = sph->GetParticleIPointer(i);
      dt = min(part.dt,
	       sphint->Timestep(part,sph));
      if (dt < timestep) timestep = dt;
      if (dt < dt_min_sph) dt_min_sph = dt;
      part.dt = dt;
    }

    // Now compute minimum timestep due to stars/systems
    for (i=0; i<nbody->Nnbody; i++) {
      dt = nbody->Timestep(nbody->nbodydata[i]);
      if (dt < timestep) timestep = dt;
      if (dt < dt_min_nbody) dt_min_nbody = dt;
      nbody->nbodydata[i]->dt = dt;
    }

    // Calculate new block timestep levels
    level_max = Nlevels - 1;
    level_step = level_max + integration_step - 1;
    dt_max = timestep*powf(2.0,level_max);

    // Calculate the maximum level occupied by all SPH particles
    level_max_sph = min((int) (invlogetwo*log(dt_max/dt_min_sph)) + 1,
			level_max);
    level_max_nbody = min((int) (invlogetwo*log(dt_max/dt_min_nbody)) + 1,
			  level_max);

    // If enforcing a single SPH timestep, set it here.  Otherwise, populate
    // the timestep levels with SPH particles.
    if (sph_single_timestep == 1)
      for (i=0; i<sph->Nsph; i++) {
        SphParticle<ndim>& part = sph->GetParticleIPointer(i);
        part.level = level_max_sph;
        part.levelneib = level_max_sph;
      }
    else {
      for (i=0; i<sph->Nsph; i++) {
        SphParticle<ndim>& part = sph->GetParticleIPointer(i);
        dt = part.dt;
        level = min((int) (invlogetwo*log(dt_max/dt)) + 1, level_max);
        level = max(level,0);
        part.level = level;
        part.levelneib = level;
        part.nlast = n;
        part.nstep = pow(2,level_step - part.level);
      }
    }

    // Populate timestep levels with N-body particles
    for (i=0; i<nbody->Nnbody; i++) {
      dt = nbody->Timestep(nbody->nbodydata[i]);
      level = min((int) (invlogetwo*log(dt_max/dt)) + 1, level_max);
      level = max(level,0);
      nbody->nbodydata[i]->level = level;
      nbody->nbodydata[i]->nlast = n;
      nbody->nbodydata[i]->nstep =
        pow(2,level_step - nbody->nbodydata[i]->level);
    }

    nresync = pow(2,level_step);
    timestep = dt_max / (DOUBLE) nresync;

  }

  // If not resynchronising, check if any SPH/N-body particles need to move
  // up or down timestep levels.
  //===========================================================================
  else {

    level_max_old = level_max;
    level_max = 0;

    // Find all SPH particles at the beginning of a new timestep
    //-------------------------------------------------------------------------
    for (i=0; i<sph->Nsph; i++) {

      SphParticle<ndim>& part = sph->GetParticleIPointer(i);

      // Skip particles that are not at end of step
      if (part.nlast == n) {
        nstep = part.nstep;
        last_level = part.level;

        // Compute new timestep value and level number
        dt = sphint->Timestep(part,sph);
        if (sph->gas_eos == "energy_eqn")
          dt = min(dt,uint->Timestep(part));
        part.dt = dt;
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

        part.nlast = n;
        part.nstep = pow(2,level_step - part.level);
      }

      // Find maximum level of all SPH particles
      level_max_sph = max(level_max_sph,part.level);
      level_max = max(level_max,part.level);
    }
    //-------------------------------------------------------------------------


    // Now find all N-body particles at the beginning of a new timestep
    //-------------------------------------------------------------------------
    for (i=0; i<nbody->Nnbody; i++) {

      // Skip particles that are not at end of step
      if (nbody->nbodydata[i]->nlast == n) {
	nstep = nbody->nbodydata[i]->nstep;
	last_level = nbody->nbodydata[i]->level;

	// Compute new timestep value and level number
	dt = nbody->Timestep(nbody->nbodydata[i]);
	nbody->nbodydata[i]->dt = dt;
	level = min((int) (invlogetwo*log(dt_max/dt)) + 1, 0);

	// Move up one level (if levels are correctly synchronised) or
	// down several levels if required
	if (level < last_level && last_level > 1 && n%(2*nstep) == 0)
	  nbody->nbodydata[i]->level = last_level - 1;
	else if (level > last_level)
	  nbody->nbodydata[i]->level = level;
	else
	  nbody->nbodydata[i]->level = last_level;

	nbody->nbodydata[i]->nlast = n;
	nbody->nbodydata[i]->nstep =
	  pow(2,level_step - nbody->nbodydata[i]->level);
      }

      // Find maximum level of all N-body particles
      level_max_nbody = max(level_max_nbody,nbody->nbodydata[i]->level);
      level_max = max(level_max,nbody->nbodydata[i]->level);
    }
    //-------------------------------------------------------------------------


    // Set fixed SPH timestep level here in case maximum has changed
    if (sph_single_timestep == 1) {
      for (i=0; i<sph->Nsph; i++) {
        SphParticle<ndim>& part = sph->GetParticleIPointer(i);
        part.level = level_max_sph;
        part.nstep = pow(2,level_step - part.level);
      }
    }

    // For now, don't allow levels to be removed
    level_max = max(level_max,level_max_old);
    level_step = level_max + integration_step - 1;

    for (i=0; i<sph->Nsph; i++) {
      SphParticle<ndim>& part = sph->GetParticleIPointer(i);
      if (part.nlast == n)
        part.nstep = pow(2,level_step - part.level);
    }
    for (i=0; i<nbody->Nnbody; i++) {
      if (nbody->nbodydata[i]->nlast == n) nbody->nbodydata[i]->nstep =
	pow(2,level_step - nbody->nbodydata[i]->level);
    }

    // Update all timestep variables if we have removed or added any levels
    //-------------------------------------------------------------------------
    if (level_max != level_max_old) {

      // Increase maximum timestep level if correctly synchronised
      istep = pow(2,level_step - level_max_old + 1);
      if (level_max <= level_max_old - 1 && level_max_old > 1 && n%istep == 0)
	level_max = level_max_old - 1;
      else if (level_max == level_max_old)
	level_max = level_max_old;

      // Adjust integer time if levels added or removed
      if (level_max > level_max_old) {
	nfactor = pow(2,level_max - level_max_old);
	n *= nfactor;
	for (i=0; i<sph->Nsph; i++) {
	  SphParticle<ndim>& part = sph->GetParticleIPointer(i);
	  part.nlast *= nfactor;
	  part.nstep *= nfactor;
	}
	for (i=0; i<nbody->Nnbody; i++) nbody->nbodydata[i]->nlast *= nfactor;
	for (i=0; i<nbody->Nnbody; i++) nbody->nbodydata[i]->nstep *= nfactor;
      }
      else if (level_max < level_max_old) {
	nfactor = pow(2,level_max_old - level_max);
	n /= nfactor;
	for (i=0; i<sph->Nsph; i++) {
	  SphParticle<ndim>& part = sph->GetParticleIPointer(i);
	  part.nlast /= nfactor;
	  part.nstep /= nfactor;
	}
	for (i=0; i<nbody->Nnbody; i++) nbody->nbodydata[i]->nlast /= nfactor;
	for (i=0; i<nbody->Nnbody; i++) nbody->nbodydata[i]->nstep /= nfactor;
      }

    }
    //-------------------------------------------------------------------------

    nresync = pow(2,level_step);
    timestep = dt_max / (DOUBLE) nresync;

  }
  //===========================================================================

  return;
}
