//=================================================================================================
//  NbodySimulation.cpp
//  Contains all main functions controlling N-body simulation work-flow.
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
#include <sstream>
#include <string>
#include <math.h>
#include <time.h>
#include <cstdio>
#include <cstring>
#include "Precision.h"
#include "Exception.h"
#include "Debug.h"
#include "Hydrodynamics.h"
#include "InlineFuncs.h"
#include "Simulation.h"
#include "Parameters.h"
#include "Nbody.h"
#include "Ghosts.h"
using namespace std;



// Create template class instances of the NbodySimulation object for each dimension (1, 2 and 3)
template class NbodySimulation<1>;
template class NbodySimulation<2>;
template class NbodySimulation<3>;



//=================================================================================================
//  NbodySimulation::ProcessParameters
/// Process all the options chosen in the parameters file, setting various
/// simulation variables and creating important simulation objects.
/// Nbody specific version
//=================================================================================================
template <int ndim>
void NbodySimulation<ndim>::ProcessParameters(void)
{
  // Local references to parameter variables for brevity
  map<string, int> &intparams = simparams->intparams;
  map<string, double> &floatparams = simparams->floatparams;
  map<string, string> &stringparams = simparams->stringparams;
  string sim = stringparams["sim"];
  string KernelName = stringparams["kernel"];

  debug2("[NbodySimulation::ProcessParameters]");


  // Common set-up
  Simulation<ndim>::ProcessParameters();



  // Boundary condition variables
  //-----------------------------------------------------------------------------------------------
  simbox.boundary_lhs[0] = setBoundaryType(stringparams["boundary_lhs[0]"]);
  simbox.boundary_rhs[0] = setBoundaryType(stringparams["boundary_rhs[0]"]);
  simbox.min[0] = floatparams["boxmin[0]"]/simunits.r.outscale;
  simbox.max[0] = floatparams["boxmax[0]"]/simunits.r.outscale;

  if (ndim > 1) {
    simbox.boundary_lhs[1] = setBoundaryType(stringparams["boundary_lhs[1]"]);
    simbox.boundary_rhs[1] = setBoundaryType(stringparams["boundary_rhs[1]"]);
    simbox.min[1] = floatparams["boxmin[1]"]/simunits.r.outscale;
    simbox.max[1] = floatparams["boxmax[1]"]/simunits.r.outscale;
  }

  if (ndim == 3) {
    simbox.boundary_lhs[2] = setBoundaryType(stringparams["boundary_lhs[2]"]);
    simbox.boundary_rhs[2] = setBoundaryType(stringparams["boundary_rhs[2]"]);
    simbox.min[2] = floatparams["boxmin[2]"]/simunits.r.outscale;
    simbox.max[2] = floatparams["boxmax[2]"]/simunits.r.outscale;
  }

  for (int k=0; k<ndim; k++) {
    simbox.size[k] = simbox.max[k] - simbox.min[k];
    simbox.half[k] = 0.5*simbox.size[k];
  }


  // Set-up dummy SPH object in order to have valid pointers in N-body object
  hydro = new NullHydrodynamics<ndim>
    (intparams["hydro_forces"], intparams["self_gravity"], floatparams["h_fac"],
     stringparams["gas_eos"], KernelName, sizeof(Particle<ndim>), simunits, simparams);


  // Process all N-body parameters and set-up main N-body objects
  this->ProcessNbodyParameters();


  // Set pointers to external potential field object
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


  // Set important variables for N-body objects
  nbody->Nstar          = intparams["Nstar"];
  nbody->Nstarmax       = intparams["Nstarmax"];
  nbody_single_timestep = intparams["nbody_single_timestep"];
  nbodytree.gpehard     = floatparams["gpehard"];
  nbodytree.gpesoft     = floatparams["gpesoft"];
  //nbody->perturbers     = intparams["perturbers"];
  //if (intparams["sub_systems"] == 1) subsystem->perturbers = intparams["perturbers"];


  // Create Ewald periodic gravity object
  periodicBoundaries = IsAnyBoundaryPeriodic(simbox);
  if (periodicBoundaries) {
    ewaldGravity = true;
    ewald = new Ewald<ndim>
      (simbox, intparams["gr_bhewaldseriesn"], intparams["in"], intparams["nEwaldGrid"],
       floatparams["ewald_mult"], floatparams["ixmin"], floatparams["ixmax"],
       floatparams["EFratio"], timing);
    simbox.PeriodicGravity = true;
  }
  else{
    simbox.PeriodicGravity = false ;
    if (IsAnyBoundaryReflecting(simbox)) {
      ExceptionHandler::getIstance().raise("Error: Reflecting boundaries not supported");
    }
  }


  // Set pointers to timing object
  nbody->timing   = timing;

  sinks = new Sinks<ndim>(0);


  // Flag that we've processed all parameters already
  ParametersProcessed = true;

  return;
}



//=================================================================================================
//  NbodySimulation::PostInitialConditionsSetup
/// ..
//=================================================================================================
template <int ndim>
void NbodySimulation<ndim>::PostInitialConditionsSetup(void)
{
  int i;                            // Particle counter
  int k;                            // Dimension counter

  debug2("[NbodySimulation::PostInitialConditionsSetup]");

  // Set time variables here (for now)
  Noutsnap = 0;
  nresync = 0;
  //tsnapnext = dt_snap;

  // Compute all initial N-body terms
  //-----------------------------------------------------------------------------------------------
  if (nbody->Nstar > 0) {

    // Zero all acceleration terms
    for (i=0; i<nbody->Nstar; i++) {
      for (k=0; k<ndim; k++) nbody->stardata[i].a[k]        = 0.0;
      for (k=0; k<ndim; k++) nbody->stardata[i].adot[k]     = 0.0;
      for (k=0; k<ndim; k++) nbody->stardata[i].a2dot[k]    = 0.0;
      for (k=0; k<ndim; k++) nbody->stardata[i].a3dot[k]    = 0.0;
      for (k=0; k<ndim; k++) nbody->stardata[i].apert[k]    = 0.0;
      for (k=0; k<ndim; k++) nbody->stardata[i].adotpert[k] = 0.0;
      nbody->stardata[i].gpot   = 0.0;
      nbody->stardata[i].gpe    = 0.0;
      nbody->stardata[i].tlast  = t;
      nbody->stardata[i].flags.set(active);
      nbody->stardata[i].level  = level_step;
      nbody->stardata[i].nstep  = 0;
      nbody->stardata[i].nlast  = 0;
      nbody->stardata[i].istar  = i;
      nbody->nbodydata[i]       = &(nbody->stardata[i]);
    }

    nbody->Nnbody = nbody->Nstar;
    if (nbody->nbody_softening == 1) {
      nbody->CalculateDirectSmoothedGravForces(nbody->Nnbody, nbody->nbodydata, simbox, ewald);
    }
    else {
      nbody->CalculateDirectGravForces(nbody->Nnbody, nbody->nbodydata, simbox, ewald);
    }
    //nbody->CalculateDirectGravForces(nbody->Nnbody,nbody->nbodydata);
    nbody->CalculateAllStartupQuantities(nbody->Nnbody, nbody->nbodydata, simbox, ewald);
    for (i=0; i<nbody->Nnbody; i++) {
      if (nbody->nbodydata[i]->flags.check(active)) {
        nbody->extpot->AddExternalPotential(nbody->nbodydata[i]->r,nbody->nbodydata[i]->v,
                                            nbody->nbodydata[i]->a,nbody->nbodydata[i]->adot,
                                            nbody->nbodydata[i]->gpot);
      }
    }

  }

  if (Nlevels == 1) this->ComputeGlobalTimestep();
  else this->ComputeBlockTimesteps();

  // Set particle values for initial step (e.g. r0, v0, a0)
  nbody->EndTimestep(n,nbody->Nstar,t,timestep,nbody->nbodydata);

  this->CalculateDiagnostics();
  this->diag0 = this->diag;

  this->setup = true;

  return;
}



//=================================================================================================
//  NbodySimulation::MainLoop
/// Main N-body simulation integration loop.
//=================================================================================================
template <int ndim>
void NbodySimulation<ndim>::MainLoop(void)
{
  int i;                            // Particle loop counter
  int it;                           // Time-symmetric iteration counter
  int k;                            // Dimension counter

  debug2("[NbodySimulation::MainLoop]");

  // If we are using sub-systems, create N-body tree here
  //-----------------------------------------------------------------------------------------------
  if (nbody->sub_systems == 1) {
    nbody->reset_tree = 1;

    // If we are obliged to re-build the tree, then recompute the grav.
    // potential for all star particles (could be optimised in the future).
    if (Nsteps%nsystembuildstep == 0 || nbody->reset_tree == 1) {

      // Zero all acceleration terms
      for (i=0; i<nbody->Nstar; i++) {
        for (k=0; k<ndim; k++) nbody->stardata[i].a[k]     = 0.0;
        for (k=0; k<ndim; k++) nbody->stardata[i].adot[k]  = 0.0;
        for (k=0; k<ndim; k++) nbody->stardata[i].a2dot[k] = 0.0;
        for (k=0; k<ndim; k++) nbody->stardata[i].a3dot[k] = 0.0;
        nbody->stardata[i].gpot = 0.0;
        nbody->stardata[i].gpe  = 0.0;
        nbody->stardata[i].flags.set(active);
        nbody->nbodydata[i] = &(nbody->stardata[i]);
      }
      nbody->Nnbody = nbody->Nstar;

      // Calculate forces for all star by direct-sum
      if (nbody->nbody_softening == 1) {
        nbody->CalculateDirectSmoothedGravForces(nbody->Nnbody, nbody->nbodydata, simbox, ewald);
      }
      else {
        nbody->CalculateDirectGravForces(nbody->Nnbody, nbody->nbodydata, simbox, ewald);
      }
      //nbody->CalculateDirectGravForces(nbody->Nnbody,nbody->nbodydata);
      nbody->CalculateAllStartupQuantities(nbody->Nnbody, nbody->nbodydata, simbox, ewald);

      // Now create nearest neighbour tree and build any sub-systems from tree
      nbodytree.CreateNbodySystemTree(nbody);
      nbodytree.BuildSubSystems(n,level_max,t,nbody);
    }

    // Find list of perturbers for all sub-systems
    if (nbody->perturbers == 1) nbodytree.FindPerturberLists(nbody);

  }
  //-----------------------------------------------------------------------------------------------



  // Advance time variables
  n = n + 1;
  Nsteps = Nsteps + 1;
  t = t + timestep;
  if (n == nresync) Nblocksteps = Nblocksteps + 1;
  if (n%integration_step == 0) Nfullsteps = Nfullsteps + 1;


  // Advance SPH particles positions and velocities
  nbody->AdvanceParticles(n,nbody->Nnbody,t,timestep,nbody->nbodydata);


  // Compute N-body forces
  //-----------------------------------------------------------------------------------------------
  if (nbody->Nnbody > 0) {

    // Iterate force calculation for time-symmetric schemes.  Otherwise compute forces once.
    //---------------------------------------------------------------------------------------------
    for (it=0; it<nbody->Npec; it++) {

      // Zero all acceleration terms
      for (i=0; i<nbody->Nnbody; i++) {
        if (nbody->nbodydata[i]->flags.check(active)) {
          for (k=0; k<ndim; k++) nbody->nbodydata[i]->a[k]     = 0.0;
          for (k=0; k<ndim; k++) nbody->nbodydata[i]->adot[k]  = 0.0;
          for (k=0; k<ndim; k++) nbody->nbodydata[i]->a2dot[k] = 0.0;
          for (k=0; k<ndim; k++) nbody->nbodydata[i]->a3dot[k] = 0.0;
          nbody->nbodydata[i]->gpot = 0.0;
        }
      }

      // Calculate forces, force derivatives etc.., for active stars/systems
      if (nbody->nbody_softening == 1) {
        nbody->CalculateDirectSmoothedGravForces(nbody->Nnbody, nbody->nbodydata, simbox, ewald);
      }
      else {
        nbody->CalculateDirectGravForces(nbody->Nnbody, nbody->nbodydata, simbox, ewald);
      }
      //nbody->CalculateDirectGravForces(nbody->Nnbody,nbody->nbodydata);

      for (i=0; i<nbody->Nnbody; i++) {
        if (nbody->nbodydata[i]->flags.check(active)) {
          nbody->extpot->AddExternalPotential(nbody->nbodydata[i]->r,nbody->nbodydata[i]->v,
                                              nbody->nbodydata[i]->a,nbody->nbodydata[i]->adot,
                                              nbody->nbodydata[i]->gpot);
        }
      }

      // Calculate correction step for all stars at end of step
      nbody->CorrectionTerms(n, nbody->Nnbody, t, timestep, nbody->nbodydata);

    }
    //---------------------------------------------------------------------------------------------

  }
  //-----------------------------------------------------------------------------------------------


  // Integrate the internal motion of all sub-systems
  if (nbody->sub_systems == 1) {
    for (i=0; i<nbody->Nnbody; i++) {
      if (nbody->nbodydata[i]->Ncomp > 1) {
        // The cast is needed because the function is defined only in SystemParticle, not in
        // NbodyParticle.  The safety of the cast relies on the correctness of the Ncomp value.
        subsystem->IntegrateInternalMotion
         (static_cast<SystemParticle<ndim>* > (nbody->nbodydata[i]), n, t - timestep, t, simbox, ewald);
      }
    }
  }


  // Correct positions of all child stars in any hierarchical systems
  //for (i=0; i<nbody->Nnbody; i++) {
  //  if (nbody->nbodydata[i]->Ncomp > 1) {
  //    // The cast is needed because the function is defined only in SystemParticle, not in
  //    // NbodyParticle.  The safety of the cast relies on the correctness of the Ncomp value.
  //    subsystem->UpdateChildStars(static_cast<SystemParticle<ndim>* > (nbody->nbodydata[i]));
  //  }
  //}

  // Calculate correction terms on perturbing stars due to sub-systems
  //if (nbody->sub_systems == 1 && nbody->perturbers == 1) {
    //nbody->PerturberCorrectionTerms(n,nbody->Nnbody,nbody->nbodydata,timestep);
    //nbody->CorrectionTerms(n,nbody->Nnbody,t,timestep,nbody->nbodydata);
  //}

  // Compute timesteps for all particles
  if (Nlevels == 1) this->ComputeGlobalTimestep();
  else this->ComputeBlockTimesteps();

  // Set all end-of-step variables
  nbody->EndTimestep(n, nbody->Nnbody, t, timestep, nbody->nbodydata);

  return;
}

/*

//=================================================================================================
//  NbodySimulation::ComputeGlobalTimestep
/// Computes global timestep for SPH simulation.  Calculates the minimum
/// timestep for all SPH and N-body particles in the simulation.
//=================================================================================================
template <int ndim>
void NbodySimulation<ndim>::ComputeGlobalTimestep(void)
{
  const FLOAT nbody_mult = simparams->floatparams["nbody_mult"];   // Local copy of multiplier
  DOUBLE dt_min = big_number_dp;                                // Local copy of minimum timestep

  debug2("[NbodySimulation::ComputeGlobalTimestep]");

  //-----------------------------------------------------------------------------------------------
  if (n == nresync) {

    n          = 0;
    level_max  = 0;
    level_step = level_max + integration_step - 1;
    nresync    = integration_step;

    // Compute minimum timestep due to stars/systems
    for (int i=0; i<nbody->Nnbody; i++) {
      if (Nlevels <= 0) {
        nbody->nbodydata[i]->dt = nbody_mult;
      }
      else {
        nbody->nbodydata[i]->dt = nbody->Timestep(nbody->nbodydata[i]);
      }
      dt_min = min(dt_min,nbody->nbodydata[i]->dt);
    }

    // Set all particles to same timestep
    timestep = dt_min;
    for (int i=0; i<nbody->Nnbody; i++) {
      nbody->nbodydata[i]->level = level_max;
      nbody->nbodydata[i]->nstep = pow(2,level_step - nbody->nbodydata[i]->level);
      nbody->nbodydata[i]->nlast = n;
      nbody->nbodydata[i]->dt = timestep;
    }

  }
  //-----------------------------------------------------------------------------------------------

  return;
}



//=================================================================================================
//  NbodySimulation::ComputeBlockTimesteps
/// Compute timesteps for all particles using hierarchical block timesteps.
//=================================================================================================
template <int ndim>
void NbodySimulation<ndim>::ComputeBlockTimesteps(void)
{
  int i;                                     // Particle counter
  int istep;                                 // Aux. variable for changing steps
  int level;                                 // Particle timestep level
  int last_level;                            // Previous timestep level
  int level_max_aux;                         // Aux. maximum level variable
  int level_max_old;                         // Old level_max
  int level_max_nbody = 0;                   // level_max for star particles only
  int level_nbody;                           // local thread var. for N-body level
  int nfactor;                               // Increase/decrease factor of n
  int nstep;                                 // Particle integer step-size
  DOUBLE dt;                                 // Aux. timestep variable
  DOUBLE dt_min = big_number_dp;             // Minimum timestep
  DOUBLE dt_min_aux;                         // Aux. minimum timestep variable
  DOUBLE dt_min_nbody = big_number_dp;       // Maximum N-body particle timestep
  DOUBLE dt_nbody;                           // Aux. minimum N-body timestep

  debug2("[SphSimulation::ComputeBlockTimesteps]");
  CodeTiming::BlockTimer timer = timing->StartNewTimer("BLOCK_TIMESTEPS");


  // Synchronise all timesteps and reconstruct block timestep structure.
  //===============================================================================================
  if (n == nresync) {
    n = 0;
    timestep = big_number_dp;

#pragma omp parallel default(none) shared(dt_min_nbody) private(dt,dt_min_aux,dt_nbody,i)
    {
      // Initialise all timestep and min/max variables
      dt_min_aux = big_number_dp;
      dt_nbody = big_number_dp;

      // Now compute minimum timestep due to stars/systems
#pragma omp for
      for (i=0; i<nbody->Nnbody; i++) {
        dt = nbody->Timestep(nbody->nbodydata[i]);
        dt_min_aux = min(dt_min_aux,dt);
        dt_nbody = min(dt_nbody,dt);
        nbody->nbodydata[i]->dt = dt;
      }

#pragma omp critical
      {
        timestep = min(timestep,dt_min_aux);
        dt_min_nbody = min(dt_min_nbody,dt_nbody);
      }
#pragma omp barrier
    }


    // Calculate new block timestep levels
    level_max = Nlevels - 1;
    level_step = level_max + integration_step - 1;
    dt_max = timestep*powf(2.0, level_max);

    // Calculate the maximum level occupied by all SPH particles
    level_max_nbody = min((int) (invlogetwo*log(dt_max/dt_min_nbody)) + 1, level_max);

    // Populate timestep levels with N-body particles.
    // Ensures that N-body particles occupy levels lower than all SPH particles
    for (i=0; i<nbody->Nnbody; i++) {
      nbody->nbodydata[i]->tlast = t;
      nbody->nbodydata[i]->nlast = n;
      if (nbody->nbodydata[i]->Ncomp > 1) {
        nbody->nbodydata[i]->level = level_max;
        nbody->nbodydata[i]->nstep = pow(2, level_step - level_max);
      }
      else {
        dt = nbody->nbodydata[i]->dt;
        level = min((int) (invlogetwo*log(dt_max/dt)) + 1, level_max);
        level = max(level, 0);
        nbody->nbodydata[i]->level = level;
        nbody->nbodydata[i]->nstep = pow(2, level_step - nbody->nbodydata[i]->level);
      }
    }

    nresync = pow(2, level_step);
    timestep = dt_max / (DOUBLE) nresync;

  }
  // If not resynchronising, check if any SPH/N-body particles need to move
  // up or down timestep levels.
  //===============================================================================================
  else {

    level_max_old = level_max;
    level_max = 0;


#pragma omp parallel default(none) shared(dt_min,dt_min_nbody,level_max_nbody)\
  private(dt,dt_min_aux,dt_nbody,i,istep,last_level,level,level_max_aux,level_nbody,nstep,nfactor)
    {
      dt_min_aux = big_number_dp;
      dt_nbody = big_number_dp;
      level_max_aux = 0;
      level_nbody = 0;

      // Now find all N-body particles at the beginning of a new timestep
      //-------------------------------------------------------------------------------------------
#pragma omp for
      for (i=0; i<nbody->Nnbody; i++) {

        // Skip particles that are not at end of step
        if (nbody->nbodydata[i]->nlast == n) {
          nstep = nbody->nbodydata[i]->nstep;
          last_level = nbody->nbodydata[i]->level;

          // Compute new timestep value and level number
          dt = nbody->Timestep(nbody->nbodydata[i]);
          nbody->nbodydata[i]->dt = dt;
          level = max((int) (invlogetwo*log(dt_max/dt)) + 1, 0);

          // Move up one level (if levels are correctly synchronised) or
          // down several levels if required
          if (level < last_level && last_level > 1 && n%(2*nstep) == 0) {
            nbody->nbodydata[i]->level = last_level - 1;
          }
          else if (level > last_level) {
            nbody->nbodydata[i]->level = level;
          }
          else {
            nbody->nbodydata[i]->level = last_level;
          }

          nbody->nbodydata[i]->tlast = t;
          nbody->nbodydata[i]->nlast = n;
          nbody->nbodydata[i]->nstep = pow(2,level_step - nbody->nbodydata[i]->level);
        }

        // Find maximum level of all N-body particles
        level_nbody = max(level_nbody,nbody->nbodydata[i]->level);
        level_max_aux = max(level_max_aux,nbody->nbodydata[i]->level);
        dt_nbody = min(dt_nbody,nbody->nbodydata[i]->dt);
      }
      //-------------------------------------------------------------------------------------------


#pragma omp critical
      {
        dt_min = min(dt_min,dt_min_aux);
        dt_min_nbody = min(dt_min_nbody,dt_nbody);
        level_max = max(level_max,level_max_aux);
        level_max_nbody = max(level_max_nbody,level_nbody);
      }
#pragma omp barrier

    }


    // For now, don't allow levels to be removed
    //level_max = max(level_max,level_max_old);
    level_step = level_max + integration_step - 1;


    // Reduce all sub-systems to lowest level (if any exist)
    for (i=0; i<nbody->Nnbody; i++) {
      if (nbody->nbodydata[i]->nlast == n && nbody->nbodydata[i]->Ncomp > 1) {
        nbody->nbodydata[i]->nlast = n;
        nbody->nbodydata[i]->level = level_max;
        nbody->nbodydata[i]->nstep = pow(2,level_step - level_max);
      }
    }


    istep = pow(2, level_step - level_max_old + 1);

    // Update all timestep variables if we have removed or added any levels
    //---------------------------------------------------------------------------------------------
    if (level_max > level_max_old) {
      nfactor = pow(2, level_max - level_max_old);
      n *= nfactor;
      for (i=0; i<nbody->Nnbody; i++) nbody->nbodydata[i]->nstep *= nfactor;
      for (i=0; i<nbody->Nnbody; i++) nbody->nbodydata[i]->nlast *= nfactor;
    }
    else if (level_max <= level_max_old - 1 && level_max_old > 1 && n%istep == 0) {
      nfactor = pow(2, level_max_old - level_max);
      n /= nfactor;
      for (i=0; i<nbody->Nnbody; i++) nbody->nbodydata[i]->nlast /= nfactor;
      for (i=0; i<nbody->Nnbody; i++) nbody->nbodydata[i]->nstep /= nfactor;
    }

    nresync = pow(2, level_step);
    timestep = dt_max / (DOUBLE) nresync;

    // Update values of nstep for both SPH and star particles
    for (i=0; i<nbody->Nnbody; i++) {
      if (nbody->nbodydata[i]->nlast == n) nbody->nbodydata[i]->nstep =
        pow(2,level_step - nbody->nbodydata[i]->level);
    }

  }
  //===============================================================================================


#if defined(VERIFY_ALL)
  //VerifyBlockTimesteps();
#endif

  return;
}

*/
