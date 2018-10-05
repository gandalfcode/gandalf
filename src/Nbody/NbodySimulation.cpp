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
      nbody->stardata[i].level  = level_step;
      nbody->stardata[i].istar  = i;
      nbody->stardata[i].flags.set(active);
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
  nbody->EndTimestep(level_step, n, nbody->Nstar, t, timestep, nbody->nbodydata);

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
  nbody->AdvanceParticles(level_step, n, nbody->Nnbody, t, timestep, nbody->nbodydata);


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
      nbody->CorrectionTerms(level_step, n, nbody->Nnbody, t, timestep, nbody->nbodydata);

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
         (static_cast<SystemParticle<ndim>* > (nbody->nbodydata[i]),
          level_step, n, t - timestep, t, simbox, ewald);
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
  nbody->EndTimestep(level_step, n, nbody->Nnbody, t, timestep, nbody->nbodydata);

  return;
}
