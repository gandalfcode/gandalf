//=============================================================================
//  NbodySimulation.cpp
//  Contains all main functions controlling SPH simulation work-flow.
//
//  This file is part of GANDALF :
//  Graphical Astrophysics code for N-body Dynamics and Lagrangian Fluids
//  https://github.com/gandalfcode/gandalf
//  Contact : gandalfcode@gmail.com
//
//  Copyright (C) 2013  D. A. Hubber, G Rosotti
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
#include "Debug.h"
#include "InlineFuncs.h"
#include "Simulation.h"
#include "Parameters.h"
#include "Nbody.h"
#include "Sph.h"
#include "Ghosts.h"
using namespace std;



// Create template class instances of the main NbodySimulation object for
// each dimension used (1, 2 and 3)
template class NbodySimulation<1>;
template class NbodySimulation<2>;
template class NbodySimulation<3>;



//=============================================================================
//  NbodySimulation::PostInitialConditionsSetup
/// ..
//=============================================================================
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
  // --------------------------------------------------------------------------
  if (nbody->Nstar > 0) {

    // Zero all acceleration terms
    for (i=0; i<nbody->Nstar; i++) {
      for (k=0; k<ndim; k++) nbody->stardata[i].a[k] = 0.0;
      for (k=0; k<ndim; k++) nbody->stardata[i].adot[k] = 0.0;
      for (k=0; k<ndim; k++) nbody->stardata[i].a2dot[k] = 0.0;
      for (k=0; k<ndim; k++) nbody->stardata[i].a3dot[k] = 0.0;
      for (k=0; k<ndim; k++) nbody->stardata[i].apert[k] = 0.0;
      for (k=0; k<ndim; k++) nbody->stardata[i].adotpert[k] = 0.0;
      nbody->stardata[i].gpot = 0.0;
      nbody->stardata[i].active = true;
      nbody->stardata[i].level = level_step;
      nbody->stardata[i].nstep = 0;
      nbody->stardata[i].nlast = 0;
      nbody->nbodydata[i] = &(nbody->stardata[i]);
    }

    nbody->Nnbody = nbody->Nstar;
    nbody->CalculateDirectGravForces(nbody->Nnbody,nbody->nbodydata);
    nbody->CalculateAllStartupQuantities(nbody->Nnbody,nbody->nbodydata);

  }

  // Set particle values for initial step (e.g. r0, v0, a0)
  nbody->EndTimestep(n,nbody->Nstar,nbody->nbodydata);

  this->CalculateDiagnostics();
  this->diag0 = this->diag;

  this->setup = true;

  return;
}



//=============================================================================
//  NbodySimulation::MainLoop
/// Main N-body simulation integration loop.
//=============================================================================
template <int ndim>
void NbodySimulation<ndim>::MainLoop(void)
{
  int i;                            // Particle loop counter
  int it;                           // Time-symmetric iteration counter
  int k;                            // Dimension counter

  debug2("[NbodySimulation::MainLoop]");

  // If we are using sub-systems, create N-body tree here
  // --------------------------------------------------------------------------
  if (nbody->sub_systems == 1) {
    nbody->reset_tree = 1;

    // If we are obliged to re-build the tree, then recompute the grav.
    // potential for all star particles (could be optimised in the future).
    if (nbody->reset_tree == 1) {
      
      // Zero all acceleration terms
      for (i=0; i<nbody->Nstar; i++) {
        for (k=0; k<ndim; k++) nbody->stardata[i].a[k] = 0.0;
        for (k=0; k<ndim; k++) nbody->stardata[i].adot[k] = 0.0;
        for (k=0; k<ndim; k++) nbody->stardata[i].a2dot[k] = 0.0;
        for (k=0; k<ndim; k++) nbody->stardata[i].a3dot[k] = 0.0;
        nbody->stardata[i].gpot = 0.0;
        nbody->stardata[i].active = true;
        nbody->nbodydata[i] = &(nbody->stardata[i]);
      }
      nbody->Nnbody = nbody->Nstar;
     
      nbody->CalculateDirectGravForces(nbody->Nnbody,nbody->nbodydata);
      nbody->CalculateAllStartupQuantities(nbody->Nnbody,nbody->nbodydata);

      nbodytree.CreateNbodySystemTree(nbody);
      nbodytree.BuildSubSystems(nbody);
    }

    nbodytree.FindPerturberLists(nbody);

  }
  // --------------------------------------------------------------------------


  // Compute timesteps for all particles
  if (Nlevels == 1)
    this->ComputeGlobalTimestep();
  else 
    this->ComputeBlockTimesteps();

  // Advance time variables
  n = n + 1;
  Nsteps = Nsteps + 1;
  t = t + timestep;

  // Advance SPH particles positions and velocities
  nbody->AdvanceParticles(n,nbody->Nnbody,nbody->nbodydata,timestep);

  // Compute N-body forces
  // --------------------------------------------------------------------------
  if (nbody->Nstar > 0) {

    // Iterate end-of-step
    // ------------------------------------------------------------------------
    for (it=0; it<nbody->Npec; it++) {

      // Zero all acceleration terms
      for (i=0; i<nbody->Nnbody; i++) {
        if (nbody->nbodydata[i]->active) {
          for (k=0; k<ndim; k++) nbody->nbodydata[i]->a[k] = 0.0;
          for (k=0; k<ndim; k++) nbody->nbodydata[i]->adot[k] = 0.0;
          for (k=0; k<ndim; k++) nbody->nbodydata[i]->a2dot[k] = 0.0;
          for (k=0; k<ndim; k++) nbody->nbodydata[i]->a3dot[k] = 0.0;
          nbody->nbodydata[i]->gpot = 0.0;
        }
      }

      // Calculate forces, force derivatives etc.., for active stars/systems
      nbody->CalculateDirectGravForces(nbody->Nnbody,nbody->nbodydata);
      
      // Calculate correction step for all stars at end of step
      nbody->CorrectionTerms(n,nbody->Nnbody,nbody->nbodydata,timestep);

    }
    // ------------------------------------------------------------------------

  }
  // --------------------------------------------------------------------------


  // Now loop over children and, if they are systems, integrate
  // their internal motion
  // --------------------------------------------------------------------------
  for (i=0; i<nbody->Nnbody; i++) {
    if (nbody->nbodydata[i]->Ncomp > 1) {
      // The cast is needed because the function is defined only in
      // SystemParticle, not in NbodyParticle.  
      // The safety of the cast relies on the correctness of the Ncomp value
      subsystem->IntegrateInternalMotion(static_cast<SystemParticle<ndim>* > 
                                         (nbody->nbodydata[i]), n, 
					 timestep, timestep);
    }
  }

  // Calculate correction terms on perturbing stars due to sub-systems
  if (nbody->sub_systems == 1 && nbody->perturbers == 1) {
    nbody->PerturberCorrectionTerms(n,nbody->Nnbody,nbody->nbodydata,timestep);
    nbody->CorrectionTerms(n,nbody->Nnbody,nbody->nbodydata,timestep);
  }

  // Update properties of child stars in sub-systems to correctly 
  // match updates to the parent system particle
  if (nbody->sub_systems == 1) {
    for (i=0; i<nbody->Nnbody; i++) {
      if (nbody->nbodydata[i]->Ncomp > 1)
	subsystem->UpdateChildStars(static_cast<SystemParticle<ndim>* > 
				    (nbody->nbodydata[i]),  
				    n, timestep, timestep);
    }
  }

  // Set all end-of-step variables
  nbody->EndTimestep(n,nbody->Nnbody,nbody->nbodydata);

  return;
}



//=============================================================================
//  NbodySimulation::ComputeGlobalTimestep
/// Computes global timestep for SPH simulation.  Calculates the minimum 
/// timestep for all SPH and N-body particles in the simulation.
//=============================================================================
template <int ndim>
void NbodySimulation<ndim>::ComputeGlobalTimestep(void)
{
  int i;                            // Particle counter
  DOUBLE dt = big_number_dp;        // Particle timestep
  DOUBLE dt_min = big_number_dp;    // Local copy of minimum timestep

  debug2("[NbodySimulation::ComputeGlobalTimestep]");

  // --------------------------------------------------------------------------
  if (n == nresync) {

    n = 0;
    level_max = 0;
    level_step = level_max + integration_step - 1;
    nresync = integration_step;

    // Compute minimum timestep due to stars/systems
    for (i=0; i<nbody->Nnbody; i++)
      dt_min = min(dt_min,nbody->Timestep(nbody->nbodydata[i]));

    // Set all particles to same timestep
    timestep = dt_min;
    for (i=0; i<nbody->Nnbody; i++) {
      nbody->nbodydata[i]->level = 0;
      nbody->nbodydata[i]->nstep = 
        pow(2,level_step - nbody->nbodydata[i]->level);
      nbody->nbodydata[i]->nlast = n;
      nbody->nbodydata[i]->dt = timestep;
    }

  }
  // --------------------------------------------------------------------------

  return;
}



//=============================================================================
//  NbodySimulation::ComputeBlockTimesteps
/// Compute timesteps for all particles using hierarchical block timesteps.
//=============================================================================
template <int ndim>
void NbodySimulation<ndim>::ComputeBlockTimesteps(void)
{
  int i;                                // Particle counter
  int istep;                            // ??
  int level;                            // Particle timestep level
  int last_level;                       // Previous timestep level
  int level_max_old;                    // Old level_max
  int level_max_nbody = 0;              // level_max for SPH particles only
  int nstep;                            // ??
  int nfactor;                          // ??
  DOUBLE dt;                            // Aux. timestep variable
  DOUBLE dt_min_nbody = big_number_dp;  // Maximum N-body particle timestep

  int *ninlevel;

  debug2("[NbodySimulation::ComputeBlockTimesteps]");

  // Synchronise all timesteps and reconstruct block timestep structure.
  // ==========================================================================
  if (n == nresync) {

    n = 0;
    timestep = big_number_dp;
    for (i=0; i<nbody->Nnbody; i++) nbody->nbodydata[i]->dt = big_number_dp;

    // Compute minimum timestep due to stars/systems
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
    level_max_nbody = min((int) (invlogetwo*log(dt_max/dt_min_nbody)) + 1, 
			  level_max);

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
  // ==========================================================================
  else {

    level_max_old = level_max;
    level_max = 0;


    // Find all N-body particles at the beginning of a new timestep
    // ------------------------------------------------------------------------
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
    // ------------------------------------------------------------------------
      

    // For now, don't allow levels to be removed
    level_max = max(level_max,level_max_old);
    level_step = level_max + integration_step - 1;

    for (i=0; i<nbody->Nnbody; i++) {
      if (nbody->nbodydata[i]->nlast == n) nbody->nbodydata[i]->nstep = 
	pow(2,level_step - nbody->nbodydata[i]->level);
    }

    // Update all timestep variables if we have removed or added any levels
    // ------------------------------------------------------------------------
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
	for (i=0; i<nbody->Nnbody; i++) nbody->nbodydata[i]->nlast *= nfactor;
	for (i=0; i<nbody->Nnbody; i++) nbody->nbodydata[i]->nstep *= nfactor;
      }
      else if (level_max < level_max_old) {
	nfactor = pow(2,level_max_old - level_max);
	n /= nfactor;
	for (i=0; i<nbody->Nnbody; i++) nbody->nbodydata[i]->nlast /= nfactor;
	for (i=0; i<nbody->Nnbody; i++) nbody->nbodydata[i]->nstep /= nfactor;
      }

    }
    // ------------------------------------------------------------------------

    nresync = pow(2,level_step);
    timestep = dt_max / (DOUBLE) nresync;

  }
  // ==========================================================================


#if defined(VERIFY_ALL)
  //VerifyBlockTimesteps();
#endif

  return;

  // Some validations
  // --------------------------------------------------------------------------
  ninlevel = new int[level_max+1];
  for (int l=0; l<level_max+1; l++) ninlevel[l] = 0;
  for (i=0; i<nbody->Nnbody; i++) ninlevel[nbody->nbodydata[i]->level]++;
  cout << "Level occupancy; timestep : " << timestep << endl;
  for (int l=0; l<level_max+1; l++) 
    cout << "level : " << l << "     N : " << ninlevel[l] << endl;
  delete[] ninlevel;
  cin >> i;
  return;
}
