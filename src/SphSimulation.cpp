//=============================================================================
//  SphSimulation.cpp
//  Contains all main functions controlling SPH simulation work-flow.
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
#include "RiemannSolver.h"
#include "SimGhostParticles.cpp"
using namespace std;



// Create template class instances of the main SphSimulation object for
// each dimension used (1, 2 and 3)
template class SphSimulation<1>;
template class SphSimulation<2>;
template class SphSimulation<3>;



//TODO: make this mess more modular (note: initial h computation
//should be done inside the neighbour search)
//=============================================================================
//  SphSimulation::PostGeneration
/// ..
//=============================================================================
template <int ndim>
void SphSimulation<ndim>::PostGeneration(void)
{
  int i;                            // Particle counter
  int k;                            // Dimension counter

  debug2("[SphSimulation::PostGeneration]");

  // Set time variables here (for now)
  Noutsnap = 0;
  tsnapnext = dt_snap;

  // Set initial smoothing lengths and create initial ghost particles
  // --------------------------------------------------------------------------
  if (sph->Nsph > 0) {

    // Set all relevant particle counters
    sph->Nghost = 0;
    sph->Nghostmax = sph->Nsphmax - sph->Nsph;
    sph->Ntot = sph->Nsph;
    for (int i=0; i<sph->Nsph; i++) sph->sphdata[i].active = true;

    sph->InitialSmoothingLengthGuess();
    sphneib->UpdateTree(sph,*simparams);

    sphneib->neibcheck = false;
    sphneib->UpdateAllSphProperties(sph);

    // Search ghost particles
    this->SearchGhostParticles();

    // Update neighbour tree
    sphneib->UpdateTree(sph,*simparams);
  }

  // Compute all initial SPH particle properties (if SPH particles exist)
  // --------------------------------------------------------------------------
  if (sph->Nsph > 0) {

    cout << "Ntot : " << sph->Ntot << endl;
    level_step = 1;

    // Zero accelerations (perhaps here)
    for (i=0; i<sph->Ntot; i++) sph->sphdata[i].active = true;

    // Calculate all SPH properties
    sphneib->UpdateAllSphProperties(sph);

    // Search ghost particles
    this->SearchGhostParticles();

    // Update neighbour tre
    sphneib->UpdateTree(sph,*simparams);
    sphneib->neibcheck = true;
    sphneib->UpdateAllSphProperties(sph);

  }


  // Compute all initial N-body terms
  // --------------------------------------------------------------------------
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
      nbody->stardata[i].nstep = 1;
      nbody->nbodydata[i] = &(nbody->stardata[i]);
    }

    nbody->CalculateDirectGravForces(nbody->Nnbody,nbody->nbodydata);
    nbody->CalculateAllStartupQuantities(nbody->Nnbody,nbody->nbodydata);

  }


  // Compute all initial SPH force terms
  // --------------------------------------------------------------------------
  if (sph->Nsph > 0) {

    // Zero accelerations (here for now)
    for (i=0; i<sph->Ntot; i++) {
      for (k=0; k<ndim; k++) sph->sphdata[i].a[k] = (FLOAT) 0.0;
      for (k=0; k<ndim; k++) sph->sphdata[i].agrav[k] = (FLOAT) 0.0;
      sph->sphdata[i].gpot = (FLOAT) 0.0;
      sph->sphdata[i].dudt = (FLOAT) 0.0;
      sph->sphdata[i].active = true;
      sph->sphdata[i].level = level_step;
      sph->sphdata[i].nstep = 1;
    }

    this->CopySphDataToGhosts();
    sphneib->UpdateTree(sph,*simparams);

    // Calculate SPH gravity and hydro forces, depending on which are activated
    if (sph->hydro_forces == 1 && sph->self_gravity == 1)
      sphneib->UpdateAllSphForces(sph);
    else if (sph->hydro_forces == 1)
      sphneib->UpdateAllSphHydroForces(sph);
    else if (sph->self_gravity == 1)
      sphneib->UpdateAllSphGravForces(sph);

    // Add accelerations
    for (i=0; i<sph->Nsph; i++) {
      sph->sphdata[i].active = false;
      for (k=0; k<ndim; k++)
        sph->sphdata[i].a[k] += sph->sphdata[i].agrav[k];
    }

    this->CopySphDataToGhosts();

  }

  // Set particle values for initial step (e.g. r0, v0, a0)
  sphint->EndTimestep(n,sph->Nsph,sph->sphdata);
  if (simparams->stringparams["gas_eos"] == "energy_eqn")
    uint->EndTimestep(n,sph->Nsph,sph->sphdata);
  nbody->EndTimestep(n,nbody->Nstar,nbody->nbodydata);

  this->CalculateDiagnostics();
  this->diag0 = this->diag;

  this->setup = true;

  return;
}



//=============================================================================
//  SphSimulation::MainLoop
/// Main SPH simulation integration loop.
//=============================================================================
template <int ndim>
void SphSimulation<ndim>::MainLoop(void)
{
  int i;                            // Particle loop counter
  int it;                           // Time-symmetric iteration counter
  int k;                            // Dimension counter

  debug2("[SphSimulation::MainLoop]");

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
  sphint->AdvanceParticles(n,sph->Nsph,sph->sphdata,(FLOAT) timestep);
  if (simparams->stringparams["gas_eos"] == "energy_eqn")
    uint->EnergyIntegration(n,sph->Nsph,sph->sphdata,(FLOAT) timestep);
  nbody->AdvanceParticles(n,nbody->Nstar,nbody->nbodydata,timestep);

  // Check all boundary conditions
  this->CheckBoundaries();

  // --------------------------------------------------------------------------
  if (sph->Nsph > 0) {
    
    // Reorder particles

    // Search ghost particles
    this->SearchGhostParticles();

    // Update neighbour tree
    sphneib->UpdateTree(sph,*simparams);

    // Calculate all SPH properties
    sphneib->UpdateAllSphProperties(sph);

    // Copy properties from original particles to ghost particles
    this->CopySphDataToGhosts();

    // Zero accelerations (perhaps)
    for (i=0; i<sph->Ntot; i++) {
      if (sph->sphdata[i].active) {
        for (k=0; k<ndim; k++) sph->sphdata[i].a[k] = (FLOAT) 0.0;
        for (k=0; k<ndim; k++) sph->sphdata[i].agrav[k] = (FLOAT) 0.0;
        sph->sphdata[i].gpot = (FLOAT) 0.0;
        sph->sphdata[i].dudt = (FLOAT) 0.0;
      }
    }

    // Calculate SPH gravity and hydro forces, depending on which are activated
    if (sph->hydro_forces == 1 && sph->self_gravity == 1)
      sphneib->UpdateAllSphForces(sph);
    else if (sph->hydro_forces == 1)
      sphneib->UpdateAllSphHydroForces(sph);
    else if (sph->self_gravity == 1)
      sphneib->UpdateAllSphGravForces(sph);

    // Compute contribution to grav. accel from stars
    for (i=0; i<sph->Nsph; i++)
      if (sph->sphdata[i].active)
	sph->ComputeStarGravForces(nbody->Nnbody,nbody->nbodydata,
                                   sph->sphdata[i]);

    // Add accelerations
    for (i=0; i<sph->Nsph; i++) {
      for (k=0; k<ndim; k++) 
        sph->sphdata[i].a[k] += sph->sphdata[i].agrav[k];
    }

  }
  // --------------------------------------------------------------------------


  // Compute N-body forces
  // --------------------------------------------------------------------------
  if (nbody->Nstar > 0) {

    // Iterate end-of-step
    // ------------------------------------------------------------------------
	for (it=0; it<nbody->Npec; it++) {

	  cout << "it : " << it << "     Npec : " << nbody->Npec << endl;

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
      nbody->CalculateDirectGravForces(nbody->Nstar,nbody->nbodydata);

      // Calculate correction step for all stars at end of step
      nbody->CorrectionTerms(n,nbody->Nnbody,nbody->nbodydata,timestep);

	}
    // ------------------------------------------------------------------------

  }
  // --------------------------------------------------------------------------


  // Apply correction steps for both particle and energy integration
  sphint->CorrectionTerms(n,sph->Nsph,sph->sphdata,(FLOAT) timestep);
  if (simparams->stringparams["gas_eos"] == "energy_eqn")
    uint->EnergyCorrectionTerms(n,sph->Nsph,sph->sphdata,(FLOAT) timestep);

  // Set all end-of-step variables
  sphint->EndTimestep(n,sph->Nsph,sph->sphdata);
  if (simparams->stringparams["gas_eos"] == "energy_eqn")
    uint->EndTimestep(n,sph->Nsph,sph->sphdata);
  nbody->EndTimestep(n,nbody->Nnbody,nbody->nbodydata);

  return;
}



//=============================================================================
//  SphSimulation::ComputeGlobalTimestep
/// Computes global timestep for SPH simulation.
//=============================================================================
template <int ndim>
void SphSimulation<ndim>::ComputeGlobalTimestep(void)
{
  int i;                            // Particle counter
  DOUBLE dt = big_number_dp;        // Particle timestep
  DOUBLE dt_min = big_number_dp;    // Local copy of minimum timestep

  debug2("[SphSimulation::ComputeGlobalTimestep]");

  // --------------------------------------------------------------------------
  if (n == nresync) {

    n = 0;
    level_max = 0;
    level_step = level_max + integration_step - 1;
    nresync = integration_step;

    // Find minimum timestep from all SPH particles
    // ------------------------------------------------------------------------
#pragma omp parallel default(shared) private(i,dt)
    {
#pragma omp for
      for (i=0; i<sph->Nsph; i++)
        dt = min(dt,sphint->Timestep(sph->sphdata[i],sph->hydro_forces));

      // If integrating energy equation, include energy timestep
      if (simparams->stringparams["gas_eos"] == "energy_eqn") {
#pragma omp for
        for (i=0; i<sph->Nsph; i++)
          dt = min(dt,uint->Timestep(sph->sphdata[i]));
      }

#pragma omp critical
      if (dt < dt_min) dt_min = dt;
    }
    // ------------------------------------------------------------------------


    // Now compute minimum timestep due to stars/systems
    for (i=0; i<nbody->Nstar; i++)
      dt_min = min(dt_min,nbody->Timestep(nbody->nbodydata[i]));

    
    // Set all particles to same timestep
    timestep = dt_min;
    for (i=0; i<sph->Nsph; i++) {
      sph->sphdata[i].level = 0;
      sph->sphdata[i].nstep = pow(2,level_step - sph->sphdata[i].level);
      sph->sphdata[i].dt = timestep;
    }
    for (i=0; i<nbody->Nnbody; i++) {
      nbody->nbodydata[i]->level = 0;
      nbody->nbodydata[i]->nstep = 
        pow(2,level_step - nbody->nbodydata[i]->level);
      nbody->nbodydata[i]->dt = timestep;
    }

  }
  // --------------------------------------------------------------------------

  return;
}



//=============================================================================
//  SphSimulation::ComputeBlockTimesteps
/// ..
//=============================================================================
template <int ndim>
void SphSimulation<ndim>::ComputeBlockTimesteps(void)
{
  int i;                            // Particle counter
  int istep;                        // ??
  int level;                        // Particle timestep level
  int last_level;                   // Previous timestep level
  int level_max_old;                // Old level_max
  int level_max_sph = 0;            // level_max for SPH particles only
  int nstep;                        // ??
  DOUBLE dt;                        // Aux. timestep variable
  DOUBLE dt_max_sph = 0.0;          // Maximum SPH particle timestep

  debug2("[SphSimulation::ComputeBlockTimesteps]");

  timestep = big_number;

  // Synchronise all timesteps and reconstruct block timestep structure.
  // ==========================================================================
  if (n == nresync) {

    n = 0;

    // Find minimum timestep from all SPH particles
    for (i=0; i<sph->Nsph; i++) {
      dt = sphint->Timestep(sph->sphdata[i],sph->hydro_forces);
      if (dt < timestep) timestep = dt;
      if (dt > dt_max_sph) dt_max_sph = dt;
      sph->sphdata[i].dt = dt;
    }
    
    // If integrating energy equation, include energy timestep
    if (sph->gas_eos == "energy_eqn") {
      for (i=0; i<sph->Nsph; i++) {
	dt = uint->Timestep(sph->sphdata[i]);
	if (dt < timestep) timestep = dt;
	sph->sphdata[i].dt = min(sph->sphdata[i].dt,dt);
      }
    }

    // Calculate new block timestep levels
    level_max = Nlevels - 1;
    level_step = level_max + integration_step - 1;
    dt_max = timestep*powf(2.0,level_max);
    nresync = pow(2,level_step);
    timestep = dt_max / (DOUBLE) nresync;
    
    // Calculate the maximum level occupied by all SPH particles
    level_max_sph = max((int) (invlogetwo*log(dt_max/dt_max_sph)) + 1, 0);

    // If enforcing a single SPH timestep, set it here.  Otherwise, populate 
    // the timestep levels with SPH particles.
    if (sph_single_timestep == 1)
      for (i=0; i<sph->Nsph; i++) sph->sphdata[i].level = level_max_sph;
    else {
      for (i=0; i<sph->Nsph; i++) {
	level = min((int) (invlogetwo*log(dt_max/dt)) + 1, level_max);
	level = max(level,0);
	sph->sphdata[i].level = level;
	sph->sphdata[i].dt = 
	  (DOUBLE) pow(2,level_step - sph->sphdata[i].level)*timestep;
      }
    }

  }

  // If not resynchronising, check if any SPH particles need to move up 
  // or down timestep levels
  // ==========================================================================
  else {

    level_max_old = level_max;
    level_max = 0;

    // Find all SPH particles at the beginning of a new timestep
    // ------------------------------------------------------------------------
    for (i=0; i<sph->Nsph; i++) {
      last_level = sph->sphdata[i].level;

      nstep = pow(2,level_step - last_level);
      istep = pow(2,level_step - last_level + 1);

      // Skip particles that are not at end of step
      if (n%nstep == 0) {
	last_level = sph->sphdata[i].level;
	dt = sphint->Timestep(sph->sphdata[i],sph->hydro_forces);
	if (sph->gas_eos == "energy_eqn") 
	  dt = min(dt,uint->Timestep(sph->sphdata[i]));
	sph->sphdata[i].dt = dt;
	level = min((int) (invlogetwo*log(dt_max/dt)) + 1, level_max);
	level = max(level,0);

	// Move up one level (if levels are correctly synchronised) or 
	// down several levels if required
	if (level < last_level && last_level > 1 && n%istep == 0) 
	  sph->sphdata[i].level--;
	else if (level > last_level) {
	  sph->sphdata[i].level = level;
	}
      }

      // Find maximum level of all SPH particles
      level_max_sph = max(level_max_sph,sph->sphdata[i].level);
      level_max = max(level_max,sph->sphdata[i].level);
    }
    // ------------------------------------------------------------------------
      

    // Set fixed SPH timestep level here in case maximum has changed
    if (sph_single_timestep == 1)
      for (i=0; i<sph->Nsph; i++) sph->sphdata[i].level = level_max_sph;


    // Update all timestep variables if we have removed or added any levels
    // ------------------------------------------------------------------------
    if (level_max != level_max_old) {

      // Increase maximum timestep level if correctly synchronised
      istep = pow(2,level_step - level_max_old + 1);
      if (level_max <= level_max_old - 1 && level_max_old > 1 && n%istep == 0)
	level_max = level_max_old - 1;
      else if (level_max == level_max_old)
	level_max = level_max_old;
      level_step = level_max + integration_step - 1;

      // Adjust integer time if levels added or removed
      if (level_max > level_max_old)
	n *= pow(2,level_max - level_max_old);
      else if (level_max < level_max_old)
	n /= pow(2,level_max_old - level_max);

    }
    // ------------------------------------------------------------------------

    nresync = pow(2,level_step);
    timestep = dt_max / (DOUBLE) nresync;

    for (i=0; i<sph->Nsph; i++) sph->sphdata[i].dt = 
      (DOUBLE) pow(2,level_step - sph->sphdata[i].level)*timestep;

  }
  // ==========================================================================

#if defined(VERIFY_ALL)
  //VerifyBlockTimesteps();
#endif

  return;
}




