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
#include "Ghosts.h"
#include "Sinks.h"
using namespace std;



// Create template class instances of the main SphSimulation object for
// each dimension used (1, 2 and 3)
template class SphSimulation<1>;
template class SphSimulation<2>;
template class SphSimulation<3>;



//TODO: make this mess more modular (note: initial h computation
//should be done inside the neighbour search)
//=============================================================================
//  SphSimulation::PostInitialConditionsSetup
/// ..
//=============================================================================
template <int ndim>
void SphSimulation<ndim>::PostInitialConditionsSetup(void)
{
  int i;                            // Particle counter
  int k;                            // Dimension counter

  debug2("[SphSimulation::PostInitialConditionsSetup]");


  // Perform initial MPI decomposition
  //---------------------------------------------------------------------------
#ifdef MPI_PARALLEL
  mpicontrol.CreateInitialDomainDecomposition(sph,nbody,simparams,simbox);
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Abort(MPI_COMM_WORLD,0);
#endif


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
    for (i=0; i<sph->Nsph; i++) sph->sphdata[i].active = true;
    for (i=0; i<sph->Nsph; i++) sph->sphintdata[i].part = &(sph->sphdata[i]);

    // Set initial artificial viscosity alpha values
    if (sph->time_dependent_avisc == 1)
      for (i=0; i<sph->Nsph; i++) {
        sph->sphdata[i].alpha = sph->alpha_visc_min;
      }
    else
      for (i=0; i<sph->Nsph; i++) sph->sphdata[i].alpha = sph->alpha_visc;

    // Compute mean mass
    sph->mmean = 0.0;
    for (i=0; i<sph->Nsph; i++) sph->mmean += sph->sphdata[i].m;
    sph->mmean /= (FLOAT) sph->Nsph;

    sph->InitialSmoothingLengthGuess();
    sphneib->BuildTree(sph,*simparams);

    sphneib->neibcheck = false;
    sphneib->UpdateAllSphProperties(sph,nbody);

    // Search ghost particles
    ghosts->SearchGhostParticles(simbox,sph);

    // Update neighbour tree
    sphneib->BuildTree(sph,*simparams);

    level_step = 1;

    // Zero accelerations
    for (i=0; i<sph->Ntot; i++) sph->sphdata[i].active = true;

    // Calculate all SPH properties
    sphneib->UpdateAllSphProperties(sph,nbody);

    // Search ghost particles
    ghosts->SearchGhostParticles(simbox,sph);

    // Update neighbour tre
    sphneib->BuildTree(sph,*simparams);
    sphneib->neibcheck = true;
    sphneib->UpdateAllSphProperties(sph,nbody);

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
      nbody->stardata[i].level = 0; //level_step;
      nbody->stardata[i].nstep = 0;
      nbody->stardata[i].nlast = 0;
      nbody->nbodydata[i] = &(nbody->stardata[i]);
    }

    nbody->Nnbody = nbody->Nstar;
    nbody->CalculateDirectGravForces(nbody->Nnbody,nbody->nbodydata);
    if (sph->self_gravity == 1)
      nbody->CalculateDirectSPHForces(nbody->Nnbody,sph->Nsph,
				      sph->sphdata,nbody->nbodydata);
    nbody->CalculateAllStartupQuantities(nbody->Nnbody,nbody->nbodydata);

  }


  // Compute all initial SPH force terms
  //---------------------------------------------------------------------------
  if (sph->Nsph > 0) {

    // Zero accelerations (here for now)
    for (i=0; i<sph->Ntot; i++) {
      for (k=0; k<ndim; k++) sph->sphdata[i].a[k] = (FLOAT) 0.0;
      for (k=0; k<ndim; k++) sph->sphdata[i].agrav[k] = (FLOAT) 0.0;
      sph->sphdata[i].gpot = (FLOAT) 0.0;
      sph->sphdata[i].dudt = (FLOAT) 0.0;
      sph->sphdata[i].dalphadt = (FLOAT) 0.0;
      sph->sphdata[i].active = true;
      sph->sphdata[i].level = 0;
      sph->sphintdata[i].nstep = 0;
      sph->sphintdata[i].nlast = 0;
    }

    ghosts->CopySphDataToGhosts(sph);
    sphneib->BuildTree(sph,*simparams);

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
      sph->sphdata[i].active = false;
      for (k=0; k<ndim; k++)
        sph->sphdata[i].a[k] += sph->sphdata[i].agrav[k];
    }

    ghosts->CopySphDataToGhosts(sph);

  }

  // Set particle values for initial step (e.g. r0, v0, a0)
  if (simparams->stringparams["gas_eos"] == "energy_eqn")
    uint->EndTimestep(n,sph->Nsph,sph->sphintdata);
  sphint->EndTimestep(n,sph->Nsph,sph->sphintdata);
  nbody->EndTimestep(n,nbody->Nstar,nbody->nbodydata);

  this->CalculateDiagnostics();
  this->OutputDiagnostics();
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
  int activecount;                  // Flag if we need to recompute particles
  int i;                            // Particle loop counter
  int it;                           // Time-symmetric iteration counter
  int k;                            // Dimension counter

  debug2("[SphSimulation::MainLoop]");


  // Compute timesteps for all particles
  //---------------------------------------------------------------------------
  // MPI : Currently, MPI commands to transmit timestep information will need
  //       to be inside these routines, unless they are partly re-written.
  //---------------------------------------------------------------------------
  if (Nlevels == 1)
    this->ComputeGlobalTimestep();
  else 
    this->ComputeBlockTimesteps();

  // Advance time variables
  n = n + 1;
  Nsteps = Nsteps + 1;
  t = t + timestep;
  if (n == nresync) Nblocksteps = Nblocksteps + 1;

  // Advance SPH particles positions and velocities
  sphint->AdvanceParticles(n,sph->Nsph,sph->sphintdata,(FLOAT) timestep);
  if (simparams->stringparams["gas_eos"] == "energy_eqn")
    uint->EnergyIntegration(n,sph->Nsph,sph->sphintdata,(FLOAT) timestep);
  nbody->AdvanceParticles(n,nbody->Nnbody,nbody->nbodydata,timestep);

  // Check all boundary conditions
  // (DAVID : Move this function to sphint and create an analagous one for N-body)
  // (Also, only check this on tree-build steps)
  ghosts->CheckBoundaries(simbox,sph);


  //---------------------------------------------------------------------------
  // MPI : On tree re-build step, determine load balancing for all MPI nodes.
  //       (How is this done?  All computed on root node??)
  //       Send/receive particles to their new nodes.
  //       Compute and transmit all bounding boxes (e.g. all particles, active
  //       particles, h-extent, ghosts, etc..) to all other MPI nodes
  //---------------------------------------------------------------------------


  // Compute all SPH quantities
  //---------------------------------------------------------------------------
  if (sph->Nsph > 0) {
    
    // Search for new ghost particles and create on local processor
    ghosts->SearchGhostParticles(simbox,sph);


    // Reorder particles to tree-walk order (not implemented yet)


    // Rebuild or update local neighbour and gravity tree
    sphneib->BuildTree(sph,*simparams);


    //-------------------------------------------------------------------------
    // MPI : Walk local tree to determine minimum tree to be sent to all other
    //       MPI nodes.  Pack and send minimum sub-tree, along with the
    //       MPI-ghost particles contained in leaf cells of the sub-tree.
    //-------------------------------------------------------------------------


    // Iterate if we need to immediately change SPH particle timesteps
    // (e.g. due to feedback, or sudden change in neighbour timesteps)
    //-------------------------------------------------------------------------
    do {

      // Update cells containing active particles
      sphneib->UpdateActiveParticleCounters(sph);

      // Calculate all SPH properties
      sphneib->UpdateAllSphProperties(sph,nbody);
      

      //-----------------------------------------------------------------------
      // MPI : Transmit updated particle properties from parent node to
      //       other MPI nodes for MPI-ghost particles.
      //-----------------------------------------------------------------------


      // Copy properties from original particles to ghost particles
      ghosts->CopySphDataToGhosts(sph);
      
      // Zero accelerations
#pragma parallel for default(none) private(k) shared(sph)
      for (i=0; i<sph->Ntot; i++) {
        if (sph->sphdata[i].active) {
          for (k=0; k<ndim; k++) sph->sphdata[i].a[k] = (FLOAT) 0.0;
          for (k=0; k<ndim; k++) sph->sphdata[i].agrav[k] = (FLOAT) 0.0;
          sph->sphdata[i].gpot = (FLOAT) 0.0;
          sph->sphdata[i].gpe = (FLOAT) 0.0;
          sph->sphdata[i].dudt = (FLOAT) 0.0;
          sph->sphdata[i].levelneib = 0;
        }
      }
      
      // Compute SPH gravity and hydro forces, depending on which are activated
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
      
      // Compute additional terms now accelerations and other derivatives 
      // have been computed for active particles
      for (i=0; i<sph->Nsph; i++) {
        if (sph->sphdata[i].active) {
          for (k=0; k<ndim; k++)
            sph->sphdata[i].a[k] += sph->sphdata[i].agrav[k];
          sph->sphdata[i].dalphadt = 0.1*sph->sphdata[i].sound*
            (sph->alpha_visc_min - sph->sphdata[i].alpha)*
            sph->sphdata[i].invh + max(sph->sphdata[i].div_v,0.0)*
            (sph->alpha_visc - sph->sphdata[i].alpha);
        }
      }

      // Check if all neighbouring timesteps are acceptable
      if (Nlevels > 1)
        activecount = sphint->CheckTimesteps(level_diff_max,n,
                                             sph->Nsph,sph->sphintdata);
      else activecount = 0;
      //activecount = 0;

    } while (activecount > 0);
    //-------------------------------------------------------------------------

  }
  //---------------------------------------------------------------------------


  // Compute N-body forces
  //---------------------------------------------------------------------------
  if (nbody->Nnbody > 0) {

    // Iterate for P(EC)^n schemes
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
      if (sink_particles == 1) {
	for (i=0; i<sinks.Nsink; i++) {
          if (sinks.sink[i].star->active) {
	    for (k=0; k<ndim; k++) sinks.sink[i].fhydro[k] = 0.0;
	  }
	}
      }
	  

      // Calculate forces, force derivatives etc.., for active stars/systems
      nbody->CalculateDirectGravForces(nbody->Nnbody,nbody->nbodydata);

      if (sph->self_gravity == 1)
        nbody->CalculateDirectSPHForces(nbody->Nnbody,sph->Nsph,
                                        sph->sphdata,nbody->nbodydata);

      // Calculate correction step for all stars at end of step
      nbody->CorrectionTerms(n,nbody->Nnbody,nbody->nbodydata,timestep);

    }
    //-------------------------------------------------------------------------

  }
  //---------------------------------------------------------------------------


  // Compute correction steps for all SPH particles
  if (sph->Nsph > 0) {
    sphint->CorrectionTerms(n,sph->Nsph,sph->sphintdata,(FLOAT) timestep);
    if (simparams->stringparams["gas_eos"] == "energy_eqn")
      uint->EnergyCorrectionTerms(n,sph->Nsph,sph->sphintdata,(FLOAT) timestep);
  }


  // Search for new sink particles (if activated)
  if (sink_particles == 1) {
    if (sinks.create_sinks == 1) sinks.SearchForNewSinkParticles(n,sph,nbody);
    if (sinks.Nsink > 0) sinks.AccreteMassToSinks(sph,nbody,n,timestep);
  }


  // End-step terms for all SPH particles
  if (sph->Nsph > 0) {
    if (simparams->stringparams["gas_eos"] == "energy_eqn")
      uint->EndTimestep(n,sph->Nsph,sph->sphintdata);
    sphint->EndTimestep(n,sph->Nsph,sph->sphintdata);
  }

  // End-step terms for all star particles
  if (nbody->Nstar > 0)
    nbody->EndTimestep(n,nbody->Nnbody,nbody->nbodydata);

  return;
}



//=============================================================================
//  SphSimulation::ComputeGlobalTimestep
/// Computes global timestep for SPH simulation.  Calculates the minimum 
/// timestep for all SPH and N-body particles in the simulation.
//=============================================================================
template <int ndim>
void SphSimulation<ndim>::ComputeGlobalTimestep(void)
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
        sph->sphdata[i].dt = sphint->Timestep(sph->sphdata[i],
                                              sph->hydro_forces);
        dt = min(dt,sph->sphdata[i].dt);
      }
      
      // If integrating energy equation, include energy timestep
      if (simparams->stringparams["gas_eos"] == "energy_eqn") {
#pragma omp for
        for (i=0; i<sph->Nsph; i++) {
          sph->sphdata[i].dt = min(sph->sphdata[i].dt,
                                   uint->Timestep(sph->sphdata[i]));
          dt = min(dt,sph->sphdata[i].dt);
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
      sph->sphdata[i].level = 0;
      sph->sphdata[i].levelneib = 0;
      sph->sphdata[i].dt = timestep;
      sph->sphintdata[i].nstep = pow(2,level_step - sph->sphdata[i].level);
      sph->sphintdata[i].nlast = n;

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

  return;
}



//=============================================================================
//  SphSimulation::ComputeBlockTimesteps
/// Compute timesteps for all particles using hierarchical block timesteps.
//=============================================================================
template <int ndim>
void SphSimulation<ndim>::ComputeBlockTimesteps(void)
{
  int i;                                // Particle counter
  int imin=-1;                          // id of particle with minimum timestep
  int istep;                            // ??
  int level;                            // Particle timestep level
  int last_level;                       // Previous timestep level
  int level_max_old;                    // Old level_max
  int level_max_sph = 0;                // level_max for SPH particles only
  int level_min_sph = 9999999;          // level_min for SPH particles
  int level_max_nbody = 0;              // level_max for star particles only
  int nstep;                            // ??
  int nfactor;                          // ??
  DOUBLE dt;                            // Aux. timestep variable
  DOUBLE dt_min_sph = big_number_dp;    // Maximum SPH particle timestep
  DOUBLE dt_min_nbody = big_number_dp;  // Maximum N-body particle timestep

  int *ninlevel;

  debug2("[SphSimulation::ComputeBlockTimesteps]");

  // Synchronise all timesteps and reconstruct block timestep structure.
  //===========================================================================
  if (n == nresync) {

    n = 0;
    timestep = big_number_dp;
    for (i=0; i<sph->Nsph; i++) sph->sphdata[i].dt = big_number_dp;
    for (i=0; i<nbody->Nnbody; i++) nbody->nbodydata[i]->dt = big_number_dp;

    // If integrating energy equation, calculate the explicit energy timestep
    if (sph->gas_eos == "energy_eqn") {
      for (i=0; i<sph->Nsph; i++)
        sph->sphdata[i].dt = uint->Timestep(sph->sphdata[i]);
    }

    // Find minimum timestep from all SPH particles
    for (i=0; i<sph->Nsph; i++) {
      dt = min(sph->sphdata[i].dt,
               sphint->Timestep(sph->sphdata[i],sph->hydro_forces));
      if (dt < timestep) timestep = dt;
      if (dt < dt_min_sph) imin = i;
      if (dt < dt_min_sph) dt_min_sph = dt;
      sph->sphdata[i].dt = dt;
      
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
        sph->sphdata[i].level = level_max_sph;
        sph->sphdata[i].levelneib = level_max_sph;
        sph->sphintdata[i].nlast = n;
        sph->sphintdata[i].nstep = pow(2,level_step - sph->sphdata[i].level);
        level_min_sph = min(level_min_sph,sph->sphdata[i].level);
      }
    else {
      for (i=0; i<sph->Nsph; i++) {
        dt = sph->sphdata[i].dt;
        level = min((int) (invlogetwo*log(dt_max/dt)) + 1, level_max);
        level = max(level,0);
        sph->sphdata[i].level = level;
        sph->sphdata[i].levelneib = level;
        sph->sphintdata[i].nlast = n;
        sph->sphintdata[i].nstep = pow(2,level_step - sph->sphdata[i].level);
        level_min_sph = min(level_min_sph,sph->sphdata[i].level);
      }
    }

    // Populate timestep levels with N-body particles
    for (i=0; i<nbody->Nnbody; i++) {
      dt = nbody->nbodydata[i]->dt;
      level = min((int) (invlogetwo*log(dt_max/dt)) + 1, level_max);
      level = max(level,0);
      nbody->nbodydata[i]->level = max(level,level_max_sph);
      //nbody->nbodydata[i]->level = max(level,level_min_sph);
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

      // Skip particles that are not at end of step
      if (sph->sphintdata[i].nlast == n) {
        nstep = sph->sphintdata[i].nstep;
        last_level = sph->sphdata[i].level;
	
        // Compute new timestep value and level number
        dt = sphint->Timestep(sph->sphdata[i],sph->hydro_forces);
        if (sph->gas_eos == "energy_eqn")
          dt = min(dt,uint->Timestep(sph->sphdata[i]));
        sph->sphdata[i].dt = dt;
        level = max((int) (invlogetwo*log(dt_max/dt)) + 1, 0);
        level = max(level,sph->sphdata[i].levelneib - level_diff_max);

        // Move up one level (if levels are correctly synchronised) or
        // down several levels if required
        if (level < last_level && last_level > 1 && n%(2*nstep) == 0)
          sph->sphdata[i].level = last_level - 1;
        else if (level > last_level)
          sph->sphdata[i].level = level;
        else
          sph->sphdata[i].level = last_level;

        sph->sphintdata[i].nlast = n;
        sph->sphintdata[i].nstep = pow(2,level_step - sph->sphdata[i].level);
      }

      // Find maximum level of all SPH particles
      level_max_sph = max(level_max_sph,sph->sphdata[i].level);
      level_min_sph = min(level_min_sph,sph->sphdata[i].level);
      level_max = max(level_max,sph->sphdata[i].level);
      if (sph->sphdata[i].dt < dt_min_sph) imin = i;
      dt_min_sph = min(dt_min_sph,sph->sphdata[i].dt);
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
        level = max((int) (invlogetwo*log(dt_max/dt)) + 1, 0);
        level = max(level,level_max_sph);
        //level = max(level,level_min_sph);

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
      dt_min_nbody = min(dt_min_nbody,nbody->nbodydata[i]->dt);
    }
    //-------------------------------------------------------------------------
      

    // Set fixed SPH timestep level here in case maximum has changed
    if (sph_single_timestep == 1) {
      for (i=0; i<sph->Nsph; i++) {
        if (sph->sphintdata[i].nlast == n)  {
          sph->sphdata[i].level = level_max_sph;
          sph->sphintdata[i].nstep = pow(2,level_step - sph->sphdata[i].level);
        }
      }
    }

    // For now, don't allow levels to be removed
    level_max = max(level_max,level_max_old);
    level_step = level_max + integration_step - 1;

    for (i=0; i<sph->Nsph; i++) {
      if (sph->sphintdata[i].nlast == n)
        sph->sphintdata[i].nstep = pow(2,level_step - sph->sphdata[i].level);
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
        for (i=0; i<sph->Nsph; i++) sph->sphintdata[i].nlast *= nfactor;
        for (i=0; i<sph->Nsph; i++) sph->sphintdata[i].nstep *= nfactor;
        for (i=0; i<nbody->Nnbody; i++) nbody->nbodydata[i]->nlast *= nfactor;
        for (i=0; i<nbody->Nnbody; i++) nbody->nbodydata[i]->nstep *= nfactor;
      }
      else if (level_max < level_max_old) {
        nfactor = pow(2,level_max_old - level_max);
        n /= nfactor;
        for (i=0; i<sph->Nsph; i++) sph->sphintdata[i].nlast /= nfactor;
        for (i=0; i<sph->Nsph; i++) sph->sphintdata[i].nstep /= nfactor;
        for (i=0; i<nbody->Nnbody; i++) nbody->nbodydata[i]->nlast /= nfactor;
        for (i=0; i<nbody->Nnbody; i++) nbody->nbodydata[i]->nstep /= nfactor;
      }

    }
    //-------------------------------------------------------------------------

    nresync = pow(2,level_step);
    timestep = dt_max / (DOUBLE) nresync;

  }
  //===========================================================================


#if defined(VERIFY_ALL)
  //VerifyBlockTimesteps();
#endif

  return;

  // Some validations
  //---------------------------------------------------------------------------
  ninlevel = new int[level_max+1];

  cout << "-----------------------------------------------------" << endl;
  cout << "Checking timesteps : " << level_max << "   " << level_max_sph << "    " << level_max_nbody << endl;
  cout << "dt_min_sph : " << dt_min_sph << "    dt_min_nbody : " << dt_min_nbody << "    timestep : " << timestep << endl;
  cout << "imin : " << imin << "    " << sph->sphdata[imin].dt << "     " 
       << sph->sphdata[imin].h << "     " 
       << sqrt(DotProduct(sph->sphdata[imin].a,sph->sphdata[imin].a,ndim)) 
       << endl;
  for (int l=0; l<level_max+1; l++) ninlevel[l] = 0;
  for (i=0; i<sph->Nsph; i++) ninlevel[sph->sphdata[i].level]++;
  cout << "SPH level occupancy" << endl;
  for (int l=0; l<level_max+1; l++) 
    cout << "level : " << l << "     N : " << ninlevel[l] << endl;

  for (int l=0; l<level_max+1; l++) ninlevel[l] = 0;
  for (i=0; i<nbody->Nstar; i++) ninlevel[nbody->nbodydata[i]->level]++;
  cout << "N-body level occupancy" << endl;
  for (int l=0; l<level_max+1; l++)
    cout << "level : " << l << "     N : " << ninlevel[l] << endl;

  for (int l=0; l<level_max+1; l++) {
    if (ninlevel[l] > 0 && l < level_max_sph) {
      cout << "Something going wrong with timesteps" << endl;
      exit(0);
    }
  }

  delete[] ninlevel;

  return;
}
