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
      for (i=0; i<sph->Nsph; i++) sph->sphdata[i].alpha = sph->alpha_visc_min;
    else
      for (i=0; i<sph->Nsph; i++) sph->sphdata[i].alpha = sph->alpha_visc;

    // Compute mean mass
    sph->mmean = 0.0;
    for (i=0; i<sph->Nsph; i++) sph->mmean += sph->sphdata[i].m;
    sph->mmean /= (FLOAT) sph->Nsph;

    // If the smoothing lengths have not been provided beforehand, then
    // calculate the initial values here
    sphneib->neibcheck = false;
    if (!this->initial_h_provided) {
      sph->InitialSmoothingLengthGuess();
      sphneib->BuildTree(rebuild_tree,0,ntreebuildstep,ntreestockstep,timestep,sph);

      sphneib->UpdateAllSphProperties(sph,nbody);
    }

#ifdef MPI_PARALLEL
    mpicontrol.UpdateAllBoundingBoxes(sph->Nsph, sph->sphdata, sph->kernp);
#endif

    // Search ghost particles
    LocalGhosts->SearchGhostParticles(0.0,simbox,sph);
#ifdef MPI_PARALLEL
    MpiGhosts->SearchGhostParticles(0.0,simbox,sph);
#endif

    // Update neighbour tree
    rebuild_tree = true;
    sphneib->BuildTree(rebuild_tree,0,ntreebuildstep,ntreestockstep,timestep,sph);

    level_step = 1;

    // Zero accelerations
    for (i=0; i<sph->Nsph; i++) sph->sphdata[i].active = true;

    // Calculate all SPH properties
    sphneib->UpdateAllSphProperties(sph,nbody);

#ifdef MPI_PARALLEL
    mpicontrol.UpdateAllBoundingBoxes(sph->Nsph, sph->sphdata, sph->kernp);
#endif

    // Search ghost particles
    LocalGhosts->SearchGhostParticles(0.0,simbox,sph);
#ifdef MPI_PARALLEL
    MpiGhosts->SearchGhostParticles(0.0,simbox,sph);
#endif

    // Update neighbour tre
    rebuild_tree = true;
    sphneib->BuildTree(rebuild_tree,0,ntreebuildstep,ntreestockstep,timestep,sph);
    sphneib->neibcheck = true;
    //sphneib->UpdateAllSphProperties(sph,nbody);

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

  }



  // Compute all initial SPH force terms
  //---------------------------------------------------------------------------
  if (sph->Nsph > 0) {

    // Zero accelerations (here for now)
    for (i=0; i<sph->Ntot; i++) {
      sph->sphdata[i].level = 0;
      sph->sphintdata[i].nstep = 0;
      sph->sphintdata[i].nlast = 0;
      sph->sphdata[i].active = false;
    }
    for (i=0; i<sph->Nsph; i++) sph->sphdata[i].active = true;

    LocalGhosts->CopySphDataToGhosts(simbox,sph);
#ifdef MPI_PARALLEL
    MpiGhosts->CopySphDataToGhosts(simbox,sph);
#endif
    sphneib->BuildTree(rebuild_tree,n,ntreebuildstep,ntreestockstep,timestep,sph);

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
    //for (i=0; i<sph->Nsph; i++) {
    //  sph->sphdata[i].active = false;
    //  for (k=0; k<ndim; k++)
    //    sph->sphdata[i].a[k] += sph->sphdata[i].agrav[k];
    //}

    LocalGhosts->CopySphDataToGhosts(simbox,sph);
#ifdef MPI_PARALLEL
    MpiGhosts->CopySphDataToGhosts(simbox,sph);
#endif

  }


  // Compute initial N-body forces
  //---------------------------------------------------------------------------
  if (nbody->Nstar > 0) {

    nbody->CalculateDirectGravForces(nbody->Nnbody,nbody->nbodydata);
    if (sph->self_gravity == 1 && sph->Nsph > 0)
      sphneib->UpdateAllStarGasForces(sph,nbody);
    nbody->CalculateAllStartupQuantities(nbody->Nnbody,nbody->nbodydata);

  }


  // Set particle values for initial step (e.g. r0, v0, a0)
  if (simparams->stringparams["gas_eos"] == "energy_eqn")
    uint->EndTimestep(n,sph->Nsph,sph->sphintdata);
  sphint->EndTimestep(n,sph->Nsph,sph->sphintdata);
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
  int activecount;                  // Flag if we need to recompute particles
  int i;                            // Particle loop counter
  int it;                           // Time-symmetric iteration counter
  int k;                            // Dimension counter
  FLOAT tghost;                     // Approx. ghost particle lifetime

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
  if (n == nresync) Nblocksteps = Nblocksteps + 1;

  // Advance SPH and N-body particles' positions and velocities
  sphint->AdvanceParticles(n,sph->Nsph,sph->sphintdata,(FLOAT) timestep);
  if (simparams->stringparams["gas_eos"] == "energy_eqn")
    uint->EnergyIntegration(n,sph->Nsph,sph->sphintdata,(FLOAT) timestep);
  nbody->AdvanceParticles(n,nbody->Nnbody,nbody->nbodydata,timestep);

  // Check all boundary conditions
  // (DAVID : Move this function to sphint and create an analagous one 
  //  for N-body.  Also, only check this on tree-build steps)
  if (Nsteps%ntreebuildstep == 0 || rebuild_tree)
    LocalGhosts->CheckBoundaries(simbox,sph);


  //---------------------------------------------------------------------------
  // MPI : On tree re-build step, determine load balancing for all MPI nodes.
  //       (How is this done?  All computed on root node??)
  //       Send/receive particles to their new nodes.
  //       Compute and transmit all bounding boxes (e.g. all particles, active
  //       particles, h-extent, ghosts, etc..) to all other MPI nodes
  //---------------------------------------------------------------------------
#ifdef MPI_PARALLEL
  if (Nsteps%ntreebuildstep == 0 || rebuild_tree) {
    mpicontrol.UpdateAllBoundingBoxes(sph->Nsph, sph->sphdata, sph->kernp);
    mpicontrol.LoadBalancing(sph,nbody);
    //exit(0);
  }
#endif


  // Compute all SPH quantities
  //---------------------------------------------------------------------------
  if (sph->Nsph > 0) {
    
    // Search for new ghost particles and create on local processor
    if (Nsteps%ntreebuildstep == 0 || rebuild_tree) {
      tghost = timestep*(FLOAT)(ntreebuildstep - 1);
      LocalGhosts->SearchGhostParticles(tghost,simbox,sph);
#ifdef MPI_PARALLEL
      MpiGhosts->SearchGhostParticles(tghost,simbox,sph);
#endif
    }
    // Otherwise copy properties from original particles to ghost particles
    else {
      LocalGhosts->CopySphDataToGhosts(simbox,sph);
#ifdef MPI_PARALLEL
      MpiGhosts->CopySphDataToGhosts(simbox,sph);
#endif
    }

    // Rebuild or update local neighbour and gravity tree
    sphneib->BuildTree(rebuild_tree,Nsteps,ntreebuildstep,
		       ntreestockstep,timestep,sph);
    rebuild_tree = false;

    // Reorder particles to tree-walk order (not implemented yet)

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
      LocalGhosts->CopySphDataToGhosts(simbox,sph);
#ifdef MPI_PARALLEL
      MpiGhosts->CopySphDataToGhosts(simbox,sph);
#endif
      
      // Zero accelerations
      //for (i=0; i<sph->Ntot; i++) {
      //  if (sph->sphdata[i].active) {
      //    for (k=0; k<ndim; k++) sph->sphdata[i].a[k] = (FLOAT) 0.0;
      //    for (k=0; k<ndim; k++) sph->sphdata[i].agrav[k] = (FLOAT) 0.0;
      //    sph->sphdata[i].gpot = (FLOAT) 0.0;
      //    sph->sphdata[i].gpe = (FLOAT) 0.0;
      //    sph->sphdata[i].dudt = (FLOAT) 0.0;
      //    sph->sphdata[i].levelneib = 0;
      //  }
      //}
      
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
      //for (i=0; i<sph->Nsph; i++) {
      //  if (sph->sphdata[i].active) {
      //    for (k=0; k<ndim; k++)
      //      sph->sphdata[i].a[k] += sph->sphdata[i].agrav[k];
      //  }
      //}

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

      if (sph->self_gravity == 1 && sph->Nsph > 0)
	sphneib->UpdateAllStarGasForces(sph,nbody);

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

      // Now compute minimum timestep due to stars/systems
#pragma omp parallel for
      for (i=0; i<nbody->Nnbody; i++)
	dt_min = min(dt_min,nbody->Timestep(nbody->nbodydata[i]));


#pragma omp critical
      if (dt < dt_min) dt_min = dt;
#pragma omp barrier

    }
    //-------------------------------------------------------------------------


    // For MPI, determine the global minimum timestep over all processors
#ifdef MPI_PARALLEL
    dt = dt_min;
    MPI_Allreduce(&dt,&dt_min,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
#endif


    // Set all particles to same timestep
    timestep = dt_min;
#pragma omp parallel for default(none)
    for (i=0; i<sph->Nsph; i++) {
      sph->sphdata[i].level = 0;
      sph->sphdata[i].levelneib = 0;
      sph->sphdata[i].dt = timestep;
      sph->sphintdata[i].nstep = pow(2,level_step - sph->sphdata[i].level);
      sph->sphintdata[i].nlast = n;

    }
#pragma omp parallel for default(none)
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
  int imin;                             // ..
  int istep;                            // ??
  int level;                            // Particle timestep level
  int last_level;                       // Previous timestep level
  int level_max_aux;                    // ..
  int level_max_old;                    // Old level_max
  int level_max_sph = 0;                // level_max for SPH particles only
  int level_min_sph = 9999999;          // level_min for SPH particles
  int level_max_nbody = 0;              // level_max for star particles only
  int level_nbody;                      // ..
  int level_sph;                        // ..
  int nfactor;                          // ??
  int nstep;                            // ??
  DOUBLE dt;                            // Aux. timestep variable
  DOUBLE dt_min = big_number_dp;        // ..
  DOUBLE dt_min_aux;                    // ..
  DOUBLE dt_min_nbody = big_number_dp;  // Maximum N-body particle timestep
  DOUBLE dt_min_sph = big_number_dp;    // Minimum SPH particle timestep
  DOUBLE dt_nbody;                      // ..
  DOUBLE dt_sph;                        // Aux. dt_min_sph_aux

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

#pragma omp parallel default(none) shared(dt_min_sph,dt_min_nbody) \
  private(dt,dt_min_aux,dt_nbody,dt_sph,i,imin)
    {
      // Initialise all timestep and min/max variables
      dt_min_aux = big_number_dp;
      dt_sph = big_number_dp;
      dt_nbody = big_number_dp;

      // Find minimum timestep from all SPH particles
#pragma omp for
      for (i=0; i<sph->Nsph; i++) {
	dt = min(sph->sphdata[i].dt,
		 sphint->Timestep(sph->sphdata[i],sph->hydro_forces));
	if (dt < dt_sph) imin = i;
	dt_min_aux = min(dt_min_aux,dt);
	dt_sph = min(dt_sph,dt);
	sph->sphdata[i].dt = dt;
      }
    
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
	dt_min_sph = min(dt_min_sph,dt_sph);
	dt_min_nbody = min(dt_min_nbody,dt_nbody);
      }
#pragma barrier
    }


    // For MPI, determine the global minimum timestep over all processors
#ifdef MPI_PARALLEL
    dt = timestep;
    MPI_Allreduce(&dt,&timestep,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
    dt = dt_min_sph;
    MPI_Allreduce(&dt,&dt_min_sph,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
#endif
    // Calculate new block timestep levels
    level_max = Nlevels - 1;
    level_step = level_max + integration_step - 1;
    dt_max = timestep*powf(2.0,level_max);
    
    // Calculate the maximum level occupied by all SPH particles
    level_max_sph = 
      min((int) (invlogetwo*log(dt_max/dt_min_sph)) + 1, level_max);
    level_max_nbody = 
      min((int) (invlogetwo*log(dt_max/dt_min_nbody)) + 1, level_max);
      
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
    

#pragma omp parallel default(none) shared(dt_min,dt_min_sph,dt_min_nbody) \
  shared(level_max_nbody,level_max_sph,level_min_sph)\
  private(dt,dt_min_aux,dt_nbody,dt_sph,i,imin,istep,last_level,level)\
  private(level_max_aux,level_nbody,level_sph,nstep,nfactor)
    {
      dt_min_aux = big_number_dp;
      dt_sph = big_number_dp;
      dt_nbody = big_number_dp;
      level_max_aux = 0;
      level_nbody = 0;
      level_sph = 0;

      // Find all SPH particles at the beginning of a new timestep
      //-----------------------------------------------------------------------
#pragma omp for
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
	level_sph = max(level_sph,sph->sphdata[i].level);
	if (sph->sphdata[i].dt < dt_sph) imin = i;
	//level_min_sph = min(level_min_sph,sph->sphdata[i].level);
	level_max_aux = max(level_max_aux,sph->sphdata[i].level);
	
	dt_sph = min(dt_sph,sph->sphdata[i].dt);
      }
      //-----------------------------------------------------------------------
      

#pragma omp critical
      {
        dt_min = min(dt_min,dt_min_aux);
        dt_min_sph = min(dt_min_sph,dt_sph);
        level_max = max(level_max,level_max_aux);
	level_max_sph = max(level_max_sph,level_sph);
      }
#pragma omp barrier

      // Now find all N-body particles at the beginning of a new timestep
      //-----------------------------------------------------------------------
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
	level_nbody = max(level_nbody,nbody->nbodydata[i]->level);
	level_max_aux = max(level_max_aux,nbody->nbodydata[i]->level);
	dt_nbody = min(dt_nbody,nbody->nbodydata[i]->dt);
      }
      //-----------------------------------------------------------------------
      

#pragma omp critical
      {
        dt_min = min(dt_min,dt_min_aux);
        dt_min_nbody = min(dt_min_nbody,dt_nbody);
        level_max = max(level_max,level_max_aux);
        level_max_nbody = max(level_max_nbody,level_nbody);
      }
#pragma omp barrier

    }


    // For MPI, find the global maximum timestep levels for each processor
#ifdef MPI_PARALLEL
    level = level_max;
    MPI_Allreduce(&level,&level_max,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
    level = level_max_sph;
    MPI_Allreduce(&level,&level_max_sph,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
#endif
    // For now, don't allow levels to be removed
    //level_max = max(level_max,level_max_old);
    level_step = level_max + integration_step - 1;
  
    // Set fixed SPH timestep level here in case maximum has changed
    if (sph_single_timestep == 1) {
      for (i=0; i<sph->Nsph; i++) {
	if (sph->sphintdata[i].nlast == n)
	  sph->sphdata[i].level = level_max_sph;
      }
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
	for (i=0; i<sph->Nsph; i++) sph->sphintdata[i].nstep *= nfactor;
	for (i=0; i<sph->Nsph; i++) sph->sphintdata[i].nlast *= nfactor;
	for (i=0; i<nbody->Nnbody; i++) nbody->nbodydata[i]->nstep *= nfactor;
	for (i=0; i<nbody->Nnbody; i++) nbody->nbodydata[i]->nlast *= nfactor;
      }
      else if (level_max < level_max_old) {
	nfactor = pow(2,level_max_old - level_max);
	n /= nfactor;
	for (i=0; i<sph->Nsph; i++) sph->sphintdata[i].nlast /= nfactor;
	for (i=0; i<sph->Nsph; i++) sph->sphintdata[i].nstep /= nfactor;
	for (i=0; i<nbody->Nnbody; i++) nbody->nbodydata[i]->nlast /= nfactor;
	for (i=0; i<nbody->Nnbody; i++) nbody->nbodydata[i]->nstep /= nfactor;
      }

      // Update values of nstep for both SPH and star particles
      for (i=0; i<sph->Nsph; i++) {
	if (sph->sphintdata[i].nlast == n)
	  sph->sphintdata[i].nstep = pow(2,level_step - sph->sphdata[i].level);
      }
      for (i=0; i<nbody->Nnbody; i++) {
	if (nbody->nbodydata[i]->nlast == n) nbody->nbodydata[i]->nstep = 
	  pow(2,level_step - nbody->nbodydata[i]->level);
      }
    
      nresync = pow(2,level_step);
      timestep = dt_max / (DOUBLE) nresync;

    }
    //-------------------------------------------------------------------------

  }
  //===========================================================================


#if defined(VERIFY_ALL)
  //VerifyBlockTimesteps();
#endif

  return;

  // Some validations
  //---------------------------------------------------------------------------
  int *ninlevel;
  int Nactive=0;
  ninlevel = new int[level_max+1];

  cout << "-----------------------------------------------------" << endl;
  cout << "Checking timesteps : " << level_max << "   " << level_max_sph << "    " << level_max_nbody << "    " << level_step << "   " << level_max_old << endl;
  cout << "n : " << n << endl;
  cout << "dt_min_sph : " << dt_min_sph << "    dt_min_nbody : " << dt_min_nbody << "    timestep : " << timestep << endl;
  cout << "imin : " << imin << "    " << sph->sphdata[imin].dt << "     " 
       << sph->sphdata[imin].h << "     " 
       << sqrt(DotProduct(sph->sphdata[imin].a,sph->sphdata[imin].a,ndim)) 
       << endl;
  for (int l=0; l<=level_max; l++) ninlevel[l] = 0;
  for (i=0; i<sph->Nsph; i++) if (sph->sphdata[i].active) Nactive++;
  for (i=0; i<sph->Nsph; i++) ninlevel[sph->sphdata[i].level]++;
  cout << "No. of active SPH particles : " << Nactive << endl;
  cout << "SPH level occupancy" << endl;
  for (int l=0; l<=level_max; l++) 
    cout << "level : " << l << "     N : " << ninlevel[l] << endl;

  for (int l=0; l<=level_max; l++) ninlevel[l] = 0;
  for (i=0; i<nbody->Nstar; i++) ninlevel[nbody->nbodydata[i]->level]++;
  cout << "N-body level occupancy" << endl;
  for (int l=0; l<=level_max; l++)
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
