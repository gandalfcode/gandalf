//=================================================================================================
//  MfvLeapFrogSimulation.cpp
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

//=================================================================================================
//  MfvLeapFrogSimulation::PostInitialConditionsSetup
/// Call routines for calculating all initial hydro and N-body quantities
/// once initial conditions have been set-up.
//=================================================================================================
template <int ndim>
void MfvLeapFrogSimulation<ndim>::PostInitialConditionsSetup(void)
{
  int i;                               // Particle counter
  int k;                               // Dimension counter

  debug2("[MfvLeapFrogSimulation::PostInitialConditionsSetup]");

  // Set iorig
  if (rank == 0) {
    MeshlessFVParticle<ndim> *partdata = mfv->GetMeshlessFVParticleArray();
    for (i=0; i<mfv->Nhydro; i++) {
      partdata[i].iorig = i;
      partdata[i].flags.set_flag(active);
      partdata[i].flags.set_flag(update_density) ;
      partdata[i].gpot = big_number;
    }
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
  mpicontrol->CreateInitialDomainDecomposition(mfv, nbody, simparams, simbox,
                                               this->initial_h_provided);
  this->AllocateParticleMemory();
#endif

  // Set time variables here (for now)
  nresync = 0;   // DAVID : Need to adapt this for block timesteps
  n = 0;
  integration_step = 1;
  level_step = 1;

  // Set all relevant particle counters
  mfv->Nghost = 0;
  mfv->Ntot = mfv->Nhydro;

  // Initialise conserved variables
  for (i=0; i<mfv->Nhydro; i++) {
    MeshlessFVParticle<ndim>& part = mfv->GetMeshlessFVParticlePointer(i);
    part.Qcons0[MeshlessFV<ndim>::irho] = part.m;
    for (int k=0; k<ndim; k++) part.Qcons0[k] = part.m*part.v[k];
    FLOAT ekin = (FLOAT) 0.0;
    for (int k=0; k<ndim; k++) ekin += part.v[k]*part.v[k];
    part.Qcons0[MeshlessFV<ndim>::ietot] = part.u*part.m + (FLOAT) 0.5*part.m*ekin;
  }



  // If the smoothing lengths have not been provided beforehand, then
  // calculate the initial values here

  mfvneib->ToggleNeighbourCheck(false);
  if (!this->initial_h_provided) {
    mfv->InitialSmoothingLengthGuess();
    mfvneib->BuildTree(rebuild_tree, 0, ntreebuildstep, ntreestockstep, timestep, mfv);
    mfvneib->UpdateAllProperties(mfv, nbody, simbox);

      for (i=0; i<mfv->Nhydro; i++) {
        MeshlessFVParticle<ndim>& part = mfv->GetMeshlessFVParticlePointer(i);
        part.flags.set_flag(update_density) ;
      }
  }
  else {
    mfvneib->BuildTree(rebuild_tree, 0, ntreebuildstep, ntreestockstep, timestep, mfv);
  }

#ifdef MPI_PARALLEL
  mpicontrol->UpdateAllBoundingBoxes(mfv->Nhydro, mfv, mfv->kernp);
#endif

  // Search ghost particles
  mfvneib->SearchBoundaryGhostParticles((FLOAT) 0.0, simbox, mfv);
  mfvneib->BuildGhostTree(true, 0, ntreebuildstep, ntreestockstep, timestep, mfv);
#ifdef MPI_PARALLEL
  mpicontrol->UpdateAllBoundingBoxes(mfv->Nhydro+mfv->NPeriodicGhost, mfv, mfv->kernp);
  for (int i=0; i<mfv->Nhydro+mfv->NPeriodicGhost; i++) {
    MeshlessFVParticle<ndim>& parti =  mfv->GetMeshlessFVParticlePointer(i);
    parti.hrangesqd = mfv->kernp->kernrangesqd*parti.h*parti.h;
  }
  MpiGhosts->SearchGhostParticles((FLOAT) 0.0, simbox, mfv);
  mfvneib->BuildMpiGhostTree(true, 0, ntreebuildstep, ntreestockstep, timestep, mfv);
#endif

  // Zero accelerations
  for (i=0; i<mfv->Nhydro; i++) {
    MeshlessFVParticle<ndim>& part = mfv->GetMeshlessFVParticlePointer(i);
    part.flags.set_flag(active);
    part.flags.set_flag(update_density) ;
    part.gpot = big_number;
  }

  // Update neighbour tree
  rebuild_tree = true;
  mfvneib->BuildTree(rebuild_tree, 0, ntreebuildstep, ntreestockstep, timestep, mfv);
#ifdef MPI_PARALLEL
  mfvneib->InitialiseCellWorkCounters();
#endif


  // Update neighbour tree
  rebuild_tree = true;
  mfvneib->BuildTree(rebuild_tree, 0, ntreebuildstep, ntreestockstep, timestep, mfv);

  // Calculate all hydro properties
  mfvneib->UpdateAllProperties(mfv, nbody, simbox);

#ifdef MPI_PARALLEL
  mpicontrol->UpdateAllBoundingBoxes(mfv->Nhydro, mfv, mfv->kernp);
#endif

  // Search ghost particles
  mfvneib->SearchBoundaryGhostParticles((FLOAT) 0.0, simbox, mfv);
  mfvneib->BuildGhostTree(true, 0, ntreebuildstep, ntreestockstep, timestep, mfv);
#ifdef MPI_PARALLEL
  mpicontrol->UpdateAllBoundingBoxes(mfv->Nhydro + mfv->NPeriodicGhost, mfv, mfv->kernp);
  MpiGhosts->SearchGhostParticles((FLOAT) 0.0, simbox, mfv);
  mfvneib->BuildMpiGhostTree(true, 0, ntreebuildstep, ntreestockstep, timestep, mfv);
#endif


  // Update neighbour tree
  rebuild_tree = true;
  mfvneib->BuildTree(rebuild_tree, 0, ntreebuildstep, ntreestockstep, timestep, mfv);
#ifdef MPI_PARALLEL
  mfvneib->InitialiseCellWorkCounters();
#endif
  mfvneib->SearchBoundaryGhostParticles((FLOAT) 0.0, simbox, mfv);
  mfvneib->BuildGhostTree(true, 0, ntreebuildstep, ntreestockstep, timestep, mfv);
  // Communicate pruned trees for MPI
#ifdef MPI_PARALLEL
  mfvneib->BuildPrunedTree(rank, simbox, mpicontrol->mpinode, mfv);
#endif
  mfvneib->ToggleNeighbourCheck(true);

  // Compute mean mass (used for smooth sink accretion)
  if (!restart) {
    MeshlessFVParticle<ndim> *partdata = mfv->GetMeshlessFVParticleArray();
    mfv->mmean = (FLOAT) 0.0;
    for (i=0; i<mfv->Nhydro; i++) mfv->mmean += partdata[i].m;
    mfv->mmean /= (FLOAT) mfv->Nhydro;
  }

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

  // Read-in N-body table here
  nbody->LoadStellarPropertiesTable(&simunits);
  nbody->UpdateStellarProperties();


  // Compute all initial hydro terms
  //-----------------------------------------------------------------------------------------------
  for (i=0; i<mfv->Ntot; i++) {
    MeshlessFVParticle<ndim>& part = mfv->GetMeshlessFVParticlePointer(i);
    for (k=0; k<ndim; k++) {
        part.a[k] = (FLOAT) 0.0;
        part.rdmdt[k] = 0.0;
        part.rdmdt0[k] = 0.0;
    }
    part.level  = 0;
    part.nstep  = 0;
    part.nlast  = 0;
    part.tlast  = t;
    part.flags.unset_flag(active);
  }
  for (i=0; i<mfv->Nhydro; i++) mfv->GetMeshlessFVParticlePointer(i).flags.set_flag(active);

  // Copy all other data from real hydro particles to ghosts
  LocalGhosts->CopyHydroDataToGhosts(simbox,mfv);

  mfvneib->BuildTree(true, 0, ntreebuildstep, ntreestockstep, timestep, mfv);
#ifdef MPI_PARALLEL
  mfvneib->InitialiseCellWorkCounters();
#endif
  mfvneib->SearchBoundaryGhostParticles((FLOAT) 0.0, simbox, mfv);
  mfvneib->BuildGhostTree(true, 0, ntreebuildstep, ntreestockstep, timestep, mfv);
#ifdef MPI_PARALLEL
  mpicontrol->UpdateAllBoundingBoxes(mfv->Nhydro + mfv->NPeriodicGhost, mfv, mfv->kernp);
  MpiGhosts->SearchGhostParticles((FLOAT) 0.0, simbox, mfv);
  mfvneib->BuildMpiGhostTree(true, 0, ntreebuildstep, ntreestockstep, timestep, mfv);
#endif

  // ..
  for (i=0; i<mfv->Nhydro; i++) {
    MeshlessFVParticle<ndim>& part = mfv->GetMeshlessFVParticlePointer(i);
    for (k=0; k<ndim; k++) part.r0[k] = part.r[k];
    for (k=0; k<ndim; k++) part.v0[k] = part.v[k];
    for (k=0; k<ndim; k++) part.a0[k] = part.a[k];
    part.flags.set_flag(active);
  }

  mfvneib->UpdateAllProperties(mfv, nbody, simbox);


  LocalGhosts->CopyHydroDataToGhosts(simbox,mfv);
#ifdef MPI_PARALLEL
  MpiGhosts->CopyHydroDataToGhosts(simbox,mfv);
#endif

  // Update the primitive vectors for all particles
  for (i=0; i<mfv->Ntot; i++) {
    MeshlessFVParticle<ndim>& part = mfv->GetMeshlessFVParticlePointer(i);
    mfv->ComputeThermalProperties(part);
    mfv->UpdatePrimitiveVector(part);
    mfv->ConvertPrimitiveToConserved(part.ndens, part.Wprim, part.Qcons0);
    for (k=0; k<ndim+2; k++) part.dQ[k] = (FLOAT) 0.0;
  }

#ifdef MPI_PARALLEL
  //mfvneib->UpdateHydroExportList(rank, mfv->Nhydro, mfv->Ntot, partdata, mfv, nbody, simbox);

  //mpicontrol->ExportParticlesBeforeForceLoop(mfv);
  // Update pointer in case there has been a reallocation
  //partdata = mfv->GetMeshlessFVParticleArray();
#endif


  mfvneib->UpdateGradientMatrices(mfv, nbody, simbox);
  LocalGhosts->CopyHydroDataToGhosts(simbox,mfv);
#ifdef MPI_PARALLEL
  MpiGhosts->CopyHydroDataToGhosts(simbox,mfv);
#endif


  if (mfv->self_gravity == 1 || nbody->Nnbody > 0) {
#ifdef MPI_PARALLEL
    if (mfv->self_gravity ==1 ) {
      MeshlessFVParticle<ndim> *partdata = mfv->GetMeshlessFVParticleArray();
      for (int i=0; i< mfv->Nhydro; i++) {
        for (int k=0; k<ndim; k++) partdata[i].a[k] = 0.0;
        partdata[i].gpot=0.0;
      }
      mfvneib->UpdateGravityExportList(rank, mfv, nbody, simbox);
      mpicontrol->ExportParticlesBeforeForceLoop(mfv);
    }
#endif
    mfvneib->UpdateAllGravForces(mfv, nbody, simbox, ewald);
#ifdef MPI_PARALLEL
    if (mfv->self_gravity ==1 ) {
    mpicontrol->GetExportedParticlesAccelerations(mfv);
    }
#endif
  }


  // Compute initial N-body forces
  //-----------------------------------------------------------------------------------------------
  if (mfv->self_gravity == 1 && mfv->Nhydro > 0) {
    mfvneib->UpdateAllStarGasForces(mfv, nbody);  //, simbox, ewald);

    // We need to sum up the contributions from the different domains
#if defined MPI_PARALLEL
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

  this->CalculateDiagnostics();
  this->diag0 = this->diag;
  this->setup = true;


  // Compute timesteps for all particles
  if (Nlevels == 1) {
    this->ComputeGlobalTimestep();
  }
  else {
    this->ComputeBlockTimesteps();
  }


#ifdef MPI_PARALLEL
  // Pruned trees are used only to compute which particles to export
  // Therefore we don't need to update them at the start of the loop, and we can do it soon before we need them
  if (Nsteps%ntreebuildstep == 0 || rebuild_tree) {
      mfvneib->BuildPrunedTree(rank, simbox, mpicontrol->mpinode, mfv);
  }
  else {
      mfvneib->StockPrunedTree(rank, mfv);
  }

  if (mfv->hydro_forces) {
    mfvneib->UpdateHydroExportList(rank, mfv, nbody, simbox);

    mpicontrol->ExportParticlesBeforeForceLoop(mfv);
  }
#endif
  // Update the numerical fluxes of all active particles
  if (mfv->hydro_forces) {
    mfvneib->UpdateGodunovFluxes(timestep/2, mfv, nbody, simbox);
  }

#ifdef MPI_PARALLEL
  if (mfv->hydro_forces) mpicontrol->GetExportedParticlesAccelerations(mfv);
#endif

  for (i=0; i<mfv->Ntot; i++) {
    MeshlessFVParticle<ndim>& part = mfv->GetMeshlessFVParticlePointer(i);
    FLOAT Qcons[ndim+2];
    for (int var=0; var<ndim+2; var++) {
      Qcons[var] = part.Qcons0[var] + part.dQ[var];
      part.dQ[var]    = (FLOAT) 0.0;
      part.dQdt[var]  = (FLOAT) 0.0;
    }
    mfv->UpdateArrayVariables(part, Qcons);
    mfv->ComputeThermalProperties(part);
    mfv->UpdatePrimitiveVector(part) ;
  }

  // Call EndTimestep to set all 'beginning-of-step' variables
  //mfv->EndTimestep(0, t, timestep);
  //nbody->EndTimestep(n, nbody->Nstar, t, timestep, nbody->nbodydata);

  return;
}



//=================================================================================================
//  MfvLeapFrogSimulation::MainLoop
/// Main SPH simulation integration loop.
//=================================================================================================
template <int ndim>
void MfvLeapFrogSimulation<ndim>::MainLoop(void)
{
  //int activecount = 0;                 // Flag if we need to recompute particles
  int i;                               // Particle loop counter
  //int it;                              // Time-symmetric iteration counter
  int k;                               // Dimension counter
  FLOAT tghost;                        // Approx. ghost particle lifetime

  debug2("[MfvLeapFrogSimulation:MainLoop]");

  // Advance all global time variables
  n++;
  Nsteps++;
  t = t + timestep;
  if (n == nresync) Nblocksteps++;
  if (n%integration_step == 0) Nfullsteps++;


  // Integrate positions of particles
  mfv->IntegrateParticles(n, t, timestep, simbox);

  // Advance N-body particle positions
  nbody->AdvanceParticles(n, nbody->Nnbody, t, timestep, nbody->nbodydata);

#ifdef MPI_PARALLEL
  if (Nsteps%ntreebuildstep == 0 || rebuild_tree) {
	// Horrible hack in order NOT to trigger a full tree rebuild
	mfvneib->BuildTree(rebuild_tree,Nsteps+1,2, ntreestockstep,timestep,mfv);
	if (rebuild_tree) {
		  mfvneib->BuildPrunedTree(rank, simbox, mpicontrol->mpinode, mfv);
	}
	else {
		mfvneib->StockPrunedTree(rank, mfv);
	}
	mpicontrol->UpdateAllBoundingBoxes(mfv->Nhydro, mfv, mfv->kernp);
    mpicontrol->LoadBalancing(mfv, nbody);
  }
#endif

  // Re-build/re-stock tree now particles have moved
  mfvneib->BuildTree(rebuild_tree, Nsteps, ntreebuildstep, ntreestockstep, timestep, mfv);
#ifdef MPI_PARALLEL
  mfvneib->InitialiseCellWorkCounters();
#endif

	//tghost = timestep*(FLOAT) (ntreebuildstep - 1);
	tghost = 0;
	mfvneib->SearchBoundaryGhostParticles(tghost, simbox, mfv);
	mfvneib->BuildGhostTree(true, Nsteps, ntreebuildstep, ntreestockstep,timestep, mfv);
	#ifdef MPI_PARALLEL
	mpicontrol->UpdateAllBoundingBoxes(mfv->Nhydro + mfv->NPeriodicGhost, mfv, mfv->kernp);
	MpiGhosts->SearchGhostParticles(tghost, simbox, mfv);
	  mfvneib->BuildMpiGhostTree(true, Nsteps, ntreebuildstep, ntreestockstep,  timestep, mfv);
	#endif


  // Search for new sink particles (if activated) and accrete to existing sinks
  if (sink_particles == 1) {
    if (sinks->create_sinks == 1 && (rebuild_tree || Nfullsteps%ntreebuildstep == 0)) {
      sinks->SearchForNewSinkParticles(n, t, mfv, nbody);
    }
    if (sinks->Nsink > 0) {
      mfv->mmean = (FLOAT) 0.0;
      for (i=0; i<mfv->Nhydro; i++) mfv->mmean += mfv->GetMeshlessFVParticlePointer(i).m;
      mfv->mmean /= (FLOAT) mfv->Nhydro;
      mfv->hmin_sink = big_number;
      for (i=0; i<sinks->Nsink; i++) {
        mfv->hmin_sink = min(mfv->hmin_sink, (FLOAT) sinks->sink[i].star->h);
      }
      sinks->AccreteMassToSinks(n, timestep, mfv, nbody);
      nbody->UpdateStellarProperties();
      //if (extra_sink_output) WriteExtraSinkOutput();
    }
    // If we will output a snapshot (regular or for restarts), then delete all accreted particles
    if ((t >= tsnapnext && sinks->Nsink > 0) || n == nresync || kill_simulation ||
         timing->RunningTime()  > (FLOAT) 0.99*tmax_wallclock) {
      hydro->DeleteDeadParticles();
      rebuild_tree = true;
    }

    // Re-build/re-stock tree now particles have moved
    mfvneib->BuildTree(rebuild_tree, Nsteps, ntreebuildstep, ntreestockstep,timestep, mfv);
    LocalGhosts->CopyHydroDataToGhosts(simbox,mfv);
    mfvneib->BuildGhostTree(rebuild_tree, Nsteps, ntreebuildstep, ntreestockstep, timestep, mfv);
#ifdef MPI_PARALLEL
	MpiGhosts->CopyHydroDataToGhosts(simbox,mfv);
    mfvneib->BuildMpiGhostTree(rebuild_tree, Nsteps, ntreebuildstep, ntreestockstep, timestep, mfv);
#endif

  }
  // Update all active cell counters in the tree
  mfvneib->UpdateActiveParticleCounters(mfv);

  //Calculate all properties (and copy updated data to ghost particles)
  mfvneib->UpdateAllProperties(mfv, nbody, simbox);

#ifdef MPI_PARALLEL
  LocalGhosts->CopyHydroDataToGhosts(simbox,mfv);
  MpiGhosts->CopyHydroDataToGhosts(simbox,mfv);
#endif


  // Calculate all matrices and gradients (and copy updated data to ghost particles)
  // TODO:
  //   Compute gradients for all cells neighbouring active ones (use levelneib?).
  mfvneib->UpdateGradientMatrices(mfv, nbody, simbox);


  // Compute timesteps for all particles
  double oldtimestep=timestep;
  if (Nlevels == 1) {
    this->ComputeGlobalTimestep();
  }
  else {
    this->ComputeBlockTimesteps();
  }


#ifdef MPI_PARALLEL
  // Pruned trees are used only to compute which particles to export
  // Therefore we don't need to update them at the start of the loop, and we can do it soon before we need them
  if (Nsteps%ntreebuildstep == 0 || rebuild_tree) {
      mfvneib->BuildPrunedTree(rank, simbox, mpicontrol->mpinode, mfv);
  }
  else {
      mfvneib->StockPrunedTree(rank, mfv);
  }

  if (mfv->hydro_forces) {
    mfvneib->UpdateHydroExportList(rank, mfv, nbody, simbox);

    mpicontrol->ExportParticlesBeforeForceLoop(mfv);
  }
#endif
  // Update the numerical fluxes of all active particles
  if (mfv->hydro_forces) {
    mfvneib->UpdateGodunovFluxes(0.5*(timestep+oldtimestep), mfv, nbody, simbox);
  }
#ifdef MPI_PARALLEL
  if (mfv->hydro_forces) mpicontrol->GetExportedParticlesAccelerations(mfv);
#endif

  // Calculate terms due to self-gravity / stars
  if (mfv->self_gravity == 1 || nbody->Nnbody > 0) {
    // Update the density to get the correct softening & grad-h terms.
    mfvneib->UpdateAllProperties(mfv, nbody, simbox);
    LocalGhosts->CopyHydroDataToGhosts(simbox,mfv);
#ifdef MPI_PARALLEL
    if (mfv->self_gravity ==1 ) {
      if (Nsteps%ntreebuildstep == 0 || rebuild_tree) {
    	     mfvneib->BuildPrunedTree(rank, simbox, mpicontrol->mpinode, mfv);
    	}
      else {
    		 mfvneib->StockPrunedTree(rank, mfv);
      }
      mfvneib->UpdateGravityExportList(rank, mfv, nbody, simbox);
      mpicontrol->ExportParticlesBeforeForceLoop(mfv);
    }
#endif
    // Does only the star forces in mfv->self_gravity != 1
    mfvneib->UpdateAllGravForces(mfv, nbody, simbox, ewald);
#ifdef MPI_PARALLEL
    if (mfv->self_gravity ==1 ) {
    mpicontrol->GetExportedParticlesAccelerations(mfv);
    }
#endif
  }

  // Compute N-body forces
  //-----------------------------------------------------------------------------------------------
  if (nbody->Nnbody > 0) {

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

    if (mfv->self_gravity == 1 && mfv->Nhydro > 0) {
      mfvneib->UpdateAllStarGasForces(mfv, nbody);
#if defined MPI_PARALLEL
      // We need to sum up the contributions from the different domains
      mpicontrol->ComputeTotalStarGasForces(nbody);
#endif
    }

    // Calculate forces, force derivatives etc.., for active stars/systems
    if (nbody->nbody_softening == 1) {
      nbody->CalculateDirectSmoothedGravForces(nbody->Nnbody, nbody->nbodydata);
    }
    else {
      nbody->CalculateDirectGravForces(nbody->Nnbody, nbody->nbodydata);
    }

    for (i=0; i<nbody->Nnbody; i++) {
      if (nbody->nbodydata[i]->active) {
        nbody->extpot->AddExternalPotential(nbody->nbodydata[i]->r, nbody->nbodydata[i]->v,
                                            nbody->nbodydata[i]->a, nbody->nbodydata[i]->adot,
                                            nbody->nbodydata[i]->gpot);
      }
    }

    // Calculate correction step for all stars at end of step.
    nbody->CorrectionTerms(n, nbody->Nnbody, t, timestep, nbody->nbodydata);

  }
  //-----------------------------------------------------------------------------------------------


  // End-step terms for all hydro particles
  mfv->EndTimestep(n, t, timestep);
  nbody->EndTimestep(n, nbody->Nnbody, t, timestep, nbody->nbodydata);


  /* Check that we have sensible smoothing lengths */
  if (periodicBoundaries) {
    double hmax = mfvneib->GetMaximumSmoothingLength() ;
    hmax *= mfv->kernp->kernrange ;
    for (i=0; i < ndim; i++)
      if (simbox.half[i] < 2*hmax){
        string message = "Error: Smoothing length too large, self-interaction will occur" ;
    	ExceptionHandler::getIstance().raise(message);
      }
  }

  return;
}



// Create template class instances of the main MfvLeapFrogSimulation object for
// each dimension used (1, 2 and 3)
template class MfvLeapFrogSimulation<1>;
template class MfvLeapFrogSimulation<2>;
template class MfvLeapFrogSimulation<3>;
