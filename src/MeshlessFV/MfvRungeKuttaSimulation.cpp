//=================================================================================================
//  MfvRungeKuttaSimulation.cpp
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
//  MfvRungeKuttaSimulation::MainLoop
/// Main SPH simulation integration loop.
//=================================================================================================
template <int ndim>
void MfvRungeKuttaSimulation<ndim>::MainLoop(void)
{
  //int activecount = 0;                 // Flag if we need to recompute particles
  int i;                               // Particle loop counter
  //int it;                              // Time-symmetric iteration counter
  int k;                               // Dimension counter
  FLOAT tghost;                        // Approx. ghost particle lifetime

  debug2("[MfvRungeKuttaSimulation::MainLoop]");


  // Compute timesteps for all particles
  if (Nlevels == 1) {
    this->ComputeGlobalTimestep();
  }
  else {
    this->ComputeBlockTimesteps();
  }
  LocalGhosts->CopyHydroDataToGhosts(simbox,mfv);

  // Update the numerical fluxes of all active particles
  if (mfv->hydro_forces) {
    mfvneib->UpdateGodunovFluxes(timestep, mfv, nbody, simbox);
  }

  // Advance all global time variables
  n++;
  Nsteps++;
  t = t + timestep;
  if (n == nresync) Nblocksteps++;
  if (n%integration_step == 0) Nfullsteps++;


  // Integrate positions of particles
  hydroint->AdvanceParticles(n, t, timestep, mfv);
  nbody->AdvanceParticles(n, nbody->Nnbody, t, timestep, nbody->nbodydata);

  // Check all boundary conditions
  // (DAVID : create an analagous of this function for N-body)
  hydroint->CheckBoundaries(simbox,mfv);

  // Re-build/re-stock tree now particles have moved
  mfvneib->BuildTree(rebuild_tree, Nsteps, ntreebuildstep, ntreestockstep, timestep, mfv);
  mfvneib->BuildGhostTree(rebuild_tree, Nsteps, ntreebuildstep, ntreestockstep, timestep, mfv);


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
         timing->RunningTime() > (FLOAT) 0.99*tmax_wallclock) {
      hydro->DeleteDeadParticles();
      rebuild_tree = true;
    }

    // Re-build/re-stock tree now particles have moved
    mfvneib->BuildTree(rebuild_tree, Nsteps, ntreebuildstep, ntreestockstep, timestep, mfv);
    mfvneib->BuildGhostTree(rebuild_tree, Nsteps, ntreebuildstep, ntreestockstep, timestep, mfv);

  }

  // Calculate terms due to self-gravity
  if (mfv->self_gravity == 1) {
    mfvneib->UpdateAllGravForces(mfv, nbody, simbox, ewald);
  }

  // Compute N-body forces
  //-----------------------------------------------------------------------------------------------
  if (nbody->Nnbody > 0) {

    // Zero all acceleration terms
    for (i=0; i<nbody->Nnbody; i++) {
      if (nbody->nbodydata[i]->flags.check(active)) {
        for (k=0; k<ndim; k++) nbody->nbodydata[i]->a[k]     = (FLOAT) 0.0;
        for (k=0; k<ndim; k++) nbody->nbodydata[i]->adot[k]  = (FLOAT) 0.0;
        for (k=0; k<ndim; k++) nbody->nbodydata[i]->a2dot[k] = (FLOAT) 0.0;
        for (k=0; k<ndim; k++) nbody->nbodydata[i]->a3dot[k] = (FLOAT) 0.0;
        nbody->nbodydata[i]->gpot = (FLOAT) 0.0;
        nbody->nbodydata[i]->gpe = (FLOAT) 0.0;
      }
    }

    if (mfv->self_gravity == 1 && mfv->Nhydro > 0) {
      mfvneib->UpdateAllStarGasForces(mfv, nbody, simbox, ewald);
#if defined MPI_PARALLEL
      // We need to sum up the contributions from the different domains
      mpicontrol->ComputeTotalStarGasForces(nbody);
#endif
    }

    // Calculate forces, force derivatives etc.., for active stars/systems
    if (nbody->nbody_softening == 1) {
      nbody->CalculateDirectSmoothedGravForces(nbody->Nnbody, nbody->nbodydata, simbox, ewald);
    }
    else {
      nbody->CalculateDirectGravForces(nbody->Nnbody, nbody->nbodydata, simbox, ewald);
    }

    for (i=0; i<nbody->Nnbody; i++) {
      if (nbody->nbodydata[i]->flags.check(active)) {
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
  hydroint->EndTimestep(n, t, timestep, mfv);
  nbody->EndTimestep(n, nbody->Nnbody, t, timestep, nbody->nbodydata);


  // Rebuild or update local neighbour and gravity tree
  mfvneib->BuildTree(rebuild_tree, Nsteps, ntreebuildstep, ntreestockstep, timestep, mfv);

  // Search for new ghost particles and create on local processor
  if (Nsteps%ntreebuildstep == 0 || rebuild_tree) {
    tghost = timestep*(FLOAT) (ntreebuildstep - 1);
    mfvneib->SearchBoundaryGhostParticles(tghost, simbox, mfv);
    LocalGhosts->CopyHydroDataToGhosts(simbox,mfv);
    mfvneib->BuildGhostTree(rebuild_tree, Nsteps, ntreebuildstep, ntreestockstep, timestep, mfv);
  }

  // Update all active cell counters in the tree
  mfvneib->UpdateActiveParticleCounters(mfv);

  // Calculate all properties (and copy updated data to ghost particles)
  mfvneib->UpdateAllProperties(mfv, nbody);
  LocalGhosts->CopyHydroDataToGhosts(simbox,mfv);


  // Calculate all matrices and gradients (and copy updated data to ghost particles)
  mfvneib->UpdateGradientMatrices(mfv, nbody, simbox);
  LocalGhosts->CopyHydroDataToGhosts(simbox,mfv);




  /*//int activecount = 0;                 // Flag if we need to recompute particles
  int i;                               // Particle loop counter
  //int it;                              // Time-symmetric iteration counter
  int k;                               // Dimension counter
  FLOAT tghost;                        // Approx. ghost particle lifetime
  MeshlessFVParticle<ndim> *partdata = mfv->GetMeshlessFVParticleArray();

  debug2("[MfvRungeKuttaSimulation::MainLoop]");


  // Compute timesteps for all particles
  if (Nlevels == 1) {
    this->ComputeGlobalTimestep();
  }
  else {
    this->ComputeBlockTimesteps();
  }
  LocalGhosts->CopyHydroDataToGhosts(simbox,mfv);


  // Update the numerical fluxes of all active particles
  if (mfv->hydro_forces) {
    mfvneib->UpdateGodunovFluxes(mfv->Nhydro, mfv->Ntot, timestep, partdata, mfv, nbody, simbox);
  }

  // Advance all global time variables
  n++;
  Nsteps++;
  t = t + timestep;
  if (n == nresync) Nblocksteps++;
  if (n%integration_step == 0) Nfullsteps++;


  // Integrate positions of particles
  mfv->IntegrateParticles(n, mfv->Nhydro, t, timestep, simbox, partdata);

  // Advance N-body particle positions
  nbody->AdvanceParticles(n, nbody->Nnbody, t, timestep, nbody->nbodydata);

  // Re-build/re-stock tree now particles have moved
  mfvneib->BuildTree(rebuild_tree, Nsteps, ntreebuildstep, ntreestockstep,
                     mfv->Ntot, mfv->Nhydromax, timestep, partdata, mfv);
  mfvneib->BuildGhostTree(rebuild_tree, Nsteps, ntreebuildstep, ntreestockstep,
                          mfv->Ntot, mfv->Nhydromax, timestep, partdata, mfv);


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
      sinks->AccreteMassToSinks(n, timestep, partdata, mfv, nbody);
      nbody->UpdateStellarProperties();
      //if (extra_sink_output) WriteExtraSinkOutput();
    }
    // If we will output a snapshot (regular or for restarts), then delete all accreted particles
    if ((t >= tsnapnext && sinks->Nsink > 0) || n == nresync || kill_simulation ||
         timing->WallClockTime() - timing->tstart_wall > (FLOAT) 0.99*tmax_wallclock) {
      hydro->DeleteDeadParticles();
      rebuild_tree = true;
    }

    // Update all array variables now accretion has probably removed some mass
    for (i=0; i<mfv->Nhydro; i++) partdata[i].m = partdata[i].Qcons[FV<ndim>::irho] + partdata[i].dQ[FV<ndim>::irho];

    // Re-build/re-stock tree now particles have moved
    mfvneib->BuildTree(rebuild_tree, Nsteps, ntreebuildstep, ntreestockstep,
                       mfv->Ntot, mfv->Nhydromax, timestep, partdata, mfv);
    mfvneib->BuildGhostTree(rebuild_tree, Nsteps, ntreebuildstep, ntreestockstep,
                            mfv->Ntot, mfv->Nhydromax, timestep, partdata, mfv);

  }

  // Calculate terms due to self-gravity
  if (mfv->self_gravity == 1) {
    mfvneib->UpdateAllGravForces(mfv->Nhydro, mfv->Ntot, partdata, mfv, nbody, simbox, ewald);
  }

  // Compute N-body forces
  //-----------------------------------------------------------------------------------------------
  if (nbody->Nnbody > 0) {

    // Zero all acceleration terms
    for (i=0; i<nbody->Nnbody; i++) {
      if (nbody->nbodydata[i]->flags.check(active)) {
        for (k=0; k<ndim; k++) nbody->nbodydata[i]->a[k]     = (FLOAT) 0.0;
        for (k=0; k<ndim; k++) nbody->nbodydata[i]->adot[k]  = (FLOAT) 0.0;
        for (k=0; k<ndim; k++) nbody->nbodydata[i]->a2dot[k] = (FLOAT) 0.0;
        for (k=0; k<ndim; k++) nbody->nbodydata[i]->a3dot[k] = (FLOAT) 0.0;
        nbody->nbodydata[i]->gpot = (FLOAT) 0.0;
        nbody->nbodydata[i]->gpe = (FLOAT) 0.0;
      }
    }
    if (sink_particles == 1) {
      for (i=0; i<sinks->Nsink; i++) {
        if (sinks->sink[i].star->flags.check(active)) {
          for (k=0; k<ndim; k++) sinks->sink[i].fhydro[k] = (FLOAT) 0.0;
        }
      }
    }

    if (mfv->self_gravity == 1 && mfv->Nhydro > 0) {
      mfvneib->UpdateAllStarGasForces(mfv->Nhydro, mfv->Ntot, partdata, mfv, nbody);
#if defined MPI_PARALLEL
      // We need to sum up the contributions from the different domains
      mpicontrol->ComputeTotalStarGasForces(nbody);
#endif
    }

    // Calculate forces, force derivatives etc.., for active stars/systems
    if (nbody->nbody_softening == 1) {
      nbody->CalculateDirectSmoothedGravForces(nbody->Nnbody, nbody->nbodydata, simbox, ewald);
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




  // Advance SPH and N-body particles' positions and velocities
  //uint->EnergyIntegration(n,mfv->Nhydro,(FLOAT) t,(FLOAT) timestep,mfv->GetMeshlessFVParticleArray());
  //sphint->AdvanceParticles(n,mfv->Nhydro,(FLOAT) t,(FLOAT) timestep,mfv->GetMeshlessFVParticleArray());
  //nbody->AdvanceParticles(n,nbody->Nnbody,t,timestep,nbody->nbodydata);

  // Check all boundary conditions
  // (DAVID : Move this function to sphint and create an analagous one
  //  for N-body.  Also, only check this on tree-build steps)
  //if (Nsteps%ntreebuildstep == 0 || rebuild_tree) sphint->CheckBoundaries(simbox,sph);


  // Rebuild or update local neighbour and gravity tree
  mfvneib->BuildTree(rebuild_tree, Nsteps, ntreebuildstep, ntreestockstep,
                     mfv->Ntot, mfv->Nhydromax, timestep, partdata, mfv);


  // Search for new ghost particles and create on local processor
  //if (Nsteps%ntreebuildstep == 0 || rebuild_tree) {
  tghost = timestep*(FLOAT)(ntreebuildstep - 1);
  mfvneib->SearchBoundaryGhostParticles(tghost, simbox, mfv);
  // Update pointer in case there has been a reallocation
  partdata = sph->GetSphParticleArray();
  mfvneib->BuildGhostTree(rebuild_tree,Nsteps,ntreebuildstep,ntreestockstep,
                          mfv->Ntot,mfv->Nhydromax,timestep,partdata,mfv);
    //}
    // Otherwise copy properties from original particles to ghost particles
    //else {
    //  LocalGhosts->CopyHydroDataToGhosts(simbox, sph);
#ifdef MPI_PARALLEL
    //  MpiGhosts->CopyHydroDataToGhosts(simbox, sph);
#endif
    //}



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


  for (i=0; i<mfv->Nhydro; i++) {
    MeshlessFVParticle<ndim>& part = mfv->GetMeshlessFVParticlePointer(i);
    part.flags.set(active);
  }
  LocalGhosts->CopyHydroDataToGhosts(simbox,mfv);



  // Drift particles to the half-step position
  for (i=0; i<mfv->Nhydro; i++) {

    for (k=0; k<ndim+2; k++) partdata[i].Qcons0[k] = partdata[i].Qcons[k];

    for (k=0; k<ndim; k++) {
      partdata[i].r[k] += (FLOAT) 0.5*partdata[i].v[k]*timestep;
      //if (partdata[i].r[k] < simbox.min[k])
      //  if (simbox.boundary_lhs[k] == periodicBoundary) {
      //    partdata[i].r[k] += simbox.size[k];
      //  }
      //if (partdata[i].r[k] > simbox.max[k])
      //  if (simbox.boundary_rhs[k] == periodicBoundary) {
      //    partdata[i].r[k] -= simbox.size[k];
      //  }
    }

  }


  // Calculate all properties at the half-step
  LocalGhosts->CopyHydroDataToGhosts(simbox,mfv);
  mfvneib->UpdateAllProperties(mfv->Nhydro, mfv->Ntot, partdata, mfv, nbody, simbox);
  LocalGhosts->CopyHydroDataToGhosts(simbox,mfv);

  // Calculate all properties at the half-step
  mfvneib->UpdateGradientMatrices(mfv->Nhydro, mfv->Ntot, partdata, mfv, nbody, simbox);
  LocalGhosts->CopyHydroDataToGhosts(simbox,mfv);
  mfvneib->UpdateGodunovFluxes(mfv->Nhydro, mfv->Ntot, timestep, partdata, mfv, nbody, simbox);


  // Integrate all conserved variables to end of timestep using predicted fluxes
  for (i=0; i<mfv->Nhydro; i++) {

    for (k=0; k<ndim+2; k++) partdata[i].Qcons[k] += partdata[i].dQdt[k]*timestep;

    //mfv->IntegrateConservedVariables(partdata[i], timestep);
    mfv->ConvertConservedToPrimitive(partdata[i].volume, partdata[i].Qcons, partdata[i].Wprim);
    mfv->UpdateArrayVariables(partdata[i]);

  }



  // Calculate all properties at the half-step
  LocalGhosts->CopyHydroDataToGhosts(simbox,mfv);
  mfvneib->UpdateGradientMatrices(mfv->Nhydro, mfv->Ntot, partdata, mfv, nbody, simbox);
  LocalGhosts->CopyHydroDataToGhosts(simbox,mfv);
  mfvneib->UpdateGodunovFluxes(mfv->Nhydro, mfv->Ntot, timestep, partdata, mfv, nbody, simbox);


  // Correction step for all conserved variables to end of timestep
  for (i=0; i<mfv->Nhydro; i++) {

    for (k=0; k<ndim+2; k++) partdata[i].Qcons[k] =
      (FLOAT) 0.5*(partdata[i].Qcons[k] + partdata[i].Qcons0[k] + partdata[i].dQdt[k]*timestep);

    mfv->ConvertConservedToPrimitive(partdata[i].volume, partdata[i].Qcons, partdata[i].Wprim);
    mfv->UpdateArrayVariables(partdata[i]);

    for (k=0; k<ndim; k++) {
      partdata[i].r[k] += (FLOAT) 0.5*partdata[i].v[k]*timestep;
      if (partdata[i].r[k] < simbox.min[k])
        if (simbox.boundary_lhs[k] == periodicBoundary) {
          partdata[i].r[k] += simbox.size[k];
        }
      if (partdata[i].r[k] > simbox.max[k])
        if (simbox.boundary_rhs[k] == periodicBoundary) {
          partdata[i].r[k] -= simbox.size[k];
        }
    }


  }

  rebuild_tree = false;

*/

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



// Create template class instances of the main MfvRungeKuttaSimulation object for
// each dimension used (1, 2 and 3)
template class MfvRungeKuttaSimulation<1>;
template class MfvRungeKuttaSimulation<2>;
template class MfvRungeKuttaSimulation<3>;
