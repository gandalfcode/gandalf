//=================================================================================================
//  MfvMusclSimulation.cpp
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
//  MfvMusclSimulation::MainLoop
/// Main SPH simulation integration loop.
//=================================================================================================
template <int ndim>
void MfvMusclSimulation<ndim>::MainLoop(void)
{
  //int activecount = 0;                 // Flag if we need to recompute particles
  int i;                               // Particle loop counter
  //int it;                              // Time-symmetric iteration counter
  int k;                               // Dimension counter
  FLOAT tghost;                        // Approx. ghost particle lifetime
  MeshlessFVParticle<ndim> *partdata = mfv->GetMeshlessFVParticleArray();

  debug2("[MfvMusclSimulation::MainLoop]");

  // Update all active cell counters in the tree
  mfvneib->UpdateActiveParticleCounters(partdata, mfv);

  // Calculate all properties (and copy updated data to ghost particles)
  mfvneib->UpdateAllProperties(mfv->Nhydro, mfv->Ntot, partdata, mfv, nbody);
  mfv->CopyDataToGhosts(simbox, partdata);

  // Calculate all matrices and gradients (and copy updated data to ghost particles)
  mfvneib->UpdateGradientMatrices(mfv->Nhydro, mfv->Ntot, partdata, mfv, nbody);
  mfv->CopyDataToGhosts(simbox, partdata);


  // Compute timesteps for all particles
  if (Nlevels == 1) {
    this->ComputeGlobalTimestep();
  }
  else {
    this->ComputeBlockTimesteps();
  }
  mfv->CopyDataToGhosts(simbox, partdata);


  // Advance all global time variables
  n++;
  Nsteps++;
  t = t + timestep;
  if (n == nresync) Nblocksteps++;
  if (n%integration_step == 0) Nfullsteps++;


  // Update the numerical fluxes of all active particles
  mfvneib->UpdateGodunovFluxes(mfv->Nhydro, mfv->Ntot, timestep, partdata, mfv, nbody);


  if (Nsteps%1 == 0) {
  stringstream ss;
  string nostring = "";
  ss << setfill('0') << setw(5) << Nsteps;
  nostring = ss.str();
  string filename = run_id + ".SLOPES." + nostring;
  ss.str(std::string());
  ofstream outfile;
  outfile.open(filename.c_str());
  for (i=0; i<mfv->Nhydro; i++) {
    for (k=0; k<ndim; k++) outfile << partdata[i].r[k] << "    ";
    for (k=0; k<ndim+2; k++) outfile << partdata[i].Wprim[k] << "    ";
    for (k=0; k<ndim+2; k++) {
      for (int kk=0; kk<ndim; kk++) outfile << partdata[i].grad[k][kk] << "    ";
    }
    for (k=0; k<ndim+2; k++) outfile << partdata[i].dQdt[k] << "    ";
    outfile << endl;
  }
  outfile.close();
  cout << "WROTE FILE : " << filename << endl;
}


  // Integrate all conserved variables to end of timestep
  //-----------------------------------------------------------------------------------------------
  if (!mfv->staticParticles) {
    for (i=0; i<mfv->Nhydro; i++) {
      MeshlessFVParticle<ndim> &part = partdata[i];
      int dn = n - part.nlast;

      //-------------------------------------------------------------------------------------------
      for (k=0; k<ndim; k++) {
        part.r[k] = part.r0[k] + part.v0[k]*timestep*(FLOAT) dn;

        // Check if particle has crossed LHS boundary
        //-----------------------------------------------------------------------------------------
        if (part.r[k] < simbox.boxmin[k]) {

          // Check if periodic boundary
          if (simbox.boundary_lhs[k] == periodicBoundary) {
            part.r[k]  += simbox.boxsize[k];
            part.r0[k] += simbox.boxsize[k];
          }

          // Check if wall or mirror boundary
          if (simbox.boundary_lhs[k] == mirrorBoundary || simbox.boundary_lhs[k] == wallBoundary) {
            part.r[k]  = (FLOAT) 2.0*simbox.boxmin[k] - part.r[k];
            part.r0[k] = (FLOAT) 2.0*simbox.boxmin[k] - part.r0[k];
            part.v[k]  = -part.v[k];
            part.v0[k] = -part.v0[k];
            part.a[k]  = -part.a[k];
            part.a0[k] = -part.a0[k];
          }
        }

        // Check if particle has crossed RHS boundary
        //-----------------------------------------------------------------------------------------
        if (part.r[k] > simbox.boxmax[k]) {

          // Check if periodic boundary
          if (simbox.boundary_rhs[k] == periodicBoundary) {
            part.r[k]  -= simbox.boxsize[k];
            part.r0[k] -= simbox.boxsize[k];
          }

          // Check if wall or mirror boundary
          if (simbox.boundary_rhs[k] == mirrorBoundary || simbox.boundary_rhs[k] == wallBoundary) {
            part.r[k]  = (FLOAT) 2.0*simbox.boxmax[k] - part.r[k];
            part.r0[k] = (FLOAT) 2.0*simbox.boxmax[k] - part.r0[k];
            part.v[k]  = -part.v[k];
            part.v0[k] = -part.v0[k];
            part.a[k]  = -part.a[k];
            part.a0[k] = -part.a0[k];
          }

        }
        //-----------------------------------------------------------------------------------------

      }
      //-------------------------------------------------------------------------------------------

    }
  }
  //-----------------------------------------------------------------------------------------------

  // Advance N-body particle positions
  nbody->AdvanceParticles(n, nbody->Nnbody, t, timestep, nbody->nbodydata);


  // Rebuild or update local neighbour and gravity tree
  mfvneib->BuildTree(rebuild_tree, Nsteps, ntreebuildstep, ntreestockstep,
                     mfv->Ntot, mfv->Nhydromax, timestep, partdata, mfv);


  // Search for new ghost particles and create on local processor
  //if (Nsteps%ntreebuildstep == 0 || rebuild_tree) {
  tghost = timestep*(FLOAT) (ntreebuildstep - 1);
  mfvneib->SearchBoundaryGhostParticles(tghost, simbox, mfv);
  mfv->CopyDataToGhosts(simbox, partdata);
  mfvneib->BuildGhostTree(rebuild_tree, Nsteps, ntreebuildstep, ntreestockstep,
                          mfv->Ntot, mfv->Nhydromax, timestep, partdata, mfv);


  // Calculate terms due to self-gravity
  if (mfv->self_gravity == 1) {
    mfvneib->UpdateAllGravForces(mfv->Nhydro, mfv->Ntot, partdata, mfv, nbody);
  }

  // Compute N-body forces
  //-----------------------------------------------------------------------------------------------
  if (nbody->Nnbody > 0) {

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

    if (mfv->self_gravity == 1 && mfv->Nhydro > 0) {
      mfvneib->UpdateAllStarGasForces(mfv->Nhydro, mfv->Ntot, partdata, mfv, nbody);
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
  mfv->EndTimestep(n, mfv->Nhydro, t, timestep, mfv->GetMeshlessFVParticleArray());

  // End-step terms for all star particles
  if (nbody->Nstar > 0) nbody->EndTimestep(n,nbody->Nnbody,t,timestep,nbody->nbodydata);


  // Search for new sink particles (if activated) and accrete to existing sinks
  if (sink_particles == 1) {
    if (sinks.create_sinks == 1 && (rebuild_tree || Nfullsteps%ntreebuildstep == 0)) {
      sinks.SearchForNewSinkParticles(n, t, mfv, nbody);
    }
    if (sinks.Nsink > 0) {
      sinks.AccreteMassToSinks(n, timestep, mfv, nbody);
      nbody->UpdateStellarProperties();
      //if (extra_sink_output) WriteExtraSinkOutput();
    }
    // If we will output a snapshot (regular or for restarts), then delete all accreted particles
    if ((t >= tsnapnext && sinks.Nsink > 0) || n == nresync || kill_simulation ||
         timing->WallClockTime() - timing->tstart_wall > 0.99*tmax_wallclock) {
      hydro->DeleteDeadParticles();
      rebuild_tree = true;
    }
  }

  rebuild_tree = true;

  return;
}



// Create template class instances of the main MfvMusclSimulation object for
// each dimension used (1, 2 and 3)
template class MfvMusclSimulation<1>;
template class MfvMusclSimulation<2>;
template class MfvMusclSimulation<3>;
