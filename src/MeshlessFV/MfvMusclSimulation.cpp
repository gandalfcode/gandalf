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
  int activecount = 0;                 // Flag if we need to recompute particles
  int i;                               // Particle loop counter
  int it;                              // Time-symmetric iteration counter
  int k;                               // Dimension counter
  FLOAT tghost;                        // Approx. ghost particle lifetime
  MeshlessFVParticle<ndim> *partdata = mfv->GetMeshlessFVParticleArray();

  debug2("[MfvMusclSimulation::MainLoop]");


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

  // Advance time variables
  n++;
  Nsteps++;
  t = t + timestep;
  if (n == nresync) Nblocksteps++;
  if (n%integration_step == 0) Nfullsteps++;


  // Update the numerical fluxes of all active particles
  mfvneib->UpdateGodunovFluxes(mfv->Nhydro, mfv->Ntot, timestep, partdata, mfv, nbody);


  // Integrate all conserved variables to end of timestep
  //-----------------------------------------------------------------------------------------------
  for (i=0; i<mfv->Nhydro; i++) {
    MeshlessFVParticle<ndim> &part = partdata[i];
    int dn = n - part.nlast;

    for (k=0; k<ndim; k++) {
      part.r[k] = part.r0[k] + part.v0[k]*timestep*(FLOAT) dn;

      if (part.r[k] < simbox.boxmin[k]) {
        if (simbox.boundary_lhs[k] == periodicBoundary) {
          part.r[k] += simbox.boxsize[k];
          part.r0[k] += simbox.boxsize[k];
        }
      }
      if (part.r[k] > simbox.boxmax[k]) {
        if (simbox.boundary_rhs[k] == periodicBoundary) {
          part.r[k] -= simbox.boxsize[k];
          part.r0[k] -= simbox.boxsize[k];
        }
      }
    }

  }
  //-----------------------------------------------------------------------------------------------


  // Rebuild or update local neighbour and gravity tree
  mfvneib->BuildTree(rebuild_tree, Nsteps, ntreebuildstep, ntreestockstep,
                     mfv->Ntot, mfv->Nhydromax, timestep, partdata, mfv);


  // Search for new ghost particles and create on local processor
  //if (Nsteps%ntreebuildstep == 0 || rebuild_tree) {
  tghost = timestep*(FLOAT)(ntreebuildstep - 1);
  mfvneib->SearchBoundaryGhostParticles(tghost, simbox, mfv);
  mfv->CopyDataToGhosts(simbox, partdata);
  mfvneib->BuildGhostTree(rebuild_tree, Nsteps, ntreebuildstep, ntreestockstep,
                          mfv->Ntot, mfv->Nhydromax, timestep, partdata, mfv);


  // Calculate terms due to self-gravity
  if (mfv->self_gravity == 1) {
    mfvneib->UpdateAllGravForces(mfv->Nhydro, mfv->Ntot, partdata, mfv, nbody);
  }


  // End-step terms for all SPH particles
  mfv->EndTimestep(n, mfv->Nhydro, t, timestep, mfv->GetMeshlessFVParticleArray());

  /*this->CalculateDiagnostics();
  this->OutputDiagnostics();
  this->UpdateDiagnostics();*/


  rebuild_tree = true;


  return;
}



// Create template class instances of the main MfvMusclSimulation object for
// each dimension used (1, 2 and 3)
template class MfvMusclSimulation<1>;
template class MfvMusclSimulation<2>;
template class MfvMusclSimulation<3>;
