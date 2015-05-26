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
  MeshlessFVParticle<ndim> *partdata;         // Pointer to main SPH data array

  debug2("[MfvMusclSimulation::MainLoop]");


  // Set pointer for SPH data array
  partdata = mfv->GetMeshlessFVParticleArray();


  // Advance SPH and N-body particles' positions and velocities
  /*uint->EnergyIntegration(n,mfv->Nhydro,(FLOAT) t,(FLOAT) timestep,mfv->GetMeshlessFVParticleArray());
  sphint->AdvanceParticles(n,mfv->Nhydro,(FLOAT) t,(FLOAT) timestep,mfv->GetMeshlessFVParticleArray());
  nbody->AdvanceParticles(n,nbody->Nnbody,t,timestep,nbody->nbodydata);*/

  // Check all boundary conditions
  // (DAVID : Move this function to sphint and create an analagous one
  //  for N-body.  Also, only check this on tree-build steps)
  //if (Nsteps%ntreebuildstep == 0 || rebuild_tree) sphint->CheckBoundaries(simbox,sph);


    //}
    // Otherwise copy properties from original particles to ghost particles
    /*else {
      LocalGhosts->CopySphDataToGhosts(simbox, sph);
#ifdef MPI_PARALLEL
      MpiGhosts->CopySphDataToGhosts(simbox, sph);
#endif
    }*/

  for (i=0; i<mfv->Nhydro; i++) {
    MeshlessFVParticle<ndim>& part = mfv->GetMeshlessFVParticlePointer(i);
    part.active = true;
    for (k=0; k<ndim+2; k++) partdata[i].Qcons0[k] = partdata[i].Qcons[k];
    for (k=0; k<ndim; k++) partdata[i].rdmdt0[k] = partdata[i].rdmdt[k];
  }

  // Calculate all properties
  mfvneib->UpdateAllProperties(mfv->Nhydro, mfv->Ntot, partdata, mfv, nbody);

  mfv->CopyDataToGhosts(simbox, partdata);

  mfvneib->UpdateGradientMatrices(mfv->Nhydro, mfv->Ntot, partdata, mfv, nbody);

  mfv->CopyDataToGhosts(simbox, partdata);


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


  mfvneib->UpdateGodunovFluxes(mfv->Nhydro, mfv->Ntot, timestep, partdata, mfv, nbody);


  // Integrate all conserved variables to end of timestep
  //-----------------------------------------------------------------------------------------------
  for (i=0; i<mfv->Nhydro; i++) {

    for (k=0; k<ndim; k++) {
      //partdata[i].r[k] += (FLOAT) 0.5*(partdata[i].v0[k] + partdata[i].v[k])*timestep;
      partdata[i].a0[k] = partdata[i].a[k];
      partdata[i].v0[k] = partdata[i].v[k];
      partdata[i].r[k] += partdata[i].v[k]*timestep;
      if (partdata[i].r[k] < simbox.boxmin[k]) {
        if (simbox.boundary_lhs[k] == periodicBoundary) {
          partdata[i].r[k] += simbox.boxsize[k];
        }
      }
      if (partdata[i].r[k] > simbox.boxmax[k]) {
        if (simbox.boundary_rhs[k] == periodicBoundary) {
          partdata[i].r[k] -= simbox.boxsize[k];
        }
      }
    }

  }
  //-----------------------------------------------------------------------------------------------



  // Integrate all conserved variables to end of timestep
  for (i=0; i<mfv->Nhydro; i++) {
    mfv->IntegrateConservedVariables(partdata[i], timestep);
    mfv->ConvertQToConserved(partdata[i].volume, partdata[i].Qcons, partdata[i].Ucons);
    mfv->ConvertConservedToPrimitive(partdata[i].Ucons, partdata[i].Wprim);
    mfv->UpdateArrayVariables(partdata[i]);
  }


  // Rebuild or update local neighbour and gravity tree
  mfvneib->BuildTree(rebuild_tree, Nsteps, ntreebuildstep, ntreestockstep,
                     mfv->Ntot, mfv->Nhydromax, timestep, partdata, mfv);


  // Search for new ghost particles and create on local processor
  //if (Nsteps%ntreebuildstep == 0 || rebuild_tree) {
  tghost = timestep*(FLOAT)(ntreebuildstep - 1);
  mfvneib->SearchBoundaryGhostParticles(tghost, simbox, mfv);
  mfvneib->BuildGhostTree(rebuild_tree,Nsteps,ntreebuildstep,ntreestockstep,
                          mfv->Ntot,mfv->Nhydromax,timestep,partdata,mfv);


  //-----------------------------------------------------------------------------------------------
  if (mfv->self_gravity == 1) {
    mfvneib->UpdateAllGravForces(mfv->Nhydro, mfv->Ntot, partdata, mfv, nbody);

    // Integrate all conserved variables to end of timestep
    for (i=0; i<mfv->Nhydro; i++) {

      for (k=0; k<ndim; k++) {
        partdata[i].Qcons[k] += (FLOAT) 0.5*timestep*
          (partdata[i].Qcons0[MeshlessFV<ndim>::irho]*partdata[i].a0[k] +
           partdata[i].Qcons[MeshlessFV<ndim>::irho]*partdata[i].a[k]);
      }
      partdata[i].Qcons[MeshlessFV<ndim>::ietot] += (FLOAT) 0.5*timestep*
        (partdata[i].Qcons0[MeshlessFV<ndim>::irho]*DotProduct(partdata[i].v0, partdata[i].a0, ndim) +
         partdata[i].Qcons[MeshlessFV<ndim>::irho]*DotProduct(partdata[i].v, partdata[i].a, ndim)); //+
         //DotProduct(partdata[i].a0, partdata[i].rdmdt0, ndim) +
         //DotProduct(partdata[i].a, partdata[i].rdmdt, ndim));

      mfv->ConvertQToConserved(partdata[i].volume, partdata[i].Qcons, partdata[i].Ucons);
      mfv->ConvertConservedToPrimitive(partdata[i].Ucons, partdata[i].Wprim);
      mfv->UpdateArrayVariables(partdata[i]);
    }

  }
  //-----------------------------------------------------------------------------------------------


  /*this->CalculateDiagnostics();
  this->OutputDiagnostics();
  this->UpdateDiagnostics();*/


  rebuild_tree = false;


  // End-step terms for all SPH particles
  /*if (mfv->Nhydro > 0) {
    uint->EndTimestep(n,mfv->Nhydro,(FLOAT) t,(FLOAT) timestep,mfv->GetMeshlessFVParticleArray());
    sphint->EndTimestep(n,mfv->Nhydro,(FLOAT) t,(FLOAT) timestep,mfv->GetMeshlessFVParticleArray());
  }*/

  return;
}



// Create template class instances of the main MfvMusclSimulation object for
// each dimension used (1, 2 and 3)
template class MfvMusclSimulation<1>;
template class MfvMusclSimulation<2>;
template class MfvMusclSimulation<3>;
