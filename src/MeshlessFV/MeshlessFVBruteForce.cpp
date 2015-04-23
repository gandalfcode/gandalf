//=================================================================================================
//  MeshlessFVBruteForce.cpp
//  Contains all routines for generating SPH neighbour lists using
//  brute-force (i.e. direct summation over all particles).
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
#include <math.h>
#include "MfvNeighbourSearch.h"
#include "NeighbourSearch.h"
#include "Sph.h"
#include "Parameters.h"
#include "Particle.h"
#include "Debug.h"
#include "InlineFuncs.h"
#include "SmoothingKernel.h"
#if defined MPI_PARALLEL
#include "MpiNode.h"
#endif
using namespace std;



//=================================================================================================
//  MeshlessFVBruteForce::MeshlessFVBruteForce
/// MeshlessFVBruteForce class constructor
//=================================================================================================
template <int ndim, template<int> class ParticleType>
MeshlessFVBruteForce<ndim,ParticleType>::MeshlessFVBruteForce
(FLOAT kernrangeaux,
 DomainBox<ndim> *boxaux,
 SmoothingKernel<ndim> *kernaux,
 CodeTiming *timingaux):
 NeighbourSearch<ndim>(kernrangeaux, boxaux, kernaux, timingaux),
 MeshlessFVNeighbourSearch<ndim>(kernrangeaux, boxaux, kernaux, timingaux),
 BruteForceSearch<ndim,ParticleType>(kernrangeaux, boxaux, kernaux, timingaux)
{
}



//=================================================================================================
//  MeshlessFVBruteForce::~MeshlessFVBruteForce
/// MeshlessFVBruteForce class destructor
//=================================================================================================
template <int ndim, template<int> class ParticleType>
MeshlessFVBruteForce<ndim,ParticleType>::~MeshlessFVBruteForce()
{
}



//=================================================================================================
//  MeshlessFVBruteForce::UpdateAllSphProperties
/// Routine for computing SPH properties (smoothing lengths, densities and forces) for all active
/// SPH particle using neighbour lists generated using brute force (i.e. direct summation).
//=================================================================================================
template <int ndim, template<int> class ParticleType>
void MeshlessFVBruteForce<ndim,ParticleType>::UpdateAllProperties
 (int Nhydro,                          ///< [in] No. of SPH particles
  int Ntot,                            ///< [in] No. of SPH + ghost particles
  MeshlessFVParticle<ndim> *mfvdata,   ///< [inout] Pointer to SPH ptcl array
  MeshlessFV<ndim> *mfv,               ///< [in] Pointer to SPH object
  Nbody<ndim> *nbody)                  ///< [in] Pointer to N-body object
{
  int i,j,jj,k;                        // Particle and dimension counters
  int Nneib = 0;                       // No. of (non-dead) neighbours
  //int okflag;                          // Flag valid smoothing length
  int *neiblist;                       // List of neighbours
  FLOAT dr[ndim];                      // Relative distance vector
  FLOAT rp[ndim];                      // Position of current particle
  FLOAT *drsqd;                        // Distance squared
  FLOAT *gpot;                         // Array of neib. grav. potentials
  FLOAT *m;                            // Array of neib. position vectors
  FLOAT *mu = 0;                       // Array of neib. mass*u values
  //ParticleType<ndim>* mfvdata = static_cast<ParticleType<ndim>* > (mfv_gen);

  debug2("[MeshlessFVBruteForce::UpdateAllProperties]");
  //timing->StartTimingSection("MFV_COMPUTE_H");

  // Store masses in separate array
  gpot = new FLOAT[Ntot];
  m = new FLOAT[Ntot];
  neiblist = new int[Ntot];
  for (i=0; i<Ntot; i++) {
    if (mfvdata[i].itype == dead) continue;
    neiblist[Nneib] = i;
    gpot[Nneib] = mfvdata[i].gpot;
    m[Nneib] = mfvdata[i].m;
    Nneib++;
  }

  // Create parallel threads
  //===============================================================================================
#pragma omp parallel default(none) private(dr,drsqd,i,j,jj,k,rp) \
  shared(gpot,m,mu,nbody,neiblist,Nneib,Nhydro,Ntot,mfv,mfvdata)
  {
    drsqd = new FLOAT[Ntot];

    // Compute smoothing lengths of all SPH particles
    //---------------------------------------------------------------------------------------------
#pragma omp for
    for (i=0; i<Nhydro; i++) {

      // Skip over inactive particles
      if (!mfvdata[i].active || mfvdata[i].itype == dead) continue;

      for (k=0; k<ndim; k++) rp[k] = mfvdata[i].r[k];

      // Compute distances and the reciprical between the current particle and all neighbours here
      //-------------------------------------------------------------------------------------------
      for (jj=0; jj<Nneib; jj++) {
        j = neiblist[jj];
        for (k=0; k<ndim; k++) dr[k] = mfvdata[j].r[k] - rp[k];
        drsqd[jj] = DotProduct(dr,dr,ndim);
      }
      //-------------------------------------------------------------------------------------------

      // Compute all SPH gather properties
      //okflag =
      mfv->ComputeH(i, Nneib, big_number, m, mu, drsqd, gpot, mfvdata[i], nbody);

    }
    //---------------------------------------------------------------------------------------------

    delete[] drsqd;

  }
  //===============================================================================================

  delete[] neiblist;
  delete[] m;
  delete[] gpot;

  //timing->EndTimingSection("MFV_COMPUTE_H");

  return;
}



//=================================================================================================
//  MeshlessFVBruteForceSearch::UpdateGradientMatrices
/// Routine for computing SPH properties (smoothing lengths, densities and
/// forces) for all active SPH particle using neighbour lists generated
/// using brute force (i.e. direct summation).
//=================================================================================================
template <int ndim, template<int> class ParticleType>
void MeshlessFVBruteForce<ndim,ParticleType>::UpdateGradientMatrices
 (int Nhydro,                          ///< [in] No. of SPH particles
  int Ntot,                            ///< [in] No. of SPH + ghost particles
  MeshlessFVParticle<ndim> *mfvdata,   ///< [inout] Pointer to SPH ptcl array
  MeshlessFV<ndim> *mfv,               ///< [in] Pointer to SPH object
  Nbody<ndim> *nbody)                  ///< [in] Pointer to N-body object
{
  debug2("[MeshlessFVBruteForce::UpdateGradientMatrices]");
  timing->StartTimingSection("MFV_UPDATE_GRADIENTS");

  // Compute forces of real and imported particles
  //-----------------------------------------------------------------------------------------------
#pragma omp parallel default(none) shared(mfv,mfvdata,Nhydro,Ntot)
  {
    int i,j,k;                           // Particle and dimension counters
    int Nneib;                           // No. of neighbours
    int *neiblist;                       // List of neighbour ids
    FLOAT draux[ndim];                   // Relative distance vector
    FLOAT drsqd;                         // Distance squared
    FLOAT hrangesqdi;                    // Gather kernel extent (squared)
    FLOAT hrangesqdj;                    // Scatter kernel extent (squared)
    FLOAT rp[ndim];                      // Position of current particle
    FLOAT *dr;                           // Array of neib. position vectors
    FLOAT *drmag;                        // Array of neib. distances
    FLOAT *invdrmag;                     // Array of neib. inverse distances

    // Allocate memory for storing neighbour ids and position data
    neiblist = new int[Ntot];
    dr       = new FLOAT[ndim*Ntot];
    drmag    = new FLOAT[Ntot];
    invdrmag = new FLOAT[Ntot];

    const int offset_imported = 0; //mfv->Nghost;

    // Compute forces of real and imported particles
    //---------------------------------------------------------------------------------------------
#pragma omp for
    for (int ipart=0; ipart<Nhydro; ipart++) {

      if (ipart < Nhydro) i = ipart;
      else i = ipart + offset_imported;

      // Skip over inactive particles
      if (!mfvdata[i].active || mfvdata[i].itype == dead) continue;

      for (k=0; k<ndim; k++) rp[k] = mfvdata[i].r[k];
      hrangesqdi = mfvdata[i].hrangesqd; //pow(kernfac*kernp->kernrange*mfvdata[i].h,2);
      Nneib = 0;

      //cout << "kern : " << kernfac << "   " << kernp->kernrange << endl;
      //cin >> j;

      // Compute distances and the reciprical between the current particle
      // and all neighbours here
      //-------------------------------------------------------------------------------------------
      for (j=0; j<mfv->Nhydro + mfv->NPeriodicGhost; j++) {
        if (mfvdata[j].itype == dead) continue;
        hrangesqdj = mfvdata[j].hrangesqd; //pow(kernfac*kernp->kernrange*mfvdata[j].h,2);
        for (k=0; k<ndim; k++) draux[k] = mfvdata[j].r[k] - rp[k];
        drsqd = DotProduct(draux,draux,ndim);

        //cout << "drsqd : " << drsqd << "   " << hrangesqdi << "   " << hrangesqdj << endl;

        if ((drsqd < hrangesqdi || drsqd < hrangesqdj) && i != j) {
          neiblist[Nneib] = j;
          drmag[Nneib] = sqrt(drsqd);
          invdrmag[Nneib] = (FLOAT) 1.0/(drmag[Nneib] + small_number);
          for (k=0; k<ndim; k++) dr[Nneib*ndim + k] = draux[k]*invdrmag[Nneib];
          Nneib++;
        }
      }
      //-------------------------------------------------------------------------------------------

      // Compute all SPH hydro forces
      mfv->ComputePsiFactors(i, Nneib, neiblist, drmag, invdrmag, dr, mfvdata[i], mfvdata);
      mfv->ComputeGradients(i, Nneib, neiblist, drmag, invdrmag, dr, mfvdata[i], mfvdata);

      //sphdata[i].active = false;

    }
    //---------------------------------------------------------------------------------------------


    // Free all allocated memory
    delete[] invdrmag;
    delete[] drmag;
    delete[] dr;
    delete[] neiblist;

  }
  //-----------------------------------------------------------------------------------------------

  timing->EndTimingSection("MFV_UPDATE_GRADIENTS");

  return;
}



//=================================================================================================
//  MeshlessFVBruteForce::UpdateGodunovFluxes
/// Routine for computing SPH properties (smoothing lengths, densities and
/// forces) for all active SPH particle using neighbour lists generated
/// using brute force (i.e. direct summation).
//=================================================================================================
template <int ndim, template<int> class ParticleType>
void MeshlessFVBruteForce<ndim,ParticleType>::UpdateGodunovFluxes
 (const int Nhydro,                    ///< [in] No. of SPH particles
  const int Ntot,                      ///< [in] No. of SPH + ghost particles
  const FLOAT timestep,                ///< [in] ..
  MeshlessFVParticle<ndim> *mfvdata,   ///< [inout] Pointer to SPH ptcl array
  MeshlessFV<ndim> *mfv,               ///< [in] Pointer to SPH object
  Nbody<ndim> *nbody)                  ///< [in] Pointer to N-body object
{
  debug2("[MeshlessFVBruteForce::UpdateGodunovFluxes]");
  timing->StartTimingSection("MFV_UPDATE_FLUXES");


  // Compute forces of real and imported particles
  //-----------------------------------------------------------------------------------------------
#pragma omp parallel default(none) shared(mfv,mfvdata)
  {
    int i,j,k;                             // Particle and dimension counters
    int Nneib;                             // No. of neighbours
    int *neiblist;                         // List of neighbour ids
    FLOAT draux[ndim];                     // Relative distance vector
    FLOAT drsqd;                           // Distance squared
    FLOAT hrangesqdi;                      // Gather kernel extent (squared)
    FLOAT hrangesqdj;                      // Scatter kernel extent (squared)
    FLOAT rp[ndim];                        // Position of current particle
    FLOAT *dr;                             // Array of neib. position vectors
    FLOAT *drmag;                          // Array of neib. distances
    FLOAT *invdrmag;                       // Array of neib. inverse distances
    FLOAT (*fluxBuffer)[ndim+2];           // ..
    MeshlessFVParticle<ndim> *neibdata;    // ..
    const int offset_imported = 0; //mfv->Nghost;

    // Allocate memory for storing neighbour ids and position data
    neiblist   = new int[Ntot];
    dr         = new FLOAT[ndim*Ntot];
    drmag      = new FLOAT[Ntot];
    invdrmag   = new FLOAT[Ntot];
    fluxBuffer = new FLOAT[Ntot][ndim+2];
    neibdata   = new MeshlessFVParticle<ndim>[Ntot];

    // ..
    for (i=0; i<Ntot; i++) neibdata[i] = mfvdata[i];
    for (i=0; i<Ntot; i++) {
      for (k=0; k<ndim+2; k++) fluxBuffer[i][k] = (FLOAT) 0.0;
    }

    //---------------------------------------------------------------------------------------------
#pragma omp for
    for (int ipart=0; ipart<Nhydro; ipart++) {

      if (ipart < Nhydro) i = ipart;
      else i = ipart + offset_imported;

      // Skip over inactive particles
      if (!mfvdata[i].active || mfvdata[i].itype == dead) continue;

      for (k=0; k<ndim; k++) rp[k] = mfvdata[i].r[k];
      hrangesqdi = mfvdata[i].hrangesqd; //pow(kernfac*kernp->kernrange*mfvdata[i].h,2);
      Nneib = 0;
      for (j=0; j<Ntot; j++) {
        for (k=0; k<ndim+2; k++) neibdata[j].dQdt[k] = (FLOAT) 0.0;
      }

      // Compute distances and the reciprical between the current particle
      // and all neighbours here
      //-------------------------------------------------------------------------------------------
      for (j=0; j<mfv->Nhydro + mfv->NPeriodicGhost; j++) {
        if (mfvdata[j].itype == dead) continue;
        hrangesqdj = mfvdata[j].hrangesqd; //pow(kernfac*kernp->kernrange*mfvdata[j].h,2);
        for (k=0; k<ndim; k++) draux[k] = mfvdata[j].r[k] - rp[k];
        drsqd = DotProduct(draux,draux,ndim);

        //if ((drsqd < hrangesqdi || drsqd < hrangesqdj) && i != j) {
        if ((drsqd < hrangesqdi || drsqd < hrangesqdj) && i != j  &&
            (!neibdata[j].active || j > i)) {
          neiblist[Nneib] = j;
          drmag[Nneib]    = sqrt(drsqd) + small_number;
          invdrmag[Nneib] = (FLOAT) 1.0/drmag[Nneib];
          for (k=0; k<ndim; k++) dr[Nneib*ndim + k] = draux[k]*invdrmag[Nneib];
          Nneib++;
        }
      }
      //-------------------------------------------------------------------------------------------


      // Compute all flux terms
      mfv->ComputeGodunovFlux(i, Nneib, timestep, neiblist, drmag, invdrmag, dr,
                              mfvdata[i], neibdata);

      // Accumlate fluxes
      for (int jj=0; jj<Nneib; jj++) {
        j = neiblist[jj];
        for (k=0; k<ndim+2; k++) fluxBuffer[j][k] += neibdata[j].dQdt[k];
      }

    }
    //---------------------------------------------------------------------------------------------

    // Add all buffers back to main arrays
#pragma omp barrier
#pragma omp critical
    {
      for (i=0; i<Nhydro; i++) {
        for (k=0; k<ndim+2; k++) mfvdata[i].dQdt[k] += fluxBuffer[i][k];
      }
    }

    // Free all allocated memory
    delete[] neibdata;
    delete[] fluxBuffer;
    delete[] invdrmag;
    delete[] drmag;
    delete[] dr;
    delete[] neiblist;

  }
  //-----------------------------------------------------------------------------------------------

  timing->EndTimingSection("MFV_UPDATE_FLUXES");

  return;
}



//=================================================================================================
//  MeshlessFVBruteForce::SearchBoundaryGhostParticles
/// Search domain to create any required ghost particles near any boundaries.
/// Currently only searches to create periodic or mirror ghost particles.
//=================================================================================================
/*template <int ndim, template <int> class ParticleType>
void MeshlessFVBruteForce<ndim,ParticleType>::SearchBoundaryGhostParticles
 (FLOAT tghost,                        ///< Ghost particle 'lifetime'
  DomainBox<ndim> &simbox,             ///< Simulation box structure
  MeshlessFV<ndim> *mfv)               ///< Sph object pointer
{
  int i;                               // Particle counter

  // Set all relevant particle counters
  mfv->Nghost         = 0;
  mfv->NPeriodicGhost = 0;
  mfv->Nghostmax      = mfv->Nhydromax - mfv->Nhydro;
  mfv->Ntot           = mfv->Nhydro;


  // If all boundaries are open, immediately return to main loop
  if (simbox.boundary_lhs[0] == openBoundary && simbox.boundary_rhs[0] == openBoundary &&
      simbox.boundary_lhs[1] == openBoundary && simbox.boundary_rhs[1] == openBoundary &&
      simbox.boundary_lhs[2] == openBoundary && simbox.boundary_rhs[2] == openBoundary)
    return;


  debug2("[BruteForceSearch::SearchBoundaryGhostParticles]");


  // Create ghost particles in x-dimension
  //-----------------------------------------------------------------------------------------------
  if ((simbox.boundary_lhs[0] == openBoundary && simbox.boundary_rhs[0] == openBoundary) == 0) {

    for (i=0; i<mfv->Ntot; i++) {
      mfv->CheckXBoundaryGhostParticle(i,tghost,simbox);
    }

    mfv->Ntot = mfv->Nhydro + mfv->Nghost;
  }


  // Create ghost particles in y-dimension
  //-----------------------------------------------------------------------------------------------
  if (ndim >= 2 && (simbox.boundary_lhs[1] == openBoundary &&
                    simbox.boundary_rhs[1] == openBoundary) == 0) {

    for (i=0; i<mfv->Ntot; i++) {
      mfv->CheckYBoundaryGhostParticle(i,tghost,simbox);
    }

    mfv->Ntot = mfv->Nhydro + mfv->Nghost;
  }


  // Create ghost particles in z-dimension
  //-----------------------------------------------------------------------------------------------
  if (ndim == 3 && (simbox.boundary_lhs[2] == openBoundary &&
                    simbox.boundary_rhs[2] == openBoundary) == 0) {

    for (i=0; i<mfv->Ntot; i++) {
      mfv->CheckZBoundaryGhostParticle(i,tghost,simbox);
    }

    mfv->Ntot = mfv->Nhydro + mfv->Nghost;
  }


  // Quit here if we've run out of memory for ghosts
  if (mfv->Ntot > mfv->Nhydromax) {
    string message="Not enough memory for ghost particles";
    ExceptionHandler::getIstance().raise(message);
  }

  mfv->NPeriodicGhost = mfv->Nghost;

  return;
}*/



template class MeshlessFVBruteForce<1,MeshlessFVParticle>;
template class MeshlessFVBruteForce<2,MeshlessFVParticle>;
template class MeshlessFVBruteForce<3,MeshlessFVParticle>;
