//=============================================================================
//  BruteForceSearch.cpp
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
//=============================================================================



#include <iostream>
#include <math.h>
#include "SphNeighbourSearch.h"
#include "Sph.h"
#include "Parameters.h"
#include "SphParticle.h"
#include "Debug.h"
#include "InlineFuncs.h"
#include "SphKernel.h"
#if defined MPI_PARALLEL
#include "MpiNode.h"
#endif
using namespace std;



//=============================================================================
//  BruteForceSearch::BruteForceSearch
/// BruteForceSearch class constructor
//=============================================================================
template <int ndim, template<int> class ParticleType>
BruteForceSearch<ndim,ParticleType>::BruteForceSearch
(FLOAT kernrangeaux,
 DomainBox<ndim> *boxaux,
 SphKernel<ndim> *kernaux,
 CodeTiming *timingaux):
  SphNeighbourSearch<ndim>(kernrangeaux,boxaux,kernaux,timingaux)
{
}



//=============================================================================
//  BruteForceSearch::~BruteForceSearch
/// BruteForceSearch class destructor
//=============================================================================
template <int ndim, template<int> class ParticleType>
BruteForceSearch<ndim,ParticleType>::~BruteForceSearch()
{
}



//=============================================================================
//  BruteForceSearch::BuildTree
/// For Brute Force neighbour searching, there is no tree to construct but 
/// we chose to delete any dead SPH particles here to be consistent with 
/// the tree neighbour search.
//=============================================================================
template <int ndim, template<int> class ParticleType>
void BruteForceSearch<ndim,ParticleType>::BuildTree
(bool rebuild_tree,                 ///< Flag to rebuild tree
 int n,                             ///< Integer time
 int ntreebuildstep,                ///< Tree build frequency
 int ntreestockstep,                ///< Tree stocking frequency
 int Npart,                         ///< No. of particles
 int Npartmax,                      ///< Max. no. of particles
 SphParticle<ndim> *sph_gen, Sph<ndim> *sph,        ///< Particle data array
 FLOAT timestep)                    ///< Smallest physical timestep
{
  sph->DeleteDeadParticles();
  return;
}



//=============================================================================
//  BruteForceSearch::UpdateActiveParticleCounters
/// For Brute Force neighbour searching, there are no counters.
//=============================================================================
template <int ndim, template<int> class ParticleType>
void BruteForceSearch<ndim,ParticleType>::UpdateActiveParticleCounters
(SphParticle<ndim> *sph_gen, Sph<ndim> *sph)
{
  return;
}



//=============================================================================
//  BruteForceSearch::UpdateAllSphProperties
/// ..
//=============================================================================
template <int ndim, template<int> class ParticleType>
int BruteForceSearch<ndim,ParticleType>::GetGatherNeighbourList
(FLOAT rp[ndim],                    ///< ..
 FLOAT rsearch,                     ///< ..
 SphParticle<ndim> *sph_gen,        ///< ..
 int Nsph,                          ///< ..
 int Nneibmax,                      ///< ..
 int *neiblist)                     ///< ..
{
  int i,k;                          // ..
  int Nneib = 0;                    // No. of (non-dead) neighbours
  FLOAT dr[ndim];                   // Relative distance vector
  FLOAT drsqd;                      // Distance squared
  ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sph_gen);

  debug2("[BruteForceSearch::GetGatherNeighbourList]");

  // Compute smoothing lengths of all SPH particles
  //---------------------------------------------------------------------------
  for (i=0; i<Nsph; i++) {

    // Skip over inactive particles
    if (!sphdata[i].active || sphdata[i].itype == dead) continue;

    for (k=0; k<ndim; k++) dr[k] = sphdata[i].r[k] - rp[k];
    drsqd = DotProduct(dr,dr,ndim);

    if (drsqd <= rsearch*rsearch && Nneib < Nneibmax)
      neiblist[Nneib++] = i;
    else if (drsqd <= rsearch*rsearch && Nneib == Nneibmax)
      return -1;

  }
  //---------------------------------------------------------------------------

  return Nneib;
}



//=============================================================================
//  BruteForceSearch::UpdateAllSphProperties
/// Routine for computing SPH properties (smoothing lengths, densities and 
/// forces) for all active SPH particle using neighbour lists generated 
/// using brute force (i.e. direct summation).
//=============================================================================
template <int ndim, template<int> class ParticleType>
void BruteForceSearch<ndim,ParticleType>::UpdateAllSphProperties
(int Nsph,                          ///< [in] No. of SPH particles
 int Ntot,                          ///< [in] No. of SPH + ghost particles
 SphParticle<ndim> *sph_gen,        ///< [inout] Pointer to SPH ptcl array
 Sph<ndim> *sph,                    ///< [in] Pointer to SPH object
 Nbody<ndim> *nbody)                ///< [in] Pointer to N-body object
{
  int i,j,jj,k;                     // Particle and dimension counters
  int Nneib = 0;                    // No. of (non-dead) neighbours
  int okflag;                       // Flag valid smoothing length
  int *neiblist;                    // List of neighbours
  FLOAT dr[ndim];                   // Relative distance vector
  FLOAT rp[ndim];                   // Position of current particle
  FLOAT *drsqd;                     // Distance squared
  FLOAT *gpot;                      // Array of neib. grav. potentials
  FLOAT *m;                         // Array of neib. position vectors
  FLOAT *mu;                        // Array of neib. mass*u values
  ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sph_gen);

  debug2("[BruteForceSearch::UpdateAllSphProperties]");

  // Store masses in separate array
  gpot = new FLOAT[Ntot];
  m = new FLOAT[Ntot];
  mu = new FLOAT[Ntot];
  neiblist = new int[Ntot];
  for (i=0; i<Ntot; i++) {
    if (sphdata[i].itype == dead) continue;
    neiblist[Nneib] = i;
    gpot[Nneib] = sphdata[i].gpot;
    m[Nneib] = sphdata[i].m;
    mu[Nneib] = sphdata[i].m*sphdata[i].u;
    Nneib++;
  }

  // Create parallel threads
  //===========================================================================
#pragma omp parallel default(none) private(dr,drsqd,i,j,jj,k,okflag,rp)	\
  shared(gpot,m,mu,nbody,neiblist,Nneib,Nsph,Ntot,sph,sphdata)
  {
    drsqd = new FLOAT[Ntot];

    // Compute smoothing lengths of all SPH particles
    //-------------------------------------------------------------------------
#pragma omp for
    for (i=0; i<Nsph; i++) {

      // Skip over inactive particles
      if (!sphdata[i].active || sphdata[i].itype == dead) continue;

      for (k=0; k<ndim; k++) rp[k] = sphdata[i].r[k];

      // Compute distances and the reciprical between the current particle 
      // and all neighbours here
      //-----------------------------------------------------------------------
      for (jj=0; jj<Nneib; jj++) { 
	j = neiblist[jj];
    	for (k=0; k<ndim; k++) dr[k] = sphdata[j].r[k] - rp[k];
    	drsqd[jj] = DotProduct(dr,dr,ndim);
      }
      //-----------------------------------------------------------------------

      // Compute all SPH gather properties
      okflag = sph->ComputeH(i,Nneib,big_number,m,mu,drsqd,
                             gpot,sphdata[i],nbody);
  
    }
    //-------------------------------------------------------------------------

    delete[] drsqd;

  }
  //===========================================================================

  delete[] neiblist;
  delete[] mu;
  delete[] m;
  delete[] gpot;

  return;
}



//=============================================================================
//  BruteForceSearch::UpdateAllSphHydroForces
/// Routine for computing SPH properties (smoothing lengths, densities and 
/// forces) for all active SPH particle using neighbour lists generated 
/// using brute force (i.e. direct summation).
//=============================================================================
template <int ndim, template<int> class ParticleType>
void BruteForceSearch<ndim,ParticleType>::UpdateAllSphHydroForces
(int Nsph,                          ///< [in] No. of SPH particles
 int Ntot,                          ///< [in] No. of SPH + ghost particles
 SphParticle<ndim> *sph_gen,        ///< [inout] Pointer to SPH ptcl array
 Sph<ndim> *sph,                    ///< [in] Pointer to SPH object
 Nbody<ndim> *nbody)                ///< [in] Pointer to N-body object
{
  int i,j,k;                        // Particle and dimension counters
  int Nneib;                        // No. of neighbours
  int *neiblist;                    // List of neighbour ids
  FLOAT draux[ndim];                // Relative distance vector
  FLOAT drsqd;                      // Distance squared
  FLOAT hrangesqdi;                 // Gather kernel extent (squared)
  FLOAT hrangesqdj;                 // Scatter kernel extent (squared)
  FLOAT rp[ndim];                   // Position of current particle
  FLOAT *dr;                        // Array of neib. position vectors
  FLOAT *drmag;                     // Array of neib. distances
  FLOAT *invdrmag;                  // Array of neib. inverse distances
  ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sph_gen);

  debug2("[BruteForceSearch::UpdateAllSphHydroForces]");

  // Allocate memory for storing neighbour ids and position data
  neiblist = new int[Ntot];
  dr = new FLOAT[ndim*Ntot];
  drmag = new FLOAT[Ntot];
  invdrmag = new FLOAT[Ntot];


  // Compute smoothing lengths of all SPH particles
  //---------------------------------------------------------------------------
  for (i=0; i<Nsph; i++) {

    // Skip over inactive particles
    if (!sphdata[i].active || sphdata[i].itype == dead) continue;

    // Zero all arrays to be updated
    for (k=0; k<ndim; k++) sphdata[i].a[k] = (FLOAT) 0.0;
    for (k=0; k<ndim; k++) sphdata[i].agrav[k] = (FLOAT) 0.0;
    sphdata[i].gpot = (FLOAT) 0.0;
    sphdata[i].gpe = (FLOAT) 0.0;
    sphdata[i].dudt = (FLOAT) 0.0;
    sphdata[i].levelneib = 0;

    for (k=0; k<ndim; k++) rp[k] = sphdata[i].r[k];
    hrangesqdi = pow(kernfac*kernp->kernrange*sphdata[i].h,2);
    Nneib = 0;

    // Compute distances and the reciprical between the current particle 
    // and all neighbours here
    //-------------------------------------------------------------------------
    for (j=0; j<Ntot; j++) {
      if (sphdata[j].itype == dead) continue;
      hrangesqdj = pow(kernfac*kernp->kernrange*sphdata[j].h,2);
      for (k=0; k<ndim; k++) draux[k] = sphdata[j].r[k] - rp[k];
      drsqd = DotProduct(draux,draux,ndim);
      if ((drsqd < hrangesqdi || drsqd < hrangesqdj) && i != j) {
    	neiblist[Nneib] = j;
    	drmag[Nneib] = sqrt(drsqd);
    	invdrmag[Nneib] = (FLOAT) 1.0/(drmag[Nneib] + small_number);
    	for (k=0; k<ndim; k++) dr[Nneib*ndim + k] = draux[k]*invdrmag[Nneib];
    	Nneib++;
      }
    }
    //-------------------------------------------------------------------------

    // Compute all SPH hydro forces
    sph->ComputeSphHydroForces(i,Nneib,neiblist,drmag,invdrmag,dr,
			       sphdata[i],sphdata);

    // Compute all star forces
    sph->ComputeStarGravForces(nbody->Nnbody,nbody->nbodydata,sphdata[i]);

    sphdata[i].active = false;

  }
  //---------------------------------------------------------------------------

  
  // Free all allocated memory
  delete[] invdrmag;
  delete[] drmag;
  delete[] dr;
  delete[] neiblist;


  return;
}



//=============================================================================
//  BruteForceSearch::UpdateAllSphForces
/// Empty function for now
//=============================================================================
template <int ndim, template<int> class ParticleType>
void BruteForceSearch<ndim,ParticleType>::UpdateAllSphForces
(int Nsph,                            ///< [in] ..
 int Ntot,                            ///< [in] ..
 SphParticle<ndim> *sph_gen,          ///< [in] ..
 Sph<ndim> *sph,                      ///< [inout] Pointer to SPH object
 Nbody<ndim> *nbody)                  ///< [in] Pointer to N-body object
{
  int i,j,k;                          // Particle and dimension counters
  int Nneib;                          // No. of neighbours
  int *neiblist;                      // List of neighbour ids
  ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sph_gen);

  debug2("[BruteForceSearch::UpdateAllSphForces]");

  // Allocate memory for storing neighbour ids and position data
  neiblist = new int[Ntot];


  // Compute smoothing lengths of all SPH particles
  //---------------------------------------------------------------------------
  for (i=0; i<Nsph; i++) {

    // Skip over inactive particles
    if (!sphdata[i].active || sphdata[i].itype == dead) continue;

    // Zero all arrays to be updated
    for (k=0; k<ndim; k++) sphdata[i].a[k] = (FLOAT) 0.0;
    for (k=0; k<ndim; k++) sphdata[i].agrav[k] = (FLOAT) 0.0;
    sphdata[i].gpot = (FLOAT) 0.0;
    sphdata[i].gpe = (FLOAT) 0.0;
    sphdata[i].dudt = (FLOAT) 0.0;
    sphdata[i].levelneib = 0;

    // Add self-contribution to gravitational potential
    sphdata[i].gpot += sphdata[i].m*
      sphdata[i].invh*kernp->wpot(0.0);

    // Determine interaction list (to ensure we don't compute pair-wise
    // forces twice)
    Nneib = 0;
    for (j=0; j<Nsph; j++)
      if (i != j && sphdata[j].itype != dead) neiblist[Nneib++] = j;

    // Compute forces between SPH neighbours (hydro and gravity)
    sph->ComputeSphHydroGravForces(i,Nneib,neiblist,
                                   sphdata[i],sphdata);

    // Compute all star forces
    sph->ComputeStarGravForces(nbody->Nnbody,nbody->nbodydata,sphdata[i]);

    for (k=0; k<ndim; k++) sphdata[i].a[k] += sphdata[i].agrav[k];
    sphdata[i].active = false;

  }
  //---------------------------------------------------------------------------

  delete[] neiblist;

  return;
}



//=============================================================================
//  BruteForceSearch::UpdateAllSphGravForces
/// Routine for computing SPH properties (smoothing lengths, densities and 
/// forces) for all active SPH particle using neighbour lists generated 
/// using brute force (i.e. direct summation).
//=============================================================================
template <int ndim, template<int> class ParticleType>
void BruteForceSearch<ndim,ParticleType>::UpdateAllSphGravForces
(int Nsph,                            ///< [in] ..
 int Ntot,                            ///< [in] ..
 SphParticle<ndim> *sph_gen,          ///< [in] ..
 Sph<ndim> *sph,                      ///< [inout] Pointer to SPH object
 Nbody<ndim> *nbody)                  ///< [in] Pointer to N-body object
{
  int i,j,k;                          // Particle and dimension counters
  int Nneib;                          // No. of neighbours
  int *neiblist;                      // List of neighbour ids
  ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sph_gen);

  debug2("[BruteForceSearch::UpdateAllSphGravForces]");

  // Allocate memory for storing neighbour ids and position data
  neiblist = new int[Ntot];

  // Compute smoothing lengths of all SPH particles
  //---------------------------------------------------------------------------
  for (i=0; i<Nsph; i++) {

    // Skip over inactive particles
    if (!sphdata[i].active || sphdata[i].itype == dead) continue;

    // Zero all arrays to be updated
    for (k=0; k<ndim; k++) sphdata[i].a[k] = (FLOAT) 0.0;
    for (k=0; k<ndim; k++) sphdata[i].agrav[k] = (FLOAT) 0.0;
    sphdata[i].gpot = (FLOAT) 0.0;
    sphdata[i].gpe = (FLOAT) 0.0;
    sphdata[i].dudt = (FLOAT) 0.0;
    sphdata[i].levelneib = 0;

    // Add self-contribution to gravitational potential
    sphdata[i].gpot += sphdata[i].m*
      sphdata[i].invh*kernp->wpot(0.0);

    // Determine interaction list (to ensure we don't compute pair-wise
    // forces twice)
    Nneib = 0;
    for (j=0; j<Nsph; j++)
      if (i != j && sphdata[j].itype != dead) neiblist[Nneib++] = j;

    // Compute forces between SPH neighbours (hydro and gravity)
    sph->ComputeSphGravForces(i,Nneib,neiblist,sphdata[i],sphdata);

    // Compute all star forces
    sph->ComputeStarGravForces(nbody->Nnbody,nbody->nbodydata,sphdata[i]);

    for (k=0; k<ndim; k++) sphdata[i].a[k] += sphdata[i].agrav[k];
    sphdata[i].active = false;

  }
  //---------------------------------------------------------------------------

  delete[] neiblist;

  return;
}



//=============================================================================
//  BruteForceSearch::UpdateAllSphDerivatives
/// Compute all SPH derivatives required for 2nd-order Riemann solver in 
/// Godunov SPH method.
//=============================================================================
template <int ndim, template<int> class ParticleType>
void BruteForceSearch<ndim,ParticleType>::UpdateAllSphDerivatives
(int Nsph,                            ///< [in] ..
 int Ntot,                            ///< [in] ..
 SphParticle<ndim> *sph_gen, Sph<ndim> *sph)          ///< [inout] Pointer to SPH object
{
  int i,j,k;                          // Particle and dimension counters
  int Nneib;                          // No. of neighbours
  int *neiblist;                      // List of neighbour ids
  FLOAT draux[ndim];                  // Relative distance vector
  FLOAT drsqd;                        // Distance squared
  FLOAT hrangesqd;                    // Kernel extent (squared)
  FLOAT rp[ndim];                     // Position of current particle
  FLOAT *dr;                          // Array of neib. position vectors
  FLOAT *drmag;                       // Array of neib. distances
  FLOAT *invdrmag;                    // Array of neib. inverse distances
  struct SphParticle<ndim> *neibpart; // Local copies of neib. particles
  ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sph_gen);

  debug2("[BruteForceSearch::UpdateAllSphForces]");

  // The potential number of neighbours is given by ALL the particles
  Nneib = Ntot;

  // Allocate memory for storing neighbour ids and position data
  neiblist = new int[Ntot];
  dr = new FLOAT[ndim*Ntot];
  drmag = new FLOAT[Ntot];
  invdrmag = new FLOAT[Ntot];
  neibpart = new SphParticle<ndim>[Ntot];

  // Record local copies of (all) neighbour properties
  for (j=0; j<Ntot; j++) neibpart[j] = sphdata[j];


  // Compute smoothing lengths of all SPH particles
  //---------------------------------------------------------------------------
  for (i=0; i<Nsph; i++) {
    for (k=0; k<ndim; k++) rp[k] = sphdata[i].r[k];
    hrangesqd = pow(kernp->kernrange*sphdata[i].h,2);
    Nneib = 0;

    // Compute distances and the reciprical between the current particle 
    // and all neighbours here
    //-------------------------------------------------------------------------
    for (j=0; j<Ntot; j++) {
      for (k=0; k<ndim; k++) draux[k] = sphdata[j].r[k] - rp[k];
      drsqd = DotProduct(draux,draux,ndim);
      if (drsqd < hrangesqd) {
    	neiblist[Nneib] = j;
    	drmag[Nneib] = sqrt(drsqd);
    	invdrmag[Nneib] = (FLOAT) 1.0/(drmag[Nneib] + small_number);
    	for (k=0; k<ndim; k++) dr[Nneib*ndim + k] = draux[k]*invdrmag[Nneib];
    	Nneib++;
      }
    }
    //-------------------------------------------------------------------------

    // Compute all SPH hydro forces
    sph->ComputeSphDerivatives(i,Nneib,neiblist,drmag,invdrmag,dr,
			       sphdata[i],neibpart);

  }
  //---------------------------------------------------------------------------

  delete[] neibpart;
  delete[] invdrmag;
  delete[] drmag;
  delete[] dr;
  delete[] neiblist;

  return;
}



//=============================================================================
//  BruteForceSearch::UpdateAllSphDudt
/// Compute the compressional heating rate (dudt) for all active particles.
//=============================================================================
template <int ndim, template<int> class ParticleType>
void BruteForceSearch<ndim,ParticleType>::UpdateAllSphDudt
(int Nsph,                            ///< [in] ..
 int Ntot,                            ///< [in] ..
 SphParticle<ndim> *sph_gen, Sph<ndim> *sph)          ///< [inout] Pointer to SPH object
{
  int i,j,k;                          // Particle and dimension counters
  int Nneib;                          // No. of neighbours
  int *neiblist;                      // List of neighbour ids
  FLOAT draux[ndim];                  // Relative distance vector
  FLOAT drsqd;                        // Distance squared
  FLOAT hrangesqdi;                   // Gather kernel range (squared)
  FLOAT hrangesqdj;                   // Scatter kernel range (squared)
  FLOAT rp[ndim];                     // Position of current particle
  FLOAT *dr;                          // Array of neib. position vectors
  FLOAT *drmag;                       // Array of neib. distances
  FLOAT *invdrmag;                    // Array of neib. inverse distances
  struct SphParticle<ndim> *neibpart; // Local copies of neighbour particles
  ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sph_gen);

  debug2("[BruteForceSearch::UpdateAllSphDudt]");

  // The potential number of neighbours is given by ALL the particles
  Nneib = Ntot;

  // Allocate memory for storing neighbour ids and position data
  neiblist = new int[Ntot];
  dr = new FLOAT[ndim*Ntot];
  drmag = new FLOAT[Ntot];
  invdrmag = new FLOAT[Ntot];
  neibpart = new SphParticle<ndim>[Ntot];

  for (j=0; j<Ntot; j++) neibpart[j] = sphdata[j];
  for (j=0; j<Ntot; j++) neibpart[j].dudt = (FLOAT) 0.0;


  // Compute smoothing lengths of all SPH particles
  //---------------------------------------------------------------------------
  for (i=0; i<Nsph; i++) {
    for (k=0; k<ndim; k++) rp[k] = sphdata[i].r[k];
    hrangesqdi = pow(kernfac*kernp->kernrange*sphdata[i].h,2);
    Nneib = 0;

    // Compute distances and the reciprical between the current particle 
    // and all neighbours here
    //-------------------------------------------------------------------------
    for (j=0; j<Ntot; j++) {
      hrangesqdj = pow(kernfac*kernp->kernrange*sphdata[j].h,2);
      for (k=0; k<ndim; k++) draux[k] = sphdata[j].r[k] - rp[k];
      drsqd = DotProduct(draux,draux,ndim);
      if ((drsqd < hrangesqdi || drsqd < hrangesqdj) &&
	  ((j < i && !sphdata[j].active) || j > i)) {
    	neiblist[Nneib] = j;
    	drmag[Nneib] = sqrt(drsqd);
    	invdrmag[Nneib] = (FLOAT) 1.0/(drmag[Nneib] + small_number);
    	for (k=0; k<ndim; k++) dr[Nneib*ndim + k] = draux[k]*invdrmag[Nneib];
    	Nneib++;
      }
    }
    //-------------------------------------------------------------------------

    // Compute all SPH hydro forces
    sph->ComputeSphNeibDudt(i,Nneib,neiblist,drmag,invdrmag,dr,
			    sphdata[i],neibpart);

  }
  //---------------------------------------------------------------------------


  // Now add all active neighbour contributions to the main arrays
  for (j=0; j<Ntot; j++) {
    if (neibpart[j].active) {
      sphdata[j].dudt += neibpart[j].dudt;
    }
  }

  delete[] neibpart;
  delete[] invdrmag;
  delete[] drmag;
  delete[] dr;
  delete[] neiblist;


  return;
}



//=============================================================================
//  BruteForceSearch::UpdateAllStarGasForces
/// ..
//=============================================================================
template <int ndim, template<int> class ParticleType>
void BruteForceSearch<ndim,ParticleType>::UpdateAllStarGasForces
(int Nsph,                            ///< [in] ..
 int Ntot,                            ///< [in] No. of SPH particles
 SphParticle<ndim> *sph_gen, Sph<ndim> *sph,          ///< [inout] Pointer to SPH ptcl array
 Nbody<ndim> *nbody)                  ///< [inout] Pointer to N-body object
{
  int i;                              // Particle and dimension counters
  int Nneib = 0;                      // No. of (non-dead) neighbours
  int *dummy;                         // Dummy var to satisfy function argument
  int *neiblist;                      // Array of (all) particle ids
  ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sph_gen);

  debug2("[BruteForceSearch::UpdateAllSphForces]");

  // Allocate memory for storing neighbour ids and position data
  neiblist = new int[Nsph];
  for (i=0; i<Nsph; i++) 
    if (sphdata[i].itype != dead) neiblist[Nneib++] = i;

  // Compute smoothing lengths of all SPH particles
  //---------------------------------------------------------------------------
  for (i=0; i<nbody->Nnbody; i++) {

    // Skip over inactive particles
    if (!nbody->nbodydata[i]->active) continue;

    // Compute forces between SPH neighbours (hydro and gravity)
    nbody->CalculateDirectSPHForces(nbody->nbodydata[i],Nneib,
                                    0,neiblist,dummy,sph);

  }
  //---------------------------------------------------------------------------

  delete[] neiblist;

  return;
}



#if defined MPI_PARALLEL
//=============================================================================
//  BruteForceSearch::FindGhostParticlesToExport
/// Compute on behalf of the MpiControl class the ghost particles we need 
/// to export to other nodes.
//=============================================================================
template <int ndim, template<int> class ParticleType>
void BruteForceSearch<ndim,ParticleType>::FindGhostParticlesToExport
(Sph<ndim>* sph,                            ///< [in] Pointer to sph class
 vector<vector<ParticleType<ndim>*>>& particles_to_export_per_node, ///< [inout] Vector that will be filled with values
 const std::vector<int>& overlapping_nodes, ///< [in] Vector containing which 
                                            ///<      nodes overlap our hbox
 MpiNode<ndim>* mpinodes)                   ///< [in] Array of other mpi nodes
{
  int i;                                    // ..
  int inode;                                // ..
  int node_number;                          // ..
  ParticleType<ndim> *sphdata = 
    static_cast<ParticleType<ndim>*> sph->GetParticlesArray();

  // Loop over particles and prepare the ones to export
  for (i=0; i<sph->Ntot; i++) {
    ParticleType<ndim>& part = sphdata[i];
    
    // Loop over potential domains and find particles to export to them
    for (inode=0; inode<overlapping_nodes.size(); inode++) {
      node_number = overlapping_nodes[inode];
      if (ParticleBoxOverlap(part,mpinodes[node_number].hbox)) {
        particles_to_export_per_node[node_number].push_back(&part);
      }
    }
  }
  
  return;
}



//=============================================================================
//  BruteForceSearch::FindParticlesToTransfer
/// Compute on behalf of the MpiControl class the particles that are outside 
/// the domain after a load balancing and need to be transferred to other nodes
//=============================================================================
template <int ndim, template<int> class ParticleType>
void BruteForceSearch<ndim,ParticleType>::FindParticlesToTransfer
(Sph<ndim>* sph,    ///< [in] Pointer to sph class
 std::vector<std::vector<int> >& particles_to_export, ///< [inout] Vector that for each node gives the list of particles to export
 std::vector<int>& all_particles_to_export,  ///< [inout] Vector containing all the particles that will be exported by this processor
 const std::vector<int>& potential_nodes, ///< [in] Vector containing the potential nodes we might be sending particles to
 MpiNode<ndim>* mpinodes) ///< [in] Array of other mpi nodes
{

  ParticleType<ndim> *sphdata = static_cast<ParticleType<ndim>* > 
    (sph->GetParticlesArray());

  //Loop over particles and prepare the ones to export
  for (int i=0; i<sph->Nsph; i++) {
    ParticleType<ndim>& part = sphdata[i];

    //Loop over potential domains and see if we need to transfer this particle to them
    for (int inode=0; inode<potential_nodes.size(); inode++) {
      int node_number = potential_nodes[inode];
      if (ParticleInBox(part,mpinodes[node_number].domain)) {
        particles_to_export[node_number].push_back(i);
        all_particles_to_export.push_back(i);
        // The particle can belong only to one domain, so we can break from this loop
        break;
      }
    }
  }
}


//=============================================================================
//  BruteForceSearch::GetExportInfo
/// Get the array with the information that needs to be exported to the given
/// processor (note: Nproc is ignored at the moment, as we always need to export
/// all particles to the other processors)
//=============================================================================
template <int ndim, template<int> class ParticleType>
vector<ExportParticleInfo<ndim> >& BruteForceSearch<ndim,ParticleType>::GetExportInfo (
    int Nproc,        ///< [in] Number of processor we want to send the information to
    Sph<ndim>*  sph ) ///< [in] Pointer to sph object
    {

  //Find number of active particles
  ids_active_particles.clear();
  ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sph->GetParticlesArray() );
  for (int i=0; i< sph->Nsph; i++) {
    if (sphdata[i].active) {
      ids_active_particles.push_back(i);
    }
  }

  const int Nactive = ids_active_particles.size();

  //Copy positions of active particles inside arrays
  int j=0;
  for (int i=0; i< sph->Nsph; i++) {
    if (sphdata[i].active) {
      for (int k=0; k<ndim; k++)
        particles_to_export[j].r[k] = sphdata[i].r[k];
      j++;
    }
  }

  assert(j==Nactive);

  return particles_to_export;

}


//=============================================================================
//  BruteForceSearch::UnpackExported
/// Unpack the information exported from the other processors, contaning the positions
/// of the particles that were exported
//=============================================================================
template <int ndim, template<int> class ParticleType>
void BruteForceSearch<ndim,ParticleType>::UnpackExported (
    vector<ExportParticleInfo<ndim> >& received_array,
    vector<int>& N_received_particles_from_proc,
    Sph<ndim>* sph,
    int rank) {


  int offset = N_received_particles_from_proc[0];

  ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sph->GetParticlesArray() );

  for (int Nproc = 0; Nproc<N_received_particles_from_proc.size(); Nproc++) {

    int N_received_particles = N_received_particles_from_proc[Nproc];

    offset += N_received_particles;

    //Avoid inserting my own particles
    if (Nproc==rank)
      continue;

    //Ensure there is enough memory
    if (sph->Nsph + N_received_particles > sph->Nsphmax) {
      ExceptionHandler::getIstance().raise("Error while receiving imported particles: not enough memory!");
    }

    //Insert particles inside SPH main arrays
    for (int i=0; i<N_received_particles; i++) {
      for (int k=0; k<ndim; k++)
        sphdata[i+sph->Nsph].r[k] = received_array[i+offset].r[k];
    }

    //Update the SPH counter (note: there is probably more to do here!!!!)
    sph->Nsph += N_received_particles;

  }

}


//=============================================================================
//  BruteForceSearch::GetBackExportInfo
/// Return the data to transmit back to the other processors (particle acceleration etc.)
//=============================================================================
template <int ndim, template<int> class ParticleType>
void BruteForceSearch<ndim,ParticleType>::GetBackExportInfo (
    vector<ExportBackParticleInfo<ndim> >& send_buffer, ///< [inout] These arrays will be overwritten with the information to send
    vector<int>& N_exported_particles_from_proc,
    Sph<ndim>* sph,   ///< [in] Pointer to the SPH object
    int rank
    ) {

  //loop BACKWARDS over the processors, removing particles as we go
  int removed_particles=0;
  for (int Nproc=N_exported_particles_from_proc.size()-1; Nproc>=0; Nproc-- ) {

    //My local particles were not inserted
    if (rank==Nproc)
      continue;

    const int N_received_particles = N_exported_particles_from_proc[Nproc];

    //Copy the accelerations and gravitational potential of the particles
    ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sph->GetParticlesArray() );
    int j=0;
    for (int i=sph->Nsph - N_received_particles; i<sph->Nsph; i++) {
      for (int k=0; k<ndim; k++)
        send_buffer[removed_particles+j].a[k] = sphdata[i].a[k];
      send_buffer[removed_particles+j].gpot = sphdata[i].gpot;
      j++;
    }

    assert(j==N_received_particles);

    removed_particles += j;

    //Decrease the particle counter (note: there is probably more to do here!!!!)
    sph->Nsph -= N_received_particles;

  }

}


//=============================================================================
//  BruteForceSearch::UnpackReturnedExportInfo
/// Unpack the data that was returned by the other processors, summing the accelerations to the particles
//=============================================================================
template <int ndim, template<int> class ParticleType>
void BruteForceSearch<ndim,ParticleType>::UnpackReturnedExportInfo (
    vector<ExportBackParticleInfo<ndim> >& received_information,
    vector<int>& recv_displs,
    Sph<ndim>* sph,
    int rank
    ) {

  ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sph->GetParticlesArray() );

  //For each particle, sum up the accelerations returned by other processors
  for (int i=0; i< ids_active_particles.size(); i++ ) {
    const int j = ids_active_particles[i];

    for(int Nproc=0; Nproc<recv_displs.size(); Nproc++ ) {

      if (rank==Nproc)
        continue;

      for (int k=0; k<ndim; k++)
        sphdata[j].a[k] += received_information[i+recv_displs[Nproc] ].a[k];
      sphdata[j].gpot += received_information[i+recv_displs[Nproc] ].gpot;
    }

  }

}



#endif



template class BruteForceSearch<1,GradhSphParticle>;
template class BruteForceSearch<2,GradhSphParticle>;
template class BruteForceSearch<3,GradhSphParticle>;
template class BruteForceSearch<1,SM2012SphParticle>;
template class BruteForceSearch<2,SM2012SphParticle>;
template class BruteForceSearch<3,SM2012SphParticle>;
template class BruteForceSearch<1,GodunovSphParticle>;
template class BruteForceSearch<2,GodunovSphParticle>;
template class BruteForceSearch<3,GodunovSphParticle>;
