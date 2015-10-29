//=================================================================================================
//  SphBruteForceSearch.cpp
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



#include <assert.h>
#include <iostream>
#include <math.h>
#include <numeric>
#include "MirrorNeighbours.hpp"
#include "NeighbourSearch.h"
#include "SphNeighbourSearch.h"
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
//  SphBruteForceSearch::SphBruteForceSearch
/// SphBruteForceSearch class constructor
//=================================================================================================
template <int ndim, template<int> class ParticleType>
SphBruteForceSearch<ndim,ParticleType>::SphBruteForceSearch
 (FLOAT kernrangeaux,
  DomainBox<ndim> *boxaux,
  SmoothingKernel<ndim> *kernaux,
  CodeTiming *timingaux):
  NeighbourSearch<ndim>(kernrangeaux, boxaux, kernaux, timingaux),
  SphNeighbourSearch<ndim>(kernrangeaux, boxaux, kernaux, timingaux),
  BruteForceSearch<ndim,ParticleType>(kernrangeaux, boxaux, kernaux, timingaux)
{
}



//=================================================================================================
//  SphBruteForceSearch::~SphBruteForceSearch
/// SphBruteForceSearch class destructor
//=================================================================================================
template <int ndim, template<int> class ParticleType>
SphBruteForceSearch<ndim,ParticleType>::~SphBruteForceSearch()
{
}



//=================================================================================================
//  SphBruteForceSearch::UpdateAllSphProperties
/// Routine for computing SPH properties (smoothing lengths, densities and
/// forces) for all active SPH particle using neighbour lists generated
/// using brute force (i.e. direct summation).
//=================================================================================================
template <int ndim, template<int> class ParticleType>
void SphBruteForceSearch<ndim,ParticleType>::UpdateAllSphProperties
 (int Nhydro,                          ///< [in] No. of SPH particles
  int Ntot,                            ///< [in] No. of SPH + ghost particles
  SphParticle<ndim> *part_gen,         ///< [inout] Pointer to SPH ptcl array
  Sph<ndim> *sph,                      ///< [in] Pointer to SPH object
  Nbody<ndim> *nbody)                  ///< [in] Pointer to N-body object
{
  int i,j,jj,k;                        // Particle and dimension counters
  int Nneib = 0;                       // No. of (non-dead) neighbours
  int *neiblist;                       // List of neighbours
  FLOAT dr[ndim];                      // Relative distance vector
  FLOAT rp[ndim];                      // Position of current particle
  FLOAT *drsqd;                        // Distance squared
  FLOAT *gpot;                         // Array of neib. grav. potentials
  FLOAT *m;                            // Array of neib. position vectors
  FLOAT *mu;                           // Array of neib. mass*u values
  ParticleType<ndim>* partdata = static_cast<ParticleType<ndim>* > (part_gen);

  debug2("[SphBruteForceSearch::UpdateAllSphProperties]");

  // Store masses in separate array
  gpot = new FLOAT[Ntot];
  m = new FLOAT[Ntot];
  mu = new FLOAT[Ntot];
  neiblist = new int[Ntot];
  for (i=0; i<Ntot; i++) {
    if (partdata[i].itype == dead) continue;
    neiblist[Nneib] = i;
    gpot[Nneib] = partdata[i].gpot;
    m[Nneib]    = partdata[i].m;
    mu[Nneib]   = partdata[i].m*partdata[i].u;
    Nneib++;
  }

  // Create parallel threads
  //===============================================================================================
#pragma omp parallel default(none) private(dr,drsqd,i,j,jj,k,rp)	\
  shared(gpot,m,mu,nbody,neiblist,Nneib,Nhydro,Ntot,sph,partdata)
  {
    drsqd = new FLOAT[Ntot];

    // Compute smoothing lengths of all SPH particles
    //---------------------------------------------------------------------------------------------
#pragma omp for
    for (i=0; i<Nhydro; i++) {

      // Skip over inactive particles
      if (!partdata[i].active || partdata[i].itype == dead) continue;

      for (k=0; k<ndim; k++) rp[k] = partdata[i].r[k];

      // Compute distances and the reciprical between the current particle
      // and all neighbours here
      //-------------------------------------------------------------------------------------------
      for (jj=0; jj<Nneib; jj++) {
        j = neiblist[jj];
        for (k=0; k<ndim; k++) dr[k] = partdata[j].r[k] - rp[k];
        drsqd[jj] = DotProduct(dr, dr, ndim) + small_number;
      }
      //-------------------------------------------------------------------------------------------

      // Compute all SPH gather properties
      sph->ComputeH(i, Nneib, big_number, m, mu, drsqd, gpot, partdata[i], nbody);

    }
    //---------------------------------------------------------------------------------------------

    delete[] drsqd;

  }
  //===============================================================================================

  delete[] neiblist;
  delete[] mu;
  delete[] m;
  delete[] gpot;

  return;
}



//=================================================================================================
//  SphBruteForceSearch::UpdateAllSphHydroForces
/// Routine for computing SPH properties (smoothing lengths, densities and
/// forces) for all active SPH particle using neighbour lists generated
/// using brute force (i.e. direct summation).
//=================================================================================================
template <int ndim, template<int> class ParticleType>
void SphBruteForceSearch<ndim,ParticleType>::UpdateAllSphHydroForces
 (int Nhydro,                          ///< [in] No. of SPH particles
  int Ntot,                            ///< [in] No. of SPH + ghost particles
  SphParticle<ndim> *part_gen,         ///< [inout] Pointer to SPH ptcl array
  Sph<ndim> *sph,                      ///< [in] Pointer to SPH object
  Nbody<ndim> *nbody)                  ///< [in] Pointer to N-body object
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
  ParticleType<ndim>* partdata = static_cast<ParticleType<ndim>* > (part_gen);
  const int offset_imported = sph->Nghost;

  debug2("[SphBruteForceSearch::UpdateAllSphHydroForces]");

  // Allocate memory for storing neighbour ids and position data
  neiblist = new int[Ntot];
  dr       = new FLOAT[ndim*Ntot];
  drmag    = new FLOAT[Ntot];
  invdrmag = new FLOAT[Ntot];


  // Compute forces of real and imported particles
  //-----------------------------------------------------------------------------------------------
  for (int ipart=0; ipart<Nhydro+sph->NImportedParticles; ipart++) {

    if (ipart < Nhydro) i = ipart;
    else i = ipart + offset_imported;

    // Skip over inactive particles
    if (!partdata[i].active || partdata[i].itype == dead) continue;

    // Zero all arrays to be updated
    for (k=0; k<ndim; k++) partdata[i].a[k] = (FLOAT) 0.0;
    for (k=0; k<ndim; k++) partdata[i].agrav[k] = (FLOAT) 0.0;
    partdata[i].gpot      = (FLOAT) 0.0;
    partdata[i].dudt      = (FLOAT) 0.0;
    partdata[i].levelneib = 0;

    for (k=0; k<ndim; k++) rp[k] = partdata[i].r[k];
    hrangesqdi = partdata[i].hrangesqd;
    //hrangesqdi = pow(kernp->kernrange*partdata[i].h,2); //pow(kernfac*kernp->kernrange*partdata[i].h,2);
    Nneib = 0;

    // Compute distances and the reciprical between the current particle
    // and all neighbours here
    //---------------------------------------------------------------------------------------------
    for (j=0; j<sph->Nhydro + sph->NPeriodicGhost; j++) {
      if (partdata[j].itype == dead) continue;
      hrangesqdj = partdata[j].hrangesqd;
      //hrangesqdj = pow(kernfac*kernp->kernrange*partdata[j].h,2);
      for (k=0; k<ndim; k++) draux[k] = partdata[j].r[k] - rp[k];
      drsqd = DotProduct(draux,draux,ndim);

      if ((drsqd < hrangesqdi || drsqd < hrangesqdj) && i != j) {
        neiblist[Nneib] = j;
        drmag[Nneib] = sqrt(drsqd);
        invdrmag[Nneib] = (FLOAT) 1.0/(drmag[Nneib] + small_number);
        for (k=0; k<ndim; k++) dr[Nneib*ndim + k] = draux[k]*invdrmag[Nneib];
        Nneib++;
      }
    }
    //---------------------------------------------------------------------------------------------

    // Compute all SPH hydro forces
    sph->ComputeSphHydroForces(i,Nneib,neiblist,drmag,invdrmag,dr,partdata[i],partdata);

    // Compute all star forces
    if (ipart < Nhydro) {
      sph->ComputeStarGravForces(nbody->Nnbody,nbody->nbodydata,partdata[i]);
    }

    partdata[i].active = false;

  }
  //-----------------------------------------------------------------------------------------------


  // Free all allocated memory
  delete[] invdrmag;
  delete[] drmag;
  delete[] dr;
  delete[] neiblist;


  return;
}



//=================================================================================================
//  SphBruteForceSearch::UpdateAllSphForces
/// Update all SPH forces (both hydro and gravity).
//=================================================================================================
template <int ndim, template<int> class ParticleType>
void SphBruteForceSearch<ndim,ParticleType>::UpdateAllSphForces
 (int Nhydro,                          ///< [in] No. of SPH particles
  int Ntot,                            ///< [in] Total no. of particles
  SphParticle<ndim> *part_gen,         ///< [inout] Pointer to SPH data
  Sph<ndim> *sph,                      ///< [inout] Pointer to SPH object
  Nbody<ndim> *nbody)                  ///< [in] Pointer to N-body object
{
  int i,j,k;                           // Particle and dimension counters
  int Nneib;                           // No. of neighbours
  int *neiblist;                       // List of neighbour ids
  ParticleType<ndim>* partdata = static_cast<ParticleType<ndim>* > (part_gen);
  const int offset_imported = sph->Nghost;

  debug2("[SphBruteForceSearch::UpdateAllSphForces]");

  // Allocate memory for storing neighbour ids and position data
  neiblist = new int[Ntot];


  // Compute forces for real and imported particles
  //-----------------------------------------------------------------------------------------------
  for (int ipart=0; ipart<Nhydro+sph->NImportedParticles; ipart++) {

    if (ipart < Nhydro) i = ipart;
    else i = ipart + offset_imported;

    // Skip over inactive particles
    if (!partdata[i].active || partdata[i].itype == dead) continue;

    // Zero all arrays to be updated
    for (k=0; k<ndim; k++) partdata[i].a[k] = (FLOAT) 0.0;
    for (k=0; k<ndim; k++) partdata[i].agrav[k] = (FLOAT) 0.0;
    partdata[i].gpot      = (FLOAT) 0.0;
    partdata[i].dudt      = (FLOAT) 0.0;
    partdata[i].levelneib = 0;

    // Add self-contribution to gravitational potential
    partdata[i].gpot += partdata[i].m*partdata[i].invh*kernp->wpot(0.0);

    // Determine interaction list (to ensure we don't compute pair-wise forces twice)
    Nneib = 0;
    for (j=0; j<sph->Nhydro + sph->NPeriodicGhost; j++) {
      if (i != j && partdata[j].itype != dead) neiblist[Nneib++] = j;
    }

    // Compute forces between SPH neighbours (hydro and gravity)
    sph->ComputeSphHydroGravForces(i, Nneib, neiblist, partdata[i], partdata);

    // Compute all star forces
    if (ipart < Nhydro) sph->ComputeStarGravForces(nbody->Nnbody, nbody->nbodydata, partdata[i]);

    for (k=0; k<ndim; k++) partdata[i].a[k] += partdata[i].agrav[k];
    partdata[i].active = false;

  }
  //-----------------------------------------------------------------------------------------------

  delete[] neiblist;

  return;
}



//=================================================================================================
//  SphBruteForceSearch::UpdateAllSphGravForces
/// Update all SPH gravity forces.
//=================================================================================================
template <int ndim, template<int> class ParticleType>
void SphBruteForceSearch<ndim,ParticleType>::UpdateAllSphGravForces
 (int Nhydro,                          ///< [in] No. of hydro particles
  int Ntot,                            ///< [in] Total no. of particles (inc. ghosts)
  SphParticle<ndim> *part_gen,         ///< [in] Pointer to main SPH particle array
  Sph<ndim> *sph,                      ///< [inout] Pointer to SPH object
  Nbody<ndim> *nbody)                  ///< [in] Pointer to N-body object
{
  int i,j,k;                           // Particle and dimension counters
  int Nneib;                           // No. of neighbours
  int *neiblist;                       // List of neighbour ids
  ParticleType<ndim>* partdata = static_cast<ParticleType<ndim>* > (part_gen);
  const int offset_imported = sph->Nghost;

  debug2("[SphBruteForceSearch::UpdateAllSphGravForces]");

  // Allocate memory for storing neighbour ids and position data
  neiblist = new int[Ntot];


  // Compute forces for real and imported particles
  //-----------------------------------------------------------------------------------------------
  for (int iparticle=0; iparticle<Nhydro+sph->NImportedParticles; iparticle++) {

    if (iparticle < Nhydro) i = iparticle;
    else i = iparticle + offset_imported;

    // Skip over inactive particles
    if (!partdata[i].active || partdata[i].itype == dead) continue;

    // Zero all arrays to be updated
    for (k=0; k<ndim; k++) partdata[i].a[k] = (FLOAT) 0.0;
    for (k=0; k<ndim; k++) partdata[i].agrav[k] = (FLOAT) 0.0;
    partdata[i].gpot      = (FLOAT) 0.0;
    partdata[i].dudt      = (FLOAT) 0.0;
    partdata[i].levelneib = 0;

    // Add self-contribution to gravitational potential
    partdata[i].gpot += partdata[i].m*partdata[i].invh*kernp->wpot(0.0);

    // Determine interaction list (to ensure we don't compute pair-wise forces twice)
    Nneib = 0;
    for (j=0; j<Nhydro; j++) {
      if (i != j && partdata[j].itype != dead) neiblist[Nneib++] = j;
    }

    // Compute forces between SPH neighbours (hydro and gravity)
    sph->ComputeSphGravForces(i,Nneib,neiblist,partdata[i],partdata);

    // Compute all star forces
    if (iparticle < Nhydro) sph->ComputeStarGravForces(nbody->Nnbody,nbody->nbodydata,partdata[i]);

    for (k=0; k<ndim; k++) partdata[i].a[k] += partdata[i].agrav[k];
    partdata[i].active = false;

  }
  //-----------------------------------------------------------------------------------------------

  delete[] neiblist;

  return;
}


//=================================================================================================
//  ParticleCellProxy
/// Proxy class that we can use to create the mirrors for the brute force search
//=================================================================================================
template<int ndim>
struct ParticleCellProxy
{
	ParticleCellProxy(const Particle<ndim>& p)
	{
		FLOAT dr = 0.5 * sqrt(p.hrangesqd) ;
		for (int k=0; k<ndim;k++)
		{
			rcell[k] = p.r[k] ;

			hboxmax[k] = rcell[k] + dr ;
			hboxmax[k] = rcell[k] - dr ;
		}
	}
	FLOAT hboxmin[ndim], hboxmax[ndim] ;
	FLOAT rcell[ndim] ;
} ;


//=================================================================================================
//  DomainCellProxy
/// Proxy class that we can use to create the mirrors for the brute force search
//=================================================================================================
template<int ndim>
struct DomainCellProxy
{
	DomainCellProxy(const DomainBox<ndim>& box)
	{
		for (int k=0; k<ndim;k++)
		{
			rcell[k] = 0.5 * (box.boxmax[k] + box.boxmin[k]) ;
			FLOAT dr = 0.6 * (box.boxmax[k] - box.boxmin[k]) ;

			hboxmax[k] = rcell[k] + dr ;
			hboxmax[k] = rcell[k] - dr ;
		}
	}
	FLOAT hboxmin[ndim], hboxmax[ndim] ;
	FLOAT rcell[ndim] ;
} ;


//=================================================================================================
//  SphBruteForceSearch::UpdateAllSphPeriodicForces
/// Update all SPH forces (both hydro and gravity) for periodic boundary conditions.
//=================================================================================================
template <int ndim, template<int> class ParticleType>
void SphBruteForceSearch<ndim,ParticleType>::UpdateAllSphPeriodicForces
 (int Nhydro,                          ///< [in] No. of SPH particles
  int Ntot,                            ///< [in] Total no. of particles
  SphParticle<ndim> *part_gen,         ///< [in] Pointer to SPH particle array
  Sph<ndim> *sph,                      ///< [inout] Pointer to SPH object
  Nbody<ndim> *nbody,                  ///< [in] Pointer to N-body object
  DomainBox<ndim> &simbox,             ///< [in] Simulation box with periodic information
  Ewald<ndim> *ewald)                  ///< [in] Pointer to Ewald object for periodic forces
{
  int i,j,k;                           // Particle and dimension counters
  int Nneib;                           // No. of neighbours
  int *neiblist;                       // List of neighbour ids
  FLOAT potperiodic;                   // Periodic potential correction
  FLOAT aperiodic[ndim];               // Periodic acceleration correction
  FLOAT dr[ndim];                      // Relative displacement vector
  FLOAT dr_corr[ndim];                 // Periodic correction displacement
  ParticleType<ndim>* neibdata;        // Local copy of neighbouring particles
  ParticleType<ndim>* partdata = static_cast<ParticleType<ndim>* > (part_gen);


  debug2("[SphBruteForceSearch::UpdateAllSphPeriodicForces]");

  // Construct the mirrors
  MirrorNeighbourFinder<ndim> NgbFinder(simbox) ;
  NgbFinder.SetTargetCell(DomainCellProxy<ndim>(simbox)) ;

  FLOAT r_ghost[ndim] ;
  int   sign[ndim] ;

  assert(NgbFinder.MaxNumMirrors() == 1) ;

  // Allocate memory for storing neighbour ids and position data
  neiblist = new int[Ntot];
  neibdata = new ParticleType<ndim>[Ntot];

  // Compute smoothing lengths of all SPH particles
  //-----------------------------------------------------------------------------------------------
  for (i=0; i<Nhydro; i++) {
    // Skip over inactive particles
    if (!partdata[i].active || partdata[i].itype == dead) continue;

    // Zero all arrays to be updated
    for (k=0; k<ndim; k++) partdata[i].a[k] = (FLOAT) 0.0;
    for (k=0; k<ndim; k++) partdata[i].agrav[k] = (FLOAT) 0.0;

    partdata[i].gpot = (FLOAT) 0.0;
    partdata[i].dudt = (FLOAT) 0.0;
    partdata[i].levelneib = 0;

    // Add self-contribution to gravitational potential
    partdata[i].gpot += partdata[i].m*partdata[i].invh*kernp->wpot(0.0);

    // Determine interaction list (to ensure we don't compute pair-wise forces twice).
    // Also make sure that only the closest periodic replica is considered.
    Nneib = 0;
    for (j=0; j<Nhydro; j++) {
      neibdata[j] = partdata[j];
      if (i != j && partdata[j].itype != dead) {
        neiblist[Nneib++] = j;
        for (k=0; k<ndim; k++) dr[k] = neibdata[j].r[k] - partdata[i].r[k];
        NgbFinder.NearestPeriodicVector(dr);
        for (k=0; k<ndim; k++) neibdata[j].r[k] = partdata[i].r[k] + dr[k];
      };
    }


    // Compute forces between SPH neighbours (hydro and gravity)
    sph->ComputeSphHydroGravForces(i,Nneib,neiblist,partdata[i],partdata);

    if (simbox.PeriodicGravity){
      // Now add the periodic correction force
      for (j=0; j<Nneib; j++) {
        for (k=0; k<ndim; k++) dr[k] = neibdata[j].r[k] - partdata[i].r[k];
        ewald->CalculatePeriodicCorrection(neibdata[j].m,dr,aperiodic,potperiodic);
        for (k=0; k<ndim; k++) partdata[i].agrav[k] += aperiodic[k];
        partdata[i].gpot += potperiodic;
      }
    }

    // Compute all star forces
    sph->ComputeStarGravForces(nbody->Nnbody,nbody->nbodydata,partdata[i]);

    for (k=0; k<ndim; k++) partdata[i].a[k] += partdata[i].agrav[k];
    partdata[i].active = false;

  }
  //-----------------------------------------------------------------------------------------------

  delete[] neibdata;
  delete[] neiblist;

  return;
}



//=================================================================================================
//  SphBruteForceSearch::UpdateAllSphPeriodicHydroForces
/// Routine for computing SPH properties (smoothing lengths, densities and
/// forces) for all active SPH particle using neighbour lists generated
/// using brute force (i.e. direct summation).
//=================================================================================================
template <int ndim, template<int> class ParticleType>
void SphBruteForceSearch<ndim,ParticleType>::UpdateAllSphPeriodicHydroForces
 (int Nhydro,                          ///< [in] No. of SPH particles
  int Ntot,                            ///< [in] No. of SPH + ghost particles
  SphParticle<ndim> *part_gen,         ///< [inout] Pointer to SPH ptcl array
  Sph<ndim> *sph,                      ///< [in] Pointer to SPH object
  Nbody<ndim> *nbody,                  ///< [in] Pointer to N-body object
  DomainBox<ndim> &simbox)             ///< [in] Simulation box with periodic information
{
  int i,j,k;                           // Particle and dimension counters
  int Nneib;                           // No. of neighbours
  int *neiblist;                       // List of neighbour ids
  FLOAT draux[ndim];                   // Relative distance vector
  FLOAT drsqd;                         // Distance squared
  FLOAT dr_corr[ndim];                 // Periodic correction vector
  FLOAT hrangesqdi;                    // Gather kernel extent (squared)
  FLOAT hrangesqdj;                    // Scatter kernel extent (squared)
  FLOAT rp[ndim];                      // Position of current particle
  FLOAT *dr;                           // Array of neib. position vectors
  FLOAT *drmag;                        // Array of neib. distances
  FLOAT *invdrmag;                     // Array of neib. inverse distances
  ParticleType<ndim>* neibdata;        // Local copy of neighbouring particles
  ParticleType<ndim>* partdata = static_cast<ParticleType<ndim>* > (part_gen);
  const int offset_imported = sph->Nghost;

  debug2("[SphBruteForceSearch::UpdateAllSphPeriodicHydroForces]");

  // Construct the mirrors
  MirrorNeighbourFinder<ndim> NgbFinder(simbox) ;
  NgbFinder.SetTargetCell(DomainCellProxy<ndim>(simbox)) ;

  // Allocate memory for storing neighbour ids and position data
  neibdata = new ParticleType<ndim>[Ntot];
  neiblist = new int[Ntot];
  dr       = new FLOAT[ndim*Ntot];
  drmag    = new FLOAT[Ntot];
  invdrmag = new FLOAT[Ntot];


  // Compute forces of real and imported particles
  //-----------------------------------------------------------------------------------------------
  for (int ipart=0; ipart<Nhydro+sph->NImportedParticles; ipart++) {

    if (ipart < Nhydro) i = ipart;
    else i = ipart + offset_imported;

    // Skip over inactive particles
    if (!partdata[i].active || partdata[i].itype == dead) continue;

    // Zero all arrays to be updated
    for (k=0; k<ndim; k++) partdata[i].a[k] = (FLOAT) 0.0;
    for (k=0; k<ndim; k++) partdata[i].agrav[k] = (FLOAT) 0.0;
    partdata[i].gpot      = (FLOAT) 0.0;
    partdata[i].dudt      = (FLOAT) 0.0;
    partdata[i].levelneib = 0;

    for (k=0; k<ndim; k++) rp[k] = partdata[i].r[k];
    hrangesqdi = partdata[i].hrangesqd;

    NgbFinder.SetTargetCell(ParticleCellProxy<ndim>(partdata[i])) ;

    // Determine interaction list (to ensure we don't compute pair-wise forces twice).
    // Also make sure that only the closest periodic replica is considered.
    Nneib = 0;
    for (j=0; j<Nhydro; j++) {
      if (i != j && partdata[j].itype != dead) {
        neiblist[Nneib] = j;

        int NumMirrors = NgbFinder.ConstructGhostsScatterGather(partdata[j], neibdata + Nneib) ;

        for (int n=0; n < NumMirrors; n++){
          for (k=0; k < ndim; k++) draux[k] = neibdata[Nneib].r[k] - partdata[i].r[k] ;
          drsqd = DotProduct(draux, draux, ndim);
          drmag[Nneib] = sqrt(drsqd + small_number);

          invdrmag[Nneib] = (FLOAT) 1.0/drmag[Nneib];
          for (k=0; k<ndim; k++) dr[Nneib*ndim + k] = draux[k]*invdrmag[Nneib];
          Nneib++;
        }
      }
    }
    // Compute all SPH hydro forces
    sph->ComputeSphHydroForces(i, Nneib, neiblist, drmag, invdrmag, dr, partdata[i], partdata);

    // Compute all star forces
    if (ipart < Nhydro) {
      sph->ComputeStarGravForces(nbody->Nnbody, nbody->nbodydata, partdata[i]);
    }

    partdata[i].active = false;

  }
  //-----------------------------------------------------------------------------------------------


  // Free all allocated memory
  delete[] invdrmag;
  delete[] drmag;
  delete[] dr;
  delete[] neiblist;
  delete[] neibdata;

  return;
}



//=================================================================================================
//  SphBruteForceSearch::UpdateAllSphPeriodicGravForces
/// Update all SPH gravity forces, including periodic Ewald corrections.
//=================================================================================================
template <int ndim, template<int> class ParticleType>
void SphBruteForceSearch<ndim,ParticleType>::UpdateAllSphPeriodicGravForces
 (int Nhydro,                          ///< [in] No. of SPH particles
  int Ntot,                            ///< [in] Total no. of particles
  SphParticle<ndim> *part_gen,         ///< [in] Pointer to SPH particle array
  Sph<ndim> *sph,                      ///< [inout] Pointer to SPH object
  Nbody<ndim> *nbody,                  ///< [in] Pointer to N-body object
  DomainBox<ndim> &simbox,             ///< [in] Simulation box with periodic information
  Ewald<ndim> *ewald)                  ///< [in] Pointer to Ewald object for periodic forces
{
  int i,j,k;                           // Particle and dimension counters
  int Nneib;                           // No. of neighbours
  int *neiblist;                       // List of neighbour ids
  FLOAT potperiodic;                   // Periodic potential correction
  FLOAT aperiodic[ndim];               // Periodic acceleration correction
  FLOAT dr[ndim];                      // Relative displacement vector
  FLOAT dr_corr[ndim];                 // Periodic correction vector
  ParticleType<ndim>* neibdata;        // Local copy of neighbouring particles
  ParticleType<ndim>* partdata = static_cast<ParticleType<ndim>* > (part_gen);

  debug2("[SphBruteForceSearch::UpdateAllSphPeriodicGravForces]");

  // Construct the mirrors
  MirrorNeighbourFinder<ndim> NgbFinder(simbox) ;
  NgbFinder.SetTargetCell(DomainCellProxy<ndim>(simbox)) ;

  FLOAT r_ghost[ndim] ;
  int   sign[ndim] ;

  assert(NgbFinder.MaxNumMirrors() == 1) ;

  // Allocate memory for storing neighbour ids and position data
  neiblist = new int[Ntot];
  neibdata = new ParticleType<ndim>[Ntot];

  // Compute smoothing lengths of all SPH particles
  //-----------------------------------------------------------------------------------------------
  for (i=0; i<Nhydro; i++) {

    // Skip over inactive particles
    if (!partdata[i].active || partdata[i].itype == dead) continue;

    // Zero all arrays to be updated
    for (k=0; k<ndim; k++) partdata[i].a[k] = (FLOAT) 0.0;
    for (k=0; k<ndim; k++) partdata[i].agrav[k] = (FLOAT) 0.0;
    partdata[i].gpot = (FLOAT) 0.0;
    partdata[i].dudt = (FLOAT) 0.0;
    partdata[i].levelneib = 0;

    // Add self-contribution to gravitational potential
    partdata[i].gpot += partdata[i].m*partdata[i].invh*kernp->wpot(0.0);

    // Determine interaction list (to ensure we don't compute pair-wise forces twice).
    // Also make sure that only the closest periodic replica is considered.
    Nneib = 0;
    for (j=0; j<Nhydro; j++) {
      neibdata[j] = partdata[j];
      if (i != j && partdata[j].itype != dead) {
        neiblist[Nneib++] = j;
        for (k=0; k<ndim; k++) dr[k] = neibdata[j].r[k] - partdata[i].r[k];
        NgbFinder.NearestPeriodicVector(dr);
        for (k=0; k<ndim; k++) neibdata[j].r[k] = partdata[i].r[k] + dr[k];
      };
    }

    // Compute forces between SPH neighbours (hydro and gravity)
    sph->ComputeSphGravForces(i,Nneib,neiblist,partdata[i],neibdata);

    // Now add the periodic correction force
    if (simbox.PeriodicGravity) {
      for (j=0; j<Nneib; j++) {
        for (k=0; k<ndim; k++) dr[k] = neibdata[j].r[k] - partdata[i].r[k];
        ewald->CalculatePeriodicCorrection(neibdata[j].m,dr,aperiodic,potperiodic);
        for (k=0; k<ndim; k++) partdata[i].agrav[k] += aperiodic[k];
        partdata[i].gpot += potperiodic;
      }
    }
    // Compute all star forces
    sph->ComputeStarGravForces(nbody->Nnbody,nbody->nbodydata,partdata[i]);

    for (k=0; k<ndim; k++) partdata[i].a[k] += partdata[i].agrav[k];
    partdata[i].active = false;

  }
  //-----------------------------------------------------------------------------------------------

  delete[] neibdata;
  delete[] neiblist;

  return;
}



template class SphBruteForceSearch<1,SphParticle>;
template class SphBruteForceSearch<2,SphParticle>;
template class SphBruteForceSearch<3,SphParticle>;
template class SphBruteForceSearch<1,GradhSphParticle>;
template class SphBruteForceSearch<2,GradhSphParticle>;
template class SphBruteForceSearch<3,GradhSphParticle>;
template class SphBruteForceSearch<1,SM2012SphParticle>;
template class SphBruteForceSearch<2,SM2012SphParticle>;
template class SphBruteForceSearch<3,SM2012SphParticle>;
