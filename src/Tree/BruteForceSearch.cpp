//=================================================================================================
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
//=================================================================================================



#include <assert.h>
#include <iostream>
#include <math.h>
#include <numeric>
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
//  BruteForceSearch::BruteForceSearch
/// BruteForceSearch class constructor
//=================================================================================================
template <int ndim, template<int> class ParticleType>
BruteForceSearch<ndim,ParticleType>::BruteForceSearch
 (FLOAT kernrangeaux,
  DomainBox<ndim> *boxaux,
  SmoothingKernel<ndim> *kernaux,
  CodeTiming *timingaux):
  SphNeighbourSearch<ndim>(kernrangeaux,boxaux,kernaux,timingaux)
{
}



//=================================================================================================
//  BruteForceSearch::~BruteForceSearch
/// BruteForceSearch class destructor
//=================================================================================================
template <int ndim, template<int> class ParticleType>
BruteForceSearch<ndim,ParticleType>::~BruteForceSearch()
{
}



//=================================================================================================
//  BruteForceSearch::BuildTree
/// For Brute Force neighbour searching, there is no tree to construct but we delete any
/// dead SPH particles here to be consistent with the tree neighbour search.
//=================================================================================================
template <int ndim, template<int> class ParticleType>
void BruteForceSearch<ndim,ParticleType>::BuildTree
 (bool rebuild_tree,                   ///< Flag to rebuild tree
  int n,                               ///< Integer time
  int ntreebuildstep,                  ///< Tree build frequency
  int ntreestockstep,                  ///< Tree stocking frequency
  int Npart,                           ///< No. of particles
  int Npartmax,                        ///< Max. no. of particles
  SphParticle<ndim> *sph_gen,          ///< Generic pointer to SPH data array
  Sph<ndim> *sph,                      ///< Particle data array
  FLOAT timestep)                      ///< Smallest physical timestep
{
  sph->DeleteDeadParticles();
  return;
}



//=================================================================================================
//  BruteForceSearch::GetGatherNeighbourList
/// ..
//=================================================================================================
template <int ndim, template<int> class ParticleType>
int BruteForceSearch<ndim,ParticleType>::GetGatherNeighbourList
 (FLOAT rp[ndim],                      ///< [in] Position
  FLOAT rsearch,                       ///< [in] Search radius
  SphParticle<ndim> *sph_gen,          ///< [in] Pointer to SPH data array
  int Nhydro,                            ///< [in] No. of SPH particles
  int Nneibmax,                        ///< [in] Max. size of neighbour list
  int *neiblist)                       ///< [out] List of neighbour ids
{
  int i,k;                             // Particle and dimension counters
  int Nneib = 0;                       // No. of (non-dead) neighbours
  FLOAT dr[ndim];                      // Relative distance vector
  FLOAT drsqd;                         // Distance squared
  ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sph_gen);

  debug2("[BruteForceSearch::GetGatherNeighbourList]");

  // Compute smoothing lengths of all SPH particles
  //-----------------------------------------------------------------------------------------------
  for (i=0; i<Nhydro; i++) {

    // Skip over inactive particles
    if (!sphdata[i].active || sphdata[i].itype == dead) continue;

    for (k=0; k<ndim; k++) dr[k] = sphdata[i].r[k] - rp[k];
    drsqd = DotProduct(dr,dr,ndim);

    if (drsqd <= rsearch*rsearch && Nneib < Nneibmax) neiblist[Nneib++] = i;
    else if (drsqd <= rsearch*rsearch && Nneib == Nneibmax) return -1;

  }
  //-----------------------------------------------------------------------------------------------

  return Nneib;
}



//=================================================================================================
//  BruteForceSearch::UpdateAllSphProperties
/// Routine for computing SPH properties (smoothing lengths, densities and
/// forces) for all active SPH particle using neighbour lists generated
/// using brute force (i.e. direct summation).
//=================================================================================================
template <int ndim, template<int> class ParticleType>
void BruteForceSearch<ndim,ParticleType>::UpdateAllSphProperties
 (int Nhydro,                            ///< [in] No. of SPH particles
  int Ntot,                            ///< [in] No. of SPH + ghost particles
  SphParticle<ndim> *sph_gen,          ///< [inout] Pointer to SPH ptcl array
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
  //===============================================================================================
#pragma omp parallel default(none) private(dr,drsqd,i,j,jj,k,okflag,rp)	\
  shared(gpot,m,mu,nbody,neiblist,Nneib,Nhydro,Ntot,sph,sphdata)
  {
    drsqd = new FLOAT[Ntot];

    // Compute smoothing lengths of all SPH particles
    //---------------------------------------------------------------------------------------------
#pragma omp for
    for (i=0; i<Nhydro; i++) {

      // Skip over inactive particles
      if (!sphdata[i].active || sphdata[i].itype == dead) continue;

      for (k=0; k<ndim; k++) rp[k] = sphdata[i].r[k];

      // Compute distances and the reciprical between the current particle
      // and all neighbours here
      //-------------------------------------------------------------------------------------------
      for (jj=0; jj<Nneib; jj++) {
        j = neiblist[jj];
        for (k=0; k<ndim; k++) dr[k] = sphdata[j].r[k] - rp[k];
        drsqd[jj] = DotProduct(dr,dr,ndim) + small_number;
      }
      //-------------------------------------------------------------------------------------------

      // Compute all SPH gather properties
      sph->ComputeH(i,Nneib,big_number,m,mu,drsqd,gpot,sphdata[i],nbody);
      //assert(okflag);

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
//  BruteForceSearch::UpdateAllSphHydroForces
/// Routine for computing SPH properties (smoothing lengths, densities and
/// forces) for all active SPH particle using neighbour lists generated
/// using brute force (i.e. direct summation).
//=================================================================================================
template <int ndim, template<int> class ParticleType>
void BruteForceSearch<ndim,ParticleType>::UpdateAllSphHydroForces
 (int Nhydro,                            ///< [in] No. of SPH particles
  int Ntot,                            ///< [in] No. of SPH + ghost particles
  SphParticle<ndim> *sph_gen,          ///< [inout] Pointer to SPH ptcl array
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
  ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sph_gen);

  debug2("[BruteForceSearch::UpdateAllSphHydroForces]");

  // Allocate memory for storing neighbour ids and position data
  neiblist = new int[Ntot];
  dr       = new FLOAT[ndim*Ntot];
  drmag    = new FLOAT[Ntot];
  invdrmag = new FLOAT[Ntot];

  const int offset_imported = sph->Nghost;

  // Compute forces of real and imported particles
  //-----------------------------------------------------------------------------------------------
  for (int ipart=0; ipart<Nhydro+sph->NImportedParticles; ipart++) {

    if (ipart < Nhydro) i = ipart;
    else i = ipart + offset_imported;

    // Skip over inactive particles
    if (!sphdata[i].active || sphdata[i].itype == dead) continue;

    // Zero all arrays to be updated
    for (k=0; k<ndim; k++) sphdata[i].a[k] = (FLOAT) 0.0;
    for (k=0; k<ndim; k++) sphdata[i].agrav[k] = (FLOAT) 0.0;
    sphdata[i].gpot      = (FLOAT) 0.0;
    sphdata[i].gpe       = (FLOAT) 0.0;
    sphdata[i].dudt      = (FLOAT) 0.0;
    sphdata[i].levelneib = 0;

    for (k=0; k<ndim; k++) rp[k] = sphdata[i].r[k];
    hrangesqdi = pow(kernfac*kernp->kernrange*sphdata[i].h,2);
    Nneib = 0;

    // Compute distances and the reciprical between the current particle
    // and all neighbours here
    //---------------------------------------------------------------------------------------------
    for (j=0; j<sph->Nhydro + sph->NPeriodicGhost; j++) {
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
    //---------------------------------------------------------------------------------------------

    // Compute all SPH hydro forces
    sph->ComputeSphHydroForces(i,Nneib,neiblist,drmag,invdrmag,dr,sphdata[i],sphdata);

    // Compute all star forces
    if ( ipart < Nhydro)
      sph->ComputeStarGravForces(nbody->Nnbody,nbody->nbodydata,sphdata[i]);

    sphdata[i].active = false;

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
//  BruteForceSearch::UpdateAllSphForces
/// Update all SPH forces (both hydro and gravity)
//=================================================================================================
template <int ndim, template<int> class ParticleType>
void BruteForceSearch<ndim,ParticleType>::UpdateAllSphForces
 (int Nhydro,                            ///< [in] No. of SPH particles
  int Ntot,                            ///< [in] Total no. of particles
  SphParticle<ndim> *sph_gen,          ///< [inout] Pointer to SPH data
  Sph<ndim> *sph,                      ///< [inout] Pointer to SPH object
  Nbody<ndim> *nbody)                  ///< [in] Pointer to N-body object
{
  int i,j,k;                           // Particle and dimension counters
  int Nneib;                           // No. of neighbours
  int *neiblist;                       // List of neighbour ids
  ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sph_gen);

  debug2("[BruteForceSearch::UpdateAllSphForces]");

  // Allocate memory for storing neighbour ids and position data
  neiblist = new int[Ntot];

  const int offset_imported = sph->Nghost;

  // Compute forces for real and imported particles
  //-----------------------------------------------------------------------------------------------
  for (int ipart=0; ipart<Nhydro+sph->NImportedParticles; ipart++) {

    if (ipart < Nhydro) i = ipart;
    else i = ipart + offset_imported;

    // Skip over inactive particles
    if (!sphdata[i].active || sphdata[i].itype == dead) continue;

    // Zero all arrays to be updated
    for (k=0; k<ndim; k++) sphdata[i].a[k] = (FLOAT) 0.0;
    for (k=0; k<ndim; k++) sphdata[i].agrav[k] = (FLOAT) 0.0;
    sphdata[i].gpot      = (FLOAT) 0.0;
    sphdata[i].gpe       = (FLOAT) 0.0;
    sphdata[i].dudt      = (FLOAT) 0.0;
    sphdata[i].levelneib = 0;

    // Add self-contribution to gravitational potential
    sphdata[i].gpot += sphdata[i].m*sphdata[i].invh*kernp->wpot(0.0);

    // Determine interaction list (to ensure we don't compute pair-wise forces twice)
    Nneib = 0;
    for (j=0; j<sph->Nhydro + sph->NPeriodicGhost; j++) {
      if (i != j && sphdata[j].itype != dead) neiblist[Nneib++] = j;
    }

    // Compute forces between SPH neighbours (hydro and gravity)
    sph->ComputeSphHydroGravForces(i,Nneib,neiblist,sphdata[i],sphdata);

    // Compute all star forces
    if ( ipart < Nhydro)
      sph->ComputeStarGravForces(nbody->Nnbody,nbody->nbodydata,sphdata[i]);

    for (k=0; k<ndim; k++) sphdata[i].a[k] += sphdata[i].agrav[k];
    sphdata[i].active = false;

  }
  //-----------------------------------------------------------------------------------------------

  delete[] neiblist;

  return;
}



//=================================================================================================
//  BruteForceSearch::UpdateAllSphGravForces
/// ...
//=================================================================================================
template <int ndim, template<int> class ParticleType>
void BruteForceSearch<ndim,ParticleType>::UpdateAllSphGravForces
 (int Nhydro,                            ///< [in] ..
  int Ntot,                            ///< [in] ..
  SphParticle<ndim> *sph_gen,          ///< [in] ..
  Sph<ndim> *sph,                      ///< [inout] Pointer to SPH object
  Nbody<ndim> *nbody)                  ///< [in] Pointer to N-body object
{
  int i,j,k;                           // Particle and dimension counters
  int Nneib;                           // No. of neighbours
  int *neiblist;                       // List of neighbour ids
  ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sph_gen);
  const int offset_imported = sph->Nghost;

  debug2("[BruteForceSearch::UpdateAllSphGravForces]");

  // Allocate memory for storing neighbour ids and position data
  neiblist = new int[Ntot];


  // Compute forces for real and imported particles
  //-----------------------------------------------------------------------------------------------
  for (int iparticle=0; iparticle<Nhydro+sph->NImportedParticles; iparticle++) {

    if (iparticle < Nhydro) i = iparticle;
    else i = iparticle + offset_imported;

    // Skip over inactive particles
    if (!sphdata[i].active || sphdata[i].itype == dead) continue;

    // Zero all arrays to be updated
    for (k=0; k<ndim; k++) sphdata[i].a[k] = (FLOAT) 0.0;
    for (k=0; k<ndim; k++) sphdata[i].agrav[k] = (FLOAT) 0.0;
    sphdata[i].gpot      = (FLOAT) 0.0;
    sphdata[i].gpe       = (FLOAT) 0.0;
    sphdata[i].dudt      = (FLOAT) 0.0;
    sphdata[i].levelneib = 0;

    // Add self-contribution to gravitational potential
    sphdata[i].gpot += sphdata[i].m*sphdata[i].invh*kernp->wpot(0.0);

    // Determine interaction list (to ensure we don't compute pair-wise
    // forces twice)
    Nneib = 0;
    for (j=0; j<Nhydro; j++) {
      if (i != j && sphdata[j].itype != dead) neiblist[Nneib++] = j;
    }

    // Compute forces between SPH neighbours (hydro and gravity)
    sph->ComputeSphGravForces(i,Nneib,neiblist,sphdata[i],sphdata);

    // Compute all star forces
    if ( iparticle < Nhydro)
      sph->ComputeStarGravForces(nbody->Nnbody,nbody->nbodydata,sphdata[i]);

    for (k=0; k<ndim; k++) sphdata[i].a[k] += sphdata[i].agrav[k];
    sphdata[i].active = false;

  }
  //-----------------------------------------------------------------------------------------------

  delete[] neiblist;

  return;
}



//=================================================================================================
//  BruteForceSearch::UpdateAllSphPeriodicForces
/// ...
//=================================================================================================
template <int ndim, template<int> class ParticleType>
void BruteForceSearch<ndim,ParticleType>::UpdateAllSphPeriodicForces
 (int Nhydro,                            ///< [in] No. of SPH particles
  int Ntot,                            ///< [in] Total no. of particles
  SphParticle<ndim> *sph_gen,          ///< [in] Pointer to SPH particle array
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
  ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sph_gen);

  debug2("[BruteForceSearch::UpdateAllSphPeriodicForces]");

  // Allocate memory for storing neighbour ids and position data
  neiblist = new int[Ntot];
  neibdata = new ParticleType<ndim>[Ntot];

  // Compute smoothing lengths of all SPH particles
  //-----------------------------------------------------------------------------------------------
  for (i=0; i<Nhydro; i++) {
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
    sphdata[i].gpot += sphdata[i].m*sphdata[i].invh*kernp->wpot(0.0);

    // Determine interaction list (to ensure we don't compute pair-wise forces twice).
    // Also make sure that only the closest periodic replica is considered.
    Nneib = 0;
    for (j=0; j<Nhydro; j++) {
      neibdata[j] = sphdata[j];
      if (i != j && sphdata[j].itype != dead) {
        neiblist[Nneib++] = j;
        for (k=0; k<ndim; k++) dr[k] = neibdata[j].r[k] - sphdata[i].r[k];
        NearestPeriodicVector(simbox,dr,dr_corr);
        for (k=0; k<ndim; k++) neibdata[j].r[k] = sphdata[i].r[k] + dr[k];
      };
    }

    // Compute forces between SPH neighbours (hydro and gravity)
    sph->ComputeSphHydroGravForces(i,Nneib,neiblist,sphdata[i],sphdata);

    // Now add the periodic correction force
    for (j=0; j<Nneib; j++) {
      for (k=0; k<ndim; k++) dr[k] = neibdata[j].r[k] - sphdata[i].r[k];
      ewald->CalculatePeriodicCorrection(neibdata[j].m,dr,aperiodic,potperiodic);
      for (k=0; k<ndim; k++) sphdata[i].agrav[k] += aperiodic[k];
      sphdata[i].gpot += potperiodic;
    }


    // Compute all star forces
    sph->ComputeStarGravForces(nbody->Nnbody,nbody->nbodydata,sphdata[i]);

    for (k=0; k<ndim; k++) sphdata[i].a[k] += sphdata[i].agrav[k];
    sphdata[i].active = false;

  }
  //-----------------------------------------------------------------------------------------------

  delete[] neibdata;
  delete[] neiblist;

  return;
}



//=================================================================================================
//  BruteForceSearch::UpdateAllSphPeriodicGravForces
/// ...
//=================================================================================================
template <int ndim, template<int> class ParticleType>
void BruteForceSearch<ndim,ParticleType>::UpdateAllSphPeriodicGravForces
 (int Nhydro,                            ///< [in] No. of SPH particles
  int Ntot,                            ///< [in] Total no. of particles
  SphParticle<ndim> *sph_gen,          ///< [in] Pointer to SPH particle array
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
  ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sph_gen);

  debug2("[BruteForceSearch::UpdateAllSphPeriodicGravForces]");

  // Allocate memory for storing neighbour ids and position data
  neiblist = new int[Ntot];
  neibdata = new ParticleType<ndim>[Ntot];

  // Compute smoothing lengths of all SPH particles
  //-----------------------------------------------------------------------------------------------
  for (i=0; i<Nhydro; i++) {

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
    sphdata[i].gpot += sphdata[i].m*sphdata[i].invh*kernp->wpot(0.0);

    // Determine interaction list (to ensure we don't compute pair-wise forces twice).
    // Also make sure that only the closest periodic replica is considered.
    Nneib = 0;
    for (j=0; j<Nhydro; j++) {
      neibdata[j] = sphdata[j];
      if (i != j && sphdata[j].itype != dead) {
        neiblist[Nneib++] = j;
        for (k=0; k<ndim; k++) dr[k] = neibdata[j].r[k] - sphdata[i].r[k];
        NearestPeriodicVector(simbox,dr,dr_corr);
        for (k=0; k<ndim; k++) neibdata[j].r[k] = sphdata[i].r[k] + dr[k];
      };
    }

    // Compute forces between SPH neighbours (hydro and gravity)
    sph->ComputeSphGravForces(i,Nneib,neiblist,sphdata[i],neibdata);

    // Now add the periodic correction force
    for (j=0; j<Nneib; j++) {
      for (k=0; k<ndim; k++) dr[k] = neibdata[j].r[k] - sphdata[i].r[k];
      ewald->CalculatePeriodicCorrection(neibdata[j].m,dr,aperiodic,potperiodic);
      for (k=0; k<ndim; k++) sphdata[i].agrav[k] += aperiodic[k];
      sphdata[i].gpot += potperiodic;
    }

    // Compute all star forces
    sph->ComputeStarGravForces(nbody->Nnbody,nbody->nbodydata,sphdata[i]);

    for (k=0; k<ndim; k++) sphdata[i].a[k] += sphdata[i].agrav[k];
    sphdata[i].active = false;

  }
  //-----------------------------------------------------------------------------------------------

  delete[] neibdata;
  delete[] neiblist;

  return;
}



//=================================================================================================
//  BruteForceSearch::UpdateAllSphDerivatives
/// Compute all SPH derivatives required for 2nd-order Riemann solver in
/// Godunov SPH method.
//==================================================================================================
template <int ndim, template<int> class ParticleType>
void BruteForceSearch<ndim,ParticleType>::UpdateAllSphDerivatives
 (int Nhydro,                          ///< [in] ..
  int Ntot,                            ///< [in] ..
  SphParticle<ndim> *sph_gen,          ///< ..
  Sph<ndim> *sph)                      ///< [inout] Pointer to SPH object
{
  int i,j,k;                           // Particle and dimension counters
  int Nneib;                           // No. of neighbours
  int *neiblist;                       // List of neighbour ids
  FLOAT draux[ndim];                   // Relative distance vector
  FLOAT drsqd;                         // Distance squared
  FLOAT hrangesqd;                     // Kernel extent (squared)
  FLOAT rp[ndim];                      // Position of current particle
  FLOAT *dr;                           // Array of neib. position vectors
  FLOAT *drmag;                        // Array of neib. distances
  FLOAT *invdrmag;                     // Array of neib. inverse distances
  struct SphParticle<ndim> *neibpart;  // Local copies of neib. particles
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
  //-----------------------------------------------------------------------------------------------
  for (i=0; i<Nhydro; i++) {
    for (k=0; k<ndim; k++) rp[k] = sphdata[i].r[k];
    hrangesqd = pow(kernp->kernrange*sphdata[i].h,2);
    Nneib = 0;

    // Compute distances and the reciprical between the current particle and all neighbours here
    //---------------------------------------------------------------------------------------------
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
    //---------------------------------------------------------------------------------------------

    // Compute all SPH hydro forces
    sph->ComputeSphDerivatives(i, Nneib, neiblist, drmag, invdrmag, dr, sphdata[i], neibpart);

  }
  //-----------------------------------------------------------------------------------------------

  delete[] neibpart;
  delete[] invdrmag;
  delete[] drmag;
  delete[] dr;
  delete[] neiblist;

  return;
}



//=================================================================================================
//  BruteForceSearch::UpdateAllSphDudt
/// Compute the compressional heating rate (dudt) for all active particles.
//=================================================================================================
template <int ndim, template<int> class ParticleType>
void BruteForceSearch<ndim,ParticleType>::UpdateAllSphDudt
 (int Nhydro,                          ///< [in] ..
  int Ntot,                            ///< [in] ..
  SphParticle<ndim> *sph_gen,          ///< [inout] ..
  Sph<ndim> *sph)                      ///< [inout] Pointer to SPH object
{
  int i,j,k;                           // Particle and dimension counters
  int Nneib;                           // No. of neighbours
  int *neiblist;                       // List of neighbour ids
  FLOAT draux[ndim];                   // Relative distance vector
  FLOAT drsqd;                         // Distance squared
  FLOAT hrangesqdi;                    // Gather kernel range (squared)
  FLOAT hrangesqdj;                    // Scatter kernel range (squared)
  FLOAT rp[ndim];                      // Position of current particle
  FLOAT *dr;                           // Array of neib. position vectors
  FLOAT *drmag;                        // Array of neib. distances
  FLOAT *invdrmag;                     // Array of neib. inverse distances
  struct SphParticle<ndim> *neibpart;  // Local copies of neighbour particles
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
  //-----------------------------------------------------------------------------------------------
  for (i=0; i<Nhydro; i++) {
    for (k=0; k<ndim; k++) rp[k] = sphdata[i].r[k];
    hrangesqdi = pow(kernfac*kernp->kernrange*sphdata[i].h,2);
    Nneib = 0;

    // Compute distances and the reciprical between the current particle
    // and all neighbours here
    //---------------------------------------------------------------------------------------------
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
    //---------------------------------------------------------------------------------------------

    // Compute all SPH hydro forces
    sph->ComputeSphNeibDudt(i,Nneib,neiblist,drmag,invdrmag,dr,sphdata[i],neibpart);

  }
  //-----------------------------------------------------------------------------------------------


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



//=================================================================================================
//  BruteForceSearch::UpdateAllStarGasForces
/// ..
//=================================================================================================
template <int ndim, template<int> class ParticleType>
void BruteForceSearch<ndim,ParticleType>::UpdateAllStarGasForces
 (int Nhydro,                          ///< [in] ..
  int Ntot,                            ///< [in] No. of SPH particles
  SphParticle<ndim> *sph_gen,          ///< [in] ..
  Sph<ndim> *sph,                      ///< [inout] Pointer to SPH ptcl array
  Nbody<ndim> *nbody)                  ///< [inout] Pointer to N-body object
{
  int i;                               // Particle and dimension counters
  int Nneib = 0;                       // No. of (non-dead) neighbours
  int *dummy = 0;                      // Dummy var to satisfy function argument
  int *neiblist;                       // Array of (all) particle ids
  ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sph_gen);

  debug2("[BruteForceSearch::UpdateAllSphForces]");

  // Allocate memory for storing neighbour ids and position data
  neiblist = new int[Nhydro];
  for (i=0; i<Nhydro; i++) {
    if (sphdata[i].itype != dead) neiblist[Nneib++] = i;
  }

  // Compute smoothing lengths of all SPH particles
  //-----------------------------------------------------------------------------------------------
  for (i=0; i<nbody->Nnbody; i++) {

    // Skip over inactive particles
    if (!nbody->nbodydata[i]->active) continue;

    // Compute forces between SPH neighbours (hydro and gravity)
    nbody->CalculateDirectSPHForces(nbody->nbodydata[i],Nneib,0,neiblist,dummy,sph);

  }
  //-----------------------------------------------------------------------------------------------

  delete[] neiblist;

  return;
}



//=================================================================================================
//  BruteForceSearch::SearchBoundaryGhostParticles
/// Search domain to create any required ghost particles near any boundaries.
/// Currently only searches to create periodic or mirror ghost particles.
//=================================================================================================
template <int ndim, template <int> class ParticleType>
void BruteForceSearch<ndim,ParticleType>::SearchBoundaryGhostParticles
(FLOAT tghost,                      ///< Ghost particle 'lifetime'
 DomainBox<ndim> simbox,            ///< Simulation box structure
 Sph<ndim> *sph)                    ///< Sph object pointer
{
  int i;                            // Particle counter

  // Set all relevant particle counters
  sph->Nghost         = 0;
  sph->NPeriodicGhost = 0;
  sph->Nghostmax      = sph->Nhydromax - sph->Nhydro;
  sph->Ntot           = sph->Nhydro;


  // If all boundaries are open, immediately return to main loop
  if (simbox.boundary_lhs[0] == openBoundary && simbox.boundary_rhs[0] == openBoundary &&
      simbox.boundary_lhs[1] == openBoundary && simbox.boundary_rhs[1] == openBoundary &&
      simbox.boundary_lhs[2] == openBoundary && simbox.boundary_rhs[2] == openBoundary)
    return;


  debug2("[BruteForceSearch::SearchBoundaryGhostParticles]");


  // Create ghost particles in x-dimension
  //-----------------------------------------------------------------------------------------------
  if ((simbox.boundary_lhs[0] == openBoundary && simbox.boundary_rhs[0] == openBoundary) == 0) {

    for (i=0; i<sph->Ntot; i++) {
      sph->CheckXBoundaryGhostParticle(i,tghost,simbox);
    }

    sph->Ntot = sph->Nhydro + sph->Nghost;
  }


  // Create ghost particles in y-dimension
  //-----------------------------------------------------------------------------------------------
  if (ndim >= 2 && (simbox.boundary_lhs[1] == openBoundary &&
                    simbox.boundary_rhs[1] == openBoundary) == 0) {

    for (i=0; i<sph->Ntot; i++) {
      sph->CheckYBoundaryGhostParticle(i,tghost,simbox);
    }

    sph->Ntot = sph->Nhydro + sph->Nghost;
  }


  // Create ghost particles in z-dimension
  //-----------------------------------------------------------------------------------------------
  if (ndim == 3 && (simbox.boundary_lhs[2] == openBoundary &&
                    simbox.boundary_rhs[2] == openBoundary) == 0) {

    for (i=0; i<sph->Ntot; i++) {
      sph->CheckZBoundaryGhostParticle(i,tghost,simbox);
    }

    sph->Ntot = sph->Nhydro + sph->Nghost;
  }


  // Quit here if we've run out of memory for ghosts
  if (sph->Ntot > sph->Nhydromax) {
    string message="Not enough memory for ghost particles";
    ExceptionHandler::getIstance().raise(message);
  }

  sph->NPeriodicGhost = sph->Nghost;

  return;
}



#if defined MPI_PARALLEL
//=================================================================================================
//  BruteForceSearch::UpdateGravityExportList
/// ..
//=================================================================================================
template <int ndim, template<int> class ParticleType>
void BruteForceSearch<ndim,ParticleType>::UpdateGravityExportList
 (int rank,                            ///< [in] ..
  int Nhydro,                          ///< [in] No. of SPH particles
  int Ntot,                            ///< [in] No. of SPH + ghost particles
  SphParticle<ndim> *sph_gen,          ///< [inout] Pointer to SPH ptcl array
  Sph<ndim> *sph,                      ///< [in] Pointer to SPH object
  Nbody<ndim> *nbody)                  ///< [in] Pointer to N-body object
{
  return;
}



//=================================================================================================
//  BruteForceSearch::UpdateHydroExportList
/// ..
//=================================================================================================
template <int ndim, template<int> class ParticleType>
void BruteForceSearch<ndim,ParticleType>::UpdateHydroExportList
 (int rank,                            ///< [in] ..
  int Nhydro,                          ///< [in] No. of SPH particles
  int Ntot,                            ///< [in] No. of SPH + ghost particles
  SphParticle<ndim> *sph_gen,          ///< [inout] Pointer to SPH ptcl array
  Sph<ndim> *sph,                      ///< [in] Pointer to SPH object
  Nbody<ndim> *nbody)                  ///< [in] Pointer to N-body object
{
  return;
}



//=================================================================================================
//  BruteForceSearch::SearchMpiGhostParticles
/// Compute on behalf of the MpiControl class the ghost particles we need to export to other nodes.
//=================================================================================================
template <int ndim, template<int> class ParticleType>
int BruteForceSearch<ndim,ParticleType>::SearchMpiGhostParticles
 (const FLOAT tghost,                  ///< [in] Expected ghost life-time
  const Box<ndim> &mpibox,             ///< [in] Bounding box of MPI domain
  Sph<ndim> *sph,                      ///< [in] Pointer to SPH object
  vector<int> &export_list)            ///< [out] List of particle ids
{
  int i;                               // ..
  int k;                               // ..
  int Nexport = 0;                     // No. of MPI ghosts to export
  FLOAT scattermin[ndim];              // ..
  FLOAT scattermax[ndim];              // ..
  const FLOAT grange = ghost_range*kernrange;
  ParticleType<ndim> *sphdata = static_cast<ParticleType<ndim>* > (sph->GetSphParticleArray());

  // Loop over particles and prepare the ones to export
  //-----------------------------------------------------------------------------------------------
  for (i=0; i<sph->Nhydro; i++) {
    ParticleType<ndim>& part = sphdata[i];

    // Construct maximum cell bounding box depending on particle velocities
    for (k=0; k<ndim; k++) {
      scattermin[k] = part.r[k] + min(0.0,part.v[k]*tghost) - grange*part.h;
      scattermax[k] = part.r[k] + max(0.0,part.v[k]*tghost) + grange*part.h;
    }

    // If maximum cell scatter box overlaps MPI domain, open cell
    if (BoxOverlap(ndim,scattermin,scattermax,mpibox.boxmin,mpibox.boxmax)) {
      export_list.push_back(i);
      Nexport++;
    }

  }
  //-----------------------------------------------------------------------------------------------

  return Nexport;
}



//=================================================================================================
//  BruteForceSearch::SearchHydroExportParticles
/// Compute on behalf of the MpiControl class the ghost particles we need to export to other nodes.
//=================================================================================================
template <int ndim, template<int> class ParticleType>
int BruteForceSearch<ndim,ParticleType>::SearchHydroExportParticles
 (const Box<ndim> &mpibox,             ///< [in] Bounding box of MPI domain
  Sph<ndim> *sph,                      ///< [in] Pointer to SPH object
  vector<int> &export_list)            ///< [out] List of particle ids
{
  int i;                               // ..
  int k;                               // ..
  int Nexport = 0;                     // No. of MPI ghosts to export
  FLOAT scattermin[ndim];              // ..
  FLOAT scattermax[ndim];              // ..
  const FLOAT grange = ghost_range*kernrange;
  ParticleType<ndim> *sphdata = static_cast<ParticleType<ndim>* > (sph->GetSphParticleArray());

  // Loop over particles and prepare the ones to export
  //-----------------------------------------------------------------------------------------------
  for (i=0; i<sph->Nhydro; i++) {
    ParticleType<ndim>& part = sphdata[i];

    // Construct maximum cell bounding box depending on particle velocities
    for (k=0; k<ndim; k++) {
      scattermin[k] = part.r[k] - grange*part.h;
      scattermax[k] = part.r[k] + grange*part.h;
    }

    // If maximum cell scatter box overlaps MPI domain, open cell
    if (BoxOverlap(ndim,scattermin,scattermax,mpibox.boxmin,mpibox.boxmax) &&
        part.active) {
      export_list.push_back(i);
      Nexport++;
    }

  }
  //-----------------------------------------------------------------------------------------------

  return Nexport;
}



//=================================================================================================
//  BruteForceSearch::FindMpiTransferParticles
/// Compute on behalf of the MpiControl class the particles that are outside
/// the domain after a load balancing and need to be transferred to other nodes
//=================================================================================================
template <int ndim, template<int> class ParticleType>
void BruteForceSearch<ndim,ParticleType>::FindMpiTransferParticles
 (Sph<ndim>* sph,                          ///< [in] Pointer to sph class
  vector<vector<int> >& id_export_buffers, ///< [inout] Vector that for each node
                                           ///<         gives the list of particles to export
  vector<int>& all_ids_export_buffer,      ///< [inout] Vector containing all the particles that
                                           ///<         will be exported by this processor
  const vector<int>& potential_nodes,      ///< [in] Potential nodes we might send particles to
  MpiNode<ndim>* mpinodes)                 ///< [in] Array of other mpi nodes
{
  int i;
  int inode;
  int node_number;
  ParticleType<ndim> *sphdata = static_cast<ParticleType<ndim>* > (sph->GetSphParticleArray());


  // Loop over particles and prepare the ones to export
  //-----------------------------------------------------------------------------------------------
  for (i=0; i<sph->Nhydro; i++) {
    ParticleType<ndim>& part = sphdata[i];

    // Loop over potential domains and see if we need to transfer
    // this particle to them
    //---------------------------------------------------------------------------------------------
    for (inode=0; inode<potential_nodes.size(); inode++) {
      node_number = potential_nodes[inode];

      // If particle belongs to this domain, add to vector and break from loop
      if (ParticleInBox(part,mpinodes[node_number].domain)) {
        id_export_buffers[node_number].push_back(i);
        all_ids_export_buffer.push_back(i);
        break;
      }

    }
    //---------------------------------------------------------------------------------------------

  }
  //-----------------------------------------------------------------------------------------------

  return;
}



//=================================================================================================
//  BruteForceSearch::FindGhostParticlesToExport
/// Compute on behalf of the MpiControl class the ghost particles we need to export to other nodes.
//=================================================================================================
template <int ndim, template<int> class ParticleType>
void BruteForceSearch<ndim,ParticleType>::FindGhostParticlesToExport
 (Sph<ndim>* sph,                                            ///< [in] Pointer to sph class
  vector<vector<ParticleType<ndim>*> >& ptcl_export_buffers, ///< [inout] Buffers with ptcls to
                                                             ///<         export to each node
  const vector<int>& overlapping_nodes,                      ///< [in] Vector containing which
                                                             ///<      nodes overlap our hbox
  MpiNode<ndim>* mpinodes)                                   ///< [in] Array of other mpi nodes
{
  int i;                                                     // ..
  int inode;                                                 // ..
  int node_number;                                           // ..
  ParticleType<ndim> *sphdata = static_cast<ParticleType<ndim>*> (sph->GetSphParticleArray());

  // Loop over particles and prepare the ones to export
  for (i=0; i<sph->Ntot; i++) {
    ParticleType<ndim>& part = sphdata[i];

    // Loop over potential domains and find particles to export to them
    for (inode=0; inode<overlapping_nodes.size(); inode++) {
      node_number = overlapping_nodes[inode];
      if (ParticleBoxOverlap(part,mpinodes[node_number].hbox)) {
        ptcl_export_buffers[node_number].push_back(&part);
      }
    }
  }

  return;
}



//=================================================================================================
//  BruteForceSearch::FindParticlesToTransfer
/// Compute on behalf of the MpiControl class the particles that are outside
/// the domain after a load balancing and need to be transferred to other nodes
//=================================================================================================
template <int ndim, template<int> class ParticleType>
void BruteForceSearch<ndim,ParticleType>::FindParticlesToTransfer
 (Sph<ndim>* sph,                            ///< [in] Pointer to sph class
  vector<vector<int> >& id_export_buffers,   ///< [inout] List of ids to export for each node
  vector<int>& all_ids_export_buffer,        ///< [inout] List of all ids to export from proc
  const vector<int>& potential_nodes,        ///< [in] Potential nodes we might send particles to
  MpiNode<ndim>* mpinodes)                   ///< [in] Array of other mpi nodes
{
  ParticleType<ndim> *sphdata = static_cast<ParticleType<ndim>* > (sph->GetSphParticleArray());

  // Loop over particles and prepare the ones to export
  for (int i=0; i<sph->Nhydro; i++) {
    ParticleType<ndim>& part = sphdata[i];

    //Loop over potential domains and see if we need to transfer this particle to them
    for (int inode=0; inode<potential_nodes.size(); inode++) {
      int node_number = potential_nodes[inode];
      if (ParticleInBox(part,mpinodes[node_number].domain)) {
        id_export_buffers[node_number].push_back(i);
        all_ids_export_buffer.push_back(i);
        // The particle can belong only to one domain, so we can break from this loop
        break;
      }
    }
  }

  return;
}



//=================================================================================================
//  BruteForceSearch::GetExportInfo
/// Get the array with the information that needs to be exported to the given processor
/// (NB: Nproc is ignored at the moment, as we must always export all ptcls to other processors).
//=================================================================================================
template <int ndim, template<int> class ParticleType>
int BruteForceSearch<ndim,ParticleType>::GetExportInfo
 (int Nproc,                           ///< [in] No of processor we want to send the information to
  Sph<ndim> *sph,                      ///< [in] Pointer to sph object
  vector<char >& ptcl_export_buffer,   ///< [inout] Buffer containing particles to be exported
  MpiNode<ndim>& mpinode,              ///< ..
  int rank,                            ///< ..
  int Nmpi)                            ///< [in]  Array with information for the other mpi nodes
{

  const bool first_proc = (Nproc==0) || (rank==0 && Nproc==1);

  ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sph->GetSphParticleArray() );

  //Find number of active particles
  if (first_proc) {
    ids_active_particles.clear();
    for (int i=0; i<sph->Nhydro; i++) {
      if (sphdata[i].active) {
        ids_active_particles.push_back(i);
      }
    }
  }

  const int Nactive = ids_active_particles.size();
  const int size_export = Nactive*sizeof(ParticleType<ndim>);

  if (first_proc) {
    ptcl_export_buffer.clear();
    ptcl_export_buffer.reserve((Nmpi-1)*size_export);
    ptcl_export_buffer.resize(size_export);
  }
  else {
    ptcl_export_buffer.resize(ptcl_export_buffer.size()+size_export);
  }

//  //Copy positions of active particles inside arrays
//  int j=0;
//  for (int i=0; i< sph->Nhydro; i++) {
//    if (sphdata[i].active) {
//      for (int k=0; k<ndim; k++)
//        ptcl_export_buffer[j].r[k] = sphdata[i].r[k];
//      j++;
//    }
//  }

  //Copy particles to export inside arrays
  int j = (ptcl_export_buffer.size() - size_export)/sizeof(ParticleType<ndim>);
  for (int i=0; i<sph->Nhydro; i++) {
    if (sphdata[i].active) {
      copy(&ptcl_export_buffer[j*sizeof(ParticleType<ndim>)] , &sphdata[i]);
      j++;
    }
  }

  return size_export;
}


//=================================================================================================
//  BruteForceSearch::UnpackExported
/// Unpack the information exported from the other processors, contaning the positions
/// of the particles that were exported
//=================================================================================================
template <int ndim, template<int> class ParticleType>
void BruteForceSearch<ndim,ParticleType>::UnpackExported
 (vector<char > &received_array,
  vector<int> &Nbytes_from_proc,
  Sph<ndim>* sph)
{
  int offset = 0;
  ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sph->GetSphParticleArray() );

  assert(sph->NImportedParticles==0);


  //-----------------------------------------------------------------------------------------------
  for (int Nproc = 0; Nproc<Nbytes_from_proc.size(); Nproc++) {

    int N_received_bytes = Nbytes_from_proc[Nproc];
    int N_received_particles = N_received_bytes/sizeof(ParticleType<ndim>);

    //Ensure there is enough memory
    if (sph->Ntot + N_received_particles > sph->Nhydromax) {
      ExceptionHandler::getIstance().raise("Error while receiving imported particles: not enough memory!");
    }

//    //Copy particle positions inside SPH main arrays
//    for (int i=0; i<N_received_particles; i++) {
//      for (int k=0; k<ndim; k++)
//        sphdata[i+sph->Ntot].r[k] = received_array[i+offset].r[k];
//    }

    //Copy received particles inside SPH main arrays
    for (int i=0; i<N_received_particles; i++) {
      copy( &sphdata[i+sph->Ntot] , &received_array[i*sizeof(ParticleType<ndim>)+offset]);
    }

    //Update the SPH counters
    sph->Ntot += N_received_particles;
    sph->NImportedParticles += N_received_particles;

    //Update the offset
    offset += N_received_bytes;

  }
  //-----------------------------------------------------------------------------------------------

  return;
}



//=================================================================================================
//  BruteForceSearch::GetBackExportInfo
/// Return the data to transmit back to the other processors (particle acceleration etc.)
//=================================================================================================
template <int ndim, template<int> class ParticleType>
void BruteForceSearch<ndim,ParticleType>::GetBackExportInfo
 (vector<char >& send_buffer,          ///< [inout] Buffer containing information to send
  vector<int>& Nbytes_from_proc,       ///< ..
  vector<int>& Nbytes_to_proc,         ///< ..
  Sph<ndim>* sph,                      ///< [in] Pointer to the SPH object
  int rank)                            ///< ..
{
  int removed_particles=0;
  int Nbytes_received_exported = std::accumulate(Nbytes_from_proc.begin(), Nbytes_from_proc.end(), 0);

  send_buffer.resize(Nbytes_received_exported);

  // Loop over the processors, removing particles as we go
  //-----------------------------------------------------------------------------------------------
  for (int Nproc=0 ; Nproc < Nbytes_from_proc.size(); Nproc++ ) {

    const int N_received_particles = Nbytes_from_proc[Nproc]/sizeof(ParticleType<ndim>);

//    //Copy the accelerations and gravitational potential of the particles
//    ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sph->GetSphParticleArray() );
//    int j=0;
//    for (int i=sph->Nhydro - N_received_particles; i<sph->Nhydro; i++) {
//      for (int k=0; k<ndim; k++)
//        send_buffer[removed_particles+j].a[k] = sphdata[i].a[k];
//      send_buffer[removed_particles+j].gpot = sphdata[i].gpot;
//      j++;
//    }

    //Copy the particles inside the send buffer
    ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sph->GetSphParticleArray() );
    int j=0;
    const int start_index = sph->Nhydro + sph->Nghost + removed_particles;
    for (int i=start_index; i<start_index + N_received_particles; i++) {
      copy (&send_buffer[(removed_particles+j)*sizeof(ParticleType<ndim>)],&sphdata[i]);
      j++;
    }

    assert(j==N_received_particles);

    removed_particles += j;

    //Decrease the particle counter
    sph->Ntot -= N_received_particles;
    sph->NImportedParticles -= N_received_particles;

  }
  //-----------------------------------------------------------------------------------------------

  assert(sph->NImportedParticles == 0);
  assert(sph->Ntot == sph->Nhydro + sph->Nghost);
  assert(send_buffer.size() == removed_particles*sizeof(ParticleType<ndim>));

  return;
}



//=================================================================================================
//  BruteForceSearch::UnpackReturnedExportInfo
/// Unpack the data that was returned by the other processors, summing the accelerations to the particles
//=================================================================================================
template <int ndim, template<int> class ParticleType>
void BruteForceSearch<ndim,ParticleType>::UnpackReturnedExportInfo
 (vector<char >& received_information,   ///< ..
  vector<int>& recv_displs,              ///< ..
  Sph<ndim>* sph,                        ///< ..
  int rank)                              ///< ..
{
  ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sph->GetSphParticleArray() );

//  //For each particle, sum up the accelerations returned by other processors
//  for (int i=0; i< ids_active_particles.size(); i++ ) {
//    const int j = ids_active_particles[i];
//
//    for(int Nproc=0; Nproc<recv_displs.size(); Nproc++ ) {
//
//      if (rank==Nproc)
//        continue;
//
//      for (int k=0; k<ndim; k++)
//        sphdata[j].a[k] += received_information[i+recv_displs[Nproc] ].a[k];
//      sphdata[j].gpot += received_information[i+recv_displs[Nproc] ].gpot;
//    }
//
//  }

  // For each particle, sum up some of the quantities returned by other processors
  //-----------------------------------------------------------------------------------------------
  for (int i=0; i<ids_active_particles.size(); i++) {
    const int j = ids_active_particles[i];

    for (int Nproc=0; Nproc<recv_displs.size(); Nproc++) {

      if (rank == Nproc) continue;

      ParticleType<ndim>* received_particle = reinterpret_cast<ParticleType<ndim>*>
        (&received_information[i * sizeof(ParticleType<ndim>) + recv_displs[Nproc] ]);


      assert(sphdata[j].iorig == received_particle->iorig);

      for (int k=0; k<ndim; k++) {
        sphdata[j].a[k] += received_particle->a[k];
        sphdata[j].agrav[k] += received_particle->agrav[k];
      }
      sphdata[j].gpot += received_particle->gpot;
      sphdata[j].gpe += received_particle->gpe;
      sphdata[j].dudt += received_particle->dudt;
      sphdata[j].div_v += received_particle->div_v;
      sphdata[j].levelneib = max(sphdata[j].levelneib, received_particle->levelneib);

    }

  }
  //-----------------------------------------------------------------------------------------------

}
#endif



template class BruteForceSearch<1,SphParticle>;
template class BruteForceSearch<2,SphParticle>;
template class BruteForceSearch<3,SphParticle>;
template class BruteForceSearch<1,GradhSphParticle>;
template class BruteForceSearch<2,GradhSphParticle>;
template class BruteForceSearch<3,GradhSphParticle>;
template class BruteForceSearch<1,SM2012SphParticle>;
template class BruteForceSearch<2,SM2012SphParticle>;
template class BruteForceSearch<3,SM2012SphParticle>;
template class BruteForceSearch<1,GodunovSphParticle>;
template class BruteForceSearch<2,GodunovSphParticle>;
template class BruteForceSearch<3,GodunovSphParticle>;
