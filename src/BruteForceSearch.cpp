// ============================================================================
// BruteForceSearch.cpp
// Contains all routines for generating SPH neighbour lists using 
// brute-force (i.e. direct summation over all particles).
// ============================================================================



#include <iostream>
#include <math.h>
#include "SphNeighbourSearch.h"
#include "Sph.h"
#include "Parameters.h"
#include "SphParticle.h"
#include "Debug.h"
#include "InlineFuncs.h"
using namespace std;



// ============================================================================
// BruteForceSearch::BruteForceSearch
// ============================================================================
BruteForceSearch::BruteForceSearch(int ndimaux)
{
#if !defined(FIXED_DIMENSIONS)
  ndim = ndimaux;
#endif
}



// ============================================================================
// BruteForceSearch::~BruteForceSearch
// ============================================================================
BruteForceSearch::~BruteForceSearch()
{
}



// ============================================================================
// BruteForceSearch::UpdateTree
// ============================================================================
void BruteForceSearch::UpdateTree(Sph *sph, Parameters &simparams)
{
  return;
}



// ============================================================================
// BruteForceSearch::UpdateAllSphProperties
// Routine for computing SPH properties (smoothing lengths, densities and 
// forces) for all active SPH particle using neighbour lists generated 
// using brute force (i.e. direct summation).
// ============================================================================
void BruteForceSearch::UpdateAllSphProperties(Sph *sph)
{
  int i,j,k;                            // Particle and dimension counters
  int okflag;                           // Flag valid smoothing length
  int Nneib;                            // No. of neighbours
  int Nfar;                             // No. of 'far' neighbours
  int *neiblist;                        // List of neighbour ids
  FLOAT draux[ndimmax];                 // Relative distance vector
  FLOAT drsqd;                          // Distance squared
  FLOAT rp[ndimmax];                    // Position of current particle
  FLOAT *dr;                            // Array of neib. position vectors
  FLOAT *drmag;                         // Array of neib. distances
  FLOAT *invdrmag;                      // Array of neib. inverse distances

  debug2("[BruteForceSearch::UpdateAllSphProperties]");

  // The potential number of neighbours is given by ALL the particles
  Nneib = sph->Ntot;

  // Allocate memory for storing neighbour ids and position data
  neiblist = new int[sph->Ntot];
  dr = new FLOAT[ndim*sph->Ntot];
  drmag = new FLOAT[sph->Ntot];
  invdrmag = new FLOAT[sph->Ntot];

  // Compute smoothing lengths of all SPH particles
  // --------------------------------------------------------------------------
  for (i=0; i<sph->Ntot; i++) {

    // Compute distances and the reciprical between the current particle 
    // and all neighbours here
    for (k=0; k<ndim; k++) rp[k] = sph->sphdata[i].r[k];
    for (j=0; j<Nneib; j++) { 
      neiblist[j] = j;
      for (k=0; k<ndim; k++) draux[k] = sph->sphdata[j].r[k] - rp[k];
      drsqd = DotProduct(draux,draux,ndim);
      drmag[j] = sqrt(drsqd);
      invdrmag[j] = 1.0/(drmag[j] + small_number);
      for (k=0; k<ndim; k++) dr[j*ndim + k] = draux[k]*invdrmag[j];
    }

    // Compute smoothing length for current particle
    okflag = sph->ComputeH(i,sph->sphdata[i],Nneib,Nneib,neiblist,
			   sph->sphdata,drmag,invdrmag,dr);
    
    // Compute all SPH hydro forces
    sph->ComputeHydroForces(i,sph->sphdata[i],Nneib,Nneib,neiblist,
			    sph->sphdata,drmag,invdrmag,dr);

  }
  // --------------------------------------------------------------------------

  delete[] invdrmag;
  delete[] drmag;
  delete[] dr;
  delete[] neiblist;

  return;
}



// ============================================================================
// BruteForceSearch::UpdateAllSphGravityProperties
// Routine for computing SPH properties for all active SPH particle using 
// neighbour lists generated using brute force (i.e. direct summation).
// ============================================================================
void BruteForceSearch::UpdateAllSphGravityProperties(Sph *sph)
{
  int i,j,k;                            // Particle and dimension counters
  int okflag;                           // Flag valid smoothing length
  int Nneib;                            // No. of neighbours
  int Nfar;                             // No. of 'far' neighbours
  int *neiblist;                        // List of neighbour ids
  FLOAT draux[ndimmax];                 // Relative distance vector
  FLOAT drsqd;                          // Distance squared
  FLOAT rp[ndimmax];                    // Position of current particle
  FLOAT *dr;                            // Array of neib. position vectors
  FLOAT *drmag;                         // Array of neib. distances
  FLOAT *invdrmag;                      // Array of neib. inverse distances

  debug2("[BruteForceSearch::UpdateAllSphGravityProperties]");

  // The potential number of neighbours is given by ALL the particles
  Nneib = sph->Ntot;
  neiblist = new int[sph->Ntot];
  drmag = new FLOAT[sph->Ntot];
  invdrmag = new FLOAT[sph->Ntot];
  dr = new FLOAT[ndim*sph->Ntot];

  // Compute smoothing lengths of all SPH particles
  // --------------------------------------------------------------------------
  for (i=0; i<sph->Ntot; i++) {

    // Compute distances and the reciprical between the current particle 
    // and all neighbours here
    for (k=0; k<ndim; k++) rp[k] = sph->sphdata[i].r[k];
    for (j=0; j<Nneib; j++) { 
      neiblist[j] = j;
      for (k=0; k<ndim; k++) draux[k] = sph->sphdata[j].r[k] - rp[k];
      drsqd = DotProduct(draux,draux,ndim);
      drmag[j] = sqrt(drsqd);
      invdrmag[j] = 1.0/(drmag[j] + small_number);
      for (k=0; k<ndim; k++) dr[j*ndim + k] = draux[k]*invdrmag[j];
    }

    // Compute smoothing length for current particle
    okflag = sph->ComputeH(i,sph->sphdata[i],Nneib,Nneib,neiblist,
			   sph->sphdata,drmag,invdrmag,dr);
    
    // Compute all SPH forces
    sph->ComputeHydroForces(i,sph->sphdata[i],Nneib,Nneib,neiblist,
			    sph->sphdata,drmag,invdrmag,dr);

  }
  // --------------------------------------------------------------------------

  delete[] dr;
  delete[] invdrmag;
  delete[] drmag;
  delete[] neiblist;

  return;
}
