// ============================================================================
// BruteForceSearch.cpp
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
// ============================================================================
void BruteForceSearch::UpdateAllSphProperties(Sph *sph)
{
  int i,j,k;
  int okflag;
  int Nneib;
  int *neiblist;
  FLOAT draux[ndimmax];
  FLOAT drsqd;
  FLOAT rp[ndimmax];
  FLOAT *dr;
  FLOAT *drmag;
  FLOAT *invdrmag;

  debug2("[BruteForceSearch::UpdateAllSphProperties]");

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

    okflag = sph->ComputeH(i,sph->sphdata[i],Nneib,Nneib,neiblist,
			   sph->sphdata,drmag,invdrmag,dr);
    
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



// ============================================================================
// BruteForceSearch::UpdateAllSphProperties
// ============================================================================
void BruteForceSearch::UpdateAllSphGravityProperties(Sph *sph)
{
  return;
}
