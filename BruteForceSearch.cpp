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
  int jj;
  int okflag;
  int *neiblist;
  int Nneib;
  FLOAT draux[ndimmax];
  FLOAT drsqd;
  FLOAT rp[ndimmax];
  FLOAT *dr;
  FLOAT *drmag;
  FLOAT *invdrmag;

  debug2("[BruteForceSearch::UpdateAllSphProperties]");

  // The potential number of neighbours is given by ALL the particles
  Nneib = sph->Ntot;

  drmag = new FLOAT[sph->Ntot];
  invdrmag = new FLOAT[sph->Ntot];
  dr = new FLOAT[ndim*sph->Ntot];

  // Compute smoothing lengths of all SPH particles
  // --------------------------------------------------------------------------
  for (int i=0; i<sph->Ntot; i++) {

    // Compute distances and the reciprical between the current particle 
    // and all neighbours here
    for (int k=0; k<ndim; k++) rp[k] = sph->sphdata[i].r[k];
    for (int jj=0; jj<Nneib; jj++) { 
      for (int k=0; k<ndim; k++) draux[k] = sph->sphdata[jj].r[k] - rp[k];
      drsqd = DotProduct(draux,draux,ndim);
      drmag[jj] = sqrt(drsqd);
      invdrmag[jj] = 1.0/(drmag[jj] + small_number);
      for (int k=0; k<ndim; k++) dr[jj*ndim + k] = draux[k]*invdrmag[jj];
    }

    okflag = sph->ComputeH(i,sph->sphdata[i],Nneib,
			   sph->sphdata,drmag,invdrmag,dr);
    
    sph->ComputeHydroForces(i,sph->sphdata[i],Nneib,
			    sph->sphdata,drmag,invdrmag,dr);

  }
  // --------------------------------------------------------------------------

  delete[] dr;
  delete[] invdrmag;
  delete[] drmag;

  // Compute SPH hydro forces for all particles
  //for (int i=0; i<sph->Nsph; i++)
  //  sph->ComputeGravForces(i,Nneib,sph->sphdata);

  return;
}

