// ============================================================================
// BruteForceSearch.cpp
// ============================================================================



#include <iostream>
#include "SphNeighbourSearch.h"
#include "Sph.h"
#include "Parameters.h"
#include "SphParticle.h"
#include "Debug.h"
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
void BruteForceSearch::UpdateAllSphProperties(Sph *sph, Parameters &simparams)
{
  int okflag;
  int *neiblist;
  int Nneib;

  debug2("[BruteForceSearch::UpdateAllSphProperties]");

  // The potential number of neighbours is given by ALL the particles
  Nneib = sph->Ntot;

  // Compute smoothing lengths of all SPH particles
  for (int i=0; i<sph->Nsph; i++) 
    okflag = sph->ComputeH(i,Nneib,sph->sphdata,simparams);

  // Compute all other SPH properties
  for (int i=0; i<sph->Nsph; i++)
    sph->ComputeSphProperties(i,Nneib,sph->sphdata,simparams);

  return;
}



// ============================================================================
// BruteForceSearch::UpdateAllSphForces
// ============================================================================
void BruteForceSearch::UpdateAllSphForces(Sph *sph, Parameters &params)
{
  int *neiblist;
  int Nneib;

  debug2("[BruteForceSearch::UpdateAllSphForces]");

  // Allocate array to store local copy of potential neighbour ids
  Nneib = sph->Ntot;

  // Compute SPH hydro forces for all particles
  if (params.intparams["hydro_forces"] == 1)
    for (int i=0; i<sph->Nsph; i++) 
      sph->ComputeHydroForces(i,Nneib,sph->sphdata,params);

  // Compute SPH hydro forces for all particles
  if (params.intparams["self_gravity"] == 1)
    for (int i=0; i<sph->Nsph; i++)
      sph->ComputeGravForces(i,Nneib,sph->sphdata);


  return;
}

