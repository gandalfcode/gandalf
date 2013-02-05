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

  // Allocate array to store local copy of potential neighbour ids
  Nneib = sph->Ntot;
  neiblist = new int[sph->Ntot];
  for (int i=0; i<sph->Ntot; i++) neiblist[i] = i;

  // Compute smoothing lengths of all SPH particles
  for (int i=0; i<sph->Nsph; i++) 
    okflag = sph->ComputeH(i,Nneib,neiblist,simparams);

  // Compute all other SPH properties
  for (int i=0; i<sph->Nsph; i++)
    sph->ComputeSphProperties(i,Nneib,neiblist,simparams);

  // Free up memory from local array
  delete[] neiblist;

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
  neiblist = new int[sph->Ntot];
  for (int i=0; i<sph->Ntot; i++) neiblist[i] = i;
  cout << "self_gravity : " << params.intparams["self_gravity"] << endl;
  // Compute SPH hydro forces for all particles
  if (params.intparams["hydro_forces"] == 1)
    for (int i=0; i<sph->Nsph; i++) 
      sph->ComputeHydroForces(i,Nneib,neiblist,params);

  // Compute SPH hydro forces for all particles
  if (params.intparams["self_gravity"] == 1)
    for (int i=0; i<sph->Nsph; i++)
      sph->ComputeGravForces(i,Nneib,neiblist);

  // Free up memory from local array
  delete[] neiblist;

  return;
}

