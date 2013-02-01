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
BruteForceSearch::BruteForceSearch()
{
}



// ============================================================================
// BruteForceSearch::~BruteForceSearch
// ============================================================================
BruteForceSearch::~BruteForceSearch()
{
}



// ============================================================================
// BruteForceSearch::UpdateAllSphProperties
// ============================================================================
void BruteForceSearch::UpdateAllSphProperties(Sph *sph, Parameters &simparams)
{
  int okflag;
  int *neiblist;
  int Nneib;

  debug2("[BruteForceSearch::UpdateAllSphProperties]\n");

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
void BruteForceSearch::UpdateAllSphForces(Sph *sph, Parameters &simparams)
{
  int *neiblist;
  int Nneib;

  debug2("[BruteForceSearch::UpdateAllSphForces]\n");

  // Allocate array to store local copy of potential neighbour ids
  Nneib = sph->Ntot;
  neiblist = new int[sph->Ntot];
  for (int i=0; i<sph->Ntot; i++) neiblist[i] = i;

  // Compute SPH hydro forces for all particles
  for (int i=0; i<sph->Nsph; i++)
    sph->ComputeHydroForces(i,Nneib,neiblist,simparams);

  // Free up memory from local array
  delete[] neiblist;

  return;
}



// ============================================================================
// BruteForceSearch::UpdateAllGravityForces
// ============================================================================
void BruteForceSearch::UpdateAllGravityForces(Sph *sph, Parameters &simparams)
{
  int *neiblist;
  int Nneib;
  float agrav[ndim];
  float gpot;

  debug2("[BruteForceSearch::UpdateAllGravityForces]\n");

  // Allocate array to store local copy of potential neighbour ids
  Nneib = sph->Nsph;
  neiblist = new int[sph->Nsph];
  for (int i=0; i<sph->Nsph; i++) neiblist[i] = i;

  // Compute SPH hydro forces for all particles
  for (int i=0; i<sph->Nsph; i++)
    sph->ComputeGravForces(i,Nneib,neiblist);

  // Free up memory from local array
  delete[] neiblist;

  return;
}
