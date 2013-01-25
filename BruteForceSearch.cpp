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
  int *neiblist;
  int Nneib;

  debug2("[BruteForceSearch::UpdateAllSphProperties]");

  Nneib = sph->Nsph;
  neiblist = new int[sph->Nsph];

  for (int i=0; i<sph->Nsph; i++) neiblist[i] = i;

  for (int i=0; i<sph->Nsph; i++) sph->ComputeH(i,Nneib,neiblist,simparams);

  for (int i=0; i<sph->Nsph; i++) sph->ComputeSphProperties(i,Nneib,neiblist,simparams);

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

  debug2("[BruteForceSearch::UpdateAllSphForces]");

  Nneib = sph->Nsph;
  neiblist = new int[sph->Nsph];

  for (int i=0; i<sph->Nsph; i++) neiblist[i] = i;

  for (int i=0; i<sph->Nsph; i++) sph->ComputeHydroForces(i,Nneib,neiblist,simparams);

  delete[] neiblist;

  return;
}
