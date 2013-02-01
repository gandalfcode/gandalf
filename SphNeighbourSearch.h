// ============================================================================
// SphNeighbourSearch.h
// ============================================================================


#ifndef _SPH_NEIGHBOUR_SEARCH_H_
#define _SPH_NEIGHBOUR_SEARCH_H_


#include "Constants.h"
#include "Dimensions.h"
#include "SphKernel.h"
#include "SphParticle.h"
#include "Sph.h"
#include "Parameters.h"



// ============================================================================
// Class SphNeighbourSearch
// ============================================================================
class SphNeighbourSearch
{
 public:

  SphNeighbourSearch();
  ~SphNeighbourSearch();

  virtual void UpdateAllSphProperties(Sph *, Parameters &) = 0;
  virtual void UpdateAllSphForces(Sph *, Parameters &) = 0;
  virtual void UpdateAllGravityForces(Sph *, Parameters &) = 0;

};



// ============================================================================
// Class BruteForceSearch
// ============================================================================
class BruteForceSearch: public SphNeighbourSearch
{
 public:

  BruteForceSearch();
  ~BruteForceSearch();

  void UpdateAllSphProperties(Sph *, Parameters &);
  void UpdateAllSphForces(Sph *, Parameters &);
  void UpdateAllGravityForces(Sph *, Parameters &);

};


#endif
