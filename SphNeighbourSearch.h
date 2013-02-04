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


struct GridCell {
  int Nactive;
  int Nptcls;
  int ifirst;
  int ilast;
};


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

};



// ============================================================================
// Class GridSearch
// ============================================================================
class GridSearch: public SphNeighbourSearch
{
 public:

  GridSearch();
  ~GridSearch();

  void UpdateAllSphProperties(Sph *, Parameters &);
  void UpdateAllSphForces(Sph *, Parameters &);

  // Additional functions for grid neighbour search
  // --------------------------------------------------------------------------
  void AllocateGridMemory(int);
  void DeallocateGridMemory(void);
  void CreateGrid(Sph *);
  int ComputeParticleGridCell(float *);
  void ComputeCellCoordinate(int, int *);
  int ComputeActiveCellList(int *);
  int ComputeActiveParticleList(int, int *);
  int ComputeNeighbourList(int, int *);

  // Additional variables for grid
  // --------------------------------------------------------------------------
  bool allocated_grid;                      // Are grid arrays allocated?
  int Ncell;                                // Current no. of grid cells
  int Ncellmax;                             // Max. allowed no. of grid cells
  int Ngrid[ndimmax];                       // No. of cells in each dimension
  int Noccupymax;                           // Max. occupancy of all cells
  int Nlistmax;                             // Max. length of neighbour list
  int Ntot;                                 // No. of current points in list
  int Ntotmax;                              // Max. no. of points in list
  int *inext;                               // Linked list for grid search
  float dx_grid;                            // Grid spacing
  float rmin[ndimmax];                      // Minimum extent of bounding box
  float rmax[ndimmax];                      // Maximum extent of bounding box
  GridCell *grid;                           // Main grid array

};


#endif
