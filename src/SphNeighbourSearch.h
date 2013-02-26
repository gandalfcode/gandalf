// ============================================================================
// SphNeighbourSearch.h
// ============================================================================


#ifndef _SPH_NEIGHBOUR_SEARCH_H_
#define _SPH_NEIGHBOUR_SEARCH_H_


#include <iostream>
#include <string>
#include "Precision.h"
#include "Constants.h"
#include "Dimensions.h"
#include "SphKernel.h"
#include "SphParticle.h"
#include "Sph.h"
#include "Parameters.h"
using namespace std;

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

  virtual void UpdateAllSphProperties(Sph *) = 0;
  virtual void UpdateAllSphGravityProperties(Sph *) = 0;
  virtual void UpdateTree(Sph *, Parameters &) = 0;

  bool neibcheck;

};



// ============================================================================
// Class BruteForceSearch
// ============================================================================
class BruteForceSearch: public SphNeighbourSearch
{
 public:

  BruteForceSearch(int);
  ~BruteForceSearch();

  void UpdateAllSphProperties(Sph *);
  void UpdateAllSphGravityProperties(Sph *);
  void UpdateTree(Sph *, Parameters &);

#if !defined(FIXED_DIMENSIONS)
  int ndim;
#endif

};



// ============================================================================
// Class GridSearch
// ============================================================================
class GridSearch: public SphNeighbourSearch
{
 public:

  GridSearch(int);
  ~GridSearch();

  void UpdateAllSphProperties(Sph *);
  void UpdateAllSphGravityProperties(Sph *);
  void UpdateTree(Sph *, Parameters &);

  // Additional functions for grid neighbour search
  // --------------------------------------------------------------------------
  void AllocateGridMemory(int);
  void DeallocateGridMemory(void);
  void CreateGrid(Sph *);
  int ComputeParticleGridCell(FLOAT *);
  void ComputeCellCoordinate(int, int *);
  int ComputeActiveCellList(int *);
  int ComputeActiveParticleList(int, int *, Sph *);
  int ComputeNeighbourList(int, int *);
  void CheckValidNeighbourList(Sph *,int,int,int *,string);
  void ValidateGrid(void);

  // Additional variables for grid
  // --------------------------------------------------------------------------
  bool allocated_grid;                      // Are grid arrays allocated?
  int Ncell;                                // Current no. of grid cells
  int Ncellmax;                             // Max. allowed no. of grid cells
  int Ngrid[ndimmax];                       // No. of cells in each dimension
  int Noccupymax;                           // Max. occupancy of all cells
  int Nlistmax;                             // Max. length of neighbour list
  int Nsph;                                 // ..
  int Ntot;                                 // No. of current points in list
  int Ntotmax;                              // Max. no. of points in list
  int *inext;                               // Linked list for grid search
  FLOAT dx_grid;                            // Grid spacing
  FLOAT rmin[ndimmax];                      // Minimum extent of bounding box
  FLOAT rmax[ndimmax];                      // Maximum extent of bounding box
  GridCell *grid;                           // Main grid array
#if !defined(FIXED_DIMENSIONS)
  int ndim;
#endif

};


#endif
