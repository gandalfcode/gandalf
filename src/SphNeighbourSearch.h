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
template <int ndim>
class SphNeighbourSearch
{
 public:

  virtual void UpdateAllSphProperties(Sph<ndim> *) = 0;
  virtual void UpdateAllSphForces(Sph<ndim> *) = 0;
  virtual void UpdateAllSphDudt(Sph<ndim> *) = 0;
  virtual void UpdateAllSphDerivatives(Sph<ndim> *) = 0;
  virtual void UpdateAllSphGravityProperties(Sph<ndim> *) = 0;
  virtual void UpdateTree(Sph<ndim>*, Parameters &) = 0;

  bool neibcheck;

};



// ============================================================================
// Class BruteForceSearch
// ============================================================================
template <int ndim>
class BruteForceSearch: public SphNeighbourSearch<ndim>
{
 public:

  BruteForceSearch();
  ~BruteForceSearch();

  void UpdateAllSphProperties(Sph<ndim> *);
  void UpdateAllSphForces(Sph<ndim> *);
  void UpdateAllSphDudt(Sph<ndim> *);
  void UpdateAllSphDerivatives(Sph<ndim> *);
  void UpdateAllSphGravityProperties(Sph<ndim> *);
  void UpdateTree(Sph<ndim> *, Parameters &);

//#if !defined(FIXED_DIMENSIONS)
//  int ndim;
//#endif

};



// ============================================================================
// Class GridSearch
// ============================================================================
template <int ndim>
class GridSearch: public SphNeighbourSearch<ndim>
{
 public:

  GridSearch();
  ~GridSearch();

  void UpdateAllSphProperties(Sph<ndim> *);
  void UpdateAllSphForces(Sph<ndim> *);
  void UpdateAllSphDudt(Sph<ndim> *);
  void UpdateAllSphDerivatives(Sph<ndim> *);
  void UpdateAllSphGravityProperties(Sph<ndim> *);
  void UpdateTree(Sph<ndim> *, Parameters &);

  // Additional functions for grid neighbour search
  // --------------------------------------------------------------------------
  void AllocateGridMemory(int);
  void DeallocateGridMemory(void);
  void CreateGrid(Sph<ndim> *);
  int ComputeParticleGridCell(FLOAT *);
  void ComputeCellCoordinate(int, int *);
  int ComputeActiveCellList(int *);
  int ComputeActiveParticleList(int, int *, Sph<ndim> *);
  int ComputeNeighbourList(int, int *);
#if defined(VERIFY_ALL)
  void CheckValidNeighbourList(Sph<ndim> *,int,int,int *,string);
  void ValidateGrid(void);
#endif

  // Additional variables for grid
  // --------------------------------------------------------------------------
  bool allocated_grid;                      // Are grid arrays allocated?
  int Ncell;                                // Current no. of grid cells
  int Ncellmax;                             // Max. allowed no. of grid cells
  int Ngrid[ndim];                       // No. of cells in each dimension
  int Noccupymax;                           // Max. occupancy of all cells
  int Nlistmax;                             // Max. length of neighbour list
  int Nsph;                                 // ..
  int Ntot;                                 // No. of current points in list
  int Ntotmax;                              // Max. no. of points in list
  int *inext;                               // Linked list for grid search
  FLOAT dx_grid;                            // Grid spacing
  FLOAT rmin[ndim];                      // Minimum extent of bounding box
  FLOAT rmax[ndim];                      // Maximum extent of bounding box
  GridCell *grid;                           // Main grid array
//#if !defined(FIXED_DIMENSIONS)
//  int ndim;
//#endif

};


#endif
