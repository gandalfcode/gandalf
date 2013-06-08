//=============================================================================
//  SphNeighbourSearch.h
//  Header file containing class definitions for all SPH neighbour searching 
//  data structures and algorithms.
//=============================================================================


#ifndef _SPH_NEIGHBOUR_SEARCH_H_
#define _SPH_NEIGHBOUR_SEARCH_H_


#include <iostream>
#include <string>
#include "Precision.h"
#include "Constants.h"
#include "SphKernel.h"
#include "SphParticle.h"
#include "Sph.h"
#include "Parameters.h"
using namespace std;


//=============================================================================
//  Structure GridCell
/// Neighbour grid cell data structure
//=============================================================================
struct GridCell {
  int Nactive;                      ///< No. of active particles in grid cell
  int Nptcls;                       ///< Total no. of particles in grid cell
  int ifirst;                       ///< i.d. of first particle in cell
  int ilast;                        ///< i.d. of last particle in cell
};



//=============================================================================
//  Structure BinaryTreeCell
/// Neighbour grid cell data structure
//=============================================================================
template <int ndim>
struct BinaryTreeCell {
  int c2;                           ///< i.d. of 2nd child cell
  int cnext;                        ///< i.d. of next cell if not opened
  int c2g;                          ///< i.d. of tree-cell c/grid-cell g
  int clevel;                       ///< Level occupied by tree-cell
  int ifirst;                       ///< i.d. of 1st particle allocated to cell
  int ilast;                        ///< i.d. of last particle in cell
  int Nactive;                      ///< No. of active particles in cell
  int N;                            ///< No. of particles in cell
  FLOAT cdistsqd;                   ///< Opening distances squared
  FLOAT r[ndim];                    ///< Position of centre of mass
  FLOAT m;                          ///< Total mass of cell
  FLOAT rmax;                       ///< Max. dist. of ptcl from COM
  FLOAT hmax;                       ///< Maximum smoothing length inside cell
};



//=============================================================================
//  Class SphNeighbourSearch
/// \brief   SphNeighbourSearch class definition.  
/// \details Contains routines for creating the SPH neighbour search data
///          structure, and for computing local neighbour lists and calling 
///          SPH functions (e.g. computing h, SPH forces, etc..).
/// \author  D. A. Hubber, G. Rosotti
/// \date    03/04/2013
//=============================================================================
template <int ndim>
class SphNeighbourSearch
{
 public:

  virtual void UpdateAllSphProperties(Sph<ndim> *) = 0;
  virtual void UpdateAllSphForces(Sph<ndim> *) = 0;
  virtual void UpdateAllSphHydroForces(Sph<ndim> *) = 0;
  virtual void UpdateAllSphGravForces(Sph<ndim> *) = 0;
  virtual void UpdateAllSphDudt(Sph<ndim> *) = 0;
  virtual void UpdateAllSphDerivatives(Sph<ndim> *) = 0;
  virtual void UpdateTree(Sph<ndim>*, Parameters &) = 0;

  bool neibcheck;

};



//=============================================================================
//  Class BruteForceSearch
/// Class for computing SPH neighbour lists using brute force only 
/// (i.e. direct summation over all particles).
//=============================================================================
template <int ndim>
class BruteForceSearch: public SphNeighbourSearch<ndim>
{
  using SphNeighbourSearch<ndim>::neibcheck;

 public:

  BruteForceSearch();
  ~BruteForceSearch();

  void UpdateAllSphProperties(Sph<ndim> *);
  void UpdateAllSphForces(Sph<ndim> *);
  void UpdateAllSphHydroForces(Sph<ndim> *);
  void UpdateAllSphGravForces(Sph<ndim> *);
  void UpdateAllSphDudt(Sph<ndim> *);
  void UpdateAllSphDerivatives(Sph<ndim> *);
  void UpdateTree(Sph<ndim> *, Parameters &);

};



//=============================================================================
//  Class GridSearch
/// Class for computing SPH neighbour lists using a uniform grid.  The size 
/// of the grid is the maximum kernel extent (e.g. 2*h_max for the M4 kernel)
/// multiplied by some tolerance.
//=============================================================================
template <int ndim>
class GridSearch: public SphNeighbourSearch<ndim>
{
  using SphNeighbourSearch<ndim>::neibcheck;

 public:

  GridSearch();
  ~GridSearch();

  void UpdateAllSphProperties(Sph<ndim> *);
  void UpdateAllSphForces(Sph<ndim> *);
  void UpdateAllSphHydroForces(Sph<ndim> *);
  void UpdateAllSphGravForces(Sph<ndim> *);
  void UpdateAllSphDudt(Sph<ndim> *);
  void UpdateAllSphDerivatives(Sph<ndim> *);
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
  bool allocated_grid;              ///< Are grid arrays allocated?
  int Ncell;                        ///< Current no. of grid cells
  int Ncellmax;                     ///< Max. allowed no. of grid cells
  int Ngrid[ndim];                  ///< No. of cells in each dimension
  int Noccupymax;                   ///< Max. occupancy of all cells
  int Nlistmax;                     ///< Max. length of neighbour list
  int Nsph;                         ///< Total no. of points/ptcls in grid
  int Ntot;                         ///< No. of current points in list
  int Ntotmax;                      ///< Max. no. of points in list
  int *inext;                       ///< Linked list for grid search
  FLOAT dx_grid;                    ///< Grid spacing
  FLOAT rmin[ndim];                 ///< Minimum extent of bounding box
  FLOAT rmax[ndim];                 ///< Maximum extent of bounding box
  GridCell *grid;                   ///< Main grid array

};



//=============================================================================
//  Class BinaryTree
/// \brief   Class containing binary tree
/// \details Binary tree data structure used for efficient neighbour searching 
///          and computation of gravitational forces
/// \author  D. A. Hubber, A. P. Whitworth
/// \date    27/05/2013
//=============================================================================
template <int ndim>
class BinaryTree: public SphNeighbourSearch<ndim>
{
  using SphNeighbourSearch<ndim>::neibcheck;

 public:

  BinaryTree(int, FLOAT, string);
  ~BinaryTree();

  void UpdateAllSphProperties(Sph<ndim> *);
  void UpdateAllSphForces(Sph<ndim> *);
  void UpdateAllSphHydroForces(Sph<ndim> *);
  void UpdateAllSphGravForces(Sph<ndim> *);
  void UpdateAllSphDudt(Sph<ndim> *);
  void UpdateAllSphDerivatives(Sph<ndim> *);
  void UpdateTree(Sph<ndim> *, Parameters &);

  // Additional functions for grid neighbour search
  // --------------------------------------------------------------------------
  void AllocateTreeMemory(int);
  void DeallocateTreeMemory(void);
  void ComputeTreeSize(int, int);
  void CreateTreeStructure(int);
  void OrderParticlesByCartCoord(int, SphParticle<ndim> *);
  void LoadParticlesToTree(int);
  void StockCellProperties(int, SphParticle<ndim> *);
  void UpdateHmaxValues(int, SphParticle<ndim> *);
  void CreateGridCells(int, SphParticle<ndim> *);
  int ComputeActiveCellList(int *);
  int ComputeActiveParticleList(int, int *, Sph<ndim> *);
  int ComputeGatherNeighbourList(int, int, int *, FLOAT);
  int ComputeNeighbourList(int, int, int *);
  int ComputeGravityInteractionList(int, int, int, int, int &, int &, int &,
                                    int *, int *, int *);
  void ComputeCellForces(int, int, int *, SphParticle<ndim> &);
#if defined(VERIFY_ALL)
  void CheckValidNeighbourList(Sph<ndim> *,int,int,int *,string);
  void ValidateTree(Sph<ndim> *);
#endif

  // Additional variables for grid
  // --------------------------------------------------------------------------
  string gravity_mac;               ///< ..
  bool allocated_tree;              ///< Are grid arrays allocated?
  int Ncell;                        ///< Current no. of grid cells
  int Ncellmax;                     ///< Max. allowed no. of grid cells
  int gtot;                         ///< Total number of grid/leaf cells
  int Ngridcells;                   ///< ""
  int ltot;                         ///< Total number of levels in tree
  int Nlevel;                       ///< ""
  int Nleafmax;                     ///< Max. number of particles per leaf cell
  int Nlistmax;                     ///< Max. length of neighbour list
  int Nsph;                         ///< Total no. of points/ptcls in grid
  int Ntot;                         ///< No. of current points in list
  int Ntotmax;                      ///< Max. no. of points in list
  int Ntotmaxold;                   ///< Old value of Ntotmax
  int *pc;                          ///< i.d. of leaf cell occupied by ptcl
  int *g2c;                         ///< i.d. of leaf(grid) cells
  int *gactivelist;                 ///< List of active cells
  int *inext;                       ///< Linked list for grid search
  int *porder[ndim];                ///< Ordered ids of Cartesian coordinates
  FLOAT *pw;                        ///< Particle weights
  FLOAT *r[ndim];                   ///< Particle Cartesian coordinates
  FLOAT theta;                      ///< ..
  FLOAT thetamaxsqd;                ///< ..
  BinaryTreeCell<ndim> *tree;       ///< Main tree array

};


#endif
