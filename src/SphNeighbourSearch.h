//=============================================================================
//  SphNeighbourSearch.h
//  Header file containing class definitions for all SPH neighbour searching 
//  data structures and algorithms.
//
//  This file is part of GANDALF :
//  Graphical Astrophysics code for N-body Dynamics And Lagrangian Fluids
//  https://github.com/gandalfcode/gandalf
//  Contact : gandalfcode@gmail.com
//
//  Copyright (C) 2013  D. A. Hubber, G. Rosotti
//
//  GANDALF is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 2 of the License, or
//  (at your option) any later version.
//
//  GANDALF is distributed in the hope that it will be useful, but
//  WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  General Public License (http://www.gnu.org/licenses) for more details.
//=============================================================================


#ifndef _SPH_NEIGHBOUR_SEARCH_H_
#define _SPH_NEIGHBOUR_SEARCH_H_

#include <assert.h>
#include <iostream>
#include <string>
#include <vector>
#include "Precision.h"
#include "Constants.h"
#include "SphKernel.h"
#include "SphParticle.h"
#include "Sph.h"
#include "Nbody.h"
#include "Sinks.h"
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
  FLOAT q[5];                       ///< Quadrupole moment tensor
  FLOAT bbmin[ndim];                ///< Minimum extent of bounding box
  FLOAT bbmax[ndim];                ///< Maximum extent of bounding box
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

  virtual void BuildTree(Sph<ndim> *, Parameters &) = 0;
  virtual void UpdateAllSphProperties(Sph<ndim> *, Nbody<ndim> *) = 0;
  virtual void UpdateAllSphForces(Sph<ndim> *) = 0;
  virtual void UpdateAllSphHydroForces(Sph<ndim> *) = 0;
  virtual void UpdateAllSphGravForces(Sph<ndim> *) = 0;
  virtual void UpdateAllSphDudt(Sph<ndim> *) = 0;
  virtual void UpdateAllSphDerivatives(Sph<ndim> *) = 0;
  virtual void UpdateTree(Sph<ndim> *, Parameters &) = 0;
  virtual void UpdateActiveParticleCounters(Sph<ndim> *) = 0;

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

  void BuildTree(Sph<ndim> *, Parameters &);
  void UpdateAllSphProperties(Sph<ndim> *, Nbody<ndim> *);
  void UpdateAllSphForces(Sph<ndim> *);
  void UpdateAllSphHydroForces(Sph<ndim> *);
  void UpdateAllSphGravForces(Sph<ndim> *);
  void UpdateAllSphDudt(Sph<ndim> *);
  void UpdateAllSphDerivatives(Sph<ndim> *);
  void UpdateTree(Sph<ndim> *, Parameters &);
  void UpdateActiveParticleCounters(Sph<ndim> *);

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

  void BuildTree(Sph<ndim> *, Parameters &);
  void UpdateAllSphProperties(Sph<ndim> *, Nbody<ndim> *);
  void UpdateAllSphForces(Sph<ndim> *);
  void UpdateAllSphHydroForces(Sph<ndim> *);
  void UpdateAllSphGravForces(Sph<ndim> *);
  void UpdateAllSphDudt(Sph<ndim> *);
  void UpdateAllSphDerivatives(Sph<ndim> *);
  void UpdateTree(Sph<ndim> *, Parameters &);
  void UpdateActiveParticleCounters(Sph<ndim> *);

  // Additional functions for grid neighbour search
  //---------------------------------------------------------------------------
  void AllocateGridMemory(int);
  void DeallocateGridMemory(void);
  void CreateGrid(Sph<ndim> *);
  int ComputeParticleGridCell(FLOAT *);
  void ComputeCellCoordinate(int, int *);
  int ComputeActiveCellList(int *);
  int ComputeActiveParticleList(int, int *, Sph<ndim> *);
  int ComputeNeighbourList(int, int *);
  int FindSplitAxis(int);
#if defined(VERIFY_ALL)
  void CheckValidNeighbourList(Sph<ndim> *,int,int,int *,string);
  void ValidateGrid(void);
#endif

  // Additional variables for grid
  //---------------------------------------------------------------------------
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
//  Class BinarySubTree
/// \brief   Class containing binary sub-tree
/// \details ..
/// \author  D. A. Hubber, A. P. Whitworth
/// \date    27/05/2013
//=============================================================================
template <int ndim>
class BinarySubTree
{
 public:

  BinarySubTree(int, FLOAT, FLOAT, string, string);
  ~BinarySubTree();

  // Additional functions for binary tree neighbour search
  //---------------------------------------------------------------------------
  void AllocateSubTreeMemory(void);
  void DeallocateSubTreeMemory(void);
  int ComputeActiveCellList(int ,BinaryTreeCell<ndim> **);
  void ComputeSubTreeSize(void);
  void CreateSubTreeStructure(void);
  void OrderParticlesByCartCoord(SphParticle<ndim> *);
  void LoadParticlesToSubTree(void);
  void StockCellProperties(SphParticle<ndim> *);
  FLOAT UpdateHmaxValues(SphParticle<ndim> *);
  void UpdateActiveParticleCounters(Sph<ndim> *);
  void BuildSubTree(Sph<ndim> *, Parameters &);
  int ComputeGatherNeighbourList(BinaryTreeCell<ndim> *, int, int, int *, FLOAT, SphParticle<ndim> *);
  int ComputeNeighbourList(BinaryTreeCell<ndim> *, int, int, int *, SphParticle<ndim> *);
  int ComputeGravityInteractionList(BinaryTreeCell<ndim> *, int, int, int, int &, int &, int &,
                                    int *, int *, BinaryTreeCell<ndim> **, SphParticle<ndim> *);
  int GlobalId(int local_id) {
	if (local_id < 0) cout << "local_id : " << local_id << endl;
    assert(local_id>=0);
    return ids[local_id];};

#if defined(VERIFY_ALL)
  void ValidateTree(Sph<ndim> *);
#endif

  // Additional variables for grid
  //---------------------------------------------------------------------------
  string gravity_mac;               ///< Multipole-acceptance criteria for tree
  string multipole;                 ///< Multipole-order for cell gravity
  bool allocated_tree;              ///< Are grid arrays allocated?
  int gtot;                         ///< Total number of grid/leaf cells
  int ifirst;                       ///< ..
  int ltot;                         ///< Total number of levels in tree
  int Ncell;                        ///< Current no. of grid cells
  int Ncellmax;                     ///< Max. allowed no. of grid cells
  int Ngridcells;                   ///< ""
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
  int *ids;                         ///< Particle ids
  int *inext;                       ///< Linked list for grid search
  int *porder[ndim];                ///< Ordered ids of Cartesian coordinates
  FLOAT kernrange;                  ///< Extent of employed kernel
  FLOAT theta;                      ///< ..
  FLOAT thetamaxsqd;                ///< ..
  FLOAT *rk[ndim];                  ///< Particle Cartesian coordinates
  BinaryTreeCell<ndim> *tree;       ///< Main tree array

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

  typedef typename vector <BinarySubTree<ndim> *>::iterator binlistiterator;

  BinaryTree(int, FLOAT, FLOAT, string, string);
  ~BinaryTree();

  //---------------------------------------------------------------------------
  void BuildTree(Sph<ndim> *, Parameters &);
  void UpdateAllSphProperties(Sph<ndim> *, Nbody<ndim> *);
  void UpdateAllSphForces(Sph<ndim> *);
  void UpdateAllSphHydroForces(Sph<ndim> *);
  void UpdateAllSphGravForces(Sph<ndim> *);
  void UpdateAllSphDudt(Sph<ndim> *);
  void UpdateAllSphDerivatives(Sph<ndim> *);
  void UpdateTree(Sph<ndim> *, Parameters &);
  void UpdateActiveParticleCounters(Sph<ndim> *);

  // Additional functions for binary tree neighbour search
  //---------------------------------------------------------------------------
  void AllocateTreeMemory(void);
  void DeallocateTreeMemory(void);
  void ComputeTreeSize(void);
  void CreateTreeStructure(void);
  void OrderParticlesByCartCoord(SphParticle<ndim> *);
  void LoadParticlesToTree(void);
  void UpdateHmaxValues(SphParticle<ndim> *);
  int ComputeActiveCellList(BinaryTreeCell<ndim> **, BinarySubTree<ndim> **);
  int ComputeActiveParticleList(BinaryTreeCell<ndim> *, BinarySubTree<ndim> *,
                                Sph<ndim> *, int *);
  int ComputeGatherNeighbourList(BinaryTreeCell<ndim> *, int, int *, 
                                 FLOAT, SphParticle<ndim> *);
  int ComputeNeighbourList(BinaryTreeCell<ndim> *, int, int *, 
                           SphParticle<ndim> *);
  int ComputeGravityInteractionList(BinaryTreeCell<ndim> *, int, int, int, 
                                    int &, int &, int &, int *, int *,
                                    BinaryTreeCell<ndim> **, 
                                    SphParticle<ndim> *);
  void ComputeCellMonopoleForces(int, int, BinaryTreeCell<ndim> **, 
                                 SphParticle<ndim> &);
  void ComputeCellQuadrupoleForces(int, int, BinaryTreeCell<ndim> **, 
                                   SphParticle<ndim> &);
#if defined(VERIFY_ALL)
  void CheckValidNeighbourList(Sph<ndim> *,int,int,int *,string);
#endif

  // Additional variables for grid
  //---------------------------------------------------------------------------
  string gravity_mac;               ///< Multipole-acceptance criteria for tree
  string multipole;                 ///< Multipole-order for cell gravity
  bool allocated_tree;              ///< Are grid arrays allocated?
  bool created_sub_trees;           ///< Have sub-tree objects been created?
  int gtot;                         ///< Total number of grid/leaf cells
  int ltot;                         ///< Total number of levels in tree
  int Ncell;                        ///< Current no. of grid cells
  int Ncellmax;                     ///< Max. allowed no. of grid cells
  int Nsubtree;                     ///< No. of sub-trees
  int Nsubtreemax;                  ///< ..
  int Nsubtreemaxold;               ///< No. of sub-trees
  int Nlevel;                       ///< ""
  int Nleafmax;                     ///< Max. number of particles per leaf cell
  int Nlistmax;                     ///< Max. length of neighbour list
  int Nsph;                         ///< Total no. of points/ptcls in grid
  int Ntot;                         ///< No. of current points in list
  int Ntotmax;                      ///< Max. no. of points in list
  int Ntotmaxold;                   ///< Old value of Ntotmax
  int *porder[ndim];                ///< Ordered ids of Cartesian coordinates
  FLOAT hmax;                       ///< Store hmax in the tree
  FLOAT kernrange;                  ///< Extent of employed kernel
  FLOAT theta;                      ///< ..
  FLOAT thetamaxsqd;                ///< ..
  int *pc;
  FLOAT *pw;                        ///< Particle weights
  FLOAT *rk[ndim];                  ///< Particle Cartesian coordinates
  BinaryTreeCell<ndim> *tree;       ///< Main tree array
  BinarySubTree<ndim> *stree;       ///< Array of sub-tree objects
  vector <BinarySubTree<ndim> *> subtrees;   ///< List containing pointers to sub-trees
};

#endif
