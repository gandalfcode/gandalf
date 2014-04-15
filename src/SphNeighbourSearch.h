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
#include "CodeTiming.h"
#include "SphKernel.h"
#include "SphParticle.h"
#include "Sph.h"
#include "Nbody.h"
#include "Sinks.h"
#include "DomainBox.h"
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
  int c1;                           ///< First child cell
  int c2;                           ///< Second child cell
  int c2g;                          ///< i.d. of tree-cell c/grid-cell g
  int cnext;                        ///< i.d. of next cell if not opened
  int id;                           ///< Cell id
  int k_divide;                     ///< Dimension along which cell is split
  int level;                        ///< Level of cell on tree
  int ifirst;                       ///< i.d. of first particle in cell
  int ilast;                        ///< i.d. of last particle in cell
  int N;                            ///< ..
  int Nactive;                      ///< ..
  FLOAT cdistsqd;                   ///< ..
  FLOAT mac;                        ///< Multipole-opening criterion value
  FLOAT bbmin[ndim];                ///< Minimum extent of bounding box
  FLOAT bbmax[ndim];                ///< Maximum extent of bounding box
  FLOAT hboxmin[ndim];              ///< Minimum extent of bounding box
  FLOAT hboxmax[ndim];              ///< Maximum extent of bounding box
  FLOAT rcell[ndim];                ///< ..
  FLOAT r[ndim];                    ///< Position of cell
  FLOAT v[ndim];                    ///< Position of cell
  FLOAT m;                          ///< Mass contained in cell
  FLOAT rmax;                       ///< Radius of bounding sphere
  FLOAT hmax;                       ///< Maximum smoothing length inside cell
  FLOAT drmaxdt;                    ///< Rate of change of bounding sphere
  FLOAT dhmaxdt;                    ///< Rate of change of maximum h
  FLOAT q[5];                       ///< Quadrupole moment tensor
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

  SphNeighbourSearch();
  ~SphNeighbourSearch();

  virtual void BuildTree(bool, int, int, int, FLOAT, Sph<ndim> *) = 0;
  virtual void UpdateAllSphProperties(Sph<ndim> *, Nbody<ndim> *) = 0;
  virtual void UpdateAllSphForces(Sph<ndim> *, Nbody<ndim> *) = 0;
  virtual void UpdateAllSphHydroForces(Sph<ndim> *, Nbody<ndim> *) = 0;
  virtual void UpdateAllSphGravForces(Sph<ndim> *, Nbody<ndim> *) = 0;
  virtual void UpdateAllSphDudt(Sph<ndim> *) = 0;
  virtual void UpdateAllSphDerivatives(Sph<ndim> *) = 0;
  virtual void UpdateActiveParticleCounters(Sph<ndim> *) = 0;
  virtual void UpdateAllStarGasForces(Sph<ndim> *, Nbody<ndim> *) = 0;
#if defined(VERIFY_ALL)
  void CheckValidNeighbourList(Sph<ndim> *,int, int, int *, string);
#endif

  bool neibcheck;                   ///< Flag to verify neighbour lists
  CodeTiming *timing;               ///< Pointer to code timing object
  DomainBox<ndim> *box;             ///< Pointer to simulation bounding box

};

#if defined MPI_PARALLEL
//Forward declare MpiNode to break circular dependency
template <int ndim>
class MpiNode;
#endif


//=============================================================================
//  Class BruteForceSearch
/// Class for computing SPH neighbour lists using brute force only 
/// (i.e. direct summation over all particles).
//=============================================================================
template <int ndim>
class BruteForceSearch: public SphNeighbourSearch<ndim>
{
  using SphNeighbourSearch<ndim>::neibcheck;
  using SphNeighbourSearch<ndim>::timing;
  
 public:

  BruteForceSearch();
  ~BruteForceSearch();

  void BuildTree(bool, int, int, int, FLOAT, Sph<ndim> *);
  void UpdateAllSphProperties(Sph<ndim> *, Nbody<ndim> *);
  void UpdateAllSphForces(Sph<ndim> *, Nbody<ndim> *);
  void UpdateAllSphHydroForces(Sph<ndim> *, Nbody<ndim> *);
  void UpdateAllSphGravForces(Sph<ndim> *, Nbody<ndim> *);
  void UpdateAllSphDudt(Sph<ndim> *);
  void UpdateAllSphDerivatives(Sph<ndim> *);
  void UpdateActiveParticleCounters(Sph<ndim> *);
  void UpdateAllStarGasForces(Sph<ndim> *, Nbody<ndim> *);
#if defined MPI_PARALLEL
  void FindGhostParticlesToExport(Sph<ndim>* sph, std::vector<std::vector<SphParticle<ndim>* > >&,
      const std::vector<int>&, MpiNode<ndim>*);
  void FindParticlesToTransfer(Sph<ndim>* sph, std::vector<std::vector<int> >& particles_to_export,
      std::vector<int>& all_particles_to_export, const std::vector<int>& potential_nodes, MpiNode<ndim>* mpinodes);
#endif
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
  using SphNeighbourSearch<ndim>::timing;

 public:

  GridSearch();
  ~GridSearch();

  void BuildTree(bool, int, int, int, FLOAT, Sph<ndim> *);
  void UpdateAllSphProperties(Sph<ndim> *, Nbody<ndim> *);
  void UpdateAllSphForces(Sph<ndim> *, Nbody<ndim> *);
  void UpdateAllSphHydroForces(Sph<ndim> *, Nbody<ndim> *);
  void UpdateAllSphGravForces(Sph<ndim> *, Nbody<ndim> *);
  void UpdateAllSphDudt(Sph<ndim> *);
  void UpdateAllSphDerivatives(Sph<ndim> *);
  void UpdateActiveParticleCounters(Sph<ndim> *);
  void UpdateAllStarGasForces(Sph<ndim> *, Nbody<ndim> *);

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
//  Class BinaryTree
/// \brief   Class containing binary tree
/// \details Binary tree data structure used for efficient neighbour searching 
///          and computation of gravitational forces
/// \author  D. A. Hubber, O. Lomax, A. P. Whitworth
/// \date    08/01/2014
//=============================================================================
template <int ndim>
class BinaryTree: public SphNeighbourSearch<ndim>
{
 public:

  using SphNeighbourSearch<ndim>::neibcheck;
  using SphNeighbourSearch<ndim>::box;
  using SphNeighbourSearch<ndim>::timing;

  BinaryTree(int, FLOAT, FLOAT, FLOAT, string, string);
  ~BinaryTree();

  //---------------------------------------------------------------------------
  void BuildTree(bool, int, int, int, FLOAT, Sph<ndim> *);
  void UpdateAllSphProperties(Sph<ndim> *, Nbody<ndim> *);
  void UpdateAllSphForces(Sph<ndim> *, Nbody<ndim> *);
  void UpdateAllSphHydroForces(Sph<ndim> *, Nbody<ndim> *);
  void UpdateAllSphGravForces(Sph<ndim> *, Nbody<ndim> *);
  void UpdateAllSphDudt(Sph<ndim> *);
  void UpdateAllSphDerivatives(Sph<ndim> *);
  void UpdateActiveParticleCounters(Sph<ndim> *);
  void UpdateAllStarGasForces(Sph<ndim> *, Nbody<ndim> *);


  // Additional functions for binary tree neighbour search
  //---------------------------------------------------------------------------
  void AllocateTreeMemory(Sph<ndim> *);
  void DeallocateTreeMemory(void);
  bool BoxOverlap(FLOAT *, FLOAT *, FLOAT *, FLOAT *);
  void ComputeTreeSize(void);
  void CreateTreeStructure(void);
  void DivideTreeCell(int, int, Sph<ndim> *, BinaryTreeCell<ndim> &);
  //void ExtrapolateCellProperties(BinaryTreeCell<ndim> &, FLOAT);
  void ExtrapolateCellProperties(FLOAT);
  FLOAT QuickSelect(int, int, int, int, Sph<ndim> *);
  void StockTree(BinaryTreeCell<ndim> &, SphParticle<ndim> *);
  void StockCellProperties(BinaryTreeCell<ndim> &, SphParticle<ndim> *);
  void UpdateHmaxValues(BinaryTreeCell<ndim> &, SphParticle<ndim> *);
  int ComputeActiveCellList(BinaryTreeCell<ndim> **);
  int ComputeActiveParticleList(BinaryTreeCell<ndim> *, Sph<ndim> *, int *);
  int ComputeGatherNeighbourList(BinaryTreeCell<ndim> *, int, int *, 
                                 FLOAT, SphParticle<ndim> *);
  int ComputeNeighbourList(BinaryTreeCell<ndim> *, int, int *, 
                           SphParticle<ndim> *);
  int ComputeGravityInteractionList(BinaryTreeCell<ndim> *, FLOAT, int, int,  
                                    int, int &, int &, int &, int *, int *,
                                    BinaryTreeCell<ndim> **, 
				    SphParticle<ndim> *);
  int ComputeStarGravityInteractionList(NbodyParticle<ndim> *, FLOAT, int, int,
					int, int &, int &, int &, int *, int *,
					BinaryTreeCell<ndim> **, 
					SphParticle<ndim> *);
  void ComputeCellMonopoleForces(FLOAT &, FLOAT *, FLOAT *, int, 
				 BinaryTreeCell<ndim> **);
  void ComputeCellQuadrupoleForces(FLOAT &, FLOAT *, FLOAT *, int, 
				   BinaryTreeCell<ndim> **);
  void ComputeFastMonopoleForces(int, int, BinaryTreeCell<ndim> **, 
				 BinaryTreeCell<ndim> *, SphParticle<ndim> *);
#if defined(VERIFY_ALL)
  void ValidateTree(Sph<ndim> *);
#endif

  // Additional variables for grid
  //---------------------------------------------------------------------------
  string gravity_mac;               ///< Multipole-acceptance criteria for tree
  string multipole;                 ///< Multipole-order for cell gravity
  bool allocated_tree;              ///< Are grid arrays allocated?
  int gmax;                         ///< Max. no. of grid/leaf cells
  int gtot;                         ///< Total number of grid/leaf cells
  int lmax;                         ///< Max. no. of levels
  int ltot;                         ///< Total number of levels in tree
  int ltot_old;                     ///< Prev. value of ltot
  int Ncell;                        ///< Current no. of grid cells
  int Ncellmax;                     ///< Max. allowed no. of grid cells
  int Nlevel;                       ///< ""
  int Nleafmax;                     ///< Max. number of particles per leaf cell
  int Nlistmax;                     ///< Max. length of neighbour list
  int Nsph;                         ///< Total no. of points/ptcls in grid
  int Ntot;                         ///< No. of current points in list
  int Ntotold;                      ///< Prev. no. of particles
  int Ntotmax;                      ///< Max. no. of points in list
  int Ntotmaxold;                   ///< Old value of Ntotmax
  int *g2c;                         ///< i.d. of leaf(grid) cells
  int *ids;                         ///< Particle ids
  int *inext;                       ///< Linked list for grid search
  FLOAT macerror;                   ///< Error tolerance for gravity tree-MAC
  FLOAT hmax;                       ///< Store hmax in the tree
  FLOAT kernrange;                  ///< Extent of employed kernel
  FLOAT theta;                      ///< Geometric opening angle
  FLOAT thetamaxsqd;                ///< Geometric opening angle squared
  FLOAT invthetamaxsqd;             ///< 1 / thetamaxsqd
  BinaryTreeCell<ndim> *tree;       ///< Main tree array


  bool allocated_buffer;
  int Nthreads;
  int *Nneibmaxbuf;
  int *Ndirectmaxbuf;
  int *Ngravcellmaxbuf;
  int **activelistbuf;
  int **levelneibbuf;
  SphParticle<ndim> **neibpartbuf;   // Local copy of neighbouring ptcls
  SphParticle<ndim> **activepartbuf; // Local copy of SPH particle  


};


#endif
