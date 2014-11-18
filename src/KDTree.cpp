//=================================================================================================
//  KDTree.cpp
//  Contains all functions for building, stocking and walking for the
//  binary KD-tree for SPH particles.
//  Based on code courtesy of O. Lomax and A. Whitworth
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
//=================================================================================================


#include <cstdlib>
#include <cassert>
#include <iostream>
#include <string>
#include <math.h>
#include "Precision.h"
#include "Exception.h"
#include "Parameters.h"
#include "InlineFuncs.h"
#include "SphParticle.h"
#include "Sph.h"
#include "KDTree.h"
#include "Debug.h"
#if defined _OPENMP
#include <omp.h>
#endif
using namespace std;



//=================================================================================================
//  KDTree::KDTree
/// KDTree constructor.  Initialises various variables.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
KDTree<ndim,ParticleType,TreeCell>::KDTree(int Nleafmaxaux, FLOAT thetamaxsqdaux,
                                           FLOAT kernrangeaux, FLOAT macerroraux,
                                           string gravity_mac_aux, string multipole_aux)
{
  allocated_tree = false;
  lactive        = 0;
  ltot           = 0;
  Ncell          = 0;
  Ncellmax       = 0;
  Ntot           = 0;
  Ntotmax        = 0;
  Ntotmaxold     = 0;
  Nleafmax       = Nleafmaxaux;
  kernrange      = kernrangeaux;
  thetamaxsqd    = thetamaxsqdaux;
  invthetamaxsqd = 1.0/thetamaxsqdaux;
  gravity_mac    = gravity_mac_aux;
  macerror       = macerroraux;
  multipole      = multipole_aux;
#if defined _OPENMP
  Nthreads       = omp_get_max_threads();
#else
  Nthreads       = 1;
#endif
#if defined MPI_PARALLEL
  Ncelltot       = 0;
  Nimportedcell  = 0;
#endif
}



//=================================================================================================
//  KDTree::~KDTree
/// KDTree destructor.  Deallocates tree memory upon object destruction.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
KDTree<ndim,ParticleType,TreeCell>::~KDTree()
{
  if (allocated_tree) DeallocateTreeMemory();
}



//=================================================================================================
//  KDTree::AllocateTreeMemory
/// Allocate memory for KD-tree as requested.  If more memory is required
/// than currently allocated, tree is deallocated and reallocated here.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void KDTree<ndim,ParticleType,TreeCell>::AllocateTreeMemory(void)
{
  int ithread;                      // Thread id number

  debug2("[KDTree::AllocateTreeMemory]");

  ComputeTreeSize();

  if (!allocated_tree || Ntotmax > Ntotmaxold) {
    if (allocated_tree) DeallocateTreeMemory();
    Ntotmax = max(Ntotmax,Ntot);
    Ntotmaxold = Ntotmax;

    g2c = new int[gmax];
    ids = new int[Ntotmax];
    inext = new int[Ntotmax];
    celldata = new struct TreeCell<ndim>[Ncellmax];

    allocated_tree = true;

    //CreateTreeStructure();

  }

  CreateTreeStructure();

  return;
}



//=================================================================================================
//  KDTree::DeallocateTreeMemory
/// Deallocates all KD-tree memory
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void KDTree<ndim,ParticleType,TreeCell>::DeallocateTreeMemory(void)
{
  debug2("[KDTree::DeallocateTreeMemory]");

  if (allocated_tree) {
    delete[] celldata;
    delete[] inext;
    delete[] ids;
    delete[] g2c;
    allocated_tree = false;
  }

  return;
}



//=================================================================================================
//  KDTree::BuildTree
/// Call all routines to build/re-build the KD-tree on the local node.
/// If OpenMP is activated, the local domain is partitioned into sub-trees
/// in order to improve the scalability of building and stocking the tree.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void KDTree<ndim,ParticleType,TreeCell>::BuildTree
 (int _ifirst,                         ///< i.d. of first particle
  int _ilast,                          ///< i.d. of last particle
  int Npart,                           ///< No. of particles
  int Npartmax,                        ///< Max. no. of particles
  ParticleType<ndim> *partdata,        ///< Particle data array
  FLOAT timestep)                      ///< Smallest physical timestep
{
  int i;                               // Particle counter
  int k;                               // Dimension counter

  debug2("[KDTree::BuildTree]");
  //timing->StartTimingSection("BUILD_TREE",2);

  // Activate nested parallelism for tree building routines
#ifdef _OPENMP
  omp_set_nested(1);
#endif

  // Set no. of tree members to total number of SPH particles (inc. ghosts)
  gtot = 0;
  ltot_old = ltot;

  // Compute the size of all tree-related arrays now we know number of points
  ComputeTreeSize();

  // Allocate (or reallocate if needed) all tree memory
  AllocateTreeMemory();

  // If the number of levels in the tree has changed (due to destruction or
  // creation of new particles) then re-create tree data structure
  // including linked lists and cell pointers
  if (ltot != ltot_old) CreateTreeStructure();


  // Set properties for root cell before constructing tree
  ifirst = _ifirst;
  ilast  = _ilast;
  celldata[0].N      = Ntot;
  celldata[0].ifirst = ifirst;
  celldata[0].ilast  = ilast;
  celldata[0].cnext  = Ncellmax;
  for (k=0; k<ndim; k++) celldata[0].bbmin[k] = -big_number;
  for (k=0; k<ndim; k++) celldata[0].bbmax[k] = big_number;
  for (k=0; k<ndim; k++) celldata[0].cexit[0][k] = -1;
  for (k=0; k<ndim; k++) celldata[0].cexit[1][k] = -1;
  for (i=ifirst; i<ilast; i++) inext[i] = i+1;
  inext[ilast] = -1;

  // If number of particles remains unchanged, use old id list
  // (nearly sorted list should be faster for quick select).
  if (Ntot > 0) {
    if (Ntot != Ntotold)
      for (i=ifirst; i<=ilast; i++) ids[i] = i;

    // Recursively build tree from root node down
    DivideTreeCell(ifirst,ilast,partdata,celldata[0]);
#if defined(VERIFY_ALL)
    ValidateTree(partdata);
#endif
  }

  return;
}



//=================================================================================================
//  KDTree::ComputeTreeSize
/// Compute the maximum size (i.e. no. of levels, cells and leaf cells) of the KD-tree.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void KDTree<ndim,ParticleType,TreeCell>::ComputeTreeSize(void)
{
  debug2("[KDTree::ComputeTreeSize]");


  // Calculate maximum level of tree that can contain max. no. of particles
  lmax = 0;
  while (Nleafmax*pow(2,lmax) < Ntotmax) {
    lmax++;
  };
  gmax = pow(2,lmax);
  Ncellmax = 2*gmax - 1;

  // Calculate level of tree that can contain all current particles
  lactive = 0;
  while (Nleafmax*pow(2,lactive) < Ntot) {
    lactive++;
  };
  gactive = pow(2,lactive);

  // Calculate level of tree that can contain all current particles
  ltot = 0;
  while (Nleafmax*pow(2,ltot) < Ntot) {
    ltot++;
  };
  gtot = pow(2,ltot);
  Ncell = 2*gtot - 1;


  // Optional output (for debugging)
#if defined(VERIFY_ALL)
  cout << "No. of ptcls in tree  : " << Ntot << "   " << Ntotmax << endl;
  cout << "No. of grid-cells     : " << gtot << "   " << gmax << endl;
  cout << "No. of levels on tree : " << ltot << "   " << lmax << "    " << lactive << endl;
  cout << "No. of cells in tree  : " << Ncell << "   " << Ncellmax << endl;
#endif

  return;
}



//=================================================================================================
//  KDTree::CreateTreeStructure
/// Create the raw tree skeleton structure once the tree size is known.
/// Sets all cell pointer variables and all cell levels.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void KDTree<ndim,ParticleType,TreeCell>::CreateTreeStructure(void)
{
  int c;                            // Dummy id of tree-level, then tree-cell
  int g;                            // Dummy id of grid-cell
  int l;                            // Dummy id of level
  int *c2L;                         // Increment to second child-cell
  int *cNL;                         // Increment to next cell if cell unopened

  debug2("[KDTree::CreateTreeStructure]");

  // Allocate memory for local arrays
  c2L = new int[ltot + 1];
  cNL = new int[ltot + 1];

  // Set pointers to second child-cell (if opened) and next cell (if unopened)
  for (l=0; l<ltot; l++) {
    c2L[l] = pow(2,ltot - l);
    cNL[l] = 2*c2L[l] - 1;
  }


  // Zero tree cell variables
  for (g=0; g<gmax; g++) g2c[g] = 0;
  for (c=0; c<Ncell; c++) {
    celldata[c].c2g = 0;
    celldata[c].c1 = -1;
    celldata[c].c2 = -1;
    celldata[c].ifirst = -1;
    celldata[c].ilast = -1;
    celldata[c].N = 0;
    celldata[c].Nactive = 0;
  }
  g = 0;
  celldata[0].level = 0;


  // Loop over all cells and set all other pointers
  //-----------------------------------------------------------------------------------------------
  for (c=0; c<Ncell; c++) {
    celldata[c].id = c;
    if (celldata[c].level == ltot) {                    // If on leaf level
      celldata[c].cnext = c + 1;                        // id of next cell
      celldata[c].c2g = g;                              // Record leaf id
      g2c[g++] = c;                                     // Record inverse id
    }
    else {
      celldata[c+1].level = celldata[c].level + 1;             // Level of 1st child
      celldata[c].c1 = c + 1;
      celldata[c].c2 = c + c2L[celldata[c].level];             // id of 2nd child
      celldata[celldata[c].c2].level = celldata[c].level + 1;  // Level of 2nd child
      celldata[c].cnext = c + cNL[celldata[c].level];          // Next cell id
    }
    if (celldata[c].level == lactive) {
      celldata[c].c2g = g;                              // Record leaf id
      g2c[g++] = c;                                     // Record inverse id
    }

    // Some assert statements (for debugging)
    assert(c >= celldata[c].level);

  }
  //-----------------------------------------------------------------------------------------------


  // Free locally allocated memory
  delete[] cNL;
  delete[] c2L;

  return;
}



//=================================================================================================
//  KDTree::DivideTreeCell
/// Recursive routine to divide a tree cell into two children cells.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void KDTree<ndim,ParticleType,TreeCell>::DivideTreeCell
 (int ifirst,                          ///< Aux. id of first particle in cell
  int ilast,                           ///< Aux. id of last particle in cell
  ParticleType<ndim> *partdata,        ///< Pointer to main SPH object
  TreeCell<ndim> &cell)                ///< Cell to be divided
{
  int i;                               // Aux. child cell counter
  int j;                               // Aux. particle counter
  int k;                               // Dimension counter
  int k_divide = 0;                    // Division dimension
  FLOAT rkmax = 0.0;                   // Max. box size of all dimensions
  FLOAT rdivide;                       // Coordinate value at division


  // If cell is a leaf cell, do not divide further and set linked lists
  if (cell.level == ltot) {
    if (cell.N > 0) {
      for (j=cell.ifirst; j<cell.ilast; j++) inext[ids[j]] = ids[j+1];
      cell.ifirst = ids[cell.ifirst];
      cell.ilast = ids[cell.ilast];
    }
    else {
      cell.ifirst = -1;
      cell.ilast = -1;
    }
    StockCellProperties(cell,partdata);
    return;
  }

#if defined(MPI_PARALLEL)
  i = cell.ifirst;
  for (k=0; k< ndim; k++) {
    cell.bbmin[k]=+big_number;
    cell.bbmax[k] = -big_number;
  }
  for (i=cell.ifirst; i<=cell.ilast; i++) {
    int j = ids[i];
    for (k=0; k< ndim; k++) {
      if (cell.bbmin[k]>partdata[j].r[k]) cell.bbmin[k]=partdata[j].r[k];
      if (cell.bbmax[k]<partdata[j].r[k]) cell.bbmax[k]=partdata[j].r[k];
    }
  }
#endif

  // Determine dimension to split the cell along.
  // For now, simply split along direction of the bounding box's longest axis
  for (k=0; k<ndim; k++) {
    if (cell.bbmax[k] - cell.bbmin[k] > rkmax) {
      rkmax = cell.bbmax[k] - cell.bbmin[k];
      k_divide = k;
    }
  }
  cell.k_divide = k_divide;


  // Find median value along selected division plane and re-order array
  // so particles reside on correct side of division
  rdivide = QuickSelect(cell.ifirst,cell.ilast,cell.ifirst+cell.N/2,k_divide,partdata);

  // Set properties of first child cell
  for (k=0; k<ndim; k++) celldata[cell.c1].bbmin[k] = cell.bbmin[k];
  for (k=0; k<ndim; k++) celldata[cell.c1].bbmax[k] = cell.bbmax[k];
  for (k=0; k<ndim; k++) celldata[cell.c1].cexit[0][k] = cell.cexit[0][k];
  for (k=0; k<ndim; k++) celldata[cell.c1].cexit[1][k] = cell.cexit[1][k];
  celldata[cell.c1].bbmax[k_divide] = rdivide;
  celldata[cell.c1].cexit[1][k_divide] = cell.c2;
  celldata[cell.c1].N = cell.N/2;
  if (celldata[cell.c1].N != 0) {
    celldata[cell.c1].ifirst = ifirst;
    celldata[cell.c1].ilast = ifirst + cell.N/2 - 1;
  }

  // Set properties of second child cell
  for (k=0; k<ndim; k++) celldata[cell.c2].bbmin[k] = cell.bbmin[k];
  for (k=0; k<ndim; k++) celldata[cell.c2].bbmax[k] = cell.bbmax[k];
  for (k=0; k<ndim; k++) celldata[cell.c2].cexit[0][k] = cell.cexit[0][k];
  for (k=0; k<ndim; k++) celldata[cell.c2].cexit[1][k] = cell.cexit[1][k];
  celldata[cell.c2].bbmin[k_divide] = rdivide;
  celldata[cell.c2].cexit[0][k_divide] = cell.c1;
  celldata[cell.c2].N = cell.N - celldata[cell.c1].N;
  if (celldata[cell.c2].N != 0) {
    celldata[cell.c2].ifirst = ifirst + cell.N/2;
    celldata[cell.c2].ilast = ilast;
  }
  assert(cell.N == celldata[cell.c1].N + celldata[cell.c2].N);


  // Now divide the new child cells as a recursive function
#if defined _OPENMP
  if (pow(2,cell.level) < Nthreads) {
#pragma omp parallel default(none) private(i) shared(cell,ifirst,ilast,partdata) num_threads(2)
    {
#pragma omp for
      for (i=0; i<2; i++) {
        if (i == 0) DivideTreeCell(ifirst,ifirst+cell.N/2-1,partdata,celldata[cell.c1]);
        else if (i == 1) DivideTreeCell(ifirst+cell.N/2,ilast,partdata,celldata[cell.c2]);
      }
#pragma omp barrier
    }
  }
  else {
    for (i=0; i<2; i++) {
      if (i == 0) DivideTreeCell(ifirst,ifirst+cell.N/2-1,partdata,celldata[cell.c1]);
      else if (i == 1) DivideTreeCell(ifirst+cell.N/2,ilast,partdata,celldata[cell.c2]);
    }
  }
#else
  for (i=0; i<2; i++) {
    if (i == 0) DivideTreeCell(ifirst,ifirst+cell.N/2-1,partdata,celldata[cell.c1]);
    else if (i == 1) DivideTreeCell(ifirst+cell.N/2,ilast,partdata,celldata[cell.c2]);
  }
#endif


  // Re-set the cell first and last particles now that child cells have been
  // re-ordered by the QuickSelect algorithm
  if (celldata[cell.c1].N > 0) {
    cell.ifirst = celldata[cell.c1].ifirst;
    inext[celldata[cell.c1].ilast] = celldata[cell.c2].ifirst;
  }
  else {
    cell.ifirst = celldata[cell.c2].ifirst;
  }
  cell.ilast = celldata[cell.c2].ilast;


  assert(cell.N == celldata[cell.c1].N + celldata[cell.c2].N);
  assert(!(cell.ifirst == -1 && cell.ilast == -1));

  // Stock all cell properties once constructed
  StockCellProperties(cell,partdata);

  return;
}



//=================================================================================================
//  KDTree::QuickSelect
/// Find median and sort particles in arrays to ensure they are the correct
/// side of the division.  Uses the QuickSelect algorithm.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
FLOAT KDTree<ndim,ParticleType,TreeCell>::QuickSelectSort
 (int left,                            ///< Left-most id of particle in array
  int right,                           ///< Right-most id of particle in array
  int jpivot,                          ///< Pivot/median point
  int k,                               ///< Dimension of sort
  ParticleType<ndim> *partdata)        ///< Pointer to main particle data array
{
  int i;                               // ..
  int j;                               // ..
  int jguess;                          // ..
  int jtemp;                           // ..
  FLOAT rpivot;                        // ..
  ParticleType<ndim> temppart;         // ..


  // Place all particles left or right of chosen pivot point.
  // Iterate until correct median pivot has been identified.
  //-----------------------------------------------------------------------------------------------
  do {

    // Make a guess of pivot value
    jguess = (left + right)/2;
    rpivot = partdata[jguess].r[k];

    // ..
    temppart = partdata[jguess];
    partdata[jguess] = partdata[right];
    partdata[right] = temppart;
    partdata[jguess] = partdata[right];

    // ..
    jguess = left;

    //---------------------------------------------------------------------------------------------
    for (j=left; j<right; j++) {

      if (partdata[j].r[k] <= rpivot) {
        temppart = partdata[j];
        partdata[j] = partdata[jguess];
        partdata[jguess] = temppart;
        partdata[j] = partdata[jguess];
        jguess++;
      }

    }
    //---------------------------------------------------------------------------------------------


    // Move ?? particle from end of array to index i
    temppart = partdata[right];
    partdata[right] = partdata[jguess];
    partdata[jguess] = temppart;
    partdata[right] = partdata[jguess];


    // jguess is lower than jpivot.
    // Only need to search between jguess+1 and right
    if (jguess < jpivot) left = jguess + 1;

    // jguess is higher than jpivot.
    // Only need to search between left and jguess-1
    else if (jguess > jpivot) right = jguess - 1;


  } while (jguess != jpivot);
  //-----------------------------------------------------------------------------------------------


  return rpivot;
}



//=================================================================================================
//  KDTree::QuickSelect
/// Find median and sort particles in arrays to ensure they are the correct
/// side of the division.  Uses the QuickSelect algorithm.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
FLOAT KDTree<ndim,ParticleType,TreeCell>::QuickSelect
(int left,                          ///< Left-most id of particle in array
 int right,                         ///< Right-most id of particle in array
 int jpivot,                        ///< Pivot/median point
 int k,                             ///< Dimension of sort
 ParticleType<ndim> *partdata)      ///< Pointer to main SPH object
{
  int i;                            // ..
  int j;                            // ..
  int jguess;                       // ..
  int jtemp;                        // ..
  FLOAT rpivot;                     // ..


  // Place all particles left or right of chosen pivot point.
  // Iterate until correct median pivot has been identified.
  //-----------------------------------------------------------------------------------------------
  do {

    // Make a guess of pivot value
    jguess = (left + right)/2;
    rpivot = partdata[ids[jguess]].r[k];

    // ..
    jtemp = ids[jguess];
    ids[jguess] = ids[right];
    ids[right] = jtemp;

    // ..
    jguess = left;

    //---------------------------------------------------------------------------------------------
    for (j=left; j<right; j++) {
      //assert(j < Ntot);
      if (partdata[ids[j]].r[k] <= rpivot) {
        jtemp = ids[j];
        ids[j] = ids[jguess];
        ids[jguess] = jtemp;
        jguess++;
      }

    }
    //---------------------------------------------------------------------------------------------


    // Move ?? particle from end of array to index i
    jtemp = ids[right];
    ids[right] = ids[jguess];
    ids[jguess] = jtemp;

    //assert(left < Ntot);
    //assert(right < Ntot);
    //assert(jguess < Ntot);
    //assert(jpivot < Ntot);


    // jguess is lower than jpivot.
    // Only need to search between jguess+1 and right
    if (jguess < jpivot) left = jguess + 1;

    // jguess is higher than jpivot.
    // Only need to search between left and jguess-1
    else if (jguess > jpivot) right = jguess - 1;

  } while (jguess != jpivot);
  //-----------------------------------------------------------------------------------------------


  return rpivot;
}



//=================================================================================================
//  KDTree::StockTree
/// Stock given tree cell in KD-tree.  If cell is not a leaf-cell, recursively
/// calls itself for its two child cells.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void KDTree<ndim,ParticleType,TreeCell>::StockTree
 (TreeCell<ndim> &cell,                ///< Reference to cell to be stocked
  ParticleType<ndim> *partdata)        ///< SPH particle data array
{
  int i;                               // Aux. child cell counter

  // If cell is not leaf, stock child cells
  if (cell.level != ltot) {
#if defined _OPENMP
    if (pow(2,cell.level) < Nthreads) {
#pragma omp parallel for default(none) private(i) shared(cell,partdata) num_threads(2)
      for (i=0; i<2; i++) {
        if (i == 0) StockTree(celldata[cell.c1],partdata);
        else if (i == 1) StockTree(celldata[cell.c2],partdata);
      }
    }
    else {
      for (i=0; i<2; i++) {
        if (i == 0) StockTree(celldata[cell.c1],partdata);
        else if (i == 1) StockTree(celldata[cell.c2],partdata);
      }
    }
#else
    for (i=0; i<2; i++) {
      if (i == 0) StockTree(celldata[cell.c1],partdata);
      else if (i == 1) StockTree(celldata[cell.c2],partdata);
    }
#endif
  }

  // Stock node once all children are stocked
  StockCellProperties(cell,partdata);

  return;
}



//=================================================================================================
//  KDTree::StockCellProperties
/// Calculate the physical properties (e.g. total mass, centre-of-mass,
/// opening-distance, etc..) of all cells in the tree.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void KDTree<ndim,ParticleType,TreeCell>::StockCellProperties
 (TreeCell<ndim> &cell,                ///< Reference to current tree cell
  ParticleType<ndim> *partdata)        ///< Particle data array
{
  int cc,ccc;                          // Cell counters
  int i;                               // Particle counter
  int iaux;                            // Aux. particle i.d. variable
  int j;                               // ..
  int k;                               // Dimension counter
  FLOAT dr[ndim];                      // Relative position vector
  FLOAT drsqd;                         // Distance squared
  FLOAT mi;                            // Mass of particle i
  FLOAT p = 0.0;                       // ..
  FLOAT lambda = 0.0;                  // ..
  TreeCell<ndim> &child1 = celldata[cell.c1];
  TreeCell<ndim> &child2 = celldata[cell.c2];


  // Zero all summation variables for all cells
  cell.Nactive  = 0;
  cell.N        = 0;
  cell.m        = 0.0;
  cell.hmax     = 0.0;
  cell.rmax     = 0.0;
  cell.dhmaxdt  = 0.0;
  cell.drmaxdt  = 0.0;
  cell.mac      = 0.0;
  cell.cdistsqd = big_number;
  for (k=0; k<5; k++) cell.q[k] = 0.0;
  for (k=0; k<ndim; k++) cell.r[k] = 0.0;
  for (k=0; k<ndim; k++) cell.v[k] = 0.0;
  for (k=0; k<ndim; k++) cell.rcell[k] = 0.0;
  for (k=0; k<ndim; k++) cell.bbmin[k] = big_number;
  for (k=0; k<ndim; k++) cell.bbmax[k] = -big_number;
  for (k=0; k<ndim; k++) cell.hboxmin[k] = big_number;
  for (k=0; k<ndim; k++) cell.hboxmax[k] = -big_number;


  // If this is a leaf cell, sum over all particles
  //-----------------------------------------------------------------------------------------------
  if (cell.level == ltot) {

    // First, check if any particles have been accreted and remove them
    // from the linked list.  If cell no longer contains any live particles,
    // then set N = 0 to ensure cell is not included in future tree-walks.
    i = cell.ifirst;
    cell.ifirst = -1;
    iaux = -1;
    while (i != -1) {
      if (partdata[i].itype != dead) {
        if (iaux == -1) cell.ifirst = i;
        else inext[iaux] = i;
        iaux = i;
      }
      if (i == cell.ilast) break;
      i = inext[i];
    };
    cell.ilast = iaux;


    // Loop over all particles in cell summing their contributions
    i = cell.ifirst;
    while (i != -1) {
      if (partdata[i].itype != dead) {
        cell.N++;
        if (partdata[i].active) cell.Nactive++;
        cell.hmax = max(cell.hmax,partdata[i].h);
        cell.m += partdata[i].m;
        for (k=0; k<ndim; k++) cell.r[k] += partdata[i].m*partdata[i].r[k];
        for (k=0; k<ndim; k++) cell.v[k] += partdata[i].m*partdata[i].v[k];
        for (k=0; k<ndim; k++) {
          if (partdata[i].r[k] < cell.bbmin[k]) cell.bbmin[k] = partdata[i].r[k];
          if (partdata[i].r[k] > cell.bbmax[k]) cell.bbmax[k] = partdata[i].r[k];
          if (partdata[i].r[k] - kernrange*partdata[i].h < cell.hboxmin[k])
            cell.hboxmin[k] = partdata[i].r[k] - kernrange*partdata[i].h;
          if (partdata[i].r[k] + kernrange*partdata[i].h > cell.hboxmax[k])
            cell.hboxmax[k] = partdata[i].r[k] + kernrange*partdata[i].h;
        }
      }
      if (i == cell.ilast) break;
      i = inext[i];
    };

    // Normalise all cell values
    if (cell.N > 0) {
      for (k=0; k<ndim; k++) cell.r[k] /= cell.m;
      for (k=0; k<ndim; k++) cell.v[k] /= cell.m;
      for (k=0; k<ndim; k++) cell.rcell[k] = 0.5*(cell.bbmin[k] + cell.bbmax[k]);
      for (k=0; k<ndim; k++) dr[k] = 0.5*(cell.bbmax[k] - cell.bbmin[k]);
      cell.cdistsqd = max(DotProduct(dr,dr,ndim),cell.hmax*cell.hmax)/thetamaxsqd;
      cell.rmax = sqrt(DotProduct(dr,dr,ndim));
    }

    // Compute quadrupole moment terms if selected
    if (multipole == "quadrupole") {
      i = cell.ifirst;

      while (i != -1) {
        if (partdata[i].itype != dead) {
          mi = partdata[i].m;
          for (k=0; k<ndim; k++) dr[k] = partdata[i].r[k] - cell.r[k];
          drsqd = DotProduct(dr,dr,ndim);
          if (ndim == 3) {
            cell.q[0] += mi*(3.0*dr[0]*dr[0] - drsqd);
            cell.q[1] += mi*3.0*dr[0]*dr[1];
            cell.q[2] += mi*(3.0*dr[1]*dr[1] - drsqd);
            cell.q[3] += mi*3.0*dr[2]*dr[0];
            cell.q[4] += mi*3.0*dr[2]*dr[1];
          }
          else if (ndim == 2) {
            cell.q[0] += mi*(3.0*dr[0]*dr[0] - drsqd);
            cell.q[1] += mi*3.0*dr[0]*dr[1];
            cell.q[2] += mi*(3.0*dr[1]*dr[1] - drsqd);
          }
        }
        if (i == cell.ilast) break;
        i = inext[i];
      }
    }

  }
  // For non-leaf cells, sum together two children cells
  //-----------------------------------------------------------------------------------------------
  else {

    if (child1.N > 0) {
      for (k=0; k<ndim; k++) cell.bbmin[k] = min(child1.bbmin[k],cell.bbmin[k]);
      for (k=0; k<ndim; k++) cell.bbmax[k] = max(child1.bbmax[k],cell.bbmax[k]);
      for (k=0; k<ndim; k++) cell.hboxmin[k] = min(child1.hboxmin[k],cell.hboxmin[k]);
      for (k=0; k<ndim; k++) cell.hboxmax[k] = max(child1.hboxmax[k],cell.hboxmax[k]);
      cell.hmax = max(cell.hmax,child1.hmax);
    }
    if (child2.N > 0) {
      for (k=0; k<ndim; k++) cell.bbmin[k] = min(child2.bbmin[k],cell.bbmin[k]);
      for (k=0; k<ndim; k++) cell.bbmax[k] = max(child2.bbmax[k],cell.bbmax[k]);
      for (k=0; k<ndim; k++) cell.hboxmin[k] = min(child2.hboxmin[k],cell.hboxmin[k]);
      for (k=0; k<ndim; k++) cell.hboxmax[k] = max(child2.hboxmax[k],cell.hboxmax[k]);
      cell.hmax = max(cell.hmax,child2.hmax);
    }

    cell.N = child1.N + child2.N;
    cell.Nactive = child1.Nactive + child2.Nactive;
    if (cell.N > 0) {
      cell.m = child1.m + child2.m;
      for (k=0; k<ndim; k++) cell.r[k] = (child1.m*child1.r[k] + child2.m*child2.r[k])/cell.m;
      for (k=0; k<ndim; k++) cell.v[k] = (child1.m*child1.v[k] + child2.m*child2.v[k])/cell.m;
      for (k=0; k<ndim; k++) cell.rcell[k] = 0.5*(cell.bbmin[k] + cell.bbmax[k]);
      for (k=0; k<ndim; k++) dr[k] = 0.5*(cell.bbmax[k] - cell.bbmin[k]);
      cell.cdistsqd = max(DotProduct(dr,dr,ndim),cell.hmax*cell.hmax)/thetamaxsqd;
      cell.rmax = sqrt(DotProduct(dr,dr,ndim));
#ifdef MPI_PARALLEL
      cell.worktot = child1.worktot + child2.worktot;
#endif
    }

    // Now add individual quadrupole moment terms
    if (multipole == "quadrupole" && child1.N > 0) {
      mi = child1.m;
      for (k=0; k<ndim; k++) dr[k] = child1.r[k] - cell.r[k];
      drsqd = DotProduct(dr,dr,ndim);
      if (ndim == 3) {
        cell.q[0] += mi*(3.0*dr[0]*dr[0] - drsqd);
        cell.q[1] += mi*3.0*dr[0]*dr[1];
        cell.q[2] += mi*(3.0*dr[1]*dr[1] - drsqd);
        cell.q[3] += mi*3.0*dr[2]*dr[0];
        cell.q[4] += mi*3.0*dr[2]*dr[1];
      }
      else if (ndim == 2) {
        cell.q[0] += mi*(3.0*dr[0]*dr[0] - drsqd);
        cell.q[1] += mi*3.0*dr[0]*dr[1];
        cell.q[2] += mi*(3.0*dr[1]*dr[1] - drsqd);
      }
    }

    if (multipole == "quadrupole" && child2.N > 0) {
      mi = child2.m;
      for (k=0; k<ndim; k++) dr[k] = child2.r[k] - cell.r[k];
      drsqd = DotProduct(dr,dr,ndim);
      if (ndim == 3) {
        cell.q[0] += mi*(3.0*dr[0]*dr[0] - drsqd);
        cell.q[1] += mi*3.0*dr[0]*dr[1];
        cell.q[2] += mi*(3.0*dr[1]*dr[1] - drsqd);
        cell.q[3] += mi*3.0*dr[2]*dr[0];
        cell.q[4] += mi*3.0*dr[2]*dr[1];
      }
      else if (ndim == 2) {
        cell.q[0] += mi*(3.0*dr[0]*dr[0] - drsqd);
        cell.q[1] += mi*3.0*dr[0]*dr[1];
        cell.q[2] += mi*(3.0*dr[1]*dr[1] - drsqd);
      }
    }

  }
  //-----------------------------------------------------------------------------------------------


  // Calculate eigenvalue MAC criteria
  if (gravity_mac == "eigenmac") {
    if (ndim == 3)
      p = cell.q[0]*cell.q[2] - (cell.q[0] + cell.q[2])*(cell.q[0] + cell.q[2]) -
        cell.q[1]*cell.q[1] - cell.q[3]*cell.q[3] - cell.q[4]*cell.q[4];
    if (p >= 0.0) cell.mac = 0.0;
    else {
      lambda = 2.0*sqrt(-p/3.0);
      cell.mac = pow(0.5*lambda/macerror,0.66666666666666);
    }
  }
  else {
    cell.mac = 0.0;
  }


  return;
}



//=================================================================================================
//  KDTree::ExtrapolateCellProperties
/// Extrapolate important physical properties of all cells in the tree.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void KDTree<ndim,ParticleType,TreeCell>::ExtrapolateCellProperties
 (FLOAT dt)                            ///< Smallest timestep size
{
  int c;                               // Cell counter
  int k;                               // Dimension counter

  debug2("[KDTree::ExtrapolateCellProperties]");


  // ...
  //-----------------------------------------------------------------------------------------------
  for (c=0; c<Ncell; c++) {

    for (k=0; k<ndim; k++) celldata[c].r[k] += celldata[c].v[k]*dt;
    for (k=0; k<ndim; k++) celldata[c].rcell[k] += celldata[c].v[k]*dt;
    for (k=0; k<ndim; k++) celldata[c].bbmin[k] += celldata[c].v[k]*dt;
    for (k=0; k<ndim; k++) celldata[c].bbmax[k] += celldata[c].v[k]*dt;
    for (k=0; k<ndim; k++) celldata[c].hboxmin[k] += celldata[c].v[k]*dt;
    for (k=0; k<ndim; k++) celldata[c].hboxmax[k] += celldata[c].v[k]*dt;
    //celldata[c].rmax += celldata[c].drmaxdt*dt;
    //celldata[c].hmax += celldata[c].dhmaxdt*dt;

  }
  //-----------------------------------------------------------------------------------------------

  return;
}



//=================================================================================================
//  KDTree::UpdateHmaxValues
/// Calculate the physical properties (e.g. total mass, centre-of-mass,
/// opening-distance, etc..) of all cells in the tree.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void KDTree<ndim,ParticleType,TreeCell>::UpdateHmaxValues
 (TreeCell<ndim> &cell,                ///< KD-tree cell
  ParticleType<ndim> *partdata)        ///< SPH particle data array
{
  int c,cc,ccc;                        // Cell counters
  int i;                               // Particle counter
  int k;                               // Dimension counter

  // If cell is not leaf, stock child cells
  if (cell.level != ltot) {
#if defined _OPENMP
    if (pow(2,cell.level) < Nthreads) {
#pragma omp parallel for default(none) private(i) shared(cell,partdata) num_threads(2)
      for (i=0; i<2; i++) {
        if (i == 0) UpdateHmaxValues(celldata[cell.c1],partdata);
        else if (i == 1) UpdateHmaxValues(celldata[cell.c2],partdata);
      }
    }
    else {
      for (i=0; i<2; i++) {
        if (i == 0) UpdateHmaxValues(celldata[cell.c1],partdata);
        else if (i == 1) UpdateHmaxValues(celldata[cell.c2],partdata);
      }
    }
#else
    for (i=0; i<2; i++) {
      if (i == 0) UpdateHmaxValues(celldata[cell.c1],partdata);
      else if (i == 1) UpdateHmaxValues(celldata[cell.c2],partdata);
    }
#endif
  }


  // Zero all summation variables for all cells
  cell.hmax = 0.0;
  for (k=0; k<ndim; k++) cell.hboxmin[k] = big_number;
  for (k=0; k<ndim; k++) cell.hboxmax[k] = -big_number;


  // If this is a leaf cell, sum over all particles
  //-----------------------------------------------------------------------------------------------
  if (cell.level == ltot) {
    i = cell.ifirst;

    // Loop over all particles in cell summing their contributions
    while (i != -1) {
      cell.hmax = max(cell.hmax,partdata[i].h);
      for (k=0; k<ndim; k++) {
        if (partdata[i].r[k] - kernrange*partdata[i].h < cell.hboxmin[k]) {
          cell.hboxmin[k] = partdata[i].r[k] - kernrange*partdata[i].h;
        }
        if (partdata[i].r[k] + kernrange*partdata[i].h > cell.hboxmax[k]) {
          cell.hboxmax[k] = partdata[i].r[k] + kernrange*partdata[i].h;
        }
      }
      if (i == cell.ilast) break;
      i = inext[i];
    };

  }
  // For non-leaf cells, sum together two children cells
  //-----------------------------------------------------------------------------------------------
  else {
    cc = cell.c1;
    ccc = cell.c2;
    if (celldata[cc].N > 0) {
      cell.hmax = max(cell.hmax,celldata[cc].hmax);
      for (k=0; k<ndim; k++) cell.hboxmin[k] = min(celldata[cc].hboxmin[k],cell.hboxmin[k]);
      for (k=0; k<ndim; k++) cell.hboxmax[k] = max(celldata[cc].hboxmax[k],cell.hboxmax[k]);
    }
    if (celldata[ccc].N > 0) {
      cell.hmax = max(cell.hmax,celldata[ccc].hmax);
      for (k=0; k<ndim; k++) cell.hboxmin[k] = min(celldata[ccc].hboxmin[k],cell.hboxmin[k]);
      for (k=0; k<ndim; k++) cell.hboxmax[k] = max(celldata[ccc].hboxmax[k],cell.hboxmax[k]);
    }

  }
  //-----------------------------------------------------------------------------------------------

  return;
}



//=================================================================================================
//  KDTree::UpdateActiveParticleCounters
/// Loop through all leaf cells in KD-tree and update all active particle counters.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void KDTree<ndim,ParticleType,TreeCell>::UpdateActiveParticleCounters
 (ParticleType<ndim> *partdata)        ///< [in] Main particle array
{
  int c;                               // Cell counter
  int i;                               // SPH particle index
  int ilast;                           // Last particle in linked list

  debug2("[KDTree::UpdateActiveParticleCounters]");
  //timing->StartTimingSection("TREE_UPDATE_COUNTERS",2);


  // Loop through all grid cells in turn
  //-----------------------------------------------------------------------------------------------
#pragma omp parallel for default(none) private(c,i,ilast) shared(partdata)
  for (c=0; c<Ncell; c++) {
    celldata[c].Nactive = 0;

    //if (celldata[c].level != ltot) continue;
    if (celldata[c].level != lactive) continue;
    i = celldata[c].ifirst;
    ilast = celldata[c].ilast;

    // Else walk through linked list to obtain list and number of active ptcls.
    while (i != -1) {
      if (i < Ntot && partdata[i].active && partdata[i].itype != dead) celldata[c].Nactive++;
      if (i == ilast) break;
      i = inext[i];
    };

  }
  //-----------------------------------------------------------------------------------------------

  //timing->EndTimingSection("TREE_UPDATE_COUNTERS");

  return;
}



//=================================================================================================
//  KDTree::ComputeActiveParticleList
/// Returns the number (Nactive) and list of ids (activelist) of all active
/// SPH particles in the given cell.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
int KDTree<ndim,ParticleType,TreeCell>::ComputeActiveParticleList
 (TreeCell<ndim> &cell,                ///< [in] Pointer to cell
  ParticleType<ndim> *partdata,        ///< [in] Pointer to particle data array
  int *activelist)                     ///< [out] List of active particles in cell
{
  int i = cell.ifirst;                 // Local particle id (set to first ptcl id)
  int ilast = cell.ilast;              // i.d. of last particle in cell c
  int Nactive = 0;                     // No. of active particles in cell

  // Walk through linked list to obtain list and number of active ptcls.
  while (i != -1) {
    if (i < Ntot && partdata[i].active && partdata[i].itype != dead) activelist[Nactive++] = i;
    if (i == ilast) break;
    i = inext[i];
  };

  return Nactive;
}



//=================================================================================================
//  KDTree::BoxOverlap
/// Check if two bounding boxes overlap.  If yes, then returns true.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
bool KDTree<ndim,ParticleType,TreeCell>::BoxOverlap
 (const FLOAT box1min[ndim],           ///< Minimum extent of box 1
  const FLOAT box1max[ndim],           ///< Maximum extent of box 1
  const FLOAT box2min[ndim],           ///< Minimum extent of box 2
  const FLOAT box2max[ndim])           ///< Maximum extent of box 2
{
  if (ndim == 1) {
    if (box1min[0] > box2max[0]) return false;
    if (box2min[0] > box1max[0]) return false;
    return true;
  }
  else if (ndim == 2) {
    if (box1min[0] > box2max[0]) return false;
    if (box2min[0] > box1max[0]) return false;
    if (box1min[1] > box2max[1]) return false;
    if (box2min[1] > box1max[1]) return false;
    return true;
  }
  else if (ndim == 3) {
    if (box1min[0] > box2max[0]) return false;
    if (box2min[0] > box1max[0]) return false;
    if (box1min[1] > box2max[1]) return false;
    if (box2min[1] > box1max[1]) return false;
    if (box1min[2] > box2max[2]) return false;
    if (box2min[2] > box1max[2]) return false;
    return true;
  }

}



//=================================================================================================
//  KDTree::ComputeGatherNeighbourList
/// Computes and returns number of neighbour, 'Nneib', and the list
/// of neighbour ids, 'neiblist', for all particles inside cell 'c'.
/// Includes all particles in the selected cell, plus all particles
/// contained in adjacent cells (including diagonal cells).
/// Wrapper around the true implementation inside KDTree
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
int KDTree<ndim,ParticleType,TreeCell>::ComputeGatherNeighbourList
 (const ParticleType<ndim> *partdata,  ///< [in] Particle data array
  const FLOAT rp[ndim],                ///< [in] Search position
  const FLOAT rsearch,                 ///< [in] Maximum smoothing length
  const int Nneibmax,                  ///< [in] Max. no. of neighbours
  int *neiblist)                       ///< [out] List of neighbour i.d.s
{
  int cc;                              // Cell counter
  int i;                               // Particle id
  int j;                               // Aux. particle counter
  int k;                               // Neighbour counter
  int Nneib = 0;                       // Neighbour counter
  int Ntemp = 0;                       // Temporary neighbour counter
  FLOAT dr[ndim];                      // Relative position vector
  FLOAT drsqd;                         // Distance squared
  FLOAT rsearchsqd;                    // Search radius squared

  // Start with root cell and walk through entire tree
  cc = 0;
  rsearchsqd = rsearch*rsearch;

  //===============================================================================================
  while (cc < Ncell) {

    for (k=0; k<ndim; k++) dr[k] = celldata[cc].rcell[k] - rp[k];
    drsqd = DotProduct(dr,dr,ndim);


    // Check if bounding boxes overlap with each other
    //---------------------------------------------------------------------------------------------
    if (drsqd < (rsearch + celldata[cc].rmax)*(rsearch + celldata[cc].rmax)) {

      // If not a leaf-cell, then open cell to first child cell
      if (celldata[cc].level != ltot)
        cc++;

      else if (celldata[cc].N == 0)
        cc = celldata[cc].cnext;

      // If leaf-cell, add particles to list
      else if (celldata[cc].level == ltot && Nneib + Nleafmax < Nneibmax) {
        i = celldata[cc].ifirst;
        while (i != -1) {
          for (k=0; k<ndim; k++) dr[k] = partdata[i].r[k] - rp[k];
          drsqd = DotProduct(dr,dr,ndim);
          if (drsqd < rsearchsqd && partdata[i].itype != dead) neiblist[Nneib++] = i;
          if (i == celldata[cc].ilast) break;
          i = inext[i];
        };
        cc = celldata[cc].cnext;
      }

      // If leaf-cell, but we've run out of memory, return with error-code (-1)
      else if (celldata[cc].level == ltot && Nneib + Nleafmax >= Nneibmax)
        return -1;

    }

    // If not in range, then open next cell
    //---------------------------------------------------------------------------------------------
    else
      cc = celldata[cc].cnext;

  };
  //===============================================================================================


  return Nneib;
}



//=================================================================================================
//  KDTree::ComputeGatherNeighbourList
/// Computes and returns number of neighbour, 'Nneib', and the list
/// of neighbour ids, 'neiblist', for all particles inside cell 'c'.
/// Includes all particles in the selected cell, plus all particles
/// contained in adjacent cells (including diagonal cells).
/// Wrapper around the true implementation inside KDTree
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
int KDTree<ndim,ParticleType,TreeCell>::ComputeGatherNeighbourList
 (const TreeCell<ndim> &cell,          ///< [in] Pointer to current cell
  const ParticleType<ndim> *partdata,  ///< [in] Particle data array
  const FLOAT hmax,                    ///< [in] Maximum smoothing length
  const int Nneibmax,                  ///< [in] Max. no. of neighbours
  int &Nneib,                          ///< [inout] No. of neighbours
  int *neiblist)                       ///< [out] List of neighbour i.d.s
{
  int cc;                              // Cell counter
  int i;                               // Particle id
  int j;                               // Aux. particle counter
  int k;                               // Neighbour counter
  int Ntemp = Nneib;                   // Temporary neighbour counter
  FLOAT dr[ndim];                      // Relative position vector
  FLOAT drsqd;                         // Distance squared
  FLOAT gatherboxmin[ndim];            // Minimum gather neighbour box
  FLOAT gatherboxmax[ndim];            // Maximum gather neighbour box
  FLOAT rc[ndim];                      // Position of cell

  // Exit immediately if we have overflowed the neighbour list buffer
  if (Nneib == -1) return -1;

  const FLOAT hrangemaxsqd = pow(cell.rmax + kernrange*hmax,2);
  for (k=0; k<ndim; k++) rc[k] = cell.rcell[k];
  for (k=0; k<ndim; k++) gatherboxmin[k] = cell.bbmin[k] - kernrange*hmax;
  for (k=0; k<ndim; k++) gatherboxmax[k] = cell.bbmax[k] + kernrange*hmax;


  // Start with root cell and walk through entire tree
  cc = 0;

  //===============================================================================================
  while (cc < Ncell) {

    // Check if bounding boxes overlap with each other
    //---------------------------------------------------------------------------------------------
    if (BoxOverlap(gatherboxmin,gatherboxmax,celldata[cc].bbmin,celldata[cc].bbmax)) {

      // If not a leaf-cell, then open cell to first child cell
      if (celldata[cc].level != ltot)
        cc++;

      else if (celldata[cc].N == 0)
        cc = celldata[cc].cnext;

      // If leaf-cell, add particles to list
      else if (celldata[cc].level == ltot && Ntemp + Nleafmax < Nneibmax) {
        i = celldata[cc].ifirst;
        while (i != -1) {
          neiblist[Ntemp++] = i;
          if (i == celldata[cc].ilast) break;
          i = inext[i];
        };
        cc = celldata[cc].cnext;
      }

      // If leaf-cell, but we've run out of memory, return with error-code (-1)
      else if (celldata[cc].level == ltot && Ntemp + Nleafmax >= Nneibmax)
    	return -1;

    }

    // If not in range, then open next cell
    //---------------------------------------------------------------------------------------------
    else
      cc = celldata[cc].cnext;

  };
  //===============================================================================================


  // Now, trim the list to remove particles that are definitely not neighbours
  for (j=Nneib; j<Ntemp; j++) {
    i = neiblist[j];
    if (partdata[i].itype == dead) continue;
    for (k=0; k<ndim; k++) dr[k] = partdata[i].r[k] - rc[k];
    drsqd = DotProduct(dr,dr,ndim);
    if (drsqd < hrangemaxsqd) neiblist[Nneib++] = i;
  }


  return Nneib;
}



//=================================================================================================
//  KDTree::ComputeNeighbourList
/// Computes and returns number of neighbour, 'Nneib', and the list
/// of neighbour ids, 'neiblist', for all particles inside cell 'c'.
/// Includes all particles in the selected cell, plus all particles
/// contained in adjacent cells (including diagonal cells).
/// Wrapper around the true implementation inside KDTree.
/// If allocated memory array containing neighbour ids (neiblist) overflows,
/// return with error code (-1) in order to reallocate more memory.
//================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
int KDTree<ndim,ParticleType,TreeCell>::ComputeNeighbourList
 (const TreeCell<ndim> &cell,          ///< [in] Cell pointer
  const ParticleType<ndim> *partdata,  ///< [in] Particle data array
  const int Nneibmax,                  ///< [in] Max. no. of neighbours
  int &Nneib,                          ///< [inout] No. of neighbours
  int *neiblist,                       ///< [out] List of neighbour i.d.s
  ParticleType<ndim> *neibpart)        ///< [out] Array of local copies of neighbours
{
  int cc;                              // Cell counter
  int i;                               // Particle id
  int j;                               // Aux. particle counter
  int k;                               // Neighbour counter
  int Ntemp = Nneib;                   // Temp. no. of neighbouts
  FLOAT dr[ndim];                      // Relative position vector
  FLOAT drsqd;                         // Distance squared
  FLOAT rc[ndim];                      // Position of cell

  // Local (const) copies of cell properties
  const FLOAT hrangemaxsqd = pow(cell.rmax + kernrange*cell.hmax,2);
  const FLOAT rmax = cell.rmax;
  for (k=0; k<ndim; k++) rc[k] = cell.rcell[k];

  // Exit immediately if we have overflowed the neighbour list buffer
  if (Nneib == -1) return -1;

  // Start with root cell and walk through entire tree
  cc = 0;


  //===============================================================================================
  while (cc < Ncell) {

    for (k=0; k<ndim; k++) dr[k] = celldata[cc].rcell[k] - rc[k];
    drsqd = DotProduct(dr,dr,ndim);

    // Check if bounding boxes overlap with each other
    //---------------------------------------------------------------------------------------------
    if (BoxOverlap(cell.bbmin,cell.bbmax,celldata[cc].hboxmin,celldata[cc].hboxmax) ||
        BoxOverlap(cell.hboxmin,cell.hboxmax,celldata[cc].bbmin,celldata[cc].bbmax)) {

      // If not a leaf-cell, then open cell to first child cell
      if (celldata[cc].level != ltot)
        cc++;

      else if (celldata[cc].N == 0)
        cc = celldata[cc].cnext;

      // If leaf-cell, add particles to list
      else if (celldata[cc].level == ltot && Ntemp + Nleafmax < Nneibmax) {
        i = celldata[cc].ifirst;
        while (i != -1) {
          //neibpart[Ntemp] = partdata[i];
          neiblist[Ntemp++] = i;
          if (i == celldata[cc].ilast) break;
          i = inext[i];
        };
        cc = celldata[cc].cnext;
      }

      // If leaf-cell, but we've run out of memory, return with error-code (-1)
      else if (celldata[cc].level == ltot && Ntemp + Nleafmax >= Nneibmax)
        return -1;

    }

    // If not in range, then open next cell
    //---------------------------------------------------------------------------------------------
    else
      cc = celldata[cc].cnext;

  };
  //===============================================================================================


  // Now, trim the list to remove particles that are definitely not neighbours
  for (j=Nneib; j<Ntemp; j++) {
    i = neiblist[j];
    if (partdata[i].itype == dead) continue;

    for (k=0; k<ndim; k++) dr[k] = partdata[i].r[k] - rc[k];
    drsqd = DotProduct(dr,dr,ndim);
    if (drsqd < hrangemaxsqd || drsqd <
        (rmax + kernrange*partdata[i].h)*(rmax + kernrange*partdata[i].h)) {
      neibpart[Nneib] = partdata[i];
      neiblist[Nneib] = i;
      Nneib++;
    }
  }


  return Nneib;
}



//=================================================================================================
//  KDTree::ComputeGravityInteractionList
/// Computes and returns number of SPH neighbours (Nneib), direct sum particles
/// (Ndirect) and number of cells (Ngravcell), including lists of ids, from
/// the gravity tree walk for active particles inside cell c.
/// Currently defaults to the geometric opening criteria.
/// If any of the interactions list arrays (neiblist,directlist,gravcelllist)
/// overflow, return with error code (-1) to reallocate more memory.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
int KDTree<ndim,ParticleType,TreeCell>::ComputeGravityInteractionList
 (const TreeCell<ndim> &cell,          ///< [in] Pointer to cell
  const ParticleType<ndim> *partdata,  ///< [in] Particle data array
  const FLOAT macfactor,               ///< [in] Gravity MAC particle factor
  const int Nneibmax,                  ///< [in] Max. no. of SPH neighbours
  const int Ngravcellmax,              ///< [in] Max. no. of cell interactions
  int &Nneib,                          ///< [out] Total no. of neighbours
  int &Nsphneib,                       ///< [out] No. of SPH neighbours
  int &Ndirect,                        ///< [out] No. of direct-sum neighbours
  int &Ngravcell,                      ///< [out] No. of cell interactions
  int *neiblist,                       ///< [out] List of all particle ids
  int *sphneiblist,                    ///< [out] List of SPH neibpart ids
  int *directlist,                     ///< [out] List of direct-sum neibpart ids
  TreeCell<ndim> *gravcell,            ///< [out] Array of local copies of tree cells
  ParticleType<ndim> *neibpart)        ///< [out] Array of local copies of neighbour particles
{
  int cc;                              // Cell counter
  int i;                               // Particle id
  int j;                               // Aux. particle counter
  int k;                               // Neighbour counter
  int Nsphneibtemp = 0;                // Aux. counter
  FLOAT dr[ndim];                      // Relative position vector
  FLOAT drsqd;                         // Distance squared
  FLOAT rc[ndim];                      // Position of cell

  // Make local copies of important cell properties
  const FLOAT hrangemaxsqd = pow(cell.rmax + kernrange*cell.hmax,2);
  const FLOAT rmax = cell.rmax;
  for (k=0; k<ndim; k++) rc[k] = cell.rcell[k];

  // Start with root cell and walk through entire tree
  cc        = 0;
  Nneib     = 0;
  Nsphneib  = 0;
  Ndirect   = 0;
  Ngravcell = 0;


  // Walk through all cells in tree to determine particle and cell interaction lists
  //===============================================================================================
  while (cc < Ncell) {

    for (k=0; k<ndim; k++) dr[k] = celldata[cc].rcell[k] - rc[k];
    drsqd = DotProduct(dr,dr,ndim);


    // Check if bounding boxes overlap with each other
    //---------------------------------------------------------------------------------------------
    if (BoxOverlap(cell.bbmin,cell.bbmax,celldata[cc].hboxmin,celldata[cc].hboxmax) ||
        BoxOverlap(cell.hboxmin,cell.hboxmax,celldata[cc].bbmin,celldata[cc].bbmax)) {

      // If not a leaf-cell, then open cell to first child cell
      if (celldata[cc].level != ltot)
        cc++;

      else if (celldata[cc].N == 0)
        cc = celldata[cc].cnext;

      // If leaf-cell, add particles to list
      else if (celldata[cc].level == ltot && Nneib + Nleafmax <= Nneibmax) {
        i = celldata[cc].ifirst;
        while (i != -1) {
          sphneiblist[Nsphneib++] = Nneib;
          neiblist[Nneib] = i;
          neibpart[Nneib++] = partdata[i];
          if (i == celldata[cc].ilast) break;
          i = inext[i];
        };
        cc = celldata[cc].cnext;
      }

      // If leaf-cell, but we've run out of memory, return with error-code (-1)
      else if (celldata[cc].level == ltot && Nneib + Nleafmax > Nneibmax)
        return -1;

    }

    // Check if cell is far enough away to use the COM approximation
    //---------------------------------------------------------------------------------------------
    else if (drsqd > celldata[cc].cdistsqd && drsqd > celldata[cc].mac*macfactor
              && celldata[cc].N > 0) {

      // If cell is a leaf-cell with only one particle, more efficient to
      // compute the gravitational contribution from the particle than the cell
      if (celldata[cc].level == ltot && celldata[cc].N == 1 && Ndirect + Nneib < Nneibmax) {
        i = celldata[cc].ifirst;
        directlist[Ndirect++] = Nneib;
        neiblist[Nneib] = i;
        neibpart[Nneib++] = partdata[i];
      }
      else if (Ngravcell < Ngravcellmax) {
        gravcell[Ngravcell++] = celldata[cc];
      }
      else {
        return -2;
      }
      cc = celldata[cc].cnext;

    }

    // If cell is too close, open cell to interogate children cells.
    // If cell is too close and a leaf cell, then add particles to direct list.
    //---------------------------------------------------------------------------------------------
    else if (!(drsqd > celldata[cc].cdistsqd && drsqd > celldata[cc].mac*macfactor )
	     && celldata[cc].N > 0) {

      // If not a leaf-cell, then open cell to first child cell
      if (celldata[cc].level != ltot)
         cc++;

      // If leaf-cell, add particles to list
      else if (celldata[cc].level == ltot && Ndirect + Nleafmax <= Nneibmax) {
        i = celldata[cc].ifirst;
        while (i != -1) {
          directlist[Ndirect++] = Nneib;
          neiblist[Nneib] = i;
          neibpart[Nneib++] = partdata[i];
          if (i == celldata[cc].ilast) break;
          i = inext[i];
        };
        cc = celldata[cc].cnext;
      }

      // If leaf-cell, but we've run out of memory, return with error-code (-1)
      else if (celldata[cc].level == ltot && Ndirect + Nleafmax > Nneibmax)
        return -3;

    }

    // If not in range, then open next cell
    //---------------------------------------------------------------------------------------------
    else
      cc = celldata[cc].cnext;

  };
  //===============================================================================================


  // Now, trim the list to remove particles that are definitely not SPH neighbours.
  // If not an SPH neighbour, then add to direct gravity sum list.
  for (j=Nsphneibtemp; j<Nsphneib; j++) {
    i = sphneiblist[j];
    if (neibpart[i].itype == dead) continue;
    for (k=0; k<ndim; k++) dr[k] = neibpart[i].r[k] - rc[k];
    drsqd = DotProduct(dr,dr,ndim);
    if (drsqd < hrangemaxsqd || drsqd <
        (rmax + kernrange*neibpart[i].h)*(rmax + kernrange*neibpart[i].h))
      sphneiblist[Nsphneibtemp++] = i;
    else if (Ndirect + Nsphneibtemp < Nneibmax)
      directlist[Ndirect++] = i;
    else
     return -3;
  }
  Nsphneib = Nsphneibtemp;


  return 1;
}



//=================================================================================================
//  KDTree::ComputePeriodicGravityInteractionList
/// Computes and returns number of SPH neighbours (Nneib), direct sum particles
/// (Ndirect) and number of cells (Ngravcell), including lists of ids, from
/// the gravity tree walk for active particles inside cell c.
/// Currently defaults to the geometric opening criteria.
/// If any of the interactions list arrays (neiblist,directlist,gravcelllist)
/// overflow, return with error code (-1) to reallocate more memory.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
int KDTree<ndim,ParticleType,TreeCell>::ComputePeriodicGravityInteractionList
 (const TreeCell<ndim> &cell,          ///< [in] Pointer to cell
  const ParticleType<ndim> *partdata,  ///< [in] Particle data array
  const DomainBox<ndim> &simbox,       ///< [in] Simulation domain box object
  const FLOAT macfactor,               ///< [in] Gravity MAC particle factor
  const int Nneibmax,                  ///< [in] Max. no. of SPH neighbours
  const int Ngravcellmax,              ///< [in] Max. no. of cell interactions
  int &Nneib,                          ///< [out] Total no. of neighbours
  int &Nsphneib,                       ///< [out] No. of SPH neighbours
  int &Ndirect,                        ///< [out] No. of direct-sum neighbours
  int &Ngravcell,                      ///< [out] No. of cell interactions
  int *neiblist,                       ///< [out] List of all particle ids
  int *sphneiblist,                    ///< [out] List of SPH neibpart ids
  int *directlist,                     ///< [out] List of direct-sum neibpart ids
  TreeCell<ndim> *gravcell,            ///< [out] Array of local copies of tree cells
  ParticleType<ndim> *neibpart)        ///< [out] Array of local copies of neighbour particles
{
  int cc;                              // Cell counter
  int i;                               // Particle id
  int j;                               // Aux. particle counter
  int k;                               // Neighbour counter
  int Nsphneibtemp = 0;                // Aux. counter
  FLOAT dr[ndim];                      // Relative position vector
  FLOAT drsqd;                         // Distance squared
  FLOAT rc[ndim];                      // Position of cell

  // Make local copies of important cell properties
  const FLOAT hrangemaxsqd = pow(cell.rmax + kernrange*cell.hmax,2);
  const FLOAT rmax = cell.rmax;
  for (k=0; k<ndim; k++) rc[k] = cell.rcell[k];

  // Start with root cell and walk through entire tree
  cc        = 0;
  Nneib     = 0;
  Nsphneib  = 0;
  Ndirect   = 0;
  Ngravcell = 0;


  // Walk through all cells in tree to determine particle and cell interaction lists
  //===============================================================================================
  while (cc < Ncell) {

    for (k=0; k<ndim; k++) dr[k] = celldata[cc].rcell[k] - rc[k];
    drsqd = DotProduct(dr,dr,ndim);


    // Check if bounding boxes overlap with each other
    //---------------------------------------------------------------------------------------------
    if (BoxOverlap(cell.bbmin,cell.bbmax,celldata[cc].hboxmin,celldata[cc].hboxmax) ||
        BoxOverlap(cell.hboxmin,cell.hboxmax,celldata[cc].bbmin,celldata[cc].bbmax)) {

      // If not a leaf-cell, then open cell to first child cell
      if (celldata[cc].level != ltot)
        cc++;

      else if (celldata[cc].N == 0)
        cc = celldata[cc].cnext;

      // If leaf-cell, add particles to list
      else if (celldata[cc].level == ltot && Nneib + Nleafmax <= Nneibmax) {
        i = celldata[cc].ifirst;
        while (i != -1) {
          sphneiblist[Nsphneib++] = Nneib;
          neiblist[Nneib] = i;
          neibpart[Nneib++] = partdata[i];
          if (i == celldata[cc].ilast) break;
          i = inext[i];
        };
        cc = celldata[cc].cnext;
      }

      // If leaf-cell, but we've run out of memory, return with error-code (-1)
      else if (celldata[cc].level == ltot && Nneib + Nleafmax > Nneibmax)
        return -1;

    }

    // Check if cell is far enough away to use the COM approximation
    //---------------------------------------------------------------------------------------------
    else if (drsqd > celldata[cc].cdistsqd && drsqd > celldata[cc].mac*macfactor
              && celldata[cc].N > 0) {

      // If cell is a leaf-cell with only one particle, more efficient to
      // compute the gravitational contribution from the particle than the cell
      if (celldata[cc].level == ltot && celldata[cc].N == 1 && Ndirect + Nneib < Nneibmax) {
        i = celldata[cc].ifirst;
        directlist[Ndirect++] = Nneib;
        neiblist[Nneib] = i;
        neibpart[Nneib++] = partdata[i];
      }
      else if (Ngravcell < Ngravcellmax) {
        gravcell[Ngravcell++] = celldata[cc];
      }
      else {
        return -2;
      }
      cc = celldata[cc].cnext;

    }

    // If cell is too close, open cell to interogate children cells.
    // If cell is too close and a leaf cell, then add particles to direct list.
    //---------------------------------------------------------------------------------------------
    else if (!(drsqd > celldata[cc].cdistsqd && drsqd > celldata[cc].mac*macfactor )
       && celldata[cc].N > 0) {

      // If not a leaf-cell, then open cell to first child cell
      if (celldata[cc].level != ltot)
         cc++;

      // If leaf-cell, add particles to list
      else if (celldata[cc].level == ltot && Ndirect + Nleafmax <= Nneibmax) {
        i = celldata[cc].ifirst;
        while (i != -1) {
          directlist[Ndirect++] = Nneib;
          neiblist[Nneib] = i;
          neibpart[Nneib++] = partdata[i];
          if (i == celldata[cc].ilast) break;
          i = inext[i];
        };
        cc = celldata[cc].cnext;
      }

      // If leaf-cell, but we've run out of memory, return with error-code (-1)
      else if (celldata[cc].level == ltot && Ndirect + Nleafmax > Nneibmax)
        return -3;

    }

    // If not in range, then open next cell
    //---------------------------------------------------------------------------------------------
    else
      cc = celldata[cc].cnext;

  };
  //===============================================================================================


  // Now, trim the list to remove particles that are definitely not SPH neighbours.
  // If not an SPH neighbour, then add to direct gravity sum list.
  for (j=Nsphneibtemp; j<Nsphneib; j++) {
    i = sphneiblist[j];
    if (neibpart[i].itype == dead) continue;
    for (k=0; k<ndim; k++) dr[k] = neibpart[i].r[k] - rc[k];
    drsqd = DotProduct(dr,dr,ndim);
    if (drsqd < hrangemaxsqd || drsqd <
        (rmax + kernrange*neibpart[i].h)*(rmax + kernrange*neibpart[i].h))
      sphneiblist[Nsphneibtemp++] = i;
    else if (Ndirect + Nsphneibtemp < Nneibmax)
      directlist[Ndirect++] = i;
    else
     return -3;
  }
  Nsphneib = Nsphneibtemp;


  return 1;
}




//=================================================================================================
//  KDTree::ComputeStarGravityInteractionList
/// Computes and returns number of SPH neighbours (Nneib), direct sum particles
/// (Ndirect) and number of cells (Ngravcell), including lists of ids, from
/// the gravity tree walk for active particles inside cell c.
/// Currently defaults to the geometric opening criteria.
/// If any of the interactions list arrays (neiblist,directlist,gravcelllist)
/// overflow, return with error code (-1) to reallocate more memory.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
int KDTree<ndim,ParticleType,TreeCell>::ComputeStarGravityInteractionList
 (const NbodyParticle<ndim> *star,     ///< [in] Nbody particle
  const FLOAT macfactor,               ///< [in] Gravity MAC factor
  const int Nneibmax,                  ///< [in] Max. no. of SPH neighbours
  const int Ndirectmax,                ///< [in] Max. no. of direct-sum neighbours
  const int Ngravcellmax,              ///< [in] Max. no. of cell interactions
  int &Nneib,                          ///< [out] No. of SPH neighbours
  int &Ndirect,                        ///< [out] No. of direct-sum neighbours
  int &Ngravcell,                      ///< [out] No. of cell interactions
  int *neiblist,                       ///< [out] List of SPH neighbour ids
  int *directlist,                     ///< [out] List of direct-sum neighbour ids
  TreeCell<ndim> *gravcell,            ///< [out] List of cell ids
  ParticleType<ndim> *partdata)        ///< [in] Particle data array
{
  int cc;                              // Cell counter
  int i;                               // Particle id
  int k;                               // Neighbour counter
  FLOAT dr[ndim];                      // Relative position vector
  FLOAT drsqd;                         // Distance squared
  FLOAT hrangemax;                     // Maximum kernel extent
  FLOAT rs[ndim];                      // Position of star

  // Make local copies of important cell properties
  for (k=0; k<ndim; k++) rs[k] = star->r[k];
  hrangemax = kernrange*star->h;

  // Start with root cell and walk through entire tree
  cc = 0;
  Nneib = 0;
  Ndirect = 0;
  Ngravcell = 0;


  // Walk through all cells in tree to determine particle and cell
  // interaction lists
  //===============================================================================================
  while (cc < Ncell) {

    for (k=0; k<ndim; k++) dr[k] = celldata[cc].rcell[k] - rs[k];
    drsqd = DotProduct(dr,dr,ndim);


    // Check if cells contain SPH neighbours
    //---------------------------------------------------------------------------------------------
    if (drsqd < pow(0.5*hrangemax + celldata[cc].rmax + 0.5*kernrange*celldata[cc].hmax,2)) {

      // If not a leaf-cell, then open cell to first child cell
      if (celldata[cc].level != ltot)
        cc++;

      // If leaf-cell, add particles to list
      else if (celldata[cc].level == ltot && Nneib + Nleafmax <= Nneibmax) {
        i = celldata[cc].ifirst;
        while (i != -1) {
          if (partdata[i].itype != dead) neiblist[Nneib++] = i;
          if (i == celldata[cc].ilast) break;
          i = inext[i];
        };
        cc = celldata[cc].cnext;
      }

      // If leaf-cell, but we've run out of memory, return with error-code (-1)
      else if (celldata[cc].level == ltot && Nneib + Nleafmax > Nneibmax)
        return -1;

    }

    // Check if cell is far enough away to use the COM approximation
    //---------------------------------------------------------------------------------------------
    else if (drsqd > celldata[cc].cdistsqd && celldata[cc].N > 0) {

      // If cell is a leaf-cell with only one particle, more efficient to
      // compute the gravitational contribution from the particle than the cell
      if (celldata[cc].level == ltot && celldata[cc].N == 1 &&
          Ndirect + Nneib < Ndirectmax) {
        if (partdata[celldata[cc].ifirst].itype != dead)
            directlist[Ndirect++] = celldata[cc].ifirst;
      }
      else if (Ngravcell < Ngravcellmax && celldata[cc].N > 0)
        gravcell[Ngravcell++] = celldata[cc];
      else
        return -1;
      cc = celldata[cc].cnext;

    }

    // If cell is too close, open cell to interogate children cells.
    // If cell is too close and a leaf cell, then add particles to direct list.
    //---------------------------------------------------------------------------------------------
    else if (drsqd <= celldata[cc].cdistsqd && celldata[cc].N > 0) {

      // If not a leaf-cell, then open cell to first child cell
      if (celldata[cc].level != ltot)
         cc++;

      // If leaf-cell, add particles to list
      else if (celldata[cc].level == ltot && Ndirect + Nleafmax <= Ndirectmax) {
        i = celldata[cc].ifirst;
        while (i != -1) {
          if (partdata[i].itype != dead) directlist[Ndirect++] = i;
          if (i == celldata[cc].ilast) break;
          i = inext[i];
        };
        cc = celldata[cc].cnext;
      }

      // If leaf-cell, but we've run out of memory, return with error-code (-1)
      else if (celldata[cc].level == ltot && Ndirect + Nleafmax > Ndirectmax)
        return -1;

    }

    // If not in range, then open next cell
    //---------------------------------------------------------------------------------------------
    else
      cc = celldata[cc].cnext;

  };
  //===============================================================================================


  return 1;
}



//=================================================================================================
//  KDTree::ComputeCellMonopoleForces
/// Compute the force on particle 'parti' due to all cells obtained in the
/// gravity tree walk.  Uses only monopole moments (i.e. COM) of the cell.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void KDTree<ndim,ParticleType,TreeCell>::ComputeCellMonopoleForces
 (FLOAT &gpot,                         ///< [inout] Grav. potential
  FLOAT agrav[ndim],                   ///< [inout] Acceleration array
  FLOAT rp[ndim],                      ///< [in] Position of point
  int Ngravcell,                       ///< [in] No. of tree cells in list
  TreeCell<ndim> *gravcell)            ///< [in] List of tree cell ids
{
  int cc;                              // Aux. cell counter
  int k;                               // Dimension counter
  FLOAT dr[ndim];                      // Relative position vector
  FLOAT drsqd;                         // Distance squared
  FLOAT invdrmag;                      // 1 / distance
  FLOAT invdrsqd;                      // 1 / drsqd
  FLOAT invdr3;                        // 1 / dist^3
  FLOAT mc;                            // Mass of cell
  TreeCell<ndim> *cellptr;             // Pointer to gravity tree cell

  // Loop over all neighbouring particles in list
  //-----------------------------------------------------------------------------------------------
  for (cc=0; cc<Ngravcell; cc++) {
    cellptr = &(gravcell[cc]);

    mc = cellptr->m;
    for (k=0; k<ndim; k++) dr[k] = cellptr->r[k] - rp[k];
    drsqd    = DotProduct(dr,dr,ndim) + small_number;
    invdrsqd = 1.0/drsqd;
    invdrmag = sqrt(invdrsqd);
    invdr3   = invdrsqd*invdrmag;

    gpot += mc*invdrmag;
    for (k=0; k<ndim; k++) agrav[k] += mc*dr[k]*invdr3;

  }
  //-----------------------------------------------------------------------------------------------

  return;
}



//=================================================================================================
//  KDTree::ComputeCellQuadrupoleForces
/// Compute the force on particle 'parti' due to all cells obtained in the
/// gravity tree walk including the quadrupole moment correction term.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void KDTree<ndim,ParticleType,TreeCell>::ComputeCellQuadrupoleForces
 (FLOAT &gpot,                         ///< [inout] Grav. potential
  FLOAT agrav[ndim],                   ///< [inout] Acceleration array
  FLOAT rp[ndim],                      ///< [in] Position of point
  int Ngravcell,                       ///< [in] No. of tree cells in list
  TreeCell<ndim> *gravcell)            ///< [in] List of tree cell ids
{
  int cc;                              // Aux. cell counter
  int k;                               // Dimension counter
  FLOAT dr[ndim];                      // Relative position vector
  FLOAT drsqd;                         // Distance squared
  FLOAT invdrsqd;                      // 1 / drsqd
  FLOAT invdrmag;                      // 1 / distance
  FLOAT invdr5;                        // 1 / distance^5
  FLOAT qfactor;                       // Constant factor for optimisation
  FLOAT qscalar;                       // Quadrupole moment scalar quantity
  TreeCell<ndim> *cellptr;             // Pointer to gravity tree cell


  // Loop over all neighbouring particles in list
  //-----------------------------------------------------------------------------------------------
  for (cc=0; cc<Ngravcell; cc++) {
    cellptr = &(gravcell[cc]);

    for (k=0; k<ndim; k++) dr[k] = cellptr->r[k] - rp[k];
    drsqd = DotProduct(dr,dr,ndim) + small_number;
    invdrsqd = 1.0/drsqd;
    invdrmag = sqrt(invdrsqd);
    invdr5 = invdrsqd*invdrsqd*invdrmag;

    // First add monopole term for acceleration
    for (k=0; k<ndim; k++) agrav[k] += cellptr->m*dr[k]*invdrsqd*invdrmag;

    // Now add quadrupole moment terms depending on dimensionality
    if (ndim == 3) {
      qscalar = cellptr->q[0]*dr[0]*dr[0] + cellptr->q[2]*dr[1]*dr[1] -
        (cellptr->q[0] + cellptr->q[2])*dr[2]*dr[2] +
         2.0*(cellptr->q[1]*dr[0]*dr[1] + cellptr->q[3]*dr[0]*dr[2] +
         cellptr->q[4]*dr[1]*dr[2]);
      qfactor = 2.5*qscalar*invdr5*invdrsqd;
      agrav[0] +=
        (cellptr->q[0]*dr[0] + cellptr->q[1]*dr[1] + cellptr->q[3]*dr[2])*invdr5 - qfactor*dr[0];
      agrav[1] +=
        (cellptr->q[1]*dr[0] + cellptr->q[2]*dr[1] + cellptr->q[4]*dr[2])*invdr5 - qfactor*dr[1];
      agrav[2] +=
        (cellptr->q[3]*dr[0] + cellptr->q[4]*dr[1] -
         (cellptr->q[0] + cellptr->q[2])*dr[2])*invdr5 - qfactor*dr[2];
      gpot += cellptr->m*invdrmag + 0.5*qscalar*invdr5;
    }
    else if (ndim == 2) {
      qscalar = cellptr->q[0]*dr[0]*dr[0] + cellptr->q[2]*dr[1]*dr[1] +
        2.0*cellptr->q[1]*dr[0]*dr[1];
      qfactor = 2.5*qscalar*invdr5*invdrsqd;
      agrav[0] += (cellptr->q[0]*dr[0] + cellptr->q[1]*dr[1])*invdr5 - qfactor*dr[0];
      agrav[1] += (cellptr->q[1]*dr[0] + cellptr->q[2]*dr[1])*invdr5 - qfactor*dr[1];
      gpot += cellptr->m*invdrmag + 0.5*qscalar*invdr5;
    }

  }
  //-----------------------------------------------------------------------------------------------


  return;
}



//=================================================================================================
//  KDTree::ComputeFastMonopoleForces
/// Compute the force on particle 'parti' due to all cells obtained in the
/// gravity tree walk.  Uses only monopole moments (i.e. COM) of the cell.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void KDTree<ndim,ParticleType,TreeCell>::ComputeFastMonopoleForces
 (int Nactive,                         ///< [in] No. of active particles
  int Ngravcell,                       ///< [in] No. of tree cells in list
  TreeCell<ndim> *gravcell,            ///< [in] List of tree cell ids
  TreeCell<ndim> &cell,                ///< [in] Current cell pointer
  ParticleType<ndim> *activepart)      ///< [inout] Active SPH particle array
{
  int cc;                              // Aux. cell counter
  int j;                               // ..
  int k;                               // Dimension counter
  FLOAT ac[ndim];                      // ..
  FLOAT dr[ndim];                      // Relative position vector
  FLOAT drsqd;                         // Distance squared
  FLOAT invdrmag;                      // 1 / distance
  FLOAT invdrsqd;                      // 1 / drsqd
  FLOAT invdr3;                        // 1 / dist^3
  FLOAT mc;                            // Mass of cell
  FLOAT q[6];                          // Local copy of quadrupole moment
  FLOAT dphi[3];                       // ..
  FLOAT cellpot;                       // ..
  FLOAT rc[ndim];                      // ..

  for (k=0; k<ndim; k++) rc[k] = cell.r[k];
  for (k=0; k<ndim; k++) ac[k] = 0.0;
  for (k=0; k<ndim; k++) dphi[k] = 0.0;
  for (k=0; k<6; k++) q[k] = 0;
  cellpot = 0.0;


  //-----------------------------------------------------------------------------------------------
  if (ndim == 3) {

    for (cc=0; cc<Ngravcell; cc++) {
#ifndef MPI_PARALLEL
      assert(cell.id != gravcell[cc].id);
#endif
      mc = gravcell[cc].m;
      for (k=0; k<ndim; k++) dr[k] = gravcell[cc].r[k] - rc[k];
      drsqd    = DotProduct(dr,dr,ndim);
      invdrsqd = 1.0/drsqd;
      invdrmag = sqrt(invdrsqd);
      invdr3   = invdrsqd*invdrmag;
      cellpot  += mc*invdrmag;
      for (k=0; k<ndim; k++) ac[k] += mc*dr[k]*invdr3;
      for (k=0; k<ndim; k++) dphi[k] += mc*dr[k]*invdr3;
      q[0] += mc*(3.0*dr[0]*dr[0]*invdr3*invdrsqd - invdrsqd*invdrmag);
      q[1] += mc*(3.0*dr[0]*dr[1]*invdr3*invdrsqd);
      q[2] += mc*(3.0*dr[1]*dr[1]*invdr3*invdrsqd - invdrsqd*invdrmag);
      q[3] += mc*(3.0*dr[2]*dr[0]*invdr3*invdrsqd);
      q[4] += mc*(3.0*dr[2]*dr[1]*invdr3*invdrsqd);
      q[5] += mc*(3.0*dr[2]*dr[2]*invdr3*invdrsqd - invdrsqd*invdrmag);
    }

    for (j=0; j<Nactive; j++) {
      for (k=0; k<ndim; k++) dr[k] = activepart[j].r[k] - rc[k];
      activepart[j].agrav[0] += ac[0] + q[0]*dr[0] + q[1]*dr[1] + q[3]*dr[2];
      activepart[j].agrav[1] += ac[1] + q[1]*dr[0] + q[2]*dr[1] + q[4]*dr[2];
      activepart[j].agrav[2] += ac[2] + q[3]*dr[0] + q[4]*dr[1] + q[5]*dr[2];
      activepart[j].gpot += cellpot + dphi[0]*dr[0] + dphi[1]*dr[1] + dphi[2]*dr[2];
    }

  }
  //-----------------------------------------------------------------------------------------------

  return;
}



//=================================================================================================
//  KDTree::ComputeActiveCellList
/// Returns the number of cells containing active particles, 'Nactive', and
/// the i.d. list of cells contains active particles, 'celllist'
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
int KDTree<ndim,ParticleType,TreeCell>::ComputeActiveCellList
 (TreeCell<ndim> *celllist)            ///< Array containing copies of cells with active ptcls
{
  int c;                               // Cell counter
  int Nactive = 0;                     // No. of active leaf cells in tree

  for (c=0; c<Ncell; c++) {
    if (celldata[c].level == lactive && celldata[c].Nactive > 0) {
      celllist[Nactive++] = celldata[c];
    }
  }

#ifdef MPI_PARALLEL
    for (c=Ncell; c<Ncell+Nimportedcell; c++) {
      if (celldata[c].Nactive > 0) celllist[Nactive++] = celldata[c];
    }
#endif

  return Nactive;
}



#ifdef MPI_PARALLEL
//=================================================================================================
//  KDTree::ComputeDistantGravityInteractionList
/// Computes and returns number of SPH neighbours (Nneib), direct sum particles
/// (Ndirect) and number of cells (Ngravcell), including lists of ids, from
/// the gravity tree walk for active particles inside cell c.
/// Currently defaults to the geometric opening criteria.
/// If any of the interactions list arrays (neiblist,directlist,gravcelllist)
/// overflow, return with error code (-1) to reallocate more memory.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
int KDTree<ndim,ParticleType,TreeCell>::ComputeDistantGravityInteractionList
(const TreeCell<ndim> *cellptr,       ///< [in] Pointer to cell
 const FLOAT macfactor,               ///< [in] Gravity MAC particle factor
 const int Ngravcellmax,              ///< [in] Max. no. of cell interactions
 int Ngravcell,                       ///< [in] Current no. of cells in array
 TreeCell<ndim> **gravcelllist)      ///< [out] Array of pointers to cells
{
  int cc;                             // Cell counter
  int i;                              // Particle id
  int j;                              // Aux. particle counter
  int k;                              // Neighbour counter
  int Ngravcelltemp = Ngravcell;      // ..
  FLOAT dr[ndim];                     // Relative position vector
  FLOAT drsqd;                        // Distance squared
  FLOAT rc[ndim];                     // Position of cell
  FLOAT hrangemax;                    // Maximum kernel extent
  FLOAT rmax;                         // Radius of sphere containing particles

  // Make local copies of important cell properties
  for (k=0; k<ndim; k++) rc[k] = cellptr->rcell[k];
  hrangemax = cellptr->rmax + kernrange*cellptr->hmax;
  rmax = cellptr->rmax;

  // Start with root cell and walk through entire tree
  cc = 0;


  // Walk through all cells in tree to determine particle and cell
  // interaction lists
  //===============================================================================================
  while (cc < Ncell) {

    for (k=0; k<ndim; k++) dr[k] = celldata[cc].rcell[k] - rc[k];
    drsqd = DotProduct(dr,dr,ndim);


    // Check if bounding boxes overlap with each other
    //---------------------------------------------------------------------------------------------
    if (BoxOverlap(cellptr->bbmin,cellptr->bbmax,
                   celldata[cc].hboxmin,celldata[cc].hboxmax) ||
        BoxOverlap(cellptr->hboxmin,cellptr->hboxmax,
                   celldata[cc].bbmin,celldata[cc].bbmax)) {

      // If not a leaf-cell, then open cell to first child cell
      if (celldata[cc].level != ltot)
        cc++;

      else if (celldata[cc].N == 0)
        cc = celldata[cc].cnext;

      // If leaf-cell, add particles to list
      else if (celldata[cc].level == ltot) {
        return -1;
      }

    }

    // Check if cell is far enough away to use the COM approximation
    //---------------------------------------------------------------------------------------------
    else if (drsqd > celldata[cc].cdistsqd && drsqd > celldata[cc].mac*macfactor
	     && celldata[cc].N > 0) {

      gravcelllist[Ngravcelltemp++] = &(celldata[cc]);
      if (Ngravcelltemp>=Ngravcellmax) {
        ExceptionHandler::getIstance().raise("Too many interaction cells in disant gravity interaction list!!!");
      }
      cc = celldata[cc].cnext;

    }

    // If cell is too close, open cell to interogate children cells.
    // If cell is too close and a leaf cell, then add particles to direct list.
    //---------------------------------------------------------------------------------------------
    else if (!(drsqd > celldata[cc].cdistsqd &&
	       drsqd > celldata[cc].mac*macfactor) && celldata[cc].N > 0) {

      // If not a leaf-cell, then open cell to first child cell
      if (celldata[cc].level != ltot)
         cc++;

      // If leaf-cell, add particles to list
      else
       	return -1;

    }

    // If not in range, then open next cell
    //---------------------------------------------------------------------------------------------
    else
      cc = celldata[cc].cnext;

  };
  //===============================================================================================

  return Ngravcelltemp;
}



//=================================================================================================
//  KDTree::ComputeHydroTreeCellOverlap
/// ...
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
bool KDTree<ndim,ParticleType,TreeCell>::ComputeHydroTreeCellOverlap
 (const TreeCell<ndim> *cellptr)       ///< [in] Pointer to cell
{
  int cc = 0;                          // Cell counter

  // Walk through all cells in tree to determine particle and cell interaction lists
  //===============================================================================================
  while (cc < Ncell) {

    // Check if bounding boxes overlap with each other
    //---------------------------------------------------------------------------------------------
    if (BoxOverlap(cellptr->bbmin,cellptr->bbmax,celldata[cc].hboxmin,celldata[cc].hboxmax) ||
        BoxOverlap(cellptr->hboxmin,cellptr->hboxmax,celldata[cc].bbmin,celldata[cc].bbmax)) {

      // If not a leaf-cell, then open cell to first child cell
      if (celldata[cc].level != ltot) cc++;

      // If cell contains no particle (dead leaf?) then move to next cell
      //else if (celldata[cc].N == 0) cc = celldata[cc].cnext;

      // If cell is overlapping with any leaf, then flag overlap on return
      else if (celldata[cc].level == ltot) return true;

    }

    // If not in range, then open next cell
    //---------------------------------------------------------------------------------------------
    else
      cc = celldata[cc].cnext;

  };
  //===============================================================================================

  // If we've walked the entire tree wihout any leaf overlaps, flag no overlap
  return false;
}
#endif



#if defined(VERIFY_ALL)
//=================================================================================================
//  KDTree::ValidateTree
/// Performs various sanity and validation checks on KD-tree structure.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void KDTree<ndim,ParticleType,TreeCell>::ValidateTree
(ParticleType<ndim> *partdata)      ///< Pointer to SPH class
{
  bool overlap_flag = false;        // Flag if cell bounding boxes overlap
  bool hmax_flag = false;           // Flag if ptcls have larger h than hmax
  bool kill_flag = false;           // ..
  int activecount;                  // Active particles in leaf cell
  int c;                            // Cell counter
  int cc;                           // Aux. cell counter
  int i;                            // Particle counter
  int j;                            // Aux. particle counter
  int l;                            // Tree level
  int leafcount;                    // Leaf cell counter
  int Nactivecount=0;               // Counter for total no. of active ptcls
  int Ncount=0;                     // Total particle counter
  int *ccount;                      // Array for counting cells
  int *lcount;                      // Array for counting ptcls on each level
  int *pcount;                      // Array for counting particles in tree
  TreeCell<ndim> cell;              // Local copy of KD-tree cell

  debug2("[KDTree::ValidateTree]");

  ccount = new int[Ncellmax];
  pcount = new int[Ntotmax];
  lcount = new int[ltot+1];
  for (i=0; i<Ntotmax; i++) pcount[i] = 0;
  for (c=0; c<Ncellmax; c++) ccount[c] = 0;
  for (l=0; l<ltot; l++) lcount[l] = 0;
  Ncount = 0;
  Nactivecount = 0;

  // Count how many times we enter a cell in a full tree walk
  c = 0;
  while (c < Ncell) {
    ccount[c]++;
    if (celldata[c].c1 != -1) c = celldata[c].c1;
    else c = celldata[c].cnext;
  }

  // Now check we enter all cells once and once only
  for (c=0; c<Ncell; c++) {
    if (ccount[c] != 1) {
      cout << "Error in cell walk count : " << ccount[c] << endl;
      PrintArray("ccount     : ",Ncell,ccount);
      exit(0);
    }
  }

  // Check inext linked list values and ids array are all valid
  for (i=ifirst; i<=ilast; i++) {
    if (!(ids[i] >= ifirst && ids[i] <= ilast)) {
      cout << "Problem with ids array : " << i << "   " << ids[i] << endl;
      exit(0);
    }
    if (!(inext[i] >= -1)) {
      cout << "Problem with inext linked lists : " << i << "   " << inext[i] << endl;
      exit(0);
    }
  }


  // Verify linked lists are valid for all levels of tree
  //-----------------------------------------------------------------------------------------------
  for (l=0; l<=ltot; l++) {
    for (i=0; i<Ntotmax; i++) pcount[i] = 0;

    for (c=0; c<Ncell; c++) {
      cell = celldata[c];

      // Check that particles are not in linked lists more than once
      if (cell.level == l) {
        i = cell.ifirst;
        while (i != -1) {
          pcount[i]++;
          if (i == cell.ilast) break;
          i = inext[i];
        }

      }
    }

    // Check particles are included in the tree once and once only
    for (i=ifirst; i<=ilast; i++) {
      if (pcount[i] != 1) {
        cout << "Problem with linked lists on level : " << l
             << " for particle : " << i << "   " << pcount[i] << "    " << ltot << endl;
        kill_flag = true;
      }
    }

  }
  for (i=0; i<Ntotmax; i++) pcount[i] = 0;


  // Loop over all cells in tree
  //-----------------------------------------------------------------------------------------------
  for (c=0; c<Ncell; c++) {
    cell = celldata[c];
    activecount = 0;
    leafcount = 0;

    // Add particles from current level
    lcount[cell.level] += cell.N;

    // Check that particles are not in linked lists more than once
    if (cell.level == ltot) {
      i = cell.ifirst;
      while (i != -1) {
        pcount[i]++;
        leafcount++;
        Ncount++;
        if (partdata[i].active) activecount++;
        if (partdata[i].active) Nactivecount++;
        if (partdata[i].h > cell.hmax) {
          cout << "hmax flag error : " << c << "    "
               << partdata[i].h << "   " << cell.hmax << endl;
          exit(0);
        }
        if (i == cell.ilast) break;
        i = inext[i];
      }
      if (leafcount > Nleafmax) {
        cout << "Leaf particle count error : " << leafcount << "   " << Nleafmax << endl;
        exit(0);
      }
      if (activecount > leafcount) {
        cout << "Leaf particle count error : " << leafcount << "   " << Nleafmax << endl;
        exit(0);
      }
    }


    // Check that child cell boxes fit inside parent box
    if (cell.level != ltot) {
      KDTreeCell<ndim>& c1 = celldata[cell.c1];
      KDTreeCell<ndim>& c2 = celldata[cell.c2];

      for (int k=0; k<ndim; k++) {
        if (c1.hboxmin[k] < cell.hboxmin[k] || c1.hboxmax[k] > cell.hboxmax[k] ||
            c2.hboxmin[k] < cell.hboxmin[k] || c2.hboxmax[k] > cell.hboxmax[k])
            overlap_flag = true;
      }
      if (overlap_flag) {
        cout << "Parent/daughter cell overlap error!! : "
             << c << "   " << cell.c1 << "    " << cell.c2 << endl;
        exit(0);
      }
    }

    // Check that bounding boxes of cells on each level do not overlap each other
    for (cc=0; cc<Ncell; cc++) {
      if (c != cc && celldata[cc].level == cell.level) {
        if (ndim == 2) {
          if (cell.bbmin[0] < celldata[cc].bbmax[0] &&
              cell.bbmax[0] > celldata[cc].bbmin[0] &&
              cell.bbmin[1] < celldata[cc].bbmax[1] &&
              cell.bbmax[1] > celldata[cc].bbmin[1]) {
            overlap_flag = true;
          }
          else if (ndim == 3) {
            if (cell.bbmin[0] < celldata[cc].bbmax[0] &&
                cell.bbmax[0] > celldata[cc].bbmin[0] &&
                cell.bbmin[1] < celldata[cc].bbmax[1] &&
                cell.bbmax[1] > celldata[cc].bbmin[1] &&
                cell.bbmin[2] < celldata[cc].bbmax[2] &&
                cell.bbmax[2] > celldata[cc].bbmin[2])
              overlap_flag = true;
          }
        }
      }
      if (overlap_flag) {
        cout << "Brother/sister cell overlap error!! : " << c << "   " << cc << endl;
        exit(0);
      }
    }
  }
  //-----------------------------------------------------------------------------------------------


  // Check particles are included in the tree once and once only
  for (i=ifirst; i<=ilast; i++) {
    if (pcount[i] != 1) {
      cout << "Problem with child cell ptcl counter : " << i << "   " << pcount[i] << endl;
      kill_flag = true;
    }
  }

  // Check all particles accounted for
  if (Ncount != Ntot) {
    cout << "Ncount problem with KD-tree : " << Ncount << "   " << Ntot << endl;
    kill_flag = true;
  }

  // Check active particles don't exceed total number of particles
  if (Nactivecount > Ntot) {
    cout << "Nactivecount problem with KD-tree : " << Nactivecount << "   " << Ntot << endl;
    kill_flag = true;
  }

  // Check number of particles on all levels is consistent
  for (l=0; l<ltot; l++) {
    if (lcount[l] != Ntot) {
      cout << "Problem with SPH particles on level : " << l << "    "
           << lcount[l] << "    " << Ntot << endl;
      kill_flag = true;
    }
  }

  delete[] lcount;
  delete[] pcount;
  delete[] ccount;

  if (kill_flag) exit(0);

  return;
}
#endif



template class KDTree<1,SphParticle,KDTreeCell>;
template class KDTree<2,SphParticle,KDTreeCell>;
template class KDTree<3,SphParticle,KDTreeCell>;
template class KDTree<1,GradhSphParticle,KDTreeCell>;
template class KDTree<2,GradhSphParticle,KDTreeCell>;
template class KDTree<3,GradhSphParticle,KDTreeCell>;
template class KDTree<1,SM2012SphParticle,KDTreeCell>;
template class KDTree<2,SM2012SphParticle,KDTreeCell>;
template class KDTree<3,SM2012SphParticle,KDTreeCell>;
template class KDTree<1,GodunovSphParticle,KDTreeCell>;
template class KDTree<2,GodunovSphParticle,KDTreeCell>;
template class KDTree<3,GodunovSphParticle,KDTreeCell>;
