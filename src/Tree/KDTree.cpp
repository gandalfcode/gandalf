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


#include <algorithm>
#include <cstdlib>
#include <cassert>
#include <iostream>
#include <string>
#include <math.h>
#include "Precision.h"
#include "Exception.h"
#include "DomainBox.h"
#include "Parameters.h"
#include "InlineFuncs.h"
#include "Particle.h"
#include "Sph.h"
#include "KDTree.h"
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
                                           string gravity_mac_aux, multipole_method multipole_aux,
                                           const DomainBox<ndim>& domain,
                                  		   const ParticleTypeRegister& reg,
										   const bool IAmPruned):
  Tree<ndim,ParticleType,TreeCell>(Nleafmaxaux, thetamaxsqdaux, kernrangeaux,
                                   macerroraux, gravity_mac_aux, multipole_aux, domain, reg,
                                   IAmPruned)
{
  allocated_tree = false;
  gmax           = 0;
  gtot           = 0;
  ifirst         = -1;
  ilast          = -1;
  lmax           = 0;
  ltot           = 0;
  ltot_old       = -1;
  Ncell          = 0;
  Ncellmax       = 0;
  Ntot           = 0;
  Ntotmax        = 0;
  hmax           = 0.0;
#if defined _OPENMP
  Nthreads       = omp_get_max_threads();
#else
  Nthreads       = 1;
#endif
#if defined MPI_PARALLEL
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
void KDTree<ndim,ParticleType,TreeCell>::AllocateTreeMemory(int Nparticles, int Ncells, bool force_realloc)
{
  debug2("[KDTree::AllocateTreeMemory]");

  //ComputeTreeSize();

  if (!allocated_tree || Nparticles > Ntotmax ||  Ncells > Ncellmax || force_realloc) {
    if (allocated_tree) DeallocateTreeMemory();

    Ncells = max(Ncells,Ncellmax);
    Nparticles = max(Nparticles,Ntotmax);

    g2c      = new int[gmax];
    celldata = new struct TreeCell<ndim>[Ncells];

    allocated_tree = true;

    Ntotmax = Nparticles;
    Ncellmax = Ncells;

    //CreateTreeStructure();

  }

  assert(Ntotmax>=Ntot);
  assert(Ncell<=Ncellmax);

  //CreateTreeStructure();

  return;
}

//=================================================================================================
//  KDTree::ReallocateMemory
/// Reallocate memory for KD-tree (when we need to grow the tree) as requested. Preserves the existing
/// information in the tree, differently from the previous function
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void KDTree<ndim,ParticleType,TreeCell>::ReallocateMemory(int Nparticles, int Ncells)
{
  debug2("[KDTree::ReallocateMemory]");

  if (!allocated_tree) {
	  ExceptionHandler::getIstance().raise("This function should not be called if the tree has not been allocated yet!");
  }


  if (Nparticles > Ntotmax ) {
	  Ntotmax = Nparticles;
  }

  if (Ncells > Ncellmax) {


    TreeCell<ndim>* celldataold = celldata;

    celldata = new struct TreeCell<ndim>[Ncells];

    std::copy(celldataold,celldataold+Ncellmax,celldata);

    delete[] celldataold;

	Ncellmax=Ncells;


  }


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
 (const int _ifirst,                   ///< i.d. of first particle
  const int _ilast,                    ///< i.d. of last particle
  const int Npart,                     ///< No. of particles
  const int Npartmax,                  ///< Max. no. of particles
  const FLOAT timestep,                ///< Smallest physical timestep
  Particle<ndim> *part_gen)            ///< Particle data array
{
  int i;                               // Particle counter
  int k;                               // Dimension counter
  FLOAT bbmin[ndim];                   // Minimum extent of local bounding box
  FLOAT bbmax[ndim];                   // Maximum extent of local bounding box

  debug2("[KDTree::BuildTree]");
  //timing->StartTimingSection("BUILD_TREE");

  ParticleType<ndim>* partdata = reinterpret_cast<ParticleType<ndim>*>(part_gen) ;

  // Activate nested parallelism for tree building routines
#ifdef _OPENMP
  omp_set_nested(1);
#endif

  // Set no. of tree members to total number of SPH particles (inc. ghosts)
  gtot = 0;
  ltot_old = ltot;

  bool force_realloc=false;
  if (Npartmax > Ntotmax) force_realloc=true;

  Ntotmax = Npartmax;

  // Compute the size of all tree-related arrays now we know number of points
  const int gmax_old = gmax;
  const int Ncellmax_old = Ncellmax;
  ComputeTreeSize();

  // Allocate (or reallocate if needed) all tree memory
  if (gmax_old<gmax || Ncellmax>Ncellmax_old) force_realloc=true;
  AllocateTreeMemory(Npartmax,0,force_realloc);


  // If the number of levels in the tree has changed (due to destruction or
  // creation of new particles) then re-create tree data structure
  // including linked lists and cell pointers
  CreateTreeStructure();
  //if (ltot != ltot_old) CreateTreeStructure();


  // Create bounding box of SPH particles
  for (k=0; k<ndim; k++) bbmin[k] = big_number;
  for (k=0; k<ndim; k++) bbmax[k] = -big_number;

  if (Npart > 0) {
    ifirst = _ifirst;
    ilast  = _ifirst + Npart - 1;
    for (i=ifirst; i<=ilast; i++) {
      for (k=0; k<ndim; k++) {
        bbmax[k] = max(bbmax[k], partdata[i].r[k] + kernrange*partdata[i].h);
        bbmin[k] = min(bbmin[k], partdata[i].r[k] - kernrange*partdata[i].h);
      }
    }
  }
  else {
    ifirst = -1;
    ilast  = -1;
  }


  // Set properties for root cell before constructing tree
  celldata[0].N      = Ntot;
  celldata[0].ifirst = ifirst;
  celldata[0].ilast  = ilast;
  celldata[0].cnext  = Ncellmax;
  celldata[0].hmax = 0;
  for (k=0; k<ndim; k++) celldata[0].bb.min[k] = bbmin[k];
  for (k=0; k<ndim; k++) celldata[0].bb.max[k] = bbmax[k];
  for (k=0; k<ndim; k++) celldata[0].v[k]= (FLOAT) 0.0;


  // If there are particles in the tree, recursively build tree from root node down
  if (Ntot > 0) {
    DivideTreeCell(ifirst, ilast, partdata, celldata[0]);
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
  cout << "No. of levels on tree : " << ltot << "   " << lmax << endl;
  cout << "No. of cells in tree  : " << Ncell << "   " << Ncellmax << endl;
  cout << "Max No part in leaf cell: " << Nleafmax << "\n";
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
    celldata[c].c2g     = 0;
    celldata[c].copen   = -1;
    celldata[c].cnext   = -1;
    celldata[c].c1      = -1;
    celldata[c].c2      = -1;
    celldata[c].ifirst  = -1;
    celldata[c].ilast   = -1;
    celldata[c].N       = 0;
    celldata[c].Nactive = 0;
  }
  g = 0;
  celldata[0].level = 0;


  // Loop over all cells and set all other pointers
  //-----------------------------------------------------------------------------------------------
  for (c=0; c<Ncell; c++) {
    celldata[c].id = c;

    if (celldata[c].level == ltot) {                           // If on leaf level
      celldata[c].cnext = c + 1;                               // id of next cell
      celldata[c].c2g = g;                                     // Record leaf id
      assert (g<gmax);
      g2c[g++] = c;                                            // Record inverse id
    }
    else {
      celldata[c+1].level            = celldata[c].level + 1;        // Level of 1st child
      celldata[c].copen              = c + 1;                        // id of 1st child
      celldata[c].c1                 = c + 1;                        // id of 1st child
      celldata[c].c2                 = c + c2L[celldata[c].level];   // id of 2nd child
      celldata[celldata[c].c2].level = celldata[c].level + 1;        // Level of 2nd child
      celldata[c].cnext              = c + cNL[celldata[c].level];   // Next cell id
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
#ifdef MPI_PARALLEL
    cell.worktot = 0.0;
#endif
    if (cell.N > 0) {
      assert(cell.ilast - cell.ifirst == cell.N - 1);
    }
    else {
      cell.ifirst = -1;
      cell.ilast = -1;
    }
    StockCellProperties(cell,partdata,true);
    return;
  }

#if defined(MPI_PARALLEL)
  i = cell.ifirst;
  for (k=0; k< ndim; k++) {
    cell.bb.min[k] = +big_number;
    cell.bb.max[k] = -big_number;
  }
  for (i=cell.ifirst; i<=cell.ilast; i++) {
    int j = ids[i];
    for (k=0; k<ndim; k++) {
      if (cell.bb.min[k] > partdata[j].r[k]) cell.bb.min[k] = partdata[j].r[k];
      if (cell.bb.max[k] < partdata[j].r[k]) cell.bb.max[k] = partdata[j].r[k];
    }
  }
#endif

  // Determine dimension to split the cell along.
  // For now, simply split along direction of the bounding box's longest axis
  for (k=0; k<ndim; k++) {
    if (cell.bb.max[k] - cell.bb.min[k] > rkmax) {
      rkmax = cell.bb.max[k] - cell.bb.min[k];
      k_divide = k;
      //cout << "Division? : " << k_divide << "    " << rkmax << endl;
    }
  }
  //cell.k_divide = k_divide;
  //cin >> i;


  // Find median value along selected division plane and re-order array
  // so particles reside on correct side of division
  rdivide = QuickSelect(cell.ifirst, cell.ilast, cell.ifirst+cell.N/2, k_divide, partdata);

  // Set properties of first child cell
  for (k=0; k<ndim; k++) celldata[cell.c1].bb.min[k] = cell.bb.min[k];
  for (k=0; k<ndim; k++) celldata[cell.c1].bb.max[k] = cell.bb.max[k];
  celldata[cell.c1].bb.max[k_divide] = rdivide;
  celldata[cell.c1].N = cell.N/2;
  if (celldata[cell.c1].N != 0) {
    celldata[cell.c1].ifirst = ifirst;
    celldata[cell.c1].ilast = ifirst + cell.N/2 - 1;
    assert(celldata[cell.c1].ilast - celldata[cell.c1].ifirst == celldata[cell.c1].N - 1);
  }

  // Set properties of second child cell
  for (k=0; k<ndim; k++) celldata[cell.c2].bb.min[k] = cell.bb.min[k];
  for (k=0; k<ndim; k++) celldata[cell.c2].bb.max[k] = cell.bb.max[k];
  celldata[cell.c2].bb.min[k_divide] = rdivide;
  celldata[cell.c2].N = cell.N - celldata[cell.c1].N;
  if (celldata[cell.c2].N != 0) {
    celldata[cell.c2].ifirst = ifirst + cell.N/2;
    celldata[cell.c2].ilast = ilast;
    assert(celldata[cell.c2].ilast - celldata[cell.c2].ifirst == celldata[cell.c2].N - 1);
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
  assert(cell.N == celldata[cell.c1].N + celldata[cell.c2].N);
  assert(!(cell.ifirst == -1 && cell.ilast == -1));

  // Stock all cell properties once constructed
  StockCellProperties(cell,partdata,true);

  return;
}



//=================================================================================================
//  KDTree::QuickSelectSort
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
  int j;                               // ..
  int jguess;                          // ..
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


// Predicate for comparing particles and sorting particles
template<class ParticleType>
class ParticleSorter {

public:
  ParticleSorter(int jdir)
    : _j(jdir) { } ;

  bool operator()(const ParticleType& l, const ParticleType& r) const {
    return l.r[_j] < r.r[_j] ;
  }
private:
  int _j ;
};



//=================================================================================================
//  KDTree::QuickSelect
/// Find median and sort particles in arrays to ensure they are the correct
/// side of the division.  Uses the QuickSelect algorithm.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
FLOAT KDTree<ndim,ParticleType,TreeCell>::QuickSelect
 (int left,                            ///< Left-most id of particle in array
  int right,                           ///< Right-most id of particle in array
  int jpivot,                          ///< Pivot/median point
  int k,                               ///< Dimension of sort
  ParticleType<ndim> *partdata)        ///< Pointer to main SPH object
{
  ParticleSorter<ParticleType<ndim> > pred(k) ;
  std::nth_element(partdata+left, partdata+jpivot, partdata+right+1, pred);
  return partdata[jpivot].r[k];
}



//=================================================================================================
//  KDTree::StockTree
/// Stock given tree cell in KD-tree.  If cell is not a leaf-cell, recursively
/// calls itself for its two child cells.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void KDTree<ndim,ParticleType,TreeCell>::StockTree
(TreeCell<ndim> &cell,                ///< Reference to cell to be stocked
 ParticleType<ndim> *partdata,        ///< SPH particle data array
 bool stock_leaf)					  ///< Whether to stock leaf cells
{
  int i;                               // Aux. child cell counter

  // If cell is not leaf, stock child cells
  if (cell.copen != -1) {
	  TreeCell<ndim>& child1 = celldata[cell.copen];
	  TreeCell<ndim>& child2 = celldata[child1.cnext];
#if defined _OPENMP
    if (pow(2,cell.level) < Nthreads) {
#pragma omp parallel for default(none) private(i) shared(cell,partdata, stock_leaf,child1,child2) num_threads(2)
      for (i=0; i<2; i++) {
        if (i == 0) StockTree(child1,partdata, stock_leaf);
        else if (i == 1) StockTree(child2,partdata, stock_leaf);
      }
    }
    else {
      for (i=0; i<2; i++) {
        if (i == 0) StockTree(child1,partdata, stock_leaf);
        else if (i == 1) StockTree(child2,partdata, stock_leaf);
      }
    }
#else
    for (i=0; i<2; i++) {
      if (i == 0) StockTree(child1,partdata, stock_leaf);
      else if (i == 1) StockTree(child2,partdata, stock_leaf);
    }
#endif
  }

  // Stock node once all children are stocked
  StockCellProperties(cell,partdata,stock_leaf);

  return;
}



//=================================================================================================
//  KDTree::StockCellProperties
/// Calculate the physical properties (e.g. total mass, centre-of-mass,
/// MAC opening-distance, etc..) of all given cell in the tree.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void KDTree<ndim,ParticleType,TreeCell>::StockCellProperties
 (TreeCell<ndim> &cell,                ///< Reference to current tree cell
  ParticleType<ndim> *partdata,        ///< Particle data array
  bool stock_leaf)					   ///< Whether to stock leaf cells
{
  int i;                               // Particle counter
  int iaux;                            // Aux. particle i.d. variable
  int k;                               // Dimension counter
  FLOAT dr[ndim];                      // Relative position vector
  FLOAT drsqd;                         // Distance squared
  FLOAT mi;                            // Mass of particle i
  FLOAT p = (FLOAT) 0.0;               // ..
  FLOAT lambda = (FLOAT) 0.0;          // ..


  const bool need_quadrupole_moments =
      multipole == quadrupole || multipole == fast_quadrupole || gravity_mac == eigenmac ;

  // Zero all summation variables for all cells
  if ((cell.level==ltot&&stock_leaf) || cell.copen != -1 ) {
	  cell.Nactive  = 0;
	  cell.N        = 0;
	  cell.maxsound = 0.0f;
	  cell.m        = (FLOAT) 0.0;
	  cell.hmax     = (FLOAT) 0.0;
	  cell.rmax     = (FLOAT) 0.0;
	  cell.mac      = (FLOAT) 0.0;
	  cell.cdistsqd = big_number;
      if (gravity_mac == gadget2)
        cell.amin = big_number ;
      else if (gravity_mac == eigenmac)
        cell.macfactor = 0 ;
      for (k=0; k<5; k++) cell.q[k]          = (FLOAT) 0.0;
	  for (k=0; k<ndim; k++) cell.r[k]       = (FLOAT) 0.0;
	  for (k=0; k<ndim; k++) cell.v[k]       = (FLOAT) 0.0;
	  //for (k=0; k<ndim; k++) cell.rcell[k]   = (FLOAT) 0.0;
	  for (k=0; k<ndim; k++) cell.bb.min[k]   = big_number;
	  for (k=0; k<ndim; k++) cell.bb.max[k]   = -big_number;
	  for (k=0; k<ndim; k++) cell.hbox.min[k] = big_number;
	  for (k=0; k<ndim; k++) cell.hbox.max[k] = -big_number;
	  for (k=0; k<ndim; k++) cell.vbox.min[k] = big_number;
	  for (k=0; k<ndim; k++) cell.vbox.max[k] = -big_number;
  }

  // If this is a leaf cell, sum over all particles
  //-----------------------------------------------------------------------------------------------
  if (cell.level == ltot && stock_leaf) {



    // Loop over all particles in cell summing their contributions
    for (i = cell.ifirst; i <= cell.ilast; ++i) {
      if (!partdata[i].flags.is_dead()) {
        cell.N++;
        if (partdata[i].flags.check(active)) cell.Nactive++;
        cell.hmax = max(cell.hmax,partdata[i].h);
        cell.maxsound = max(cell.maxsound, partdata[i].sound);
        if (gravmask[partdata[i].ptype]) {
          cell.m += partdata[i].m;
          for (k=0; k<ndim; k++) cell.r[k] += partdata[i].m*partdata[i].r[k];
          for (k=0; k<ndim; k++) cell.v[k] += partdata[i].m*partdata[i].v[k];
        }
        for (k=0; k<ndim; k++) {
         if (partdata[i].r[k] < cell.bb.min[k]) cell.bb.min[k] = partdata[i].r[k];
         if (partdata[i].r[k] > cell.bb.max[k]) cell.bb.max[k] = partdata[i].r[k];
         if (partdata[i].r[k] - kernrange*partdata[i].h < cell.hbox.min[k])
           cell.hbox.min[k] = partdata[i].r[k] - kernrange*partdata[i].h;
         if (partdata[i].r[k] + kernrange*partdata[i].h > cell.hbox.max[k])
           cell.hbox.max[k] = partdata[i].r[k] + kernrange*partdata[i].h;
         if (partdata[i].v[k] > cell.vbox.max[k]) cell.vbox.max[k] = partdata[i].v[k];
         if (partdata[i].v[k] < cell.vbox.min[k]) cell.vbox.min[k] = partdata[i].v[k];
        }
        if (gravity_mac == gadget2)
          cell.amin = min(cell.amin,
                          sqrt(DotProduct(partdata[i].atree,partdata[i].atree,ndim)));
        else if (gravity_mac == eigenmac)
          cell.macfactor = max(cell.macfactor,pow(partdata[i].gpot,-twothirds));
      }
    }

    // Normalise all cell values
    if (cell.N > 0) {
      for (k=0; k<ndim; k++) cell.r[k] /= cell.m;
      for (k=0; k<ndim; k++) cell.v[k] /= cell.m;
      //for (k=0; k<ndim; k++) cell.rcell[k] = (FLOAT) 0.5*(cell.bb.min[k] + cell.bb.max[k]);
      for (k=0; k<ndim; k++) dr[k] = (FLOAT) 0.5*(cell.bb.max[k] - cell.bb.min[k]);
      cell.cdistsqd = max(DotProduct(dr,dr,ndim),cell.hmax*cell.hmax)/thetamaxsqd;
      cell.rmax = sqrt(DotProduct(dr,dr,ndim));
    }

    // Compute quadrupole moment terms if selected
    if (need_quadrupole_moments) {

      for (i = cell.ifirst; i <= cell.ilast; ++i) {
        if (!partdata[i].flags.is_dead() && gravmask[partdata[i].ptype]) {
          mi = partdata[i].m;
          for (k=0; k<ndim; k++) dr[k] = partdata[i].r[k] - cell.r[k];
          drsqd = DotProduct(dr,dr,ndim);
          if (ndim == 3) {
            cell.q[0] += mi*((FLOAT) 3.0*dr[0]*dr[0] - drsqd);
            cell.q[1] += mi*(FLOAT) 3.0*dr[0]*dr[1];
            cell.q[2] += mi*((FLOAT) 3.0*dr[1]*dr[1] - drsqd);
            cell.q[3] += mi*(FLOAT) 3.0*dr[2]*dr[0];
            cell.q[4] += mi*(FLOAT) 3.0*dr[2]*dr[1];
          }
          else if (ndim == 2) {
            cell.q[0] += mi*((FLOAT) 3.0*dr[0]*dr[0] - drsqd);
            cell.q[1] += mi*(FLOAT) 3.0*dr[0]*dr[1];
            cell.q[2] += mi*((FLOAT) 3.0*dr[1]*dr[1] - drsqd);
          } else if (ndim == 1) {
            cell.q[0] += mi*((FLOAT) 3.0*dr[0]*dr[0] - drsqd);
          }
        }
      }
    }

  }
  // For non-leaf cells, sum together two children cells
  //-----------------------------------------------------------------------------------------------
  else if (cell.copen != -1) {

	  TreeCell<ndim> &child1 = celldata[cell.copen];
	  TreeCell<ndim> &child2 = celldata[child1.cnext];

    if (child1.N > 0) {
      for (k=0; k<ndim; k++) cell.bb.min[k] = min(child1.bb.min[k],cell.bb.min[k]);
      for (k=0; k<ndim; k++) cell.bb.max[k] = max(child1.bb.max[k],cell.bb.max[k]);
      for (k=0; k<ndim; k++) cell.hbox.min[k] = min(child1.hbox.min[k],cell.hbox.min[k]);
      for (k=0; k<ndim; k++) cell.hbox.max[k] = max(child1.hbox.max[k],cell.hbox.max[k]);
      for (k=0; k<ndim; k++) cell.vbox.min[k] = min(child1.vbox.min[k],cell.vbox.min[k]);
      for (k=0; k<ndim; k++) cell.vbox.max[k] = max(child1.vbox.max[k],cell.vbox.max[k]);
      cell.hmax = max(cell.hmax,child1.hmax);
      cell.maxsound = max(cell.maxsound,child1.maxsound);
      if (gravity_mac == gadget2)
        cell.amin = min(cell.amin, child1.amin);
      else if (gravity_mac == eigenmac)
        cell.macfactor = max(cell.macfactor, child1.macfactor) ;
    }
    if (child2.N > 0) {
      for (k=0; k<ndim; k++) cell.bb.min[k] = min(child2.bb.min[k],cell.bb.min[k]);
      for (k=0; k<ndim; k++) cell.bb.max[k] = max(child2.bb.max[k],cell.bb.max[k]);
      for (k=0; k<ndim; k++) cell.hbox.min[k] = min(child2.hbox.min[k],cell.hbox.min[k]);
      for (k=0; k<ndim; k++) cell.hbox.max[k] = max(child2.hbox.max[k],cell.hbox.max[k]);
      for (k=0; k<ndim; k++) cell.vbox.min[k] = min(child2.vbox.min[k],cell.vbox.min[k]);
      for (k=0; k<ndim; k++) cell.vbox.max[k] = max(child2.vbox.max[k],cell.vbox.max[k]);
      cell.hmax = max(cell.hmax,child2.hmax);
      cell.maxsound = max(cell.maxsound,child2.maxsound);
      if (gravity_mac == gadget2)
        cell.amin = min(cell.amin, child2.amin);
      else if (gravity_mac == eigenmac)
        cell.macfactor = max(cell.macfactor, child2.macfactor) ;
    }

    cell.N = child1.N + child2.N;
    cell.Nactive = child1.Nactive + child2.Nactive;
    if (cell.N > 0) {
      cell.m = child1.m + child2.m;
      for (k=0; k<ndim; k++) cell.r[k] = (child1.m*child1.r[k] + child2.m*child2.r[k])/cell.m;
      for (k=0; k<ndim; k++) cell.v[k] = (child1.m*child1.v[k] + child2.m*child2.v[k])/cell.m;
      //for (k=0; k<ndim; k++) cell.rcell[k] = (FLOAT) 0.5*(cell.bb.min[k] + cell.bb.max[k]);
      for (k=0; k<ndim; k++) dr[k] = (FLOAT) 0.5*(cell.bb.max[k] - cell.bb.min[k]);
      cell.cdistsqd = max(DotProduct(dr,dr,ndim),cell.hmax*cell.hmax)/thetamaxsqd;
      cell.rmax = sqrt(DotProduct(dr,dr,ndim));
#ifdef MPI_PARALLEL
      cell.worktot = child1.worktot + child2.worktot;
#endif
    }

    // Now add individual quadrupole moment terms
    if (need_quadrupole_moments && child1.N > 0) {
      mi = child1.m;
      for (k=0; k<ndim; k++) dr[k] = child1.r[k] - cell.r[k];
      drsqd = DotProduct(dr,dr,ndim);
      if (ndim == 3) {
        for (k=0; k<5; k++) cell.q[k] += child1.q[k] ;
        cell.q[0] += mi*((FLOAT) 3.0*dr[0]*dr[0] - drsqd);
        cell.q[1] += mi*(FLOAT) 3.0*dr[0]*dr[1];
        cell.q[2] += mi*((FLOAT) 3.0*dr[1]*dr[1] - drsqd);
        cell.q[3] += mi*(FLOAT) 3.0*dr[2]*dr[0];
        cell.q[4] += mi*(FLOAT) 3.0*dr[2]*dr[1];
      }
      else if (ndim == 2) {
        for (k=0; k<3; k++) cell.q[k] += child1.q[k] ;
        cell.q[0] += mi*((FLOAT) 3.0*dr[0]*dr[0] - drsqd);
        cell.q[1] += mi*(FLOAT) 3.0*dr[0]*dr[1];
        cell.q[2] += mi*((FLOAT) 3.0*dr[1]*dr[1] - drsqd);
      }
      else if (ndim == 1) {
        cell.q[0] += child1.q[0] ;
        cell.q[0] += mi*((FLOAT) 3.0*dr[0]*dr[0] - drsqd);
      }
    }

    if (need_quadrupole_moments && child2.N > 0) {
      mi = child2.m;
      for (k=0; k<ndim; k++) dr[k] = child2.r[k] - cell.r[k];
      drsqd = DotProduct(dr,dr,ndim);
      if (ndim == 3) {
        for (k=0; k<5; k++) cell.q[k] += child2.q[k] ;
        cell.q[0] += mi*((FLOAT) 3.0*dr[0]*dr[0] - drsqd);
        cell.q[1] += mi*(FLOAT) 3.0*dr[0]*dr[1];
        cell.q[2] += mi*((FLOAT) 3.0*dr[1]*dr[1] - drsqd);
        cell.q[3] += mi*(FLOAT) 3.0*dr[2]*dr[0];
        cell.q[4] += mi*(FLOAT) 3.0*dr[2]*dr[1];
      }
      else if (ndim == 2) {
        for (k=0; k<3; k++) cell.q[k] += child2.q[k] ;
        cell.q[0] += mi*((FLOAT) 3.0*dr[0]*dr[0] - drsqd);
        cell.q[1] += mi*(FLOAT) 3.0*dr[0]*dr[1];
        cell.q[2] += mi*((FLOAT) 3.0*dr[1]*dr[1] - drsqd);
      }
      else if (ndim == 1) {
        cell.q[0] += child2.q[0] ;
        cell.q[0] += mi*((FLOAT) 3.0*dr[0]*dr[0] - drsqd);
      }
    }

  }
  //-----------------------------------------------------------------------------------------------


  // Calculate eigenvalue MAC criteria
  if (gravity_mac == eigenmac) {
    if (ndim == 3) {
      p = cell.q[0]*cell.q[2] - (cell.q[0] + cell.q[2])*(cell.q[0] + cell.q[2]) -
          cell.q[1]*cell.q[1] - cell.q[3]*cell.q[3] - cell.q[4]*cell.q[4];
      if (p >= (FLOAT) 0.0) {
        lambda = 0;
      } else {
        lambda = (FLOAT) 2.0*sqrt(-p/(FLOAT) 3.0);
      }
    } else if (ndim == 2) {
      p = (cell.q[0]-cell.q[2])*(cell.q[0]-cell.q[2]) + 4*cell.q[1]*cell.q[1];
      lambda = (FLOAT) 0.5*max(cell.q[0] + cell.q[2] + sqrt(p), (FLOAT) 0.0);
    } else {
      lambda = fabs(cell.q[0]) ;
    }

    cell.mac = pow((FLOAT) 0.5*lambda/macerror,(FLOAT) 0.66666666666666);

  }
  else {
    cell.mac = (FLOAT) 0.0;
  }


  return;
}



//=================================================================================================
//  KDTree::ExtrapolateCellProperties
/// Extrapolate important physical properties of all cells in the tree.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void KDTree<ndim,ParticleType,TreeCell>::ExtrapolateCellProperties
 (const FLOAT dt)                            ///< Smallest timestep size
{
  int c;                               // Cell counter
  int k;                               // Dimension counter

  debug2("[KDTree::ExtrapolateCellProperties]");


  // ..
  //-----------------------------------------------------------------------------------------------
  for (c=0; c<Ncell; c++) {

    for (k=0; k<ndim; k++) celldata[c].r[k] += celldata[c].v[k]*dt;
    //for (k=0; k<ndim; k++) celldata[c].rcell[k] += celldata[c].v[k]*dt;
    for (k=0; k<ndim; k++) celldata[c].bb.min[k] += celldata[c].v[k]*dt;
    for (k=0; k<ndim; k++) celldata[c].bb.max[k] += celldata[c].v[k]*dt;
    for (k=0; k<ndim; k++) celldata[c].hbox.min[k] += celldata[c].v[k]*dt;
    for (k=0; k<ndim; k++) celldata[c].hbox.max[k] += celldata[c].v[k]*dt;

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
  int cc,ccc;                          // Cell counters
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
  cell.hmax = (FLOAT) 0.0;
  for (k=0; k<ndim; k++) cell.hbox.min[k] = big_number;
  for (k=0; k<ndim; k++) cell.hbox.max[k] = -big_number;


  // If this is a leaf cell, sum over all particles
  //-----------------------------------------------------------------------------------------------
  if (cell.level == ltot) {

    // Loop over all particles in cell summing their contributions
    for (i = cell.ifirst; i <= cell.ilast; ++i) {
      cell.hmax = max(cell.hmax,partdata[i].h);
      for (k=0; k<ndim; k++) {
        if (partdata[i].r[k] - kernrange*partdata[i].h < cell.hbox.min[k]) {
          cell.hbox.min[k] = partdata[i].r[k] - kernrange*partdata[i].h;
        }
        if (partdata[i].r[k] + kernrange*partdata[i].h > cell.hbox.max[k]) {
          cell.hbox.max[k] = partdata[i].r[k] + kernrange*partdata[i].h;
        }
      }
    };

  }
  // For non-leaf cells, sum together two children cells
  //-----------------------------------------------------------------------------------------------
  else {
    cc = cell.c1;
    ccc = cell.c2;
    if (celldata[cc].N > 0) {
      cell.hmax = max(cell.hmax,celldata[cc].hmax);
      for (k=0; k<ndim; k++) cell.hbox.min[k] = min(celldata[cc].hbox.min[k],cell.hbox.min[k]);
      for (k=0; k<ndim; k++) cell.hbox.max[k] = max(celldata[cc].hbox.max[k],cell.hbox.max[k]);
    }
    if (celldata[ccc].N > 0) {
      cell.hmax = max(cell.hmax,celldata[ccc].hmax);
      for (k=0; k<ndim; k++) cell.hbox.min[k] = min(celldata[ccc].hbox.min[k],cell.hbox.min[k]);
      for (k=0; k<ndim; k++) cell.hbox.max[k] = max(celldata[ccc].hbox.max[k],cell.hbox.max[k]);
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
 (Particle<ndim> *part_gen)        ///< [in] Main particle array
{
  int c;                               // Cell counter
  int i;                               // SPH particle index
  int ilast;                           // Last particle in linked list

  debug2("[KDTree::UpdateActiveParticleCounters]");
  //timing->StartTimingSection("TREE_UPDATE_COUNTERS");

  ParticleType<ndim>* partdata = reinterpret_cast<ParticleType<ndim>*>(part_gen) ;


  // Loop through all grid cells in turn
  //-----------------------------------------------------------------------------------------------
#pragma omp parallel for default(none) private(c,i,ilast) shared(partdata)
  for (c=0; c<Ncell; c++) {
    celldata[c].Nactive = 0;

    if (celldata[c].level != ltot) continue;
    i = celldata[c].ifirst;
    ilast = celldata[c].ilast;

    // Else walk through linked list to obtain list and number of active ptcls.
    for (; i <= ilast; ++i) {
      if (i < Ntot && partdata[i].flags.check(active) && !partdata[i].flags.is_dead())
        celldata[c].Nactive++;
    };

  }
  //-----------------------------------------------------------------------------------------------

  //timing->EndTimingSection("TREE_UPDATE_COUNTERS");

  return;
}



#ifdef MPI_PARALLEL
//=================================================================================================
//  KDTree::UpdateWorkCounters
/// Calculate the physical properties (e.g. total mass, centre-of-mass,
/// opening-distance, etc..) of all cells in the tree.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void KDTree<ndim,ParticleType,TreeCell>::UpdateWorkCounters
 (TreeCell<ndim> &cell)                ///< KD-tree cell
{
  int cc,ccc;                          // Cell counters
  int i;                               // Particle counter

  // If cell is not leaf, stock child cells
  //-----------------------------------------------------------------------------------------------
  if (cell.level != ltot && cell.c1 >= 0) {
#if defined _OPENMP
    if (pow(2,cell.level) < Nthreads) {
#pragma omp parallel for default(none) private(i) shared(cell) num_threads(2)
      for (i=0; i<2; i++) {
        if (i == 0) UpdateWorkCounters(celldata[cell.c1]);
        else if (i == 1) UpdateWorkCounters(celldata[cell.c2]);
      }
    }
    else {
      for (i=0; i<2; i++) {
        if (i == 0) UpdateWorkCounters(celldata[cell.c1]);
        else if (i == 1) UpdateWorkCounters(celldata[cell.c2]);
      }
    }
#else
    for (i=0; i<2; i++) {
      if (i == 0) UpdateWorkCounters(celldata[cell.c1]);
      else if (i == 1) UpdateWorkCounters(celldata[cell.c2]);
    }
#endif
  }
  //-----------------------------------------------------------------------------------------------


  // If this is a leaf cell, sum over all particles
  //-----------------------------------------------------------------------------------------------
  if (cell.level < ltot) {
    cell.worktot = (FLOAT) 0.0;
    cc = cell.c1;
    ccc = cell.c2;

    if (celldata[cc].N > 0) {
      cell.worktot += celldata[cc].worktot;
    }
    if (celldata[ccc].N > 0) {
      cell.worktot += celldata[ccc].worktot;
    }

  }
  //-----------------------------------------------------------------------------------------------

  //cout << "Work in cell " << cell.id << " : " << cell.worktot << "       level : " << cell.level << endl;

  return;
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
  bool kill_flag = false;           // ..
  int activecount;                  // Active particles in leaf cell
  int c;                            // Cell counter
  int cc;                           // Aux. cell counter
  int i;                            // Particle counter
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
    if (celldata[c].copen != -1) c = celldata[c].copen;
    else c = celldata[c].cnext;
  }

  // Now check we enter all cells once and once only
  for (c=0; c<Ncell; c++) {
    if (ccount[c] != 1) {
      cout << "Error in cell walk count : " << ccount[c] << endl;
      PrintArray("ccount     : ",Ncell,ccount);
      ExceptionHandler::getIstance().raise("Error in cell walk count in KDTree");
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
        for (i = cell.ifirst; i <= cell.ilast; ++i) {
          pcount[i]++;
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
      for (i = cell.ifirst; i <= cell.ilast; ++i) {
        pcount[i]++;
        leafcount++;
        Ncount++;
        if (partdata[i].flags.check(active)) activecount++;
        if (partdata[i].flags.check(active)) Nactivecount++;
        if (partdata[i].h > cell.hmax) {
          cout << "hmax flag error : " << c << "    "
               << partdata[i].h << "   " << cell.hmax << endl;
          ExceptionHandler::getIstance().raise("hmax flag error in KDTree");
        }
      }
      if (leafcount > Nleafmax) {
        ExceptionHandler::getIstance().raise("Error : leaf particle error in KDTree");
      }
      if (activecount > leafcount) {
        ExceptionHandler::getIstance().raise("Error : leaf particle error in KDTree");
      }
    }


    // Check that child cell boxes fit inside parent box
    if (cell.level != ltot) {
      KDTreeCell<ndim>& c1 = celldata[cell.c1];
      KDTreeCell<ndim>& c2 = celldata[cell.c2];

      for (int k=0; k<ndim; k++) {
        if (c1.hbox.min[k] < cell.hbox.min[k] || c1.hbox.max[k] > cell.hbox.max[k] ||
            c2.hbox.min[k] < cell.hbox.min[k] || c2.hbox.max[k] > cell.hbox.max[k])
            overlap_flag = true;
      }
      if (overlap_flag) {
        cout << "Parent/daughter cell overlap error!! : "
             << c << "   " << cell.c1 << "    " << cell.c2 << endl;
        ExceptionHandler::getIstance().raise("Parent/daughter cell overlap error in KDTree");
      }
    }

    // Check that bounding boxes of cells on each level do not overlap each other
    for (cc=0; cc<Ncell; cc++) {
      if (c != cc && celldata[cc].level == cell.level) {
        if (ndim == 1) {
          if (cell.bb.min[0] < celldata[cc].bb.max[0] &&
              cell.bb.max[0] > celldata[cc].bb.min[0]) {
            overlap_flag = true;
          }
        }
        else if (ndim == 2) {
          if (cell.bb.min[0] < celldata[cc].bb.max[0] &&
              cell.bb.max[0] > celldata[cc].bb.min[0] &&
              cell.bb.min[1] < celldata[cc].bb.max[1] &&
              cell.bb.max[1] > celldata[cc].bb.min[1]) {
            overlap_flag = true;
          }
        }
        else if (ndim == 3) {
          if (cell.bb.min[0] < celldata[cc].bb.max[0] &&
              cell.bb.max[0] > celldata[cc].bb.min[0] &&
              cell.bb.min[1] < celldata[cc].bb.max[1] &&
              cell.bb.max[1] > celldata[cc].bb.min[1] &&
              cell.bb.min[2] < celldata[cc].bb.max[2] &&
              cell.bb.max[2] > celldata[cc].bb.min[2]) {
            overlap_flag = true;
          }
        }
      }
      if (overlap_flag) {
        cout << "Brother/sister cell overlap error!! : " << c << "   " << cc << endl;
        ExceptionHandler::getIstance().raise("Brother/sister cell overlap error in KDTree");
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

  if (kill_flag) {
    ExceptionHandler::getIstance().raise("Kill flag set in KDTree");
  }

  return;
}
#endif




template class KDTree<1,GradhSphParticle,KDTreeCell>;
template class KDTree<2,GradhSphParticle,KDTreeCell>;
template class KDTree<3,GradhSphParticle,KDTreeCell>;
template class KDTree<1,SM2012SphParticle,KDTreeCell>;
template class KDTree<2,SM2012SphParticle,KDTreeCell>;
template class KDTree<3,SM2012SphParticle,KDTreeCell>;
template class KDTree<1,MeshlessFVParticle,KDTreeCell>;
template class KDTree<2,MeshlessFVParticle,KDTreeCell>;
template class KDTree<3,MeshlessFVParticle,KDTreeCell>;
