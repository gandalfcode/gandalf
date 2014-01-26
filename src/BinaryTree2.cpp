//=============================================================================
//  BinaryTree2.cpp
//  Contains all functions for building, stocking and walking for the 
//  binary tree for SPH particles.  
//  Based on Fortran code courtesy of O. Lomax.
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


#include <cstdlib>
#include <cassert>
#include <iostream>
#include <string>
#include <math.h>
#include "Precision.h"
#include "Exception.h"
#include "SphNeighbourSearch.h"
#include "Sph.h"
#include "Parameters.h"
#include "InlineFuncs.h"
#include "SphParticle.h"
#include "Debug.h"
#if defined _OPENMP
#include <omp.h>
#endif
using namespace std;



//=============================================================================
//  BinaryTree::BinaryTree
/// BinaryTree constructor.  Initialises various variables.
//=============================================================================
template <int ndim>
BinaryTree2<ndim>::BinaryTree2(int Nleafmaxaux, FLOAT thetamaxsqdaux, 
                               FLOAT kernrangeaux, string gravity_mac_aux,
                               string multipole_aux, int Nthreads, int Nmpi)
{
  allocated_tree = false;
  created_sub_trees = false;
  neibcheck = true;
  Ntot = 0;
  Ntotmax = 0;
  Ntotmaxold = 0;
  Nleafmax = Nleafmaxaux;
  kernrange = kernrangeaux;
  thetamaxsqd = thetamaxsqdaux;
  gravity_mac = gravity_mac_aux;
  multipole = multipole_aux;
}



//=============================================================================
//  BinaryTree2::~BinaryTree2
/// BinaryTree2 destructor.  Deallocates tree memory upon object destruction.
//=============================================================================
template <int ndim>
BinaryTree2<ndim>::~BinaryTree2()
{
  if (allocated_tree) DeallocateTreeMemory();
}



//=============================================================================
//  BinaryTree2::AllocateTreeMemory
/// Allocate memory for binary tree as requested.  If more memory is required 
/// than currently allocated, tree is deallocated and reallocated here.
//=============================================================================
template <int ndim>
void BinaryTree2<ndim>::AllocateTreeMemory(void)
{
  debug2("[BinaryTree2::AllocateTreeMemory]");

  if (!allocated_tree || Ntotmax > Ntotmaxold) {
    if (allocated_tree) DeallocateTreeMemory();
    Ntotmax = max(Ntotmax,Ntot);
    Ntotmaxold = Ntotmax;

    g2c = new int[gtot];
    ids = new int[Ntotmax];
    inext = new int[Ntotmax];
    tree = new struct BinaryTree2Cell<ndim>[Ncellmax];
    allocated_tree = true;
  }

  return;
}



//=============================================================================
//  BinaryTree2::DeallocateTreeMemory
/// Deallocates all binary tree memory
//=============================================================================
template <int ndim>
void BinaryTree2<ndim>::DeallocateTreeMemory(void)
{
  debug2("[BinaryTree2::DeallocateTreeMemory]");

  if (allocated_tree) {
    delete[] tree;
    delete[] inext;
    delete[] ids;
    delete[] g2c;
    allocated_tree = false;
  }

  return;
}



//=============================================================================
//  BinaryTree2::BuildTree
/// Call all routines to build/re-build the binary tree on the local node.
/// If OpenMP is activated, the local domain is partitioned into sub-trees 
/// in order to improve the scalability of building and stocking the tree.
//=============================================================================
template <int ndim>
void BinaryTree2<ndim>::BuildTree
(bool rebuild_tree,                 ///< Flag to rebuild tree
 int n,                             ///< Integer time
 int ntreebuildstep,                ///< Tree build frequency
 int ntreestockstep,                ///< Tree stocking frequency
 FLOAT timestep,                    ///< Smallest physical timestep
 Sph<ndim> *sph)                    ///< Pointer to SPH object
{
  debug2("[BinaryTree2::BuildTree]");

  // For tree rebuild steps
  //---------------------------------------------------------------------------
  if (n%ntreebuildstep == 0 || rebuild_tree) {

    // Set no. of tree members to total number of SPH particles (inc. ghosts)
    Nsph = sph->Nsph;
    Ntot = sph->Ntot;
    Ntotmax = max(Ntot,Ntotmax);
    gtot = 0;

    // Compute the size of all tree-related arrays now we know number of points
    ComputeTreeSize();

    // Allocate (or reallocate if needed) all tree memory
    AllocateTreeMemory();

    // Create tree data structure including linked lists and cell pointers
    CreateTreeStructure();

    // Now add particles to tree depending on Cartesian coordinates
    LoadParticlesToTree(sph);

    // Calculate all cell quantities (e.g. COM, opening distance)
    StockTree(tree[0],sph->sphdata);

#if defined(VERIFY_ALL)
    ValidateTree(sph);
#endif
    
  }

  // Else stock the tree
  //---------------------------------------------------------------------------
  else if (n%ntreestockstep == 0) {

    StockTree(tree[0],sph->sphdata);

  }

  // Otherwise simply extrapolate tree cell properties
  //---------------------------------------------------------------------------
  else {

    //ExtrapolateCellProperties(tree[0],timestep);
    ExtrapolateCellProperties(timestep);

  }
  //---------------------------------------------------------------------------


  return;
}



//=============================================================================
//  BinaryTree2::UpdateActiveParticleCounters
/// Loop through all leaf cells in binary tree and update all active 
/// particle counters.
//=============================================================================
template <int ndim>
void BinaryTree2<ndim>::UpdateActiveParticleCounters
(Sph<ndim> *sph)                    ///< Pointer to main SPH object
{
  return;
}



//=============================================================================
//  BinaryTree2::ComputeTreeSize
/// Compute the maximum size (i.e. no. of levels, cells and leaf cells) of 
/// the binary tree.
//=============================================================================
template <int ndim>
void BinaryTree2<ndim>::ComputeTreeSize(void)
{
  debug2("[BinaryTree2::ComputeTreeSize]");

  // Increase level until tree can contain all particles
  ltot = 0;
  while (Nleafmax*pow(2,ltot) < Ntotmax) {
    ltot++;
  };

  // Set total number of leaf/grid cells and tree cells
  gtot = pow(2,ltot);
  Ncell = 2*gtot - 1;
  Ncellmax = Ncell;

  // Optional output (for debugging)
#if defined(VERIFY_ALL)
  cout << "Calculating tree size variables" << endl;
  cout << "No. of grid-cells     : " << gtot << endl;
  cout << "No. of levels on tree : " << ltot << endl;
  cout << "No. of cells in tree  : " << Ncell << endl; 
#endif

  return;
}



//=============================================================================
//  BinaryTree2::CreateTreeStructure
/// Create the raw tree skeleton structure once the tree size is known.
/// Sets all cell pointer variables and all cell levels.
//=============================================================================
template <int ndim>
void BinaryTree2<ndim>::CreateTreeStructure(void)
{
  int c;                            // Dummy id of tree-level, then tree-cell
  int g;                            // Dummy id of grid-cell
  int l;                            // Dummy id of level
  int *c2L;                         // Increment to second child-cell
  int *cNL;                         // Increment to next cell if cell unopened

  debug2("[BinaryTree2::CreateTreeStructure]");

  // Allocate memory for local arrays
  c2L = new int[ltot + 1];
  cNL = new int[ltot + 1];

  // Set pointers to second child-cell (if opened) and next cell (if unopened)
  for (l=0; l<ltot; l++) {
    c2L[l] = pow(2,ltot - l);
    cNL[l] = 2*c2L[l] - 1;
  }

  // Zero tree cell variables
  for (g=0; g<gtot; g++) g2c[g] = 0;
  for (c=0; c<Ncell; c++) {
    tree[c].c2g = 0;
    tree[c].c1 = -1;
    tree[c].c2 = -1;
    tree[c].ifirst = -1;
    tree[c].ilast = -1;
    tree[c].N = 0;
  }
  g = 0;
  tree[0].level = 0;

  // Loop over all cells and set all other pointers
  //---------------------------------------------------------------------------
  for (c=0; c<Ncell; c++) {
    tree[c].id = c;
    if (tree[c].level == ltot) {                    // If on leaf level
      tree[c].cnext = c + 1;                        // id of next cell
      tree[c].c2g = g;                              // Record leaf id
      g2c[g++] = c;                                 // Record inverse id
    }
    else {
      tree[c+1].level = tree[c].level + 1;          // Level of 1st child
      tree[c].c1 = c + 1;
      tree[c].c2 = c + c2L[tree[c].level];          // id of 2nd child
      tree[tree[c].c2].level = tree[c].level + 1;   // Level of 2nd child
      tree[c].cnext = c + cNL[tree[c].level];       // Next cell id
    }

  }
  //---------------------------------------------------------------------------


  // Free locally allocated memory
  delete[] cNL;
  delete[] c2L;

  return;
}



//=============================================================================
//  BinaryTree2::LoadParticlesToTree
/// Create tree structure by adding particles to leaf cells.
//=============================================================================
template <int ndim>
void BinaryTree2<ndim>::LoadParticlesToTree
(Sph<ndim> *sph)                   ///< Pointer to main SPH object
{
  int i;                           // SPH particle counter
  int k;                           // Dimensionality counter

  debug2("[BinaryTree2::LoadParticleToTree]");

  // Set root cell counters and properties
  tree[0].N = Ntot;
  tree[0].ifirst = 0;
  tree[0].ilast = Ntot - 1;
  for (k=0; k<ndim; k++) tree[0].bbmin[k] = box->boxmin[k];
  for (k=0; k<ndim; k++) tree[0].bbmax[k] = box->boxmax[k];

  for (i=0; i<Ntot; i++) ids[i] = i;
  for (i=0; i<Ntot-1; i++) inext[i] = -1;
  
  // Recursively build tree from root node down
  DivideTreeCell(0,Ntot-1,sph,tree[0]);

  return;
}



//=============================================================================
//  BinaryTree2::DivideTreeCell
/// Recursive routine to divide a tree cell into two children cells.
//=============================================================================
template <int ndim>
void BinaryTree2<ndim>::DivideTreeCell
(int ifirst,                        ///< Aux. id of first particle in cell
 int ilast,                         ///< Aux. id of last particle in cell
 Sph<ndim> *sph,                    ///< Pointer to main SPH object
 BinaryTree2Cell<ndim> &cell)       ///< Cell to be divided
{
  int i;                            // Aux. child cell counter
  int j;                            // Aux. particle counter
  int k;                            // Dimension counter
  int k_divide = 0;                 // Division dimension
  FLOAT rkmax = 0.0;                // Max. box size of all dimensions
  FLOAT rdivide;                    // Coordinate value at division


  // If cell is a leaf cell, do not divide further and set linked lists
  if (cell.level == ltot) {
    if (cell.N > 0) {
      for (j=cell.ifirst; j<cell.ilast; j++) inext[ids[j]] = ids[j+1];
      inext[ids[cell.ilast]] = -1;
      cell.ifirst = ids[cell.ifirst];
      cell.ilast = ids[cell.ilast];
    }
    return;
  }

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
  rdivide = QuickSelect(cell.ifirst,cell.ilast,
			cell.ifirst+cell.N/2,k_divide,sph);


  // Set properties of first child cell
  for (k=0; k<ndim; k++) tree[cell.c1].bbmin[k] = cell.bbmin[k];
  for (k=0; k<ndim; k++) tree[cell.c1].bbmax[k] = cell.bbmax[k];
  tree[cell.c1].bbmax[k_divide] = rdivide;
  tree[cell.c1].N = cell.N/2;
  if (tree[cell.c1].N != 0) {
    tree[cell.c1].ifirst = ifirst;
    tree[cell.c1].ilast = ifirst + cell.N/2 - 1;
  }

  // Set properties of second child cell
  for (k=0; k<ndim; k++) tree[cell.c2].bbmin[k] = cell.bbmin[k];
  for (k=0; k<ndim; k++) tree[cell.c2].bbmax[k] = cell.bbmax[k];
  tree[cell.c2].bbmin[k_divide] = rdivide;
  tree[cell.c2].N = cell.N - tree[cell.c1].N;
  if (tree[cell.c2].N != 0) {
    tree[cell.c2].ifirst = ifirst + cell.N/2;
    tree[cell.c2].ilast = ilast;
  }


  // Now divide the new child cells as a recursive function
#if defined _OPENMP
  if (pow(2,cell.level) <= omp_get_max_threads()) {
#pragma omp parallel for default(none) private(i) shared(ifirst,ilast,sph) num_threads(2)
    for (i=0; i<2; i++) {
      if (i == 0) DivideTreeCell(ifirst,ifirst+cell.N/2-1,sph,tree[cell.c1]);
      else if (i == 1) DivideTreeCell(ifirst+cell.N/2,ilast,sph,tree[cell.c2]);
    }
  }
  else {
    for (i=0; i<2; i++) {
      if (i == 0) DivideTreeCell(ifirst,ifirst+cell.N/2-1,sph,tree[cell.c1]);
      else if (i == 1) DivideTreeCell(ifirst+cell.N/2,ilast,sph,tree[cell.c2]);
    }
  }
#else
  for (i=0; i<2; i++) {
    if (i == 0) DivideTreeCell(ifirst,ifirst+cell.N/2-1,sph,tree[cell.c1]);
    else if (i == 1) DivideTreeCell(ifirst+cell.N/2,ilast,sph,tree[cell.c2]);
  }
#endif

  return;
}



//=============================================================================
//  BinaryTree2::QuickSelect
/// Find median and sort particles in arrays to ensure they are the correct 
/// side of the division.  Uses the QuickSelect algorithm.
//=============================================================================
template <int ndim>
FLOAT BinaryTree2<ndim>::QuickSelect
(int left,                          ///< ..
 int right,                         ///< ..
 int jpivot,                        ///< ..
 int k,                             ///< ..
 Sph<ndim> *sph)                    ///< ..
{
  int i;                            // ..
  int j;                            // ..
  int jguess;                       // ..
  int jtemp;                        // ..
  FLOAT rpivot;                     // ..


  // ..
  //---------------------------------------------------------------------------
  do {

    // Make a guess of pivot value
    jguess = (left + right)/2;
    rpivot = sph->sphdata[ids[jguess]].r[k];

    // ..
    jtemp = ids[jguess];
    ids[jguess] = ids[right];
    ids[right] = jtemp;

    // ..
    jguess = left;

    //-------------------------------------------------------------------------
    for (j=left; j<right; j++) {

      if (sph->sphdata[ids[j]].r[k] <= rpivot) {
	jtemp = ids[j];
	ids[j] = ids[jguess];
	ids[jguess] = jtemp;
        jguess++;
      }

    }
    //-------------------------------------------------------------------------


    // Move ?? particle from end of array to index i
    jtemp = ids[right];
    ids[right] = ids[jguess];
    ids[jguess] = jtemp;

    // jguess is lower than jpivot.  
    // Only need to search between jguess+1 and right
    if (jguess < jpivot) left = jguess + 1;

    // jguess is higher than jpivot.  
    // Only need to search between left and jguess-1
    else if (jguess > jpivot) right = jguess - 1;

  } while (jguess != jpivot);
  //---------------------------------------------------------------------------


  return rpivot;
}



//=============================================================================
//  BinaryTree2::StockTree
/// ..
//=============================================================================
template <int ndim>
void BinaryTree2<ndim>::StockTree
(BinaryTree2Cell<ndim> &cell,       ///< Reference to cell to be stocked
 SphParticle<ndim> *sphdata)        ///< SPH particle data array
{
  int i;                            // Aux. child cell counter

  // If cell is not leaf, stock child cells
  if (cell.level != ltot) {
#if defined _OPENMP
    if (pow(2,cell.level) <= omp_get_max_threads()) {
#pragma omp parallel for default(none) private(i) shared(sphdata) num_threads(2)
      for (i=0; i<2; i++) {
	if (i == 0) StockTree(tree[cell.c1],sphdata);
	else if (i == 1) StockTree(tree[cell.c2],sphdata);
      }
    }
    else {
      for (i=0; i<2; i++) {
	if (i == 0) StockTree(tree[cell.c1],sphdata);
	else if (i == 1) StockTree(tree[cell.c2],sphdata);
      }
    }
#else
    for (i=0; i<2; i++) {
      if (i == 0) StockTree(tree[cell.c1],sphdata);
      else if (i == 1) StockTree(tree[cell.c2],sphdata);
    }  
#endif
  }

  // Stock node once all children are stocked
  StockCellProperties(cell,sphdata);

  return;
}



//=============================================================================
//  BinaryTree2::StockCellProperties
/// Calculate the physical properties (e.g. total mass, centre-of-mass, 
/// opening-distance, etc..) of all cells in the tree.
//=============================================================================
template <int ndim>
void BinaryTree2<ndim>::StockCellProperties
(BinaryTree2Cell<ndim> &cell,       ///< ..
 SphParticle<ndim> *sphdata)        ///< SPH particle data array
{
  int c,cc,ccc;                     // Cell counters
  int i;                            // Particle counter
  int j;                            // ..
  int k;                            // Dimension counter
  FLOAT dr[ndim];                   // Relative position vector
  FLOAT drsqd;                      // Distance squared
  FLOAT factor = 1.0/thetamaxsqd;   // Geometric MAC aux. variable
  FLOAT mi;                         // Mass of particle i
  BinaryTree2Cell<ndim> &child1 = tree[cell.c1];
  BinaryTree2Cell<ndim> &child2 = tree[cell.c2];


  // Zero all summation variables for all cells
  cell.Nactive = 0;
  cell.N = 0;
  cell.m = 0.0;
  cell.hmax = 0.0;
  cell.rmax = 0.0;
  cell.dhmaxdt = 0.0;
  cell.drmaxdt = 0.0;
  cell.cdistsqd = big_number;
  for (k=0; k<ndim; k++) cell.r[k] = 0.0;
  for (k=0; k<ndim; k++) cell.v[k] = 0.0;
  for (k=0; k<ndim; k++) cell.bbmin[k] = big_number;
  for (k=0; k<ndim; k++) cell.bbmax[k] = -big_number;
  for (k=0; k<5; k++) cell.q[k] = 0.0;


  // If this is a leaf cell, sum over all particles
  //---------------------------------------------------------------------------
  if (cell.level == ltot) {
    i = cell.ifirst;
    
    // Loop over all particles in cell summing their contributions
    while (i != -1) {
      cell.N++;
      if (sphdata[i].active) cell.Nactive++;
      cell.hmax = max(cell.hmax,sphdata[i].h);
      cell.m += sphdata[i].m;
      for (k=0; k<ndim; k++) cell.r[k] += sphdata[i].m*sphdata[i].r[k];
      for (k=0; k<ndim; k++) cell.v[k] += sphdata[i].m*sphdata[i].v[k];
      for (k=0; k<ndim; k++) {
	if (sphdata[i].r[k] < cell.bbmin[k]) cell.bbmin[k] = sphdata[i].r[k];
	if (sphdata[i].r[k] > cell.bbmax[k]) cell.bbmax[k] = sphdata[i].r[k];
      }
      if (i == cell.ilast) break;
      i = inext[i];
    };
    
    // Normalise all cell values
    if (cell.N > 0) {
      for (k=0; k<ndim; k++) cell.r[k] /= cell.m;
      for (k=0; k<ndim; k++) cell.v[k] /= cell.m;
      for (k=0; k<ndim; k++) dr[k] = 0.5*(cell.bbmax[k] - cell.bbmin[k]);
      cell.cdistsqd = factor*DotProduct(dr,dr,ndim);
      for (k=0; k<ndim; k++) 
	dr[k] = max(cell.bbmax[k] - cell.r[k],cell.r[k] - cell.bbmin[k]);
      cell.rmax = sqrt(DotProduct(dr,dr,ndim));
    }
    
    // Compute quadrupole moment terms if selected
    if (multipole == "quadrupole") {
      i = cell.ifirst;
      
      while (i != -1) {
	mi = sphdata[i].m;
	for (k=0; k<ndim; k++) dr[k] = sphdata[i].r[k] - cell.r[k];
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
	if (i == cell.ilast) break;
	i = inext[i];
      }
    }

  }
  // For non-leaf cells, sum together two children cells
  //---------------------------------------------------------------------------
  else {
    cell.N = child1.N + child2.N;
    
    if (child1.N > 0) {
      for (k=0; k<ndim; k++)
	cell.bbmin[k] = min(child1.bbmin[k],cell.bbmin[k]);
      for (k=0; k<ndim; k++)
	cell.bbmax[k] = max(child1.bbmax[k],cell.bbmax[k]);
    }
    if (child2.N > 0) {
      for (k=0; k<ndim; k++)
	cell.bbmin[k] = min(child2.bbmin[k],cell.bbmin[k]);
      for (k=0; k<ndim; k++)
	cell.bbmax[k] = max(child2.bbmax[k],cell.bbmax[k]);
    }


    if (cell.N > 0) {
      cell.hmax = max(child1.hmax,child2.hmax);
      cell.m = child1.m + child2.m;
      for (k=0; k<ndim; k++) cell.r[k] =
	(child1.m*child1.r[k] + child2.m*child2.r[k])/cell.m;
      for (k=0; k<ndim; k++) cell.v[k] =
	(child1.m*child1.v[k] + child2.m*child2.v[k])/cell.m;
      for (k=0; k<ndim; k++) 
	dr[k] = 0.5*(cell.bbmax[k] - cell.bbmin[k]);
      cell.cdistsqd = factor*DotProduct(dr,dr,ndim);
      for (k=0; k<ndim; k++) dr[k] = max(cell.bbmax[k] - cell.r[k],
					 cell.r[k] - cell.bbmin[k]);
      cell.rmax = sqrt(DotProduct(dr,dr,ndim));
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
  //---------------------------------------------------------------------------


  return;
}



//=============================================================================
//  BinarySubTree::ExtrapolateCellProperties
/// Extrapolate important physical properties of all cells in the tree.
//=============================================================================
template <int ndim>
void BinaryTree2<ndim>::ExtrapolateCellProperties
(FLOAT dt)                          ///< Smallest timestep size
{
  int c;                            // ..
  int k;                            // ..

  debug2("[BinaryTree2::ExtrapolateCellProperties]");


  // Loop backwards over all tree cells to ensure child cells are always
  // computed first before being summed in parent cells.
  //===========================================================================
  for (c=Ncell-1; c>=0; c--) {

    for (k=0; k<ndim; k++) tree[c].r[k] += tree[c].v[k]*dt;
    //tree[c].rmax += tree[c].drmaxdt*dt;
    //tree[c].hmax += tree[c].dhmaxdt*dt;

  }
  //===========================================================================

  return;
}



//=============================================================================
//  BinarySubTree::UpdateHmaxValues
/// Calculate the physical properties (e.g. total mass, centre-of-mass,
/// opening-distance, etc..) of all cells in the tree.
//=============================================================================
template <int ndim>
void BinaryTree2<ndim>::UpdateHmaxValues
(SphParticle<ndim> *sphdata)        ///< SPH particle data array
{
  int c,cc,ccc;                     // Cell counters
  int i;                            // Particle counter
  int j;

  debug2("[BinarySubTree::UpdateHmaxValues]");

  // Zero all summation variables for all cells
  for (c=0; c<Ncellmax; c++) tree[c].hmax = 0.0;

  // Loop backwards over all tree cells to ensure child cells are always 
  // computed first before being summed in parent cells.
  //===========================================================================
  for (c=Ncell-1; c>=0; c--) {

    // If this is a leaf cell, sum over all particles
    // ------------------------------------------------------------------------
    if (tree[c].level == ltot) {
      i = tree[c].ifirst;

      // Loop over all particles in cell summing their contributions
      while (i != -1) {
        tree[c].hmax = max(tree[c].hmax,sphdata[i].h);
	if (i == tree[c].ilast) break;
        i = inext[i];
      };

    }
    // For non-leaf cells, sum together two children cells
    //-------------------------------------------------------------------------
    else {
      cc = tree[c].c1;
      ccc = tree[c].c2;
      if (tree[cc].N > 0) tree[c].hmax = max(tree[c].hmax,tree[cc].hmax);
      if (tree[ccc].N > 0) tree[c].hmax = max(tree[c].hmax,tree[ccc].hmax);

    }
    //-------------------------------------------------------------------------

  }
  //===========================================================================

  return;
}



//=============================================================================
//  BinaryTree2::ComputeActiveParticleList
/// Returns the number (Nactive) and list of ids (activelist) of all active
/// SPH particles in the given cell.
//=============================================================================
template <int ndim>
int BinaryTree2<ndim>::ComputeActiveParticleList
(BinaryTree2Cell<ndim> *cell,       ///< [in] Pointer to cell
 Sph<ndim> *sph,                    ///< [in] SPH object pointer
 int *activelist)                   ///< [out] List of active particles in cell
{
  int i = cell->ifirst;             // Local particle id (set to first ptcl id)
  int ilast = cell->ilast;          // i.d. of last particle in cell c
  int Nactive = 0;                  // No. of active particles in cell

  // Walk through linked list to obtain list and number of active ptcls.
  while (i != -1) {
    if (i < sph->Nsph && sph->sphdata[i].active) activelist[Nactive++] = i;
    if (i == ilast) break;
    i = inext[i];
  };

  return Nactive;
}



//=============================================================================
//  BinaryTree2::ComputeGatherNeighbourList
/// Computes and returns number of neighbour, 'Nneib', and the list
/// of neighbour ids, 'neiblist', for all particles inside cell 'c'.
/// Includes all particles in the selected cell, plus all particles
/// contained in adjacent cells (including diagonal cells).
/// Wrapper around the true implementation inside BinarySubTree
//=============================================================================
template <int ndim>
int BinaryTree2<ndim>::ComputeGatherNeighbourList
(BinaryTree2Cell<ndim> *cell,       ///< [in] Pointer to current cell
 int Nneibmax,                      ///< [in] Max. no. of neighbours
 int *neiblist,                     ///< [out] List of neighbour i.d.s
 FLOAT hmax,                        ///< [in] Maximum smoothing length
 SphParticle<ndim> *sphdata)        ///< [in] SPH particle data
{
  int cc;                           // Cell counter
  int i;                            // Particle id
  int j;                            // Aux. particle counter
  int k;                            // Neighbour counter
  int Nneib = 0;                    // Neighbour counter
  int Ntemp = 0;                    // Temporary neighbour counter
  FLOAT dr[ndim];                   // Relative position vector
  FLOAT drsqd;                      // Distance squared
  FLOAT rc[ndim];                   // Position of cell
  FLOAT hrangemax;                  // Maximum SPH kernel extent

  for (k=0; k<ndim; k++) rc[k] = cell->r[k];
  hrangemax = cell->rmax + kernrange*hmax;

  // Start with root cell and walk through entire tree
  cc = 0;

  //===========================================================================
  while (cc < Ncell) {

    for (k=0; k<ndim; k++) dr[k] = tree[cc].r[k] - rc[k];
    drsqd = DotProduct(dr,dr,ndim);

    // Check if circular range overlaps cell bounding sphere
    //-------------------------------------------------------------------------
    if (drsqd < pow(tree[cc].rmax + hrangemax,2)) {

      // If not a leaf-cell, then open cell to first child cell
      if (tree[cc].level != ltot)
        cc++;

      // If leaf-cell, add particles to list
      else if (tree[cc].level == ltot && Nneib + Nleafmax < Nneibmax) {
        i = tree[cc].ifirst;
    	while (i != -1) {
          neiblist[Nneib++] = i;
          if (i == tree[cc].ilast) break;
    	  i = inext[i];
        };
        cc = tree[cc].cnext;
      }

      // If leaf-cell, but we've run out of memory, return with error-code (-1)
      else if (tree[cc].level == ltot && Nneib + Nleafmax >= Nneibmax)
    	return -1;

    }

    // If not in range, then open next cell
    //-------------------------------------------------------------------------
    else
      cc = tree[cc].cnext;

  };
  //===========================================================================


  // Now, trim the list to remove particles that are definitely not neighbours
  hrangemax = hrangemax*hrangemax;
  for (j=0; j<Nneib; j++) {
    i = neiblist[j];
    for (k=0; k<ndim; k++) dr[k] = sphdata[i].r[k] - rc[k];
    drsqd = DotProduct(dr,dr,ndim);
    if (drsqd < hrangemax) neiblist[Ntemp++] = i;
  }
  Nneib = Ntemp;


  return Nneib;
}



//=============================================================================
//  BinaryTree2::ComputeNeighbourList
/// Computes and returns number of neighbour, 'Nneib', and the list
/// of neighbour ids, 'neiblist', for all particles inside cell 'c'.
/// Includes all particles in the selected cell, plus all particles
/// contained in adjacent cells (including diagonal cells).
/// Wrapper around the true implementation inside BinarySubTree.
/// If allocated memory array containing neighbour ids (neiblist) overflows, 
/// return with error code (-1) in order to reallocate more memory.
//=============================================================================
template <int ndim>
int BinaryTree2<ndim>::ComputeNeighbourList
(BinaryTree2Cell<ndim> *cell,       ///< [in] Cell pointer
 int Nneibmax,                      ///< [in] Max. no. of neighbours
 int *neiblist,                     ///< [out] List of neighbour i.d.s
 SphParticle<ndim> *sphdata)        ///< [in] SPH particle data
{
  int cc;                           // Cell counter
  int i;                            // Particle id
  int j;                            // Aux. particle counter
  int k;                            // Neighbour counter
  int Nneib = 0;                    // No. of SPH neighbours
  int Ntemp = 0;                    // Temp. no. of neighbouts
  FLOAT dr[ndim];                   // Relative position vector
  FLOAT drsqd;                      // Distance squared
  FLOAT rc[ndim];                   // Position of cell
  FLOAT hrangemax;                  // Maximum kernel extent
  FLOAT rmax;                       // Max. extent of particles from cell COM


  for (k=0; k<ndim; k++) rc[k] = cell->r[k];
  hrangemax = cell->rmax + kernrange*cell->hmax;
  rmax = cell->rmax;

  // Start with root cell and walk through entire tree
  cc = 0;


  //===========================================================================
  while (cc < Ncell) {

    for (k=0; k<ndim; k++) dr[k] = tree[cc].r[k] - rc[k];
    drsqd = DotProduct(dr,dr,ndim);


    // Check if circular range overlaps cell bounding sphere
    //-------------------------------------------------------------------------
    if (drsqd < pow(tree[cc].rmax + hrangemax,2) ||
        drsqd < pow(rmax + tree[cc].rmax + kernrange*tree[cc].hmax,2)) {

      // If not a leaf-cell, then open cell to first child cell
      if (tree[cc].level != ltot)
        cc++;

      // If leaf-cell, add particles to list
      else if (tree[cc].level == ltot && Nneib + Nleafmax < Nneibmax) {
        i = tree[cc].ifirst;
    	while (i != -1) {
          neiblist[Nneib++] = i;
          if (i == tree[cc].ilast) break;
    	  i = inext[i];
        };
        cc = tree[cc].cnext;
      }

      // If leaf-cell, but we've run out of memory, return with error-code (-1)
      else if (tree[cc].level == ltot && Nneib + Nleafmax >= Nneibmax)
    	return -1;

    }

    // If not in range, then open next cell
    //-------------------------------------------------------------------------
    else
      cc = tree[cc].cnext;
  };
  //===========================================================================


  // Now, trim the list to remove particles that are definitely not neighbours
  hrangemax = hrangemax*hrangemax;
  for (j=0; j<Nneib; j++) {
    i = neiblist[j];
    for (k=0; k<ndim; k++) dr[k] = sphdata[i].r[k] - rc[k];
    drsqd = DotProduct(dr,dr,ndim);
    if (drsqd < hrangemax || drsqd <
        (rmax + kernrange*sphdata[i].h)*(rmax + kernrange*sphdata[i].h));
      neiblist[Ntemp++] = i;
  }
  Nneib = Ntemp;


  return Nneib;
}



//=============================================================================
//  BinarySubTree::ComputeGravityInteractionList
/// Computes and returns number of SPH neighbours (Nneib), direct sum particles
/// (Ndirect) and number of cells (Ngravcell), including lists of ids, from
/// the gravity tree walk for active particles inside cell c.
/// Currently defaults to the geometric opening criteria.
/// If any of the interactions list arrays (neiblist,directlist,gravcelllist) 
/// overflow, return with error code (-1) to reallocate more memory.
//=============================================================================
template <int ndim>
int BinaryTree2<ndim>::ComputeGravityInteractionList
(BinaryTree2Cell<ndim> *cell,       ///< [in] Pointer to cell
 int Nneibmax,                      ///< [in] Max. no. of SPH neighbours
 int Ndirectmax,                    ///< [in] Max. no. of direct-sum neighbours
 int Ngravcellmax,                  ///< [in] Max. no. of cell interactions
 int &Nneib,                        ///< [out] No. of SPH neighbours
 int &Ndirect,                      ///< [out] No. of direct-sum neighbours
 int &Ngravcell,                    ///< [out] No. of cell interactions
 int *neiblist,                     ///< [out] List of SPH neighbour ids
 int *directlist,                   ///< [out] List of direct-sum neighbour ids
 BinaryTree2Cell<ndim> **gravcelllist,  ///< [out] List of cell ids
 SphParticle<ndim> *sphdata)        ///< [in] SPH particle data
{
  int cc;                           // Cell counter
  int i;                            // Particle id
  int j;                            // Aux. particle counter
  int k;                            // Neighbour counter
  int Nneibtemp = Nneib;            // Aux. counter
  FLOAT cdistsqd;                   // ..
  FLOAT dr[ndim];                   // Relative position vector
  FLOAT drsqd;                      // Distance squared
  FLOAT rc[ndim];                   // Position of cell
  FLOAT hrangemax;                  // Maximum kernel extent 
  FLOAT rmax;                       // Radius of sphere containing particles

  // Make local copies of important cell properties
  for (k=0; k<ndim; k++) rc[k] = cell->r[k];
  hrangemax = cell->rmax + kernrange*cell->hmax;
  cdistsqd = cell->cdistsqd;
  rmax = cell->rmax;

  // Start with root cell and walk through entire tree
  cc = 0;
  Nneib = 0;
  Ndirect = 0;
  Ngravcell = 0;


  // Walk through all cells in tree to determine particle and cell 
  // interaction lists
  //===========================================================================
  while (cc < Ncell) {
    for (k=0; k<ndim; k++) dr[k] = tree[cc].r[k] - rc[k];
    drsqd = DotProduct(dr,dr,ndim);

    // Check if cells contain SPH neighbours
    //-------------------------------------------------------------------------
    if (drsqd < (tree[cc].rmax + hrangemax)*(tree[cc].rmax + hrangemax) ||
        drsqd < pow(tree[cc].rmax + rmax + kernrange*tree[cc].hmax,2)) {

      // If not a leaf-cell, then open cell to first child cell
      if (tree[cc].level != ltot)
        cc++;

      // If leaf-cell, add particles to list
      else if (tree[cc].level == ltot && Nneib + Nleafmax <= Nneibmax) {
        i = tree[cc].ifirst;
    	while (i != -1) {
          neiblist[Nneib++] = i;
          if (i == tree[cc].ilast) break;
    	  i = inext[i];
        };
        cc = tree[cc].cnext;
      }

      // If leaf-cell, but we've run out of memory, return with error-code (-1)
      else if (tree[cc].level == ltot && Nneib + Nleafmax > Nneibmax)
    	return -1;

    }

    // Check if cell is far enough away to use the COM approximation
    //-------------------------------------------------------------------------
    //else if (drsqd >= tree[c].cdistsqd + tree[cc].cdistsqd) {
    else if (drsqd > cdistsqd || drsqd > tree[cc].cdistsqd) {

      // If cell is a leaf-cell with only one particle, more efficient to
      // compute the gravitational contribution from the particle than the cell
      if (tree[cc].level == ltot && tree[cc].N == 1 && Ndirect < Ndirectmax)
        directlist[Ndirect++] = tree[cc].ifirst;
      else if (Ngravcell < Ngravcellmax)
        gravcelllist[Ngravcell++] = &(tree[cc]);
      else
        return -1;
      cc = tree[cc].cnext;

    }

    // If cell is too close, open cell to interogate children cells.
    // If cell is too close and a leaf cell, then add particles to direct list.
    //-------------------------------------------------------------------------
    //else if (drsqd < tree[c].cdistsqd + tree[cc].cdistsqd) {
    else if (drsqd <= cdistsqd && drsqd <= tree[cc].cdistsqd) {

      // If not a leaf-cell, then open cell to first child cell
      if (tree[cc].level != ltot)
         cc++;

      // If leaf-cell, add particles to list
      else if (tree[cc].level == ltot && Ndirect + Nleafmax <= Ndirectmax) {
        i = tree[cc].ifirst;
        while (i != -1) {
          directlist[Ndirect++] = i;
          if (i == tree[cc].ilast) break;
       	  i = inext[i];
        };
        cc = tree[cc].cnext;
      }

      // If leaf-cell, but we've run out of memory, return with error-code (-1)
      else if (tree[cc].level == ltot && Ndirect + Nleafmax > Ndirectmax)
       	return -1;

    }

    // If not in range, then open next cell
    //-------------------------------------------------------------------------
    else
      cc = tree[cc].cnext;

  };
  //===========================================================================

  // Now, trim the list to remove particles that are definitely not neighbours.
  // If not an SPH neighbour, then add to direct gravity sum list.
  hrangemax = hrangemax*hrangemax;
  for (j=Nneibtemp; j<Nneib; j++) {
    i = neiblist[j];
    for (k=0; k<ndim; k++) dr[k] = sphdata[i].r[k] - rc[k];
    drsqd = DotProduct(dr,dr,ndim);
    if (drsqd < hrangemax || drsqd < 
        (rmax + kernrange*sphdata[i].h)*(rmax + kernrange*sphdata[i].h))
      neiblist[Nneibtemp++] = i;
    else if (Ndirect < Ndirectmax)
      directlist[Ndirect++] = i;
    else
     return -1;
  }
  Nneib = Nneibtemp;

  return 1;
}



//=============================================================================
//  BinaryTree2::ComputeCellMonopoleForces
/// Compute the force on particle 'parti' due to all cells obtained in the 
/// gravity tree walk.  Uses only monopole moments (i.e. COM) of the cell.
//=============================================================================
template <int ndim>
void BinaryTree2<ndim>::ComputeCellMonopoleForces
(int i,                               ///< [in] i.d. of particle
 int Ngravcell,                       ///< [in] No. of tree cells in list
 BinaryTree2Cell<ndim> **gravcelllist, ///< [in] List of tree cell ids
 SphParticle<ndim> &parti)            ///< [inout] SPH particle
{
  int cc;                           // Aux. cell counter
  int k;                            // Dimension counter
  FLOAT dr[ndim];                   // Relative position vector
  FLOAT drsqd;                      // Distance squared
  FLOAT invdrmag;                   // 1 / distance
  FLOAT invdrsqd;                   // 1 / drsqd
  FLOAT invdr3;                     // 1 / dist^3
  FLOAT mc;                         // Mass of cell
  FLOAT rp[ndim];                   // Position of particle
  BinaryTree2Cell<ndim> *cell;       // Pointer to gravity tree cell

  for (k=0; k<ndim; k++) rp[k] = parti.r[k];

  // Loop over all neighbouring particles in list
  //---------------------------------------------------------------------------
  for (cc=0; cc<Ngravcell; cc++) {
    cell = gravcelllist[cc];

    mc = cell->m;
    for (k=0; k<ndim; k++) dr[k] = cell->r[k] - rp[k];
    drsqd = DotProduct(dr,dr,ndim) + small_number;
    invdrsqd = 1.0/drsqd;
    invdrmag = sqrt(invdrsqd);
    invdr3 = invdrsqd*invdrmag;

    parti.gpot += mc*invdrmag;
    for (k=0; k<ndim; k++) parti.agrav[k] += mc*dr[k]*invdr3;

  }
  //---------------------------------------------------------------------------

  return;
}



//=============================================================================
//  BinaryTree2::ComputeCellQuadrupoleForces
/// Compute the force on particle 'parti' due to all cells obtained in the 
/// gravity tree walk including the quadrupole moment correction term.
//=============================================================================
template <int ndim>
void BinaryTree2<ndim>::ComputeCellQuadrupoleForces
(int i,                               ///< [in] i.d. of particle
 int Ngravcell,                       ///< [in] No. of tree cells in list
 BinaryTree2Cell<ndim> **gravcelllist, ///< [in] List of tree cell ids
 SphParticle<ndim> &parti)            ///< [inout] SPH particle
{
  int cc;                           // Aux. cell counter
  int k;                            // Dimension counter
  FLOAT dr[ndim];                   // Relative position vector
  FLOAT drsqd;                      // Distance squared
  FLOAT invdrsqd;                   // 1 / drsqd
  FLOAT invdrmag;                   // 1 / distance
  FLOAT invdr5;                     // 1 / distance^5
  FLOAT qscalar;                    // Quadrupole moment scalar quantity
  BinaryTree2Cell<ndim> *cell;       // Pointer to gravity tree cell

  // Loop over all neighbouring particles in list
  //---------------------------------------------------------------------------
  for (cc=0; cc<Ngravcell; cc++) {
    cell = gravcelllist[cc];

    for (k=0; k<ndim; k++) dr[k] = cell->r[k] - parti.r[k];
    drsqd = DotProduct(dr,dr,ndim);
    invdrsqd = 1.0/(drsqd + small_number);
    invdrmag = sqrt(invdrsqd);
    invdr5 = invdrsqd*invdrsqd*invdrmag;

    // First add monopole terms
    parti.gpot += cell->m*invdrmag;
    for (k=0; k<ndim; k++) 
      parti.agrav[k] += cell->m*dr[k]*invdrsqd*invdrmag;

    // Now add quadrupole moment terms depending on dimensionality
    if (ndim == 3) {
      qscalar = cell->q[0]*dr[0]*dr[0] + cell->q[2]*dr[1]*dr[1] -
        (cell->q[0] + cell->q[2])*dr[2]*dr[2] +
         2.0*(cell->q[1]*dr[0]*dr[1] + cell->q[3]*dr[0]*dr[2] +
         cell->q[4]*dr[1]*dr[2]);
      parti.agrav[0] += 
        (cell->q[0]*dr[0] + cell->q[1]*dr[1] + cell->q[3]*dr[2])*invdr5
         - 2.5*qscalar*dr[0]*invdr5*invdrsqd;
      parti.agrav[1] += 
        (cell->q[1]*dr[0] + cell->q[2]*dr[1] + cell->q[4]*dr[2])*invdr5
         - 2.5*qscalar*dr[1]*invdr5*invdrsqd;
      parti.agrav[2] += 
        (cell->q[3]*dr[0] + cell->q[4]*dr[1] -
         (cell->q[0] + cell->q[2])*dr[2])*invdr5
          - 2.5*qscalar*dr[1]*invdr5*invdrsqd;
    }
    else if (ndim == 2) {
      qscalar = cell->q[0]*dr[0]*dr[0] + cell->q[2]*dr[1]*dr[1] +
        2.0*cell->q[1]*dr[0]*dr[1];
      parti.agrav[0] += (cell->q[0]*dr[0] + cell->q[1]*dr[1])*invdr5 -
        2.5*qscalar*dr[0]*invdr5*invdrsqd;
      parti.agrav[1] += (cell->q[1]*dr[0] + cell->q[2]*dr[1])*invdr5 -
        2.5*qscalar*dr[1]*invdr5*invdrsqd;
    }
    parti.gpot += 0.5*qscalar*invdr5;

  }
  //---------------------------------------------------------------------------

  return;
}



//=============================================================================
//  BinarySubTree::ComputeActiveCellList
/// Returns the number of cells containing active particles, 'Nactive', and
/// the i.d. list of cells contains active particles, 'celllist'
//=============================================================================
template <int ndim>
int BinaryTree2<ndim>::ComputeActiveCellList
(BinaryTree2Cell<ndim> **celllist) ///< Cells id array containing active ptcls
{
  int c;                           // Cell counter
  int Nactive = 0;                 // ..

  for (c=0; c<Ncell; c++)
    if (tree[c].Nactive > 0) celllist[Nactive++] = &tree[c];

  return Nactive;
}



//=============================================================================
//  BinaryTree2::UpdateAllSphProperties
/// Compute all local 'gather' properties of currently active particles, and 
/// then compute each particle's contribution to its (active) neighbour 
/// neighbour hydro forces.  Optimises the algorithm by using grid-cells to 
/// construct local neighbour lists for all particles  inside the cell.
//=============================================================================
template <int ndim>
void BinaryTree2<ndim>::UpdateAllSphProperties
(Sph<ndim> *sph,                   ///< [inout] Pointer to main SPH object
 Nbody<ndim> *nbody)               ///> [in] ..
{
  int celldone;                    // ..
  int okflag;                      // ..
  int cc;                          // Aux. cell counter
  int cactive;                     // No. of active
  int i;                           // Particle id
  int j;                           // Aux. particle counter
  int jj;                          // Aux. particle counter
  int k;                           // Dimension counter
  int Nactive;                     // No. of active particles in cell
  int Ngather;                     // No. of near gather neighbours
  int Nneib;                       // No. of neighbours
  int Nneibmax;                    // Max. no. of neighbours
  int *activelist;                 // List of active particle ids
  int *neiblist;                   // List of neighbour ids
  FLOAT draux[ndim];               // Aux. relative position vector var
  FLOAT drsqdaux;                  // Distance squared
  FLOAT hrangesqd;                 // Kernel extent
  FLOAT hmax;                      // Maximum smoothing length
  FLOAT rp[ndim];                  // Local copy of particle position
  FLOAT *drsqd;                    // Position vectors to gather neibs
  FLOAT *gpot;                     // Potential for particles
  FLOAT *gpot2;                    // ..
  FLOAT *m;                        // Distances to gather neibs
  FLOAT *m2;                       // ..
  FLOAT *mu;                       // mass*u for gather neibs
  FLOAT *mu2;                      // ..
  FLOAT *r;                        // Positions of neibs
  BinaryTree2Cell<ndim> *cell;      // Pointer to binary tree cell
  BinaryTree2Cell<ndim> **celllist; // List of binary cell pointers
  SphParticle<ndim> *data = sph->sphdata;  // Pointer to SPH particle data

  debug2("[BinaryTree2::UpdateAllSphProperties]");

  // Find list of all cells that contain active particles
  celllist = new BinaryTree2Cell<ndim>*[gtot];
  cactive = ComputeActiveCellList(celllist);


  // Set-up all OMP threads
  //===========================================================================
#pragma omp parallel default(none) private(activelist,cc,cell,celldone,draux)\
  private(drsqd,drsqdaux,hmax,hrangesqd,i,j,jj,k,okflag,m,mu,Nactive,neiblist)\
  private(Nneib,Nneibmax,r,rp,gpot,gpot2,m2,mu2,Ngather)\
  shared(sph,celllist,cactive,data,nbody)
  {
    Nneibmax = 2*sph->Ngather;
    activelist = new int[Nleafmax];
    neiblist = new int[Nneibmax];
    gpot = new FLOAT[Nneibmax];
    gpot2 = new FLOAT[Nneibmax];
    drsqd = new FLOAT[Nneibmax];
    m = new FLOAT[Nneibmax];
    m2 = new FLOAT[Nneibmax];
    mu = new FLOAT[Nneibmax];
    mu2 = new FLOAT[Nneibmax];
    r = new FLOAT[Nneibmax*ndim];

    // Loop over all active cells
    //=========================================================================
#pragma omp for schedule(dynamic)
    for (cc=0; cc<cactive; cc++) {
      cell = celllist[cc];
      celldone = 1;
      hmax = cell->hmax;

      // If hmax is too small so the neighbour lists are invalid, make hmax
      // larger and then recompute for the current active cell.
      //-----------------------------------------------------------------------
      do {
        hmax = 1.05*hmax;
        celldone = 1;

        // Find list of active particles in current cell
        Nactive = ComputeActiveParticleList(cell,sph,activelist);

        // Compute neighbour list for cell depending on physics options
        Nneib = ComputeGatherNeighbourList(cell,Nneibmax,neiblist,
                                           hmax,sph->sphdata);

        // If there are too many neighbours, reallocate the arrays and
        // recompute the neighbour lists.
        while (Nneib == -1) {
          delete[] r;
          delete[] mu2;
          delete[] mu;
          delete[] m2;
          delete[] m;
          delete[] drsqd;
          delete[] gpot2;
          delete[] gpot;
          delete[] neiblist;
          Nneibmax = 2*Nneibmax;
          neiblist = new int[Nneibmax];
          gpot = new FLOAT[Nneibmax];
          gpot2 = new FLOAT[Nneibmax];
          drsqd = new FLOAT[Nneibmax];
          m = new FLOAT[Nneibmax];
          m2 = new FLOAT[Nneibmax];
          mu = new FLOAT[Nneibmax];
          mu2 = new FLOAT[Nneibmax];
          r = new FLOAT[Nneibmax*ndim];
          Nneib = ComputeGatherNeighbourList(cell,Nneibmax,neiblist,
                                             hmax,sph->sphdata);
        };

        // Make local copies of important neib information (mass and position)
        for (jj=0; jj<Nneib; jj++) {
          j = neiblist[jj];
          assert(j >= 0 && j < sph->Ntot);
          gpot[jj] = data[j].gpot;
          m[jj] = data[j].m;
          mu[jj] = data[j].m*data[j].u;
          for (k=0; k<ndim; k++) r[ndim*jj + k] = (FLOAT) data[j].r[k];
        }


        // Loop over all active particles in the cell
        //---------------------------------------------------------------------
        for (j=0; j<Nactive; j++) {
          i = activelist[j];
          assert(i >= 0 && i < sph->Nsph);
          for (k=0; k<ndim; k++) rp[k] = data[i].r[k];

          // Set gather range as current h multiplied by some tolerance factor
          hrangesqd = sph->kernp->kernrangesqd*hmax*hmax;
          Ngather = 0;

          // Compute distance (squared) to all
          //-------------------------------------------------------------------
          for (jj=0; jj<Nneib; jj++) {
            for (k=0; k<ndim; k++) draux[k] = r[ndim*jj + k] - rp[k];
            drsqdaux = DotProduct(draux,draux,ndim);

            // Record distance squared for all potential gather neighbours
            if (drsqdaux <= hrangesqd) {
              gpot[Ngather] = gpot[jj];
              drsqd[Ngather] = drsqdaux;
              m2[Ngather] = m[jj];
              mu2[Ngather] = mu[jj];
              Ngather++;
            }

          }
          //-------------------------------------------------------------------

          // Validate that gather neighbour list is correct
#if defined(VERIFY_ALL)
          if (neibcheck) CheckValidNeighbourList(sph,i,Nneib,neiblist,"gather");
#endif

          // Compute smoothing length and other gather properties for ptcl i
          okflag = sph->ComputeH(i,Ngather,hmax,m2,mu2,
                                 drsqd,gpot,data[i],nbody);

          // If h-computation is invalid, then break from loop and recompute
          // larger neighbour lists
          if (okflag == 0) {
            celldone = 0;
            break;
          }

        }
        //---------------------------------------------------------------------

      } while (celldone == 0);
      //-----------------------------------------------------------------------

    }
    //=========================================================================

    // Free-up all memory
    delete[] r;
    delete[] mu2;
    delete[] mu;
    delete[] m2;
    delete[] m;
    delete[] drsqd;
    delete[] gpot2;
    delete[] gpot;
    delete[] neiblist;
    delete[] activelist;

  }
  //===========================================================================

  delete[] celllist;

  // Update all tree smoothing length values
  //UpdateHmaxValues(sph->sphdata);

  return;
}



//=============================================================================
//  BinaryTree2::UpdateAllSphHydroForces
/// Compute all local 'gather' properties of currently active particles, and 
/// then compute each particle's contribution to its (active) neighbour 
/// neighbour hydro forces.  Optimises the algorithm by using grid-cells to 
/// construct local neighbour lists for all particles  inside the cell.
//=============================================================================
template <int ndim>
void BinaryTree2<ndim>::UpdateAllSphHydroForces
(Sph<ndim> *sph)                   ///< Pointer to SPH object
{
  int cactive;                     // No. of active cells
  int cc;                          // Aux. cell counter
  int i;                           // Particle id
  int j;                           // Aux. particle counter
  int jj;                          // Aux. particle counter
  int k;                           // Dimension counter
  int Nactive;                     // No. of active particles in cell
  int Ninteract;                   // No. of near gather neighbours
  int Nneib;                       // No. of neighbours
  int Nneibmax;                    // Max. no. of neighbours
  int *activelist;                 // List of active particle ids
  int *interactlist;               // List of interactng SPH neighbours
  int *neiblist;                   // List of neighbour ids
  FLOAT draux[ndim];               // Aux. relative position vector
  FLOAT drsqd;                     // Distance squared
  FLOAT hrangesqdi;                // Kernel gather extent
  FLOAT rp[ndim];                  // Local copy of particle position
  FLOAT *dr;                       // Array of relative position vectors
  FLOAT *drmag;                    // Array of neighbour distances
  FLOAT *invdrmag;                 // Array of 1/drmag between particles
  BinaryTree2Cell<ndim> *cell;      // Pointer to binary tree cell
  BinaryTree2Cell<ndim> **celllist; // List of binary tree pointers
  SphParticle<ndim> *activepart;   // Local copy of active particles
  SphParticle<ndim> *neibpart;     // Local copy of neighbouring ptcls
  SphParticle<ndim> *data = sph->sphdata;   // Pointer to SPH particle data

  debug2("[BinaryTree2::UpdateAllSphHydroForces]");


  // Update tree smoothing length values here
  UpdateHmaxValues(sph->sphdata);

  // Find list of all cells that contain active particles
  celllist = new BinaryTree2Cell<ndim>*[gtot];
  cactive = ComputeActiveCellList(celllist);


  // Set-up all OMP threads
  //===========================================================================
#pragma omp parallel default(none) private(activelist,activepart,cc,cell,dr)\
  private(draux,drmag,drsqd,hrangesqdi,i,interactlist,invdrmag,j,jj,k) \
  private(Nactive,neiblist,neibpart,Ninteract,Nneib,Nneibmax,rp)\
  shared(cactive,celllist,data,sph)
  {
    Nneibmax = 4*sph->Ngather;
    activelist = new int[Nleafmax];
    activepart = new SphParticle<ndim>[Nleafmax];
    neiblist = new int[Nneibmax];
    interactlist = new int[Nneibmax];
    dr = new FLOAT[Nneibmax*ndim];
    drmag = new FLOAT[Nneibmax];
    invdrmag = new FLOAT[Nneibmax];
    neibpart = new SphParticle<ndim>[Nneibmax];

    // Loop over all active cells
    //=========================================================================
#pragma omp for schedule(dynamic)
    for (cc=0; cc<cactive; cc++) {
      cell = celllist[cc];

      // Find list of active particles in current cell
      Nactive = ComputeActiveParticleList(cell,sph,activelist);

      // Make local copies of active particles
      for (j=0; j<Nactive; j++) {
        assert(activelist[j] >= 0 && activelist[j] < sph->Nsph);
        activepart[j] = data[activelist[j]];
        activepart[j].div_v = (FLOAT) 0.0;
        activepart[j].dudt = (FLOAT) 0.0;
        activepart[j].levelneib = 0;
        for (k=0; k<ndim; k++) activepart[j].a[k] = (FLOAT) 0.0;
      }

      // Compute neighbour list for cell depending on physics options
      Nneib = ComputeNeighbourList(cell,Nneibmax,neiblist,sph->sphdata);

      // If there are too many neighbours, reallocate the arrays and
      // recompute the neighbour list.
      while (Nneib == -1) {
        delete[] neibpart;
        delete[] invdrmag;
        delete[] drmag;
        delete[] dr;
        delete[] interactlist;
        delete[] neiblist;
        Nneibmax = 2*Nneibmax;
        neiblist = new int[Nneibmax];
        interactlist = new int[Nneibmax];
        dr = new FLOAT[Nneibmax*ndim];
        drmag = new FLOAT[Nneibmax];
        invdrmag = new FLOAT[Nneibmax];
        neibpart = new SphParticle<ndim>[Nneibmax];
        Nneib = ComputeNeighbourList(cell,Nneibmax,neiblist,sph->sphdata);
      };

      // Make local copies of all potential neighbours
      for (j=0; j<Nneib; j++) {
        assert(neiblist[j] >= 0 && neiblist[j] < sph->Ntot);
        neibpart[j] = data[neiblist[j]];
        neibpart[j].div_v = (FLOAT) 0.0;
        neibpart[j].dudt = (FLOAT) 0.0;
        neibpart[j].levelneib=0;
        for (k=0; k<ndim; k++) neibpart[j].a[k] = (FLOAT) 0.0;
      }

      //cout << "Nactive2 : " << Nactive << "     " << Nneib << endl;

      // Loop over all active particles in the cell
      //-----------------------------------------------------------------------
      for (j=0; j<Nactive; j++) {
        i = activelist[j];

        for (k=0; k<ndim; k++) rp[k] = activepart[j].r[k]; //data[i].r[k];
        hrangesqdi = activepart[j].hrangesqd;
        Ninteract = 0;

        // Validate that gather neighbour list is correct
#if defined(VERIFY_ALL)
        if (neibcheck) CheckValidNeighbourList(sph,i,Nneib,neiblist,"all");
#endif

        // Compute distances and the inverse between the current particle
        // and all neighbours here, for both gather and inactive scatter neibs.
        // Only consider particles with j > i to compute pair forces once
        // unless particle j is inactive.
        //---------------------------------------------------------------------
        for (jj=0; jj<Nneib; jj++) {

          // Skip neighbour if it's not the correct part of an active pair
          if (neiblist[jj] <= i && neibpart[jj].active) continue;

          for (k=0; k<ndim; k++) draux[k] = neibpart[jj].r[k] - rp[k];
          drsqd = DotProduct(draux,draux,ndim) + small_number;

          // Compute relative position and distance quantities for pair
          if (drsqd <= hrangesqdi || drsqd <= neibpart[jj].hrangesqd) {
            drmag[Ninteract] = sqrt(drsqd);
            invdrmag[Ninteract] = (FLOAT) 1.0/drmag[Ninteract];
            for (k=0; k<ndim; k++)
              dr[Ninteract*ndim + k] = draux[k]*invdrmag[Ninteract];
            interactlist[Ninteract] = jj;
            Ninteract++;
          }

        }
        //---------------------------------------------------------------------

        // Compute all gather neighbour contributions to hydro forces
        sph->ComputeSphHydroForces(i,Ninteract,interactlist,
				   drmag,invdrmag,dr,activepart[j],neibpart);

      }
      //-----------------------------------------------------------------------


      // Add all active particles contributions to main array
      for (j=0; j<Nactive; j++) {
    	i = activelist[j];
#if defined _OPENMP
    	omp_lock_t& lock = sph->GetParticleILock(i);
    	omp_set_lock(&lock);
#endif
        for (k=0; k<ndim; k++) {
          data[i].a[k] += activepart[j].a[k];
        }
        data[i].dudt += activepart[j].dudt;
        data[i].div_v += activepart[j].div_v;
        data[i].levelneib = max(data[i].levelneib,activepart[j].levelneib);
#if defined _OPENMP
        omp_unset_lock(&lock);
#endif
      }

      // Now add all active neighbour contributions to main array
      for (jj=0; jj<Nneib; jj++) {
        j = neiblist[jj];
#if defined _OPENMP
        omp_lock_t& lock = sph->GetParticleILock(j);
        omp_set_lock(&lock);
#endif
        if (neibpart[jj].active) {
          for (k=0; k<ndim; k++) {
            data[j].a[k] += neibpart[jj].a[k];
          }
          data[j].dudt += neibpart[jj].dudt;
          data[j].div_v += neibpart[jj].div_v;
        }
        data[j].levelneib = max(data[j].levelneib,neibpart[jj].levelneib);
#if defined _OPENMP
        omp_unset_lock(&lock);
#endif
      }

    }
    //=========================================================================

    // Free-up local memory for OpenMP thread
    delete[] neibpart;
    delete[] invdrmag;
    delete[] drmag;
    delete[] dr;
    delete[] interactlist;
    delete[] neiblist;
    delete[] activepart;
    delete[] activelist;

  }
  //===========================================================================

  delete[] celllist;


  // Compute other important SPH quantities after hydro forces are computed
  if (sph->hydro_forces == 1) {
    for (i=0; i<sph->Nsph; i++) {
      if (sph->sphdata[i].active)
    	  sph->ComputePostHydroQuantities(sph->sphdata[i]);
    }
  }

  return;
}



//=============================================================================
//  BinaryTree2::UpdateAllSphForces
/// Compute all local 'gather' properties of currently active particles, and 
/// then compute each particle's contribution to its (active) neighbour 
/// neighbour hydro forces.  Optimises the algorithm by using grid-cells to 
/// construct local neighbour lists for all particles  inside the cell.
//=============================================================================
template <int ndim>
void BinaryTree2<ndim>::UpdateAllSphForces
(Sph<ndim> *sph)                    ///< Pointer to SPH object
{
  int cactive;                      // No. of active cells
  int cc;                           // Aux. cell counter
  int i;                            // Particle id
  int j;                            // Aux. particle counter
  int jj;                           // Aux. particle counter
  int k;                            // Dimension counter
  int okflag;                       // Flag if h-rho iteration is valid
  int Nactive;                      // No. of active particles in cell
  int Ndirect;                      // No. of direct-sum gravity particles
  int Ndirectmax;                   // Max. no. of direct sum particles
  int Ngravcell;                    // No. of gravity cells
  int Ngravcellmax;                 // Max. no. of gravity cells
  int Ninteract;                    // No. of interactions with hydro neibs
  int Nneib;                        // No. of neighbours
  int Nneibmax;                     // Max. no. of neighbours
  int *activelist;                  // List of active particle ids
  int *directlist;                  // List of direct sum particle ids
  int *interactlist;                // List of interacting neighbour ids
  int *neiblist;                    // List of neighbour ids
  FLOAT *agrav;                     // Local copy of gravitational accel.
  FLOAT *gpot;                      // Local copy of gravitational pot.
  BinaryTree2Cell<ndim> *cell;      // Pointer to binary tree cell
  BinaryTree2Cell<ndim> **celllist; // List of pointers to binary tree cells
  BinaryTree2Cell<ndim> **gravcelllist; // List of pointers to grav. cells
  SphParticle<ndim> *neibpart;      // Local copy of neighbouring ptcls
  SphParticle<ndim> *activepart;    // Local copy of SPH particle


  debug2("[BinaryTree2::UpdateAllSphForces]");


  // Update tree smoothing length values here
  UpdateHmaxValues(sph->sphdata);


  // Find list of all cells that contain active particles
  celllist = new BinaryTree2Cell<ndim>*[gtot];
  cactive = ComputeActiveCellList(celllist);


  // Set-up all OMP threads
  //===========================================================================
#pragma omp parallel default(none) private(activelist,agrav,cc,cell)\
  private(gpot,i,interactlist,j,jj,activepart)\
  private(k,okflag,Nactive,neiblist,neibpart,Ninteract,Nneib,directlist)\
  private(gravcelllist,Ngravcell,Ndirect,Nneibmax,Ndirectmax,Ngravcellmax)\
  shared(celllist,cactive,sph,cout)
  {
    Nneibmax = 4*sph->Ngather;
    Ndirectmax = 2*Nneibmax;
    Ngravcellmax = 2*Nneibmax;
    agrav = new FLOAT[ndim*sph->Nsph];
    gpot = new FLOAT[ndim*sph->Nsph];
    activelist = new int[Nleafmax];
    activepart = new SphParticle<ndim>[Nleafmax];
    neiblist = new int[Nneibmax];
    interactlist = new int[Nneibmax];
    directlist = new int[Ndirectmax];
    gravcelllist = new BinaryTree2Cell<ndim>*[Ngravcellmax];
    neibpart = new SphParticle<ndim>[Nneibmax];

    // Zero temporary grav. accel array
    for (i=0; i<ndim*sph->Nsph; i++) agrav[i] = 0.0;
    for (i=0; i<sph->Nsph; i++) gpot[i] = 0.0;


    // Loop over all active cells
    //=========================================================================
#pragma omp for schedule(dynamic)
    for (cc=0; cc<cactive; cc++) {
      cell = celllist[cc];

      // Find list of active particles in current cell
      Nactive = ComputeActiveParticleList(cell,sph,activelist);

      // Make local copies of active particles
      for (j=0; j<Nactive; j++) {
        assert(activelist[j] >= 0 && activelist[j] < sph->Nsph);
        activepart[j] = sph->sphdata[activelist[j]];
        activepart[j].div_v = (FLOAT) 0.0;
        activepart[j].dudt = (FLOAT) 0.0;
        activepart[j].gpot = activepart[j].m*activepart[j].invh*sph->kernp->wpot(0.0);
        for (k=0; k<ndim; k++) activepart[j].a[k] = (FLOAT) 0.0;
        for (k=0; k<ndim; k++) activepart[j].agrav[k] = (FLOAT) 0.0;
      }

      // Compute neighbour list for cell depending on physics options
      okflag = ComputeGravityInteractionList(cell,Nneibmax,Ndirectmax,
                                             Ngravcellmax,Nneib,Ndirect,
                                             Ngravcell,neiblist,directlist,
                                             gravcelllist,sph->sphdata);

      // If there are too many neighbours, reallocate the arrays and
      // recompute the neighbour lists.
      while (okflag == -1) {
        delete[] neibpart;
        delete[] gravcelllist;
        delete[] directlist;
        delete[] interactlist;
        delete[] neiblist;
        Nneibmax = 2*Nneibmax;
        Ndirectmax = 2*Ndirectmax;
        Ngravcellmax = 2*Ngravcellmax;
        neiblist = new int[Nneibmax];
        interactlist = new int[Nneibmax];
        directlist = new int[Ndirectmax];
        gravcelllist = new BinaryTree2Cell<ndim>*[Ngravcellmax];
        neibpart = new SphParticle<ndim>[Nneibmax];
        okflag = ComputeGravityInteractionList(cell,Nneibmax,Ndirectmax,
                                               Ngravcellmax,Nneib,Ndirect,
                                               Ngravcell,neiblist,directlist,
                                               gravcelllist,sph->sphdata);
      };


      // Make local copies of all potential neighbours
      for (j=0; j<Nneib; j++) {
        assert(neiblist[j] >= 0 && neiblist[j] < sph->Ntot);
        neibpart[j] = sph->sphdata[neiblist[j]];
        neibpart[j].div_v = (FLOAT) 0.0;
        neibpart[j].dudt = (FLOAT) 0.0;
        neibpart[j].gpot = (FLOAT) 0.0;
        for (k=0; k<ndim; k++) neibpart[j].a[k] = (FLOAT) 0.0;
        for (k=0; k<ndim; k++) neibpart[j].agrav[k] = (FLOAT) 0.0;
      }

      // Loop over all active particles in the cell
      //-----------------------------------------------------------------------
      for (j=0; j<Nactive; j++) {
        i = activelist[j];

        // Determine SPH neighbour interaction list 
        // (to ensure we don't compute pair-wise forces twice)
        Ninteract = 0;
        for (jj=0; jj<Nneib; jj++) {
          if ((neiblist[jj] < i && !neibpart[jj].active) || neiblist[jj] > i)
            interactlist[Ninteract++] = jj;
        }

        // Compute forces between SPH neighbours (hydro and gravity)
        sph->ComputeSphHydroGravForces(i,Ninteract,interactlist,
                                       activepart[j],neibpart);

        // Compute direct gravity forces between distant particles
        sph->ComputeDirectGravForces(i,Ndirect,directlist,
                                     agrav,gpot,activepart[j],sph->sphdata);

        // Compute gravitational force due to distant cells
        if (multipole == "monopole")
          ComputeCellMonopoleForces(i,Ngravcell,gravcelllist,activepart[j]);
        else if (multipole == "quadrupole")
          ComputeCellQuadrupoleForces(i,Ngravcell,gravcelllist,activepart[j]);

      }
      //-----------------------------------------------------------------------


      // Add all active particles contributions to main array
      for (j=0; j<Nactive; j++) {
    	i = activelist[j];
#if defined _OPENMP
        omp_lock_t& lock = sph->GetParticleILock(i);
        omp_set_lock(&lock);
#endif
        for (k=0; k<ndim; k++) {
          sph->sphdata[i].a[k] += activepart[j].a[k];
          sph->sphdata[i].agrav[k] += activepart[j].agrav[k];
        }
        sph->sphdata[i].gpot += activepart[j].gpot;
        sph->sphdata[i].dudt += activepart[j].dudt;
        sph->sphdata[i].div_v += activepart[j].div_v;
        sph->sphdata[i].levelneib = max(sph->sphdata[i].levelneib,activepart[j].levelneib);
#if defined _OPENMP
        omp_unset_lock(&lock);
#endif
      }

      // Now add all active neighbour contributions to the main arrays
      for (jj=0; jj<Nneib; jj++) {
        j = neiblist[jj];
#if defined _OPENMP
        omp_lock_t& lock = sph->GetParticleILock(j);
        omp_set_lock(&lock);
#endif
        if (neibpart[jj].active) {
          for (k=0; k<ndim; k++) {
            sph->sphdata[j].a[k] += neibpart[jj].a[k];
            sph->sphdata[j].agrav[k] += neibpart[jj].agrav[k];
          }
          sph->sphdata[j].gpot += neibpart[jj].gpot;
          sph->sphdata[j].dudt += neibpart[jj].dudt;
          sph->sphdata[j].div_v += neibpart[jj].div_v;
        }
        sph->sphdata[j].levelneib = max(sph->sphdata[j].levelneib,neibpart[jj].levelneib);
#if defined _OPENMP
        omp_unset_lock(&lock);
#endif
      }

    }
    //=========================================================================


    // Finally, add all contributions from distant pair-wise forces to arrays
    for (i=0; i<sph->Nsph; i++) {
      if (sph->sphdata[i].active) {
        for (k=0; k<ndim; k++) {
#pragma omp atomic
          sph->sphdata[i].agrav[k] += agrav[ndim*i + k];
        }
#pragma omp atomic
        sph->sphdata[i].gpot += gpot[i];
      }
    }


    // Free-up local memory for OpenMP thread
    delete[] neibpart;
    delete[] gravcelllist;
    delete[] directlist;
    delete[] interactlist;
    delete[] neiblist;
    delete[] activepart;
    delete[] activelist;
    delete[] gpot;
    delete[] agrav;

  }
  //===========================================================================

  delete[] celllist;


  // Compute other important SPH quantities after hydro forces are computed
  if (sph->hydro_forces == 1) {
    for (i=0; i<sph->Nsph; i++) {
      if (sph->sphdata[i].active)
        sph->ComputePostHydroQuantities(sph->sphdata[i]);
    }
  }

  return;
}



//=============================================================================
//  BinaryTree2::UpdateAllSphGravityProperties
/// Compute all local 'gather' properties of currently active particles, and 
/// then compute each particle's contribution to its (active) neighbour 
/// neighbour hydro forces.  Optimises the algorithm by using grid-cells to 
/// construct local neighbour lists for all particles  inside the cell.
//=============================================================================
template <int ndim>
void BinaryTree2<ndim>::UpdateAllSphGravForces
(Sph<ndim> *sph)                    ///< Pointer to SPH object
{
  return;
}



//=============================================================================
//  BinaryTree2::UpdateAllSphDerivatives
/// ..
//=============================================================================
template <int ndim>
void BinaryTree2<ndim>::UpdateAllSphDerivatives(Sph<ndim> *sph)
{
  return;
}



//=============================================================================
//  BinaryTree2::UpdateAllSphDudt
/// Compute all local 'gather' properties of currently active particles, and 
/// then compute each particle's contribution to its (active) neighbour 
/// neighbour hydro forces.  Optimises the algorithm by using grid-cells to 
/// construct local neighbour lists for all particles  inside the cell.
//=============================================================================
template <int ndim>
void BinaryTree2<ndim>::UpdateAllSphDudt(Sph<ndim> *sph)
{
  return;
}



#if defined(VERIFY_ALL)
//=============================================================================
//  BinaryTree2::ValidateTree
/// ..
//=============================================================================
template <int ndim>
void BinaryTree2<ndim>::ValidateTree
(Sph<ndim> *sph)
{
  bool overlap_flag = false;
  int c;
  int cc;
  int i;
  int j;
  int l;
  int *pcount;
  BinaryTree2Cell<ndim> cell;

  debug2("[BinaryTree2::ValidateTree]");

  pcount = new int[Ntot];
  for (i=0; i<Ntot; i++) pcount[i] = 0;


  //---------------------------------------------------------------------------
  for (c=0; c<Ncell; c++) {
    cell = tree[c];

    // Check that particles are not in linked lists more than once
    if (cell.level == ltot) {
      i = cell.ifirst;
      while (i != -1) {
	pcount[i]++;
        if (i == cell.ilast) break;
	i = inext[i];
      }
    }

    // Check that bounding boxes of cells on each level do not overlap 
    // each other
    for (cc=0; cc<Ncell; cc++) {
      if (c != cc && tree[cc].level == cell.level) {
	if (ndim == 2) {
	  if (cell.bbmin[0] < tree[cc].bbmax[0] &&
	      cell.bbmax[0] > tree[cc].bbmin[0] &&
	      cell.bbmin[1] < tree[cc].bbmax[1] &&
	      cell.bbmax[1] > tree[cc].bbmin[1])
	    overlap_flag = true;
	}
      }
      if (overlap_flag) {
	cout << "Brother/sister cell overlap error!! : " << c << "   " 
             << cc << endl;
	exit(0);
      }
    }
  }
  //---------------------------------------------------------------------------


  for (i=0; i<Ntot; i++) {
    if (pcount[i] != 1) {
      cout << "Problem with child cell ptcl counter : " << i << "   " 
	   << pcount[i] << endl;
      exit(0);
    }
  }


  delete[] pcount;

  return;
}



//=============================================================================
//  BinaryTree2::CheckValidNeighbourList
/// Checks that the neighbour list generated by the grid is valid in that it 
/// (i) does include all true neighbours, and 
/// (ii) all true neigbours are only included once and once only.
//=============================================================================
template <int ndim>
void BinaryTree2<ndim>::CheckValidNeighbourList
(Sph<ndim> *sph,                    ///< [in] SPH object pointer
 int i,                             ///< [in] Particle i.d.
 int Nneib,                         ///< [in] No. of potential neighbours
 int *neiblist,                     ///< [in] List of potential neighbour i.d.s
 string neibtype)                   ///< [in] Neighbour search type
{
  int count;                        // Valid neighbour counter
  int j;                            // Neighbour particle counter
  int k;                            // Dimension counter
  int Ntrueneib = 0;                // No. of 'true' neighbours
  int *trueneiblist;                // List of true neighbour ids
  FLOAT drsqd;                      // Distance squared
  FLOAT dr[ndim];                   // Relative position vector

  // Allocate array to store local copy of potential neighbour ids
  trueneiblist = new int[sph->Ntot];

  // First, create list of 'true' neighbours by looping over all particles
  if (neibtype == "gather") {
    for (j=0; j<sph->Ntot; j++) {
      for (k=0; k<ndim; k++)
	dr[k] = sph->sphdata[j].r[k] - sph->sphdata[i].r[k];
      drsqd = DotProduct(dr,dr,ndim);
      if (drsqd <
	  sph->kernp->kernrangesqd*sph->sphdata[i].h*sph->sphdata[i].h)
	trueneiblist[Ntrueneib++] = j;
    }
  }
  else if (neibtype == "all") {
    for (j=0; j<sph->Ntot; j++) {
      for (k=0; k<ndim; k++) 
        dr[k] = sph->sphdata[j].r[k] - sph->sphdata[i].r[k];
      drsqd = DotProduct(dr,dr,ndim);
      if (drsqd < sph->kernp->kernrangesqd*sph->sphdata[i].h*sph->sphdata[i].h ||
    	  drsqd < sph->kernp->kernrangesqd*sph->sphdata[j].h*sph->sphdata[j].h)
	trueneiblist[Ntrueneib++] = j;
    }
  }

  // Now compare each given neighbour with true neighbour list for validation
  for (j=0; j<Ntrueneib; j++) {
    count = 0;
    int jj = trueneiblist[j];
    for (k=0; k<Nneib; k++) if (neiblist[k] == trueneiblist[j]) count++;

    // If the true neighbour is not in the list, or included multiple times, 
    // then output to screen and terminate program
    if (count != 1) {
      cout << "Problem with neighbour lists : " << neibtype << "   " 
	   << i << "   " << jj << "   " << count << "   "
	   << sph->sphdata[i].r[0] << "   " << sph->sphdata[i].h << endl;
      cout << "Nsph  : " << Nsph << "    Ntot : " << Ntot << "    Ntotmax : " << Ntotmax << endl;
      cout << "True Nsph  : " << sph->Nsph << "    True Ntot : " << sph->Ntot << endl;
      cout << "Nneib : " << Nneib << "   Ntrueneib : " << Ntrueneib << endl;
   	  for (k=0; k<ndim; k++) dr[k] = sph->sphdata[jj].r[k] - sph->sphdata[i].r[k];
      drsqd = DotProduct(dr,dr,ndim);
      cout << "drmag : " << sqrt(drsqd) << "     "
	   << sph->kernp->kernrange*sph->sphdata[i].h
	   << "    " << sph->kernp->kernrange*sph->sphdata[jj].h 
	   << "     " << sph->sphdata[jj].h << endl;
      PrintArray("neiblist     : ",Nneib,neiblist);
      PrintArray("trueneiblist : ",Ntrueneib,trueneiblist);
      string message = "Problem with neighbour lists in Binary tree";
      ExceptionHandler::getIstance().raise(message);
    }
  }

  delete[] trueneiblist;

  //cout << "List okay!! : " << Nneib << "   " << Ntrueneib << endl;

  return;
}
#endif



template class BinaryTree2<1>;
template class BinaryTree2<2>;
template class BinaryTree2<3>;
