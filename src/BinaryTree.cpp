//=============================================================================
//  BinaryTree.cpp
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
BinaryTree<ndim>::BinaryTree(int Nleafmaxaux, FLOAT thetamaxsqdaux, 
                               FLOAT kernrangeaux, string gravity_mac_aux,
                               string multipole_aux)
{
  allocated_buffer = false;
  allocated_tree = false;
  neibcheck = true;
  Ntot = 0;
  Ntotmax = 0;
  Ntotmaxold = 0;
  Nleafmax = Nleafmaxaux;
  kernrange = kernrangeaux;
  thetamaxsqd = thetamaxsqdaux;
  gravity_mac = gravity_mac_aux;
  multipole = multipole_aux;
#if defined _OPENMP
  Nthreads = omp_get_max_threads();
#else
  Nthreads = 1;
#endif
}



//=============================================================================
//  BinaryTree::~BinaryTree
/// BinaryTree destructor.  Deallocates tree memory upon object destruction.
//=============================================================================
template <int ndim>
BinaryTree<ndim>::~BinaryTree()
{
  if (allocated_tree) DeallocateTreeMemory();
}



//=============================================================================
//  BinaryTree::AllocateTreeMemory
/// Allocate memory for binary tree as requested.  If more memory is required 
/// than currently allocated, tree is deallocated and reallocated here.
//=============================================================================
template <int ndim>
void BinaryTree<ndim>::AllocateTreeMemory(Sph<ndim> *sph)
{
  int ithread;                      // Thread id number

  debug2("[BinaryTree::AllocateTreeMemory]");

  if (!allocated_tree || Ntotmax > Ntotmaxold) {
    if (allocated_tree) DeallocateTreeMemory();
    Ntotmax = max(Ntotmax,Ntot);
    Ntotmaxold = Ntotmax;

    g2c = new int[gtot];
    ids = new int[Ntotmax];
    inext = new int[Ntotmax];
    tree = new struct BinaryTreeCell<ndim>[Ncellmax];


    if (!allocated_buffer) {
      Nneibmaxbuf = new int[Nthreads];
      Ndirectmaxbuf = new int[Nthreads];
      Ngravcellmaxbuf = new int[Nthreads];
      levelneibbuf = new int*[Nthreads];
      gpotbuf = new FLOAT*[Nthreads];
      divvbuf = new FLOAT*[Nthreads];
      dudtbuf = new FLOAT*[Nthreads];
      abuf = new FLOAT*[Nthreads];
      agravbuf = new FLOAT*[Nthreads];
      activepartbuf = new SphParticle<ndim>*[Nthreads];
      neibpartbuf = new SphParticle<ndim>*[Nthreads];

      for (ithread=0; ithread<Nthreads; ithread++) {
	Nneibmaxbuf[ithread] = max(1,2*sph->Ngather);
	Ndirectmaxbuf[ithread] = max(1,4*sph->Ngather);
	Ngravcellmaxbuf[ithread] = max(1,2*sph->Ngather);
	cout << "Nneibmaxbuf     : " << Nneibmaxbuf[ithread] << endl;
	cout << "Ndirectmaxbuf   : " << Ndirectmaxbuf[ithread] << endl;
	cout << "Ngravcellmaxbuf : " << Ngravcellmaxbuf[ithread] << endl;

	levelneibbuf[ithread] = new int[Ntotmax];
	gpotbuf[ithread] = new FLOAT[Ntotmax];
	divvbuf[ithread] = new FLOAT[Ntotmax];
	dudtbuf[ithread] = new FLOAT[Ntotmax];
	abuf[ithread] = new FLOAT[ndim*Ntotmax];
	agravbuf[ithread] = new FLOAT[ndim*Ntotmax];
	activepartbuf[ithread] = new SphParticle<ndim>[Nleafmax];
	neibpartbuf[ithread] = new SphParticle<ndim>[Nneibmaxbuf[ithread]];
      }
      allocated_buffer = true;
    }



    allocated_tree = true;
  }

  return;
}



//=============================================================================
//  BinaryTree::DeallocateTreeMemory
/// Deallocates all binary tree memory
//=============================================================================
template <int ndim>
void BinaryTree<ndim>::DeallocateTreeMemory(void)
{
  int ithread;                      // Thread id number

  debug2("[BinaryTree::DeallocateTreeMemory]");

  if (allocated_tree) {

    for (ithread=0; ithread<Nthreads; ithread++) {
      delete[] neibpartbuf[ithread];
      delete[] activepartbuf[ithread];
      delete[] agravbuf[ithread];
      delete[] abuf[ithread];
      delete[] dudtbuf[ithread];
      delete[] divvbuf[ithread];
      delete[] gpotbuf[ithread];
      delete[] levelneibbuf[ithread];
    }
    delete[] neibpartbuf;
    delete[] activepartbuf;
    delete[] agravbuf;
    delete[] abuf;
    delete[] dudtbuf;
    delete[] divvbuf;
    delete[] gpotbuf;
    delete[] levelneibbuf;
    delete[] Ngravcellmaxbuf;
    delete[] Ndirectmaxbuf;
    delete[] Nneibmaxbuf;

    delete[] tree;
    delete[] inext;
    delete[] ids;
    delete[] g2c;
    allocated_tree = false;
  }

  return;
}



//=============================================================================
//  BinaryTree::BuildTree
/// Call all routines to build/re-build the binary tree on the local node.
/// If OpenMP is activated, the local domain is partitioned into sub-trees 
/// in order to improve the scalability of building and stocking the tree.
//=============================================================================
template <int ndim>
void BinaryTree<ndim>::BuildTree
(bool rebuild_tree,                 ///< Flag to rebuild tree
 int n,                             ///< Integer time
 int ntreebuildstep,                ///< Tree build frequency
 int ntreestockstep,                ///< Tree stocking frequency
 FLOAT timestep,                    ///< Smallest physical timestep
 Sph<ndim> *sph)                    ///< Pointer to SPH object
{
  debug2("[BinaryTree::BuildTree]");

  // For tree rebuild steps
  //---------------------------------------------------------------------------
  if (n%ntreebuildstep == 0 || rebuild_tree) {

    // Set no. of tree members to total number of SPH particles (inc. ghosts)
    Nsph = sph->Nsph;
    Ntot = sph->Ntot;
    Ntotmaxold = Ntotmax;
    Ntotmax = max(Ntot,Ntotmax);
    Ntotmax = max(Ntotmax,sph->Nsphmax);
    gtot = 0;

    // Compute the size of all tree-related arrays now we know number of points
    ComputeTreeSize();

    // Allocate (or reallocate if needed) all tree memory
    AllocateTreeMemory(sph);

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
//  BinaryTree::ComputeTreeSize
/// Compute the maximum size (i.e. no. of levels, cells and leaf cells) of 
/// the binary tree.
//=============================================================================
template <int ndim>
void BinaryTree<ndim>::ComputeTreeSize(void)
{
  debug2("[BinaryTree::ComputeTreeSize]");

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
//  BinaryTree::CreateTreeStructure
/// Create the raw tree skeleton structure once the tree size is known.
/// Sets all cell pointer variables and all cell levels.
//=============================================================================
template <int ndim>
void BinaryTree<ndim>::CreateTreeStructure(void)
{
  int c;                            // Dummy id of tree-level, then tree-cell
  int g;                            // Dummy id of grid-cell
  int l;                            // Dummy id of level
  int *c2L;                         // Increment to second child-cell
  int *cNL;                         // Increment to next cell if cell unopened

  debug2("[BinaryTree::CreateTreeStructure]");

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
    tree[c].Nactive = 0;
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
//  BinaryTree::LoadParticlesToTree
/// Create tree structure by adding particles to leaf cells.
//=============================================================================
template <int ndim>
void BinaryTree<ndim>::LoadParticlesToTree
(Sph<ndim> *sph)                   ///< Pointer to main SPH object
{
  int i;                           // SPH particle counter
  int k;                           // Dimensionality counter

  debug2("[BinaryTree::LoadParticlesToTree]");

  // Set root cell counters and properties
  tree[0].N = Ntot;
  tree[0].ifirst = 0;
  tree[0].ilast = Ntot - 1;
  for (k=0; k<ndim; k++) tree[0].bbmin[k] = -big_number;
  for (k=0; k<ndim; k++) tree[0].bbmax[k] = big_number;

  for (i=0; i<Ntot; i++) ids[i] = i;
  for (i=0; i<Ntot-1; i++) inext[i] = -1;
  
  // Recursively build tree from root node down
  DivideTreeCell(0,Ntot-1,sph,tree[0]);

  return;
}



//=============================================================================
//  BinaryTree::DivideTreeCell
/// Recursive routine to divide a tree cell into two children cells.
//=============================================================================
template <int ndim>
void BinaryTree<ndim>::DivideTreeCell
(int ifirst,                        ///< Aux. id of first particle in cell
 int ilast,                         ///< Aux. id of last particle in cell
 Sph<ndim> *sph,                    ///< Pointer to main SPH object
 BinaryTreeCell<ndim> &cell)        ///< Cell to be divided
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
  if (pow(2,cell.level) <= Nthreads) {
#pragma omp parallel for default(none) private(i) shared(cell,ifirst,ilast,sph) num_threads(2)
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
//  BinaryTree::QuickSelect
/// Find median and sort particles in arrays to ensure they are the correct 
/// side of the division.  Uses the QuickSelect algorithm.
//=============================================================================
template <int ndim>
FLOAT BinaryTree<ndim>::QuickSelect
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
//  BinaryTree::StockTree
/// ..
//=============================================================================
template <int ndim>
void BinaryTree<ndim>::StockTree
(BinaryTreeCell<ndim> &cell,       ///< Reference to cell to be stocked
 SphParticle<ndim> *sphdata)        ///< SPH particle data array
{
  int i;                            // Aux. child cell counter

  // If cell is not leaf, stock child cells
  if (cell.level != ltot) {
#if defined _OPENMP
    if (pow(2,cell.level) <= Nthreads) {
#pragma omp parallel for default(none) private(i) shared(cell,sphdata) num_threads(2)
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
//  BinaryTree::StockCellProperties
/// Calculate the physical properties (e.g. total mass, centre-of-mass, 
/// opening-distance, etc..) of all cells in the tree.
//=============================================================================
template <int ndim>
void BinaryTree<ndim>::StockCellProperties
(BinaryTreeCell<ndim> &cell,       ///< ..
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
  FLOAT p = 0.0;
  FLOAT lambda = 0.0;
  BinaryTreeCell<ndim> &child1 = tree[cell.c1];
  BinaryTreeCell<ndim> &child2 = tree[cell.c2];


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
  //for (k=0; k<ndim; k++) cell.hboxmin[k] = big_number;
  //for (k=0; k<ndim; k++) cell.hboxmax[k] = -big_number;
  for (k=0; k<5; k++) cell.q[k] = 0.0;
  for (k=0; k<ndim; k++) cell.rcell[k] = 0.0;

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
	//if (sphdata[i].r[k] - kernrange*sphdata[i].h < cell.hboxmin[k]) 
	//  cell.hboxmin[k] = sphdata[i].r[k] - kernrange*sphdata[i].h;
	//if (sphdata[i].r[k] + kernrange*sphdata[i].h > cell.hboxmax[k])
	//  cell.hboxmax[k] = sphdata[i].r[k] + kernrange*sphdata[i].h;
      }
      if (i == cell.ilast) break;
      i = inext[i];
    };
    
    // Normalise all cell values
    if (cell.N > 0) {
      for (k=0; k<ndim; k++) cell.r[k] /= cell.m;
      for (k=0; k<ndim; k++) cell.v[k] /= cell.m;
      for (k=0; k<ndim; k++) dr[k] = 0.5*(cell.bbmax[k] - cell.bbmin[k]);
      cell.cdistsqd = max(factor*DotProduct(dr,dr,ndim),cell.hmax*cell.hmax);
      for (k=0; k<ndim; k++) cell.rcell[k] = cell.bbmin[k] + dr[k];
      //for (k=0; k<ndim; k++) 
      //dr[k] = max(cell.bbmax[k] - cell.r[k],cell.r[k] - cell.bbmin[k]);
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
      //for (k=0; k<ndim; k++)
      //cell.hboxmin[k] = min(child1.hboxmin[k],cell.hboxmin[k]);
      //for (k=0; k<ndim; k++)
      //cell.hboxmax[k] = max(child1.hboxmax[k],cell.hboxmax[k]);
      cell.hmax = max(cell.hmax,child1.hmax);
    }
    if (child2.N > 0) {
      for (k=0; k<ndim; k++)
	cell.bbmin[k] = min(child2.bbmin[k],cell.bbmin[k]);
      for (k=0; k<ndim; k++)
	cell.bbmax[k] = max(child2.bbmax[k],cell.bbmax[k]);
      //for (k=0; k<ndim; k++)
      //cell.hboxmin[k] = min(child2.hboxmin[k],cell.hboxmin[k]);
      //for (k=0; k<ndim; k++)
      //cell.hboxmax[k] = max(child2.hboxmax[k],cell.hboxmax[k]);
      cell.hmax = max(cell.hmax,child2.hmax);
    }


    if (cell.N > 0) {
      cell.m = child1.m + child2.m;
      for (k=0; k<ndim; k++) cell.r[k] =
	(child1.m*child1.r[k] + child2.m*child2.r[k])/cell.m;
      for (k=0; k<ndim; k++) cell.v[k] =
	(child1.m*child1.v[k] + child2.m*child2.v[k])/cell.m;
      for (k=0; k<ndim; k++) 
	dr[k] = 0.5*(cell.bbmax[k] - cell.bbmin[k]);
      for (k=0; k<ndim; k++) cell.rcell[k] = cell.bbmin[k] + dr[k];
      cell.cdistsqd = max(factor*DotProduct(dr,dr,ndim),cell.hmax*cell.hmax);
      //for (k=0; k<ndim; k++) dr[k] = max(cell.bbmax[k] - cell.r[k],
      //				 cell.r[k] - cell.bbmin[k]);
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


    // Calculate eigenvalue MAC criteria
    //abserror = 0.00001;
    //if (ndim == 3)
    //  p = cell.q[0]*cell.q[2] - (cell.q[0] + cell.q[2])*(cell.q[0] + cell.q[2])
    //	- cell.q[1]*cell.q[1] - cell.q[3]*cell.q[3] - cell.q[4]*cell.q[4];
    //if (p >= 0.0) cell.mac = 0.0;
    //else {
    //  lambda = 2.0*sqrt(-p/3.0);
    //  cell.mac = pow(0.5*lambda/abserror,0.66666666666666);
    //}

    
  }
  //---------------------------------------------------------------------------


  return;
}



//=============================================================================
//  BinarySubTree::ExtrapolateCellProperties
/// Extrapolate important physical properties of all cells in the tree.
//=============================================================================
template <int ndim>
void BinaryTree<ndim>::ExtrapolateCellProperties
(FLOAT dt)                          ///< Smallest timestep size
{
  int c;                            // ..
  int k;                            // ..

  debug2("[BinaryTree::ExtrapolateCellProperties]");


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
void BinaryTree<ndim>::UpdateHmaxValues
(SphParticle<ndim> *sphdata)        ///< SPH particle data array
{
  int c,cc,ccc;                     // Cell counters
  int i;                            // Particle counter
  int j;
  int k;

  debug2("[BinaryTree::UpdateHmaxValues]");

  // Zero all summation variables for all cells
  for (c=0; c<Ncellmax; c++) tree[c].hmax = 0.0;
  //for (c=0; c<Ncellmax; c++) {
  //  for (k=0; k<ndim; k++) tree[c].hboxmin[k] = big_number;
  // for (k=0; k<ndim; k++) tree[c].hboxmax[k] = -big_number;
  //}

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
	//for (k=0; k<ndim; k++) {
	//  if (sphdata[i].r[k] - kernrange*sphdata[i].h < tree[c].hboxmin[k]) 
	//    tree[c].hboxmin[k] = sphdata[i].r[k] - kernrange*sphdata[i].h;
	//  if (sphdata[i].r[k] + kernrange*sphdata[i].h > tree[c].hboxmax[k])
	//    tree[c].hboxmax[k] = sphdata[i].r[k] + kernrange*sphdata[i].h;
	//}
	if (i == tree[c].ilast) break;
        i = inext[i];
      };

    }
    // For non-leaf cells, sum together two children cells
    //-------------------------------------------------------------------------
    else {
      cc = tree[c].c1;
      ccc = tree[c].c2;
      if (tree[cc].N > 0) {
	tree[c].hmax = max(tree[c].hmax,tree[cc].hmax);
	//for (k=0; k<ndim; k++)
	//  tree[c].hboxmin[k] = min(tree[c].hboxmin[k],tree[cc].hboxmin[k]);
	//for (k=0; k<ndim; k++)
	//  tree[c].hboxmax[k] = max(tree[c].hboxmax[k],tree[cc].hboxmax[k]);
      }
      if (tree[ccc].N > 0) {
	tree[c].hmax = max(tree[c].hmax,tree[ccc].hmax);
	//for (k=0; k<ndim; k++)
	//  tree[c].hboxmin[k] = min(tree[c].hboxmin[k],tree[ccc].hboxmin[k]);
	//for (k=0; k<ndim; k++)
	//  tree[c].hboxmax[k] = max(tree[c].hboxmax[k],tree[ccc].hboxmax[k]);
      }

    }
    //-------------------------------------------------------------------------

  }
  //===========================================================================

  return;
}



//=============================================================================
//  BinaryTree::UpdateActiveParticleCounters
/// Loop through all leaf cells in binary tree and update all active 
/// particle counters.
//=============================================================================
template <int ndim>
void BinaryTree<ndim>::UpdateActiveParticleCounters
(Sph<ndim> *sph)                    ///< Pointer to main SPH object
{
  int c;
  int i;
  int ilast;
  int cactive=0;
  int Nptcl=0;
  //int *pcount;

  //pcount = new int[sph->Ntot];
  //for (i=0; i<sph->Ntot; i++) pcount[i] = 0;

  debug2("[BinaryTree::UpdateActiveParticleCounters]");

#if defined(VERIFY_ALL)
  ValidateTree(sph);
#endif



  if (!(inext[0] >= -1 && inext[0] < sph->Ntot)) {
    cout << "inext : " << inext[0] << endl;
    cin >> i;
  }

  // Loop through all grid cells in turn
  //---------------------------------------------------------------------------
  for (c=0; c<Ncell; c++) {
    tree[c].Nactive = 0;

    if (tree[c].level != ltot) continue;
    i = tree[c].ifirst;
    ilast = tree[c].ilast;
    //cout << "c : " << c << "    " << tree[c].N << "    " << i << "    " << ilast << endl;

    // Else walk through linked list to obtain list and number of active ptcls.
    while (i != -1) {
      //cout << "i : " << i << "   " << inext[i] << "    " << ilast << endl;
      //pcount[i]++;
      if (sph->sphdata[i].active) tree[c].Nactive++;
      if (i == ilast) break;
      i = inext[i];
    };
    //cout << "Next cell? : " << i << "    " << ilast << endl;

  }
  //---------------------------------------------------------------------------

  //for (c=0; c<Ncell; c++) {
  //  if (tree[c].level == ltot) Nptcl += tree[c].N;
  //  if (tree[c].Nactive > 0) cactive += tree[c].Nactive;
  //}
  //cout << "CACTIVE2? : " << cactive << "    " << Nptcl << endl;

  //for (i=0; i<sph->Ntot; i++) {
  //  if (pcount[i] != 1) {
  //    cout << "PCOUNT PROBLEM : " << i << "   " << pcount[i] << endl;
  //    exit(0);
  //  }
  // }

  return;
}



//=============================================================================
//  BinaryTree::ComputeActiveParticleList
/// Returns the number (Nactive) and list of ids (activelist) of all active
/// SPH particles in the given cell.
//=============================================================================
template <int ndim>
int BinaryTree<ndim>::ComputeActiveParticleList
(BinaryTreeCell<ndim> *cell,       ///< [in] Pointer to cell
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
//  BinaryTree::BoxOverlap
/// ..
//=============================================================================
template <int ndim>
bool BinaryTree<ndim>::BoxOverlap
(FLOAT box1min[ndim],
 FLOAT box1max[ndim],
 FLOAT box2min[ndim],
 FLOAT box2max[ndim])
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



//=============================================================================
//  BinaryTree::ComputeGatherNeighbourList
/// Computes and returns number of neighbour, 'Nneib', and the list
/// of neighbour ids, 'neiblist', for all particles inside cell 'c'.
/// Includes all particles in the selected cell, plus all particles
/// contained in adjacent cells (including diagonal cells).
/// Wrapper around the true implementation inside BinarySubTree
//=============================================================================
template <int ndim>
int BinaryTree<ndim>::ComputeGatherNeighbourList
(BinaryTreeCell<ndim> *cell,       ///< [in] Pointer to current cell
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
  FLOAT gatherboxmin[ndim];
  FLOAT gatherboxmax[ndim];

  for (k=0; k<ndim; k++) rc[k] = cell->rcell[k];
  hrangemax = cell->rmax + kernrange*hmax;

  for (k=0; k<ndim; k++) gatherboxmin[k] = cell->bbmin[k] - kernrange*hmax;
  for (k=0; k<ndim; k++) gatherboxmax[k] = cell->bbmax[k] + kernrange*hmax;

  // Start with root cell and walk through entire tree
  cc = 0;

  //===========================================================================
  while (cc < Ncell) {

    for (k=0; k<ndim; k++) dr[k] = tree[cc].rcell[k] - rc[k];
    drsqd = DotProduct(dr,dr,ndim);

    // Check if circular range overlaps cell bounding sphere
    //-------------------------------------------------------------------------
    if (drsqd < pow(tree[cc].rmax + hrangemax,2)) {

    
    // Check if bounding boxes overlap with each other
    //if (BoxOverlap(gatherboxmin,gatherboxmax,tree[cc].bbmin,tree[cc].bbmax)) {


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
//  BinaryTree::ComputeNeighbourList
/// Computes and returns number of neighbour, 'Nneib', and the list
/// of neighbour ids, 'neiblist', for all particles inside cell 'c'.
/// Includes all particles in the selected cell, plus all particles
/// contained in adjacent cells (including diagonal cells).
/// Wrapper around the true implementation inside BinarySubTree.
/// If allocated memory array containing neighbour ids (neiblist) overflows, 
/// return with error code (-1) in order to reallocate more memory.
//=============================================================================
template <int ndim>
int BinaryTree<ndim>::ComputeNeighbourList
(BinaryTreeCell<ndim> *cell,       ///< [in] Cell pointer
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

    // Check if bounding boxes overlap with each other
    //if (BoxOverlap(cell->bbmin,cell->bbmax,tree[cc].hboxmin,tree[cc].hboxmax)
    //|| BoxOverlap(cell->hboxmin,cell->hboxmax,tree[cc].bbmin,tree[cc].bbmax)){

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
int BinaryTree<ndim>::ComputeGravityInteractionList
(BinaryTreeCell<ndim> *cell,       ///< [in] Pointer to cell
 int Nneibmax,                      ///< [in] Max. no. of SPH neighbours
 int Ndirectmax,                    ///< [in] Max. no. of direct-sum neighbours
 int Ngravcellmax,                  ///< [in] Max. no. of cell interactions
 int &Nneib,                        ///< [out] No. of SPH neighbours
 int &Ndirect,                      ///< [out] No. of direct-sum neighbours
 int &Ngravcell,                    ///< [out] No. of cell interactions
 int *neiblist,                     ///< [out] List of SPH neighbour ids
 int *directlist,                   ///< [out] List of direct-sum neighbour ids
 BinaryTreeCell<ndim> **gravcelllist,  ///< [out] List of cell ids
 SphParticle<ndim> *sphdata)        ///< [in] SPH particle data
{
  int cc;                           // Cell counter
  int i;                            // Particle id
  int j;                            // Aux. particle counter
  int k;                            // Neighbour counter
  int Nneibtemp = 0;                // Aux. counter
  FLOAT cdistsqd;                   // ..
  FLOAT dr[ndim];                   // Relative position vector
  FLOAT drsqd;                      // Distance squared
  FLOAT rc[ndim];                   // Position of cell
  FLOAT hrangemax;                  // Maximum kernel extent 
  FLOAT rmax;                       // Radius of sphere containing particles

  // Make local copies of important cell properties
  for (k=0; k<ndim; k++) rc[k] = cell->rcell[k];
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


    for (k=0; k<ndim; k++) dr[k] = tree[cc].rcell[k] - rc[k];
    drsqd = DotProduct(dr,dr,ndim);

    // Check if cells contain SPH neighbours
    //-------------------------------------------------------------------------
    if (drsqd < (tree[cc].rmax + hrangemax)*(tree[cc].rmax + hrangemax) ||
	drsqd < pow(tree[cc].rmax + rmax + kernrange*tree[cc].hmax,2)) {


    // Check if bounding boxes overlap with each other
    //if (BoxOverlap(cell->bbmin,cell->bbmax,tree[cc].hboxmin,tree[cc].hboxmax) 
	//|| BoxOverlap(cell->hboxmin,cell->hboxmax,tree[cc].bbmin,tree[cc].bbmax)){

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
//  BinaryTree::ComputeCellMonopoleForces
/// Compute the force on particle 'parti' due to all cells obtained in the 
/// gravity tree walk.  Uses only monopole moments (i.e. COM) of the cell.
//=============================================================================
template <int ndim>
void BinaryTree<ndim>::ComputeCellMonopoleForces
(int i,                               ///< [in] i.d. of particle
 int Ngravcell,                       ///< [in] No. of tree cells in list
 BinaryTreeCell<ndim> **gravcelllist, ///< [in] List of tree cell ids
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
  BinaryTreeCell<ndim> *cell;       // Pointer to gravity tree cell

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
//  BinaryTree::ComputeCellQuadrupoleForces
/// Compute the force on particle 'parti' due to all cells obtained in the 
/// gravity tree walk including the quadrupole moment correction term.
//=============================================================================
template <int ndim>
void BinaryTree<ndim>::ComputeCellQuadrupoleForces
(int i,                               ///< [in] i.d. of particle
 int Ngravcell,                       ///< [in] No. of tree cells in list
 BinaryTreeCell<ndim> **gravcelllist, ///< [in] List of tree cell ids
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
  BinaryTreeCell<ndim> *cell;       // Pointer to gravity tree cell

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
int BinaryTree<ndim>::ComputeActiveCellList
(BinaryTreeCell<ndim> **celllist) ///< Cells id array containing active ptcls
{
  int c;                           // Cell counter
  int Nactive = 0;                 // ..

  for (c=0; c<Ncell; c++)
    if (tree[c].Nactive > 0) celllist[Nactive++] = &tree[c];

  return Nactive;
}



//=============================================================================
//  BinaryTree::UpdateAllSphProperties
/// Compute all local 'gather' properties of currently active particles, and 
/// then compute each particle's contribution to its (active) neighbour 
/// neighbour hydro forces.  Optimises the algorithm by using grid-cells to 
/// construct local neighbour lists for all particles  inside the cell.
//=============================================================================
template <int ndim>
void BinaryTree<ndim>::UpdateAllSphProperties
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
  BinaryTreeCell<ndim> *cell;      // Pointer to binary tree cell
  BinaryTreeCell<ndim> **celllist; // List of binary cell pointers
  SphParticle<ndim> *data = sph->sphdata;  // Pointer to SPH particle data

  int Nneibcount = 0;
  int ithread;

  debug2("[BinaryTree::UpdateAllSphProperties]");

  // Find list of all cells that contain active particles
  celllist = new BinaryTreeCell<ndim>*[gtot];
  cactive = ComputeActiveCellList(celllist);

  //cout << "CACTIVE : " << cactive << "    " << Ncell << "    " << gtot << endl;
  //cout << "NSPH : " << sph->Nsph << "    " << sph->Ntot << endl;
  //if (cactive > sph->Nsph) exit(0);


  // Set-up all OMP threads
  //===========================================================================
#pragma omp parallel default(none) private(activelist,cc,cell,celldone,draux)\
  private(drsqd,drsqdaux,hmax,hrangesqd,i,ithread,j,jj,k,okflag,m,mu)\
  private(Nactive,neiblist,Nneib,Nneibmax,r,rp,gpot,gpot2,m2,mu2,Ngather) \
  shared(sph,celllist,cactive,cout,data,nbody) reduction(+:Nneibcount)
  {
#if defined _OPENMP
    ithread = omp_get_thread_num();
#else
    ithread = 0;
#endif
    Nneibmax = Nneibmaxbuf[ithread];
    //cout << "Nneibmax : " << Nneibmax << endl;

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
	  //cout << "Resizing arrays? : " << Nneib << "   " << Nneibmax 
	  //     << "    " << cactive << "    " << Nactive << endl;
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
          Nneibmaxbuf[ithread] = Nneibmax;
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
	  //if (Nneib != -1) cout << "Done!! : " << Nneib << "   " << Nneibmax << endl;
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
	  Nneibcount += Nneib;
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

	    //cout << "MEH?? : " << Ngather << "   " << Nneib << "    " << drsqdaux << "     " << hrangesqd << endl;

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

	  if (Ngather == 0) {
	    cout << "WTF?? : " << Ngather << "    " << Nneib << "    " 
		 << hrangesqd << endl;
	    exit(0);
	  }


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

  //cout << "Average Ngather : " << Nneibcount/sph->Nsph << endl;


  // Update all tree smoothing length values
  UpdateHmaxValues(sph->sphdata);


#if defined(VERIFY_ALL)
  ValidateTree(sph);
#endif

  return;
}



//=============================================================================
//  BinaryTree::UpdateAllSphHydroForces
/// Compute all local 'gather' properties of currently active particles, and 
/// then compute each particle's contribution to its (active) neighbour 
/// neighbour hydro forces.  Optimises the algorithm by using grid-cells to 
/// construct local neighbour lists for all particles  inside the cell.
//=============================================================================
template <int ndim>
void BinaryTree<ndim>::UpdateAllSphHydroForces
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
  BinaryTreeCell<ndim> *cell;      // Pointer to binary tree cell
  BinaryTreeCell<ndim> **celllist; // List of binary tree pointers
  SphParticle<ndim> *data = sph->sphdata;   // Pointer to SPH particle data

  int ithread;
  int Nneibcount = 0;
  SphParticle<ndim> *activepart;
  SphParticle<ndim> *neibpart;

  debug2("[BinaryTree::UpdateAllSphHydroForces]");

  // Update tree smoothing length values here
  UpdateHmaxValues(sph->sphdata);


  // Find list of all cells that contain active particles
  celllist = new BinaryTreeCell<ndim>*[gtot];
  cactive = ComputeActiveCellList(celllist);


  // Set-up all OMP threads
  //===========================================================================
#pragma omp parallel default(none) private(activelist,activepart,cc,cell,dr)\
  private(draux,drmag,drsqd,hrangesqdi,i,interactlist,ithread,invdrmag,j,jj,k)\
  private(Nactive,neiblist,neibpart,Ninteract,Nneib,Nneibmax,rp)\
  shared(cactive,celllist,data,sph) reduction(+:Nneibcount)
  {
#if defined _OPENMP
    ithread = omp_get_thread_num();
#else
    ithread = 0;
#endif
    Nneibmax = Nneibmaxbuf[ithread];
    activepart = activepartbuf[ithread];
    //neibpart = neibpartbuf[ithread];
    neibpart = new SphParticle<ndim>[Nneibmax];

    activelist = new int[Nleafmax];
    neiblist = new int[Nneibmax];
    interactlist = new int[Nneibmax];
    dr = new FLOAT[Nneibmax*ndim];
    drmag = new FLOAT[Nneibmax];
    invdrmag = new FLOAT[Nneibmax];

    // Loop over all active cells
    //=========================================================================
#pragma omp for schedule(dynamic)
    for (cc=0; cc<cactive; cc++) {
      cell = celllist[cc];

     
      //#if defined(VERIFY_ALL)
      //cout << "Cell 0 : " << cc << "    " << cactive << endl;
      //ValidateTree(sph);
      //#endif

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

      //#if defined(VERIFY_ALL)
      //cout << "Cell 0.6 : " << cc << "    " << celllist[cc]->id << endl;
      //ValidateTree(sph);
      //#endif

      // If there are too many neighbours, reallocate the arrays and
      // recompute the neighbour list.
      while (Nneib == -1) {
	cout << "Resizing arrays : " << Nneibmax << endl;
        delete[] neibpartbuf[ithread];
        delete[] invdrmag;
        delete[] drmag;
        delete[] dr;
        delete[] interactlist;
        delete[] neiblist;
        Nneibmax = 2*Nneibmax;
	Nneibmaxbuf[ithread] = Nneibmax;
        neiblist = new int[Nneibmax];
        interactlist = new int[Nneibmax];
        dr = new FLOAT[Nneibmax*ndim];
        drmag = new FLOAT[Nneibmax];
        invdrmag = new FLOAT[Nneibmax];
        //neibpartbuf[ithread] = new SphParticle<ndim>[Nneibmax]; neibpart = neibpartbuf[ithread];
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
	//#if defined(VERIFY_ALL)
	//cout << "Cell 0.8 : " << cc << "    " << j << "   " << Nneib << "    " << Nneibmax << "    " << celllist[cc]->id << endl;
	//cout << &(neibpart[j]) << "    " << (&ids[0]) << "    " << &(inext[0]) << endl;
	//ValidateTree(sph);
	//#endif
      }


      //#if defined(VERIFY_ALL)
      //cout << "Cell 1 : " << cc << "    " << celllist[cc]->id << endl;
      //ValidateTree(sph);
      //#endif

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

	Nneibcount += Ninteract;

      }
      //-----------------------------------------------------------------------

      //#if defined(VERIFY_ALL)
      //cout << "Cell 2 : " << cc << "    " << celllist[cc]->id << endl;
      //ValidateTree(sph);
      //#endif


      // Add all active particles contributions to main array
      for (j=0; j<Nactive; j++) {
    	i = activelist[j];
        for (k=0; k<ndim; k++) {
          data[i].a[k] += activepart[j].a[k];
        }
        data[i].dudt += activepart[j].dudt;
        data[i].div_v += activepart[j].div_v;
        data[i].levelneib = max(data[i].levelneib,activepart[j].levelneib);
      }

      // Now add all active neighbour contributions to main array
      for (jj=0; jj<Nneib; jj++) {
        j = neiblist[jj];
        if (neibpart[jj].active) {
          for (k=0; k<ndim; k++) {
            data[j].a[k] += neibpart[jj].a[k];
          }
          data[j].dudt += neibpart[jj].dudt;
          data[j].div_v += neibpart[jj].div_v;
        }
        data[j].levelneib = max(data[j].levelneib,neibpart[jj].levelneib);
      }

    }
    //=========================================================================


#if defined(VERIFY_ALL)
    cout << "Cell 3 : " << cc << "    " << cactive << endl;
    ValidateTree(sph);
#endif
    
    //cout << "Addresses : " << invdrmag << "    " << drmag << "     " << dr << endl;
    
    // Free-up local memory for OpenMP thread
    delete[] neibpart;
    delete[] invdrmag;
    delete[] drmag;
    delete[] dr;
    delete[] interactlist;
    delete[] neiblist;
    delete[] activelist;

  }
  //===========================================================================

  delete[] celllist;

  //cout << "Average Nneib     : " << Nneibcount/sph->Nsph << endl;

  // Compute other important SPH quantities after hydro forces are computed
  if (sph->hydro_forces == 1) {
    for (i=0; i<sph->Nsph; i++) {
      if (sph->sphdata[i].active)
    	  sph->ComputePostHydroQuantities(sph->sphdata[i]);
    }
  }

  // Update tree smoothing length values here
  //UpdateHmaxValues(sph->sphdata);


#if defined(VERIFY_ALL)
  ValidateTree(sph);
#endif

  return;
}



//=============================================================================
//  BinaryTree::UpdateAllSphForces
/// Compute all local 'gather' properties of currently active particles, and 
/// then compute each particle's contribution to its (active) neighbour 
/// neighbour hydro forces.  Optimises the algorithm by using grid-cells to 
/// construct local neighbour lists for all particles  inside the cell.
//=============================================================================
template <int ndim>
void BinaryTree<ndim>::UpdateAllSphForces
(Sph<ndim> *sph)                    ///< Pointer to SPH object
{
  int cactive;                      // No. of active cells
  int cc;                           // Aux. cell counter
  int i;                            // Particle id
  int ithread;
  int j;                            // Aux. particle counter
  int jj;                           // Aux. particle counter
  int k;                            // Dimension counter
  int okflag;                       // Flag if h-rho iteration is valid
  int Nactive;                      // No. of active particles in cell
  int Ndirect;                      // No. of direct-sum gravity particles
  int Ndirectaux;                   // ..
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
  BinaryTreeCell<ndim> *cell;      // Pointer to binary tree cell
  BinaryTreeCell<ndim> **celllist; // List of pointers to binary tree cells
  BinaryTreeCell<ndim> **gravcelllist; // List of pointers to grav. cells

  FLOAT draux[ndim];               // Aux. relative position vector
  FLOAT drsqd;                     // Distance squared
  FLOAT hrangesqdi;                // Kernel gather extent
  FLOAT rp[ndim];                  //

  SphParticle<ndim> *activepart;
  SphParticle<ndim> *neibpart;
  int *levelneib;
  FLOAT *gpot;
  FLOAT *div_v;
  FLOAT *dudt;
  FLOAT *a;
  FLOAT *agrav;

  int Nneibcount = 0;
  int Ndirectcount = 0;
  int Ncellcount = 0;

  debug2("[BinaryTree::UpdateAllSphForces]");


  // Update tree smoothing length values here
  UpdateHmaxValues(sph->sphdata);


  // Find list of all cells that contain active particles
  celllist = new BinaryTreeCell<ndim>*[gtot];
  cactive = ComputeActiveCellList(celllist);


  // Set-up all OMP threads
  //===========================================================================
#pragma omp parallel default(none) shared(celllist,cactive,sph,cout)\
  private(activelist,agrav,cc,cell,draux,drsqd,hrangesqdi,Ndirectaux,rp) \
  private(gpot,i,interactlist,ithread,j,jj,activepart,a,div_v,dudt,levelneib)\
  private(k,okflag,Nactive,neiblist,neibpart,Ninteract,Nneib,directlist)\
  private(gravcelllist,Ngravcell,Ndirect,Nneibmax,Ndirectmax,Ngravcellmax)\
  reduction(+:Ncellcount,Ndirectcount,Nneibcount)
  {
#if defined _OPENMP
    ithread = omp_get_thread_num();
#else
    ithread = 0;
#endif
    Nneibmax = Nneibmaxbuf[ithread];
    Ndirectmax = Ndirectmaxbuf[ithread];
    Ngravcellmax = Ngravcellmaxbuf[ithread];

    levelneib = levelneibbuf[ithread];
    gpot = gpotbuf[ithread];
    div_v = divvbuf[ithread];
    dudt = dudtbuf[ithread];
    a = abuf[ithread];
    agrav = agravbuf[ithread];
    activepart = activepartbuf[ithread];
    neibpart = neibpartbuf[ithread];

    activelist = new int[Nleafmax];
    neiblist = new int[Nneibmax];
    interactlist = new int[Nneibmax];
    directlist = new int[Ndirectmax];
    gravcelllist = new BinaryTreeCell<ndim>*[Ngravcellmax];


    // Zero temporary grav. accel array
    for (i=0; i<sph->Nsph; i++) gpot[i] = 0.0;
    for (i=0; i<sph->Nsph; i++) div_v[i] = 0.0;
    for (i=0; i<sph->Nsph; i++) dudt[i] = 0.0;
    for (i=0; i<ndim*sph->Nsph; i++) a[i] = 0.0;
    for (i=0; i<ndim*sph->Nsph; i++) agrav[i] = 0.0;


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
        delete[] neibpartbuf[ithread];
        delete[] gravcelllist;
        delete[] directlist;
        delete[] interactlist;
        delete[] neiblist;
        Nneibmax = 2*Nneibmax;
        Ndirectmax = 2*Ndirectmax;
        Ngravcellmax = 2*Ngravcellmax;
	Nneibmaxbuf[ithread] = Nneibmax;
	Ndirectmaxbuf[ithread] = Ndirectmax;
	Ngravcellmaxbuf[ithread] = Ngravcellmax;
        neiblist = new int[Nneibmax];
        interactlist = new int[Nneibmax];
        directlist = new int[Ndirectmax];
        gravcelllist = new BinaryTreeCell<ndim>*[Ngravcellmax];
        neibpartbuf[ithread] = new SphParticle<ndim>[Nneibmax]; neibpart = neibpartbuf[ithread];
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
	Ndirectaux = Ndirect;
        for (k=0; k<ndim; k++) rp[k] = activepart[j].r[k]; //data[i].r[k];
        hrangesqdi = activepart[j].hrangesqd;

        for (jj=0; jj<Nneib; jj++) {

          for (k=0; k<ndim; k++) draux[k] = neibpart[jj].r[k] - rp[k];
          drsqd = DotProduct(draux,draux,ndim) + small_number;

          // Compute relative position and distance quantities for pair
          if (drsqd > hrangesqdi && drsqd >= neibpart[jj].hrangesqd)
	    directlist[Ndirectaux++] = neiblist[jj];
	  else if ((neiblist[jj] < i && !neibpart[jj].active) || neiblist[jj] > i)
            interactlist[Ninteract++] = jj;
        }

        // Compute forces between SPH neighbours (hydro and gravity)
        sph->ComputeSphHydroGravForces(i,Ninteract,interactlist,
                                       activepart[j],neibpart);

        // Compute direct gravity forces between distant particles
        sph->ComputeDirectGravForces(i,Ndirectaux,directlist,
                                     agrav,gpot,activepart[j],sph->sphdata);

        // Compute gravitational force due to distant cells
        if (multipole == "monopole")
          ComputeCellMonopoleForces(i,Ngravcell,gravcelllist,activepart[j]);
        else if (multipole == "quadrupole")
          ComputeCellQuadrupoleForces(i,Ngravcell,gravcelllist,activepart[j]);

	Nneibcount += Ninteract;
	Ndirectcount += Ndirectaux;
	Ncellcount += Ngravcell;

      }
      //-----------------------------------------------------------------------


      // Add all active particles contributions to main array
      for (j=0; j<Nactive; j++) {
    	i = activelist[j];
        for (k=0; k<ndim; k++) a[ndim*i + k] += activepart[j].a[k];
        for (k=0; k<ndim; k++) agrav[ndim*i + k] += activepart[j].agrav[k];
        gpot[i] += activepart[j].gpot;
        dudt[i] += activepart[j].dudt;
        div_v[i] += activepart[j].div_v;
        levelneib[i] = max(levelneib[i],activepart[j].levelneib);
      }

      // Now add all active neighbour contributions to the main arrays
      for (jj=0; jj<Nneib; jj++) {
        j = neiblist[jj];
        if (neibpart[jj].active) {
          for (k=0; k<ndim; k++) a[ndim*j + k] += neibpart[jj].a[k];
          for (k=0; k<ndim; k++) agrav[ndim*j + k] += neibpart[jj].agrav[k];
          gpot[j] += neibpart[jj].gpot;
          dudt[j] += neibpart[jj].dudt;
          div_v[j] += neibpart[jj].div_v;
        }
        levelneib[j] = max(levelneib[j],neibpart[jj].levelneib);
      }

    }
    //=========================================================================


    // Finally, add all contributions from distant pair-wise forces to arrays
#pragma omp critical
    {
      for (i=0; i<sph->Nsph; i++) {
	if (sph->sphdata[i].active) {
	  sph->sphdata[i].levelneib = 
	    max(sph->sphdata[i].levelneib,levelneib[i]);
	  sph->sphdata[i].gpot += gpot[i];
	  sph->sphdata[i].div_v += div_v[i];
	  sph->sphdata[i].dudt += dudt[i];
	  for (k=0; k<ndim; k++) sph->sphdata[i].a[k] += a[ndim*i + k];
	  for (k=0; k<ndim; k++) sph->sphdata[i].agrav[k] += agrav[ndim*i + k];
	}
      }
    }


    // Free-up local memory for OpenMP thread
    delete[] gravcelllist;
    delete[] directlist;
    delete[] interactlist;
    delete[] neiblist;
    delete[] activelist;

  }
  //===========================================================================

  delete[] celllist;

  //cout << "Average Nneib     : " << Nneibcount/sph->Nsph << endl;
  //cout << "Average Ndirect   : " << Ndirectcount/sph->Nsph << endl;
  //cout << "Average Ngravcell : " << Ncellcount/sph->Nsph << endl;

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
//  BinaryTree::UpdateAllSphGravityProperties
/// Compute all local 'gather' properties of currently active particles, and 
/// then compute each particle's contribution to its (active) neighbour 
/// neighbour hydro forces.  Optimises the algorithm by using grid-cells to 
/// construct local neighbour lists for all particles  inside the cell.
//=============================================================================
template <int ndim>
void BinaryTree<ndim>::UpdateAllSphGravForces
(Sph<ndim> *sph)                    ///< Pointer to SPH object
{
  return;
}



//=============================================================================
//  BinaryTree::UpdateAllSphDerivatives
/// ..
//=============================================================================
template <int ndim>
void BinaryTree<ndim>::UpdateAllSphDerivatives(Sph<ndim> *sph)
{
  return;
}



//=============================================================================
//  BinaryTree::UpdateAllSphDudt
/// Compute all local 'gather' properties of currently active particles, and 
/// then compute each particle's contribution to its (active) neighbour 
/// neighbour hydro forces.  Optimises the algorithm by using grid-cells to 
/// construct local neighbour lists for all particles  inside the cell.
//=============================================================================
template <int ndim>
void BinaryTree<ndim>::UpdateAllSphDudt(Sph<ndim> *sph)
{
  return;
}



#if defined(VERIFY_ALL)
//=============================================================================
//  BinaryTree::ValidateTree
/// ..
//=============================================================================
template <int ndim>
void BinaryTree<ndim>::ValidateTree
(Sph<ndim> *sph)
{
  bool overlap_flag = false;
  bool hmax_flag = false;
  int activecount;
  int c;
  int cc;
  int i;
  int j;
  int l;
  int leafcount;
  int Nactivecount;
  int Ncount;
  int *ccount;
  int *pcount;
  BinaryTreeCell<ndim> cell;

  debug2("[BinaryTree::ValidateTree]");

  ccount = new int[Ncellmax];
  pcount = new int[Ntot];
  for (i=0; i<Ntot; i++) pcount[i] = 0;
  for (c=0; c<Ncellmax; c++) ccount[c] = 0;
  Ncount = 0;
  Nactivecount = 0;

  // Count how many times we enter a cell in a full tree walk
  c = 0;
  while (c < Ncell) {
    ccount[c]++;
    if (tree[c].c1 != -1) c = tree[c].c1;
    else c = tree[c].cnext;
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
  for (i=0; i<sph->Ntot; i++) {
    if (!(ids[i] >= 0 && ids[i] < sph->Ntot)) {
      cout << "Problem with ids array : " 
	   << i << "   " << ids[i] << endl;
      exit(0);
    }
    if (!(inext[i] >= -1 && inext[i] < sph->Ntot)) {
      cout << "Problem with inext linked lists : " 
	   << i << "   " << inext[i] << endl;
      exit(0);
    }
  }


  //---------------------------------------------------------------------------
  for (c=0; c<Ncell; c++) {
    cell = tree[c];
    activecount = 0;
    leafcount = 0;

    // Check that particles are not in linked lists more than once
    if (cell.level == ltot) {
      i = cell.ifirst;
      while (i != -1) {
	pcount[i]++;
	leafcount++;
	Ncount++;
	if (sph->sphdata[i].active) activecount++;
	if (sph->sphdata[i].active) Nactivecount++;
        if (sph->sphdata[i].h > cell.hmax) hmax_flag = true;
        if (i == cell.ilast) break;
	i = inext[i];
      }
      if (hmax_flag) {
	cout << "hmax flag error : " << c << endl;
	exit(0);
      }
      if (leafcount > Nleafmax) {
	cout << "Leaf particle count error : " 
	     << leafcount << "   " << Nleafmax << endl;
	exit(0);
      }
      if (activecount > leafcount) {
	cout << "Leaf particle count error : " 
	     << leafcount << "   " << Nleafmax << endl;
	exit(0);
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


  if (Ncount != sph->Ntot) {
    cout << "Ncount problem with binary tree : " 
	 << Ncount << "   " << sph->Ntot << endl;
    exit(0);
  }
  if (Nactivecount > sph->Nsph) {
    cout << "Nactivecount problem with binary tree : " 
	 << Nactivecount << "   " << sph->Nsph << "   " << sph->Ntot << endl;
    exit(0);
  }


  cout << "Tree all fine (apparently)" << endl;

  delete[] pcount;
  delete[] ccount;

  return;
}



//=============================================================================
//  BinaryTree::CheckValidNeighbourList
/// Checks that the neighbour list generated by the grid is valid in that it 
/// (i) does include all true neighbours, and 
/// (ii) all true neigbours are only included once and once only.
//=============================================================================
template <int ndim>
void BinaryTree<ndim>::CheckValidNeighbourList
(Sph<ndim> *sph,                    ///< [in] SPH object pointer
 int i,                             ///< [in] Particle i.d.
 int Nneib,                         ///< [in] No. of potential neighbours
 int *neiblist,                     ///< [in] List of potential neighbour i.d.s
 string neibtype)                   ///< [in] Neighbour search type
{
  int c;                            // ..
  int count;                        // Valid neighbour counter
  int ii;                           // ..
  int j;                            // Neighbour particle counter
  int k;                            // Dimension counter
  int Ntrueneib = 0;                // No. of 'true' neighbours
  int Nleaflist = 0;                // ..
  int *trueneiblist;                // List of true neighbour ids
  int *leaflist;                    // ..
  FLOAT drsqd;                      // Distance squared
  FLOAT dr[ndim];                   // Relative position vector
  return;
  // Allocate array to store local copy of potential neighbour ids
  trueneiblist = new int[sph->Ntot];
  leaflist = new int[sph->Ntot];

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

      // Compute list from leaf cells
      for (c=0; c<Ncell; c++) {
	if (tree[c].level == ltot && tree[c].N > 0) {
	  ii = tree[c].ifirst;
	  while (ii != -1) {
	    for (k=0; k<ndim; k++) 
	      dr[k] = sph->sphdata[ii].r[k] - sph->sphdata[i].r[k];
	    drsqd = DotProduct(dr,dr,ndim);
	    if (neibtype == "all" && 
		(drsqd < sph->kernp->kernrangesqd*sph->sphdata[i].h
		 *sph->sphdata[i].h || drsqd < sph->kernp->kernrangesqd*
		 sph->sphdata[ii].h*sph->sphdata[ii].h))
	      leaflist[Nleaflist++] = ii;
	    else if (neibtype == "gather" && 
		     drsqd < sph->kernp->kernrangesqd*sph->sphdata[i].h
		     *sph->sphdata[i].h) {
	      leaflist[Nleaflist++] = ii;
	      cout << "Found ptcl " << ii << "  on leaf : " << c << endl;
	    }
	    if (ii == tree[c].ilast) break;
	    ii = inext[ii];
	  };
	}
      }

      cout << "Problem with neighbour lists : " << neibtype << "   " 
	   << i << "   " << jj << "   " << count << "   "
	   << sph->sphdata[i].r[0] << "   " << sph->sphdata[i].h << endl;
      cout << "Nsph  : " << Nsph << "    Ntot : " << Ntot << "    Ntotmax : " << Ntotmax << endl;
      cout << "True Nsph  : " << sph->Nsph << "    True Ntot : " << sph->Ntot << endl;
      cout << "Nneib : " << Nneib << "   Ntrueneib : " << Ntrueneib << "    Nleaflist : " << Nleaflist << endl;
   	  for (k=0; k<ndim; k++) dr[k] = sph->sphdata[jj].r[k] - sph->sphdata[i].r[k];
      drsqd = DotProduct(dr,dr,ndim);
      cout << "drmag : " << sqrt(drsqd) << "     "
	   << sph->kernp->kernrange*sph->sphdata[i].h
	   << "    " << sph->kernp->kernrange*sph->sphdata[jj].h 
	   << "     " << sph->sphdata[jj].h << endl;
      PrintArray("neiblist     : ",Nneib,neiblist);
      PrintArray("trueneiblist : ",Ntrueneib,trueneiblist);
      PrintArray("leaflist     : ",Nleaflist,leaflist);
      string message = "Problem with neighbour lists in Binary tree";
      ExceptionHandler::getIstance().raise(message);
    }
  }

  delete[] trueneiblist;
  delete[] leaflist;

  //cout << "List okay!! : " << Nneib << "   " << Ntrueneib << endl;

  return;
}
#endif



template class BinaryTree<1>;
template class BinaryTree<2>;
template class BinaryTree<3>;
