//=============================================================================
//  MpiTree.cpp
//  ...
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
#include "MpiTree.h"
#include "Precision.h"
#include "Exception.h"
#include "Parameters.h"
#include "InlineFuncs.h"
#include "Particle.h"
#include "Debug.h"
#if defined _OPENMP
#include <omp.h>
#endif
using namespace std;



//=============================================================================
//  MpiTree::MpiTree()
/// Constructor for KD-tree radiation class
//=============================================================================
template <int ndim, template<int> class ParticleType>
MpiTree<ndim,ParticleType>::MpiTree(int Nmpiaux)
{
  allocated_tree = false;
  ltot     = 0;
  Ntot     = 0;
  Ntotmax  = 0;
  Nmpi     = Nmpiaux;
#if defined _OPENMP
  Nthreads = omp_get_max_threads();
#else
  Nthreads = 1;
#endif
}



//=============================================================================
//  MpiTree::~MpiTree()
/// Destructor for KD-tree radiation class
//=============================================================================
template <int ndim, template<int> class ParticleType>
MpiTree<ndim,ParticleType>::~MpiTree()
{
}



//=============================================================================
//  MpiTree::AllocateTreeMemory
/// Allocate memory for KD-tree as requested.  If more memory is required
/// than currently allocated, tree is deallocated and reallocated here.
//=============================================================================
template <int ndim, template<int> class ParticleType>
void MpiTree<ndim,ParticleType>::AllocateMemory(void)
{
  debug2("[MpiTree::AllocateMemory]");

  if (!allocated_tree) {
    if (allocated_tree) DeallocateMemory();
    Ntotmax = max(Ntotmax,Ntot);

    klevel = new int[lmax];
    g2c    = new int[gmax];
    ids    = new int[Ntotmax];
    inext  = new int[Ntotmax];
    tree   = new struct MpiTreeCell<ndim>[Ncellmax];

    allocated_tree = true;
  }

  return;
}



//=============================================================================
//  MpiTree::DeallocateMemory
/// Deallocates all KD-tree memory
//=============================================================================
template <int ndim, template<int> class ParticleType>
void MpiTree<ndim,ParticleType>::DeallocateMemory(void)
{
  debug2("[MpiTree::DeallocateMemory]");

  if (allocated_tree) {
    delete[] klevel;
    delete[] tree;
    delete[] inext;
    delete[] ids;
    delete[] g2c;
    allocated_tree = false;
  }

  return;
}



//=============================================================================
//  MpiTree::ComputeTreeSize
/// Compute the maximum size (i.e. no. of levels, cells and leaf cells) of
/// the MPI KD-tree.
//=============================================================================
template <int ndim, template<int> class ParticleType>
void MpiTree<ndim,ParticleType>::ComputeTreeSize(void)
{
  debug2("[MpiTree::ComputeTreeSize]");

  // Calculate level of tree that can contain required number of MPI
  // nodes on lowest level
  ltot = 0;
  while (pow(2,ltot) < Nmpi) {
    ltot++;
  };
  lmax = ltot;
  gtot = pow(2,ltot);
  gmax = pow(2,ltot);
  Ncell = 2*gtot - 1;
  Ncellmax = Ncell;

#ifdef VERIFY_ALL
  cout << "No. of levels on MPI tree : " << ltot << "   " << ltot << endl;
  cout << "No. of MPI cells in tree  : " << Nmpi << "   " << gtot << endl;
#endif

  return;
}



//=============================================================================
//  MpiTree::CreateTreeStructure
/// Create the raw tree skeleton structure once the tree size is known.
/// Sets all cell pointer variables and all cell levels.
//=============================================================================
template <int ndim, template<int> class ParticleType>
void MpiTree<ndim,ParticleType>::CreateTreeStructure
 (MpiNode<ndim> *mpinode)              ///< ..
{
  int c;                               // Dummy id of tree-level, then tree-cell
  int g;                               // Dummy id of grid-cell
  int l;                               // Dummy id of level
  int *c2L;                            // Increment to second child-cell
  int *cNL;                            // Increment to next cell if cell unopened

  debug2("[MpiTree::CreateTreeStructure]");

  // Allocate memory for local arrays
  c2L = new int[ltot + 1];
  cNL = new int[ltot + 1];

  // Set pointers to second child-cell (if opened) and next cell (if unopened)
  for (l=0; l<ltot; l++) {
    c2L[l] = pow(2,ltot - l);
    cNL[l] = 2*c2L[l] - 1;
  }

  // Zero tree cell variables
  for (c=0; c<Ncell; c++) {
    tree[c].id     = c;
    tree[c].c2g    = 0;
    tree[c].c1     = -1;
    tree[c].c2     = -1;
    tree[c].ifirst = -1;
    tree[c].ilast  = -1;
    tree[c].N      = 0;
    tree[c].nodes.clear();
  }
  g = 0;
  tree[0].level = 0;

  // Loop over all cells and set all other pointers
  //---------------------------------------------------------------------------
  for (c=0; c<Ncell; c++) {
    if (tree[c].level == ltot) {                    // If on leaf level
      tree[c].cnext = c + 1;                        // id of next cell
      tree[c].c2g = g;                              // Record leaf id
      g2c[g++] = c;                                 // Record inverse id
    }
    else {
      tree[c+1].level = tree[c].level + 1;          // Level of 1st child
      tree[c].c1 = c + 1;                           // ..
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
//  MpiTree::DivideTreeCell
/// Recursive routine to divide a tree cell into two children cells.
//=============================================================================
template <int ndim, template<int> class ParticleType>
void MpiTree<ndim,ParticleType>::DivideTreeCell
 (int ifirst,                          ///< [in] Aux. id of first particle in cell
  int ilast,                           ///< [in] Aux. id of last particle in cell
  ParticleType<ndim> *partdata,        ///< [in] Pointer to main SPH object
  MpiTreeCell<ndim> &cell)             ///< [inout] Cell to be divided
{
  int i;                               // Aux. child cell counter
  int j;                               // Aux. particle counter
  int k;                               // Dimension counter
  int k_divide = 0;                    // Division dimension
  FLOAT rkmax = 0.0;                   // Max. box size of all dimensions
  FLOAT r_divide;                      // Coordinate value at division


  // If cell is a leaf cell, do not divide further and set linked lists
  if (cell.level == ltot) {
    if (cell.N > 0) {
      for (j=cell.ifirst; j<cell.ilast; j++) inext[ids[j]] = ids[j+1];
      inext[ids[cell.ilast]] = -1;
      cell.ifirst = ids[cell.ifirst];
      cell.ilast = ids[cell.ilast];
    }
    else {
      cell.ifirst = -1;
      cell.ilast = -1;
    }
    cell.nodes.clear();
    cell.nodes.push_back(cell.c2g);
    return;
  }

  i = cell.ifirst;

  // Determine dimension to split the cell along.
  // For now, simply split along direction of the bounding box's longest axis
  for (k=0; k<ndim; k++) {
    if (cell.box.max[k] - cell.box.min[k] > rkmax) {
      rkmax = cell.box.max[k] - cell.box.min[k];
      k_divide = k;
    }
  }
  cell.k_divide = k_divide;


  // Find median value along selected division plane and re-order array
  // so particles reside on correct side of division
  r_divide = QuickSelect(cell.ifirst,cell.ilast,cell.ifirst+cell.N/2,k_divide,partdata);
  cell.r_divide = r_divide;

  // Set properties of first child cell
  for (k=0; k<ndim; k++) tree[cell.c1].box.min[k] = cell.box.min[k];
  for (k=0; k<ndim; k++) tree[cell.c1].box.max[k] = cell.box.max[k];
  tree[cell.c1].N = cell.N/2;
  if (tree[cell.c1].N != 0) {
    tree[cell.c1].ifirst = ifirst;
    tree[cell.c1].ilast = ifirst + cell.N/2 - 1;
  }

  // Set properties of second child cell
  for (k=0; k<ndim; k++) tree[cell.c2].box.min[k] = cell.box.min[k];
  for (k=0; k<ndim; k++) tree[cell.c2].box.max[k] = cell.box.max[k];
  tree[cell.c2].N = cell.N - tree[cell.c1].N;
  if (tree[cell.c2].N != 0) {
    tree[cell.c2].ifirst = ifirst + cell.N/2;
    tree[cell.c2].ilast = ilast;
  }
  assert(cell.N == tree[cell.c1].N + tree[cell.c2].N);


  // Set new cell boundaries depending on number of particles in cells
  if (tree[cell.c1].N > 0 && tree[cell.c2].N > 0) {
    tree[cell.c1].box.max[k_divide] = r_divide;
    tree[cell.c2].box.min[k_divide] = r_divide;
  }
  else if (tree[cell.c2].N > 0) {
    tree[cell.c1].box.max[k_divide] = -big_number;
  }


  // Now divide the new child cells as a recursive function
#if defined _OPENMP
  if (pow(2,cell.level) < Nthreads) {
#pragma omp parallel default(none) private(i) \
  shared(cell,ifirst,ilast,partdata) num_threads(2)
    {
#pragma omp for
      for (i=0; i<2; i++) {
        if (i == 0) DivideTreeCell(ifirst,ifirst+cell.N/2-1,partdata,tree[cell.c1]);
        else if (i == 1) DivideTreeCell(ifirst+cell.N/2,ilast,partdata,tree[cell.c2]);
      }
#pragma omp barrier
    }
  }
  else {
    for (i=0; i<2; i++) {
      if (i == 0) DivideTreeCell(ifirst,ifirst+cell.N/2-1,partdata,tree[cell.c1]);
      else if (i == 1) DivideTreeCell(ifirst+cell.N/2,ilast,partdata,tree[cell.c2]);
    }
  }
#else
  for (i=0; i<2; i++) {
    if (i == 0) DivideTreeCell(ifirst,ifirst+cell.N/2-1,partdata,tree[cell.c1]);
    else if (i == 1) DivideTreeCell(ifirst+cell.N/2,ilast,partdata,tree[cell.c2]);
  }
#endif

  // Join node lists of two cells together
  cell.nodes = tree[cell.c1].nodes;
  cell.nodes.insert(cell.nodes.end(),tree[cell.c2].nodes.begin(),tree[cell.c2].nodes.end());

  if (cell.N != tree[cell.c1].N + tree[cell.c2].N) {
    cout << "Checking : " << cell.N << "   " << tree[cell.c1].N
         << "    " << tree[cell.c2].N << endl;
  }
  assert(cell.N == tree[cell.c1].N + tree[cell.c2].N);

  return;
}



//=============================================================================
//  MpiTree::QuickSelect
/// Find median and sort particles in arrays to ensure they are the correct
/// side of the division.  Uses the QuickSelect algorithm.
//=============================================================================
template <int ndim, template<int> class ParticleType>
FLOAT MpiTree<ndim,ParticleType>::QuickSelect
(int left,                          ///< Left-most id of particle in array
 int right,                         ///< Right-most id of particle in array
 int jpivot,                        ///< Pivot/median point
 int k,                             ///< Dimension of sort
 ParticleType<ndim> *partdata)      ///< Pointer to main SPH object
{
  int j;                            // Aux.
  int jguess;                       // ..
  int jtemp;                        // ..
  FLOAT rpivot;                     // Position pivot for quick-select


  // Place all particles left or right of chosen pivot point.
  // Iterate until correct median pivot has been identified.
  //---------------------------------------------------------------------------
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


    //-------------------------------------------------------------------------
    for (j=left; j<right; j++) {
      assert(j < Ntot);
      if (partdata[ids[j]].r[k] < rpivot) {
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
//  MpiTree::FindCell
/// ...
//=============================================================================
template <int ndim, template<int> class ParticleType>
int MpiTree<ndim,ParticleType>::FindCell
(int cparent,                       ///< [in] i.d. of larger parent cell
 FLOAT rp[ndim])                    ///< [in] Position of point/ray
{
  int c = cparent;                  // Cell i.d.
  int c1;                           // i.d. of 1st cell child
  int k_divide;                     // Dimension of cell division

  // Walk back down through tree to bottom level
  //---------------------------------------------------------------------------
  while (tree[c].level < ltot) {

#ifdef OUTPUT_ALL
    cout << "Searching for cell : " << tree[c].level << "   " << ltot << endl;
#endif

    c1 = c + 1;
    k_divide = tree[c].k_divide;

    // If point is left of divide, pick 1st child cell.  Else pick 2nd child.
    if (rp[k_divide] < tree[c1].box.max[k_divide])
      c = c1;
    else
      c = tree[c].c2;

  };
  //---------------------------------------------------------------------------

  return c;
}



//=================================================================================================
//  MpiTree::UpdateBoundingBoxes
/// Update the bounding boxes of all nodes on lower levels of the tree based on adjustment of
/// higher level division (such as for load-balancing).
//=================================================================================================
template <int ndim, template<int> class ParticleType>
void MpiTree<ndim,ParticleType>::UpdateBoundingBoxes(void)
{
  int c;                               // Cell counter
  int c1;                              // 1st child cell id
  int c2;                              // 2nd child cell id
  int k;                               // Dimension counter
  int k_divide;                        // Dimension of cell division (between children)
  FLOAT r_divide;                      // Position of division

  //-----------------------------------------------------------------------------------------------
  for (c=0; c<Ncell; c++) {
    if (tree[c].level == ltot) continue;

    c1 = c + 1;
    c2 = tree[c].c2;
    k_divide = tree[c].k_divide;
    r_divide = tree[c].r_divide;

    // Adjust bounding box of 1st child
    for (k=0; k<ndim; k++) tree[c1].box.min[k] = tree[c].box.min[k];
    for (k=0; k<ndim; k++) tree[c1].box.max[k] = tree[c].box.max[k];
    tree[c1].box.max[k_divide] = r_divide;

    // Adjust bounding box of 2nd child
    for (k=0; k<ndim; k++) tree[c2].box.min[k] = tree[c].box.min[k];
    for (k=0; k<ndim; k++) tree[c2].box.max[k] = tree[c].box.max[k];
    tree[c2].box.min[k_divide] = r_divide;

  }
  //-----------------------------------------------------------------------------------------------

  return;
}




template class MpiTree<1,GradhSphParticle>;
template class MpiTree<2,GradhSphParticle>;
template class MpiTree<3,GradhSphParticle>;
template class MpiTree<1,SM2012SphParticle>;
template class MpiTree<2,SM2012SphParticle>;
template class MpiTree<3,SM2012SphParticle>;
template class MpiTree<1,MeshlessFVParticle>;
template class MpiTree<2,MeshlessFVParticle>;
template class MpiTree<3,MeshlessFVParticle>;
