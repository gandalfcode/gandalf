//=============================================================================
//  BinaryTree.cpp
//  Contains all functions for building, stocking and walking for the 
//  binary tree for SPH particles.
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
                             string multipole_aux, int Nthreads, int Nmpi)
{
  Nlocalsubtrees = Nthreads;
  Nmpisubtrees = max(Nmpi - 1,0);
  Nsubtreemax = Nlocalsubtrees + Nmpisubtrees;
  Nsubtree = Nsubtreemax;
  Ntot = 0;
  Ntotmax = 0;
  Ntotmaxold = 0;
  Nleafmax = Nleafmaxaux;
  kernrange = kernrangeaux;
  thetamaxsqd = thetamaxsqdaux;
  gravity_mac = gravity_mac_aux;
  multipole = multipole_aux;
  allocated_tree = false;
  created_sub_trees = false;
#if defined _OPENMP
  // Check that no. of threads is valid
  int ltot = 0;
  while (pow(2,ltot) < Nthreads) {
    ltot++;
  };
  if (Nthreads != pow(2,ltot)) {
    Nlocalsubtrees = pow(2,ltot-1);
    Nsubtreemax = Nlocalsubtrees + Nmpisubtrees;
    cout << "Warning: the number of OpenMP threads is not a power of two. This is sub-optimal for the binary tree parallelization" << endl;
  }
#endif
  assert(Nlocalsubtrees > 0);
  assert(Nsubtreemax > 0);
  }



//=============================================================================
//  BinaryTree::~BinaryTree
/// BinaryTree destructor.  Deallocates tree memory upon object destruction.
//=============================================================================
template <int ndim>
BinaryTree<ndim>::~BinaryTree()
{
  debug2("[BinaryTree::DeallocateTreeMemory]");
  if (allocated_tree) DeallocateTreeMemory();
}



//=============================================================================
//  BinaryTree::AllocateTreeMemory
/// Allocate memory for binary tree as requested.  If more memory is required 
/// than currently allocated, tree is deallocated and reallocated here.
//=============================================================================
template <int ndim>
void BinaryTree<ndim>::AllocateTreeMemory(void)
{
  debug2("[BinaryTree::AllocateTreeMemory]");

  if (!allocated_tree || Ntotmax > Ntotmaxold || Nsubtree > Nsubtreemax) {
    if (allocated_tree) DeallocateTreeMemory();
    Ntotmax = max(Ntotmax,Ntot);
    Ntotmaxold = Ntotmax;
    Nsubtreemax = max(Nsubtreemax,Nsubtree);
    Nsubtreemaxold = Nsubtreemax;

    pc = new int[Ntotmax];
    pw = new FLOAT[Ntotmax];
    tree = new struct BinaryTreeCell<ndim>[Ncellmax];
    for (int k=0; k<ndim; k++) porder[k] = new int[Ntotmax]; 
    for (int k=0; k<ndim; k++) rk[k] = new FLOAT[Ntotmax];
    allocated_tree = true;

    if (!created_sub_trees) {
      for (int i=0; i<Nsubtree; i++) {
        subtrees.push_back(new BinarySubTree<ndim>(Nleafmax, thetamaxsqd,
						   kernrange, gravity_mac, multipole));
      }
      created_sub_trees = true;
    }


    // ..
    for (int i=0; i<Nsubtree; i++) {
      subtrees[i]->Nsph = 0;
      subtrees[i]->Ntot = 0;
      subtrees[i]->Ntotmax = max(subtrees[i]->Ntotmax,Ntotmax/Nsubtree + 1);

      // ..
      subtrees[i]->ComputeSubTreeSize();
      subtrees[i]->AllocateSubTreeMemory();
    }

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
  debug2("[BinaryTree::DeallocateTreeMemory]");

  if (allocated_tree) {
    for (int i=Nsubtree-1; i>=0; i--) subtrees[i]->DeallocateSubTreeMemory();
    for (int k=ndim-1; k>=0; k--) delete[] rk[k];
    for (int k=ndim-1; k>=0; k--) delete[] porder[k];
    delete[] stree;
    delete[] tree;
    delete[] pw;
    delete[] pc;
    allocated_tree = false;
  }

  return;
}



//=============================================================================
//  BinaryTree::BuildTree
/// Call all routines to build/re-build the binary tree.
//=============================================================================
template <int ndim>
void BinaryTree<ndim>::BuildTree
(Sph<ndim> *sph,                    ///< Pointer to main SPH object
 Parameters &simparams)             ///< Simulation parameters
{
  int i;
  int Ncheck = 0;
  int localgtot = 0;

  debug2("[BinaryTree::BuildTree]");

  // Set number of tree members to total number of SPH particles (inc. ghosts)
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

  // Find ordered list of particle positions ready for adding particles to tree
  OrderParticlesByCartCoord(sph->sphdata);

  // Now add particles to tree depending on Cartesian coordinates
  LoadParticlesToTree(sph->rsph);

  // Build and stock all local sub-trees
  //---------------------------------------------------------------------------
#pragma omp parallel for default(none) private(i) shared(sph) \
  reduction(+:localgtot,Ncheck)
  for (i = 0; i < Nsubtree; i++) {

    BinarySubTree<ndim>* subtree = subtrees[i];

    // Build individual sub-trees
    subtree->BuildSubTree(sph,simparams);

    // Calculate all cell quantities (e.g. COM, opening distance)
    subtree->StockCellProperties(sph->sphdata);

    // Calculate total number of leaf cells in trees
    localgtot += subtree->gtot;
    Ncheck += subtree->Ntot;

  }
  //---------------------------------------------------------------------------

  gtot = localgtot;

  assert(Ncheck == sph->Ntot);

  return;
}



//=============================================================================
//  BinaryTree::UpdateTree
/// Call all routines to build/re-build the binary tree.
//=============================================================================
template <int ndim>
void BinaryTree<ndim>::UpdateTree
(Sph<ndim> *sph,                    ///< Pointer to main SPH object
 Parameters &simparams)             ///< Simulation parameters
{
  int i;

  debug2("[BinaryTree::UpdateTree]");

  //---------------------------------------------------------------------------
#pragma omp parallel for default(none) private(i) shared(sph)
  for (i=0; i<Nsubtree; i++) {

    // Calculate all cell quantities (e.g. COM, opening distance)
    subtrees[i]->StockCellProperties(sph->sphdata);

  }
  //---------------------------------------------------------------------------

  return;
}



//=============================================================================
//  BinaryTree::UpdateActiveParticleCounters
/// Loop through all leaf cells in binary tree and update all active 
/// particle counters.
//=============================================================================
template <int ndim>
void BinaryTree<ndim>::UpdateActiveParticleCounters(Sph<ndim> *sph)
{
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
  while (pow(2,ltot) < Nsubtree) {
    ltot++;
  };

  // Set total number of leaf/grid cells and tree cells
  gtot = pow(2,ltot);
  Ncell = 2*gtot - 1;
  Ncellmax = Ncell;

  // Optional output (for debugging)
#if defined(VERIFY_ALL)
  cout << "Calculating tree size variables" << endl;
  cout << "No. of sub-trees      : " << Nsubtree << endl;
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
  debug2("[BinaryTree::CreateTreeStructure]");

  int c;                            // Dummy id of tree-level, then tree-cell
  int g;                            // Dummy id of grid-cell
  int l;                            // Dummy id of level
  int *c2L;                         // Increment to second child-cell
  int *cNL;                         // Increment to next cell if cell unopened

  // Allocate memory for local arrays
  c2L = new int[ltot + 1];
  cNL = new int[ltot + 1];

  // Set pointers to second child-cell (if opened) and next cell (if unopened)
  for (l=0; l<ltot; l++) {
    c2L[l] = pow(2,ltot - l);
    cNL[l] = 2*c2L[l] - 1;
  }

  // Zero tree cell variables
  //for (g=0; g<gtot; g++) g2c[g] = 0;
  for (c=0; c<Ncell; c++) {
    tree[c].c2 = 0;
    tree[c].c2g = 0;
  }
  g = 0;
  tree[0].clevel = 0;

  // Loop over all cells and set all other pointers
  //---------------------------------------------------------------------------
  for (c=0; c<Ncell; c++) {
    if (tree[c].clevel == ltot) {                   // If on leaf level
      tree[c].cnext = c + 1;                        // id of next cell
      tree[c].c2g = g;                              // Record leaf id
      g++; //g2c[g++] = c;                          // Record inverse id
    }
    else {
      tree[c+1].clevel = tree[c].clevel + 1;        // Level of 1st child
      tree[c].c2 = c + c2L[tree[c].clevel];         // id of 2nd child
      tree[tree[c].c2].clevel = tree[c].clevel + 1; // Level of 2nd child
      tree[c].cnext = c + cNL[tree[c].clevel];      // Next cell id
    }

  }
  //---------------------------------------------------------------------------


  // Free locally allocated memory
  delete[] cNL;
  delete[] c2L;

  return;
}



//=============================================================================
//  BinaryTree::OrderParticlesByCartCoord
/// Compute lists of particle ids in order of x, y and z positions.
/// Used for efficiently creating the binary tree data structure.
//=============================================================================
template <int ndim>
void BinaryTree<ndim>::OrderParticlesByCartCoord
(SphParticle<ndim> *sphdata)        ///< SPH particle data array
{
  int i;                            // Particle counter
  int k;                            // Dimension counter

  debug2("[BinaryTree::OrderParticlesByCartCoord]");

  // First copy all values to local arrays
  for (k=0; k<ndim; k++)
    for (i=0; i<Ntot; i++) {
      rk[k][i] = sphdata[i].r[k];
    }

  // Now copy list of particle ids
  for (k=0; k<ndim; k++)
    for (i=0; i<Ntot; i++)
      porder[k][i] = i;

  // Divide sorting routines amongst (up to 3) threads
  //---------------------------------------------------------------------------
#pragma omp sections
  {
    // Sort x-values
#pragma omp section
    Heapsort(Ntot,porder[0],rk[0]);

    // Sort y-values
#pragma omp section
    if (ndim >= 2) Heapsort(Ntot,porder[1],rk[1]);

    // Sort z-values
#pragma omp section
    if (ndim == 3) Heapsort(Ntot,porder[2],rk[2]);
  }

  // Check that particles are ordered correctly
#if defined(VERIFY_ALL)
  for (k=0; k<ndim; k++) {
    for (int j=1; j<Ntot; j++) {
      int i1 = porder[k][j-1];
      int i2 = porder[k][j];
      assert(i1 >= 0 && i1 < Ntot);
      assert(i2 >= 0 && i2 < Ntot);
      if (sphdata[i2].r[k] < sphdata[i1].r[k]) {
        cout << "Problem with particle ordering : " << k << "   " << j 
	     << "    " << i1 << "    " << i2 << "    " 
	     << sphdata[i1].r[k] << "    " << sphdata[i2].r[k] << endl;
        exit(0);
      }
    }
  }
#endif

  return;
}



//=============================================================================
//  BinaryTree::LoadParticlesToTree
/// Create tree structure by adding particles to leaf cells.
//=============================================================================
template <int ndim>
void BinaryTree<ndim>::LoadParticlesToTree
(FLOAT *r)                         ///< Positions of particles
{
  int c;                           // Cell counter
  int cc;                          // Secondary cell counter
  int g;                           // ..
  int k;                           // Dimensionality counter
  int kk;                          // ..
  int i;                           // Particle counter
  int j;                           // Dummy particle id
  int l;                           // Level counter
  FLOAT *ccap;                     // Maximum capacity of cell
  FLOAT *ccon;                     // Current contents of cell

  debug2("[BinaryTree::LoadParticleToTree]");

  // Allocate memory for local arrays
  ccap = new FLOAT[Ncellmax];
  ccon = new FLOAT[Ncellmax];

  // Set capacity of root-cell using particle weights
  for (i=0; i<Ntot; i++) pw[i] = 1.0/(FLOAT) Ntot;
  for (c=0; c<Ncell; c++) ccap[c] = 0.0;
  for (i=0; i<Ntot; i++) ccap[0] += pw[i];

  // Initialise all particle and cell values before building tree structure
  for (i=0; i<Ntot; i++) pc[i] = 0;
  for (c=0; c<Ncell; c++) ccon[c] = 0.0;
  for (c=0; c<Ncell; c++) tree[c].ifirst = -1;
  for (c=0; c<Ncell; c++) tree[c].ilast = -1;

  // Set bounding box of root cell to main bounding box
  for (k=0; k<ndim; k++) tree[0].bbmin[k] = box->boxmin[k];
  for (k=0; k<ndim; k++) tree[0].bbmax[k] = box->boxmax[k];

  // Zero counters for each subtree
  for (i=0; i<Nsubtree; i++) subtrees[i]->Ntot = 0;
  for (i=0; i<Nsubtree; i++) subtrees[i]->Nsph = 0;

  // Start at top level (l = 0) dividing the cell along the x-axis (k = 0)
  l = 0;
  k = 0;


  // Loop through each level of the tree
  //---------------------------------------------------------------------------
  while (l < ltot) {

	//cout << "LEVEL : " << l << "    " << ltot << "    " << Ntot << endl;

    // Loop over all particles (in order of current split)
    //-------------------------------------------------------------------------
    for (i=0; i<Ntot; i++) {
      j = porder[k][i];
      cc = pc[j];                            // Cell currently occupied by j
      ccon[cc] += pw[j];                     // Add particle weighting to cell

      // If cell contains less than maximum allowed capacity, then add 
      // particle to the first child cell.  Otherwise add to second child.
      if (ccon[cc] < 0.5000000000001*ccap[cc]) {
        pc[j]++;
        ccap[pc[j]] += pw[j];
        if (tree[pc[j]].ifirst == -1) tree[pc[j]].ifirst = j;
        tree[pc[j]].ilast = j;
      }
      else {
        pc[j] = tree[cc].c2;
        ccap[tree[cc].c2] += pw[j];
        if (tree[tree[cc].c2].ifirst == -1) {
          tree[tree[cc].c2].ifirst = j;
          for (kk=0; kk<ndim; kk++) tree[cc+1].bbmin[kk] = tree[cc].bbmin[kk];
          for (kk=0; kk<ndim; kk++) tree[cc+1].bbmax[kk] = tree[cc].bbmax[kk];
          for (kk=0; kk<ndim; kk++) tree[tree[cc].c2].bbmin[kk] = tree[cc].bbmin[kk];
          for (kk=0; kk<ndim; kk++) tree[tree[cc].c2].bbmax[kk] = tree[cc].bbmax[kk];
          tree[cc+1].bbmax[k] = 0.5*(r[ndim*tree[cc+1].ifirst + k] + r[ndim*j + k]);
          tree[tree[cc].c2].bbmin[k] = 0.5*(r[ndim*tree[cc+1].ifirst + k] + r[ndim*j + k]);
          //cout << "FOUND SPLIT : " << l << "    " << cc << "     " << k << "     "
        	//	  << tree[cc].bbmin[k] << "    " << tree[cc+1].bbmax[k] << "     " << tree[cc].bbmax[k] << endl;
        }
        tree[pc[j]].ilast = j;
      }
    }
    //-------------------------------------------------------------------------

    // Move to next level and cycle through each dimension in turn
    // (Need more sophisticated algorithm here in future)
    l++;
    k = (k + 1)%ndim;

  }
  //---------------------------------------------------------------------------


  // Compute capacities of leaf cells here
  for (i=0; i<Ntot; i++) {
    cc = pc[i];
    ccon[cc] += pw[i];
  }


  // Loop over all particles and add the particle to the list of ids
  for (i=0; i<Ntot; i++) {
    c = pc[i];
    g = tree[c].c2g;
    assert(g >= 0 && g < Nsubtree);
    subtrees[g]->ids[subtrees[g]->Ntot++] = i;
    subtrees[g]->Nsph++;
  }


  // Set bounding boxes for sub-trees
  for (c=0; c<Ncell; c++) {
	//cout << "CELL BB : " << tree[c].bbmin[0] << "   " << tree[c].bbmax[0] << "     " << tree[c].bbmin[1] << "     " << tree[c].bbmax[1] << endl;

	if (tree[c].c2 == 0) {
      g = tree[c].c2g;
      for (k=0; k<ndim; k++) subtrees[g]->box.boxmin[k] = tree[c].bbmin[k];
      for (k=0; k<ndim; k++) subtrees[g]->box.boxmax[k] = tree[c].bbmax[k];
	}

  }


  // Free all locally allocated memory
  delete[] ccon;
  delete[] ccap;

  return;
}



//=============================================================================
//  BinaryTree::UpdateHmaxValues
/// Calculate the physical properties (e.g. total mass, centre-of-mass, 
/// opening-distance, etc..) of all cells in the tree.
//=============================================================================
template <int ndim>
void BinaryTree<ndim>::UpdateHmaxValues
(SphParticle<ndim> *sphdata)        ///< SPH particle data array
{
  int i;
  FLOAT hmax_aux=0;


#pragma omp parallel default(none) shared(sphdata,hmax_aux) private(i)
  {
    FLOAT hmax_local = 0;
#pragma omp for
    for (i = 0; i < Nsubtree; i++) {
      FLOAT subhmax = subtrees[i]->UpdateHmaxValues(sphdata);
      if (subhmax > hmax_local)
        hmax_local = subhmax;
    }
#pragma omp critical (hmax_critical)
    {
      if (hmax_local>hmax_aux) hmax_aux=hmax_local;
    }
  }

  hmax = hmax_aux;
  return;
}



//=============================================================================
//  BinaryTree::ComputeActiveCellList
/// Returns the number of cells containing active particles, 'Nactive', and
/// the i.d. list of cells contains active particles, 'celllist'
//=============================================================================
template <int ndim>
int BinaryTree<ndim>::ComputeActiveCellList
(BinaryTreeCell<ndim> **celllist,   ///< Cells id array containing active ptcls
 BinarySubTree<ndim> **treelist)    ///< 
{
  binlistiterator it;               // ..
  int Nactive = 0;                  // No. of cells containing active ptcls
  int Nfirst = 0;                   // ..

  debug2("[BinaryTree::ComputeActiveCellList]");

  // Iterate/loop over all sub-trees to find all cells containing active ptcls
  for (it = subtrees.begin(); it != subtrees.end(); it++) {
    Nactive = (*it)->ComputeActiveCellList(Nactive,celllist);
    for (int i=Nfirst; i<Nactive; i++) treelist[i] = (*it);
    Nfirst = Nactive;
  }

  return Nactive;
}



//=============================================================================
//  BinaryTree::ComputeActiveParticleList
/// Returns the number (Nactive) and list of ids (activelist) of all active
/// SPH particles in the given cell.
//=============================================================================
template <int ndim>
int BinaryTree<ndim>::ComputeActiveParticleList
(BinaryTreeCell<ndim> *cell,        ///< [in] Pointer to cell
 BinarySubTree<ndim> *treeptr,      ///< [in] List of ptrs to sub-trees
 Sph<ndim> *sph,                    ///< [in] SPH object pointer
 int *activelist)                   ///< [out] List of active particles in cell
{
  int Nactive = 0; // No. of active particles in cell
  int i = cell->ifirst; // Particle id (set to first ptcl id)
  int ilast = cell->ilast; // i.d. of last particle in cell c
  int j;

  // Walk through linked list to obtain list and number of active ptcls.
  while (i != -1) {
    j = treeptr->GlobalId(i);
    if (j < sph->Nsph && sph->sphdata[j].active)
      activelist[Nactive++] = j;
    if (i == ilast) break;
    i = treeptr->inext[i];
  };

  return Nactive;
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
(BinaryTreeCell<ndim> *cell,        ///< [in] Pointer
 int Nneibmax,                      ///< [in] Max. no. of neighbours
 int *neiblist,                     ///< [out] List of neighbour i.d.s
 FLOAT hmax,                        ///< [in] Maximum smoothing length
 SphParticle<ndim> *sphdata)        ///< [in] SPH particle data
{
  binlistiterator it;               // ..
  int Nneib = 0;                    // Total number of gather neighbours
                                    // (summed over all sub-trees)

  // Iterate/loop over all sub-trees to find all gather neighbours
  for (it = subtrees.begin(); it != subtrees.end(); it++) {
    Nneib = (*it)->ComputeGatherNeighbourList(cell,Nneib,Nneibmax,
                                              neiblist,hmax,sphdata);
    if (Nneib == -1) return Nneib;
  }

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
(BinaryTreeCell<ndim> *cell,        ///< [in] Pointer
 int Nneibmax,                      ///< [in] Max. no. of neighbours
 int *neiblist,                     ///< [out] List of neighbour i.d.s
 SphParticle<ndim> *sphdata)        ///< [in] SPH particle data
{
  binlistiterator it;               // ..
  int Nneib = 0;                    // Total number of gather neighbours
                                    // (summed over all sub-trees)

  for (it = subtrees.begin(); it != subtrees.end(); it++) {
    Nneib = (*it)->ComputeNeighbourList(cell,Nneib,Nneibmax,
                                        neiblist,sphdata);
    if (Nneib == -1) return Nneib;
  }

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
(BinaryTreeCell<ndim> *cell,        ///< [in] Pointer to cell
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
  binlistiterator it;               // ..
  int okflag;

  // Zero counters here (ComputeGravityInteractionList routine appends data)
  Nneib = 0;
  Ndirect = 0;
  Ngravcell = 0;

  // Iterate/loop over all sub-trees to find all SPH neighbours, direct-sum 
  // gravity and distant cell interaction lists
  for (it = subtrees.begin(); it != subtrees.end(); it++) {
    okflag = (*it)->ComputeGravityInteractionList(cell,Nneibmax,Ndirectmax,
						  Ngravcellmax,Nneib,Ndirect,
						  Ngravcell,neiblist,
						  directlist,gravcelllist,
						  sphdata);
    if (okflag == -1) return -1;
  }

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
    invdrmag = 1.0/sqrt(drsqd);
    invdr3 = invdrmag*invdrmag*invdrmag; //pow(invdrmag,3);

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
  BinarySubTree<ndim> **treelist;  // ..
  BinaryTreeCell<ndim> *cell;      // Pointer to binary tree cell
  BinaryTreeCell<ndim> **celllist; // List of binary cell pointers
  SphParticle<ndim> *data = sph->sphdata;  // Pointer to SPH particle data

  debug2("[BinaryTree::UpdateAllSphProperties]");

  // Find list of all cells that contain active particles
  celllist = new BinaryTreeCell<ndim>*[gtot];
  treelist = new BinarySubTree<ndim>*[gtot];
  cactive = ComputeActiveCellList(celllist,treelist);


  // Set-up all OMP threads
  //===========================================================================
#pragma omp parallel default(none) private(activelist,cc,cell,celldone,draux)\
  private(drsqd,drsqdaux,hmax,hrangesqd,i,j,jj,k,okflag,m,mu,Nactive,neiblist)\
  private(Nneib,Nneibmax,r,rp,gpot,gpot2,m2,mu2,Ngather)\
  shared(sph,celllist,cactive,data,nbody,treelist)
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
        Nactive = ComputeActiveParticleList(cell,treelist[cc],sph,activelist);

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

  delete[] treelist;
  delete[] celllist;

  // Update all tree smoothing length values
  UpdateHmaxValues(sph->sphdata);

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
  int *interactlist;               // ..
  int *neiblist;                   // List of neighbour ids
  FLOAT draux[ndim];               // Aux. relative position vector var
  FLOAT drsqd;                     // Distance squared
  FLOAT hrangesqdi;                // Kernel gather extent
  FLOAT rp[ndim];                  // Local copy of particle position
  FLOAT *dr;                       // Array of relative position vectors
  FLOAT *drmag;                    // Array of neighbour distances
  FLOAT *invdrmag;                 // Array of 1/drmag between particles
  BinarySubTree<ndim> **treelist;  // ..
  BinaryTreeCell<ndim> *cell;      // Pointer to binary tree cell
  BinaryTreeCell<ndim> **celllist; // List of binary tree pointers
  SphParticle<ndim> *activepart;   // Local copy of active particles
  SphParticle<ndim> *neibpart;     // Local copy of neighbouring ptcls
  SphParticle<ndim> *data = sph->sphdata;   // Pointer to SPH particle data

  debug2("[BinaryTree::UpdateAllSphHydroForces]");


  // Find list of all cells that contain active particles
  celllist = new BinaryTreeCell<ndim>*[gtot];
  treelist = new BinarySubTree<ndim>*[gtot];
  cactive = ComputeActiveCellList(celllist,treelist);


  // Set-up all OMP threads
  //===========================================================================
#pragma omp parallel default(none) private(activelist,activepart,cc,cell,dr)\
  private(draux,drmag,drsqd,hrangesqdi,i,interactlist,invdrmag,j,jj,k) \
  private(Nactive,neiblist,neibpart,Ninteract,Nneib,Nneibmax,rp)\
  shared(cactive,celllist,data,sph,treelist)
  {
    Nneibmax = 2*sph->Ngather;
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
      Nactive = ComputeActiveParticleList(cell,treelist[cc],sph,activelist);

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

  delete[] treelist;
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
  int j;                            // Aux. particle counter
  int jj;                           // Aux. particle counter
  int k;                            // Dimension counter
  int okflag;                       // Flag if h-rho iteration is valid
  int Nactive;                      // No. of active particles in cell
  int Ndirect;                      // ..
  int Ndirectmax;                   // ..
  int Ngravcell;                    // ..
  int Ngravcellmax;                 // ..
  int Ninteract;                    // ..
  int Nneib;                        // No. of neighbours
  int Nneibmax;                     // Max. no. of neighbours
  int *activelist;                  // List of active particle ids
  int *directlist;                  // ..
  int *interactlist;                // ..
  int *neiblist;                    // List of neighbour ids
  //FLOAT draux[ndim];                // Aux. relative position vector var
  //FLOAT drsqd;                      // Distance squared
  FLOAT *agrav;                     // Local copy of gravitational accel.
  FLOAT *gpot;                      // ..
  BinarySubTree<ndim> **treelist;   // ..
  BinaryTreeCell<ndim> *cell;       // Pointer to binary tree cell
  BinaryTreeCell<ndim> **celllist;  // ..
  BinaryTreeCell<ndim> **gravcelllist; // ..
  SphParticle<ndim> *neibpart;      // Local copy of neighbouring ptcls
  SphParticle<ndim> *activepart;    // Local copy of SPH particle
  SphParticle<ndim> *data = sph->sphdata;   // Pointer to SPH particle data

  debug2("[BinaryTree::UpdateAllSphForces]");


  // Find list of all cells that contain active particles
  celllist = new BinaryTreeCell<ndim>*[gtot];
  treelist = new BinarySubTree<ndim>*[gtot];
  cactive = ComputeActiveCellList(celllist,treelist);


  // Set-up all OMP threads
  //===========================================================================
#pragma omp parallel default(none) private(activelist,agrav,cc,cell)\
  private(gpot,i,interactlist,j,jj,activepart)\
  private(k,okflag,Nactive,neiblist,neibpart,Ninteract,Nneib,directlist)\
  private(gravcelllist,Ngravcell,Ndirect,Nneibmax,Ndirectmax,Ngravcellmax) \
  shared(celllist,cactive,sph,data,treelist,cout)
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
    gravcelllist = new BinaryTreeCell<ndim>*[Ngravcellmax];
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
      Nactive = ComputeActiveParticleList(cell,treelist[cc],sph,activelist);

      // Make local copies of active particles
      for (j=0; j<Nactive; j++) {
        assert(activelist[j] >= 0 && activelist[j] < sph->Nsph);
        activepart[j] = data[activelist[j]];
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
        gravcelllist = new BinaryTreeCell<ndim>*[Ngravcellmax];
        neibpart = new SphParticle<ndim>[Nneibmax];
        okflag = ComputeGravityInteractionList(cell,Nneibmax,Ndirectmax,
                                               Ngravcellmax,Nneib,Ndirect,
                                               Ngravcell,neiblist,directlist,
                                               gravcelllist,sph->sphdata);
      };


      // Make local copies of all potential neighbours
      for (j=0; j<Nneib; j++) {
        assert(neiblist[j] >= 0 && neiblist[j] < sph->Ntot);
        neibpart[j] = data[neiblist[j]];
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

        // Determine interaction list (to ensure we don't compute pair-wise
        // forces twice)
        Ninteract = 0;
        for (jj=0; jj<Nneib; jj++) {
          if ((neiblist[jj] < i && !neibpart[jj].active) || neiblist[jj] > i) {
            interactlist[Ninteract++] = jj;
          }
        }

        // Compute forces between SPH neighbours (hydro and gravity)
        sph->ComputeSphHydroGravForces(i,Ninteract,interactlist,
                                                          activepart[j],neibpart);

        // Compute direct gravity forces between distant particles
        sph->ComputeDirectGravForces(i,Ndirect,directlist,
                                                      agrav,gpot,activepart[j],data);

        // Compute gravitational force dues to distant cells
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
          data[i].a[k] += activepart[j].a[k];
          data[i].agrav[k] += activepart[j].agrav[k];
        }
        data[i].gpot += activepart[j].gpot;
        data[i].dudt += activepart[j].dudt;
        data[i].div_v += activepart[j].div_v;
        data[i].levelneib = max(data[i].levelneib,activepart[j].levelneib);
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
            data[j].a[k] += neibpart[jj].a[k];
            data[j].agrav[k] += neibpart[jj].agrav[k];
          }
          data[j].gpot += neibpart[jj].gpot;
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


    // Finally, add all contributions from distant pair-wise forces to arrays
    for (i=0; i<sph->Nsph; i++) {
      if (data[i].active) {
        for (k=0; k<ndim; k++) {
#pragma omp atomic
          data[i].agrav[k] += agrav[ndim*i + k];
        }
#pragma omp atomic
        data[i].gpot += gpot[i];
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

  delete[] treelist;
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
      cout << "Problem with neighbour lists : " << i << "  " << j << "   "
	   << count << "   "
	   << sph->sphdata[i].r[0] << "   " << sph->sphdata[i].h << endl;
      cout << "Nneib : " << Nneib << "   Ntrueneib : " << Ntrueneib << endl;
   	  for (k=0; k<ndim; k++) dr[k] = sph->sphdata[jj].r[k] - sph->sphdata[i].r[k];
      drsqd = DotProduct(dr,dr,ndim);
      cout << "drmag : " << sqrt(drsqd) << "     "
    		  << sph->kernp->kernrange*sph->sphdata[i].h
    		  << "    " << sph->kernp->kernrange*sph->sphdata[jj].h << endl;
      PrintArray("neiblist     : ",Nneib,neiblist);
      PrintArray("trueneiblist : ",Ntrueneib,trueneiblist);
      string message = "Problem with neighbour lists in Binary tree";
      ExceptionHandler::getIstance().raise(message);
    }
  }

  delete[] trueneiblist;

  return;
}
#endif



template class BinaryTree<1>;
template class BinaryTree<2>;
template class BinaryTree<3>;
