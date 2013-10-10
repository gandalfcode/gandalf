//=============================================================================
//  BinarySubTree.cpp
//  Contains all functions for building, stocking and walking for the 
//  binary tree for SPH particles.
//
//  This file is part of GANDALF :
//  Graphical Astrophysics code for N-body Dynamics and Lagrangian Fluids
//  https://github.com/gandalfcode/gandalf
//  Contact : gandalfcode@gmail.com
//
//  Copyright (C) 2013  D. A. Hubber, G Rosotti
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
using namespace std;



//=============================================================================
//  BinarySubTree::BinarySubTree
/// BinarySubTree constructor.  Initialises various variables.
//=============================================================================
template <int ndim>
BinarySubTree<ndim>::BinarySubTree(int Nleafmaxaux, FLOAT thetamaxsqdaux,
			     FLOAT kernrangeaux, string gravity_mac_aux,
                             string multipole_aux)
{
  allocated_tree = false;
  Ncell = 0;
  Ncellmax = 0;
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
//  BinarySubTree::~BinarySubTree
/// BinarySubTree destructor.  Deallocates tree memory upon object destruction.
//=============================================================================
template <int ndim>
BinarySubTree<ndim>::~BinarySubTree()
{
  if (allocated_tree) DeallocateSubTreeMemory();
}



//=============================================================================
//  BinarySubTree::AllocateSubTreeMemory
/// Allocate memory for binary tree as requested.  If more memory is required 
/// than currently allocated, tree is deallocated and reallocated here.
//=============================================================================
template <int ndim>
void BinarySubTree<ndim>::AllocateSubTreeMemory(void)
{
  debug2("[BinarySubTree::AllocateTreeMemory]");

  if (Ntotmax > Ntotmaxold || (!allocated_tree)) {
	//cout << "CHECKING : " << Ntotmax << "   " << Ntotmaxold << "    " << allocated_tree << endl;
    if (allocated_tree) DeallocateSubTreeMemory();
    ids = new int[Ntotmax];
    inext = new int[Ntotmax];
    pc = new int[Ntotmax];
    g2c = new int[gtot];
    tree = new struct BinaryTreeCell<ndim>[Ncellmax];
    for (int k=0; k<ndim; k++) porder[k] = new int[Ntotmax]; 
    for (int k=0; k<ndim; k++) rk[k] = new FLOAT[Ntotmax];
    allocated_tree = true;
    Ntotmaxold = Ntotmax;
    //cout << "Allocated " << Ntotmax << " for subtree" << endl;
  }

  return;
}



//=============================================================================
//  BinarySubTree::DeallocateSubTreeMemory
/// Deallocates all binary tree memory
//=============================================================================
template <int ndim>
void BinarySubTree<ndim>::DeallocateSubTreeMemory(void)
{
  debug2("[BinarySubTree::DeallocateTreeMemory]");

  if (allocated_tree) {
    for (int k=ndim-1; k>=0; k--) delete[] rk[k];
    for (int k=ndim-1; k>=0; k--) delete[] porder[k];
    delete[] tree;
    delete[] g2c;
    delete[] pc;
    delete[] inext;
    delete[] ids;
    allocated_tree = false;
  }

  return;
}



//=============================================================================
//  BinarySubTree::UpdateTree
/// Call all routines to build/re-build the binary tree.
//=============================================================================
template <int ndim>
void BinarySubTree<ndim>::BuildSubTree
(Sph<ndim> *sph,                    ///< Pointer to main SPH object
 Parameters &simparams)             ///< Simulation parameters
{
  int output = 0;

  debug2("[BinarySubTree::BuildSubTree]");

  //cout << "MEMORY : " << Ntot << "    " << Ntotmax << "    " << Ntotmaxold << "     " << allocated_tree << endl;

  // Allocate (or reallocate if needed) all tree memory
  AllocateSubTreeMemory();

  // Create tree data structure including linked lists and cell pointers
  CreateSubTreeStructure();

  // Find ordered list of particle positions ready for adding particles to tree
  OrderParticlesByCartCoord(sph->sphdata);

  // Now add particles to tree depending on Cartesian coordinates
  LoadParticlesToSubTree();

  // Calculate all cell quantities (e.g. COM, opening distance)
  StockCellProperties(sph->sphdata);

  // Validate tree structure
#if defined(VERIFY_ALL)
  ValidateTree(sph);
#endif

  return;
}



//=============================================================================
//  BinarySubTree::UpdateActiveParticleCounters
/// Loop through all leaf cells in binary tree and update all active 
/// particle counters.
//=============================================================================
template <int ndim>
void BinarySubTree<ndim>::UpdateActiveParticleCounters(Sph<ndim> *sph)
{
  return;
}



//=============================================================================
//  BinarySubTree::ComputeSubTreeSize
/// Compute the maximum size (i.e. no. of levels, cells and leaf cells) of 
/// the binary tree.
//=============================================================================
template <int ndim>
void BinarySubTree<ndim>::ComputeSubTreeSize(void)
{
  debug2("[BinarySubTree::ComputeTreeSize]");

  // Increase level until tree can contain all particles
  ltot = 0;
  while (Nleafmax*pow(2,ltot) < Ntotmax) {
    ltot++;
  };

  // Set total number of leaf/grid cells and tree cells
  gtot = pow(2,ltot);
  Ncell = 2*gtot - 1;
  Ncellmax = Ncell;
  Ntotmax = Nleafmax*pow(2,ltot);

  // Optional output (for debugging)
#if defined(VERIFY_ALL)
  cout << "Calculating tree size variables" << endl;
  cout << "Max. no of particles in leaf-cell : " << Nleafmax << endl;
  cout << "Max. no of particles in tree : " << Ntotmax << endl;
  cout << "No. of particles in tree : " << Ntot << endl;
  cout << "No. of levels on tree : " << ltot << endl;
  cout << "No. of tree cells : " << Ncell << endl;
  cout << "No. of grid cells : " << gtot << endl;
#endif

  return;
}



//=============================================================================
//  BinarySubTree::CreateSubTreeStructure
/// Create the raw tree skeleton structure once the tree size is known.
/// Sets all cell pointer variables and all cell levels.
//=============================================================================
template <int ndim>
void BinarySubTree<ndim>::CreateSubTreeStructure(void)
{
  debug2("[BinarySubTree::CreateTreeStructure]");

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
  for (g=0; g<gtot; g++) g2c[g] = 0;
  for (c=0; c<Ncell; c++) {
    tree[c].c2 = 0;
    tree[c].c2g = 0;
  }
  g = 0;
  tree[0].clevel = 0;

  // Loop over all cells and set all other pointers
  // --------------------------------------------------------------------------
  for (c=0; c<Ncell; c++) {
    if (tree[c].clevel == ltot) {                   // If on leaf level
      tree[c].cnext = c + 1;                        // id of next cell
      tree[c].c2g = g;                              // Record leaf id
      g2c[g++] = c;                                 // Record inverse id
    }
    else {
      tree[c+1].clevel = tree[c].clevel + 1;        // Level of 1st child
      tree[c].c2 = c + c2L[tree[c].clevel];         // id of 2nd child
      tree[tree[c].c2].clevel = tree[c].clevel + 1; // Level of 2nd child
      tree[c].cnext = c + cNL[tree[c].clevel];      // Next cell id
    }

  }
  // --------------------------------------------------------------------------


  // Free locally allocated memory
  delete[] cNL;
  delete[] c2L;

  return;
}



//=============================================================================
//  BinarySubTree::OrderParticlesByCartCoord
/// Compute lists of particle ids in order of x, y and z positions.
/// Used for efficiently creating the binary tree data structure.
//=============================================================================
template <int ndim>
void BinarySubTree<ndim>::OrderParticlesByCartCoord
(SphParticle<ndim> *sphdata)        ///< SPH particle data array
{
  int i;                            // Particle counter
  int j;                            // ..
  int k;                            // Dimension counter

  debug2("[BinarySubTree::OrderParticlesByCartCoord]");

  // First copy all values to local arrays
  for (k=0; k<ndim; k++) {
    for (j=0; j<Ntot; j++) {
      i = ids[j];
      //cout << "WTF?? : " << k << "   " << j << "    " << i << endl;
      rk[k][j] = sphdata[i].r[k];
    }
  }

  // Now copy list of particle ids
  for (k=0; k<ndim; k++)
    for (i=0; i<Ntot; i++)
      porder[k][i] = i;

  // Divide sorting routines amongst (up to 3) threads
  // --------------------------------------------------------------------------
  {
    // Sort x-values
    Heapsort(Ntot,porder[0],rk[0]);

    // Sort y-values
    if (ndim >= 2) Heapsort(Ntot,porder[1],rk[1]);

    // Sort z-values
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
      i1 = GlobalId(i1);
      i2 = GlobalId(i2);
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
//  BinarySubTree::LoadParticlesToTree
/// Create tree structure by adding particles to leaf cells.
//=============================================================================
template <int ndim>
void BinarySubTree<ndim>::LoadParticlesToSubTree(void)
{
  int c;                            // Cell counter
  int cc;                           // Secondary cell counter
  int k;                            // Dimensionality counter
  int i;                            // Particle counter
  int j;                            // Dummy particle id
  int l;                            // Level counter
  FLOAT *ccap;                      // Maximum capacity of cell
  FLOAT *ccon;                      // Current contents of cell

  debug2("[BinarySubTree::LoadParticleToTree]");

  // Allocate memory for local arrays
  ccap = new FLOAT[Ncellmax];
  ccon = new FLOAT[Ncellmax];

  // Set capacity of root-cell using particle weights
  //for (i=0; i<Ntot; i++) pw[i] = 1.0/(FLOAT) Ntot;
  for (c=0; c<Ncell; c++) ccap[c] = 0.0;
  for (i=0; i<Ntot; i++) ccap[0] += 1.0; //pw[i];

  // Initialise all particle and cell values before building tree structure
  for (i=0; i<Ntot; i++) pc[i] = 0;
  for (i=0; i<Ntot; i++) inext[i] = -1;
  for (c=0; c<Ncell; c++) ccon[c] = 0.0;
  for (c=0; c<Ncell; c++) tree[c].ifirst = -1;
  for (c=0; c<Ncell; c++) tree[c].ilast = -1;

  for (i=0; i<Ntot; i++) assert(ids[i] >= 0);

  // Start at top level (l = 0) dividing the cell along the x-axis (k = 0)
  l = 0;
  k = 0;


  // Loop through each level of the tree
  // --------------------------------------------------------------------------
  while (l < ltot) {

    // Loop over all particles (in order of current split)
    // ------------------------------------------------------------------------
    for (i=0; i<Ntot; i++) {
      j = porder[k][i];
      cc = pc[j];                            // Cell currently occupied by j
      ccon[cc] += 1.0; //pw[j];              // Add particle weighting to cell

      // If cell contains less than maximum allowed capacity, then add 
      // particle to the first child cell.  Otherwise add to second child.
      if (ccon[cc] < 0.5000000000001*ccap[cc]) {
        pc[j]++;
        ccap[pc[j]] += 1.0; //pw[j];
      }
      else {
        pc[j] = tree[cc].c2;
        ccap[pc[j]] += 1.0; //pw[j];
      }
    }
    // ------------------------------------------------------------------------

    // Move to next level and cycle through each dimension in turn
    // (Need more sophisticated algorithm here in future)
    l++;
    k = (k + 1)%ndim;

  }
  // --------------------------------------------------------------------------


  // Compute capacities of leaf cells here
  for (i=0; i<Ntot; i++) {
    cc = pc[i];
    ccon[cc] += 1.0; //pw[i];
  }


  // Loop over all particles and set id of first particle in each cell, plus 
  // the linked list values
  for (i=0; i<Ntot; i++) {
    c = pc[i];
    if (tree[c].ifirst == -1)
      tree[c].ifirst = i;
    else
      inext[tree[c].ilast] = i;
    tree[c].ilast = i;
  }


  // Free all locally allocated memory
  delete[] ccon;
  delete[] ccap;

  return;
}



//=============================================================================
//  BinarySubTree::StockCellProperties
/// Calculate the physical properties (e.g. total mass, centre-of-mass, 
/// opening-distance, etc..) of all cells in the tree.
//=============================================================================
template <int ndim>
void BinarySubTree<ndim>::StockCellProperties
(SphParticle<ndim> *sphdata)        ///< SPH particle data array
{
  int c,cc,ccc;                     // Cell counters
  int i;                            // Particle counter
  int j;                            // ..
  int k;                            // Dimension counter
  FLOAT dr[ndim];                   // Relative position vector
  FLOAT drsqd;                      // Distance squared
  FLOAT factor = 1.0/thetamaxsqd;   // Geometric MAC aux. variable
  FLOAT mi;                         // Mass of particle i
  FLOAT *crmax;                     // Max. extent of cell bounding boxes
  FLOAT *crmin;                     // Min. extent of cell bounding boxes

  debug2("[BinarySubTree::StockCellProperties]");

  // Allocate local memory
  crmax = new FLOAT[ndim*Ncellmax];
  crmin = new FLOAT[ndim*Ncellmax];

  // Zero all summation variables for all cells
  for (c=0; c<Ncell; c++) {
    tree[c].Nactive = 0;
    tree[c].N = 0;
    tree[c].m = 0.0;
    tree[c].hmax = 0.0;
    tree[c].rmax = 0.0;
    tree[c].cdistsqd = big_number;
    for (k=0; k<ndim; k++) tree[c].r[k] = 0.0;
    for (k=0; k<5; k++) tree[c].q[k] = 0.0;
  }

  for (c=0; c<Ncell; c++)
    for (k=0; k<ndim; k++) crmin[c*ndim + k] = big_number;

  for (c=0; c<Ncell; c++)
    for (k=0; k<ndim; k++) crmax[c*ndim + k] = -big_number;


  // Loop backwards over all tree cells to ensure child cells are always 
  // computed first before being summed in parent cells.
  // ==========================================================================
  for (c=Ncell-1; c>=0; c--) {

    // If this is a leaf cell, sum over all particles
    // ------------------------------------------------------------------------
    if (tree[c].c2 == 0) {
      j = tree[c].ifirst;

      // Loop over all particles in cell summing their contributions
      while (j != -1) {
        i = GlobalId(j);
        tree[c].N++;
        if (sphdata[i].active) tree[c].Nactive++;
        tree[c].hmax = max(tree[c].hmax,sphdata[i].h);
        tree[c].m += sphdata[i].m;
        for (k=0; k<ndim; k++) tree[c].r[k] += sphdata[i].m*sphdata[i].r[k];
        for (k=0; k<ndim; k++) {
          if (sphdata[i].r[k] < crmin[c*ndim + k])
            crmin[c*ndim + k] = sphdata[i].r[k];
          if (sphdata[i].r[k] > crmax[c*ndim + k])
            crmax[c*ndim + k] = sphdata[i].r[k];
        }
        j = inext[j];
      };

      // Normalise all cell values
      if (tree[c].N > 0) {
        for (k=0; k<ndim; k++) tree[c].r[k] /= tree[c].m;
        for (k=0; k<ndim; k++) dr[k] = 0.5*(crmax[c*ndim + k] - crmin[c*ndim + k]);
        tree[c].cdistsqd = factor*DotProduct(dr,dr,ndim);
        for (k=0; k<ndim; k++) dr[k] = max(crmax[c*ndim + k] - tree[c].r[k],
		  			 tree[c].r[k] - crmin[c*ndim + k]);
        tree[c].rmax = sqrt(DotProduct(dr,dr,ndim));
      }

      // Compute quadrupole moment terms if selected
      if (multipole == "quadrupole") {
        j = tree[c].ifirst;

        while (j != -1) {
          i = GlobalId(j);
          mi = sphdata[i].m;
          for (k=0; k<ndim; k++) dr[k] = sphdata[i].r[k] - tree[c].r[k];
          drsqd = DotProduct(dr,dr,ndim);
          if (ndim == 3) {
            tree[c].q[0] += mi*(3.0*dr[0]*dr[0] - drsqd);
            tree[c].q[1] += mi*3.0*dr[0]*dr[1];
            tree[c].q[2] += mi*(3.0*dr[1]*dr[1] - drsqd);
            tree[c].q[3] += mi*3.0*dr[2]*dr[0];
            tree[c].q[4] += mi*3.0*dr[2]*dr[1];
          }
          else if (ndim == 2) {
            tree[c].q[0] += mi*(3.0*dr[0]*dr[0] - drsqd);
            tree[c].q[1] += mi*3.0*dr[0]*dr[1];
            tree[c].q[2] += mi*(3.0*dr[1]*dr[1] - drsqd);
          }
          j = inext[j];
        }
      }

    }
    // For non-leaf cells, sum together two children cells
    // ------------------------------------------------------------------------
    else {
      cc = c + 1;
      ccc = tree[c].c2;
      tree[c].N = tree[cc].N + tree[ccc].N;

      if (tree[c].N > 0) {
        tree[c].hmax = max(tree[cc].hmax,tree[ccc].hmax);
        tree[c].m = tree[cc].m + tree[ccc].m;
        for (k=0; k<ndim; k++) tree[c].r[k] =
          (tree[cc].m*tree[cc].r[k] + tree[ccc].m*tree[ccc].r[k])/tree[c].m;
        for (k=0; k<ndim; k++)
          crmin[ndim*c + k] = min(crmin[ndim*cc+k],crmin[ndim*ccc+k]);
        for (k=0; k<ndim; k++)
          crmax[ndim*c + k] = max(crmax[ndim*cc+k],crmax[ndim*ccc+k]);
        for (k=0; k<ndim; k++) dr[k] = 0.5*(crmax[c*ndim + k] - crmin[c*ndim + k]);
        tree[c].cdistsqd = factor*DotProduct(dr,dr,ndim);
        for (k=0; k<ndim; k++) dr[k] = max(crmax[c*ndim + k] - tree[c].r[k],
                                           tree[c].r[k] - crmin[c*ndim + k]);
        tree[c].rmax = sqrt(DotProduct(dr,dr,ndim));
      }

      // Now add individual quadrupole moment terms
      if (multipole == "quadrupole" && tree[cc].N > 0) {
        mi = tree[cc].m;
        for (k=0; k<ndim; k++) dr[k] = tree[cc].r[k] - tree[c].r[k];
        drsqd = DotProduct(dr,dr,ndim);
        if (ndim == 3) {
          tree[c].q[0] += mi*(3.0*dr[0]*dr[0] - drsqd);
          tree[c].q[1] += mi*3.0*dr[0]*dr[1];
          tree[c].q[2] += mi*(3.0*dr[1]*dr[1] - drsqd);
          tree[c].q[3] += mi*3.0*dr[2]*dr[0];
          tree[c].q[4] += mi*3.0*dr[2]*dr[1];
        }
        else if (ndim == 2) {
          tree[c].q[0] += mi*(3.0*dr[0]*dr[0] - drsqd);
          tree[c].q[1] += mi*3.0*dr[0]*dr[1];
          tree[c].q[2] += mi*(3.0*dr[1]*dr[1] - drsqd);
        }
      }

      if (multipole == "quadrupole" && tree[ccc].N > 0) {
        mi = tree[ccc].m;
        for (k=0; k<ndim; k++) dr[k] = tree[ccc].r[k] - tree[c].r[k];
        drsqd = DotProduct(dr,dr,ndim);
        if (ndim == 3) {
          tree[c].q[0] += mi*(3.0*dr[0]*dr[0] - drsqd);
          tree[c].q[1] += mi*3.0*dr[0]*dr[1];
          tree[c].q[2] += mi*(3.0*dr[1]*dr[1] - drsqd);
          tree[c].q[3] += mi*3.0*dr[2]*dr[0];
          tree[c].q[4] += mi*3.0*dr[2]*dr[1];
        }
        else if (ndim == 2) {
          tree[c].q[0] += mi*(3.0*dr[0]*dr[0] - drsqd);
          tree[c].q[1] += mi*3.0*dr[0]*dr[1];
          tree[c].q[2] += mi*(3.0*dr[1]*dr[1] - drsqd);
        }
      }

    }
    // ------------------------------------------------------------------------

  }
  // ==========================================================================


#if defined(VERIFY_ALL)
  cout << "Root cell position1 : " << tree[0].r[0] << "   " 
       << tree[0].r[1] << endl;
  cout << "Mass of root cell1 : " << tree[0].m << endl;
  cout << "Bounding box : " << crmin[0] << "   " << crmax[0] 
       << "   " << crmin[1] << "   " << crmax[1] << endl;

  // Go through each tree level in turn and print info
  /*for (int l=0; l<ltot+1; l++) {
    cout << "LEVEL : " << l << endl;
    cout << "----------------------" << endl;
    for (c=0; c<Ncell; c++) {
      if (tree[c].clevel == l) {
	    cout << "c : " << c << "   N : " 
      	   << tree[c].N << "   " << Nleafmax << endl;
      }
    }
  }*/
#endif


  // Free all locally allocated memory
  delete[] crmin;
  delete[] crmax;

  return;
}



//=============================================================================
//  BinarySubTree::UpdateHmaxValues
/// Calculate the physical properties (e.g. total mass, centre-of-mass, 
/// opening-distance, etc..) of all cells in the tree.
//=============================================================================
template <int ndim>
FLOAT BinarySubTree<ndim>::UpdateHmaxValues
(SphParticle<ndim> *sphdata)        ///< SPH particle data array
{
  int c,cc,ccc;                     // Cell counters
  int i;                            // Particle counter
  int j;

  debug2("[BinarySubTree::UpdateHmaxValues]");

  // Zero all summation variables for all cells
  for (c=0; c<Ncell; c++) tree[c].hmax = 0.0;

  // Loop backwards over all tree cells to ensure child cells are always 
  // computed first before being summed in parent cells.
  // ==========================================================================
  for (c=Ncell-1; c>=0; c--) {

    // If this is a leaf cell, sum over all particles
    // ------------------------------------------------------------------------
    if (tree[c].c2 == 0) {
      i = tree[c].ifirst;

      // Loop over all particles in cell summing their contributions
      while (i != -1) {
    	j = GlobalId(i);
        tree[c].hmax = max(tree[c].hmax,sphdata[j].h);
        i = inext[i];
      };

    }
    // For non-leaf cells, sum together two children cells
    // ------------------------------------------------------------------------
    else {
      cc = c + 1;
      ccc = tree[c].c2;
      if (tree[c].N > 0) tree[c].hmax = max(tree[cc].hmax,tree[ccc].hmax);

    }
    // ------------------------------------------------------------------------

  }
  // ==========================================================================

  return tree[0].hmax;
}



//=============================================================================
//  BinarySubTree::ComputeActiveCellList
/// Returns the number of cells containing active particles, 'Nactive', and
/// the i.d. list of cells contains active particles, 'celllist'
//=============================================================================
template <int ndim>
int BinarySubTree<ndim>::ComputeActiveCellList
(int Nactive,
 BinaryTreeCell<ndim> **celllist)  ///< Cells id array containing active ptcls
{
  int c;                           // Cell counter

  for (c=0; c<Ncell; c++)
    if (tree[c].Nactive > 0) celllist[Nactive++] = &tree[c];

  return Nactive;
}



//=============================================================================
//  BinarySubTree::ComputeGatherNeighbourList
/// Computes and returns number of neighbour, 'Nneib', and the list
/// of neighbour ids, 'neiblist', for all particles inside cell 'c'.
/// Includes all particles in the selected cell, plus all particles
/// contained in adjacent cells (including diagonal cells).
//=============================================================================
template <int ndim>
int BinarySubTree<ndim>::ComputeGatherNeighbourList
(BinaryTreeCell<ndim> *cell,        ///< [in] Pointer
 int Nneib,                         ///< [in] ..
 int Nneibmax,                      ///< [in] Max. no. of neighbours
 int *neiblist,                     ///< [out] List of neighbour i.d.s
 FLOAT hmax,                        ///< [in] Maximum smoothing length
 SphParticle<ndim> *sphdata)        ///< [in] SPH particle data
{
  int cc;                           // Cell counter
  int i;                            // Particle id
  int j;                            // Aux. particle counter
  int k;                            // Neighbour counter
  int Ntemp = Nneib;                // Aux. neighbour counter
  FLOAT dr[ndim];                   // Relative position vector
  FLOAT drsqd;                      // Distance squared
  FLOAT rc[ndim];                   // Position of cell
  FLOAT hrangemax;                  // Maximum SPH kernel extent
  //FLOAT neibrange;                  // Max. range of gather neighbours

  for (k=0; k<ndim; k++) rc[k] = cell->r[k];
  hrangemax = cell->rmax + kernrange*hmax;

  // Start with root cell and walk through entire tree
  cc = 0;

  // ==========================================================================
  while (cc < Ncell) {
    for (k=0; k<ndim; k++) dr[k] = tree[cc].r[k] - rc[k];
    drsqd = DotProduct(dr,dr,ndim);
    //neibrange = tree[cc].rmax + hrangemax;

    // Check if circular range overlaps cell bounding sphere
    // ------------------------------------------------------------------------
    if (drsqd < pow(tree[cc].rmax + hrangemax,2)) {

      // If not a leaf-cell, then open cell to first child cell
      if (tree[cc].c2 != 0)
        cc++;

      // If leaf-cell, add particles to list
      else if (tree[cc].c2 == 0 && Nneib + Nleafmax < Nneibmax) {
        i = tree[cc].ifirst;
    	while (i != -1) {
          neiblist[Nneib++] = i;
          if (i == tree[cc].ilast);
    	  i = inext[i];
        };
        cc = tree[cc].cnext;
      }

      // If leaf-cell, but we've run out of memory, return with error-code (-1)
      else if (tree[cc].c2 == 0 && Nneib + Nleafmax >= Nneibmax)
    	return -1;

    }

    // If not in range, then open next cell
    // ------------------------------------------------------------------------
    else
      cc = tree[cc].cnext;

  };
  // ==========================================================================


  // Now, trim the list to remove particles that are definitely not neighbours
  hrangemax = hrangemax*hrangemax;
  for (j=Ntemp; j<Nneib; j++) {
    i = GlobalId(neiblist[j]);
    //cout << "i : " << i << "    " << j << "     " << Nneibmax << "    " << neiblist[j] << endl;
    for (k=0; k<ndim; k++) dr[k] = sphdata[i].r[k] - rc[k];
    drsqd = DotProduct(dr,dr,ndim);
    if (drsqd < hrangemax) neiblist[Ntemp++] = i;
  }
  Nneib = Ntemp;

  return Nneib;
}



//=============================================================================
//  BinarySubTree::ComputeNeighbourList
/// Computes and returns number of neighbour, 'Nneib', and the list
/// of neighbour ids, 'neiblist', for all particles inside cell 'c'.
/// Includes all particles in the selected cell, plus all particles
/// contained in adjacent cells (including diagonal cells).
//=============================================================================
template <int ndim>
int BinarySubTree<ndim>::ComputeNeighbourList
(BinaryTreeCell<ndim> *cell,        ///< [in] Cell pointer
 int Nneib,                         ///< [in] ..
 int Nneibmax,                      ///< [in] Max. no. of neighbours
 int *neiblist,                     ///< [out] List of neighbour i.d.s
 SphParticle<ndim> *sphdata)        ///< [in] SPH particle data
{
  int cc;                           // Cell counter
  int i;                            // Particle id
  int j;                            // Aux. particle counter
  int k;                            // Neighbour counter
  int Ntemp = Nneib;                // ..
  FLOAT dr[ndim];                   // Relative position vector
  FLOAT drsqd;                      // Distance squared
  FLOAT rc[ndim];                   // Position of cell
  FLOAT hrangemax;                  // Maximum kernel extent
  FLOAT rmax;                       // ..
  //FLOAT gatherrange;                // ..
  //FLOAT scatterrange;               // ..

  for (k=0; k<ndim; k++) rc[k] = cell->r[k];
  hrangemax = cell->rmax + kernrange*cell->hmax;
  rmax = cell->rmax;

  // Start with root cell and walk through entire tree
  cc = 0;

  // ==========================================================================
  while (cc < Ncell) {
    for (k=0; k<ndim; k++) dr[k] = tree[cc].r[k] - rc[k];
    drsqd = DotProduct(dr,dr,ndim);
    //scatterrange = tree[cc].rmax + rmax + kernrange*tree[cc].hmax;
    //gatherrange = tree[cc].rmax + hrangemax;

    // Check if circular range overlaps cell bounding sphere
    // ------------------------------------------------------------------------
    //if (drsqd < gatherrange*gatherrange ||
      //  drsqd < scatterrange*scatterrange) {
    if (drsqd < pow(tree[cc].rmax + hrangemax,2) ||
        drsqd < pow(rmax + tree[cc].rmax + kernrange*tree[cc].hmax,2)) {

      // If not a leaf-cell, then open cell to first child cell
      if (tree[cc].c2 != 0)
        cc++;

      // If leaf-cell, add particles to list
      else if (tree[cc].c2 == 0 && Nneib + Nleafmax < Nneibmax) {
        i = tree[cc].ifirst;
    	while (i != -1) {
          neiblist[Nneib++] = i;
    	  i = inext[i];
        };
        cc = tree[cc].cnext;
      }

      // If leaf-cell, but we've run out of memory, return with error-code (-1)
      else if (tree[cc].c2 == 0 && Nneib + Nleafmax >= Nneibmax)
    	return -1;

    }

    // If not in range, then open next cell
    // ------------------------------------------------------------------------
    else
      cc = tree[cc].cnext;
  };
  // ==========================================================================


  // Now, trim the list to remove particles that are definitely not neighbours
  hrangemax = hrangemax*hrangemax;
  for (j=Ntemp; j<Nneib; j++) {
    i = GlobalId(neiblist[j]);
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
//=============================================================================
template <int ndim>
int BinarySubTree<ndim>::ComputeGravityInteractionList
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
  int cc;                           // Cell counter
  int i;                            // Particle id
  int j;                            // Aux. particle counter
  int ilast;                        // id of last particle in current cell
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
  //Ndirect = 0;
  //Ngravcell = 0;


  // Walk through all cells in tree to determine particle and cell 
  // interaction lists
  // ==========================================================================
  while (cc < Ncell) {
    for (k=0; k<ndim; k++) dr[k] = tree[cc].r[k] - rc[k];
    drsqd = DotProduct(dr,dr,ndim);

    // Check if cells contain SPH neighbours
    // ------------------------------------------------------------------------
    if (drsqd < (tree[cc].rmax + hrangemax)*(tree[cc].rmax + hrangemax) ||
        drsqd < pow(tree[cc].rmax + rmax + kernrange*tree[cc].hmax,2)) {

      // If not a leaf-cell, then open cell to first child cell
      if (tree[cc].c2 != 0)
        cc++;

      // If leaf-cell, add particles to list
      else if (tree[cc].c2 == 0 && Nneib + Nleafmax <= Nneibmax) {
        i = tree[cc].ifirst;
    	while (i != -1) {
          neiblist[Nneib++] = i;
    	  i = inext[i];
        };
        cc = tree[cc].cnext;
      }

      // If leaf-cell, but we've run out of memory, return with error-code (-1)
      else if (tree[cc].c2 == 0 && Nneib + Nleafmax > Nneibmax)
    	return -1;

    }

    // Check if cell is far enough away to use the COM approximation
    // ------------------------------------------------------------------------
    //else if (drsqd >= tree[c].cdistsqd + tree[cc].cdistsqd) {
    else if (drsqd > cdistsqd || drsqd > tree[cc].cdistsqd) {

      // If cell is a leaf-cell with only one particle, more efficient to
      // compute the gravitational contribution from the particle than the cell
      if (tree[cc].c2 == 0 && tree[cc].N == 1 && Ndirect < Ndirectmax)
        directlist[Ndirect++] = tree[cc].ifirst;
      else if (Ngravcell < Ngravcellmax)
        gravcelllist[Ngravcell++] = &(tree[cc]);
      else
        return -1;
      cc = tree[cc].cnext;

    }

    // If cell is too close, open cell to interogate children cells.
    // If cell is too close and a leaf cell, then add particles to direct list.
    // ------------------------------------------------------------------------
    //else if (drsqd < tree[c].cdistsqd + tree[cc].cdistsqd) {
    else if (drsqd <= cdistsqd && drsqd <= tree[cc].cdistsqd) {

      // If not a leaf-cell, then open cell to first child cell
      if (tree[cc].c2 != 0)
         cc++;

      // If leaf-cell, add particles to list
      else if (tree[cc].c2 == 0 && Ndirect + Nleafmax <= Ndirectmax) {
        i = tree[cc].ifirst;
        while (i != -1) {
          directlist[Ndirect++] = i;
       	  i = inext[i];
        };
        cc = tree[cc].cnext;
      }

      // If leaf-cell, but we've run out of memory, return with error-code (-1)
      else if (tree[cc].c2 == 0 && Ndirect + Nleafmax > Ndirectmax)
       	return -1;

    }

    // If not in range, then open next cell
    // ------------------------------------------------------------------------
    else
      cc = tree[cc].cnext;

  };
  // ==========================================================================

  // Now, trim the list to remove particles that are definitely not neighbours.
  // If not an SPH neighbour, then add to direct gravity sum list.
  hrangemax = hrangemax*hrangemax;
  for (j=0; j<Nneib; j++) {
    i = GlobalId(neiblist[j]);
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



#if defined(VERIFY_ALL)
//=============================================================================
//  BinarySubTree::ValidateTree
/// Perform various consistency checks to ensure the tree structure and all 
/// cell values are valid.
//=============================================================================
template <int ndim>
void BinarySubTree<ndim>::ValidateTree
(Sph<ndim> *sph)                    ///< [in] SPH object pointer
{
  bool treeflag;                    // ..
  int c;                            // .. 
  int cc;                           // ..
  int i;                            // ..
  int j;                            // ..
  int k;                            // ..
  int N;                            // ..
  FLOAT dr[ndim];                   // ..
  FLOAT drmag;                      // ..
  FLOAT mc;                         // ..

  debug2("[BinarySubTree::ValidateTree]");


  // Check all tree cells are on valid levels
  // --------------------------------------------------------------------------
  for (c=0; c<Ncell; c++) {
    if (tree[c].clevel > ltot) {
      cout << "Problem with tree levels : " << cc << "   " << tree[cc].clevel
           << "    " << ltot << endl;
      exit(0);
    }
  }


  // Check all leaf cells the correct number of particles
  // --------------------------------------------------------------------------
  for (c=0; c<Ncell; c++) {
    if (tree[c].c2 == 0) {
      N = 0;
      i = tree[c].ifirst;
      while (i != -1) {
        N++;
        assert(i >= 0 && i < Ntot);
        assert(ids[i] >= 0);
        i = inext[i];
      };
      if (N > Nleafmax) {
        cout << "Problem with leaf cells : " << N
             << "   " << Nleafmax << "   " << c << endl;
        exit(0);
      }
    }
  }


  // Walk all cells in tree to compute quantities
  // --------------------------------------------------------------------------
  for (c=0; c<Ncell; c++) {

    treeflag = true;
    cc = c;
    mc = 0;

    // Now loop over all child cells below cell to sum all properties
    // ------------------------------------------------------------------------
    while (cc < tree[c].cnext) {

      if (tree[cc].c2 == 0) {

    	i = tree[cc].ifirst;
    	while (i != -1) {
    	  j = GlobalId(i);
          for (k=0; k<ndim; k++) dr[k] = tree[c].r[k] - sph->sphdata[j].r[k];
          mc += sph->sphdata[j].m;
          drmag = sqrt(DotProduct(dr,dr,ndim));
          if (drmag > 1.00001*tree[c].rmax) treeflag = false;
          if (sph->sphdata[j].h > 1.00001*tree[c].hmax) treeflag = false;

          if (!treeflag) {
    	    cout << "Problem with tree : " << c << "   " 
                 << cc << "   " << i << endl;
    	    cout << "rc : " << tree[c].r[0] << "   " << tree[c].r[1] << endl;
    	    cout << "hmax : " << tree[c].hmax << "   rmax : " 
                 << tree[c].rmax << endl;
    	    cout << "rp : " << sph->sphdata[j].r[0] << "    "
                 << sph->sphdata[j].r[1] << endl;
    	    cout << "h : " << sph->sphdata[j].h << endl;
    	    cout << "drmag : " << drmag << "    rmax : " 
                 << tree[c].rmax << endl;

    	    exit(0);
          }
          i = inext[i];
    	};

      }

      cc++;
    }
    // ------------------------------------------------------------------------

    if (fabs(mc - tree[c].m) > 0.00001*mc) {
      cout << "mass : " << mc << "    " << tree[c].m << endl;
      exit(0);
    }
  }
  // --------------------------------------------------------------------------

  return;
}
#endif




template class BinarySubTree<1>;
template class BinarySubTree<2>;
template class BinarySubTree<3>;
