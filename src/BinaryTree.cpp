//=============================================================================
//  BinaryTree.cpp
//  ..
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
//  BinaryTree::BinaryTree
/// BinaryTree constructor.  Initialises various variables.
//=============================================================================
template <int ndim>
BinaryTree<ndim>::BinaryTree(int Nleafmaxaux, FLOAT thetamaxsqdaux, 
			     string gravity_mac_aux)
{
  allocated_tree = false;
  Ncell = 0;
  Ncellmax = 0;
  Ntot = 0;
  Ntotmax = 0;
  Ntotmaxold = 0;
  Nleafmax = Nleafmaxaux;
  thetamaxsqd = thetamaxsqdaux;
  gravity_mac = gravity_mac_aux;
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
/// than currently allocated, grid is deallocated and reallocated here.
//=============================================================================
template <int ndim>
void BinaryTree<ndim>::AllocateTreeMemory(int Npart)
{
  debug2("[BinaryTree::AllocateTreeMemory]");

  if (Ntotmax > Ntotmaxold || (!allocated_tree)) {
    if (allocated_tree) DeallocateTreeMemory();
    inext = new int[Ntotmax];
    pw = new FLOAT[Ntotmax];
    pc = new int[Ntotmax];
    g2c = new int[gtot];
    tree = new struct BinaryTreeCell<ndim>[Ncellmax];
    for (int k=0; k<ndim; k++) porder[k] = new int[Ntotmax]; 
    for (int k=0; k<ndim; k++) r[k] = new FLOAT[Ntotmax];
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
  debug2("[BinaryTree::DeallocateTreeMemory]");

  if (allocated_tree) {
    for (int k=ndim-1; k>=0; k--) delete[] r[k];
    for (int k=ndim-1; k>=0; k--) delete[] porder[k];
    delete[] tree;
    delete[] g2c;
    delete[] pc;
    delete[] pw;
    delete[] inext;
    allocated_tree = false;
  }

  return;
}



//=============================================================================
//  BinaryTree::UpdateTree
/// Call all routines to build/re-build the binary tree.
//=============================================================================
template <int ndim>
void BinaryTree<ndim>::UpdateTree(Sph<ndim> *sph, Parameters &simparams)
{
  int output = 0;

  // Set number of tree members to total number of SPH particles (inc. ghosts)
  Ntotmaxold = Ntotmax;
  Nsph = sph->Nsph;
  Ntot = sph->Ntot;

  // Compute the size of all tree-related arrays now we know number of points
  ComputeTreeSize(output,Ntot);

  // Allocate (or reallocate if needed) all tree memory
  AllocateTreeMemory(Ntot);

  // Create tree data structure including linked lists and cell pointers
  CreateTreeStructure(output);

  // Find ordered list of particle positions ready for adding particles to tree
  OrderParticlesByCartCoord(output,sph->sphdata);

  // Now add particles to tree depending on Cartesian coordinates
  LoadParticlesToTree(output);

  // Calculate all cell quantities (e.g. COM, opening distance)
  StockCellProperties(output,sph->sphdata);

  // Validate tree structure
#if defined(VERIFY_ALL)
  ValidateTree(sph);
#endif

  return;
}



//=============================================================================
//  BinaryTree::ComputeTreeSize
/// Compute the maximum size (i.e. no. of levels, cells and leaf cells) of 
/// the binary tree.
//=============================================================================
template <int ndim>
void BinaryTree<ndim>::ComputeTreeSize(int output, int Npart)
{
  debug2("[BinaryTree::ComputeTreeSize]");

  // Increase level until tree can contain all particles
  ltot = 0;
  while (Nleafmax*pow(2,ltot) < Ntot) {
    ltot++;
  };
  gtot = pow(2,ltot);
  Ncell = 2*gtot - 1;
  Ncellmax = Ncell;
  Ntotmax = Nleafmax*pow(2,ltot);

  // Optional output (for debugging)
  if (output == 1) {
    cout << "Calculating tree size variables" << endl;
    cout << "Max. no of particles in leaf-cell : " << Nleafmax << endl;
    cout << "Max. no of particles in tree : " << Ntotmax << endl;
    cout << "No. of particles in tree : " << Ntot << endl;
    cout << "No. of levels on tree : " << ltot << endl;
    cout << "No. of tree cells : " << Ncell << endl;
    cout << "No. of grid cells : " << gtot << endl;
  }

  return;
}



//=============================================================================
//  BinaryTree::CreateTreeStructure
/// Create the raw tree skeleton structure once the tree size is known.
/// Sets all cell pointer variables and all cell levels.
//=============================================================================
template <int ndim>
void BinaryTree<ndim>::CreateTreeStructure(int output)
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


#if defined(VERIFY_ALL)
  if (output == 1) {
    cout << "Diagnostic output from BinaryTree::CreateTreeStructure" << endl;
    for (c=0; c<min(20,Ncell-20); c++) {
      cout << "c : " << c << "   " << tree[c].clevel << "   " << tree[c].c2
	   << "   " << tree[c].cnext << "   " << tree[c].c2g << "   "
	   << g2c[tree[c].c2g] << endl;
    }
    for (c=max(20,Ncell-20); c<Ncell; c++) {
      cout << "c : " << c << "   " << tree[c].clevel << "   " << tree[c].c2
	   << "   " << tree[c].cnext << "   " << tree[c].c2g << "   "
	   << g2c[tree[c].c2g] << endl;
    }
  }
#endif


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
void BinaryTree<ndim>::OrderParticlesByCartCoord(int output,
						 SphParticle<ndim> *sphdata)
{
  int i;                            // Particle counter
  int k;                            // Dimension counter

  debug2("[BinaryTree::OrderParticlesByCartCoord]");

  // First copy all values to local arrays
  for (k=0; k<ndim; k++)
    for (i=0; i<Ntot; i++)
      r[k][i] = sphdata[i].r[k];

  // Now copy list of particle ids
  for (k=0; k<ndim; k++)
    for (i=0; i<Ntot; i++)
      porder[k][i] = i;

  // Sort x-values
  Heapsort(output,Ntot,porder[0],r[0]);

  // Sort y-values
  if (ndim >= 2) Heapsort(output,Ntot,porder[1],r[1]);

  // Sort z-values
  if (ndim == 3) Heapsort(output,Ntot,porder[2],r[2]);

  // Check that particles are ordered correctly
#if defined(VERIFY_ALL)
  for (k=0; k<ndim; k++) {
    for (int j=1; j<Ntot; j++) {
      int i1 = porder[k][j-1];
      int i2 = porder[k][j];
      if (sphdata[i2].r[k] < sphdata[i1].r[k]) {
        cout << "Problem with particle ordering : "
	      << k << "   " << j << endl;
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
void BinaryTree<ndim>::LoadParticlesToTree(int output)
{
  debug2("[BinaryTree::LoadParticleToTree]");

  int c;                            // Cell counter
  int cc;                           // Secondary cell counter
  int k;                            // Dimensionality counter
  int i;                            // Particle counter
  int j;                            // Dummy particle id
  int l;                            // Level counter
  FLOAT *ccap;                      // Maximum capacity of cell
  FLOAT *ccon;                      // Current contents of cell

  // Allocate memory for local arrays
  ccap = new FLOAT[Ncellmax];
  ccon = new FLOAT[Ncellmax];

  // Set capacity of root-cell using particle weights
  for (i=0; i<Ntot; i++) pw[i] = 1.0/(FLOAT) Ntot;
  for (c=0; c<Ncell; c++) ccap[c] = 0.0;
  for (i=0; i<Ntot; i++) ccap[0] += pw[i];

  // Initialise all particle and cell values before building tree structure
  for (i=0; i<Ntot; i++) pc[i] = 0;
  for (i=0; i<Ntot; i++) inext[i] = -1;
  for (c=0; c<Ncell; c++) ccon[c] = 0.0;
  for (c=0; c<Ncell; c++) tree[c].ifirst = -1;
  for (c=0; c<Ncell; c++) tree[c].ilast = -1;

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
      ccon[cc] += pw[j];                     // Add particle weighting to cell

      // If cell contains less than maximum allowed capacity, then add 
      // particle to the first child cell.  Otherwise add to second child.
      if (ccon[cc] < 0.5000000000001*ccap[cc]) {
        pc[j]++;
        ccap[pc[j]] += pw[j];
      }
      else {
        pc[j] = tree[cc].c2;
        ccap[pc[j]] += pw[j];
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
    ccon[cc] += pw[i];
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
//  BinaryTree::StockCellProperties
/// Calculate the physical properties (e.g. total mass, centre-of-mass, 
/// opening-distance, etc..) of all cells in the tree.
//=============================================================================
template <int ndim>
void BinaryTree<ndim>::StockCellProperties
(int output,
 SphParticle<ndim> *sphdata)
{
  int c,cc,ccc;                     // Cell counters
  int i;                            // Particle counter
  int k;                            // Dimension counter
  FLOAT dr[ndim];                   // Relative position vector
  FLOAT factor = 1.0/thetamaxsqd;   // ??
  FLOAT *crmax;                     // Max. extent of cell bounding boxes
  FLOAT *crmin;                     // Min. extent of cell bounding boxes

  debug2("[BinaryTree::StockCellProperties]");

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
  }

  for (c=0; c<Ncell; c++)
    for (k=0; k<ndim; k++) crmin[c*ndim + k] = big_number;

  for (c=0; c<Ncell; c++)
    for (k=0; k<ndim; k++) crmax[c*ndim + k] = -big_number;


  // Loop backwards over all tree cells to ensure child cells are always 
  // computed first before being summed in parent cells.
  // ==========================================================================
  for (c=Ncell-1; c>=0; c--) {

    // If this is a leaf cell, some over all particles
    // ------------------------------------------------------------------------
    if (tree[c].c2 == 0) {
      i = tree[c].ifirst;

      // Loop over all particles in cell summing their contributions
      while (i != -1) {
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
        i = inext[i];
      };

      // Normalise all cell values
      if (tree[c].N > 0) {
        for (k=0; k<ndim; k++) tree[c].r[k] /= tree[c].m;
        for (k=0; k<ndim; k++) dr[k] = crmax[c*ndim + k] - crmin[c*ndim + k];
        tree[c].cdistsqd = factor*DotProduct(dr,dr,ndim);
        for (k=0; k<ndim; k++) dr[k] = max(crmax[c*ndim + k] - tree[c].r[k],
		  			 tree[c].r[k] - crmin[c*ndim + k]);
        tree[c].rmax = sqrt(DotProduct(dr,dr,ndim));
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
        for (k=0; k<ndim; k++) dr[k] = crmax[c*ndim + k] - crmin[c*ndim + k];
        tree[c].cdistsqd = factor*DotProduct(dr,dr,ndim);
        for (k=0; k<ndim; k++) dr[k] = max(crmax[c*ndim + k] - tree[c].r[k],
                                           tree[c].r[k] - crmin[c*ndim + k]);
        tree[c].rmax = sqrt(DotProduct(dr,dr,ndim));
      }

    }
    // ------------------------------------------------------------------------

  }
  // ==========================================================================


#if defined(VERIFY_ALL)
  if (output == 1) {
    cout << "Root cell position1 : " << tree[0].r[0] << "   " 
	 << tree[0].r[1] << endl;
    cout << "Mass of root cell1 : " << tree[0].m << endl;
    cout << "Bounding box : " << crmin[0] << "   " << crmax[0] 
	 << "   " << crmin[1] << "   " << crmax[1] << endl;
  }

  // Go through each tree level in turn and print info
  if (output == 1) {
    for (int l=0; l<ltot+1; l++) {
      cout << "LEVEL : " << l << endl;
      cout << "----------------------" << endl;
      for (c=0; c<Ncell; c++) {
        if (tree[c].clevel == l) {
	      cout << "c : " << c << "   N : " 
		   << tree[c].N << "   " << Nleafmax << endl;
        }
      }
    }
  }
#endif


  // Free all locally allocated memory
  delete[] crmin;
  delete[] crmax;

  return;
}



//=============================================================================
//  BinaryTree::UpdateHmaxValues
/// Calculate the physical properties (e.g. total mass, centre-of-mass, 
/// opening-distance, etc..) of all cells in the tree.
//=============================================================================
template <int ndim>
void BinaryTree<ndim>::UpdateHmaxValues
(int output,
 SphParticle<ndim> *sphdata)
{
  int c,cc,ccc;                     // Cell counters
  int i;                            // Particle counter
  int k;                            // Dimension counter

  debug2("[BinaryTree::StockCellProperties]");

  // Zero all summation variables for all cells
  for (c=0; c<Ncell; c++) tree[c].hmax = 0.0;

  // Loop backwards over all tree cells to ensure child cells are always 
  // computed first before being summed in parent cells.
  // ==========================================================================
  for (c=Ncell-1; c>=0; c--) {

    // If this is a leaf cell, some over all particles
    // ------------------------------------------------------------------------
    if (tree[c].c2 == 0) {
      i = tree[c].ifirst;

      // Loop over all particles in cell summing their contributions
      while (i != -1) {
        tree[c].hmax = max(tree[c].hmax,sphdata[i].h);
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

  return;
}



//=============================================================================
//  BinaryTree::ComputeActiveCellList
/// Returns the number of cells containing active particles, 'Nactive', and
/// the i.d. list of cells contains active particles, 'celllist'
//=============================================================================
template <int ndim>
int BinaryTree<ndim>::ComputeActiveCellList(int *celllist)
{
  int c;                           // ..
  int Nactive = 0;                 // No. of cells containing active ptcls

  debug2("[BinaryTree::ComputeActiveCellList]");

  for (c=0; c<Ncell; c++) {
    if (tree[c].Nactive > 0) celllist[Nactive++] = c;
    //cout << "Nactive : " << c << "    " << tree[c].Nactive << "   " << Nactive << endl;
  }

  return Nactive;
}



//=============================================================================
//  GridSearch::ComputeActiveParticleList
/// Returns the number (Nactive) and list of ids (activelist) of all active
/// SPH particles in the given cell 'c'.
//=============================================================================
template <int ndim>
int BinaryTree<ndim>::ComputeActiveParticleList
(int c,                             ///< [in] Cell i.d.
 int *activelist,                   ///< [out] List of active particles in cell
 Sph<ndim> *sph)                    ///< [in] SPH object pointer
{
  int Nactive = 0;                  // No. of active particles in cell
  int i = tree[c].ifirst;           // Particle id (set to first ptcl id)
  int ilast = tree[c].ilast;        // i.d. of last particle in cell c

  // Else walk through linked list to obtain list and number of active ptcls.
  while (i != -1) {
    if (i < sph->Nsph && sph->sphdata[i].active) activelist[Nactive++] = i;
    if (i == ilast) break;
    i = inext[i];
  };

  return Nactive;
}



//=============================================================================
//  BinaryTree::ComputeGatherNeighbourList
/// Computes and returns number of neighbour, 'Nneib', and the list
/// of neighbour ids, 'neiblist', for all particles inside cell 'c'.
/// Includes all particles in the selected cell, plus all particles
/// contained in adjacent cells (including diagonal cells).
//=============================================================================
template <int ndim>
int BinaryTree<ndim>::ComputeGatherNeighbourList
(int c,                             ///< [in] i.d. of cell
 int Nneibmax,                      ///< [in] Max. no. of neighbours
 int *neiblist,                     ///< [out] List of neighbour i.d.s
 FLOAT hmax)                        ///< [in] Maximum smoothing length
{
  int cc;                           // Cell counter
  int i;                            // Particle id
  int ilast;                        // id of last particle in current cell
  int k;                            // Neighbour counter
  int Nneib = 0;                    // No. of neighbours
  FLOAT dr[ndim];                   // ..
  FLOAT drsqd;                      // ..
  FLOAT rc[ndim];                   // ..
  FLOAT hrangemax;                  // ..
  FLOAT kernrange = 3.0;            // ..

  for (k=0; k<ndim; k++) rc[k] = tree[c].r[k];
  hrangemax = tree[c].rmax + kernrange*hmax;

  // Start with root cell and walk through entire tree
  cc = 0;

  // ==========================================================================
  while (cc < Ncell) {
    for (k=0; k<ndim; k++) dr[k] = tree[cc].r[k] - rc[k];
    drsqd = DotProduct(dr,dr,ndim);

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

  return Nneib;
}



//=============================================================================
//  BinaryTree::ComputeNeighbourList
/// Computes and returns number of neighbour, 'Nneib', and the list
/// of neighbour ids, 'neiblist', for all particles inside cell 'c'.
/// Includes all particles in the selected cell, plus all particles
/// contained in adjacent cells (including diagonal cells).
//=============================================================================
template <int ndim>
int BinaryTree<ndim>::ComputeNeighbourList
(int c,                             ///< [in] i.d. of cell
 int Nneibmax,                      ///< [in] Max. no. of neighbours
 int *neiblist)                     ///< [out] List of neighbour i.d.s
{
  int cc;                           // Cell counter
  int i;                            // Particle id
  int ilast;                        // id of last particle in current cell
  int k;                            // Neighbour counter
  int Nneib = 0;                    // No. of neighbours
  FLOAT dr[ndim];                   // ..
  FLOAT drsqd;                      // ..
  FLOAT rc[ndim];                   // ..
  FLOAT hrangemax;                  // ..
  FLOAT kernrange = 3.0;            // ..
  FLOAT rmax;                       // ..

  for (k=0; k<ndim; k++) rc[k] = tree[c].r[k];
  hrangemax = tree[c].rmax + kernrange*tree[c].hmax;
  rmax = tree[c].rmax;

  // Start with root cell and walk through entire tree
  cc = 0;

  // ==========================================================================
  while (cc < Ncell) {
    for (k=0; k<ndim; k++) dr[k] = tree[cc].r[k] - rc[k];
    drsqd = DotProduct(dr,dr,ndim);

    // Check if circular range overlaps cell bounding sphere
    // ------------------------------------------------------------------------
    if (drsqd < pow(tree[cc].rmax + hrangemax,2) ||
        drsqd < pow(tree[cc].rmax + rmax + kernrange*tree[cc].hmax,2)) {

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

  return Nneib;
}



//=============================================================================
//  BinaryTree::ComputeGravityInteractionList
/// ..
//=============================================================================
template <int ndim>
int BinaryTree<ndim>::ComputeGravityInteractionList
(int c,                             ///< [in] i.d. of cell
 int Nneibmax,                      ///< [in] Max. no. of neighbours
 int Ndirectmax,                    ///< [in] ..
 int Ngravcellmax,                  ///< [in] ..
 int &Nneib,                        ///< [out] ..
 int &Ndirect,                      ///< [out] ..
 int &Ngravcell,                    ///< [out] ..
 int *neiblist,                     ///< [out] List of neighbour i.d.s
 int *directlist,                   ///< [out] ..
 int *gravcelllist)                 ///< [out] ..
{
  int cc;                           // Cell counter
  int i;                            // Particle id
  int ilast;                        // id of last particle in current cell
  int k;                            // Neighbour counter
  FLOAT dr[ndim];                   // ..
  FLOAT drsqd;                      // ..
  FLOAT rc[ndim];                   // ..
  FLOAT hrangemax;                  // ..
  FLOAT kernrange = 2.0;            // ..
  FLOAT rmax;                       // ..

  for (k=0; k<ndim; k++) rc[k] = tree[c].r[k];
  hrangemax = tree[c].rmax + kernrange*tree[c].hmax;
  rmax = tree[c].rmax;

  // Start with root cell and walk through entire tree
  cc = 0;
  Nneib = 0;
  Ndirect = 0;
  Ngravcell = 0;

  // ==========================================================================
  while (cc < Ncell) {
    for (k=0; k<ndim; k++) dr[k] = tree[cc].r[k] - rc[k];
    drsqd = DotProduct(dr,dr,ndim);

    // Check if cells contain SPH neighbours
    // ------------------------------------------------------------------------
    if (drsqd < pow(tree[cc].rmax + hrangemax,2) ||
        drsqd < pow(tree[cc].rmax + rmax + kernrange*tree[cc].hmax,2)) {

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
    // ------------------------------------------------------------------------
    else if (drsqd >= tree[c].cdistsqd + tree[cc].cdistsqd) {

      // If cell is a leaf-cell with only one particles, more efficient to
      // compute the particle than the cell
      if (tree[cc].c2 == 0 && tree[cc].N == 1 && Ndirect < Ndirectmax)
        directlist[Ndirect++] = tree[cc].ifirst;
      else if (Ngravcell < Ngravcellmax)
        gravcelllist[Ngravcell++] = cc;
      else
        return -1;
      cc = tree[cc].cnext;

    }
    // ------------------------------------------------------------------------
    else if (drsqd < tree[c].cdistsqd + tree[cc].cdistsqd) {

      // If not a leaf-cell, then open cell to first child cell
      if (tree[cc].c2 != 0)
         cc++;

      // If leaf-cell, add particles to list
      else if (tree[cc].c2 == 0 && Ndirect + Nleafmax < Ndirectmax) {
        i = tree[cc].ifirst;
        while (i != -1) {
          directlist[Ndirect++] = i;
       	  i = inext[i];
        };
        cc = tree[cc].cnext;
      }

      // If leaf-cell, but we've run out of memory, return with error-code (-1)
      else if (tree[cc].c2 == 0 && Ndirect + Nleafmax >= Ndirectmax)
       	return -1;

    }

    // If not in range, then open next cell
    // ------------------------------------------------------------------------
    else
      cc = tree[cc].cnext;
  };
  // ==========================================================================

  return 1;
}



//=============================================================================
//  BinaryTree::ComputeCellForces
/// ...
//=============================================================================
template <int ndim>
void BinaryTree<ndim>::ComputeCellForces
(int i,                             ///< [in] i.d. of particle
 int Ngravcell,                     ///< [in] ..
 int *gravcelllist,                 ///< [in] ..
 SphParticle<ndim> &parti)          ///< [inout] ..
{
  int c;                            // ..
  int cc;                           // ..
  int k;                            // ..
  FLOAT dr[ndim];                   // ..
  FLOAT drsqd;                      // ..
  FLOAT invdrmag;                   // ..

  // Loop over all neighbouring particles in list
  // --------------------------------------------------------------------------
  for (cc=0; cc<Ngravcell; cc++) {
    c = gravcelllist[cc];

    for (k=0; k<ndim; k++) dr[k] = tree[c].r[k] - parti.r[k];
    drsqd = DotProduct(dr,dr,ndim);
    invdrmag = 1.0/(sqrt(drsqd) + small_number);

    parti.gpot += tree[c].m*invdrmag;
    for (k=0; k<ndim; k++) 
      parti.agrav[k] += tree[c].m*dr[k]*pow(invdrmag,3);

  }
  // --------------------------------------------------------------------------

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
  int celldone;                   // ..
  int okflag;                     // ..
  int c;                           // Cell id
  int cc;                          // Aux. cell counter
  int g;                           // Leaf cell counter
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
  int *celllist;                   // List of active cells
  int *gatherlist;                 // ..
  int *neiblist;                   // List of neighbour ids
  FLOAT draux[ndim];               // Aux. relative position vector var
  FLOAT drsqdaux;                  // Distance squared
  FLOAT hrangesqd;                 // Kernel extent
  FLOAT hmax;                      // Maximum smoothing length
  FLOAT rp[ndim];                  // Local copy of particle position
  FLOAT *drsqd;                    // Position vectors to gather neibs
  FLOAT *m;                        // Distances to gather neibs
  FLOAT *m2;
  FLOAT *mu;                       // mass*u for gather neibs
  FLOAT *mu2;
  FLOAT *r;                        // Positions of neibs
  SphParticle<ndim> *data = sph->sphdata;  // Pointer to SPH particle data

  // Find list of all cells that contain active particles
  celllist = new int[gtot];
  cactive = ComputeActiveCellList(celllist);


  // Set-up all OMP threads
  // ==========================================================================
#pragma omp parallel default(shared) private(activelist,c,cc,celldone,draux,drsqd,drsqdaux,hmax,hrangesqd,i,j,jj,k,okflag,m,mu,Nactive,neiblist,Nneib,Nneibmax,r,rp)
  {
    Nneibmax = 2*sph->Ngather;
    activelist = new int[Nleafmax];
    gatherlist = new int[Nneibmax];
    neiblist = new int[Nneibmax];
    drsqd = new FLOAT[Nneibmax];
    m = new FLOAT[Nneibmax];
    m2 = new FLOAT[Nneibmax];
    mu = new FLOAT[Nneibmax];
    mu2 = new FLOAT[Nneibmax];
    r = new FLOAT[Nneibmax*ndim];

    // Loop over all active cells
    // ========================================================================
#pragma omp for schedule(dynamic)
    for (cc=0; cc<cactive; cc++) {
      c = celllist[cc];
      celldone = 1;
      hmax = tree[c].hmax;

      // If hmax is too small so the neighbour lists are invalid, make hmax
      // larger and then recompute for the current active cell.
      // ----------------------------------------------------------------------
      do {
        hmax = 1.05*hmax;

        // Find list of active particles in current cell
        Nactive = ComputeActiveParticleList(c,activelist,sph);

        // Compute neighbour list for cell depending on physics options
        Nneib = ComputeGatherNeighbourList(c,Nneibmax,neiblist,hmax);

        // If there are too many neighbours, reallocate the arrays and
        // recompute the neighbour lists.
        while (Nneib == -1) {
          delete[] r;
          delete[] mu2;
          delete[] mu;
          delete[] m2;
          delete[] m;
          delete[] drsqd;
          delete[] neiblist;
          delete[] gatherlist;
          Nneibmax = 2*Nneibmax;
          gatherlist = new int[Nneibmax];
          neiblist = new int[Nneibmax];
          drsqd = new FLOAT[Nneibmax];
          m = new FLOAT[Nneibmax];
          m2 = new FLOAT[Nneibmax];
          mu = new FLOAT[Nneibmax];
          mu2 = new FLOAT[Nneibmax];
          r = new FLOAT[Nneibmax*ndim];
          Nneib = ComputeGatherNeighbourList(c,Nneibmax,neiblist,hmax);
        };

        // Make local copies of important neib information (mass and position)
        for (jj=0; jj<Nneib; jj++) {
          j = neiblist[jj];
          m[jj] = data[j].m;
          mu[jj] = data[j].m*data[j].u;
          for (k=0; k<ndim; k++) r[ndim*jj + k] = (FLOAT) data[j].r[k];
        }

        // Loop over all active particles in the cell
        // --------------------------------------------------------------------
        for (j=0; j<Nactive; j++) {
          i = activelist[j];
          for (k=0; k<ndim; k++) rp[k] = data[i].r[k];

          // Set gather range as current h multiplied by some tolerance factor
          hrangesqd = pow(sph->kernp->kernrange*hmax,2);
          Ngather = 0;

          // Compute distance (squared) to all
          // ------------------------------------------------------------------
          for (jj=0; jj<Nneib; jj++) {
            for (k=0; k<ndim; k++) draux[k] = r[ndim*jj + k] - rp[k];
            drsqdaux = DotProduct(draux,draux,ndim);

            // Record distance squared for all potential gather neighbours
            if (drsqdaux <= hrangesqd) {
              gatherlist[Ngather] = jj;
              drsqd[Ngather] = drsqdaux;
              m2[Ngather] = m[jj];
              mu2[Ngather] = mu[jj];
              Ngather++;
            }

          }
          // ------------------------------------------------------------------

          // Validate that gather neighbour list is correct
#if defined(VERIFY_ALL)
          if (neibcheck) CheckValidNeighbourList(sph,i,Nneib,neiblist,"gather");
#endif

          // Compute smoothing length and other gather properties for ptcl i
          okflag = sph->ComputeH(i,Ngather,hmax,m2,mu2,drsqd,data[i],nbody);

          // If h-computation is invalid, then break from loop and recompute
          // larger neighbour lists
          if (okflag == 0) {
            celldone = 0;
            break;
          }

        }
        // --------------------------------------------------------------------

      } while (celldone == 0);
      // ----------------------------------------------------------------------

    }
    // ========================================================================

    // Free-up all memory
    delete[] r;
    delete[] mu2;
    delete[] mu;
    delete[] m2;
    delete[] m;
    delete[] drsqd;
    delete[] neiblist;
    delete[] gatherlist;
    delete[] activelist;

  }
  // ==========================================================================

  delete[] celllist;

  // Update all tree smoothing length values
  UpdateHmaxValues(0,sph->sphdata);

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
(Sph<ndim> *sph                     ///< Pointer to SPH object
)
{
  int c;                            // Cell id
  int cactive;                      // No. of active cells
  int cc;                           // Aux. cell counter
  int i;                            // Particle id
  int j;                            // Aux. particle counter
  int jj;                           // Aux. particle counter
  int k;                            // Dimension counter
  int okflag;                       // Flag if h-rho iteration is valid
  int Nactive;                      // No. of active particles in cell
  int Ngrav;                        // No. of direct sum gravity ptcls
  int Ninteract;                    // No. of near gather neighbours
  int Nneib;                        // No. of neighbours
  int Nneibmax;                     // Max. no. of neighbours
  int *activelist;                  // List of active particle ids
  int *celllist;                    // List of active cells
  int *interactlist;                // ..
  int *neiblist;                    // List of neighbour ids
  FLOAT draux[ndim];                // Aux. relative position vector var
  FLOAT drsqd;                      // Distance squared
  FLOAT hrangesqdi;                 // Kernel extent
  FLOAT hrangesqdj;                 // ..
  FLOAT rp[ndim];                   // Local copy of particle position
  FLOAT *dr;                        // Array of relative position vectors
  FLOAT *drmag;                     // Array of neighbour distances
  FLOAT *invdrmag;                  // Array of 1/drmag between particles
  SphParticle<ndim> *neibpart;      // Local copy of neighbouring ptcls
  SphParticle<ndim> parti;          // Local copy of SPH particle
  SphParticle<ndim> *data = sph->sphdata;   // Pointer to SPH particle data

  debug2("[BinaryTree::UpdateAllSphHydroForces]");


  // Find list of all cells that contain active particles
  celllist = new int[gtot];
  cactive = ComputeActiveCellList(celllist);
  Nneibmax = 2*sph->Ngather;


  // Set-up all OMP threads
  // ==========================================================================
#pragma omp parallel default(shared) private(activelist,c,cc,dr,draux,drmag,drsqd,hrangesqdi,hrangesqdj,i,interactlist,invdrmag,j,jj,k,okflag,Nactive,neiblist,neibpart,Ninteract,Nneib,parti,rp)
  {
    activelist = new int[Nleafmax];
    neiblist = new int[Nneibmax];
    interactlist = new int[Nneibmax];
    dr = new FLOAT[Nneibmax*ndim];
    drmag = new FLOAT[Nneibmax];
    invdrmag = new FLOAT[Nneibmax];
    neibpart = new SphParticle<ndim>[Nneibmax];

    // Loop over all active cells
    // ========================================================================
#pragma omp for schedule(dynamic)
    for (cc=0; cc<cactive; cc++) {
      c = celllist[cc];

      // Find list of active particles in current cell
      Nactive = ComputeActiveParticleList(c,activelist,sph);

      // Compute neighbour list for cell depending on physics options
      Nneib = ComputeNeighbourList(c,Nneibmax,neiblist);

      // If there are too many neighbours, reallocate the arrays and
      // recompute the neighbour lists.
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
        Nneib = ComputeNeighbourList(c,Nneibmax,neiblist);
      };

      // Make local copies of all potential neighbours
      for (j=0; j<Nneib; j++) {
        neibpart[j] = data[neiblist[j]];
        neibpart[j].div_v = (FLOAT) 0.0;
        neibpart[j].dudt = (FLOAT) 0.0;
        for (k=0; k<ndim; k++) neibpart[j].a[k] = (FLOAT) 0.0;
      }

      // Loop over all active particles in the cell
      // ----------------------------------------------------------------------
      for (j=0; j<Nactive; j++) {
        i = activelist[j];
        parti = data[i];
        parti.div_v = (FLOAT) 0.0;
        parti.dudt = (FLOAT) 0.0;
        for (k=0; k<ndim; k++) parti.a[k] = (FLOAT) 0.0;

        for (k=0; k<ndim; k++) rp[k] = parti.r[k]; //data[i].r[k];
        hrangesqdi = pow(sph->kernfac*sph->kernp->kernrange*parti.h,2);
        Ninteract = 0;

        // Validate that gather neighbour list is correct
#if defined(VERIFY_ALL)
        if (neibcheck) CheckValidNeighbourList(sph,i,Nneib,neiblist,"all");
#endif

        // Compute distances and the inverse between the current particle
        // and all neighbours here, for both gather and inactive scatter neibs.
        // Only consider particles with j > i to compute pair forces once
        // unless particle j is inactive.
        // --------------------------------------------------------------------
        for (jj=0; jj<Nneib; jj++) {

          hrangesqdj = pow(sph->kernfac*sph->kernp->kernrange*neibpart[jj].h,2);
          for (k=0; k<ndim; k++) draux[k] = neibpart[jj].r[k] - rp[k];
          drsqd = DotProduct(draux,draux,ndim);

          // Compute list of particle-neighbour interactions and also
          // compute ..
          if ((drsqd <= hrangesqdi || drsqd <= hrangesqdj) &&
              ((neiblist[jj] < i && !neibpart[jj].active) ||
              neiblist[jj] > i)) {
            interactlist[Ninteract] = jj;
            drmag[Ninteract] = sqrt(drsqd);
            invdrmag[Ninteract] = (FLOAT) 1.0/
              (drmag[Ninteract] + small_number);
            for (k=0; k<ndim; k++)
              dr[Ninteract*ndim + k] = draux[k]*invdrmag[Ninteract];
            Ninteract++;
          }

        }
        // --------------------------------------------------------------------

        // Compute all gather neighbour contributions to hydro forces
        sph->ComputeSphHydroForces(i,Ninteract,interactlist,
				   drmag,invdrmag,dr,parti,neibpart);

        // Add all summation variables to main arrays
        for (k=0; k<ndim; k++) {
#pragma omp atomic
          data[i].a[k] += parti.a[k];
        }
#pragma omp atomic
        data[i].dudt += parti.dudt;
#pragma omp atomic
        data[i].div_v += parti.div_v;

      }
      // ----------------------------------------------------------------------

      // Now add all active neighbour contributions to the main arrays
      for (jj=0; jj<Nneib; jj++) {
        if (neibpart[jj].active) {
          j = neiblist[jj];
          for (k=0; k<ndim; k++) {
#pragma omp atomic
            data[j].a[k] += neibpart[jj].a[k];
          }
#pragma omp atomic
          data[j].dudt += neibpart[jj].dudt;
#pragma omp atomic
          data[j].div_v += neibpart[jj].div_v;
        }
      }

    }
    // ========================================================================

    // Free-up local memory for OpenMP thread
    delete[] neibpart;
    delete[] invdrmag;
    delete[] drmag;
    delete[] dr;
    delete[] interactlist;
    delete[] neiblist;
    delete[] activelist;

  }
  // ==========================================================================

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
(Sph<ndim> *sph                     ///< Pointer to SPH object
)
{
  int c;                            // Cell id
  int cactive;                      // No. of active cells
  int cc;                           // Aux. cell counter
  int i;                            // Particle id
  int j;                            // Aux. particle counter
  int jj;                           // Aux. particle counter
  int k;                            // Dimension counter
  int okflag;                       // Flag if h-rho iteration is valid
  int Nactive;                      // No. of active particles in cell
  int Ngrav;                        // No. of direct sum gravity ptcls
  int Ndirect;                      // ..
  int Ndirectmax;                   // ..
  int Ngravcell;                    // ..
  int Ngravcellmax;                 // ..
  int Ninteract;                    // ..
  int Nneib;                        // No. of neighbours
  int Nneibmax;                     // Max. no. of neighbours
  int *activelist;                  // List of active particle ids
  int *celllist;                    // List of active cells
  int *directlist;                  // ..
  int *gravcelllist;                // ..
  int *interactlist;                // ..
  int *neiblist;                    // List of neighbour ids
  FLOAT draux[ndim];                // Aux. relative position vector var
  FLOAT drsqd;                      // Distance squared
  FLOAT hrangesqdi;                 // Kernel extent
  FLOAT hrangesqdj;                 // ..
  FLOAT rp[ndim];                   // Local copy of particle position
  FLOAT *agrav;                     // Local copy of gravitational accel.
  FLOAT *dr;                        // Array of relative position vectors
  FLOAT *drmag;                     // Array of neighbour distances
  FLOAT *gpot;                      // ..
  FLOAT *invdrmag;                  // Array of 1/drmag between particles
  SphParticle<ndim> *neibpart;      // Local copy of neighbouring ptcls
  SphParticle<ndim> parti;          // Local copy of SPH particle
  SphParticle<ndim> *data = sph->sphdata;   // Pointer to SPH particle data

  debug2("[BinaryTree::UpdateAllSphForces]");


  // Find list of all cells that contain active particles
  celllist = new int[gtot];
  cactive = ComputeActiveCellList(celllist);
  Nneibmax = 2*sph->Ngather;
  Ndirectmax = 2*Nneibmax;
  Ngravcellmax = 4*Nneibmax;


  // Set-up all OMP threads
  // ==========================================================================
#pragma omp parallel default(shared) private(activelist,agrav,c,cc,dr,draux,drmag,drsqd,hrangesqdi,hrangesqdj,i,interactlist,invdrmag,j,jj,k,okflag,Nactive,neiblist,neibpart,Ninteract,Nneib,parti,rp)
  {
    agrav = new FLOAT[ndim*sph->Nsph];
    gpot = new FLOAT[ndim*sph->Nsph];
    activelist = new int[Nleafmax];
    neiblist = new int[Nneibmax];
    interactlist = new int[Nneibmax];
    directlist = new int[Ndirectmax];
    gravcelllist = new int[Ngravcellmax];
    dr = new FLOAT[Nneibmax*ndim];
    drmag = new FLOAT[Nneibmax];
    invdrmag = new FLOAT[Nneibmax];
    neibpart = new SphParticle<ndim>[Nneibmax];

    // Zero temporary grav. accel array
    for (i=0; i<ndim*sph->Nsph; i++) agrav[i] = 0.0;

    // Loop over all active cells
    // ========================================================================
#pragma omp for schedule(dynamic)
    for (cc=0; cc<cactive; cc++) {
      c = celllist[cc];

      // Find list of active particles in current cell
      Nactive = ComputeActiveParticleList(c,activelist,sph);

      // Compute neighbour list for cell depending on physics options
      okflag = ComputeGravityInteractionList(c,Nneibmax,Ndirectmax,
					     Ngravcellmax,Nneib,Ndirect,
					     Ngravcell,neiblist,directlist,
					     gravcelllist);

      // If there are too many neighbours, reallocate the arrays and
      // recompute the neighbour lists.
      while (okflag == -1) {
        delete[] neibpart;
        delete[] invdrmag;
        delete[] drmag;
        delete[] dr;
        delete[] gravcelllist;
        delete[] interactlist;
        delete[] directlist;
        delete[] neiblist;
        Nneibmax = 2*Nneibmax;
        Ndirectmax = 2*Ndirectmax;
        Ngravcellmax = 2*Ngravcellmax;
        neiblist = new int[Nneibmax];
        interactlist = new int[Nneibmax];
        directlist = new int[Ndirectmax];
        gravcelllist = new int[Ngravcellmax];
        dr = new FLOAT[Nneibmax*ndim];
        drmag = new FLOAT[Nneibmax];
        invdrmag = new FLOAT[Nneibmax];
        neibpart = new SphParticle<ndim>[Nneibmax];
        okflag = ComputeGravityInteractionList(c,Nneibmax,Ndirectmax,
					       Ngravcellmax,Nneib,Ndirect,
					       Ngravcell,neiblist,directlist,
					       gravcelllist);
      };

      // Make local copies of all potential neighbours
      for (j=0; j<Nneib; j++) {
        neibpart[j] = data[neiblist[j]];
        neibpart[j].div_v = (FLOAT) 0.0;
        neibpart[j].dudt = (FLOAT) 0.0;
        neibpart[j].gpot = (FLOAT) 0.0;
        for (k=0; k<ndim; k++) neibpart[j].a[k] = (FLOAT) 0.0;
        for (k=0; k<ndim; k++) neibpart[j].agrav[k] = (FLOAT) 0.0;
      }

      // Loop over all active particles in the cell
      // ----------------------------------------------------------------------
      for (j=0; j<Nactive; j++) {
        i = activelist[j];
        parti = data[i];
        parti.div_v = (FLOAT) 0.0;
        parti.dudt = (FLOAT) 0.0;
        parti.gpot = (FLOAT) 0.0;
        for (k=0; k<ndim; k++) parti.a[k] = (FLOAT) 0.0;
        for (k=0; k<ndim; k++) parti.agrav[k] = (FLOAT) 0.0;

        // Determine interaction list (to ensure we don't compute pair-wise
        // forces twice)
        Ninteract = 0;
        for (jj=0; jj<Nneib; jj++) {
          if ((neiblist[jj] < i && !neibpart[jj].active) || neiblist[jj] >= i) {
            interactlist[Ninteract++] = jj;
          }
        }

        // Compute forces between SPH neighbours (hydro and gravity)
        sph->ComputeSphHydroGravForces(i,Ninteract,interactlist,
				       parti,neibpart);

        // Compute direct gravity forces between distant particles
        sph->ComputeDirectGravForces(i,Ndirect,directlist,
				     agrav,gpot,parti,data);

        // Compute gravitational force dues to distant cells
        ComputeCellForces(i,Ngravcell,gravcelllist,parti);

        // Add summed variables (acceleration, etc..) to main arrays
        for (k=0; k<ndim; k++) {
#pragma omp atomic
          data[i].a[k] += parti.a[k];
#pragma omp atomic
          data[i].agrav[k] += parti.agrav[k];
        }
#pragma omp atomic
        data[i].gpot += parti.gpot;
#pragma omp atomic
        data[i].dudt += parti.dudt;
#pragma omp atomic
        data[i].div_v += parti.div_v;

      }
      // ----------------------------------------------------------------------

      // Now add all active neighbour contributions to the main arrays
      for (jj=0; jj<Nneib; jj++) {
        if (neibpart[jj].active) {
          j = neiblist[jj];
          for (k=0; k<ndim; k++) {
#pragma omp atomic
            data[j].a[k] += neibpart[jj].a[k];
#pragma omp atomic
            data[j].agrav[k] += neibpart[jj].agrav[k];
          }
#pragma omp atomic
          data[j].dudt += neibpart[jj].dudt;
#pragma omp atomic
          data[j].div_v += neibpart[jj].div_v;
        }
      }

    }
    // ========================================================================


    // Now add all active ..
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
    delete[] invdrmag;
    delete[] drmag;
    delete[] dr;
    delete[] gravcelllist;
    delete[] directlist;
    delete[] interactlist;
    delete[] neiblist;
    delete[] activelist;
    delete[] gpot;
    delete[] agrav;

  }
  // ==========================================================================

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



//=============================================================================
//  BinaryTree::ValidateTree
/// Perform various consistency checks to ensure the tree structure and all 
/// cell values are valid.
//=============================================================================
template <int ndim>
void BinaryTree<ndim>::ValidateTree
(Sph<ndim> *sph)                    ///< [in] SPH object pointer
{
  bool treeflag;                    // ..
  int c;                            // .. 
  int cc;                           // ..
  int i;                            // ..
  int k;                            // ..
  int N;                            // ..
  FLOAT dr[ndim];                   // ..
  FLOAT drmag;                      // ..

  debug2("[BinaryTree::ValidateTree]");


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

    // Now loop over all child cells below cell to sum all properties
    // ------------------------------------------------------------------------
    while (cc < tree[c].cnext) {

      if (tree[cc].c2 == 0) {
    	i = tree[cc].ifirst;
    	while (i != -1) {
          for (k=0; k<ndim; k++) dr[k] = tree[c].r[k] - sph->sphdata[i].r[k];
          drmag = sqrt(DotProduct(dr,dr,ndim));
          if (drmag > 1.00001*tree[c].rmax) treeflag = false;
          if (sph->sphdata[i].h > 1.00001*tree[c].hmax) treeflag = false;
          if (!treeflag) {
    	    cout << "Problem with tree : " << c << "   " 
		 << cc << "   " << i << endl;
    	    cout << "rc : " << tree[c].r[0] << "   " << tree[c].r[1] << endl;
    	    cout << "hmax : " << tree[c].hmax << "   rmax : " 
		 << tree[c].rmax << endl;
    	    cout << "rp : " << sph->sphdata[i].r[0] << "    " 
		 << sph->sphdata[i].r[1] << endl;
    	    cout << "h : " << sph->sphdata[i].h << endl;
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

  }
  // --------------------------------------------------------------------------

  return;
}
#endif



template class BinaryTree<1>;
template class BinaryTree<2>;
template class BinaryTree<3>;
