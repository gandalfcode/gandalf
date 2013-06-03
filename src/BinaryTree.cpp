//=============================================================================
//  BinaryTree.cpp
//  Contains functions for grid neighbour search routines.
//  Creates a uniform grid from particle distribution where the spacing is 
//  the size of the maximum kernel range (i.e. kernrange*h_max) over all ptcls.
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
BinaryTree<ndim>::BinaryTree()
{
  allocated_tree = false;
  Ncell = 0;
  Ncellmax = 0;
  Ntot = 0;
  Ntotmax = 0;
  Ntotmaxold = 0;
  Nleafmax = 1;
}



//=============================================================================
//  BinaryTree::~BinaryTree
/// BinaryTree destructor.  Deallocates grid memory upon object destruction.
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
    grid = new struct GridCell[gtot];
    tree = new struct BinaryTreeCell<ndim>[Ncellmax];
    porder = new int* [ndim];
    for (int k=0; k<ndim; k++) porder[k] = new int[Ntotmax]; 
    r = new FLOAT* [ndim];
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
    for (int k=ndim; k>=0; k--) delete[] r[k];
    delete[] r;
    for (int k=ndim; k>=0; k--) delete[] porder[k];
    delete[] porder;
    delete[] tree;
    delete[] grid;
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
/// Creates a new grid structure each time the neighbour 'tree' needs to be 
/// updated.
//=============================================================================
template <int ndim>
void BinaryTree<ndim>::UpdateTree(Sph<ndim> *sph, Parameters &simparams)
{
  int output = 1;

  // Set number of tree members to total number of SPH particles (inc. ghosts)
  Ntotmaxold = Ntotmax;
  Nsph = sph->Nsph;
  Ntot = sph->Ntot;

  // Compute the size of all tree-related arrays now we know number of points
  ComputeTreeSize(output,sph->Nsph);

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

  exit(0);

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
  while (pow(2,ltot*Nleafmax) < Ntot) {
    ltot++;
  };
  gtot = pow(2,ltot);
  Ncell = 2*gtot - 1;
  Ncellmax = Ncell;
  Ntotmax = pow(2,ltot*Nleafmax);

  if (output == 1) {
    cout << "Calculating tree size variables" << endl;
    cout << "Max. no of particles in leaf-cell : " << Nleafmax << endl;
    cout << "Max. no of particles in tree : " << Ntotmax << endl;
    cout << "No. of particles in tree : " << Nsph << endl;
    cout << "No. of levels on tree : " << ltot << endl;
    cout << "No. of tree cells : " << Ncell << endl;
    cout << "No. of grid cells : " << gtot << endl;
  }

  return;
}



//=============================================================================
//  BinaryTree::CreateTreeStructure
/// Create the raw tree skeleton structure once the tree size is known.
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
    if (tree[c].clevel == ltot) {                   // If on leaf level then
      tree[c].cnext = c + 1;                        // id of next cell
      g++;                                          // Increment leaf cell no.
      tree[c].c2g = g;                              // Record leaf id
      g2c[g] = c;                                   // Record inverse id
    }
    else {
      tree[c+1].clevel = tree[c].clevel + 1;        // Level of 1st child
      tree[c].c2 = c + c2L[tree[c].clevel];         // id of 2nd child
      tree[tree[c].c2].clevel = tree[c].clevel + 1; // Level of 2nd child
      tree[c].cnext = c + cNL[tree[c].clevel];      // Next cell id
    }

#if defined(VERIFY_ALL)
    for (int cc=0; cc<c; cc++) { 
      if (tree[cc].clevel > ltot) {
	cout << "Problem with tree levels : " << cc << "   " << tree[cc].clevel
	     << "    " << ltot << endl;
	cout << "Other info : " << cc+1 << "   " << tree[cc].c2 
	     << "    " << tree[tree[cc].c2].clevel << endl;
	exit(0);
      }
    }
#endif
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
  debug2("[BinaryTree::LoadTree]");

  int c;                            // Cell counter
  int cc;                           // Secondary cell counter
  int k;                            // Dimensionality counter
  int i;                            // Particle counter
  int j;                            // Dummy particle id
  FLOAT *ccap;                      // Maximum capacity of cell
  FLOAT *ccon;                      // Current contents of cell

  // Allocate memory for local arrays
  ccap = new FLOAT[Ncellmax];
  ccon = new FLOAT[Ncellmax];

  // Set capacity of root-cell using particle weights
  ccap[0] = 0.0;
  for (i=0; i<Ntot; i++) pw[i] = 1.0/(FLOAT) Ntot;
  for (i=0; i<Ntot; i++) ccap[0] += pw[i];

  // Record capacity of 1st and 2nd child-cells as half that of parent cells
  for (c=0; c<Ncell; c++) {
    if (tree[c].c2 == 0) continue;
    ccap[c + 1] = 0.5*ccap[c];
    ccap[tree[c].c2] = ccap[c + 1];
  }

  // Initialise all particle and cell values before building tree structure
  for (i=0; i<Ntot; i++) pc[i] = 0;
  for (i=0; i<Ntot; i++) inext[i] = -1;
  for (c=0; c<Ncell; c++) ccon[c] = 0.0;
  for (c=0; c<Ncell; c++) ccap[c] = 0.50000000001*ccap[c];
  for (c=0; c<Ncell; c++) tree[c].ifirst = -1;
  for (c=0; i<Ncell; c++) tree[c].ilast = -1;
  c = 0;
  k = 0;


  // Loop through first cell on each level of the tree
  // --------------------------------------------------------------------------
  while (c < ltot) {

    // Loop over all particles (in order of current split)
    for (i=0; i<Ntot; i++) {
      j = porder[k][i];
      cc = pc[j];                           // Cell currently occupied by j
      ccon[cc] = ccon[cc] + pw[j];          // Add particle weighting to cell

      // If cell contains less than maximum allowed capacity, then add 
      // particle to the first child cell.  Otherwise add to second child.
      if (ccon[cc] < ccap[cc])
	pc[j]++;
      else
	pc[j] = tree[pc[j]].c2;
    }

    // Move to next level and cycle through each dimension in turn
    // (Need more sophisticated algorithm here in future)
    c++;
    k = (k + 1)%ndim;

  }
  // --------------------------------------------------------------------------

  // Loop over all particles and set id of first particle in each cell, plus 
  // the linked list values
  for (i=0; i<Ntot; i++) {
    if (tree[pc[i]].ifirst == 0)
      tree[pc[i]].ifirst = i;
    else
      inext[tree[pc[i]].ilast] = i;
    tree[pc[i]].ilast = i;
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
  FLOAT factor = 1.0/(theta*theta); // ??
  FLOAT *crmax;                     // Max. extent of tree cell bounding boxes
  FLOAT *crmin;                     // Min. extent of tree cell bounding boxes

  debug2("[BinaryTree::StockCellProperties]");

  // Allocate local memory
  crmax = new FLOAT[ndim*Ncellmax];
  crmin = new FLOAT[ndim*Ncellmax];

  // Zero all summation variables for all cells
  for (c=0; c<Ncell; c++) {
    tree[c].ifirst = -1;
    tree[c].ilast = -1;
    tree[c].Nactive = 0;
    tree[c].m = 0.0;
    tree[c].hmax = 0.0;
    for (k=0; k<ndim; k++) tree[c].r[k] = 0.0;
  }

  // Loop backwards over all tree cells to ensure child cells are always 
  // computed first before being summed in parant cells.
  // --------------------------------------------------------------------------
  for (c=Ncell-1; c>=0; c--) {

    // If this is a leaf cell, some over all particles
    // ------------------------------------------------------------------------
    if (tree[c].c2 == 0) {
      for (k=0; k<ndim; k++) crmin[c*ndim + k] = big_number;
      for (k=0; k<ndim; k++) crmax[c*ndim + k] = -big_number;
      i = tree[c].ifirst;

      // Loop over all particles in cell summing their contributions
      while (i != -1) {
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
      for (k=0; k<ndim; k++) tree[c].r[k] /= tree[c].m;
      for (k=0; k<ndim; k++) dr[k] = crmax[c*ndim + k] - crmin[c*ndim + k];
      tree[c].cdistsqd = factor*DotProduct(dr,dr,ndim);
      for (k=0; k<ndim; k++) dr[k] = max(crmax[c*ndim + k] - tree[c].r[k],
					 tree[c].r[k] - crmin[c*ndim + k]);
      tree[c].rmax = sqrt(DotProduct(dr,dr,ndim));

    }
    // For non-leaf cells, sum together two children cells
    // ------------------------------------------------------------------------
    else {
      cc = c + 1;
      ccc = tree[c].c2;
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
    // ------------------------------------------------------------------------

  }
  // --------------------------------------------------------------------------


#if defined(VERIFY_ALL)
  if (output == 1) {
    cout << "Root cell position1 : " << tree[0].r[0] << "   " 
	 << tree[0].r[1] << endl;
    cout << "Mass of root cell1 : " << tree[0].m << endl;
    cout << "Bounding box : " << crmin[0] << "   " << crmax[0] 
	 << "   " << crmin[1] << "   " << crmax[1] << endl;
  }
#endif

  // Free all locally allocated memory
  delete[] crmin;
  delete[] crmax;

  return;
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

  // If there are no active particles in this cell, return without walking list
  //if (tree[c].Nptcls == 0) return 0;

  // Else walk through linked list to obtain list and number of active ptcls.
  while (i != -1) {
    if (i < sph->Nsph && sph->sphdata[i].active) activelist[Nactive++] = i;
    if (i == ilast) break;
    i = inext[i];
  };

  return Nactive;
}



//=============================================================================
//  GridSearch::ComputeNeighbourList
/// Computes and returns number of neighbour, 'Nneib', and the list
/// of neighbour ids, 'neiblist', for all particles inside cell 'c'.
/// Includes all particles in the selected cell, plus all particles
/// contained in adjacent cells (including diagonal cells).
//=============================================================================
template <int ndim>
int GridSearch<ndim>::ComputeNeighbourList
(int c,                             ///< [in] i.d. of cell
 int *neiblist)                     ///< [out] List of neighbour i.d.s
{
  int i;                            // Particle id
  int ilast;                        // id of last particle in current cell
  int Nneib = 0;                    // No. of neighbours


  return Nneib;
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
(Sph<ndim> *sph                     ///< [inout] Pointer to main SPH object
 )
{
  int c;                           // Cell id
  int cc;                          // Aux. cell counter
  int g;                           // Leaf cell counter
  int cactive;                     // No. of active
  int i;                           // Particle id
  int j;                           // Aux. particle counter
  int jj;                          // Aux. particle counter
  int k;                           // Dimension counter
  int okflag;                      // Flag if h-rho iteration is valid
  int Nactive;                     // No. of active particles in cell
  int Ngather;                     // No. of near gather neighbours
  int Nneib;                       // No. of neighbours
  int Nneibmax;                    // Max. no. of neighbours
  int *activelist;                 // List of active particle ids
  int *celllist;                   // List of active cells
  int *neiblist;                   // List of neighbour ids
  FLOAT draux[ndim];               // Aux. relative position vector var
  FLOAT drsqdaux;                  // Distance squared
  FLOAT hrangesqd;                 // Kernel extent
  FLOAT rp[ndim];                  // Local copy of particle position
  FLOAT *drsqd;                    // Position vectors to gather neibs
  FLOAT *m;                        // Distances to gather neibs
  FLOAT *mu;                       // mass*u for gather neibs
  FLOAT *r;                        // Positions of neibs
  SphParticle<ndim> *data = sph->sphdata;  // Pointer to SPH particle data

  // Find list of all cells that contain active particles
  celllist = new int[gtot];
  //cactive = ComputeActiveCellList(celllist);

  // Set-up all OMP threads
  // ==========================================================================
#pragma omp parallel default(shared) private(activelist,c,cc,draux,drsqd,drsqdaux,hrangesqd,i,j,jj,k,okflag,m,mu,Nactive,neiblist,Nneib,r,rp)
  {
    activelist = new int[Nleafmax];
    neiblist = new int[Nneibmax];
    drsqd = new FLOAT[Nneibmax];
    m = new FLOAT[Nneibmax];
    mu = new FLOAT[Nneibmax];
    r = new FLOAT[Nneibmax*ndim];

    // Loop over all active cells
    // ========================================================================
#pragma omp for schedule(dynamic)
    for (cc=0; cc<cactive; cc++) {
      c = celllist[cc];

      // Find list of active particles in current cell
      Nactive = ComputeActiveParticleList(c,activelist,sph);

      // Compute neighbour list for cell depending on physics options
      //Nneib = ComputeGatherNeighbourList(c,neiblist);

      // Make local copies of important neib information (mass and position)
      for (jj=0; jj<Nneib; jj++) {
        j = neiblist[jj];
        m[jj] = data[j].m;
        mu[jj] = data[j].m*data[j].u;
        for (k=0; k<ndim; k++) r[ndim*jj + k] = (FLOAT) data[j].r[k];
      }

      // Loop over all active particles in the cell
      // ----------------------------------------------------------------------
      for (j=0; j<Nactive; j++) {
        i = activelist[j];
        for (k=0; k<ndim; k++) rp[k] = data[i].r[k];

        // Compute distance (squared) to all potential neighbours
        for (jj=0; jj<Nneib; jj++) {
          for (k=0; k<ndim; k++) draux[k] = r[ndim*jj + k] - rp[k];
	  drsqd[jj] = DotProduct(draux,draux,ndim);
        }

	// Validate that gather neighbour list is correct
#if defined(VERIFY_ALL)
	if (neibcheck) CheckValidNeighbourList(sph,i,Nneib,
					       neiblist,"gather");
#endif

        // Compute smoothing length and other gather properties for particle i
        okflag = sph->ComputeH(i,Nneib,m,mu,drsqd,data[i]);

      }
      // ----------------------------------------------------------------------

    }
    // ========================================================================

    // Free-up all memory
    delete[] r;
    delete[] mu;
    delete[] m;
    delete[] drsqd;
    delete[] neiblist;
    delete[] activelist;

  }
  // ==========================================================================

  delete[] celllist;

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
(Sph<ndim> *sph                     ///< Pointer to SPH object
)
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
(Sph *sph,                          ///< [in] SPH object pointer
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
      for (k=0; k<ndimmax; k++)
	dr[k] = sph->sphdata[j].r[k] - sph->sphdata[i].r[k];
      drsqd = DotProduct(dr,dr,ndim);
      if (drsqd <= 
	  sph->kernp->kernrangesqd*sph->sphdata[i].h*sph->sphdata[i].h)
	trueneiblist[Ntrueneib++] = j;
    }
  }

  // Now compare each given neighbour with true neighbour list for validation
  for (j=0; j<Ntrueneib; j++) {
    count = 0;
    for (k=0; k<Nneib; k++) if (neiblist[k] == trueneiblist[j]) count++;

    // If the true neighbour is not in the list, or included multiple times, 
    // then output to screen and terminate program
    if (count != 1) {
      cout << "Problem with neighbour lists : " << i << "  " << j << "   "
	   << count << "   "
	   << sph->sphdata[i].r[0] << "   " << sph->sphdata[i].h << endl;
      cout << "Nneib : " << Nneib << "   Ntrueneib : " << Ntrueneib << endl;
      PrintArray("neiblist     : ",Nneib,neiblist);
      PrintArray("trueneiblist : ",Ntrueneib,trueneiblist);
      string message = "Problem with neighbour lists in grid search";
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
