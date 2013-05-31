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

  Ntot = Npart;

  if (Ntot > Ntotmax || (!allocated_tree)) {
    if (allocated_tree) DeallocateTreeMemory();
    Ntotmax = Ntot;
    inext = new int[Ntotmax];
    pw = new FLOAT[Ntotmax];
    pc = new int[Ntotmax];
    g2c = new int[gtot];
    tree = new struct BinaryTreeCell<ndim>[Ncellmax];
    porder = new int* [ndim];
    r = new FLOAT* [ndim];
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
/// Creates a new grid structure each time the neighbour 'tree' needs to be 
/// updated.
//=============================================================================
template <int ndim>
void BinaryTree<ndim>::UpdateTree(Sph<ndim> *sph, Parameters &simparams)
{
  int output = 1;

  // Set number of tree members to total number of SPH particles (inc. ghosts)
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
  LoadTree(output);

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
  tree[0].cL = 0;

  // Loop over all cells and set all other pointers
  // --------------------------------------------------------------------------
  for (c=0; c<Ncell; c++) {
    if (tree[c].cL == ltot) {                 // If on leaf level then
      tree[c].cnext = c + 1;                  // ..
      g++;                                    // ..
      tree[c].c2g = g;                        // ..
      g2c[g] = c;                             // ..
      //cout << "Leaf cell : " << c << "   " << tree[c].cL << "   " << &tree[c].cnext << "   " << &tree[c].c2g << "   " << &g2c[g] << "    " << &tree[0].cL << endl;
      //cout << "g : " << g << "    " << gtot << endl;
    }
    else {
      tree[c+1].cL = tree[c].cL + 1;          // Compute level of 1st child
      tree[c].c2 = c + c2L[tree[c].cL];       // Compute id of 2nd child
      tree[tree[c].c2].cL = tree[c].cL + 1;   // Compute level of 2nd child
      tree[c].cnext = c + cNL[tree[c].cL];    // ..
      //cout << "Other cells : " << tree[c].cL << "   " << tree[c+1].cL << "    " << tree[c].c2 << "    " << tree[tree[c].c2].cL << "    " << tree[c].cnext << "   " << tree[0].cL << endl;
    }
    //cout << "Root tree levels : " << tree[0].cL << "    " << ltot << endl;
    for (int cc=0; cc<c; cc++) { 
      if (tree[cc].cL > ltot) {
	cout << "Problem with tree levels : " << cc << "   " << tree[cc].cL
	     << "    " << ltot << endl;
	cout << "Other info : " << cc+1 << "   " << tree[cc].c2 
	     << "    " << tree[tree[cc].c2].cL << endl;
	exit(0);
      }
    }
  }
  // --------------------------------------------------------------------------


  if (output == 1) {
    cout << "Diagnostic output from BinaryTree::CreateTreeStructure" << endl;
    for (c=0; c<min(20,Ncell-20); c++) {
      cout << "c : " << c << "   " << tree[c].cL << "   " << tree[c].c2
	   << "   " << tree[c].cnext << "   " << tree[c].c2g << "   "
	   << g2c[tree[c].c2g] << endl;
    }
    for (c=max(20,Ncell-20); c<Ncell; c++) {
      cout << "c : " << c << "   " << tree[c].cL << "   " << tree[c].c2
	   << "   " << tree[c].cnext << "   " << tree[c].c2g << "   "
	   << g2c[tree[c].c2g] << endl;
    }
  }

  return;
}



//=============================================================================
//  BinaryTree::OrderParticlesByCartCoord
/// ..
//=============================================================================
template <int ndim>
void BinaryTree<ndim>::OrderParticlesByCartCoord(int output,
						 SphParticle<ndim> *sphdata)
{
  int i,k;                          // Particle and dimension counters

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
  debug2("Ordering x-positions");
  Heapsort(output,Ntot,porder[0],r[0]);

  // Sort y-values
  if (ndim >= 2) debug2("Ordering y-positions");
  if (ndim >= 2) Heapsort(output,Ntot,porder[1],r[1]);

  // Sort z-values
  if (ndim == 3) debug2("Ordering z-positions");
  if (ndim == 3) Heapsort(output,Ntot,porder[2],r[2]);

  // Check that particles are ordered correctly
  for (i=0; i<Ntot; i++) {
    int j=porder[0][i];
    cout << "Order : " << i << "   " << j << "   " << sphdata[j].r[0] << endl;
  }
  for (k=0; k<ndim; k++) {
    cout << "Checking " << k << " dimenion ordering" << endl;
    for (int j=1; j<Ntot; j++) {
      int i1 = porder[k][j-1];
      int i2 = porder[k][j];
      if (sphdata[i2].r[k] < sphdata[i1].r[k]) {
	cout << "Problem with particle ordering : " << k << "   " << j << endl;
	exit(0);
      }
    }
  }

  return;
}



//=============================================================================
//  BinaryTree::LoadTree
/// ..
//=============================================================================
template <int ndim>
void BinaryTree<ndim>::LoadTree(int output)
{
  debug2("[BinaryTree::LoadTree]");

  int c;                            // Cell counter
  int cc;                           // Secondary cell counter
  int k;                            // Dimensionality counter
  int i;                            // Particle counter
  int j;                            // Dummy particle id
  int *ilast;                       // i.d. of last ptcl in leaf-cell thus far
  FLOAT *ccap;                      // Capacity of cell
  FLOAT *ccon;                      // Contents of cell

  // Allocate all local array
  ccap = new FLOAT[Ncellmax];
  ccon = new FLOAT[Ncellmax];
  ilast = new int[Ncellmax];
  //pc = new int[Ntot];

  // Set capacity of root-cell using particle weights
  ccap[0] = 0.0;
  for (i=0; i<Ntot; i++) pw[i] = 1.0/(FLOAT) Ntot;
  for (i=0; i<Ntot; i++) ccap[0] += pw[i];

  // Record capacity of 1st and 2nd child-cells
  for (c=0; c<Ncell; c++) {
    if (tree[c].c2 == 0) continue;
    ccap[c + 1] = 0.5*ccap[c];
    ccap[tree[c].c2] = ccap[c + 1];
  }


  for (i=0; i<Ntot; i++) pc[i] = 0;
  for (c=0; c<Ncell; c++) ccon[c] = 0.0;
  for (c=0; c<Ncell; c++) ccap[c] = 0.5000001*ccap[c];
  c = 0;
  k = 0;


  // --------------------------------------------------------------------------
  while (c < ltot) {

    // Loop over all particles (in order of current split)
    for (i=0; i<Ntot; i++) {
      cout << "Order2 : " << i << "   " << k << "   " << porder[k][i] << endl;
      j = porder[k][i];
      cc = pc[j];
      ccon[cc] = ccon[cc] + pw[j];
      if (ccon[cc] < ccap[cc])
	pc[j]++;
      else
	pc[j] = tree[pc[j]].c2;
    }
    c++;
    k = (k + 1)%ndim;
  }
  // --------------------------------------------------------------------------

  for (c=0; c<Ncell; c++) tree[c].cp = 0;
  for (i=0; i<Ntot; i++) inext[i] = -1;
  for (i=0; i<Ntot; i++) ilast[i] = -1;
  
  for (i=0; i<Ntot; i++) {
    if (tree[pc[i]].cp == 0)
      tree[pc[i]].cp = i;
    else
      inext[ilast[pc[i]]] = i;
    ilast[pc[i]] = i;
  }

  // Free all locally allocated memory
  delete[] ilast;
  delete[] ccon;
  delete[] ccap;

  return;
}



//=============================================================================
//  BinaryTree::StockCellProperties
/// ..
//=============================================================================
template <int ndim>
void BinaryTree<ndim>::StockCellProperties
(int output,
 SphParticle<ndim> *sphdata)
{
  int c,cc,ccc;
  int i;
  int k;
  FLOAT dr[ndim];                   
  FLOAT factor = 1.0/(theta*theta);
  FLOAT *crmax;
  FLOAT *crmin;
  FLOAT rroot[ndim];
  FLOAT mroot = 0.0;

  debug2("[BinaryTree::StockCellProperties]");

  // Allocate local memory
  crmax = new FLOAT[ndim*Ncellmax];
  crmin = new FLOAT[ndim*Ncellmax];

  // Zero all summation variables for all cells
  for (c=0; c<Ncell; c++) {
    tree[c].m = 0.0;
    for (k=0; k<ndim; k++) tree[c].r[k] = 0.0;
  }
  for (k=0; k<ndim; k++) rroot[k] = 0.0;

  // Loop backwards over all tree cells
  // --------------------------------------------------------------------------
  for (c=Ncell-1; c>=0; c--) {

    // If this is a leaf cell
    if (tree[c].c2 == 0) {
      for (k=0; k<ndim; k++) crmin[c*ndim + k] = big_number;
      for (k=0; k<ndim; k++) crmax[c*ndim + k] = -big_number;
      i = tree[c].cp;

      // Loop over all particles in cell summing their contributions
      while (i != -1) {
	tree[c].m += sphdata[i].m;
	for (k=0; k<ndim; k++) tree[c].r[k] += sphdata[i].m*sphdata[i].r[k];
	mroot += sphdata[i].m;
	for (k=0; k<ndim; k++) rroot[k] += sphdata[i].m*sphdata[i].r[k];
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
      tree[c].cd = factor*DotProduct(dr,dr,ndim);
      cout << "Leaf bbmin : " << crmin[c*ndim] << "   " << crmin[c*ndim+1] << endl;
    }
    else {
      cc = c + 1;
      ccc = tree[c].c2;
      tree[c].m = tree[cc].m + tree[ccc].m;
      for (k=0; k<ndim; k++) tree[c].r[k] = 
	(tree[cc].m*tree[cc].r[k] + tree[ccc].m*tree[ccc].r[k])/tree[c].m;
      for (k=0; k<ndim; k++)
	crmin[ndim*c + k] = min(crmin[ndim*cc+k],crmin[ndim*ccc+k]);
      for (k=0; k<ndim; k++)
	crmax[ndim*c + k] = max(crmax[ndim*cc+k],crmax[ndim*ccc+k]);
      for (k=0; k<ndim; k++) dr[k] = crmax[c*ndim + k] - crmin[c*ndim + k];
      tree[c].cd = factor*DotProduct(dr,dr,ndim);
      cout << "Cell bbmin : " << crmin[c*ndim] << "   " << crmin[c*ndim+1] << endl;
    }

  }
  // --------------------------------------------------------------------------

  for (k=0; k<ndim; k++) rroot[k] /= mroot;

  if (output == 1) {
    cout << "Root cell position1 : " << tree[0].r[0] << "   " << tree[0].r[1] << endl;
    cout << "Mass of root cell1 : " << tree[0].m << endl;
    cout << "Bounding box : " << crmin[0] << "   " << crmax[0] << "   " << crmin[1] << "   " << crmax[1] << endl;
    cout << "Root cell position2 : " << rroot[0] << "   " << rroot[1] << endl;
    cout << "Mass of root cell2 : " << mroot << endl;
  }

  // Free all locally allocated memory
  delete[] crmin;
  delete[] crmax;

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
(Sph<ndim> *sph                     ///< [inout] Pointer to main SPH object
 )
{
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
