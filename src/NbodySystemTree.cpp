//=============================================================================
//  NbodySystemTree.cpp
//  Contains all function definitions for constructing N-body nearest-neighbour
//  tree for creating systems and sub-systems.
//=============================================================================


#include <math.h>
#include <string>
#include "Precision.h"
#include "Constants.h"
#include "NbodyParticle.h"
#include "StarParticle.h"
#include "SystemParticle.h"
#include "NbodySystemTree.h"
#include "InlineFuncs.h"
#include "Debug.h"
using namespace std;


//=============================================================================
//  NbodySystemTree::NbodySystemTree()
/// NbodySystemTree constructor
//=============================================================================
template <int ndim>
NbodySystemTree<ndim>::NbodySystemTree()
{
  allocated_tree = false;
  Nnode = 0;
  Nnodemax = 0;
}



//=============================================================================
//  NbodySystemTree::~NbodySystemTree()
/// NbodySystemTree destructor
//=============================================================================
template <int ndim>
NbodySystemTree<ndim>::~NbodySystemTree()
{
  if (allocated_tree) DeallocateMemory();
}



//=============================================================================
//  NbodySystemTree::AllocateMemory
/// Allocate nearest-neighbour tree memory
//=============================================================================
template <int ndim>
void NbodySystemTree<ndim>::AllocateMemory(int N)
{
  debug2("[NbodySystemTree::AllocateMemory]");

  if (N > Nnodemax || !allocated_tree) {
    if (allocated_tree) DeallocateMemory();
    Nnodemax = N;
    NNtree = new NNTreeCell<ndim>[Nnodemax];
  }

  return;
}



//=============================================================================
//  NbodySystemTree::DeallocateMemory
/// Deallocate nearest-neighbour tree memory
//=============================================================================
template <int ndim>
void NbodySystemTree<ndim>::DeallocateMemory(void)
{
  debug2("[NbodySystemTree::DeallocateMemory]");

  if (allocated_tree) {
    delete[] NNtree;
  }
  allocated_tree = false;

  return;
}



//=============================================================================
//  NbodySystemTree::CreateNbodySystemTree
/// Creates a nearest neighbour tree based on the positions of all stars
/// contained in the nbody object that is passed.
//=============================================================================
template <int ndim>
void NbodySystemTree<ndim>::CreateNbodySystemTree
(Nbody<ndim> *nbody)               ///< [in] Nbody object containing stars
{
  int i,ii,j,jj;                   // Node ids and counters
  int k;                           // Dimension counter
  int Nfreenode;                   // No. of free (i.e. unattached) nodes
  int *nodelist;                   // List of unattached nodes
  DOUBLE dr[ndim];                 // Relative position vector
  DOUBLE drsqd;                    // Distance squared

  debug2("NbodySystemTree::CreateNbodySystemTree");

  Nnode = 0;
  Nfreenode = 0;

  // Allocate memory for tree
  AllocateMemory(2*nbody->Nstar);
  nodelist = new int[Nnodemax];

  // Initialise all node variables before adding stars
  for (i=0; i<Nnodemax; i++) {
    NNtree[i].ichild1 = -1;
    NNtree[i].ichild2 = -1;
    NNtree[i].iparent = -1;
    NNtree[i].inearest = -1;
    NNtree[i].rsqdnearest = big_number_dp;
    NNtree[i].gpe = 0.0;
    NNtree[i].gpe_internal = 0.0;
  }

  // Add one star to each lead node, recording the position and id
  for (i=0; i<nbody->Nstar; i++) {
    nodelist[i] = i;
    NNtree[i].Ncomp = 1;
    for (k=0; k<ndim; k++) NNtree[i].r[k] = nbody->stardata[i].r[k];
    Nnode++;
    Nfreenode++;
  }


  // Process all remaining unconnected nodes to find new set of mutually
  // nearest neighbours for next phase of tree construction
  // ==========================================================================
  while (Nnode < Nnodemax) {

	// Construct list of remaining nodes
    Nfreenode = 0;
    for (i=0; i<Nnode; i++) {
      if (NNtree[i].iparent == -1) nodelist[Nfreenode++] = i;
    }

    // If we have only one remaining unconnected node (i.e. the root) exit loop
    if (Nfreenode == 1) break;

    // Identify the nearest neighbouring node for each node
    // ------------------------------------------------------------------------
    for (ii=0; ii<Nfreenode; ii++) {
      i = nodelist[ii];
      NNtree[i].rsqdnearest = big_number;

      for (jj=0; jj<Nfreenode; jj++) {
        j = nodelist[jj];
        if (i == j) continue;

        for (k=0; k<ndim; k++) dr[k] = NNtree[i].r[k] - NNtree[j].r[k];
        drsqd = DotProduct(dr,dr,ndim);
        if (drsqd < NNtree[i].rsqdnearest) {
          NNtree[i].rsqdnearest = drsqd;
          NNtree[i].inearest = j;
        }
      }
    }
    // ------------------------------------------------------------------------


    // Now identify all mutually nearest neighbours to create a new generation
    // of nodes
    // ------------------------------------------------------------------------
    for (ii=0; ii<Nfreenode-1; ii++) {
      i = nodelist[ii];

      for (jj=ii+1; jj<Nfreenode; jj++) {
        j = nodelist[jj];

        // If each node is the others nearest neighbour, then create a new
        // parent node with the two original nodes as child nodes
        if (NNtree[i].inearest == j && NNtree[j].inearest == i) {
          NNtree[i].iparent = Nnode;
          NNtree[j].iparent = Nnode;
          NNtree[Nnode].ichild1 = i;
          NNtree[Nnode].ichild2 = j;
          for (k=0; k<ndim; k++)
            NNtree[Nnode].r[k] = 0.5*(NNtree[i].r[k] + NNtree[j].r[k]);
          NNtree[Nnode].Ncomp = NNtree[i].Ncomp + NNtree[j].Ncomp;
          Nnode++;
        }

      }

    }
    // ------------------------------------------------------------------------

  };
  // ==========================================================================

  // Free all locally allocated memory
  delete[] nodelist;

  return;
}



//=============================================================================
//  NbodySystemTree::StockNbodySystemTree
/// Calculate the properties of all nearest-neighbour tree nodes, starting from
/// single-star child-cells and then working up through the parent nodes and
/// finally the root node.
//=============================================================================
template <int ndim>
void NbodySystemTree<ndim>::StockNbodySystemTree
(Nbody<ndim> *nbody)               ///< [in] Nbody object containing stars
{
  int c,c1,c2;                     // Cell counter
  int i,ii,j,jj;                   // ..
  int k;                           // Dimension counter
  //DOUBLE dr[ndim];                 // ..
  //DOUBLE dv[ndim];                 // ..


  // Loop through all nodes of tree and compute all physical quantities
  // ==========================================================================
  for (c=0; c<Nnode; c++) {

	// If node contains one star, set all properties equal to star values
    // ------------------------------------------------------------------------
    if (NNtree[c].Ncomp == 1) {
      i = c;
      NNtree[c].m = nbody->stardata[i].m;
      //NNtree[c].h = nbody->stardata[i].h;
      NNtree[c].gpe = 0.5*nbody->stardata[i].m*nbody->stardata[i].gpot;
      NNtree[c].gpe_internal = 0.0;
      NNtree[c].tcross = big_number;
      for (k=0; k<ndim; k++) NNtree[c].r[k] = nbody->stardata[i].r[k];
      for (k=0; k<ndim; k++) NNtree[c].v[k] = nbody->stardata[i].v[k];
      for (k=0; k<ndim; k++) NNtree[c].a[k] = nbody->stardata[i].a[k];
      for (k=0; k<ndim; k++) NNtree[c].adot[k] = nbody->stardata[i].adot[k];

    }

	// Else, add together both child node properties
    // ------------------------------------------------------------------------
    else {
      c1 = NNtree[c].ichild1;
      c2 = NNtree[c].ichild2;
      NNtree[c].m = NNtree[c1].m + NNtree[c2].m;
      NNtree[c].gpe = NNtree[c1].gpe + NNtree[c2].gpe;
      NNtree[c].gpe_internal = 0.0;
      NNtree[c].tcross = big_number;
      for (k=0; k<ndim; k++) NNtree[c].r[k] = (NNtree[c1].m*NNtree[c1].r[k] +
    		  NNtree[c2].m*NNtree[c2].r[k])/NNtree[c].m;
      for (k=0; k<ndim; k++) NNtree[c].v[k] = (NNtree[c1].m*NNtree[c1].v[k] +
    		  NNtree[c2].m*NNtree[c2].v[k])/NNtree[c].m;
      for (k=0; k<ndim; k++) NNtree[c].a[k] = (NNtree[c1].m*NNtree[c1].a[k] +
    		  NNtree[c2].m*NNtree[c2].a[k])/NNtree[c].m;
      for (k=0; k<ndim; k++) NNtree[c].adot[k] = (NNtree[c1].m*NNtree[c1].adot[k] +
    		  NNtree[c2].m*NNtree[c2].adot[k])/NNtree[c].m;
    }
    // ------------------------------------------------------------------------


  }
  // ==========================================================================


  return;
}



//=============================================================================
//  NbodySystemTree::BuildSubSystems
//=============================================================================
template <int ndim>
void NbodySystemTree<ndim>::BuildSubSystems
(Nbody<ndim> *nbody)               ///< [in] Nbody object containing stars
{
  return;
}



//=============================================================================
//  NbodySystemTree::FindPerturberLists
//=============================================================================
template <int ndim>
void NbodySystemTree<ndim>::FindPerturberLists
(Nbody<ndim> *nbody)               ///< [in] Nbody object containing stars
{
  return;
}



template class NbodySystemTree<1>;
template class NbodySystemTree<2>;
template class NbodySystemTree<3>;
