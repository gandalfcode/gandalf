//=============================================================================
//  NbodySystemTree.cpp
//  Contains all function definitions for constructing N-body nearest-neighbour
//  tree for creating systems and sub-systems.
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


#include <math.h>
#include <string>
#include "Precision.h"
#include "Constants.h"
#include "NbodyParticle.h"
#include "StarParticle.h"
#include "SystemParticle.h"
#include "NbodySystemTree.h"
#include "MergeList.h"
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

  debug2("[NbodySystemTree::CreateNbodySystemTree]");

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
    NNtree[i].Ncomp = 0;
    NNtree[i].Nstar = 0;
    NNtree[i].Nchildlist = 0;
    NNtree[i].clist.clear();
  }

  // Add one star to each lead node, recording the position and id
  for (i=0; i<nbody->Nstar; i++) {
    nodelist[i] = i;
    NNtree[i].Ncomp = 1;
    NNtree[i].Nstar = 1;
    for (k=0; k<ndim; k++) NNtree[i].r[k] = nbody->stardata[i].r[k];
    Nnode++;
    Nfreenode++;
#if defined(VERIFY_ALL)
    cout << "Setting leaf nodes of tree : " << i << endl;
#endif
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
          NNtree[Nnode].Nstar = NNtree[i].Nstar + NNtree[j].Nstar;
	  NNtree[Nnode].Ncomp = NNtree[i].Ncomp + NNtree[j].Ncomp;
          Nnode++;
#if defined(VERIFY_ALL)
	  cout << "Adding new node to tree : " << Nnode - 1 << "   " 
	       << i << "    " << j << endl;
#endif
        }

      }

    }
    // ------------------------------------------------------------------------

  };
  // ==========================================================================

  // Free all locally allocated memory
  delete[] nodelist;

#if defined(VERIFY_ALL)
  cout << "Nnode : " << Nnode << endl;
#endif

  return;
}



//=============================================================================
//  NbodySystemTree::BuildSubSystems
/// Calculate the properties of all nearest-neighbour tree nodes, starting from
/// single-star child-cells and then working up through the parent nodes and
/// finally the root node.
//=============================================================================
template <int ndim>
void NbodySystemTree<ndim>::BuildSubSystems
(Nbody<ndim> *nbody)               ///< [inout] Nbody object containing stars
{
  int c,c1,c2;                     // Cell counter
  int i,ii,j,jj;                   // Star and system particle counters
  int k;                           // Dimension counter
  int Nsystem;                     // No. of systems
  DOUBLE dr[ndim];                 // Relative position vector
  DOUBLE drsqd;                    // Distance squared
  DOUBLE dv[ndim];                 // Relative velocity vector
  DOUBLE invdrmag;                 // 1 / drmag
  DOUBLE ketot = 0.0;              // Kinetic energy
  DOUBLE vmean;                    // Average speed of stars in sub-cluster
  NbodyParticle<ndim> *si;         // ..
  NbodyParticle<ndim> *sj;         // ..

  debug2("[NbodySystemTree::BuildSubSystems]");

  // Set all counters
  nbody->Nnbody = 0;
  nbody->Nsystem = 0;
  Nsystem = 0;

  // Loop through all nodes of tree and compute all physical quantities
  // ==========================================================================
  for (c=0; c<Nnode; c++) {

#if defined(VERIFY_ALL)
    cout << "Stocking node : " << c << "    Ncomp : " 
	 << NNtree[c].Ncomp << endl;
#endif

    // If node contains one star, set all properties equal to star values
    // ------------------------------------------------------------------------
    if (NNtree[c].Ncomp == 1) {
      i = c;
      NNtree[c].m = nbody->stardata[i].m;
      NNtree[c].h = nbody->stardata[i].h;
      for (k=0; k<ndim; k++) NNtree[c].r[k] = nbody->stardata[i].r[k];
      for (k=0; k<ndim; k++) NNtree[c].v[k] = nbody->stardata[i].v[k];
      for (k=0; k<ndim; k++) NNtree[c].a[k] = nbody->stardata[i].a[k];
      for (k=0; k<ndim; k++) NNtree[c].adot[k] = nbody->stardata[i].adot[k];
      for (k=0; k<ndim; k++) NNtree[c].a2dot[k] = nbody->stardata[i].a2dot[k];
      for (k=0; k<ndim; k++) NNtree[c].a3dot[k] = nbody->stardata[i].a3dot[k];
      NNtree[c].gpe = 0.5*nbody->stardata[i].m*nbody->stardata[i].gpot;
      NNtree[c].gpe_internal = 0.0;
      NNtree[c].tcross = big_number;
      NNtree[c].Nchildlist = 1;
      NNtree[c].childlist[0] = &nbody->stardata[i];
      //NNtree[c].clist.push_back(&(nbody->stardata[i]));

#if defined(VERIFY_ALL)
      cout << "Stocking single star data : " << i << "    " << NNtree[c].r[0] 
	   << "    " << NNtree[c].r[1] << "    " << NNtree[c].gpe << endl;
#endif

    }

    // Else, add together both child node properties
    // ------------------------------------------------------------------------
    else {

      c1 = NNtree[c].ichild1;
      c2 = NNtree[c].ichild2;
      NNtree[c].Ncomp = NNtree[c1].Ncomp + NNtree[c2].Ncomp;
      NNtree[c].m = NNtree[c1].m + NNtree[c2].m;
      NNtree[c].h = max(NNtree[c1].h,NNtree[c2].h);
      NNtree[c].gpe = NNtree[c1].gpe + NNtree[c2].gpe; 
      //- NNtree[c1].gpe_internal - NNtree[c2].gpe_internal;
      for (k=0; k<ndim; k++) {
        NNtree[c].r[k] = (NNtree[c1].m*NNtree[c1].r[k] +
			  NNtree[c2].m*NNtree[c2].r[k])/NNtree[c].m;
        NNtree[c].v[k] = (NNtree[c1].m*NNtree[c1].v[k] +
			  NNtree[c2].m*NNtree[c2].v[k])/NNtree[c].m;
        NNtree[c].a[k] = (NNtree[c1].m*NNtree[c1].a[k] +
			  NNtree[c2].m*NNtree[c2].a[k])/NNtree[c].m;
        NNtree[c].adot[k] = (NNtree[c1].m*NNtree[c1].adot[k] +
			  NNtree[c2].m*NNtree[c2].adot[k])/NNtree[c].m;
        NNtree[c].a2dot[k] = (NNtree[c1].m*NNtree[c1].a2dot[k] +
			  NNtree[c2].m*NNtree[c2].a2dot[k])/NNtree[c].m;
        NNtree[c].a3dot[k] = (NNtree[c1].m*NNtree[c1].a3dot[k] +
			  NNtree[c2].m*NNtree[c2].a3dot[k])/NNtree[c].m;
      }
      NNtree[c].gpe_internal = 0.0;
      ketot = 0.0;

#if defined(VERIFY_ALL)
      cout << "Stocking system data : " << c << "    " << NNtree[c].r[0] 
	   << "    " << NNtree[c].r[1] << endl;
      //cout << "Child data : " << c1 << "    " << c2 << "    " 
      //   << NNtree[c1].gpe << "    " << NNtree[c2].gpe << endl;
      cout << "Child energies : " << NNtree[c].gpe << "    " 
	   << NNtree[c1].gpe << "   " << NNtree[c2].gpe
	   << "     " << NNtree[c].gpe_internal << "     " 
	   << NNtree[c1].gpe_internal << "    " 
	   << NNtree[c2].gpe_internal << endl;
#endif


      // If node contains within maximum allowed number of components, then 
      // check if node is a new system or not
      // ----------------------------------------------------------------------
      if (NNtree[c].Ncomp <= Ncompmax) {

	// Merge lists of stars/systems in child nodes
	NNtree[c].Nchildlist = NNtree[c1].Nchildlist + NNtree[c2].Nchildlist;
	for (i=0; i<NNtree[c1].Nchildlist; i++)
	  NNtree[c].childlist[i] = NNtree[c1].childlist[i];
	for (i=0; i<NNtree[c2].Nchildlist; i++)
	  NNtree[c].childlist[i+NNtree[c1].Nchildlist] = NNtree[c2].childlist[i];
	NNtree[c].clist = NNtree[c1].clist + NNtree[c2].clist;

#if defined(VERIFY_ALL)
	cout << "No. of child systems for " << c << "   :   " 
	     << NNtree[c].Nchildlist << endl;
#endif
	
	// Compute internal kinetic energy
	ketot = 0.0;
	for (i=0; i<NNtree[c].Nchildlist; i++) {
	  si = NNtree[c].childlist[i];
	  for (k=0; k<ndim; k++) dv[k] = si->v[k] - NNtree[c].v[k];
	  ketot += 0.5*si->m*DotProduct(dv,dv,ndim);
	}
	vmean = sqrt(2.0*ketot/NNtree[c].m);
	
	// Compute internal gravitational potential energy
	for (i=0; i<NNtree[c].Nchildlist-1; i++) {
	  si = NNtree[c].childlist[i];
	  for (j=i+1; j<NNtree[c].Nchildlist; j++) {
	    sj = NNtree[c].childlist[j];
	    for (k=0; k<ndim; k++) dr[k] = sj->r[k] - si->r[k];
	    drsqd = DotProduct(dr,dr,ndim);
	    invdrmag = 1.0 / (sqrt(drsqd) + small_number_dp);
	    NNtree[c].gpe_internal += si->m*sj->m*invdrmag;
	  }
	}

	// Then, estimate the crossing time as tcross = Rgrav / vdisp
	// where Rgrav = sqrt(m^2/G/Egrav).  Should give similar
	// timescale to binary period in case of bound binary.
	NNtree[c].tcross = 
	  sqrt(NNtree[c].m*NNtree[c].m/NNtree[c].gpe_internal)/vmean;
	
#if defined(VERIFY_ALL)
	cout << "gpe : " << NNtree[c].gpe << "    gpe_internal : "
	     << NNtree[c].gpe_internal << "    ketot : " << ketot << endl;
	if (NNtree[c].gpe < NNtree[c].gpe_internal) {
	  cout << "Grav. energy problem" << endl;
	  exit(0);
	}
#endif
	
	// Compute and store binary properties if bound
	if (NNtree[c].Ncomp == 2) {
	  for (k=0; k<ndim; k++) dv[k] = si->v[k] - sj->v[k];
	  for (k=0; k<ndim; k++) dr[k] = sj->r[k] - si->r[k];
	  drsqd = DotProduct(dr,dr,ndim);
	  invdrmag = 1.0 / (sqrt(drsqd) + small_number_dp);
#if defined(VERIFY_ALL)
	  cout << "sma : " << -0.5*NNtree[c].m/
	    (0.5*DotProduct(dv,dv,ndim) - NNtree[c].m*invdrmag) << endl;
#endif
	}


#if defined(VERIFY_ALL)
	cout << "Checking system criteria : " 
	   << fabs(NNtree[c].gpe - NNtree[c].gpe_internal)/NNtree[c].gpe 
	   << "    " << gpefrac << "    " << NNtree[c].Ncomp << "     " 
	   << Ncompmax << endl;
#endif


	// Now check energies and decide if we should create a new sub-system
	// object.  If yes, create new system in main arrays
	// --------------------------------------------------------------------
	if (fabs(NNtree[c].gpe - NNtree[c].gpe_internal) 
	    < gpefrac*NNtree[c].gpe) {
	  
	  // Copy centre-of-mass properties of new sub-system
	  nbody->system[Nsystem].inode = c;
	  nbody->system[Nsystem].dt_internal = NNtree[c].tcross;
	  nbody->system[Nsystem].m = NNtree[c].m;
          nbody->system[Nsystem].h = NNtree[c].h;
	  for (k=0; k<ndim; k++) {
	    nbody->system[Nsystem].r[k] = NNtree[c].r[k];
	    nbody->system[Nsystem].v[k] = NNtree[c].v[k];
	    nbody->system[Nsystem].a[k] = NNtree[c].a[k];
	    nbody->system[Nsystem].adot[k] = NNtree[c].adot[k];
	    nbody->system[Nsystem].a2dot[k] = NNtree[c].a2dot[k];
	    nbody->system[Nsystem].a3dot[k] = NNtree[c].a3dot[k];
	    nbody->system[Nsystem].r0[k] = NNtree[c].r[k];
	    nbody->system[Nsystem].v0[k] = NNtree[c].v[k];
	    nbody->system[Nsystem].a0[k] = NNtree[c].a[k];
	    nbody->system[Nsystem].adot0[k] = NNtree[c].adot[k];
	  }
	  
	  // Copy list of contained N-body particles (either systems or stars)
	  // to system object, and also total no. of stars/components
	  nbody->system[Nsystem].Ncomp = NNtree[c].Ncomp;
	  nbody->system[Nsystem].Nchildren = 0;
	  nbody->system[Nsystem].Npert = 0;
	  for (i=0; i<NNtree[c].Nchildlist; i++)
	    nbody->system[Nsystem].children[nbody->system[Nsystem].Nchildren++]
	      = NNtree[c].childlist[i];

#if defined(VERIFY_ALL)
	  cout << "Adding binary system : " << c << "    " << c1 << "    " 
	       << c2 << "     " << nbody->system[Nsystem].Ncomp << endl;
	  cout << "System no. : " << Nsystem << endl;
	  cout << "r : " << NNtree[c].r[0] << "    " << NNtree[c].r[1] << endl;
#endif
	  //cout << "v : " << NNtree[c].v[0] << "    " << NNtree[c].v[1] << endl;
	  //cout << "m : " << NNtree[c].m << endl;
	  //cout << "dt_internal : " << NNtree[c].tcross << endl;

	  // Finally, clear list and append newly created system particle.
	  // Also, zero internal energy to allow detection of hierarchical 
	  // systems.
          NNtree[c].Ncomp = 1;
	  NNtree[c].Nchildlist = 1;
	  NNtree[c].childlist[0] = &(nbody->system[Nsystem]);
	  NNtree[c].gpe = NNtree[c].gpe - NNtree[c].gpe_internal;
	  NNtree[c].gpe_internal = 0.0;
	  //NNtree[c].clist.clear();
	  //NNtree[c].clist.push_back(&(nbody->system[Nsystem]));
	  Nsystem++;

	}

      }

      
      // If lists contain more than the maximum allowed number of components,
      // then flush list of systems to main arrays
      // ----------------------------------------------------------------------
      else {

	// Merge lists of stars/systems in child nodes
	for (i=0; i<NNtree[c1].Nchildlist; i++)
	  nbody->nbodydata[nbody->Nnbody++] = NNtree[c1].childlist[i];
	for (i=0; i<NNtree[c2].Nchildlist; i++)
	  nbody->nbodydata[nbody->Nnbody++] = NNtree[c2].childlist[i];
	NNtree[c].Nchildlist = 0;
	NNtree[c1].Nchildlist = 0;
	NNtree[c2].Nchildlist = 0;

      }
      // ----------------------------------------------------------------------

    }
    // ------------------------------------------------------------------------

  }
  // ==========================================================================

  
  // For root cell, push all remaining systems (or single sub-systems) to list
  c = Nnode - 1;
#if defined(VERIFY_ALL)
  cout << "No. of remaining systems in root node " << c << "  :  " 
       << NNtree[c].Nchildlist << endl;
#endif
  for (i=0; i<NNtree[c].Nchildlist; i++)
    nbody->nbodydata[nbody->Nnbody++] = NNtree[c].childlist[i];

  nbody->Nsystem = Nsystem;

#if defined(VERIFY_ALL)
  cout << "List all main N-body particles : " << nbody->Nnbody << endl;
  for (i=0; i<nbody->Nnbody; i++) {
    cout << "i : " << i << "    r : " << nbody->nbodydata[i]->r[0] 
	 << "    " << nbody->nbodydata[i]->r[1] << "    Ncomp : " 
	 << nbody->nbodydata[i]->Ncomp << endl;
  }
  cout << "List all main system particles : " << nbody->Nsystem << endl;
  for (i=0; i<nbody->Nsystem; i++) {
    cout << "i : " << i << "    r : " << nbody->system[i].r[0] 
	 << "    " << nbody->system[i].r[1] << "     Ncomp : " 
	 << nbody->system[i].Ncomp << "   Nchildren : " 
	 << nbody->system[i].Nchildren << endl;
  }
#endif

  return;
}



//=============================================================================
//  NbodySystemTree::FindPerturberLists
/// Walk the N-body tree top-down and compute all perturber lists
//=============================================================================
template <int ndim>
void NbodySystemTree<ndim>::FindPerturberLists
(Nbody<ndim> *nbody)                ///< [in] Nbody object containing stars
{
  int c;                            // Node id
  int i;                            // Particle id
  int k;                            // Dimension counter
  int Nstack;                       // ..
  int s;                            // System counter
  int *cellstack;                   // ..

  debug2("[NbodySystemTree::FindPerturberLists]");


  // Start from main root node and walk down entire tree
  c = Nnode - 1;

  // Loop over all systems to find nearest perturbers
  // --------------------------------------------------------------------------
  for (s=0; s<nbody->Nsystem; s++) {

    // Find node on NN-tree corresponding to system particle
    c = nbody->system[s].inode;

    // Add 'sister' cell to stack
    Nstack = 1;
    //cparent = NNtree[c].iparent;
    //if (NNtree[c].c1 == c) cellstack[Nstack++] = NNtree[c]
    

    // Now walk up the tree in turn finding the nearest perturbers
    do {

      

    } while (Nstack > 0);
    

  }
  // --------------------------------------------------------------------------


  return;
}



template class NbodySystemTree<1>;
template class NbodySystemTree<2>;
template class NbodySystemTree<3>;
