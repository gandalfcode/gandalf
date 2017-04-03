//=================================================================================================
//  NbodySystemTree.cpp
//  Contains all function definitions for constructing N-body nearest-neighbour
//  tree for creating systems and sub-systems.
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
//=================================================================================================


#include <math.h>
#include <string>
#include "Precision.h"
#include "Constants.h"
#include "NbodyParticle.h"
#include "StarParticle.h"
#include "SystemParticle.h"
#include "NbodySystemTree.h"
#include "BinaryOrbit.h"
#include "MergeList.h"
#include "InlineFuncs.h"
#include "Debug.h"
using namespace std;



//=================================================================================================
//  NbodySystemTree::NbodySystemTree()
/// NbodySystemTree constructor
//=================================================================================================
template <int ndim>
NbodySystemTree<ndim>::NbodySystemTree()
{
  allocated_tree = false;
  Nnode      = 0;
  Nnodemax   = 0;
  Nbinary    = 0;
  Ntriple    = 0;
  Nquadruple = 0;
  Norbit     = 0;
  Norbitmax  = 0;
}



//=================================================================================================
//  NbodySystemTree::~NbodySystemTree()
/// NbodySystemTree destructor
//=================================================================================================
template <int ndim>
NbodySystemTree<ndim>::~NbodySystemTree()
{
  if (allocated_tree) DeallocateMemory();
}



//=================================================================================================
//  NbodySystemTree::AllocateMemory
/// Allocate nearest-neighbour tree memory
//=================================================================================================
template <int ndim>
void NbodySystemTree<ndim>::AllocateMemory(int N)
{
  debug2("[NbodySystemTree::AllocateMemory]");

  if (2*N > Nnodemax || !allocated_tree) {
    if (allocated_tree) DeallocateMemory();
    Nnodemax       = 2*N;
    Norbitmax      = N;
    NNtree         = new NNTreeCell<ndim>[Nnodemax];
    orbit          = new BinaryOrbit[Norbitmax];
    allocated_tree = true;
  }

  return;
}



//=================================================================================================
//  NbodySystemTree::DeallocateMemory
/// Deallocate nearest-neighbour tree memory
//=================================================================================================
template <int ndim>
void NbodySystemTree<ndim>::DeallocateMemory(void)
{
  debug2("[NbodySystemTree::DeallocateMemory]");

  if (allocated_tree) {
    delete[] orbit;
    delete[] NNtree;
  }
  allocated_tree = false;

  return;
}



//=================================================================================================
//  NbodySystemTree::CreateNbodySystemTree
/// Creates a nearest neighbour tree based on the positions of all stars
/// contained in the nbody object that is passed.
//=================================================================================================
template <int ndim>
void NbodySystemTree<ndim>::CreateNbodySystemTree
 (Nbody<ndim> *nbody)                  ///< [in] Nbody object containing stars
{
  int i,ii,j,jj;                       // Node ids and counters
  int k;                               // Dimension counter
  int Nfreenode;                       // No. of free (i.e. unattached) nodes
  int *nodelist;                       // List of unattached nodes
  FLOAT dr[ndim];                      // Relative position vector
  FLOAT drsqd;                         // Distance squared

  debug2("[NbodySystemTree::CreateNbodySystemTree]");

  Nnode = 0;
  Nfreenode = 0;

  // Allocate memory for tree
  AllocateMemory(nbody->Nstar);
  nodelist = new int[Nnodemax];

  // Initialise all node variables before adding stars
  for (i=0; i<Nnodemax; i++) {
    NNtree[i].ichild1      = -1;
    NNtree[i].ichild2      = -1;
    NNtree[i].iparent      = -1;
    NNtree[i].inearest     = -1;
    NNtree[i].Ncomp        = 0;
    NNtree[i].Nstar        = 0;
    NNtree[i].Nchildlist   = 0;
    NNtree[i].gpe          = (FLOAT) 0.0;
    NNtree[i].gpe_internal = (FLOAT) 0.0;
    NNtree[i].radius       = (FLOAT) 0.0;
    NNtree[i].rsqdnearest  = big_number_dp;
  }

  // Add one star to each lead node, recording the position and id
  for (i=0; i<nbody->Nstar; i++) {
    nodelist[i] = i;
    NNtree[i].Ncomp = 1;
    NNtree[i].Nstar = 1;
    for (k=0; k<ndim; k++) NNtree[i].rpos[k] = nbody->stardata[i].r[k];
    Nnode++;
    Nfreenode++;
  }


  // Process all remaining unconnected nodes to find new set of mutually
  // nearest neighbours for next phase of tree construction
  //===============================================================================================
  while (Nnode < Nnodemax) {

    // Construct list of remaining nodes
    Nfreenode = 0;
    for (i=0; i<Nnode; i++) {
      if (NNtree[i].iparent == -1) nodelist[Nfreenode++] = i;
    }

    // If we have only one remaining unconnected node (i.e. the root) exit loop
    if (Nfreenode == 1) break;

    // Identify the nearest neighbouring node for each node
    //---------------------------------------------------------------------------------------------
    for (ii=0; ii<Nfreenode; ii++) {
      i = nodelist[ii];
      NNtree[i].inearest = -1;
      NNtree[i].rsqdnearest = big_number;

      for (jj=0; jj<Nfreenode; jj++) {
        j = nodelist[jj];
        if (i == j) continue;

        for (k=0; k<ndim; k++) dr[k] = NNtree[i].rpos[k] - NNtree[j].rpos[k];
        drsqd = DotProduct(dr,dr,ndim);
        if (drsqd < NNtree[i].rsqdnearest) {
          NNtree[i].rsqdnearest = drsqd;
          NNtree[i].inearest = j;
        }
      }
    }
    //---------------------------------------------------------------------------------------------


    // Now identify all mutually nearest neighbours to create a new generation of nodes
    //---------------------------------------------------------------------------------------------
    for (ii=0; ii<Nfreenode-1; ii++) {
      i = nodelist[ii];

      for (jj=ii+1; jj<Nfreenode; jj++) {
        j = nodelist[jj];

        // If each node is the others nearest neighbour, then create a new
        // parent node with the two original nodes as child nodes
        if (NNtree[i].inearest == j && NNtree[j].inearest == i) {
          for (k=0; k<ndim; k++) NNtree[Nnode].rpos[k] = (FLOAT) 0.5*(NNtree[i].rpos[k] + NNtree[j].rpos[k]);
          for (k=0; k<ndim; k++) dr[k] = NNtree[Nnode].rpos[k] - NNtree[i].rpos[k];
          drsqd = DotProduct(dr,dr,ndim);
          NNtree[Nnode].radius  = sqrt(drsqd);
          NNtree[i].iparent     = Nnode;
          NNtree[j].iparent     = Nnode;
          NNtree[Nnode].ichild1 = i;
          NNtree[Nnode].ichild2 = j;
          NNtree[Nnode].Nstar   = NNtree[i].Nstar + NNtree[j].Nstar;
          NNtree[Nnode].Ncomp   = NNtree[i].Ncomp + NNtree[j].Ncomp;
          Nnode++;
#if defined(VERIFY_ALL)
          cout << "Adding new node to tree : " << Nnode - 1 << "   " << i << "    " << j << endl;
#endif
        }

      }

    }
    //---------------------------------------------------------------------------------------------

  };
  //===============================================================================================

  // Free all locally allocated memory
  delete[] nodelist;

#if defined(VERIFY_ALL)
  cout << "Nnode : " << Nnode << endl;
#endif

  return;
}



//=================================================================================================
//  NbodySystemTree::BuildSubSystems
/// Calculate the properties of all nearest-neighbour tree nodes, starting from
/// single-star child-cells and then working up through the parent nodes and
/// finally the root node.
//=================================================================================================
template <int ndim>
void NbodySystemTree<ndim>::BuildSubSystems
 (int n,                               ///< [in] Current integer step
  int level_max,                       ///< [in] Current no. of timestep levels
  FLOAT t,                             ///< [in] Current physical time
  Nbody<ndim> *nbody)                  ///< [inout] Nbody object containing stars
{
  int c,c1,c2;                         // Cell counter
  int i,j;                             // Star and system particle counters
  int k;                               // Dimension counter
  int Nsystem;                         // No. of systems
  FLOAT dr[ndim];                      // Relative position vector
  FLOAT drsqd;                         // Distance squared
  FLOAT dv[ndim];                      // Relative velocity vector
  FLOAT invdrmag;                      // 1 / drmag
  FLOAT ketot = (FLOAT) 0.0;           // Kinetic energy
  FLOAT vmean;                         // Average speed of stars in sub-cluster
  NbodyParticle<ndim> *si;             // Pointer to star/system i
  NbodyParticle<ndim> *sj;             // Pointer to star/system j

  debug2("[NbodySystemTree::BuildSubSystems]");

  // Zero all counters
  nbody->Nnbody  = 0;
  nbody->Nsystem = 0;
  Nsystem        = 0;
  Norbit         = 0;
  Nbinary        = 0;
  Ntriple        = 0;
  Nquadruple     = 0;


  // Loop through all nodes of tree and compute all physical quantities
  //===============================================================================================
  for (c=0; c<Nnode; c++) {

#if defined(VERIFY_ALL)
    cout << "Stocking node : " << c << "    Ncomp : " << NNtree[c].Ncomp << endl;
#endif

    // If node contains a single star, set all properties equal to star values
    //---------------------------------------------------------------------------------------------
    if (c < nbody->Nstar) {
      i = c;
      NNtree[c].m            = nbody->stardata[i].m;
      NNtree[c].h            = nbody->stardata[i].h;
      NNtree[c].gpe          = (FLOAT) 0.5*nbody->stardata[i].m*nbody->stardata[i].gpot;
      NNtree[c].gpe_internal = (FLOAT) 0.0;
      NNtree[c].tcross       = big_number;
      NNtree[c].Nchildlist   = 1;
      NNtree[c].childlist[0] = &(nbody->stardata[i]);
      for (k=0; k<ndim; k++) NNtree[c].r[k]     = nbody->stardata[i].r[k];
      for (k=0; k<ndim; k++) NNtree[c].rpos[k]  = nbody->stardata[i].r[k];
      for (k=0; k<ndim; k++) NNtree[c].v[k]     = nbody->stardata[i].v[k];
      for (k=0; k<ndim; k++) NNtree[c].a[k]     = nbody->stardata[i].a[k];
      for (k=0; k<ndim; k++) NNtree[c].adot[k]  = nbody->stardata[i].adot[k];
      for (k=0; k<ndim; k++) NNtree[c].a2dot[k] = nbody->stardata[i].a2dot[k];
      for (k=0; k<ndim; k++) NNtree[c].a3dot[k] = nbody->stardata[i].a3dot[k];
    }

    // Else, add together both child node properties
    //---------------------------------------------------------------------------------------------
    else {
      c1 = NNtree[c].ichild1;
      c2 = NNtree[c].ichild2;
      NNtree[c].Ncomp = NNtree[c1].Ncomp + NNtree[c2].Ncomp;
      NNtree[c].m     = NNtree[c1].m + NNtree[c2].m;
      NNtree[c].h     = max(NNtree[c1].h,NNtree[c2].h);
      NNtree[c].gpe   = NNtree[c1].gpe + NNtree[c2].gpe;
      NNtree[c].gpe_internal = (FLOAT) 0.0;
      for (k=0; k<ndim; k++) {
        NNtree[c].rpos[k] = (FLOAT) 0.5*(NNtree[c1].rpos[k] + NNtree[c2].rpos[k]);
        NNtree[c].r[k] = (NNtree[c1].m*NNtree[c1].r[k] + NNtree[c2].m*NNtree[c2].r[k])/NNtree[c].m;
        NNtree[c].v[k] = (NNtree[c1].m*NNtree[c1].v[k] + NNtree[c2].m*NNtree[c2].v[k])/NNtree[c].m;
        NNtree[c].a[k] = (NNtree[c1].m*NNtree[c1].a[k] + NNtree[c2].m*NNtree[c2].a[k])/NNtree[c].m;
        NNtree[c].adot[k] =
          (NNtree[c1].m*NNtree[c1].adot[k] + NNtree[c2].m*NNtree[c2].adot[k])/NNtree[c].m;
        NNtree[c].a2dot[k] =
          (NNtree[c1].m*NNtree[c1].a2dot[k] + NNtree[c2].m*NNtree[c2].a2dot[k])/NNtree[c].m;
        NNtree[c].a3dot[k] =
          (NNtree[c1].m*NNtree[c1].a3dot[k] + NNtree[c2].m*NNtree[c2].a3dot[k])/NNtree[c].m;
      }

#if defined(VERIFY_ALL)
      cout << "Stocking system data : " << c << "    " << NNtree[c].r[0]
           << "    " << NNtree[c].r[1] << endl;
      cout << "Child energies : " << NNtree[c].gpe << "    " << NNtree[c1].gpe << "    "
           << NNtree[c2].gpe << "     " << NNtree[c].gpe_internal << "     "
           << NNtree[c1].gpe_internal << "    " << NNtree[c2].gpe_internal << endl;
#endif


      // If node contains less than the maximum allowed number of components, then
      // check if node can be replaced by a new system particle or not.
      //-------------------------------------------------------------------------------------------
      if (NNtree[c].Ncomp <= Ncompmax) {

        // Merge lists of stars/systems in child nodes
        NNtree[c].Nchildlist = NNtree[c1].Nchildlist + NNtree[c2].Nchildlist;
        for (i=0; i<NNtree[c1].Nchildlist; i++) NNtree[c].childlist[i] = NNtree[c1].childlist[i];
        for (i=0; i<NNtree[c2].Nchildlist; i++) {
          NNtree[c].childlist[i+NNtree[c1].Nchildlist] = NNtree[c2].childlist[i];
        }

#if defined(VERIFY_ALL)
        cout << "No. of child systems for " << c << "   :   " << NNtree[c].Nchildlist << endl;
#endif

        // Compute internal kinetic energy (relative to local centre-of-mass) and the
        // velocity dispersion, vmean.
        ketot = (FLOAT) 0.0;
        for (i=0; i<NNtree[c].Nchildlist; i++) {
          si = NNtree[c].childlist[i];
          for (k=0; k<ndim; k++) dv[k] = si->v[k] - NNtree[c].v[k];
          ketot += (FLOAT) 0.5*si->m*DotProduct(dv,dv,ndim);
        }
        vmean = sqrt((FLOAT) 2.0*ketot/NNtree[c].m);

        // Compute internal gravitational potential energy
        // (N.B. counting all pairs once so no need for factor of 1/2)
        for (i=0; i<NNtree[c].Nchildlist-1; i++) {
          si = NNtree[c].childlist[i];
          for (j=i+1; j<NNtree[c].Nchildlist; j++) {
            sj = NNtree[c].childlist[j];
            for (k=0; k<ndim; k++) dr[k] = sj->r[k] - si->r[k];
            drsqd = DotProduct(dr,dr,ndim);
            invdrmag = (FLOAT) 1.0 / (sqrt(drsqd) + small_number_dp);
            NNtree[c].gpe_internal += si->m*sj->m*invdrmag;
          }
        }

        // Estimate the crossing time as tcross = Rgrav / vdisp where Rgrav = sqrt(m^2/G/Egrav).
        // Should give similar timescale to binary period in case of bound binary.
        NNtree[c].tcross = sqrt(NNtree[c].m*NNtree[c].m/NNtree[c].gpe_internal)/vmean;

#if defined(VERIFY_ALL)
        cout << "gpe : " << NNtree[c].gpe << "    gpe_internal : "
             << NNtree[c].gpe_internal << "    ketot : " << ketot << endl;
        //if (NNtree[c].gpe < NNtree[c].gpe_internal) {
        //  cout << "Grav. energy problem" << NNtree[c].gpe - NNtree[c].gpe_internal << endl;
        //}
        cout << "Checking system criteria : "
             << fabs(NNtree[c].gpe - NNtree[c].gpe_internal)/NNtree[c].gpe
             << "    " << gpesoft << "    " << NNtree[c].Ncomp << "     " << Ncompmax << endl;
#endif


        // Now check energies and decide if we should create a new sub-system
        // object.  If yes, create new system in main arrays.
        //-----------------------------------------------------------------------------------------
        if (fabs(NNtree[c].gpe - NNtree[c].gpe_internal) < gpesoft*NNtree[c].gpe) {
//cout << "CREATING NEW SYSTEM" << "    " << Nsystem << "    " << NNtree[c].Ncomp << endl;
//cout << "A : " << NNtree[c].a[0] << "    " << NNtree[c].a[1] << "    " << NNtree[c].a[2] << endl;
          // Copy centre-of-mass properties of new sub-system
          nbody->system[Nsystem].inode        = c;
          nbody->system[Nsystem].nlast        = n;
          nbody->system[Nsystem].nstep        = 1;
          nbody->system[Nsystem].level        = level_max;
          nbody->system[Nsystem].tlast        = t;
          nbody->system[Nsystem].dt_internal  = NNtree[c].tcross;
          nbody->system[Nsystem].m            = NNtree[c].m;
          nbody->system[Nsystem].h            = NNtree[c].h;
          nbody->system[Nsystem].gpe_internal = NNtree[c].gpe_internal;
          for (k=0; k<ndim; k++) {
            nbody->system[Nsystem].r[k]     = NNtree[c].r[k];
            nbody->system[Nsystem].v[k]     = NNtree[c].v[k];
            nbody->system[Nsystem].a[k]     = NNtree[c].a[k];
            nbody->system[Nsystem].adot[k]  = NNtree[c].adot[k];
            nbody->system[Nsystem].a2dot[k] = NNtree[c].a2dot[k];
            nbody->system[Nsystem].a3dot[k] = NNtree[c].a3dot[k];
            nbody->system[Nsystem].r0[k]    = NNtree[c].r[k];
            nbody->system[Nsystem].v0[k]    = NNtree[c].v[k];
            nbody->system[Nsystem].a0[k]    = NNtree[c].a[k];
            nbody->system[Nsystem].adot0[k] = NNtree[c].adot[k];
          }

          // Compute and store binary properties if bound
          if (NNtree[c].Ncomp == 2) {
            ComputeNewBinaryOrbit(c1,c2,nbody->Nstar,NNtree[c],NNtree[c1],NNtree[c2]);
          }

          // Copy list of contained N-body particles (either systems or stars)
          // to system object, and also total no. of stars/components
          nbody->system[Nsystem].Ncomp     = NNtree[c].Ncomp;
          nbody->system[Nsystem].Nchildren = 0;
          nbody->system[Nsystem].Npert     = 0;
          for (i=0; i<NNtree[c].Nchildlist; i++) {
            nbody->system[Nsystem].children[nbody->system[Nsystem].Nchildren++]
              = NNtree[c].childlist[i];
          }

#if defined(VERIFY_ALL)
          cout << "Adding binary system : " << c << "    " << c1 << "    "
            << c2 << "     " << nbody->system[Nsystem].Ncomp << endl;
          cout << "System no. : " << Nsystem << endl;
#endif

          // Finally, clear list and append newly created system particle.
          // Also, zero internal energy to allow detection of hierarchical systems.
          NNtree[c].Ncomp        = 1;
          NNtree[c].Nchildlist   = 1;
          NNtree[c].childlist[0] = &(nbody->system[Nsystem]);
          NNtree[c].gpe          = NNtree[c].gpe - NNtree[c].gpe_internal;
          NNtree[c].gpe_internal = (FLOAT) 0.0;
          Nsystem++;

        }

      }

      // If lists contain more than the maximum allowed number of components,
      // then flush list of systems to main arrays
      //-------------------------------------------------------------------------------------------
      else {

        // Merge lists of stars/systems in child nodes
        for (i=0; i<NNtree[c1].Nchildlist; i++) {
          nbody->nbodydata[nbody->Nnbody++] = NNtree[c1].childlist[i];
        }
        for (i=0; i<NNtree[c2].Nchildlist; i++) {
          nbody->nbodydata[nbody->Nnbody++] = NNtree[c2].childlist[i];
        }
        NNtree[c].Nchildlist = 0;
        NNtree[c1].Nchildlist = 0;
        NNtree[c2].Nchildlist = 0;

      }
      //-------------------------------------------------------------------------------------------

    }
    //---------------------------------------------------------------------------------------------

  }
  //===============================================================================================


  // For root cell, push all remaining systems (or single sub-systems) to list
  c = Nnode - 1;
  for (i=0; i<NNtree[c].Nchildlist; i++) {
    nbody->nbodydata[nbody->Nnbody++] = NNtree[c].childlist[i];
  }
  nbody->Nsystem = Nsystem;

#if defined(VERIFY_ALL)
  cout << "No. of remaining systems in root node " << c << "  :  " << NNtree[c].Nchildlist << endl;
  cout << "List all main N-body particles : " << nbody->Nnbody << endl;
  for (i=0; i<nbody->Nnbody; i++) {
    cout << "i : " << i << "      Ncomp : " << nbody->nbodydata[i]->Ncomp << endl;
  }
  cout << "List all main system particles : " << nbody->Nsystem << endl;
  for (i=0; i<nbody->Nsystem; i++) {
    cout << "i : " << i << "      Ncomp : " << nbody->system[i].Ncomp << "   Nchildren : "
         << nbody->system[i].Nchildren << endl;
  }
#endif

  return;
}



//=================================================================================================
//  NbodySystemTree::RestockTreeNodes
/// Recompute the properties of all nearest-neighbour tree nodes while simulataneously checking
/// that the tree is still valid (i.e. all sister nodes are still mutually nearest neighbours).
/// If tree is found to be invalid, then flag the tree to be re-built.
//=================================================================================================
template <int ndim>
void NbodySystemTree<ndim>::RestockTreeNodes
 (Nbody<ndim> *nbody)                  ///< [in] Nbody object containing stars
{
  //bool validtrue = true;               // Flag if tree is still valid or not
  int c;                               // Node counter
  int c1;                              // Child cell 1 id
  int c2;                              // Child cell 2 id
  int i;                               // Counter
  int k;                               // Dimension counter
  //FLOAT dr[ndim];                      // Relative position vector
  //FLOAT drsqd;                         // Distance squared

  debug2("[NbodySystemTree::RestockTreeNodes]");


  //===============================================================================================
  for (c=0; c<Nnode; c++) {

    // If node contains one star, set all properties equal to star values
    //---------------------------------------------------------------------------------------------
    if (c < nbody->Nstar) {

      i = c;
      NNtree[c].m            = nbody->stardata[i].m;
      NNtree[c].h            = nbody->stardata[i].h;
      NNtree[c].radius       = (FLOAT) 0.0;
      NNtree[c].gpe          = (FLOAT) 0.5*nbody->stardata[i].m*nbody->stardata[i].gpot;
      NNtree[c].gpe_internal = (FLOAT) 0.0;
      NNtree[c].tcross       = big_number;
      NNtree[c].Nchildlist   = 1;
      NNtree[c].childlist[0] = &(nbody->stardata[i]);
      for (k=0; k<ndim; k++) NNtree[c].rpos[k]  = nbody->stardata[i].r[k];
      for (k=0; k<ndim; k++) NNtree[c].r[k]     = nbody->stardata[i].r[k];
      for (k=0; k<ndim; k++) NNtree[c].v[k]     = nbody->stardata[i].v[k];
      for (k=0; k<ndim; k++) NNtree[c].a[k]     = nbody->stardata[i].a[k];
      for (k=0; k<ndim; k++) NNtree[c].adot[k]  = nbody->stardata[i].adot[k];
      for (k=0; k<ndim; k++) NNtree[c].a2dot[k] = nbody->stardata[i].a2dot[k];
      for (k=0; k<ndim; k++) NNtree[c].a3dot[k] = nbody->stardata[i].a3dot[k];

    }

    // Else, add together both child node properties
    //---------------------------------------------------------------------------------------------
    else {

      c1 = NNtree[c].ichild1;
      c2 = NNtree[c].ichild2;
      NNtree[c].Ncomp = NNtree[c1].Ncomp + NNtree[c2].Ncomp;
      NNtree[c].m = NNtree[c1].m + NNtree[c2].m;
      NNtree[c].h = max(NNtree[c1].h,NNtree[c2].h);
      NNtree[c].gpe = NNtree[c1].gpe + NNtree[c2].gpe;
      for (k=0; k<ndim; k++) {
        NNtree[c].rpos[k] = (FLOAT) 0.5*(NNtree[c1].rpos[k] + NNtree[c2].rpos[k]);
        NNtree[c].r[k] = (NNtree[c1].m*NNtree[c1].r[k] + NNtree[c2].m*NNtree[c2].r[k])/NNtree[c].m;
        NNtree[c].v[k] = (NNtree[c1].m*NNtree[c1].v[k] + NNtree[c2].m*NNtree[c2].v[k])/NNtree[c].m;
        NNtree[c].a[k] = (NNtree[c1].m*NNtree[c1].a[k] + NNtree[c2].m*NNtree[c2].a[k])/NNtree[c].m;
        NNtree[c].adot[k] =
          (NNtree[c1].m*NNtree[c1].adot[k] + NNtree[c2].m*NNtree[c2].adot[k])/NNtree[c].m;
        NNtree[c].a2dot[k] =
          (NNtree[c1].m*NNtree[c1].a2dot[k] + NNtree[c2].m*NNtree[c2].a2dot[k])/NNtree[c].m;
        NNtree[c].a3dot[k] =
          (NNtree[c1].m*NNtree[c1].a3dot[k] + NNtree[c2].m*NNtree[c2].a3dot[k])/NNtree[c].m;
        //dr[k] = NNtree[c].rpos[k] - NNtree[c1].rpos[k];
      }
      //drsqd = DotProduct(dr,dr,ndim);

      NNtree[c].gpe_internal = (FLOAT) 0.0;

    }
    //---------------------------------------------------------------------------------------------


  }
  //===============================================================================================

  return;
}



//=================================================================================================
//  NbodySystemTree::ComputeNewBinaryOrbit
/// Compute orbital properties of binary/multiple system from nearest-neighbour tree nodes.
/// Stores all orbital elements in 'orbit' array.
//=================================================================================================
template <int ndim>
void NbodySystemTree<ndim>::ComputeNewBinaryOrbit
 (const int c1,                        ///< [in] ..
  const int c2,                        ///< [in] ..
  const int Nstar,                     ///< [in] ..
  NNTreeCell<ndim> &cell,              ///< [in] ..
  NNTreeCell<ndim> &cell1,             ///< [in] ..
  NNTreeCell<ndim> &cell2)             ///< [in] ..
{
  int k;                               // Dimension counter
  FLOAT angmomsqd = (FLOAT) 0.0;       // Angular momentum squared
  FLOAT dr[ndim];                      // Relative position vector
  FLOAT drsqd;                         // Distance squared
  FLOAT dv[ndim];                      // Relative velocity vector
  FLOAT invdrmag;                      // 1 / drmag
  FLOAT mu;                            // Reduced mass of binary system

  // Compute properties of binary relative to the system COM using reduced mass
  mu = cell1.m*cell2.m/cell.m;
  for (k=0; k<ndim; k++) dv[k] = cell1.v[k] - cell2.v[k];
  for (k=0; k<ndim; k++) dr[k] = cell1.r[k] - cell2.r[k];
  if (ndim == 2) {
    orbit[Norbit].angmom[2] = mu*(dr[0]*dv[1] - dr[1]*dv[0]);
    angmomsqd = orbit[Norbit].angmom[2]*orbit[Norbit].angmom[2];
  }
  else if (ndim == 3) {
    orbit[Norbit].angmom[0] = mu*(dr[1]*dv[2] - dr[2]*dv[1]);
    orbit[Norbit].angmom[1] = mu*(dr[2]*dv[0] - dr[0]*dv[2]);
    orbit[Norbit].angmom[2] = mu*(dr[0]*dv[1] - dr[1]*dv[0]);
    angmomsqd = DotProduct(orbit[Norbit].angmom,orbit[Norbit].angmom,ndim);
  }
  drsqd = DotProduct(dr,dr,ndim);
  invdrmag = (FLOAT) 1.0 / (sqrt(drsqd) + small_number_dp);
  orbit[Norbit].ichild1 = c1;
  orbit[Norbit].ichild2 = c2;
  orbit[Norbit].m = cell.m;
  for (k=0; k<ndim; k++) orbit[Norbit].r[k] = cell.r[k];
  for (k=0; k<ndim; k++) orbit[Norbit].v[k] = cell.v[k];
  orbit[Norbit].binen = (FLOAT) 0.5*DotProduct(dv,dv,ndim) - cell.m*invdrmag;

  // If system is unbound, skip recording orbit
  if (orbit[Norbit].binen >= (FLOAT) 0.0) return;

  // Compute orbital properties
  orbit[Norbit].sma    = -(FLOAT) 0.5*cell.m/orbit[Norbit].binen;
  orbit[Norbit].period = twopi*sqrt(pow(orbit[Norbit].sma,3)/orbit[Norbit].m);
  orbit[Norbit].ecc    = (FLOAT) 1.0 - angmomsqd/(orbit[Norbit].m*orbit[Norbit].sma*mu*mu);
  orbit[Norbit].ecc    = sqrt(max((FLOAT) 0.0, orbit[Norbit].ecc));
  if (cell1.m > cell2.m) orbit[Norbit].q = cell2.m/cell1.m;
  else orbit[Norbit].q = cell1.m/cell2.m;

  // Increment object counter depending on multiple system type
  if (c1 < Nstar && c2 < Nstar) {
    orbit[Norbit].systemtype = "binary";
    Nbinary++;
  }
  else if (c1 < Nstar || c2 < Nstar) {
    orbit[Norbit].systemtype = "triple";
    Ntriple++;
  }
  else {
    orbit[Norbit].systemtype = "quadruple";
    Nquadruple++;
  }
  Norbit++;

  return;
}



//=================================================================================================
//  NbodySystemTree::FindBinarySystems
/// Search through the nearest neighbour tree to identify binaries and higher-order hierarchical
/// systems that are both mutually nearest neighbours and bound systems.
//=================================================================================================
template <int ndim>
void NbodySystemTree<ndim>::FindBinarySystems
 (Nbody<ndim> *nbody)                  ///< [inout] Nbody object containing stars
{
  int c,c1,c2;                         // Cell counter
  int i;                               // Star and system particle counters
  int k;                               // Dimension counter
  FLOAT dr[ndim];                      // Relative position vector
  FLOAT drsqd;                         // Distance squared
  FLOAT invdrmag;                      // 1 / drmag

  debug2("[NbodySystemTree::FindBinarySystems]");

  // Set all counters
  Norbit     = 0;
  Nbinary    = 0;
  Ntriple    = 0;
  Nquadruple = 0;

  // Loop through all nodes of tree and compute all physical quantities
  //===============================================================================================
  for (c=0; c<Nnode; c++) {

    // If node contains one star, set all properties equal to star values
    //---------------------------------------------------------------------------------------------
    if (c < nbody->Nstar) {

      i = c;
      NNtree[c].m = nbody->stardata[i].m;
      NNtree[c].h = nbody->stardata[i].h;
      for (k=0; k<ndim; k++) NNtree[c].r[k] = nbody->stardata[i].r[k];
      for (k=0; k<ndim; k++) NNtree[c].v[k] = nbody->stardata[i].v[k];
      for (k=0; k<ndim; k++) NNtree[c].a[k] = nbody->stardata[i].a[k];
      for (k=0; k<ndim; k++) NNtree[c].adot[k] = nbody->stardata[i].adot[k];
      NNtree[c].gpe = (FLOAT) 0.5*nbody->stardata[i].m*nbody->stardata[i].gpot;
      NNtree[c].gpe_internal = (FLOAT) 0.0;

    }

    // Else, add together both child node properties
    //---------------------------------------------------------------------------------------------
    else {

      c1 = NNtree[c].ichild1;
      c2 = NNtree[c].ichild2;
      NNtree[c].Ncomp        = NNtree[c1].Ncomp + NNtree[c2].Ncomp;
      NNtree[c].m            = NNtree[c1].m + NNtree[c2].m;
      NNtree[c].h            = max(NNtree[c1].h,NNtree[c2].h);
      NNtree[c].gpe          = NNtree[c1].gpe + NNtree[c2].gpe;
      NNtree[c].gpe_internal = (FLOAT) 0.0;
      for (k=0; k<ndim; k++) {
        NNtree[c].r[k] = (NNtree[c1].m*NNtree[c1].r[k] + NNtree[c2].m*NNtree[c2].r[k])/NNtree[c].m;
        NNtree[c].v[k] = (NNtree[c1].m*NNtree[c1].v[k] + NNtree[c2].m*NNtree[c2].v[k])/NNtree[c].m;
        NNtree[c].a[k] = (NNtree[c1].m*NNtree[c1].a[k] + NNtree[c2].m*NNtree[c2].a[k])/NNtree[c].m;
        NNtree[c].adot[k] =
          (NNtree[c1].m*NNtree[c1].adot[k] + NNtree[c2].m*NNtree[c2].adot[k])/NNtree[c].m;
      }


      // If node contains within maximum allowed number of components, then
      // check if node is a new system or not
      //-------------------------------------------------------------------------------------------
      if (NNtree[c].Ncomp <= Ncompmax) {

        // Compute internal gravitational potential energy
        for (k=0; k<ndim; k++) dr[k] = NNtree[c1].r[k] - NNtree[c2].r[k];
        drsqd = DotProduct(dr,dr,ndim);
        invdrmag = (FLOAT) 1.0 / (sqrt(drsqd) + small_number_dp);
        NNtree[c].gpe_internal += NNtree[c1].m*NNtree[c2].m*invdrmag;

        // Now check energies and decide if we should create a new sub-system
        // object.  If yes, create new system in main arrays
        //-----------------------------------------------------------------------------------------
        if (fabs(NNtree[c].gpe - NNtree[c].gpe_internal) < gpesoft*NNtree[c].gpe) {

          // Compute and store binary properties if bound
          if (NNtree[c].Ncomp == 2) {
            ComputeNewBinaryOrbit(c1,c2,nbody->Nstar,NNtree[c],NNtree[c1],NNtree[c2]);
          }

          // Finally, clear list and append newly created system particle.
          // Also, zero internal energy to allow detection of hierarchical systems.
          NNtree[c].Ncomp = 1;
          NNtree[c].gpe = NNtree[c].gpe - NNtree[c].gpe_internal;
          NNtree[c].gpe_internal = (FLOAT) 0.0;

        }
        //-----------------------------------------------------------------------------------------

      }
      //-------------------------------------------------------------------------------------------

    }
    //---------------------------------------------------------------------------------------------

  }
  //===============================================================================================


  return;
}



//=================================================================================================
//  NbodySystemTree::FindPerturberLists
/// Walk the N-body tree top-down and compute all perturber lists
//=================================================================================================
template <int ndim>
void NbodySystemTree<ndim>::FindPerturberLists
 (Nbody<ndim> *nbody)                  ///< [in] Nbody object containing stars
{
  int c;                               // Node id
  int caux;                            // Aux. node id
  int k;                               // Dimension counter
  int Nstack;                          // No. of unprocessed nodes on stack
  int s;                               // System counter
  int *cellstack;                      // ids of unprocessed nodes
  FLOAT dr[ndim];                      // Relative position vector
  FLOAT drsqd;                         // Distance squared
  FLOAT gpeaux;                        // Grav. potential energy
  FLOAT gpe_internal;                  // Internal gpe
  FLOAT ms;                            // Mass of system
  FLOAT rs[ndim];                      // Position of system
  SystemParticle<ndim> *s1;            // Pointer to system particle

  debug2("[NbodySystemTree::FindPerturberLists]");

  cellstack = new int[nbody->Nstar];
  //cout << "Nsystem : " << nbody->Nsystem << "    Nnode : " << Nnode << endl;


  // Loop over all systems to find nearest perturbers
  //===============================================================================================
  for (s=0; s<nbody->Nsystem; s++) {
    s1 = &(nbody->system[s]);
    s1->Npert = 0;
    //continue;

    // Find node on NN-tree corresponding to system particle
    c = s1->inode;
    ms = NNtree[c].m;
    gpe_internal = s1->gpe_internal;
    for (k=0; k<ndim; k++) rs[k] = s1->r[k];

    // If system contains all N-body particles, no perturbers exist
    //if (c == Nnode - 1) cout << "System contains all stars; no perturbers : " << c << endl;
    if (c == Nnode - 1) continue;

    //cout << "Finding perturbers for system " << s << "   " << nbody->Nsystem << endl;
    //cout << "gpe_internal : " << s << "   " << gpe_internal << endl;

    // Zero stack
    Nstack = 0;


    //---------------------------------------------------------------------------------------------
    do {

      // Check if nearest node contains potential perturbers
      caux = NNtree[c].inearest;
      c = NNtree[c].iparent;
      cellstack[Nstack++] = caux;
      //cout << "Now looking at sister node " << caux << "   " << c << endl;


      // Now walk through the stack interrogating any cells to see if they are valid perturbers.
      //-------------------------------------------------------------------------------------------
      do {
        Nstack--;
        caux = cellstack[Nstack];

        // Check if node may contain potential perturbers or not
        for (k=0; k<ndim; k++) dr[k] = rs[k] - NNtree[caux].r[k];
        drsqd = DotProduct(dr, dr, ndim);
        gpeaux = ms*NNtree[caux].m/sqrt(drsqd);

        //cout << "Checking perturber;  s : " << s << "    c : " << c << "    caux : "
          //   << caux << "    drsqd : " << drsqd << "    gpeaux : " << gpeaux
          //   << "      " << gpehard*gpe_internal << endl;

        // If perturber grav. energy is high enough, then add cell as potential perturber.
        if (gpeaux > gpehard*gpe_internal) {

          // If node is star/system, add to perturber list.
          // Else, open cell and add child cells to stack.
          if (NNtree[caux].Ncomp == 1 && caux != Nnode - 1) {
            if (s1->Npert < Npertmax) s1->perturber[s1->Npert++] = NNtree[caux].childlist[0];
            //cout << "Added star " << NNtree[caux].childlist[0]->istar << " as perturber" << endl;
          }
          else {
            cellstack[Nstack++] = NNtree[caux].ichild1;
            cellstack[Nstack++] = NNtree[caux].ichild2;
          }
        }

        //cout << "Nstack : " << Nstack << "     Npert : " << s1->Npert << endl;

      } while (Nstack > 0 && s1->Npert < Npertmax);
      //-------------------------------------------------------------------------------------------

    } while (s1->Npert > 0 && s1->Npert < Npertmax && c < Nnode - 1);
    //---------------------------------------------------------------------------------------------

    //cout << "Found " << s1->Npert << " perturbers for system " << s << endl;

  }
  //===============================================================================================


  delete[] cellstack;

  return;
}



//=================================================================================================
//  NbodySystemTree::OutputBinaryProperties
/// Output all binary properties to screen.
//=================================================================================================
template <int ndim>
void NbodySystemTree<ndim>::OutputBinaryProperties
 (Nbody<ndim> *nbody)                  ///< [in] Nbody object containing stars
{
  int i;                               // Binary counter

  debug2("[NbodySystemTree::OutputBinaryProperties]");

  cout << "---------------------------------------" << endl;
  cout << "Binary statistics" << endl;
  cout << "---------------------------------------" << endl;
  cout << "No. of systems           : " << nbody->Nsystem << endl;
  cout << "No. of orbits            : " << Norbit << endl;
  cout << "No. of binary orbits     : " << Nbinary << endl;  //Norbit - 2*Ntriple - 3*Nquadruple
  cout << "No. of triple orbits     : " << Ntriple << endl;
  cout << "No. of quadruple orbits  : " << Nquadruple << endl;

  for (i=0; i<Norbit; i++) {
    cout << "---------------------------------------" << endl;
    cout << "System " << i << " : " << orbit[i].systemtype << endl;
    cout << "Components           : " << orbit[i].ichild1 << "    " << orbit[i].ichild2 << endl;
    cout << "Semi-major axis      : " << orbit[i].sma << endl;
    cout << "Orbital period       : " << orbit[i].period << endl;
    cout << "Orbital eccentricity : " << orbit[i].ecc << endl;
    if (ndim == 2) {
      cout << "Orbital ang. mom     : " << orbit[i].angmom[2] << endl;
    }
    else {
      cout << "Orbital ang. mom     : " << orbit[i].angmom[0] << "    "
           << orbit[i].angmom[1] << "    " << orbit[i].angmom[2] << endl;
    }
    cout << "Total mass           : " << orbit[i].m << "    " << orbit[i].q << endl;
  }
  cout << "---------------------------------------" << endl;

  return;
}



template class NbodySystemTree<1>;
template class NbodySystemTree<2>;
template class NbodySystemTree<3>;
