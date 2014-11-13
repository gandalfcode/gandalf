//=================================================================================================
//  GodunovSphTree.cpp
//  Contains all functions for building, stocking and walking for the
//  tree for SPH particles.
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
//  GodunovSphKDTree::GodunovSphKDTree
/// SphTree constructor.  Initialises various variables.
//=============================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
GodunovSphKDTree<ndim,ParticleType,TreeCell>::GodunovSphKDTree
 (int Nleafmaxaux,
  int Nmpiaux,
  FLOAT thetamaxsqdaux,
  FLOAT kernrangeaux,
  FLOAT macerroraux,
  string gravity_mac_aux,
  string multipole_aux,
  DomainBox<ndim> *boxaux,
  SphKernel<ndim> *kernaux,
  CodeTiming *timingaux):
 GodunovSphTree<ndim,ParticleType,TreeCell>(Nleafmaxaux,Nmpiaux,thetamaxsqdaux,kernrangeaux,
                                            macerroraux,gravity_mac_aux,multipole_aux,
                                            boxaux,kernaux,timingaux)
{
  // Set-up main tree object
  tree = new KDTree<ndim,ParticleType,TreeCell>(Nleafmaxaux, thetamaxsqdaux, kernrangeaux,
                                                macerroraux, gravity_mac_aux, multipole_aux);

  // Set-up ghost-particle tree object
  ghosttree = new KDTree<ndim,ParticleType,TreeCell>(Nleafmaxaux, thetamaxsqdaux, kernrangeaux,
                                                     macerroraux, gravity_mac_aux, multipole_aux);

#ifdef MPI_PARALLEL
  // Set-up ghost-particle tree object
  mpighosttree = new KDTree<ndim,ParticleType,TreeCell>(Nleafmaxaux, thetamaxsqdaux, kernrangeaux,
                                                       macerroraux, gravity_mac_aux, multipole_aux);

  // Set-up multiple pruned trees, one for each MPI process
  *(prunedtree) = *(new KDTree<ndim,ParticleType,TreeCell>*[Nmpi]);
  for (int j=0; j<Nmpi; j++) {
    prunedtree[j] = new KDTree<ndim,ParticleType,TreeCell>
      (Nleafmaxaux, thetamaxsqdaux, kernrangeaux, macerroraux, gravity_mac_aux, multipole_aux);
  }
#endif
}



//=============================================================================
//  GodunovSphOctTree::GodunovSphOctTree
/// SphTree constructor.  Initialises various variables.
//=============================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
GodunovSphOctTree<ndim,ParticleType,TreeCell>::GodunovSphOctTree
 (int Nleafmaxaux,
  int Nmpiaux,
  FLOAT thetamaxsqdaux,
  FLOAT kernrangeaux,
  FLOAT macerroraux,
  string gravity_mac_aux,
  string multipole_aux,
  DomainBox<ndim> *boxaux,
  SphKernel<ndim> *kernaux,
  CodeTiming *timingaux):
 GodunovSphTree<ndim,ParticleType,TreeCell>(Nleafmaxaux,Nmpiaux,thetamaxsqdaux,kernrangeaux,
                                          macerroraux,gravity_mac_aux,multipole_aux,
                                          boxaux,kernaux,timingaux)
{
  // Set-up main tree object
  tree = new OctTree<ndim,ParticleType,TreeCell>(Nleafmaxaux, thetamaxsqdaux, kernrangeaux,
                                                 macerroraux, gravity_mac_aux, multipole_aux);

  // Set-up ghost-particle tree object
  ghosttree = new OctTree<ndim,ParticleType,TreeCell>(Nleafmaxaux, thetamaxsqdaux, kernrangeaux,
                                                      macerroraux, gravity_mac_aux, multipole_aux);

#ifdef MPI_PARALLEL
  // Set-up ghost-particle tree object
  mpighosttree = new OctTree<ndim,ParticleType,TreeCell>(Nleafmaxaux, thetamaxsqdaux, kernrangeaux,
                                                         macerroraux, gravity_mac_aux, multipole_aux);

  // Set-up multiple pruned trees, one for each MPI process
  *(prunedtree) = *(new OctTree<ndim,ParticleType,TreeCell>*[Nmpi]);
  for (int j=0; j<Nmpi; j++) {
    prunedtree[j] = new OctTree<ndim,ParticleType,TreeCell>
      (Nleafmaxaux, thetamaxsqdaux, kernrangeaux, macerroraux, gravity_mac_aux, multipole_aux);
  }
#endif
}



//=============================================================================
//  GodunovSphTree::GodunovSphTree
/// GodunovSphTree constructor.  Initialises various variables.
//=============================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
GodunovSphTree<ndim,ParticleType,TreeCell>::GodunovSphTree
(int Nleafmaxaux,
 int Nmpiaux,
 FLOAT thetamaxsqdaux,
 FLOAT kernrangeaux,
 FLOAT macerroraux,
 string gravity_mac_aux,
 string multipole_aux,
 DomainBox<ndim> *boxaux,
 SphKernel<ndim> *kernaux,
 CodeTiming *timingaux):
  SphTree<ndim,ParticleType,TreeCell>(Nleafmaxaux,Nmpiaux,thetamaxsqdaux,kernrangeaux,
                                      macerroraux,gravity_mac_aux,multipole_aux,
                                      boxaux,kernaux,timingaux)
{
}



//=============================================================================
//  GodunovSphTree::~GodunovSphTree
/// GodunovSphTree destructor. Deallocates tree memory upon object destruction.
//=============================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
GodunovSphTree<ndim,ParticleType,TreeCell>::~GodunovSphTree()
{
  if (tree->allocated_tree) {
    this->DeallocateMemory();
    tree->DeallocateTreeMemory();
  }
}



template class GodunovSphTree<1,GodunovSphParticle,KDTreeCell>;
template class GodunovSphTree<2,GodunovSphParticle,KDTreeCell>;
template class GodunovSphTree<3,GodunovSphParticle,KDTreeCell>;
template class GodunovSphKDTree<1,GodunovSphParticle,KDTreeCell>;
template class GodunovSphKDTree<2,GodunovSphParticle,KDTreeCell>;
template class GodunovSphKDTree<3,GodunovSphParticle,KDTreeCell>;

template class GodunovSphTree<1,GodunovSphParticle,OctTreeCell>;
template class GodunovSphTree<2,GodunovSphParticle,OctTreeCell>;
template class GodunovSphTree<3,GodunovSphParticle,OctTreeCell>;
template class GodunovSphOctTree<1,GodunovSphParticle,OctTreeCell>;
template class GodunovSphOctTree<2,GodunovSphParticle,OctTreeCell>;
template class GodunovSphOctTree<3,GodunovSphParticle,OctTreeCell>;
