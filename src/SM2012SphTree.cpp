//=================================================================================================
//  SM2012SphTree.cpp
//  Contains all functions for building, stocking and walking for the
//  binary KD tree for SPH particles.
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



//=================================================================================================
//  SM2012SphKDTree::SM2012SphKDTree
/// SphTree constructor.  Initialises various variables.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
SM2012SphKDTree<ndim,ParticleType,TreeCell>::SM2012SphKDTree
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
 SM2012SphTree<ndim,ParticleType,TreeCell>(Nleafmaxaux,Nmpiaux,thetamaxsqdaux,kernrangeaux,
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



//=================================================================================================
//  SM2012SphOctTree::Sm2012SphOctTree
/// SphTree constructor.  Initialises various variables.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
SM2012SphOctTree<ndim,ParticleType,TreeCell>::SM2012SphOctTree
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
 SM2012SphTree<ndim,ParticleType,TreeCell>(Nleafmaxaux,Nmpiaux,thetamaxsqdaux,kernrangeaux,
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



//=================================================================================================
//  SM2012SphTree::SM2012SphTree
/// SM2012SphTree constructor.  Initialises various variables.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
SM2012SphTree<ndim,ParticleType,TreeCell>::SM2012SphTree
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



//=================================================================================================
//  SM2012SphTree::~SM2012SphTree
/// SM2012SphTree destructor.  Deallocates tree memory upon object destruction.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
SM2012SphTree<ndim,ParticleType,TreeCell>::~SM2012SphTree()
{
  if (tree->allocated_tree) {
    this->DeallocateMemory();
    tree->DeallocateTreeMemory();
  }
}



template class SM2012SphTree<1,SM2012SphParticle,KDTreeCell>;
template class SM2012SphTree<2,SM2012SphParticle,KDTreeCell>;
template class SM2012SphTree<3,SM2012SphParticle,KDTreeCell>;
template class SM2012SphKDTree<1,SM2012SphParticle,KDTreeCell>;
template class SM2012SphKDTree<2,SM2012SphParticle,KDTreeCell>;
template class SM2012SphKDTree<3,SM2012SphParticle,KDTreeCell>;

template class SM2012SphTree<1,SM2012SphParticle,OctTreeCell>;
template class SM2012SphTree<2,SM2012SphParticle,OctTreeCell>;
template class SM2012SphTree<3,SM2012SphParticle,OctTreeCell>;
template class SM2012SphOctTree<1,SM2012SphParticle,OctTreeCell>;
template class SM2012SphOctTree<2,SM2012SphParticle,OctTreeCell>;
template class SM2012SphOctTree<3,SM2012SphParticle,OctTreeCell>;
