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
#include "Particle.h"
#include "Debug.h"
#if defined _OPENMP
#include <omp.h>
#endif
using namespace std;



//=================================================================================================
//  SM2012SphKDTree::SM2012SphKDTree
/// SM2012SphKDTree constructor.  Initialises various variables and creates tree objects.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
SM2012SphKDTree<ndim,ParticleType,TreeCell>::SM2012SphKDTree
 (int _Nleafmax, int _Nmpi, int _pruning_level_min, int _pruning_level_max, FLOAT _thetamaxsqd,
  FLOAT _kernrange, FLOAT _macerror, string _gravity_mac, string _multipole,
  DomainBox<ndim>* _box, SmoothingKernel<ndim>* _kern, CodeTiming* _timing):
 NeighbourSearch<ndim>(_kernrange, _box, _kern, _timing),
 SM2012SphTree<ndim,ParticleType,TreeCell>
  (_Nleafmax, _Nmpi, _pruning_level_min, _pruning_level_max, _thetamaxsqd,
   _kernrange, _macerror, _gravity_mac, _multipole, _box, _kern, _timing)
{
  // Set-up main tree object
  tree = new KDTree<ndim,ParticleType,TreeCell>(_Nleafmax, _thetamaxsqd, _kernrange,
                                                _macerror, _gravity_mac, _multipole, *_box);

  // Set-up ghost-particle tree object
  ghosttree = new KDTree<ndim,ParticleType,TreeCell>(_Nleafmax, _thetamaxsqd, _kernrange,
                                                     _macerror, _gravity_mac, _multipole, *_box);

#ifdef MPI_PARALLEL
  // Set-up ghost-particle tree object
  mpighosttree = new KDTree<ndim,ParticleType,TreeCell>(_Nleafmax, _thetamaxsqd, _kernrange,
                                                        _macerror, _gravity_mac, _multipole, *_box);

  // Set-up multiple pruned trees, one for each MPI process
  KDTree<ndim,ParticleType,TreeCell>** prunedtree_derived = new KDTree<ndim,ParticleType,TreeCell>*[Nmpi];
  prunedtree = (Tree<ndim,ParticleType,TreeCell> **) prunedtree_derived;
  KDTree<ndim,ParticleType,TreeCell>** sendprunedtree_derived = new KDTree<ndim,ParticleType,TreeCell>*[Nmpi];
  sendprunedtree = (Tree<ndim,ParticleType,TreeCell> **) sendprunedtree_derived;

  for (int i=0; i<Nmpi; i++) {
    prunedtree[i] = new KDTree<ndim,ParticleType,TreeCell>
     (_Nleafmax, _thetamaxsqd, _kernrange, _macerror, _gravity_mac, _multipole, *_box);
  }
  for (int i=0; i<Nmpi; i++) {
    sendprunedtree[i] = new KDTree<ndim,ParticleType,TreeCell>
     (_Nleafmax, _thetamaxsqd, _kernrange, _macerror, _gravity_mac, _multipole, *_box);
  }
#endif
}



//=================================================================================================
//  SM2012SphOctTree::SM2012SphOctTree
/// SphTree constructor.  Initialises various variables.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
SM2012SphOctTree<ndim,ParticleType,TreeCell>::SM2012SphOctTree
 (int _Nleafmax, int _Nmpi, int _pruning_level_min, int _pruning_level_max, FLOAT _thetamaxsqd,
  FLOAT _kernrange, FLOAT _macerror, string _gravity_mac, string _multipole,
  DomainBox<ndim>* _box, SmoothingKernel<ndim>* _kern, CodeTiming* _timing):
 NeighbourSearch<ndim>(_kernrange, _box, _kern, _timing),
 SM2012SphTree<ndim,ParticleType,TreeCell>
  (_Nleafmax, _Nmpi, _pruning_level_min, _pruning_level_max, _thetamaxsqd,
   _kernrange, _macerror, _gravity_mac, _multipole, _box, _kern, _timing)
{
  // Set-up main tree object
  tree = new OctTree<ndim,ParticleType,TreeCell>(_Nleafmax, _thetamaxsqd, _kernrange,
                                                 _macerror, _gravity_mac, _multipole, *_box);

  // Set-up ghost-particle tree object
  ghosttree = new OctTree<ndim,ParticleType,TreeCell>(_Nleafmax, _thetamaxsqd, _kernrange,
                                                      _macerror, _gravity_mac, _multipole, *_box);

#ifdef MPI_PARALLEL
  // Set-up ghost-particle tree object
  mpighosttree = new OctTree<ndim,ParticleType,TreeCell>(_Nleafmax, _thetamaxsqd, _kernrange,
                                                         _macerror, _gravity_mac, _multipole, *_box);

  // Set-up multiple pruned trees, one for each MPI process
  //*(prunedtree) = *(new OctTree<ndim,ParticleType,TreeCell>*[Nmpi]);
  // Set-up multiple pruned trees, one for each MPI process
  OctTree<ndim,ParticleType,TreeCell>** prunedtree_derived = new OctTree<ndim,ParticleType,TreeCell>*[Nmpi];
  prunedtree = (Tree<ndim,ParticleType,TreeCell> **) prunedtree_derived;
  OctTree<ndim,ParticleType,TreeCell>** sendprunedtree_derived = new OctTree<ndim,ParticleType,TreeCell>*[Nmpi];
  sendprunedtree = (Tree<ndim,ParticleType,TreeCell> **) sendprunedtree_derived;

  for (int j=0; j<Nmpi; j++) {
    prunedtree[j] = new OctTree<ndim,ParticleType,TreeCell>
     (_Nleafmax, _thetamaxsqd, _kernrange, _macerror, _gravity_mac, _multipole, *_box);
  }
  for (int i=0; i<Nmpi; i++) {
    sendprunedtree[i] = new OctTree<ndim,ParticleType,TreeCell>
     (_Nleafmax, _thetamaxsqd, _kernrange, _macerror, _gravity_mac, _multipole, *_box);
  }
#endif
}



//=================================================================================================
//  SM2012SphTree::SM2012SphTree
/// SM2012SphTree constructor.  Initialises various variables.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
SM2012SphTree<ndim,ParticleType,TreeCell>::SM2012SphTree
 (int _Nleafmax, int _Nmpi, int _pruning_level_min, int _pruning_level_max, FLOAT _thetamaxsqd,
  FLOAT _kernrange, FLOAT _macerror, string _gravity_mac, string _multipole,
  DomainBox<ndim>* _box, SmoothingKernel<ndim>* _kern, CodeTiming* _timing):
 NeighbourSearch<ndim>(_kernrange, _box, _kern, _timing),
 SphTree<ndim,ParticleType,TreeCell>
  (_Nleafmax, _Nmpi, _pruning_level_min, _pruning_level_max, _thetamaxsqd,
   _kernrange, _macerror, _gravity_mac, _multipole, _box, _kern, _timing)
{
}



//=================================================================================================
//  SM2012SphTree::~SM2012SphTree
/// SM2012SphTree destructor.  Deallocates tree memory upon object destruction.
//=================================================================================================
/*template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
SM2012SphTree<ndim,ParticleType,TreeCell>::~SM2012SphTree()
{
  if (tree->allocated_tree) {
    this->DeallocateMemory();
    tree->DeallocateTreeMemory();
  }
}*/



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
