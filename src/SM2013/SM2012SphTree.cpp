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
//  SM2012SphTree::SM2012SphTree
/// SM2012SphTree constructor.  Initialises various variables.
//=================================================================================================
template <int ndim, template<int> class ParticleType>
SM2012SphTree<ndim,ParticleType>::SM2012SphTree
 (string tree_type,
  int _Nleafmax, int _Nmpi, int _pruning_level_min, int _pruning_level_max, FLOAT _thetamaxsqd,
  FLOAT _kernrange, FLOAT _macerror, string _gravity_mac, multipole_method _multipole,
  DomainBox<ndim>* _box, SmoothingKernel<ndim>* _kern, CodeTiming* _timing,
  ParticleTypeRegister& types):
 SphTree<ndim,ParticleType>
  (tree_type, _Nleafmax, _Nmpi, _pruning_level_min, _pruning_level_max, _thetamaxsqd,
   _kernrange, _macerror, _gravity_mac, _multipole, _box, _kern, _timing, types)
{

}



template class SM2012SphTree<1,SM2012SphParticle>;
template class SM2012SphTree<2,SM2012SphParticle>;
template class SM2012SphTree<3,SM2012SphParticle>;


