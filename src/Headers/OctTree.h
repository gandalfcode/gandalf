//=================================================================================================
//  OctTree.h
//  Header file containing class definitions for constructing and updating the octal tree.
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


#ifndef _OCT_TREE_H_
#define _OCT_TREE_H_

#include <assert.h>
#include <iostream>
#include <math.h>
#include <string>
#include <vector>
#include "Precision.h"
#include "Constants.h"
#include "CodeTiming.h"
#include "SmoothingKernel.h"
#include "Particle.h"
#include "Sph.h"
#include "Nbody.h"
#include "DomainBox.h"
#include "Parameters.h"
#include "Tree.h"
using namespace std;


static const int Noctchild = 8;
static const int nfreq=1;


//=================================================================================================
//  Struct OctTreeCell
/// Octal-tree cell data structure
//=================================================================================================
template <int ndim>
struct OctTreeCell : public TreeCellBase<ndim> {
  FLOAT rcentre[ndim];                 ///< Centre of cell when creating tree structure
#ifdef MPI_PARALLEL
  int c1;                              ///< Added just to make it compile - do it properly!!!!
  int c2;                              ///< ..
#endif
};



//=================================================================================================
//  Struct TreeRayCell
/// TreeRay cell data structure
//=================================================================================================
template <int ndim>
struct TreeRayCell : public OctTreeCell<ndim>
{
  FLOAT volume;                               ///< Volume of cell
  FLOAT srcEUV[nfreq];                        ///< Source function of EUV radiation
  FLOAT erdEUVold[nfreq];                     ///< Old radiation energy density (for iteration)
  FLOAT erdEUV[nfreq];                        ///< Radiation energy density
};



//=================================================================================================
//  Class OctTree
/// \brief   Class containing Barnes-Hut octal tree.
/// \details Class containing Barnes-Hut octal tree.
/// \author  D. A. Hubber
/// \date    17/08/2014
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
class OctTree : public Tree<ndim,ParticleType,TreeCell>
{
 public:

  using Tree<ndim,ParticleType,TreeCell>::allocated_tree;
  using Tree<ndim,ParticleType,TreeCell>::celldata;
  using Tree<ndim,ParticleType,TreeCell>::gravity_mac;
  using Tree<ndim,ParticleType,TreeCell>::gtot;
  using Tree<ndim,ParticleType,TreeCell>::ids;
  using Tree<ndim,ParticleType,TreeCell>::ifirst;
  using Tree<ndim,ParticleType,TreeCell>::ilast;
  using Tree<ndim,ParticleType,TreeCell>::inext;
  using Tree<ndim,ParticleType,TreeCell>::lmax;
  using Tree<ndim,ParticleType,TreeCell>::ltot;
  using Tree<ndim,ParticleType,TreeCell>::ltot_old;
  using Tree<ndim,ParticleType,TreeCell>::multipole;
  using Tree<ndim,ParticleType,TreeCell>::Ncell;
  using Tree<ndim,ParticleType,TreeCell>::Ncellmax;
  using Tree<ndim,ParticleType,TreeCell>::Ncellmaxold;
  using Tree<ndim,ParticleType,TreeCell>::Nleafmax;
  using Tree<ndim,ParticleType,TreeCell>::Nthreads;
  using Tree<ndim,ParticleType,TreeCell>::Ntot;
  using Tree<ndim,ParticleType,TreeCell>::Ntotold;
  using Tree<ndim,ParticleType,TreeCell>::Ntotmax;
  using Tree<ndim,ParticleType,TreeCell>::Ntotmaxold;
  using Tree<ndim,ParticleType,TreeCell>::macerror;
  using Tree<ndim,ParticleType,TreeCell>::hmax;
  using Tree<ndim,ParticleType,TreeCell>::kernrange;
  using Tree<ndim,ParticleType,TreeCell>::theta;
  using Tree<ndim,ParticleType,TreeCell>::thetamaxsqd;
  using Tree<ndim,ParticleType,TreeCell>::invthetamaxsqd;
#ifdef MPI_PARALLEL
  using Tree<ndim,ParticleType,TreeCell>::Nimportedcell;
  using Tree<ndim,ParticleType,TreeCell>::Ncelltot;
#endif


  // Constructor and destructor
  //-----------------------------------------------------------------------------------------------
  OctTree(int, FLOAT, FLOAT, FLOAT, string, string);
  ~OctTree();


  //-----------------------------------------------------------------------------------------------
  void BuildTree(const int, const int, const int, const int, const FLOAT, ParticleType<ndim> *);
  void AllocateTreeMemory(void);
  void DeallocateTreeMemory(void);
  //bool BoxOverlap(const FLOAT *, const FLOAT *, const FLOAT *, const FLOAT *);
  //void ExtrapolateCellProperties(FLOAT);
  void StockTree(TreeCell<ndim> &, ParticleType<ndim> *);
  void UpdateHmaxValues(TreeCell<ndim> &, ParticleType<ndim> *);
  void UpdateActiveParticleCounters(ParticleType<ndim> *);
#ifdef MPI_PARALLEL
  void UpdateWorkCounters(TreeCell<ndim> &) {};
  int GetMaxCellNumber(const int _level) {return pow(pow(2,ndim),_level);};
#endif
#if defined(VERIFY_ALL)
  void ValidateTree(ParticleType<ndim> *);
#endif


  // Additional variables for octal tree
  //-----------------------------------------------------------------------------------------------
  int *firstCell;                      ///< Array containing first cells in each tree level
  int *lastCell;                       ///<   "        "     last       "          "
  FLOAT rootCellSize;                  ///< Length of side of root cell cube

};
#endif
