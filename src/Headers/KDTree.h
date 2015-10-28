//=================================================================================================
//  KDTree.h
//  Header file containing class definitions for constructing and updating the KD-tree.
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


#ifndef _KD_TREE_H_
#define _KD_TREE_H_

#include <assert.h>
#include <iostream>
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



//=================================================================================================
//  Struct KDTreeCell
/// KD-tree cell data structure
//=================================================================================================
template <int ndim>
struct KDTreeCell : public TreeCellBase<ndim> {
  int c1;                           ///< First child cell
  int c2;                           ///< Second child cell
  int c2g;                          ///< i.d. of tree-cell c/grid-cell g
  int k_divide;                     ///< Dimension along which cell is split
};



//=================================================================================================
//  Class KDTree
/// \brief   Class containing binary tree
/// \details Binary tree data structure used for efficient neighbour searching
///          and computation of gravitational forces
/// \author  D. A. Hubber, O. Lomax, A. P. Whitworth
/// \date    08/01/2014
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
class KDTree : public Tree<ndim,ParticleType,TreeCell>
{
 public:

  using Tree<ndim,ParticleType,TreeCell>::allocated_tree;
  using Tree<ndim,ParticleType,TreeCell>::celldata;
  using Tree<ndim,ParticleType,TreeCell>::gmax;
  using Tree<ndim,ParticleType,TreeCell>::gravity_mac;
  using Tree<ndim,ParticleType,TreeCell>::gtot;
  using Tree<ndim,ParticleType,TreeCell>::g2c;
  using Tree<ndim,ParticleType,TreeCell>::hmax;
  using Tree<ndim,ParticleType,TreeCell>::ids;
  using Tree<ndim,ParticleType,TreeCell>::ifirst;
  using Tree<ndim,ParticleType,TreeCell>::ilast;
  using Tree<ndim,ParticleType,TreeCell>::inext;
  using Tree<ndim,ParticleType,TreeCell>::invthetamaxsqd;
  using Tree<ndim,ParticleType,TreeCell>::kernrange;
  using Tree<ndim,ParticleType,TreeCell>::lmax;
  using Tree<ndim,ParticleType,TreeCell>::ltot;
  using Tree<ndim,ParticleType,TreeCell>::ltot_old;
  using Tree<ndim,ParticleType,TreeCell>::macerror;
  using Tree<ndim,ParticleType,TreeCell>::multipole;
  using Tree<ndim,ParticleType,TreeCell>::Ncell;
  using Tree<ndim,ParticleType,TreeCell>::Ncellmax;
  using Tree<ndim,ParticleType,TreeCell>::Ncellmaxold;
  using Tree<ndim,ParticleType,TreeCell>::Nleafmax;
  using Tree<ndim,ParticleType,TreeCell>::Nthreads;
  using Tree<ndim,ParticleType,TreeCell>::Ntot;
  using Tree<ndim,ParticleType,TreeCell>::Ntotmax;
  using Tree<ndim,ParticleType,TreeCell>::Ntotmaxold;
  using Tree<ndim,ParticleType,TreeCell>::Ntotold;
  using Tree<ndim,ParticleType,TreeCell>::theta;
  using Tree<ndim,ParticleType,TreeCell>::thetamaxsqd;
#ifdef MPI_PARALLEL
  using Tree<ndim,ParticleType,TreeCell>::Ncelltot;
  using Tree<ndim,ParticleType,TreeCell>::Nimportedcell;
#endif


  // Constructor and destructor
  //-----------------------------------------------------------------------------------------------
  KDTree(int, FLOAT, FLOAT, FLOAT, string, string);
  ~KDTree();


  //-----------------------------------------------------------------------------------------------
  void BuildTree(const int, const int, const int, const int, const FLOAT, ParticleType<ndim> *);
  void AllocateTreeMemory(void);
  void DeallocateTreeMemory(void);
  void ComputeTreeSize(void);
  void CreateTreeStructure(void);
  void DivideTreeCell(int, int, ParticleType<ndim> *, TreeCell<ndim> &);
  void ExtrapolateCellProperties(const FLOAT);
  FLOAT QuickSelect(int, int, int, int, ParticleType<ndim> *);
  FLOAT QuickSelectSort(int, int, int, int, ParticleType<ndim> *);
  void StockTree(TreeCell<ndim> &, ParticleType<ndim> *);
  void StockCellProperties(TreeCell<ndim> &, ParticleType<ndim> *);
  void UpdateHmaxValues(TreeCell<ndim> &, ParticleType<ndim> *);
  void UpdateActiveParticleCounters(ParticleType<ndim> *);
#ifdef MPI_PARALLEL
  void UpdateWorkCounters(TreeCell<ndim> &);
  int GetMaxCellNumber(const int _level) {return pow(2,_level);};
#endif
#if defined(VERIFY_ALL)
  void ValidateTree(ParticleType<ndim> *);
#endif

};
#endif
