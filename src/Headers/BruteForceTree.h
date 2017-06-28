//=================================================================================================
//  BruteForceTree.h
//  Header file containing class definitions for constructing and updating the brute force tree.
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

#ifndef _BF_TREE_H_
#define _BF_TREE_H_

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


//=================================================================================================
//  Struct BFTreeCell
/// Octal-tree cell data structure
//=================================================================================================
template <int ndim>
struct BruteForceTreeCell : public TreeCellBase<ndim> {
#ifdef MPI_PARALLEL
  typedef TreeCommunicationHandler<ndim> HandlerType;
#endif
};



//=================================================================================================
//  Class BruteForceTree
/// \brief   Class containing Brute force wrapper tree.
/// \details Class containing Brute force wrapper tree.
/// \author  D. A. Hubber
/// \date    17/08/2014
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
class BruteForceTree : public Tree<ndim,ParticleType,TreeCell>
{
 public:

  using Tree<ndim,ParticleType,TreeCell>::allocated_tree;
  using Tree<ndim,ParticleType,TreeCell>::celldata;
  using Tree<ndim,ParticleType,TreeCell>::gmax;
  using Tree<ndim,ParticleType,TreeCell>::gravity_mac;
  using Tree<ndim,ParticleType,TreeCell>::gtot;
  using Tree<ndim,ParticleType,TreeCell>::g2c;
  using Tree<ndim,ParticleType,TreeCell>::hmax;
  using Tree<ndim,ParticleType,TreeCell>::ifirst;
  using Tree<ndim,ParticleType,TreeCell>::ilast;
  using Tree<ndim,ParticleType,TreeCell>::invthetamaxsqd;
  using Tree<ndim,ParticleType,TreeCell>::kernrange;
  using Tree<ndim,ParticleType,TreeCell>::lmax;
  using Tree<ndim,ParticleType,TreeCell>::ltot;
  using Tree<ndim,ParticleType,TreeCell>::ltot_old;
  using Tree<ndim,ParticleType,TreeCell>::macerror;
  using Tree<ndim,ParticleType,TreeCell>::multipole;
  using Tree<ndim,ParticleType,TreeCell>::Ncell;
  using Tree<ndim,ParticleType,TreeCell>::Ncellmax;
  using Tree<ndim,ParticleType,TreeCell>::Nleafmax;
  using Tree<ndim,ParticleType,TreeCell>::Nthreads;
  using Tree<ndim,ParticleType,TreeCell>::Ntot;
  using Tree<ndim,ParticleType,TreeCell>::Ntotmax;
  using Tree<ndim,ParticleType,TreeCell>::theta;
  using Tree<ndim,ParticleType,TreeCell>::thetamaxsqd;
  using Tree<ndim,ParticleType,TreeCell>::gravmask ;
#ifdef MPI_PARALLEL
  using Tree<ndim,ParticleType,TreeCell>::Nimportedcell;
#endif

  // Constructor and destructor
  //-----------------------------------------------------------------------------------------------
  BruteForceTree(int, FLOAT, FLOAT, FLOAT, string, multipole_method, const DomainBox<ndim>&,
		  	  	 const ParticleTypeRegister& reg, const bool);
  ~BruteForceTree() ;


  //-----------------------------------------------------------------------------------------------
  void BuildTree(const int, const int, const int, const int, const FLOAT, Particle<ndim> *);
  void AllocateTreeMemory(int,int,bool);
  void ReallocateMemory(int,int);
  void DeallocateTreeMemory(void);
  void StockTree(Particle<ndim> *part_gen, bool stock_leaf) {
    ParticleType<ndim>* partdata = reinterpret_cast<ParticleType<ndim>*>(part_gen) ;
    StockTree(celldata[0], partdata, stock_leaf) ;
  }
  void StockTree(TreeCell<ndim>&, ParticleType<ndim> *, bool stock_leaf);
  void StockTreeProperties(TreeCell<ndim> &, ParticleType<ndim> *, bool);
  void UpdateAllHmaxValues(Particle<ndim> *part_gen, bool stock_leaf) {
    ParticleType<ndim>* partdata = reinterpret_cast<ParticleType<ndim>*>(part_gen) ;
    UpdateHmaxValues(celldata[0], partdata, stock_leaf) ;
  }
  void UpdateHmaxValues(TreeCell<ndim>&, ParticleType<ndim> *, bool);
  void UpdateHmaxValuesCell(TreeCell<ndim> &, ParticleType<ndim> *, bool);

  void UpdateActiveParticleCounters(Particle<ndim> *);
#ifdef MPI_PARALLEL
  void UpdateWorkCounters() ;
  int GetMaxCellNumber(const int _level) {return _level == 0 ? 1 : -1 ;}
#endif
#if defined(VERIFY_ALL)
  void ValidateTree(ParticleType<ndim> *);
#endif

};
#endif
