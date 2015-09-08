//=================================================================================================
//  Tree.h
//  Header file containing class definition for generic spatial-tree data structure.
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


#ifndef _TREE_H_
#define _TREE_H_


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
#ifdef MPI_PARALLEL
#include "MpiNode.h"
#endif
using namespace std;



//=================================================================================================
//  Class Tree
/// \brief   Generic Tree class for partitioning particles spatially.
/// \details Generic Tree class for partitioning particles spatially.
/// \author  D. A. Hubber
/// \date    12/09/2014
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
class Tree
{
 public:

  //Tree() {};
  Tree(int _Nleafmax, FLOAT _thetamaxsqd, FLOAT _kernrange, FLOAT _macerror,
       string _gravity_mac, string _multipole):
    gravity_mac(_gravity_mac), multipole(_multipole), Nleafmax(_Nleafmax),
    invthetamaxsqd(1.0/_thetamaxsqd), kernrange(_kernrange), macerror(_macerror),
    theta(sqrt(_thetamaxsqd)), thetamaxsqd(_thetamaxsqd) {};


  //-----------------------------------------------------------------------------------------------
  void ExtrapolateCellProperties(const FLOAT);
  bool BoxOverlap(const FLOAT *, const FLOAT *, const FLOAT *, const FLOAT *);
  int ComputeActiveParticleList(TreeCell<ndim> &, ParticleType<ndim> *, int *);
  int ComputeActiveCellList(TreeCell<ndim> *);
  int ComputeActiveCellPointers(TreeCell<ndim> **celllist);
  int ComputeGatherNeighbourList(const ParticleType<ndim> *, const FLOAT *,
                                 const FLOAT, const int, int &, int *);
  int ComputeGatherNeighbourList(const TreeCell<ndim> &, const ParticleType<ndim> *,
                                 const FLOAT, const int, int &, int *);
  int ComputeNeighbourList(const TreeCell<ndim> &, const ParticleType<ndim> *,
                           const int, int &, int *, ParticleType<ndim> *);
  int ComputePeriodicNeighbourList(const TreeCell<ndim> &, const ParticleType<ndim> *,
                                   const DomainBox<ndim> &, const int, int &,
                                   int *, ParticleType<ndim> *);
  int ComputeGravityInteractionList(const TreeCell<ndim> &, const ParticleType<ndim> *,
                                    const FLOAT, const int, const int, int &, int &, int &, int &,
                                    int *, int *, int *, TreeCell<ndim> *, ParticleType<ndim> *);
  int ComputePeriodicGravityInteractionList(const TreeCell<ndim> &, const ParticleType<ndim> *,
                                            const DomainBox<ndim> &, const FLOAT, const int,
                                            const int, int &, int &, int &, int &, int *, int *,
                                            int *, TreeCell<ndim> *, ParticleType<ndim> *);
  int ComputeStarGravityInteractionList(const NbodyParticle<ndim> *, const FLOAT, const int,
                                        const int, const int, int &, int &, int &, int *, int *,
                                        TreeCell<ndim> *, ParticleType<ndim> *);
#ifdef MPI_PARALLEL
  int CreatePrunedTreeForMpiNode(const MpiNode<ndim> &, const DomainBox<ndim> &, const FLOAT,
                                 const bool, const int, const int, const int, TreeCell<ndim> *);
  int ComputeDistantGravityInteractionList(const TreeCell<ndim> *, const FLOAT,
                                           const int, int, TreeCell<ndim> *);
  bool ComputeHydroTreeCellOverlap(const TreeCell<ndim> *);
  FLOAT ComputeWorkInBox(const FLOAT *, const FLOAT *);
#endif


  //-----------------------------------------------------------------------------------------------
  virtual void BuildTree(const int, const int, const int, const int,
                         const FLOAT, ParticleType<ndim> *) = 0;
  virtual void AllocateTreeMemory(void) = 0;
  virtual void DeallocateTreeMemory(void) = 0;
  virtual void StockTree(TreeCell<ndim> &, ParticleType<ndim> *) = 0;
  virtual void UpdateHmaxValues(TreeCell<ndim> &, ParticleType<ndim> *) = 0;
  virtual void UpdateActiveParticleCounters(ParticleType<ndim> *) = 0;
#ifdef MPI_PARALLEL
  virtual void UpdateWorkCounters(TreeCell<ndim> &) = 0;
  virtual int GetMaxCellNumber(const int) = 0;
#endif
#if defined(VERIFY_ALL)
  virtual void ValidateTree(ParticleType<ndim> *) = 0;
#endif


  // Const variables for tree class
  //-----------------------------------------------------------------------------------------------
  const string gravity_mac;            ///< Multipole-acceptance criteria for tree
  const string multipole;              ///< Multipole-order for cell gravity
  const int Nleafmax;                  ///< Max. number of particles per leaf cell
  const FLOAT invthetamaxsqd;          ///< 1 / thetamaxsqd
  const FLOAT kernrange;               ///< Extent of employed kernel
  const FLOAT macerror;                ///< Error tolerance for gravity tree-MAC
  const FLOAT theta;                   ///< Geometric opening angle
  const FLOAT thetamaxsqd;             ///< Geometric opening angle squared

  // Additional variables for tree class
  //-----------------------------------------------------------------------------------------------
  bool allocated_tree;                 ///< Are grid arrays allocated?
  int gmax;                            ///< Max. no. of grid/leaf cells
  int gtot;                            ///< Total number of grid/leaf cells
  int ifirst;                          ///< i.d. of first particle in tree
  int ilast;                           ///< i.d. of last particle in tree
  int lactive;                         ///< ..
  int lmax;                            ///< Max. no. of levels
  int ltot;                            ///< Total number of levels in tree
  int ltot_old;                        ///< Prev. value of ltot
  int Ncell;                           ///< Current no. of grid cells
  int Ncellmax;                        ///< Max. allowed no. of grid cells
  int Ncellmaxold;                     ///< Old value of Ncellmax
  int Nthreads;                        ///< No. of OpenMP threads
  int Ntot;                            ///< No. of current points in list
  int Ntotmax;                         ///< Max. no. of points in list
  int Ntotmaxold;                      ///< Old value of Ntotmax
  int Ntotold;                         ///< Prev. no. of particles
  FLOAT hmax;                          ///< Store hmax in the tree
  int *g2c;                            ///< i.d. of leaf(grid) cells
  int *ids;                            ///< Particle ids
  int *inext;                          ///< Linked list for grid search
  TreeCell<ndim> *celldata;            ///< Main tree cell data array

#if defined MPI_PARALLEL
  int Nimportedcell;                   ///< No. of imported cells
  int Ncelltot;                        ///< Total number of cells
#endif

};
#endif
