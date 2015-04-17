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


//=================================================================================================
//  Struct OctTreeCell
/// KD-tree cell data structure
//=================================================================================================
template <int ndim>
struct OctTreeCell {
  int copen;                        ///< ..
  int c2g;                          ///< i.d. of tree-cell c/grid-cell g
  int cnext;                        ///< i.d. of next cell if not opened
  int id;                           ///< Cell id
  int k_divide;                     ///< Dimension along which cell is split
  int level;                        ///< Level of cell on tree
  int ifirst;                       ///< i.d. of first particle in cell
  int ilast;                        ///< i.d. of last particle in cell
  int N;                            ///< ..
  int Nactive;                      ///< ..
  int cexit[2][ndim];               ///< Left and right exit cells (per dim)
  int childof[Noctchild];           ///< ..
  FLOAT cdistsqd;                   ///< ..
  FLOAT mac;                        ///< Multipole-opening criterion value
  FLOAT bbmin[ndim];                ///< Minimum extent of bounding box
  FLOAT bbmax[ndim];                ///< Maximum extent of bounding box
  FLOAT hboxmin[ndim];              ///< Minimum extent of bounding box
  FLOAT hboxmax[ndim];              ///< Maximum extent of bounding box
  FLOAT rcell[ndim];                ///< ..
  FLOAT r[ndim];                    ///< Position of cell
  FLOAT v[ndim];                    ///< Position of cell
  FLOAT m;                          ///< Mass contained in cell
  FLOAT rmax;                       ///< Radius of bounding sphere
  FLOAT hmax;                       ///< Maximum smoothing length inside cell
  FLOAT drmaxdt;                    ///< Rate of change of bounding sphere
  FLOAT dhmaxdt;                    ///< Rate of change of maximum h
  FLOAT q[5];                       ///< Quadrupole moment tensor
#ifdef MPI_PARALLEL
  double worktot;                   ///< Total work in cell
  int c1;                           /// Added just to make it compile - do it properly!!!!
  int c2;
#endif
};



//=================================================================================================
//  Class OctTree
/// \brief   Class containing binary tree
/// \details Binary tree data structure used for efficient neighbour searching
///          and computation of gravitational forces
/// \author  D. A. Hubber, O. Lomax, A. P. Whitworth
/// \date    17/08/2014
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
class OctTree : public Tree<ndim,ParticleType,TreeCell>
{
 public:

  using Tree<ndim,ParticleType,TreeCell>::allocated_tree;
  using Tree<ndim,ParticleType,TreeCell>::celldata;
  using Tree<ndim,ParticleType,TreeCell>::gravity_mac;
  using Tree<ndim,ParticleType,TreeCell>::gmax;
  using Tree<ndim,ParticleType,TreeCell>::gtot;
  using Tree<ndim,ParticleType,TreeCell>::g2c;
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
  void BuildTree(int, int, int, int, ParticleType<ndim> *, FLOAT);
  void AllocateTreeMemory(void);
  void DeallocateTreeMemory(void);
  bool BoxOverlap(const FLOAT *, const FLOAT *, const FLOAT *, const FLOAT *);
  void ExtrapolateCellProperties(FLOAT);
  void StockTree(TreeCell<ndim> &, ParticleType<ndim> *);
  void UpdateHmaxValues(TreeCell<ndim> &, ParticleType<ndim> *);
  void UpdateActiveParticleCounters(ParticleType<ndim> *);
  int ComputeActiveCellList(TreeCell<ndim> *);
  int ComputeActiveCellPointers(TreeCell<ndim> **celllist);
  int ComputeActiveParticleList(TreeCell<ndim> &, ParticleType<ndim> *, int *);
  int ComputeGatherNeighbourList(const ParticleType<ndim> *, const FLOAT *,
                                 const FLOAT, const int, int &, int *);
  int ComputeGatherNeighbourList(const TreeCell<ndim> &, const ParticleType<ndim> *,
                                 const FLOAT, const int, int &, int *);
  int ComputeNeighbourList(const TreeCell<ndim> &, const ParticleType<ndim> *,
                           const int, int &, int *, ParticleType<ndim> *);
  int ComputeGravityInteractionList(const TreeCell<ndim> &, const ParticleType<ndim> *,
                                    const FLOAT, const int, const int, int &, int &, int &, int &,
                                    int *, int *, int *, TreeCell<ndim> *, ParticleType<ndim> *);
  int ComputePeriodicGravityInteractionList(const TreeCell<ndim> &, const ParticleType<ndim> *,
                                            const DomainBox<ndim> &, const FLOAT, const int, const int,
                                            int &, int &, int &, int &, int *, int *, int *,
                                            TreeCell<ndim> *, ParticleType<ndim> *) {return 0;}
  int ComputeStarGravityInteractionList(const NbodyParticle<ndim> *, const FLOAT, const int,
                                        const int, const int, int &, int &, int &, int *, int *,
                                        TreeCell<ndim> *, ParticleType<ndim> *);
#ifdef MPI_PARALLEL
  int ComputeDistantGravityInteractionList(const TreeCell<ndim> *, const FLOAT,
                                           const int, int, TreeCell<ndim> *);
  bool ComputeHydroTreeCellOverlap(const TreeCell<ndim> *);
  FLOAT ComputeWorkInBox(const FLOAT *, const FLOAT *) {};
  void UpdateWorkCounters(TreeCell<ndim> &) {};
#endif
#if defined(VERIFY_ALL)
  void ValidateTree(ParticleType<ndim> *);
#endif


  // Additional variables for octal tree
  //-----------------------------------------------------------------------------------------------
  int *firstCell;
  int *lastCell;

};
#endif
