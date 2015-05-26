//=================================================================================================
//  KDTree.h
//  Header file containing class definitions for constructing and updating
//  the KD-tree.
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
struct KDTreeCell {
  int c1;                           ///< First child cell
  int c2;                           ///< Second child cell
  int c2g;                          ///< i.d. of tree-cell c/grid-cell g
  int cnext;                        ///< i.d. of next cell if not opened
  int id;                           ///< Cell id
  int k_divide;                     ///< Dimension along which cell is split
  int level;                        ///< Level of cell on tree
  int ifirst;                       ///< i.d. of first particle in cell
  int ilast;                        ///< i.d. of last particle in cell
  int N;                            ///< No. of particles in cell
  int Nactive;                      ///< No. of active particles in cell
  int cexit[2][ndim];               ///< Left and right exit cells (per dim)
  FLOAT cdistsqd;                   ///< Minimum distance to use COM values
  FLOAT mac;                        ///< Multipole-opening criterion value
  FLOAT bbmin[ndim];                ///< Minimum extent of bounding box
  FLOAT bbmax[ndim];                ///< Maximum extent of bounding box
  FLOAT hboxmin[ndim];              ///< Minimum extent of bounding box
  FLOAT hboxmax[ndim];              ///< Maximum extent of bounding box
  FLOAT vboxmin[ndim];
  FLOAT vboxmax[ndim];
  FLOAT rcell[ndim];                ///< Geometric centre of cell bounding box
  FLOAT r[ndim];                    ///< Position of cell COM
  FLOAT v[ndim];                    ///< Velocity of cell COM
  FLOAT m;                          ///< Mass contained in cell
  FLOAT rmax;                       ///< Radius of bounding sphere
  FLOAT hmax;                       ///< Maximum smoothing length inside cell
  FLOAT drmaxdt;                    ///< Rate of change of bounding sphere
  FLOAT dhmaxdt;                    ///< Rate of change of maximum h
  FLOAT q[5];                       ///< Quadrupole moment tensor
#ifdef MPI_PARALLEL
  double worktot;                   ///< Total work in cell
#endif
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
  void BuildTree(int, int, int, int, ParticleType<ndim> *, FLOAT);
  void AllocateTreeMemory(void);
  void DeallocateTreeMemory(void);
  bool BoxOverlap(const FLOAT *, const FLOAT *, const FLOAT *, const FLOAT *);
  void ComputeTreeSize(void);
  void CreateTreeStructure(void);
  void DivideTreeCell(int, int, ParticleType<ndim> *, TreeCell<ndim> &);
  void ExtrapolateCellProperties(FLOAT);
  FLOAT QuickSelect(int, int, int, int, ParticleType<ndim> *);
  FLOAT QuickSelectSort(int, int, int, int, ParticleType<ndim> *);
  void StockTree(TreeCell<ndim> &, ParticleType<ndim> *);
  void StockCellProperties(TreeCell<ndim> &, ParticleType<ndim> *);
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
                           const int, int & ,int *, ParticleType<ndim> *);
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
  int ComputeDistantGravityInteractionList(const TreeCell<ndim> *, const FLOAT, const int,
                                           int, TreeCell<ndim> *);
  bool ComputeHydroTreeCellOverlap(const TreeCell<ndim> *);
  FLOAT ComputeWorkInBox(const FLOAT *, const FLOAT *);
  void UpdateWorkCounters(TreeCell<ndim> &);
#endif
#if defined(VERIFY_ALL)
  void ValidateTree(ParticleType<ndim> *);
#endif


  // Additional KD-tree variables
  //-----------------------------------------------------------------------------------------------
  int gactive;                         ///< No. of active leaf/grid cells
  int lactive;                         ///< Chosen level for computing active particle loops

};
#endif
