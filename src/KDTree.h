//=============================================================================
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
//=============================================================================


#ifndef _KD_TREE_H_
#define _KD_TREE_H_

#include <assert.h>
#include <iostream>
#include <string>
#include <vector>
#include "Precision.h"
#include "Constants.h"
#include "CodeTiming.h"
#include "SphKernel.h"
#include "SphParticle.h"
#include "Sph.h"
#include "Nbody.h"
#include "DomainBox.h"
#include "Parameters.h"
using namespace std;



//=============================================================================
//  Struct KDTreeCell
/// KD-tree cell data structure
//=============================================================================
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
  int N;                            ///< ..
  int Nactive;                      ///< ..
  int cexit[2][ndim];               ///< Left and right exit cells (per dim)
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
};



//=============================================================================
//  Class KDTree
/// \brief   Class containing binary tree
/// \details Binary tree data structure used for efficient neighbour searching 
///          and computation of gravitational forces
/// \author  D. A. Hubber, O. Lomax, A. P. Whitworth
/// \date    08/01/2014
//=============================================================================
template <int ndim, template<int> class ParticleType>
class KDTree
{
 public:

  KDTree(int, FLOAT, FLOAT, FLOAT, string, string);
  ~KDTree();


  //---------------------------------------------------------------------------
  void BuildTree(int, int, ParticleType<ndim> *, FLOAT);
  void AllocateTreeMemory(void);
  void DeallocateTreeMemory(void);
  bool BoxOverlap(const FLOAT *, const FLOAT *, const FLOAT *, const FLOAT *);
  void ComputeTreeSize(void);
  void CreateTreeStructure(void);
  void DivideTreeCell(int, int, ParticleType<ndim> *, KDTreeCell<ndim> &);
  //void ExtrapolateCellProperties(KDTreeCell<ndim> &, FLOAT);
  void ExtrapolateCellProperties(FLOAT);
  FLOAT QuickSelect(int, int, int, int, ParticleType<ndim> *);
  void StockTree(KDTreeCell<ndim> &, ParticleType<ndim> *);
  void StockCellProperties(KDTreeCell<ndim> &, ParticleType<ndim> *);
  void UpdateHmaxValues(KDTreeCell<ndim> &, ParticleType<ndim> *);
  void UpdateActiveParticleCounters(ParticleType<ndim> *);
  int ComputeActiveCellList(KDTreeCell<ndim> **);
  int ComputeActiveParticleList(KDTreeCell<ndim> *, 
                                ParticleType<ndim> *, int *);
  int ComputeGatherNeighbourList(FLOAT *, FLOAT, int, int *, 
                                 ParticleType<ndim> *);
  int ComputeGatherNeighbourList(const KDTreeCell<ndim> *, const int, int *, 
                                 const FLOAT, const ParticleType<ndim> *);
  int ComputeNeighbourList(const KDTreeCell<ndim> *, const int, int *, 
                           const ParticleType<ndim> *);
  int ComputeGravityInteractionList(KDTreeCell<ndim> *, FLOAT, int, int,  
                                    int, int &, int &, int &, int *, int *,
                                    KDTreeCell<ndim> **, 
				    ParticleType<ndim> *);
  int ComputeStarGravityInteractionList(NbodyParticle<ndim> *, FLOAT, int, int,
					int, int &, int &, int &, int *, int *,
					KDTreeCell<ndim> **, 
					ParticleType<ndim> *);
  void ComputeCellMonopoleForces(FLOAT &, FLOAT *, FLOAT *, int, 
				 KDTreeCell<ndim> **);
  void ComputeCellQuadrupoleForces(FLOAT &, FLOAT *, FLOAT *, int, 
				   KDTreeCell<ndim> **);
  void ComputeFastMonopoleForces(int, int, KDTreeCell<ndim> **, 
				 KDTreeCell<ndim> *, ParticleType<ndim> *);
#if defined(VERIFY_ALL)
  void ValidateTree(ParticleType<ndim> *);
#endif

  // Additional variables for KD-tree
  //---------------------------------------------------------------------------
  string gravity_mac;               ///< Multipole-acceptance criteria for tree
  string multipole;                 ///< Multipole-order for cell gravity
  bool allocated_tree;              ///< Are grid arrays allocated?
  int gmax;                         ///< Max. no. of grid/leaf cells
  int gtot;                         ///< Total number of grid/leaf cells
  int lmax;                         ///< Max. no. of levels
  int ltot;                         ///< Total number of levels in tree
  int ltot_old;                     ///< Prev. value of ltot
  int Ncell;                        ///< Current no. of grid cells
  int Ncellmax;                     ///< Max. allowed no. of grid cells
  int Nlevel;                       ///< ""
  int Nleafmax;                     ///< Max. number of particles per leaf cell
  int Nlistmax;                     ///< Max. length of neighbour list
  int Nthreads;                     ///< No. of OpenMP threads
  int Ntot;                         ///< No. of current points in list
  int Ntotold;                      ///< Prev. no. of particles
  int Ntotmax;                      ///< Max. no. of points in list
  int Ntotmaxold;                   ///< Old value of Ntotmax
  int *g2c;                         ///< i.d. of leaf(grid) cells
  int *ids;                         ///< Particle ids
  int *inext;                       ///< Linked list for grid search
  FLOAT macerror;                   ///< Error tolerance for gravity tree-MAC
  FLOAT hmax;                       ///< Store hmax in the tree
  FLOAT kernrange;                  ///< Extent of employed kernel
  FLOAT theta;                      ///< Geometric opening angle
  FLOAT thetamaxsqd;                ///< Geometric opening angle squared
  FLOAT invthetamaxsqd;             ///< 1 / thetamaxsqd
  KDTreeCell<ndim> *kdcell;         ///< Main tree array

};
#endif
