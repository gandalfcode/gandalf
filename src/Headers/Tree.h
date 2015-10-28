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
#include "CodeTiming.h"
#include "Constants.h"
#include "DomainBox.h"
#include "Hydrodynamics.h"
#include "InlineFuncs.h"

#include "Nbody.h"
#include "Parameters.h"
#include "Particle.h"
#include "Precision.h"
#include "SmoothingKernel.h"
#ifdef MPI_PARALLEL
#include "MpiNode.h"
#endif
using namespace std;



//=================================================================================================
//  Struct TreeCellBase
/// Base tree cell data structure which contains all data elements common to all trees.
//=================================================================================================
template <int ndim>
struct TreeCellBase {
  int cnext;                           ///< i.d. of next cell if not opened
  int copen;                           ///< i.d. of first child cell
  int id;                              ///< Cell id
  int level;                           ///< Level of cell on tree
  int ifirst;                          ///< i.d. of first particle in cell
  int ilast;                           ///< i.d. of last particle in cell
  int N;                               ///< No. of particles in cell
  int Nactive;                         ///< No. of active particles in cell
  int cexit[2][ndim];                  ///< Left and right exit cells (per dim)
  FLOAT cdistsqd;                      ///< Minimum distance to use COM values
  FLOAT mac;                           ///< Multipole-opening criterion value
  FLOAT bbmin[ndim];                   ///< Minimum extent of bounding box
  FLOAT bbmax[ndim];                   ///< Maximum extent of bounding box
  FLOAT hboxmin[ndim];                 ///< Minimum extent of bounding box
  FLOAT hboxmax[ndim];                 ///< Maximum extent of bounding box
  FLOAT vboxmin[ndim];                 ///< ..
  FLOAT vboxmax[ndim];                 ///< ..
  FLOAT rcell[ndim];                   ///< Geometric centre of cell bounding box
  FLOAT r[ndim];                       ///< Position of cell COM
  FLOAT v[ndim];                       ///< Velocity of cell COM
  FLOAT m;                             ///< Mass contained in cell
  FLOAT rmax;                          ///< Radius of bounding sphere
  FLOAT hmax;                          ///< Maximum smoothing length inside cell
  FLOAT drmaxdt;                       ///< Rate of change of bounding sphere
  FLOAT dhmaxdt;                       ///< Rate of change of maximum h
  FLOAT q[5];                          ///< Quadrupole moment tensor
#ifdef MPI_PARALLEL
  double worktot;                      ///< Total work in cell
#endif
};



//=================================================================================================
//  Class TreeBase
/// \brief   Tree base class.
/// \details Tree base class.
/// \author  D. A. Hubber
/// \date    21/09/2015
/// \modified: 21/10/2015, Richard Booth
///            Added the interface implemented by Tree class, where straightforward.
///            Gravity and ActiveCells need doing.
//=================================================================================================
template <int ndim>
class TreeBase
{
 public:
	TreeBase() {};
	virtual ~TreeBase() { } ;

	virtual int MaxNumPartInLeafCell() const = 0 ;
	virtual FLOAT MaxKernelRange() const = 0 ;
	virtual int MaxNumCells() const = 0 ;

	virtual void ExtrapolateCellProperties(const FLOAT) = 0 ;
	virtual int ComputeActiveParticleList(TreeCellBase<ndim> &, Particle<ndim> *, int *) = 0 ;
	virtual int ComputeActiveCellPointers(TreeCellBase<ndim> **celllist) = 0 ;
	virtual int ComputeGatherNeighbourList(const Particle<ndim> *, const FLOAT *,
	                                       const FLOAT, const int, int &, int *) = 0 ;
	virtual int ComputeGatherNeighbourList(const TreeCellBase<ndim> &, const Particle<ndim> *,
	                                       const FLOAT, const int, int &, int *) = 0 ;
	virtual int ComputeNeighbourList(const TreeCellBase<ndim> &, const Particle<ndim> *,
	                                 const int, int &, int *, Particle<ndim> *) = 0 ;
	virtual int ComputePeriodicNeighbourList(const TreeCellBase<ndim> &, const Particle<ndim> *,
	                                         const DomainBox<ndim> &, const int, int &,
	                                         int *, Particle<ndim> *) = 0 ;

	/* TODO: Members that need more work to be part of a common interface
	 * 		 Fix: Use a proxy class to hold multipole data, rather than returning the cells.
	 *
	virtual int ComputeActiveCellList(TreeCell<ndim> *) ;
	virtual int ComputeGravityInteractionList(const TreeCell<ndim> &, const Particle<ndim> *,
	                                    const FLOAT, const int, const int, int &, int &, int &, int &,
	                                    int *, int *, int *, TreeCell<ndim> *, Particle<ndim> *);
	virtual  int ComputePeriodicGravityInteractionList(const TreeCell<ndim> &, const Particle<ndim> *,
	                                            const DomainBox<ndim> &, const FLOAT, const int,
	                                            const int, int &, int &, int &, int &, int *, int *,
	                                            int *, TreeCell<ndim> *, Particle<ndim> *);
	virtual  int ComputeStarGravityInteractionList(const NbodyParticle<ndim> *, const FLOAT, const int,
												   const int, const int, int &, int &, int &, int *, int *,
	                                               TreeCell<ndim> *, Particle<ndim> *);
	*/

	virtual int FindLeafCell(const FLOAT *) = 0;

	//-----------------------------------------------------------------------------------------------
	// virtual void BuildTree(const int, const int, const int, const int,
	//                        const FLOAT, ParticleType<ndim> *) = 0;
	virtual void AllocateTreeMemory(void) = 0;
	virtual void DeallocateTreeMemory(void) = 0;
	//virtual void UpdateAllHmaxValues(Particle<ndim> *) = 0;
	//virtual void UpdateActiveParticleCounters(Particle<ndim> *) = 0;

#if defined(VERIFY_ALL)
    virtual void ValidateTree(Particle<ndim> *) = 0;
#endif

};



//=================================================================================================
//  Class Tree
/// \brief   Generic Tree class for partitioning particles spatially.
/// \details Generic Tree class for partitioning particles spatially.  Used for (i) computing
///          near-neighbour lists for hydro forces; (ii) computing gravitational forces using an
///          efficient MAC + multipole expansion; (iii) computing radiation properties of all
///          particles via the TreeRay algorithm.
/// \author  D. A. Hubber
/// \date    12/09/2014
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
class Tree : public TreeBase<ndim>
{
 public:


  Tree(int _Nleafmax, FLOAT _thetamaxsqd, FLOAT _kernrange, FLOAT _macerror,
       string _gravity_mac, string _multipole) :
    gravity_mac(_gravity_mac), multipole(_multipole), Nleafmax(_Nleafmax),
    invthetamaxsqd(1.0/_thetamaxsqd), kernrange(_kernrange), macerror(_macerror),
    theta(sqrt(_thetamaxsqd)), thetamaxsqd(_thetamaxsqd) {};

  virtual ~Tree() { } ;

  //-----------------------------------------------------------------------------------------------
  int MaxNumPartInLeafCell() const
  { return Nleafmax ; }
  virtual FLOAT MaxKernelRange() const
  { return kernrange ; }
  virtual int MaxNumCells() const
  { return gtot ; }

  virtual void ExtrapolateCellProperties(const FLOAT);
  bool BoxOverlap(const FLOAT box1min[ndim], const FLOAT box1max[ndim],
		          const FLOAT box2min[ndim], const FLOAT box2max[ndim])
  {
    using ::BoxOverlap ;
    return BoxOverlap(ndim, box1min, box1max, box2min, box2max) ;
  }

  int ComputeActiveParticleList(TreeCellBase<ndim> &, Particle<ndim> *, int *);
  int ComputeActiveCellList(TreeCell<ndim> *);
  int ComputeActiveCellPointers(TreeCellBase<ndim> **celllist);
  int ComputeGatherNeighbourList(const Particle<ndim> *, const FLOAT *,
                                 const FLOAT, const int, int &, int *);
  int ComputeGatherNeighbourList(const TreeCellBase<ndim> &, const Particle<ndim> *,
                                 const FLOAT, const int, int &, int *);
  int ComputeNeighbourList(const TreeCellBase<ndim> &, const Particle<ndim> *,
                           const int, int &, int *, Particle<ndim> *);
  int ComputePeriodicNeighbourList(const TreeCellBase<ndim> &, const Particle<ndim> *,
                                   const DomainBox<ndim> &, const int, int &,
                                   int *, Particle<ndim> *);
  int ComputeGravityInteractionList(const TreeCell<ndim> &, const Particle<ndim> *,
                                    const FLOAT, const int, const int, int &, int &, int &, int &,
                                    int *, int *, int *, TreeCell<ndim> *, Particle<ndim> *);
  int ComputePeriodicGravityInteractionList(const TreeCell<ndim> &, const Particle<ndim> *,
                                            const DomainBox<ndim> &, const FLOAT, const int,
                                            const int, int &, int &, int &, int &, int *, int *,
                                            int *, TreeCell<ndim> *, Particle<ndim> *);
  int ComputeStarGravityInteractionList(const NbodyParticle<ndim> *, const FLOAT, const int,
                                        const int, const int, int &, int &, int &, int *, int *,
                                        TreeCell<ndim> *, Particle<ndim> *);
  virtual int FindLeafCell(const FLOAT *);
#ifdef MPI_PARALLEL
  int CreatePrunedTreeForMpiNode(const MpiNode<ndim> &, const DomainBox<ndim> &, const FLOAT,
                                 const bool, const int, const int, const int, TreeCell<ndim> *);
  int ComputeDistantGravityInteractionList(const TreeCell<ndim> *, const DomainBox<ndim> &,
                                           const FLOAT, const int, int, TreeCell<ndim> *);
  bool ComputeHydroTreeCellOverlap(const TreeCell<ndim> *, const DomainBox<ndim> &);
  FLOAT ComputeWorkInBox(const FLOAT *, const FLOAT *);
#endif


  //-----------------------------------------------------------------------------------------------

  virtual void BuildTree(const int, const int, const int, const int,
                         const FLOAT, ParticleType<ndim> *) = 0;
  void UpdateAllHmaxValues(ParticleType<ndim>* sphdata)
  { UpdateHmaxValues(celldata[0], sphdata) ; }
  virtual void UpdateHmaxValues(TreeCell<ndim> &, ParticleType<ndim> *) = 0;
  virtual void StockTree(TreeCell<ndim> &, ParticleType<ndim> *) = 0 ;
  virtual void UpdateActiveParticleCounters(ParticleType<ndim> *) = 0;

#ifdef MPI_PARALLEL
  virtual void UpdateWorkCounters(TreeCell<ndim> &) = 0;
  virtual int GetMaxCellNumber(const int) = 0;
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
