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
#include "TreeCell.h"
#include "CodeTiming.h"
#include "Constants.h"
#include "DomainBox.h"
#include "Hydrodynamics.h"
#include "InlineFuncs.h"
#include "Multipole.h"
#include "Nbody.h"
#include "Parameters.h"
#include "Particle.h"
#include "Precision.h"
#include "SmoothingKernel.h"
#include "NeighbourManagerBase.h"
#ifdef MPI_PARALLEL
#include "MpiNode.h"
template<int ndim> class TreeCommunicationHandler;
#endif
using namespace std;



enum MAC_Type {
  geometric = 0,
  eigenmac  = 1,
  gadget2   = 2
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
protected:
    const DomainBox<ndim>& _domain ;     ///< Whole simulation domain
#ifdef MPI_PARALLEL
	vector<int> Nleaf_indices; ///< Indices of the leaf cells (only used if this is a pruned tree)
	vector<int> Nleaf_indices_inlocal; ///< Indices of the leaf cells in the full local tree
#endif
 public:
	TreeBase(  const DomainBox<ndim>& domain)
    : _domain(domain), Ntot(0), Ncell(0)
#ifdef MPI_PARALLEL
   , Nimportedcell(0), Ncelltot(0), first_stock(true)
#endif
 {};
	virtual ~TreeBase() { } ;

	virtual void ReallocateMemory (int,int) = 0;

	const DomainBox<ndim>& GetDomain() {return _domain;};
	virtual int MaxNumPartInLeafCell() const = 0 ;
	virtual FLOAT MaxKernelRange() const = 0 ;
	virtual int MaxNumCells() const = 0 ;

	virtual void BuildTree(const int, const int, const int, const int,
	                       const FLOAT, Particle<ndim> *) = 0;
	virtual void StockTree(Particle<ndim> *, bool) = 0 ;
	virtual void UpdateActiveParticleCounters(Particle<ndim> *) = 0;


  virtual MAC_Type GetMacType() const  = 0;
  virtual void SetMacType(MAC_Type) = 0;


	virtual int ComputeActiveParticleList(TreeCellBase<ndim> &, Particle<ndim> *, int *) = 0 ;
	virtual int ComputeActiveCellPointers(TreeCellBase<ndim> **celllist) = 0 ;
	virtual int ComputeActiveCellList(vector<TreeCellBase<ndim> >& ) = 0 ;
	virtual int ComputeGatherNeighbourList(const Particle<ndim> *, const FLOAT *,
	                                       const FLOAT, const int, int &, int *) = 0 ;
	virtual void ComputeGatherNeighbourList(const TreeCellBase<ndim> &, const Particle<ndim> *,
			                                const FLOAT,NeighbourManagerBase& neibmanager) = 0 ;
	virtual int ComputeNeighbourList(const TreeCellBase<ndim> &, const Particle<ndim> *,
	                                 const int, int &, int *, Particle<ndim> *) = 0 ;
	virtual void ComputeNeighbourList(const TreeCellBase<ndim> &cell,NeighbourManagerBase& neibmanager)=0;
	virtual void ComputeNeighbourAndGhostList(const TreeCellBase<ndim> &, NeighbourManagerBase&) = 0 ;
	virtual void ComputeGravityInteractionAndGhostList(const TreeCellBase<ndim> &,
	                                                   NeighbourManagerDim<ndim>& neibmanager)=0;
	virtual int ComputeStarGravityInteractionList(const NbodyParticle<ndim> *, const FLOAT, const int,
	                                              const int, const int, int &, int &, int &, int *, int *,
	                                              MultipoleMoment<ndim> *, Particle<ndim> *) = 0;
#if defined(MPI_PARALLEL)
	virtual int ComputeImportedCellList(vector<TreeCellBase<ndim> >& ) = 0;
	int GetNLeafCells() {return Nleaf_indices.size();};
	virtual int CreatePrunedTreeForMpiNode(const MpiNode<ndim> &, const DomainBox<ndim> &, const FLOAT,
	                                       const bool, const int, const int, const int, TreeBase<ndim> *) = 0;
	virtual int ComputeDistantGravityInteractionList(const TreeCellBase<ndim>&, const DomainBox<ndim> &,
	                                                 vector<MultipoleMoment<ndim> >&) = 0;
	virtual  bool ComputeHydroTreeCellOverlap(const TreeCellBase<ndim> *, const DomainBox<ndim> &) = 0;
	virtual  FLOAT ComputeWorkInBox(const FLOAT *, const FLOAT *) = 0;
	virtual void UpdateLeafCells(TreeBase<ndim>*)=0;
	enum direction { to_buffer, from_buffer };
	virtual void CopyLeafCells(vector<char>& buffer, direction)=0;
        virtual void UpdateBoundingBoxData(MpiNode<ndim>& node) = 0;
#endif

	virtual bool ComputeSignalVelocityFromDistantInteractions(const TreeCellBase<ndim>& cell,
	                                                          int Nactive, Particle<ndim>* active_gen,
	                                                          Particle<ndim>* part_gen) = 0;

	virtual int FindLeafCell(const FLOAT *) = 0;

	//-----------------------------------------------------------------------------------------------
	// virtual void BuildTree(const int, const int, const int, const int,
	//                        const FLOAT, ParticleType<ndim> *) = 0;
	virtual void AllocateTreeMemory(int,int,bool) = 0;
	virtual void DeallocateTreeMemory(void) = 0;
	virtual void UpdateAllHmaxValues(Particle<ndim> *, bool stock_leaf) = 0;
	void UpdateAllHmaxValues(Particle<ndim> *part_gen) {
      UpdateAllHmaxValues(part_gen, true);
    }

	virtual void UpdateHmaxLeaf(TreeCellBase<ndim>&, Particle<ndim> *) = 0;
	virtual double GetMaximumSmoothingLength() const = 0 ;
	//virtual void UpdateActiveParticleCounters(Particle<ndim> *) = 0;

	virtual void GenerateBoundaryGhostParticles(const FLOAT, const FLOAT, const int,
	                                            const DomainBox<ndim>&,
	                                            Hydrodynamics<ndim>*) = 0;

#if defined(VERIFY_ALL)
   // virtual void ValidateTree(Particle<ndim> *) = 0;
#endif

#if defined(MPI_PARALLEL)
    virtual int GetMaxCellNumber(const int) = 0;
	virtual void InitialiseCellWorkCounters() = 0 ;
	virtual void UpdateWorkCounters() = 0;

	virtual void AddWorkCost(vector<TreeCellBase<ndim> >&, double twork, int& Nactivetot) = 0;
	virtual int GetTreeCellSize() const = 0 ;
	virtual void FindBoxOverlapParticles(const Box<ndim>&, vector<int>&, const Particle<ndim>*) = 0;
	virtual int FindBoxGhostParticles(const FLOAT, const FLOAT, const Box<ndim> &,
	                                  vector<int> &export_list) = 0;
	virtual int PackParticlesAndCellsForMPITransfer(const vector<int>& celllist,
                                                    vector<int>& cell_ids, vector<int>& part_ids,
                                                    vector<char>& send_buffer,
                                                    Particle<ndim>*) = 0;
	virtual void PackParticlesAndCellsForMPIReturn(int start_part, int Npart,
                                                   int start_cell, int Ncall,
                                                   vector<char>& send_buffer,
                                                   Particle<ndim>* part_gen) = 0;

	virtual void UnpackParticlesAndCellsFromMPITransfer(int offset_part, int Npart,
                                                        int offset_cell, int Ncell,
                                                        const vector<char>& recv_buffer,
                                                        Hydrodynamics<ndim>*) = 0 ;

	virtual void UnpackParticlesAndCellsForMPIReturn(const vector<int>& part_ids,
                                                     const vector<int>& cell_ids,
                                                     vector<char>& recv_buffer,
                                                     Hydrodynamics<ndim>*) = 0;

	virtual int GetSizeOfExportedParticleData(int Nparticles) const = 0 ;
	virtual int GetSizeOfExportedCellData(int Ncell) const = 0 ;
	virtual int GetSizeOfReturnedParticleData(int Nparticles) const = 0 ;
	virtual int GetSizeOfReturnedCellData(int Ncell) const = 0 ;
	virtual const void* GetCellDataPointer() const = 0 ;
	virtual void* GetCellDataPointer() = 0 ;
#endif

	int Ntot;                              ///< No. of current points in list
	int Ncell;                             ///< Current no. of grid cells
	int gmax;                              ///< Max. no. of grid/leaf cells
    int gtot;                              ///< Total number of grid/leaf cells
    int lmax;                              ///< Max. no. of levels
    int ltot;                              ///< Total number of levels in tree
    int ltot_old;                          ///< Prev. value of ltot
#if defined MPI_PARALLEL
	int Nimportedcell;                     ///< No. of imported cells
	int Ncelltot;                          ///< Total number of cells
	bool first_stock;		   ///< Whether this is the first time we are stocking this tree (only relevant if this is a received pruned tree)
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
private:
#ifdef MPI_PARALLEL
	using TreeBase<ndim>::Nleaf_indices;
	using TreeBase<ndim>::Nleaf_indices_inlocal;
	using TreeBase<ndim>::first_stock;
#endif
protected:
	int Ncellmax;                        ///< Max. allowed no. of grid cells
	int Ntotmax;                         ///< Max. no. of points in list
	const bool IAmPruned;				 ///< Whether we are a pruned tree
 public:

  Tree(int _Nleafmax, FLOAT _thetamaxsqd, FLOAT _kernrange, FLOAT _macerror,
       string _gravity_mac, multipole_method _multipole, const DomainBox<ndim>& domain,
       const ParticleTypeRegister& pt_reg, const bool _IAmPruned) :
    	   TreeBase<ndim>(domain),
    gravity_mac(geometric), multipole(_multipole), Nleafmax(_Nleafmax),
    invthetamaxsqd(1.0/_thetamaxsqd), kernrange(_kernrange), macerror(_macerror),
    theta(sqrt(_thetamaxsqd)), thetamaxsqd(_thetamaxsqd),
    gravmask(pt_reg.gravmask), IAmPruned(_IAmPruned)
    {
      if (_gravity_mac == "eigenmac")
        gravity_mac = eigenmac ;
      else if(_gravity_mac == "gadget2")
        gravity_mac = gadget2;
      else if (_gravity_mac != "geometric")
        ExceptionHandler::getIstance().raise
        ("Error: gravity_mac type, " + _gravity_mac + ", not recognised");
    };

  virtual ~Tree() { } ;

  virtual void ReallocateMemory (int,int) = 0;

  //-----------------------------------------------------------------------------------------------
  int MaxNumPartInLeafCell() const
  { return Nleafmax ; }
  virtual FLOAT MaxKernelRange() const
  { return kernrange ; }
  virtual int MaxNumCells() const
  { return gtot ; }

  bool BoxOverlap(const FLOAT box1min[ndim], const FLOAT box1max[ndim],
		          const FLOAT box2min[ndim], const FLOAT box2max[ndim])
  {
    using ::BoxOverlap ;
    return BoxOverlap(ndim, box1min, box1max, box2min, box2max) ;
  }

  MAC_Type GetMacType() const {
    return gravity_mac ;
  }
  void SetMacType(MAC_Type value) {
    gravity_mac = value;
  }


  int ComputeActiveParticleList(TreeCellBase<ndim> &, Particle<ndim> *, int *);
  int ComputeActiveCellList(vector<TreeCellBase<ndim> >& );
  int ComputeActiveCellPointers(TreeCellBase<ndim> **celllist);
  int ComputeGatherNeighbourList(const Particle<ndim> *, const FLOAT *,
                                 const FLOAT, const int, int &, int *);
  void ComputeGatherNeighbourList(const TreeCellBase<ndim> &, const Particle<ndim> *,
                                  const FLOAT,NeighbourManagerBase& neibmanager);
  int ComputeNeighbourList(const TreeCellBase<ndim> &, const Particle<ndim> *,
                           const int, int &, int *, Particle<ndim> *);
  void ComputeNeighbourList(const TreeCellBase<ndim> &cell,NeighbourManagerBase& neibmanager);
  void ComputeNeighbourAndGhostList(const TreeCellBase<ndim> &, NeighbourManagerBase&);
  void ComputeGravityInteractionAndGhostList(const TreeCellBase<ndim> &, NeighbourManagerDim<ndim>& neibmanager);
  int ComputeStarGravityInteractionList(const NbodyParticle<ndim> *, const FLOAT, const int,
                                        const int, const int, int &, int &, int &, int *, int *,
                                        MultipoleMoment<ndim> *, Particle<ndim> *);

  virtual bool ComputeSignalVelocityFromDistantInteractions(const TreeCellBase<ndim>& cell,
                                                            int Nactive, Particle<ndim>* active_gen,
                                                            Particle<ndim>* part_gen);
  virtual int FindLeafCell(const FLOAT *);
#ifdef MPI_PARALLEL
  int ComputeImportedCellList(vector<TreeCellBase<ndim> >& );
  int CreatePrunedTreeForMpiNode(const MpiNode<ndim> &, const DomainBox<ndim> &, const FLOAT,
                                 const bool, const int, const int, const int, TreeBase<ndim> *);
  int ComputeDistantGravityInteractionList(const TreeCellBase<ndim>&, const DomainBox<ndim> &,
                                           vector<MultipoleMoment<ndim> >&);
  bool ComputeHydroTreeCellOverlap(const TreeCellBase<ndim> *, const DomainBox<ndim> &);
  FLOAT ComputeWorkInBox(const FLOAT *, const FLOAT *);
  virtual void UpdateLeafCells(TreeBase<ndim>*);
  virtual void CopyLeafCells(vector<char>& buffer, enum TreeBase<ndim>::direction);
  virtual void UpdateBoundingBoxData(MpiNode<ndim>& node) {
	for (int k=0; k<ndim; k++) node.rbox.min[k] = min(node.rbox.min[k],celldata[0].bb.min[k]);
	for (int k=0; k<ndim; k++) node.rbox.max[k] = max(node.rbox.max[k],celldata[0].bb.max[k]);
	for (int k=0; k<ndim; k++) node.hbox.min[k] = min(node.hbox.min[k],celldata[0].hbox.min[k]);
	for (int k=0; k<ndim; k++) node.hbox.max[k] = max(node.hbox.max[k],celldata[0].hbox.max[k]);
  }
#endif


  //-----------------------------------------------------------------------------------------------

  virtual void BuildTree(const int, const int, const int, const int,
                         const FLOAT, Particle<ndim> *) = 0;
  virtual void UpdateAllHmaxValues(Particle<ndim> *, bool stock_leaf) = 0;
  virtual void UpdateHmaxLeaf(TreeCellBase<ndim>&, Particle<ndim> *);
  virtual void StockTree(Particle<ndim> *, bool) = 0 ;
  virtual void UpdateActiveParticleCounters(Particle<ndim> *) = 0;
  virtual double GetMaximumSmoothingLength() const {
    return celldata[0].hmax ;
  }

  virtual void GenerateBoundaryGhostParticles(const FLOAT, const FLOAT, const int,
                                              const DomainBox<ndim>&,
                                              Hydrodynamics<ndim>*) ;

#ifdef MPI_PARALLEL
  virtual void InitialiseCellWorkCounters() {
    assert(Ncell > 0);
    for (int c=0; c<Ncell; c++) celldata[c].worktot = 0 ;
  }
  virtual void UpdateWorkCounters() = 0;
  virtual int GetMaxCellNumber(const int) = 0;
  virtual int GetTreeCellSize() const { return sizeof(TreeCell<ndim>) ;}
  virtual void AddWorkCost(vector<TreeCellBase<ndim> >&, double twork, int& Nactivetot) ;
  virtual void FindBoxOverlapParticles(const Box<ndim>&, vector<int>&, const Particle<ndim>*) ;
  virtual int FindBoxGhostParticles(const FLOAT, const FLOAT, const Box<ndim> &,
                                    vector<int> &export_list) ;
  virtual int PackParticlesAndCellsForMPITransfer(const vector<int>& celllist,
                                                  vector<int>& cell_ids, vector<int>& part_ids,
                                                  vector<char>& send_buffer,
                                                  Particle<ndim>*);
  virtual void PackParticlesAndCellsForMPIReturn(int start_part, int Npart,
                                                 int start_cell, int Ncell,
                                                 vector<char>& send_buffer,
                                                 Particle<ndim>*);

  virtual void UnpackParticlesAndCellsFromMPITransfer(int offset_part, int Npart,
                                                      int offset_cell, int Ncell,
                                                      const vector<char>& recv_buffer,
                                                      Hydrodynamics<ndim>*) ;

  virtual void UnpackParticlesAndCellsForMPIReturn(const vector<int>& part_ids,
                                                   const vector<int>& cell_ids,
                                                   vector<char>& recv_buffer,
                                                   Hydrodynamics<ndim>*) ;

  virtual int GetSizeOfExportedParticleData(int Nparticles) const  ;
  virtual int GetSizeOfExportedCellData(int Ncell) const ;
  virtual int GetSizeOfReturnedParticleData(int Nparticles) const ;
  virtual int GetSizeOfReturnedCellData(int Ncell) const ;
  virtual const void* GetCellDataPointer() const {
    return celldata ;
  }
  virtual void* GetCellDataPointer() {
    return celldata ;
  }
#endif


  // Const variables for tree class
  //-----------------------------------------------------------------------------------------------
  MAC_Type gravity_mac;                ///< Multipole-acceptance criteria for tree
  const multipole_method multipole;    ///< Multipole-order for cell gravity
  const int Nleafmax;                  ///< Max. number of particles per leaf cell
  const FLOAT invthetamaxsqd;          ///< 1 / thetamaxsqd
  const FLOAT kernrange;               ///< Extent of employed kernel
  const FLOAT macerror;                ///< Error tolerance for gravity tree-MAC
  const FLOAT theta;                   ///< Geometric opening angle
  const FLOAT thetamaxsqd;             ///< Geometric opening angle squared



  // Additional variables for tree class
  //-----------------------------------------------------------------------------------------------
  using TreeBase<ndim>::Ntot;
  using TreeBase<ndim>::Ncell;
  using TreeBase<ndim>::gtot;
  using TreeBase<ndim>::gmax;
  using TreeBase<ndim>::lmax;
  using TreeBase<ndim>::ltot;
  using TreeBase<ndim>::ltot_old;
  using TreeBase<ndim>::_domain;
#if defined MPI_PARALLEL
  using TreeBase<ndim>::Nimportedcell;
  using TreeBase<ndim>::Ncelltot;
#endif

  bool allocated_tree;                 ///< Are grid arrays allocated?
  int ifirst;                          ///< i.d. of first particle in tree
  int ilast;                           ///< i.d. of last particle in tree
  int Nthreads;                        ///< No. of OpenMP threads
  FLOAT hmax;                          ///< Store hmax in the tree
  int *g2c;                            ///< i.d. of leaf(grid) cells
  TreeCell<ndim> *celldata;            ///< Main tree cell data array
  Typemask gravmask ;                  ///< Particle types that contribute to gravity

 private:
  bool open_cell_for_gravity(const TreeCell<ndim>& cell,
                             double drsqd, double macfactor, double amag) const
  {
    if (drsqd < cell.cdistsqd)
      return true;

    bool open ;
    switch (gravity_mac) {
    case eigenmac:
      open = drsqd < cell.mac*macfactor;
      break  ;
    case gadget2:
      open = drsqd*drsqd*amag*macerror < cell.rmax*cell.rmax*cell.m;
      break ;
    default:
      open = false ;
      break ;
    }
    return open ;
  }

};

#endif
