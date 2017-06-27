//=================================================================================================
//  NeighbourSearch.h
//  Header file containing virtual class definitions for all hydro neighbour searching algorithms.
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


#ifndef _NEIGHBOUR_SEARCH_H_
#define _NEIGHBOUR_SEARCH_H_


#include <assert.h>
#include <iostream>
#include <string>
#include <vector>
#include "Precision.h"
#include "Constants.h"
#include "CodeTiming.h"
#include "Hydrodynamics.h"
#include "InlineFuncs.h"
#include "Nbody.h"
#include "SmoothingKernel.h"
#include "Multipole.h"
#include "Particle.h"
#include "MeshlessFV.h"
#include "DomainBox.h"
#include "Ewald.h"
#include "Parameters.h"
#include "KDTree.h"
#include "OctTree.h"
#include "BruteForceTree.h"
#include "Tree.h"
#if defined MPI_PARALLEL
#include "MpiExport.h"
#include "MpiNode.h"
#include "CommunicationHandler.h"
#endif
using namespace std;



//=================================================================================================
//  Class NeighbourSearch
/// \brief   NeighbourSearch class definition.
/// \details Class for creating the hydro neighbour search data structure, and for computing local
///          neighbour lists and calling hydro functions (e.g. computing h, forces/fluxes, etc..).
/// \author  D. A. Hubber, G. Rosotti
/// \date    20/04/2015
//=================================================================================================
template <int ndim>
class NeighbourSearch
{
#if defined MPI_PARALLEL
protected:

#endif
 public:

  //-----------------------------------------------------------------------------------------------
  virtual ~NeighbourSearch() {};


  //-----------------------------------------------------------------------------------------------
  virtual void BuildTree(const bool, const int, const int,
                         const int, const FLOAT, Hydrodynamics<ndim> *) = 0;
  virtual void BuildGhostTree(const bool, const int, const int,
                              const int, const FLOAT, Hydrodynamics<ndim> *) = 0;
  virtual int GetGatherNeighbourList(FLOAT *, FLOAT, Particle<ndim> *, int, int, int *) = 0;
  virtual void SearchBoundaryGhostParticles(FLOAT, const DomainBox<ndim> &, Hydrodynamics<ndim> *) = 0;
  virtual void UpdateActiveParticleCounters(Hydrodynamics<ndim> *) = 0;
  virtual void UpdateAllStarGasForces(Hydrodynamics<ndim> *, Nbody<ndim> *,
                                      DomainBox<ndim> &, Ewald<ndim> *) = 0;
  virtual void UpdateAllProperties(Hydrodynamics<ndim> *, Nbody<ndim> *, DomainBox<ndim> &) = 0;
  virtual double GetMaximumSmoothingLength() const = 0;
  virtual TreeBase<ndim>* GetTree() const = 0;
  virtual TreeBase<ndim>* GetGhostTree() const = 0;
#ifdef MPI_PARALLEL
  virtual TreeBase<ndim>* GetMPIGhostTree() const = 0;
#endif
  virtual void SetTimingObject(CodeTiming*) = 0 ;
  virtual void ToggleNeighbourCheck(bool do_check) = 0 ;
  virtual void UpdateTimestepsLimitsFromDistantParticles(Hydrodynamics<ndim>*,const bool) = 0 ;

  virtual MAC_Type GetOpeningCriterion() const = 0;
  virtual void SetOpeningCriterion(MAC_Type) = 0;

#ifdef MPI_PARALLEL
  virtual TreeBase<ndim>** GetPrunedTrees() const = 0;
  virtual TreeBase<ndim>*  GetPrunedTree(int i) const = 0;

  virtual void BuildPrunedTree(const int, const DomainBox<ndim> &,
                               const MpiNode<ndim> *, Hydrodynamics<ndim> *) = 0;
  virtual void StockPrunedTree(const int rank,Hydrodynamics<ndim>* hydro) = 0;
  virtual void BuildMpiGhostTree(const bool, const int, const int, const int,
                                 const FLOAT, Hydrodynamics<ndim> *) = 0;
  virtual FLOAT FindLoadBalancingDivision(int, FLOAT, FLOAT *, FLOAT *) = 0;
  virtual void FindMpiTransferParticles(Hydrodynamics<ndim> *, vector<vector<int> >&,
                                        vector<int>&, const vector<int>&, MpiNode<ndim>*) = 0;
  virtual void GetBackExportInfo(vector<char >& received_array, Hydrodynamics<ndim> *hydro,
		  	  	  	  	  	  	  const int rank, const int iproc) = 0;
  virtual vector<char> ExportSize (const int iproc, Hydrodynamics<ndim>* hydro) const =0;
  virtual int ExportInfoSize(const int i) const =0;
  virtual int GetExportInfo(int Nproc, Hydrodynamics<ndim> *, vector<char >&,
                            MpiNode<ndim>&, int, int) = 0;
  virtual void InitialiseCellWorkCounters(void) = 0;
  virtual int SearchMpiGhostParticles(const FLOAT, const Box<ndim> &,
                                      Hydrodynamics<ndim> *, vector<int> &) {return 0;};
  virtual void UnpackExported(vector<char >& arrays, Hydrodynamics<ndim> *, const int,vector< vector<char> >&,
                              const int, const bool) = 0;
  virtual void UpdateGravityExportList(int, Hydrodynamics<ndim> *,
                                       Nbody<ndim> *, const DomainBox<ndim> &,  Ewald<ndim> *) = 0;
  virtual void UpdateHydroExportList(int, Hydrodynamics<ndim> *,
                                     Nbody<ndim> *, const DomainBox<ndim> &) = 0;
  virtual void UnpackReturnedExportInfo(vector<char >& received_information,
                                        Hydrodynamics<ndim>* hydro,
                                        const int rank, const int iproc) = 0;
  virtual void FindParticlesToTransfer(Hydrodynamics<ndim> *, vector<vector<int> >& ,
                                       vector<int> &, const vector<int> &, MpiNode<ndim> *) = 0;
  virtual void ResetCountersExportInfo(Hydrodynamics<ndim>* hydro) = 0;
  virtual void UpdateBoundingBoxData(MpiNode<ndim>& node, const bool) = 0;
#endif

};



//=================================================================================================
//  Class HydroTree
/// \brief   Class containing tree for efficient neighbour searching and gravity calculations.
/// \details Class containing tree for efficient neighbour searching and gravity calculations.
/// \author  D. A. Hubber
/// \date    08/01/2014
//=================================================================================================
template <int ndim, template<int> class ParticleType>
class HydroTree : public virtual NeighbourSearch<ndim>
{
	void ReallocateMemory(void);
#if defined MPI_PARALLEL
  vector<vector<int> > ids_sent_particles;
  vector<vector<int> > ids_sent_cells;
  vector<int> N_imported_part_per_proc;
  vector<int> N_imported_cells_per_proc;
  int rank;
#endif
protected:
 public:


  //-----------------------------------------------------------------------------------------------
  HydroTree(string, int, int, int, int, FLOAT, FLOAT, FLOAT, string, multipole_method,
            DomainBox<ndim> *, SmoothingKernel<ndim> *, CodeTiming *, ParticleTypeRegister&);
  virtual ~HydroTree();


  //-----------------------------------------------------------------------------------------------
  virtual void BuildTree(const bool, const int, const int,
                         const int, const FLOAT, Hydrodynamics<ndim> *);
  virtual void BuildGhostTree(const bool, const int, const int,
                              const int, const FLOAT, Hydrodynamics<ndim> *);
  virtual int GetGatherNeighbourList(FLOAT *, FLOAT, Particle<ndim> *, int, int, int *);
  virtual void SearchBoundaryGhostParticles(FLOAT, const DomainBox<ndim> &, Hydrodynamics<ndim> *);
  virtual void UpdateActiveParticleCounters(Hydrodynamics<ndim> *);
  virtual void UpdateAllStarGasForces(Hydrodynamics<ndim> *, Nbody<ndim> *,
                                      DomainBox<ndim> &, Ewald<ndim> *);
  virtual double GetMaximumSmoothingLength() const;
  virtual TreeBase<ndim>* GetTree() const { return tree; }
  virtual TreeBase<ndim>* GetGhostTree() const { return ghosttree; }
#ifdef MPI_PARALLEL
  virtual TreeBase<ndim>* GetMPIGhostTree() const { return mpighosttree; }
#endif
  virtual void SetTimingObject(CodeTiming* timer) { timing = timer ; }
  virtual void ToggleNeighbourCheck(bool do_check) { neibcheck = do_check; }

  virtual MAC_Type GetOpeningCriterion() const ;
  virtual void SetOpeningCriterion(const MAC_Type) ;


  virtual void UpdateTimestepsLimitsFromDistantParticles(Hydrodynamics<ndim>*,const bool);

#ifdef MPI_PARALLEL
  virtual TreeBase<ndim>** GetPrunedTrees() const { return prunedtree; }
  virtual TreeBase<ndim>*  GetPrunedTree(int i) const { return prunedtree[i]; }

  virtual void BuildPrunedTree(const int, const DomainBox<ndim> &,
                               const MpiNode<ndim> *, Hydrodynamics<ndim> *);
  virtual void StockPrunedTree(const int rank,Hydrodynamics<ndim>* hydro);
  virtual void BuildMpiGhostTree(const bool, const int, const int, const int,
                                 const FLOAT, Hydrodynamics<ndim> *);
  virtual FLOAT FindLoadBalancingDivision(int, FLOAT, FLOAT *, FLOAT *);
  virtual void FindMpiTransferParticles(Hydrodynamics<ndim> *, vector<vector<int> >&,
                                        vector<int>&, const vector<int>&, MpiNode<ndim>*);
  virtual void GetBackExportInfo(vector<char > &,
                                 Hydrodynamics<ndim> *, const int, const int);
  virtual vector<char> ExportSize (const int iproc, Hydrodynamics<ndim>* hydro) const{
    int cactive = cellexportlist[iproc].size();

    int Nactive = Npartexport[iproc];

    const int size_particles  = tree->GetSizeOfExportedParticleData(Nactive);
    const int size_cells      = tree->GetSizeOfExportedCellData(cactive);

    int size = size_particles + size_cells ;

    vector<char> result(3*sizeof(int));
    copy(&result[0],&size);
    copy(&result[sizeof(int)],&Nactive);
    copy(&result[2*sizeof(int)],&cactive);

    return result;
  };
  virtual int ExportInfoSize(const int i) const {
    const int size_particles = tree->GetSizeOfReturnedParticleData(ids_sent_particles[i].size());
    const int size_cells     = tree->GetSizeOfReturnedCellData(ids_sent_cells[i].size()) ;
    return size_particles+size_cells;
  };
  virtual int GetExportInfo(int, Hydrodynamics<ndim> *, vector<char >&, MpiNode<ndim>&, int, int);
  virtual int SearchMpiGhostParticles(const FLOAT, const Box<ndim> &,
                                      Hydrodynamics<ndim> *, vector<int> &);
  virtual void UnpackExported(vector<char> &, Hydrodynamics<ndim> *,
      const int, vector< vector<char> >&, const int, const bool);
  virtual void UpdateGravityExportList(int, Hydrodynamics<ndim> *,
                                       Nbody<ndim> *, const DomainBox<ndim> &,  Ewald<ndim> *);
  virtual void UpdateHydroExportList(int, Hydrodynamics<ndim> *,
                                     Nbody<ndim> *, const DomainBox<ndim> &);
  virtual void UnpackReturnedExportInfo(vector<char > &, Hydrodynamics<ndim> *, const int, const int);
  virtual void FindParticlesToTransfer(Hydrodynamics<ndim> *, vector<vector<int> >& ,
                                       vector<int> &, const vector<int> &, MpiNode<ndim> *);
  virtual void ResetCountersExportInfo (Hydrodynamics<ndim>* hydro) {
	  tree->Ntot -= hydro->NImportedParticles;
	  hydro->Ntot -= hydro->NImportedParticles;
	  assert(hydro->Ntot == hydro->Nhydro + hydro->Nghost);
	  hydro->NImportedParticles=0;
	  tree->Nimportedcell = 0;
  };
  virtual void InitialiseCellWorkCounters() {
CodeTiming::BlockTimer timer = timing->StartNewTimer("INITIALISE_CELL_WORK");
    tree->InitialiseCellWorkCounters() ;
  }
  virtual void UpdateBoundingBoxData(MpiNode<ndim>& node, const bool UseGhosts) {
     // Initialise bounding box values
     for (int k=0; k<ndim; k++) node.rbox.min[k] = big_number;
     for (int k=0; k<ndim; k++) node.rbox.max[k] = -big_number;
     for (int k=0; k<ndim; k++) node.hbox.min[k] = big_number;
     for (int k=0; k<ndim; k++) node.hbox.max[k] = -big_number;
     tree->UpdateBoundingBoxData(node);
     if (UseGhosts)
	ghosttree->UpdateBoundingBoxData(node);
  }
#endif
#if defined(VERIFY_ALL)
  void CheckValidNeighbourList(int, int, int, int *, ParticleType<ndim> *, string);
#endif


  // Additional functions for binary tree neighbour search
  //-----------------------------------------------------------------------------------------------
  void AllocateMemory(const int);
  void DeallocateMemory(void);



  // Const variables
  //-----------------------------------------------------------------------------------------------

  const int pruning_level_min;                     ///< Minimum pruned tree level
  const int pruning_level_max;                     ///< Maximum pruned tree level
  const int Nleafmax;                              ///< Max. number of particles per leaf cell
  const int Nmpi;                                  ///< No. of MPI processes
  const FLOAT thetamaxsqd;                         ///< Geometric opening angle squared
  const FLOAT invthetamaxsqd;                      ///< 1 / thetamaxsqd
  const FLOAT macerror;                            ///< Error tolerance for gravity tree-MAC
  const string gravity_mac;                        ///< Multipole-acceptance criteria for tree
  const multipole_method multipole;                ///< Multipole-order for cell gravity


  //-----------------------------------------------------------------------------------------------
  bool neibcheck;                      ///< Flag to verify neighbour lists
  FLOAT kernrange;                     ///< Kernel extent (in units of h)
  FLOAT kernrangesqd;                  ///< Kernel extent (squared)


  // Class variables
  //-----------------------------------------------------------------------------------------------
  //int Nlistmax;                                    ///< Max. length of neighbour list
  int Ntot;                                        ///< No. of current points in list
  int Ntotold;                                     ///< Prev. no. of particles
  int Ntotmax;                                     ///< Max. no. of points in list
  int Ntotmaxold;                                  ///< Old value of Ntotmax
  //FLOAT hmax;                                      ///< Store hmax in the tree
  FLOAT theta;                                     ///< Geometric opening angle

  bool allocated_buffer;                           ///< Is buffer memory allocated?
  int Nthreads;                                    ///< No. of OpenMP threads
  int *Nneibmaxbuf;                                ///< Size of neighbour buffers (for each thread)
  int *Ngravcellmaxbuf;                            ///< Size of tree-cell buffers (for each thread)
  int **activelistbuf;                             ///< Arrays of active particle ids
  int **levelneibbuf;                              ///< Arrays of neighbour timestep levels
  int **neiblistbuf;                               ///< Array of neighbour ids
  int **ptypebuf;                                  ///< Array of neighbour ptypes
  //ParticleType<ndim> **neibpartbuf;                ///< Local copy of neighbouring ptcls
  ParticleType<ndim> **activepartbuf;              ///< Local copy of SPH particle

  TreeBase<ndim> *tree;                            ///< Pointer to main (local) tree
  TreeBase<ndim> *ghosttree;                       ///< Pointer to tree containing ghosts
                                                   ///< on local domain

  CodeTiming *timing;                              ///< Pointer to code timing object
  DomainBox<ndim> *box;                            ///< Pointer to simulation bounding box
  SmoothingKernel<ndim> *kernp;                    ///< Pointer to SPH kernel object

#ifdef MPI_PARALLEL
  int Nprunedcellmax;                              ///< Max. number of cells in pruned tree
  int *Npartexport;                                ///< No. of ptcls to be exported (per MPI node)

  vector<vector<int> > cellexportlist;             ///< List of cell ids
  TreeBase<ndim> *mpighosttree;                    ///< Pointer to tree containing
                                                   ///< ghosts from other MPI procs.
  TreeBase<ndim> **prunedtree;                     ///< 'Pruned' tree for MPI nodes.
  TreeBase<ndim> **sendprunedtree;                 ///< 'Pruned' tree for MPI nodes.
#endif

 private:
  TreeBase<ndim>* CreateTree(string tree_type,
                             int Nleafmax, FLOAT thetamaxsqd,
                             FLOAT kernrange, FLOAT macerror,
                             string gravity_mac, multipole_method multipole,
                             const DomainBox<ndim>& domain,
                             const ParticleTypeRegister& reg,
							 const bool IAmPruned) ;
};

#endif
