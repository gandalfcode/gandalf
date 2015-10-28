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
#include "Particle.h"
#include "MeshlessFV.h"
#include "DomainBox.h"
#include "Ewald.h"
#include "Parameters.h"
#include "KDTree.h"
#include "OctTree.h"
#include "Tree.h"
#if defined MPI_PARALLEL
#include "MpiExport.h"
#include "MpiNode.h"
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
  vector<int> ids_active_particles;
#endif
 public:

  //-----------------------------------------------------------------------------------------------
  NeighbourSearch(FLOAT kernrangeaux, DomainBox<ndim> *boxaux,
                  SmoothingKernel<ndim> *kernaux, CodeTiming *timingaux) :
    kernfac(1.0),
    kernrange(kernrangeaux),
    kernrangesqd(kernrangeaux*kernrangeaux),
    timing(timingaux),
    box(boxaux),
    kernp(kernaux) {};
  virtual ~NeighbourSearch() {};


  //-----------------------------------------------------------------------------------------------
  virtual void BuildTree(const bool, const int, const int, const int, const int,
                         const int, const FLOAT, Particle<ndim> *, Hydrodynamics<ndim> *) = 0;
  virtual void BuildGhostTree(const bool, const int, const int, const int, const int,
                              const int, const FLOAT, Particle<ndim> *, Hydrodynamics<ndim> *) = 0;
  virtual int GetGatherNeighbourList(FLOAT *, FLOAT, Particle<ndim> *, int, int, int *) = 0;
  virtual void SearchBoundaryGhostParticles(FLOAT, DomainBox<ndim> &, Hydrodynamics<ndim> *) = 0;
  virtual void UpdateActiveParticleCounters(Particle<ndim> *, Hydrodynamics<ndim> *) = 0;
  virtual void UpdateAllStarGasForces(int, int, Particle<ndim> *,
                                      Hydrodynamics<ndim> *, Nbody<ndim> *) = 0;
#ifdef MPI_PARALLEL
  virtual void BuildPrunedTree(const int, const int, const DomainBox<ndim> &,
                               const MpiNode<ndim> *, Particle<ndim> *) = 0;
  virtual void BuildMpiGhostTree(const bool, const int, const int, const int, const int, const int,
                                 const FLOAT, Particle<ndim> *, Hydrodynamics<ndim> *) = 0;
  virtual void CommunicatePrunedTrees(vector<int>&, int) = 0;
  virtual FLOAT FindLoadBalancingDivision(int, FLOAT, FLOAT *, FLOAT *) = 0;
  virtual void FindMpiTransferParticles(Hydrodynamics<ndim> *, vector<vector<int> >&,
                                        vector<int>&, const vector<int>&, MpiNode<ndim>*) = 0;
  virtual void GetBackExportInfo(vector<char >& received_array,
                                 vector<int>& N_exported_particles_from_proc,
                                 vector<int>&, Hydrodynamics<ndim> *hydro, int rank) = 0;
  virtual int GetExportInfo(int Nproc, Hydrodynamics<ndim> *, vector<char >&,
                            MpiNode<ndim>&, int, int) = 0;
  virtual void InitialiseCellWorkCounters(void) = 0;
  virtual int SearchMpiGhostParticles(const FLOAT, const Box<ndim> &,
                                      Hydrodynamics<ndim> *, vector<int> &) {return 0;};
  virtual void UnpackExported(vector<char >& arrays, vector<int>& N_received_particles_from_proc,
                              Hydrodynamics<ndim> *) = 0;
  virtual void UpdateGravityExportList(int, int, int, Particle<ndim> *, Hydrodynamics<ndim> *,
                                       Nbody<ndim> *, const DomainBox<ndim> &) = 0;
  virtual void UpdateHydroExportList(int, int, int, Particle<ndim> *, Hydrodynamics<ndim> *,
                                     Nbody<ndim> *, const DomainBox<ndim> &) = 0;
  virtual void UnpackReturnedExportInfo(vector<char >& received_information,
                                        vector<int>& recv_displs, Hydrodynamics<ndim>* hydro,
                                        int rank) = 0;
  virtual void FindParticlesToTransfer(Hydrodynamics<ndim> *, vector<vector<int> >& ,
                                       vector<int> &, const vector<int> &, MpiNode<ndim> *) = 0;
#endif


  //-----------------------------------------------------------------------------------------------
  bool neibcheck;                      ///< Flag to verify neighbour lists
  FLOAT kernfac;                       ///< Deprecated variable (to be removed)
  FLOAT kernrange;                     ///< Kernel extent (in units of h)
  FLOAT kernrangesqd;                  ///< Kernel extent (squared)

  CodeTiming *timing;                  ///< Pointer to code timing object
  DomainBox<ndim> *box;                ///< Pointer to simulation bounding box
  SmoothingKernel<ndim> *kernp;        ///< Pointer to SPH kernel object

  TreeBase<ndim>* _tree;               ///< Pointer to main tree
#if defined MPI_PARALLEL
  TreeBase<ndim>** _prunedtree;        ///< Pointer to pruned tree arrays
#endif

};



//=================================================================================================
//  Class BruteForceSearch
/// \brief   Class for computing hydro neighbour lists using brute force only.
/// \details Class for computing hydro neighbour lists using brute force only
///          (i.e. direct summation over all particles).
/// \author  D. A. Hubber, G. Rosotti
/// \date    03/04/2013
//=================================================================================================
template <int ndim, template<int> class ParticleType>
class BruteForceSearch : public virtual NeighbourSearch<ndim>
{
 public:

  using NeighbourSearch<ndim>::neibcheck;
  using NeighbourSearch<ndim>::timing;
  using NeighbourSearch<ndim>::kernp;
  using NeighbourSearch<ndim>::kernfac;
  using NeighbourSearch<ndim>::kernrange;
  using NeighbourSearch<ndim>::kernrangesqd;
#if defined MPI_PARALLEL
  using NeighbourSearch<ndim>::ids_active_particles;
#endif


  //-----------------------------------------------------------------------------------------------
  BruteForceSearch(FLOAT kernrangeaux, DomainBox<ndim> *boxaux,
                   SmoothingKernel<ndim> *kernaux, CodeTiming *timingaux) :
    NeighbourSearch<ndim>(kernrangeaux, boxaux, kernaux, timingaux) {};
  virtual ~BruteForceSearch() {};


  //-----------------------------------------------------------------------------------------------
  virtual void BuildTree(const bool, const int, const int, const int, const int,
                         const int, const FLOAT, Particle<ndim> *, Hydrodynamics<ndim> *);
  virtual void BuildGhostTree(const bool, const int, const int, const int, const int,
                              const int, const FLOAT, Particle<ndim> *, Hydrodynamics<ndim> *) {};
  virtual int GetGatherNeighbourList(FLOAT *, FLOAT, Particle<ndim> *, int, int, int *);
  virtual void SearchBoundaryGhostParticles(FLOAT, DomainBox<ndim> &, Hydrodynamics<ndim> *);
  virtual void UpdateActiveParticleCounters(Particle<ndim> *, Hydrodynamics<ndim> *) {};
  virtual void UpdateAllStarGasForces(int, int, Particle<ndim> *,
                                      Hydrodynamics<ndim> *, Nbody<ndim> *);
#ifdef MPI_PARALLEL
  virtual void BuildPrunedTree(const int, const int, const DomainBox<ndim> &,
                               const MpiNode<ndim> *, Particle<ndim> *) {};
  virtual void BuildMpiGhostTree(const bool, const int, const int, const int, const int, const int,
                                 const FLOAT, Particle<ndim> *, Hydrodynamics<ndim> *) {};
  virtual void CommunicatePrunedTrees(vector<int>&, int) {};
  virtual FLOAT FindLoadBalancingDivision(int, FLOAT, FLOAT *, FLOAT *) {return (FLOAT) 0.0;};
  virtual void FindMpiTransferParticles(Hydrodynamics<ndim> *, vector<vector<int> >&,
                                        vector<int>&, const vector<int>&, MpiNode<ndim>*);
  virtual void GetBackExportInfo(vector<char >& received_array,
                                 vector<int>& N_exported_particles_from_proc,
                                 vector<int>&, Hydrodynamics<ndim> *hydro, int rank);
  virtual int GetExportInfo(int Nproc, Hydrodynamics<ndim> *, vector<char >&,
                            MpiNode<ndim>&, int, int);
  virtual void InitialiseCellWorkCounters(void) {};
  virtual int SearchMpiGhostParticles(const FLOAT, const Box<ndim> &,
                                      Hydrodynamics<ndim> *, vector<int> &);
  virtual void UnpackExported(vector<char>& arrays, vector<int>& N_received_particles_from_proc,
                              Hydrodynamics<ndim> *);
  virtual void UpdateGravityExportList(int, int, int, Particle<ndim> *, Hydrodynamics<ndim> *,
                                       Nbody<ndim> *, const DomainBox<ndim> &) {};
  virtual void UpdateHydroExportList(int, int, int, Particle<ndim> *, Hydrodynamics<ndim> *,
                                     Nbody<ndim> *, const DomainBox<ndim> &) {};
  virtual void UnpackReturnedExportInfo(vector<char> &, vector<int> &, Hydrodynamics<ndim> *, int);
  virtual void FindParticlesToTransfer(Hydrodynamics<ndim> *, vector<vector<int> >& ,
                                       vector<int> &, const vector<int> &, MpiNode<ndim> *);
#endif

};



//=================================================================================================
//  Class HydroTree
/// \brief   Class containing tree for efficient neighbour searching and gravity calculations.
/// \details Class containing tree for efficient neighbour searching and gravity calculations.
/// \author  D. A. Hubber
/// \date    08/01/2014
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
class HydroTree : public virtual NeighbourSearch<ndim>
{
#if defined MPI_PARALLEL
  vector<vector<int> > ids_sent_particles;
protected:
  using NeighbourSearch<ndim>::ids_active_particles;
  vector<int> N_imported_part_per_proc;
#endif
 public:

  using NeighbourSearch<ndim>::neibcheck;
  using NeighbourSearch<ndim>::box;
  using NeighbourSearch<ndim>::timing;
  using NeighbourSearch<ndim>::kernp;
  using NeighbourSearch<ndim>::kernfac;
  using NeighbourSearch<ndim>::kernrange;
  using NeighbourSearch<ndim>::kernrangesqd;


  //-----------------------------------------------------------------------------------------------
  HydroTree(int, int, int, int, FLOAT, FLOAT, FLOAT, string, string,
            DomainBox<ndim> *, SmoothingKernel<ndim> *, CodeTiming *);
  virtual ~HydroTree();


  //-----------------------------------------------------------------------------------------------
  virtual void BuildTree(const bool, const int, const int, const int, const int,
                         const int, const FLOAT, Particle<ndim> *, Hydrodynamics<ndim> *);
  virtual void BuildGhostTree(const bool, const int, const int, const int, const int,
                              const int, const FLOAT, Particle<ndim> *, Hydrodynamics<ndim> *);
  virtual int GetGatherNeighbourList(FLOAT *, FLOAT, Particle<ndim> *, int, int, int *);
  virtual void SearchBoundaryGhostParticles(FLOAT, DomainBox<ndim> &, Hydrodynamics<ndim> *);
  virtual void UpdateActiveParticleCounters(Particle<ndim> *, Hydrodynamics<ndim> *);
  virtual void UpdateAllStarGasForces(int, int, Particle<ndim> *,
                                      Hydrodynamics<ndim> *, Nbody<ndim> *);
#ifdef MPI_PARALLEL
  virtual void BuildPrunedTree(const int, const int, const DomainBox<ndim> &,
                               const MpiNode<ndim> *, Particle<ndim> *);
  virtual void BuildMpiGhostTree(const bool, const int, const int, const int, const int, const int,
                                 const FLOAT, Particle<ndim> *, Hydrodynamics<ndim> *);
  virtual void CommunicatePrunedTrees(vector<int> &, int);
  virtual FLOAT FindLoadBalancingDivision(int, FLOAT, FLOAT *, FLOAT *);
  virtual void FindMpiTransferParticles(Hydrodynamics<ndim> *, vector<vector<int> >&,
                                        vector<int>&, const vector<int>&, MpiNode<ndim>*);
  virtual void GetBackExportInfo(vector<char > &, vector<int> &,
                                 vector<int> &, Hydrodynamics<ndim> *, int);
  virtual int GetExportInfo(int, Hydrodynamics<ndim> *, vector<char >&, MpiNode<ndim>&, int, int);
  virtual void InitialiseCellWorkCounters(void);
  virtual int SearchMpiGhostParticles(const FLOAT, const Box<ndim> &,
                                      Hydrodynamics<ndim> *, vector<int> &);
  virtual void UnpackExported(vector<char> &, vector<int> &, Hydrodynamics<ndim> *);
  virtual void UpdateGravityExportList(int, int, int, Particle<ndim> *, Hydrodynamics<ndim> *,
                                       Nbody<ndim> *, const DomainBox<ndim> &);
  virtual void UpdateHydroExportList(int, int, int, Particle<ndim> *, Hydrodynamics<ndim> *,
                                     Nbody<ndim> *, const DomainBox<ndim> &);
  virtual void UnpackReturnedExportInfo(vector<char > &, vector<int> &,
                                        Hydrodynamics<ndim> *, int);
  virtual void FindParticlesToTransfer(Hydrodynamics<ndim> *, vector<vector<int> >& ,
                                       vector<int> &, const vector<int> &, MpiNode<ndim> *);
#endif
#if defined(VERIFY_ALL)
  void CheckValidNeighbourList(int, int, int, int *, ParticleType<ndim> *, string);
#endif


  // Additional functions for binary tree neighbour search
  //-----------------------------------------------------------------------------------------------
  void AllocateMemory(const int);
  void DeallocateMemory(void);
  void ComputeCellMonopoleForces(FLOAT &, FLOAT *, FLOAT *, int, TreeCell<ndim> *);
  void ComputeCellQuadrupoleForces(FLOAT &, FLOAT *, FLOAT *, int, TreeCell<ndim> *);
  void ComputeFastMonopoleForces(int, int, TreeCell<ndim> *,
                                 TreeCell<ndim> &, ParticleType<ndim> *);


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
  const string multipole;                          ///< Multipole-order for cell gravity

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
  ParticleType<ndim> **neibpartbuf;                ///< Local copy of neighbouring ptcls
  ParticleType<ndim> **activepartbuf;              ///< Local copy of SPH particle
  TreeCell<ndim> **cellbuf;                        ///< Buffers of tree-cell copies

  Tree<ndim,ParticleType,TreeCell> *tree;          ///< Pointer to main (local) tree
  Tree<ndim,ParticleType,TreeCell> *ghosttree;     ///< Pointer to tree containing ghosts
                                                   ///< on local domain

#ifdef MPI_PARALLEL
  int Nprunedcellmax;                              ///< Max. number of cells in pruned tree
  int *Ncellexport;                                ///< No. of cells to be exported (per MPI node)
  int *Npartexport;                                ///< No. of ptcls to be exported (per MPI node)
  TreeCell<ndim> ***cellexportlist;                ///< List of cells
  Tree<ndim,ParticleType,TreeCell> *mpighosttree;  ///< Pointer to tree containing
                                                   ///< ghosts from other MPI procs.
  Tree<ndim,ParticleType,TreeCell> **prunedtree;   ///< 'Pruned' tree for MPI nodes.
  Tree<ndim,ParticleType,TreeCell> **sendprunedtree;  ///< 'Pruned' tree for MPI nodes.
#endif

};
#endif
