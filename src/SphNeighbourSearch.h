//=================================================================================================
//  SphNeighbourSearch.h
//  Header file containing class definitions for all SPH neighbour searching
//  data structures and algorithms.
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


#ifndef _SPH_NEIGHBOUR_SEARCH_H_
#define _SPH_NEIGHBOUR_SEARCH_H_


#include <assert.h>
#include <iostream>
#include <string>
#include <vector>
#include "Precision.h"
#include "Constants.h"
#include "CodeTiming.h"
#include "InlineFuncs.h"
#include "Nbody.h"
#include "SmoothingKernel.h"
#include "Particle.h"
#include "Sph.h"
#include "DomainBox.h"
#include "Ewald.h"
#include "Parameters.h"
#include "KDTree.h"
#include "OctTree.h"
#if defined MPI_PARALLEL
#include "MpiExport.h"
#include "MpiNode.h"
#endif
using namespace std;


// Forward declaration of Sph to break circular dependency
template <int ndim>
class Sph;

// Forward declare MpiNode to break circular dependency
#if defined MPI_PARALLEL
template <int ndim>
class MpiNode;
#endif



//=================================================================================================
//  Class SphNeighbourSearch
/// \brief   SphNeighbourSearch class definition.
/// \details Class for creating the SPH neighbour search data structure, and for computing local
///          neighbour lists and calling SPH functions (e.g. computing h, SPH forces, etc..).
/// \author  D. A. Hubber, G. Rosotti
/// \date    03/04/2013
//=================================================================================================
template <int ndim>
class SphNeighbourSearch
{
#if defined MPI_PARALLEL
protected:
  vector<int> ids_active_particles;
#endif
 public:

  //-----------------------------------------------------------------------------------------------
  SphNeighbourSearch(FLOAT, DomainBox<ndim> *,
                     SmoothingKernel<ndim> *, CodeTiming *);
  ~SphNeighbourSearch();


  //-----------------------------------------------------------------------------------------------
  virtual void BuildTree(bool, int, int, int, int, int,
                         SphParticle<ndim> *, Sph<ndim> *, FLOAT) = 0;
  virtual void BuildGhostTree(bool, int, int, int, int, int,
                              SphParticle<ndim> *, Sph<ndim> *, FLOAT) = 0;
  virtual int GetGatherNeighbourList(FLOAT *, FLOAT, SphParticle<ndim> *, int, int, int *) = 0;
  virtual void UpdateAllSphProperties(int, int, SphParticle<ndim> *,
                                      Sph<ndim> *, Nbody<ndim> *) = 0;
  virtual void UpdateAllSphForces(int, int, SphParticle<ndim> *, Sph<ndim> *, Nbody<ndim> *) = 0;
  virtual void UpdateAllSphHydroForces(int, int, SphParticle<ndim> *,
                                       Sph<ndim> *, Nbody<ndim> *) = 0;
  virtual void UpdateAllSphGravForces(int, int, SphParticle<ndim> *,
                                      Sph<ndim> *, Nbody<ndim> *) = 0;
  virtual void UpdateAllSphPeriodicForces(int, int, SphParticle<ndim> *, Sph<ndim> *,
                                          Nbody<ndim> *, DomainBox<ndim> &, Ewald<ndim> *) = 0;
  virtual void UpdateAllSphPeriodicGravForces(int, int, SphParticle<ndim> *, Sph<ndim> *,
                                              Nbody<ndim> *, DomainBox<ndim> &, Ewald<ndim> *) = 0;
  virtual void UpdateActiveParticleCounters(SphParticle<ndim> *, Sph<ndim> *) = 0;
  virtual void UpdateAllStarGasForces(int, int, SphParticle<ndim> *,
                                      Sph<ndim> *, Nbody<ndim> *) = 0;
  virtual void UpdateAllSphDudt(int, int, SphParticle<ndim> *, Sph<ndim> *) = 0;
  virtual void UpdateAllSphDerivatives(int, int, SphParticle<ndim> *, Sph<ndim> *) = 0;
  virtual void SearchBoundaryGhostParticles(FLOAT, DomainBox<ndim>, Sph<ndim> *) = 0;
#ifdef MPI_PARALLEL
  virtual void UpdateGravityExportList(int, int, int, SphParticle<ndim> *,
                                       Sph<ndim> *, Nbody<ndim> *) = 0;
  virtual void UpdateHydroExportList(int, int, int, SphParticle<ndim> *,
                                     Sph<ndim> *, Nbody<ndim> *) = 0;
  virtual void BuildPrunedTree(int, int) {};
  virtual void BuildGhostPrunedTree(const int, const DomainBox<ndim> &) {};
  virtual void BuildMpiGhostTree(bool, int, int, int, int, int,
                                 SphParticle<ndim> *, Sph<ndim> *, FLOAT) = 0;
  virtual int SearchMpiGhostParticles(const FLOAT, const Box<ndim> &,
                                      Sph<ndim> *, vector<int> &) = 0;
  virtual void FindMpiTransferParticles(Sph<ndim> *, vector<vector<int> >&,
                                        vector<int>&, const vector<int>&, MpiNode<ndim>*) = 0;
  virtual FLOAT FindLoadBalancingDivision(int, FLOAT, FLOAT *, FLOAT *) = 0;
  virtual int GetExportInfo(int Nproc, Sph<ndim>* sph, vector<char >&,
                            MpiNode<ndim>&, int, int) = 0;
  virtual void UnpackExported (vector<char >& arrays, vector<int>& N_received_particles_from_proc,
                               Sph<ndim>* sph) = 0;
  virtual void GetBackExportInfo(vector<char >& received_array,
                                 vector<int>& N_exported_particles_from_proc,
                                 vector<int>&, Sph<ndim>* sph, int rank) = 0;
  virtual void UnpackReturnedExportInfo(vector<char >& received_information,
                                        vector<int>& recv_displs, Sph<ndim>* sph, int rank) = 0;
  virtual void CommunicatePrunedTrees(vector<int>&, int) = 0;
#endif


  //-----------------------------------------------------------------------------------------------
  const FLOAT kernrange;            ///< Kernel extent (in units of h)
  const FLOAT kernrangesqd;         ///< Kernel extent (squared)


  //-----------------------------------------------------------------------------------------------
  CodeTiming *timing;               ///< Pointer to code timing object
  DomainBox<ndim> *box;             ///< Pointer to simulation bounding box
  SmoothingKernel<ndim> *kernp;           ///< Pointer to SPH kernel object

  bool neibcheck;                   ///< Flag to verify neighbour lists
  FLOAT kernfac;                    ///< ..

};



//=================================================================================================
//  Class BruteForceSearch
/// \brief   ..
/// \details Class for computing SPH neighbour lists using brute force only
///          (i.e. direct summation over all particles).
/// \author  D. A. Hubber, G. Rosotti
/// \date    03/04/2013
//=================================================================================================
template <int ndim, template<int> class ParticleType>
class BruteForceSearch: public SphNeighbourSearch<ndim>
{
 public:

  using SphNeighbourSearch<ndim>::neibcheck;
  using SphNeighbourSearch<ndim>::timing;
  using SphNeighbourSearch<ndim>::kernp;
  using SphNeighbourSearch<ndim>::kernfac;
  using SphNeighbourSearch<ndim>::kernrange;
  using SphNeighbourSearch<ndim>::kernrangesqd;


  //-----------------------------------------------------------------------------------------------
  BruteForceSearch(FLOAT, DomainBox<ndim> *, SmoothingKernel<ndim> *, CodeTiming *);
  ~BruteForceSearch();


  //-----------------------------------------------------------------------------------------------
  void BuildTree(bool, int, int, int, int, int, SphParticle<ndim> *, Sph<ndim> *, FLOAT);
  void BuildGhostTree(bool, int, int, int, int, int,
                      SphParticle<ndim> *, Sph<ndim> *, FLOAT) {};
  int GetGatherNeighbourList(FLOAT *, FLOAT, SphParticle<ndim> *, int, int, int *);
  void UpdateAllSphProperties(int, int, SphParticle<ndim> *, Sph<ndim> *, Nbody<ndim> *);
  void UpdateAllSphForces(int, int, SphParticle<ndim> *, Sph<ndim> *, Nbody<ndim> *);
  void UpdateAllSphHydroForces(int, int, SphParticle<ndim> *, Sph<ndim> *, Nbody<ndim> *);
  void UpdateAllSphGravForces(int, int, SphParticle<ndim> *, Sph<ndim> *, Nbody<ndim> *);
  void UpdateAllSphPeriodicForces(int, int, SphParticle<ndim> *, Sph<ndim> *,
                                  Nbody<ndim> *, DomainBox<ndim> &, Ewald<ndim> *);
  void UpdateAllSphPeriodicGravForces(int, int, SphParticle<ndim> *, Sph<ndim> *,
                                      Nbody<ndim> *, DomainBox<ndim> &, Ewald<ndim> *);
  void UpdateActiveParticleCounters(SphParticle<ndim> *, Sph<ndim> *) {};
  void UpdateAllStarGasForces(int, int, SphParticle<ndim> *, Sph<ndim> *, Nbody<ndim> *);
  void UpdateAllSphDudt(int, int, SphParticle<ndim> *, Sph<ndim> *);
  void UpdateAllSphDerivatives(int, int, SphParticle<ndim> *, Sph<ndim> *);
  void SearchBoundaryGhostParticles(FLOAT, DomainBox<ndim>, Sph<ndim> *);

#ifdef MPI_PARALLEL
  using SphNeighbourSearch<ndim>::ids_active_particles;
  void UpdateGravityExportList(int, int, int, SphParticle<ndim> *, Sph<ndim> *, Nbody<ndim> *);
  void UpdateHydroExportList(int, int, int, SphParticle<ndim> *,  Sph<ndim> *, Nbody<ndim> *);
  void BuildMpiGhostTree(bool, int, int, int, int, int,
                         SphParticle<ndim> *, Sph<ndim> *, FLOAT) {};
  int SearchMpiGhostParticles(const FLOAT, const Box<ndim> &, Sph<ndim> *, vector<int> &);
  int SearchHydroExportParticles(const Box<ndim> &, Sph<ndim> *, vector<int> &);
  void FindMpiTransferParticles(Sph<ndim> *, vector<vector<int> >&,
                                vector<int>&, const vector<int>&, MpiNode<ndim>*);
  void FindGhostParticlesToExport(Sph<ndim>* sph, vector<vector<ParticleType<ndim>* > >&,
                                  const vector<int>&, MpiNode<ndim>*);
  FLOAT FindLoadBalancingDivision(int, FLOAT, FLOAT *, FLOAT *) {};
  void FindParticlesToTransfer(Sph<ndim>* sph, std::vector<std::vector<int> >& particles_to_export,
      std::vector<int>& all_particles_to_export, const std::vector<int>& potential_nodes, MpiNode<ndim>* mpinodes);
  virtual int GetExportInfo(int Nproc, Sph<ndim>* sph, vector<char >&, MpiNode<ndim>&, int, int);
  virtual void UnpackExported (vector<char >& arrays, vector<int>& N_received_particles_from_proc,
      Sph<ndim>* sph);
  virtual void GetBackExportInfo(vector<char >& received_array, vector<int>& N_exported_particles_from_proc,
      vector<int>&, Sph<ndim>* sph, int rank);
  virtual void UnpackReturnedExportInfo(vector<char >& received_information, vector<int>& recv_displs,
      Sph<ndim>* sph, int rank);
  virtual void CommunicatePrunedTrees(vector<int>&, int) {};
#endif
};



//=================================================================================================
//  Class GradhSphBruteForce
/// \brief   ..
/// \details Class for computing SPH neighbour lists using brute force only
///          (i.e. direct summation over all particles).
/// \author  D. A. Hubber, G. Rosotti
/// \date    12/05/2014
//=================================================================================================
template <int ndim, template<int> class ParticleType>
class GradhSphBruteForce: public BruteForceSearch<ndim,ParticleType>
{
 public:

  using SphNeighbourSearch<ndim>::neibcheck;
  using SphNeighbourSearch<ndim>::timing;
  using SphNeighbourSearch<ndim>::kernp;
  using SphNeighbourSearch<ndim>::kernfac;
  using SphNeighbourSearch<ndim>::kernrange;
  using SphNeighbourSearch<ndim>::kernrangesqd;


  //-----------------------------------------------------------------------------------------------
  GradhSphBruteForce(FLOAT, DomainBox<ndim> *, SmoothingKernel<ndim> *, CodeTiming *);
  ~GradhSphBruteForce();


  //-----------------------------------------------------------------------------------------------
  void UpdateAllSphProperties(int, int, SphParticle<ndim> *, Sph<ndim> *, Nbody<ndim> *);

};



//=================================================================================================
//  Class SM2012SphBruteForce
/// \brief   ..
/// \details Class for computing SPH neighbour lists using brute force only
///          (i.e. direct summation over all particles).
/// \author  D. A. Hubber, G. Rosotti
/// \date    12/05/2014
//=================================================================================================
template <int ndim, template<int> class ParticleType>
class SM2012SphBruteForce: public BruteForceSearch<ndim,ParticleType>
{
 public:

  using SphNeighbourSearch<ndim>::neibcheck;
  using SphNeighbourSearch<ndim>::timing;
  using SphNeighbourSearch<ndim>::kernp;
  using SphNeighbourSearch<ndim>::kernfac;
  using SphNeighbourSearch<ndim>::kernrange;
  using SphNeighbourSearch<ndim>::kernrangesqd;


  //-----------------------------------------------------------------------------------------------
  SM2012SphBruteForce(FLOAT, DomainBox<ndim> *, SmoothingKernel<ndim> *, CodeTiming *);
  ~SM2012SphBruteForce();


  //-----------------------------------------------------------------------------------------------
  void UpdateAllSphProperties(int, int, SphParticle<ndim> *, Sph<ndim> *, Nbody<ndim> *);

};



//=================================================================================================
//  Class GodunovSphBruteForce
/// \brief   ..
/// \details Class for computing SPH neighbour lists using brute force only
///          (i.e. direct summation over all particles).
/// \author  D. A. Hubber, G. Rosotti
/// \date    12/05/2014
//=================================================================================================
template <int ndim, template<int> class ParticleType>
class GodunovSphBruteForce: public BruteForceSearch<ndim,ParticleType>
{
 public:

  using SphNeighbourSearch<ndim>::neibcheck;
  using SphNeighbourSearch<ndim>::timing;
  using SphNeighbourSearch<ndim>::kernp;
  using SphNeighbourSearch<ndim>::kernfac;
  using SphNeighbourSearch<ndim>::kernrange;
  using SphNeighbourSearch<ndim>::kernrangesqd;


  //-----------------------------------------------------------------------------------------------
  GodunovSphBruteForce(FLOAT, DomainBox<ndim> *, SmoothingKernel<ndim> *, CodeTiming *);
  ~GodunovSphBruteForce();


  //-----------------------------------------------------------------------------------------------
  void UpdateAllSphProperties(int, int, SphParticle<ndim> *, Sph<ndim> *, Nbody<ndim> *);
  void UpdateAllSphDudt(int, int, SphParticle<ndim> *, Sph<ndim> *);
  void UpdateAllSphDerivatives(int, int, SphParticle<ndim> *, Sph<ndim> *);

};



//=================================================================================================
//  Class SphTree
/// \brief   Class containing tree for efficient neighbour searching and gravity calculations.
/// \details ..
/// \author  D. A. Hubber
/// \date    08/01/2014
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
class SphTree: public SphNeighbourSearch<ndim>
{
#if defined MPI_PARALLEL
  vector<vector<int> > ids_sent_particles;
protected:
  using SphNeighbourSearch<ndim>::ids_active_particles;
  vector<int> N_imported_part_per_proc;
#endif
 public:

  using SphNeighbourSearch<ndim>::neibcheck;
  using SphNeighbourSearch<ndim>::box;
  using SphNeighbourSearch<ndim>::timing;
  using SphNeighbourSearch<ndim>::kernp;
  using SphNeighbourSearch<ndim>::kernfac;
  using SphNeighbourSearch<ndim>::kernrange;
  using SphNeighbourSearch<ndim>::kernrangesqd;


  //-----------------------------------------------------------------------------------------------
  SphTree(int, int, FLOAT, FLOAT, FLOAT, string, string,
          DomainBox<ndim> *, SmoothingKernel<ndim> *, CodeTiming *);
  ~SphTree();


  //-----------------------------------------------------------------------------------------------
  void BuildTree(bool, int, int, int, int, int, SphParticle<ndim> *, Sph<ndim> *, FLOAT);
  void BuildGhostTree(bool, int, int, int, int, int, SphParticle<ndim> *, Sph<ndim> *, FLOAT);
  int GetGatherNeighbourList(FLOAT *, FLOAT, SphParticle<ndim> *, int, int, int *);
  void UpdateAllSphProperties(int, int, SphParticle<ndim> *, Sph<ndim> *, Nbody<ndim> *) {};
  void UpdateAllSphForces(int, int, SphParticle<ndim> *, Sph<ndim> *, Nbody<ndim> *) {};
  void UpdateAllSphHydroForces(int, int, SphParticle<ndim> *, Sph<ndim> *, Nbody<ndim> *) {};
  void UpdateAllSphGravForces(int, int, SphParticle<ndim> *, Sph<ndim> *, Nbody<ndim> *) {};
  void UpdateAllSphPeriodicForces(int, int, SphParticle<ndim> *, Sph<ndim> *,
                                  Nbody<ndim> *, DomainBox<ndim> &, Ewald<ndim> *) {};
  void UpdateAllSphPeriodicGravForces(int, int, SphParticle<ndim> *, Sph<ndim> *,
                                      Nbody<ndim> *, DomainBox<ndim> &, Ewald<ndim> *) {};
  void UpdateAllSphDudt(int, int, SphParticle<ndim> *, Sph<ndim> *) {};
  void UpdateAllSphDerivatives(int, int, SphParticle<ndim> *, Sph<ndim> *) {};
  void UpdateActiveParticleCounters(SphParticle<ndim> *, Sph<ndim> *);
  void UpdateAllStarGasForces(int, int, SphParticle<ndim> *, Sph<ndim> *, Nbody<ndim> *) {};
  void SearchBoundaryGhostParticles(FLOAT, DomainBox<ndim>, Sph<ndim> *);
  void ComputeCellMonopoleForces(FLOAT &, FLOAT *, FLOAT *, int, TreeCell<ndim> *);
  void ComputeCellQuadrupoleForces(FLOAT &, FLOAT *, FLOAT *, int, TreeCell<ndim> *);
  void ComputeFastMonopoleForces(int, int, TreeCell<ndim> *,
                                 TreeCell<ndim> &, ParticleType<ndim> *);
#ifdef MPI_PARALLEL
  void UpdateGravityExportList(int, int, int, SphParticle<ndim> *, Sph<ndim> *, Nbody<ndim> *);
  void UpdateHydroExportList(int, int, int, SphParticle<ndim> *, Sph<ndim> *, Nbody<ndim> *);
  void BuildMpiGhostTree(bool, int, int, int, int, int, SphParticle<ndim> *, Sph<ndim> *, FLOAT);
  void BuildPrunedTree(int, int);
  void BuildGhostPrunedTree(const int, const DomainBox<ndim> &);
  int SearchMpiGhostParticles(const FLOAT, const Box<ndim> &, Sph<ndim> *, vector<int> &);
  int SearchHydroExportParticles(const Box<ndim> &, Sph<ndim> *, vector<TreeCell<ndim> *> &);
  void FindMpiTransferParticles(Sph<ndim> *, vector<vector<int> >&,
                                vector<int>&, const vector<int>&, MpiNode<ndim>*);
  FLOAT FindLoadBalancingDivision(int, FLOAT, FLOAT *, FLOAT *);
  virtual int GetExportInfo(int Nproc, Sph<ndim>* sph, vector<char >&, MpiNode<ndim>&, int, int);
  virtual void UnpackExported (vector<char >& arrays, vector<int>& N_received_particles_from_proc,
                               Sph<ndim>* sph);
  virtual void GetBackExportInfo(vector<char >& received_array,
    vector<int>& N_exported_particles_from_proc, vector<int>&, Sph<ndim>* sph, int rank);
  virtual void UnpackReturnedExportInfo(vector<char >& received_information,
                                        vector<int>& recv_displs, Sph<ndim>* sph, int rank);
  virtual void CommunicatePrunedTrees(vector<int>&, int);
#endif
#if defined(VERIFY_ALL)
  void CheckValidNeighbourList(int, int, int, int *, ParticleType<ndim> *, string);
#endif


  // Additional functions for binary tree neighbour search
  //-----------------------------------------------------------------------------------------------
  void AllocateMemory(Sph<ndim> *);
  void DeallocateMemory(void);


  // Additional variables for grid
  //-----------------------------------------------------------------------------------------------
  const int Nleafmax;                              ///< Max. number of particles per leaf cell
  const int Nmpi;                                  ///< No. of MPI processes
  const FLOAT thetamaxsqd;                         ///< Geometric opening angle squared
  const FLOAT invthetamaxsqd;                      ///< 1 / thetamaxsqd
  const FLOAT macerror;                            ///< Error tolerance for gravity tree-MAC
  const string gravity_mac;                        ///< Multipole-acceptance criteria for tree
  const string multipole;                          ///< Multipole-order for cell gravity

  // Class variables
  //-----------------------------------------------------------------------------------------------
  int Ncell;                                       ///< Current no. of grid cells
  int Ncellmax;                                    ///< Max. allowed no. of grid cells
  int Nlistmax;                                    ///< Max. length of neighbour list
  int Ntot;                                        ///< No. of current points in list
  int Ntotold;                                     ///< Prev. no. of particles
  int Ntotmax;                                     ///< Max. no. of points in list
  int Ntotmaxold;                                  ///< Old value of Ntotmax
  FLOAT hmax;                                      ///< Store hmax in the tree
  FLOAT theta;                                     ///< Geometric opening angle
  Tree<ndim,ParticleType,TreeCell> *tree;          ///< Pointer to main (local) tree
  Tree<ndim,ParticleType,TreeCell> *ghosttree;     ///< Pointer to tree containing ghosts
                                                   ///< on local domain

  bool allocated_buffer;                           ///< ..
  int Nthreads;                                    ///< ..
  int *Nneibmaxbuf;                                ///< ..
  int *Ngravcellmaxbuf;                            ///< ..
  int **activelistbuf;                             ///< ..
  int **levelneibbuf;                              ///< ..
  ParticleType<ndim> **neibpartbuf;                ///< Local copy of neighbouring ptcls
  ParticleType<ndim> **activepartbuf;              ///< Local copy of SPH particle
  TreeCell<ndim> **cellbuf;                        ///< ..

#ifdef MPI_PARALLEL
  static const int Nghostprunedmax = 27;           ///< ..
  int Nprunedcellmax;                              ///< ..
  int *Nghostpruned;                               ///< ..
  int *Ncellexport;                                ///< ..
  int *Npartexport;                                ///< ..
  TreeCell<ndim> ***cellexportlist;                ///< List of cells
  Tree<ndim,ParticleType,TreeCell> *mpighosttree;  ///< Pointer to tree containing
                                                   ///< ghosts from other MPI procs.
  Tree<ndim,ParticleType,TreeCell> **prunedtree;   ///< 'Pruned' tree for MPI nodes.
                                                   ///< i.e. only uses top levels
  Tree<ndim,ParticleType,TreeCell> ***ghostprunedtree; ///< Tree of periodic ghost cells created
                                                         ///< from pruned trees
#endif

};



//=================================================================================================
//  Class GradhSphTree
/// \brief   Class containing tree for computing grad-h SPH force loops.
/// \details ...
/// \author  D. A. Hubber
/// \date    08/01/2014
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
class GradhSphTree: public SphTree<ndim,ParticleType,TreeCell>
{
 public:

  using SphTree<ndim,ParticleType,TreeCell>::activelistbuf;
  using SphTree<ndim,ParticleType,TreeCell>::activepartbuf;
  using SphTree<ndim,ParticleType,TreeCell>::allocated_buffer;
  using SphTree<ndim,ParticleType,TreeCell>::box;
  using SphTree<ndim,ParticleType,TreeCell>::cellbuf;
  using SphTree<ndim,ParticleType,TreeCell>::gravity_mac;
  using SphTree<ndim,ParticleType,TreeCell>::kernp;
  using SphTree<ndim,ParticleType,TreeCell>::kernrange;
  using SphTree<ndim,ParticleType,TreeCell>::kernrangesqd;
  using SphTree<ndim,ParticleType,TreeCell>::levelneibbuf;
  using SphTree<ndim,ParticleType,TreeCell>::multipole;
  using SphTree<ndim,ParticleType,TreeCell>::neibcheck;
  using SphTree<ndim,ParticleType,TreeCell>::neibpartbuf;
  using SphTree<ndim,ParticleType,TreeCell>::Ngravcellmaxbuf;
  using SphTree<ndim,ParticleType,TreeCell>::Nleafmax;
  using SphTree<ndim,ParticleType,TreeCell>::Nneibmaxbuf;
  using SphTree<ndim,ParticleType,TreeCell>::Ntot;
  using SphTree<ndim,ParticleType,TreeCell>::Ntotmax;
  using SphTree<ndim,ParticleType,TreeCell>::Ntotmaxold;
  using SphTree<ndim,ParticleType,TreeCell>::Ntotold;
  using SphTree<ndim,ParticleType,TreeCell>::timing;
  using SphTree<ndim,ParticleType,TreeCell>::tree;
  using SphTree<ndim,ParticleType,TreeCell>::ghosttree;
#ifdef MPI_PARALLEL
  using SphTree<ndim,ParticleType,TreeCell>::ghostprunedtree;
  using SphTree<ndim,ParticleType,TreeCell>::mpighosttree;
  using SphTree<ndim,ParticleType,TreeCell>::Nghostpruned;
  using SphTree<ndim,ParticleType,TreeCell>::Nghostprunedmax;
  using SphTree<ndim,ParticleType,TreeCell>::Nmpi;
  using SphTree<ndim,ParticleType,TreeCell>::prunedtree;
#endif


  //-----------------------------------------------------------------------------------------------
  GradhSphTree(int, int, FLOAT, FLOAT, FLOAT, string, string,
               DomainBox<ndim> *, SmoothingKernel<ndim> *, CodeTiming *);
  ~GradhSphTree();


  //-----------------------------------------------------------------------------------------------
  void UpdateAllSphProperties(int, int, SphParticle<ndim> *, Sph<ndim> *, Nbody<ndim> *);
  void UpdateAllSphForces(int, int, SphParticle<ndim> *, Sph<ndim> *, Nbody<ndim> *);
  void UpdateAllSphHydroForces(int, int, SphParticle<ndim> *, Sph<ndim> *, Nbody<ndim> *);
  void UpdateAllSphGravForces(int, int, SphParticle<ndim> *, Sph<ndim> *, Nbody<ndim> *);
  void UpdateAllStarGasForces(int, int, SphParticle<ndim> *, Sph<ndim> *, Nbody<ndim> *);
  void UpdateAllSphPeriodicForces(int, int, SphParticle<ndim> *, Sph<ndim> *,
                                  Nbody<ndim> *, DomainBox<ndim> &, Ewald<ndim> *);

};



//=================================================================================================
//  Class GradhSphKDTree
/// \brief   ...
/// \details ...
/// \author  D. A. Hubber
/// \date    17/09/2014
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
class GradhSphKDTree: public GradhSphTree<ndim,ParticleType,TreeCell>
{
 public:

  using SphTree<ndim,ParticleType,TreeCell>::tree;
  using SphTree<ndim,ParticleType,TreeCell>::ghosttree;
#ifdef MPI_PARALLEL
  using SphTree<ndim,ParticleType,TreeCell>::ghostprunedtree;
  using SphTree<ndim,ParticleType,TreeCell>::mpighosttree;
  using SphTree<ndim,ParticleType,TreeCell>::Nghostpruned;
  using SphTree<ndim,ParticleType,TreeCell>::Nghostprunedmax;
  using SphTree<ndim,ParticleType,TreeCell>::Nmpi;
  using SphTree<ndim,ParticleType,TreeCell>::prunedtree;
#endif


  //-----------------------------------------------------------------------------------------------
  GradhSphKDTree(int, int, FLOAT, FLOAT, FLOAT, string, string,
                 DomainBox<ndim> *, SmoothingKernel<ndim> *, CodeTiming *);

};



//=================================================================================================
//  Class GradhSphOctTree
/// \brief   Class containing octal tree for computing grad-h SPH force loops.
/// \details ...
/// \author  D. A. Hubber
/// \date    17/09/2014
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
class GradhSphOctTree: public GradhSphTree<ndim,ParticleType,TreeCell>
{
 public:

  using SphTree<ndim,ParticleType,TreeCell>::tree;
  using SphTree<ndim,ParticleType,TreeCell>::ghosttree;
#ifdef MPI_PARALLEL
  using SphTree<ndim,ParticleType,TreeCell>::mpighosttree;
  using SphTree<ndim,ParticleType,TreeCell>::Nmpi;
  using SphTree<ndim,ParticleType,TreeCell>::prunedtree;
#endif


  //-----------------------------------------------------------------------------------------------
  GradhSphOctTree(int, int, FLOAT, FLOAT, FLOAT, string, string,
                  DomainBox<ndim> *, SmoothingKernel<ndim> *, CodeTiming *);

};



//=================================================================================================
//  Class SM2012SphTree
/// \brief   Class containing tree for computing Saitoh & Makino (2012) SPH force loops.
/// \details kd-tree data structure used for efficient neighbour searching
///          and computation of gravitational forces for grad-h SPH.
/// \author  D. A. Hubber
/// \date    08/01/2014
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
class SM2012SphTree: public SphTree<ndim,ParticleType,TreeCell>
{
 public:

  using SphTree<ndim,ParticleType,TreeCell>::activelistbuf;
  using SphTree<ndim,ParticleType,TreeCell>::activepartbuf;
  using SphTree<ndim,ParticleType,TreeCell>::allocated_buffer;
  using SphTree<ndim,ParticleType,TreeCell>::box;
  using SphTree<ndim,ParticleType,TreeCell>::cellbuf;
  using SphTree<ndim,ParticleType,TreeCell>::gravity_mac;
  using SphTree<ndim,ParticleType,TreeCell>::kernp;
  using SphTree<ndim,ParticleType,TreeCell>::kernrange;
  using SphTree<ndim,ParticleType,TreeCell>::kernrangesqd;
  using SphTree<ndim,ParticleType,TreeCell>::levelneibbuf;
  using SphTree<ndim,ParticleType,TreeCell>::multipole;
  using SphTree<ndim,ParticleType,TreeCell>::neibcheck;
  using SphTree<ndim,ParticleType,TreeCell>::neibpartbuf;
  using SphTree<ndim,ParticleType,TreeCell>::Ngravcellmaxbuf;
  using SphTree<ndim,ParticleType,TreeCell>::Nleafmax;
  using SphTree<ndim,ParticleType,TreeCell>::Nneibmaxbuf;
  using SphTree<ndim,ParticleType,TreeCell>::Ntot;
  using SphTree<ndim,ParticleType,TreeCell>::Ntotmax;
  using SphTree<ndim,ParticleType,TreeCell>::Ntotmaxold;
  using SphTree<ndim,ParticleType,TreeCell>::Ntotold;
  using SphTree<ndim,ParticleType,TreeCell>::timing;
  using SphTree<ndim,ParticleType,TreeCell>::tree;
  using SphTree<ndim,ParticleType,TreeCell>::ghosttree;
#ifdef MPI_PARALLEL
  using SphTree<ndim,ParticleType,TreeCell>::mpighosttree;
  using SphTree<ndim,ParticleType,TreeCell>::Nmpi;
  using SphTree<ndim,ParticleType,TreeCell>::prunedtree;
#endif


  //-----------------------------------------------------------------------------------------------
  SM2012SphTree(int, int, FLOAT, FLOAT, FLOAT, string, string,
                DomainBox<ndim> *, SmoothingKernel<ndim> *, CodeTiming *);
  ~SM2012SphTree();


  //-----------------------------------------------------------------------------------------------
  void UpdateAllSphProperties(int, int, SphParticle<ndim> *, Sph<ndim> *, Nbody<ndim> *) {};
  void UpdateAllSphForces(int, int, SphParticle<ndim> *, Sph<ndim> *, Nbody<ndim> *) {};
  void UpdateAllSphHydroForces(int, int, SphParticle<ndim> *, Sph<ndim> *, Nbody<ndim> *) {};
  void UpdateAllSphGravForces(int, int, SphParticle<ndim> *, Sph<ndim> *, Nbody<ndim> *) {};
  void UpdateAllStarGasForces(int, int, SphParticle<ndim> *, Sph<ndim> *, Nbody<ndim> *) {};

};



//=================================================================================================
//  Class SM2012SphKDTree
/// \brief   Class containing kd-tree for computing grad-h SPH force loops.
/// \details kd-tree data structure used for efficient neighbour searching
///          and computation of gravitational forces for grad-h SPH.
/// \author  D. A. Hubber
/// \date    17/09/2014
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
class SM2012SphKDTree: public SM2012SphTree<ndim,ParticleType,TreeCell>
{
 public:

  using SphTree<ndim,ParticleType,TreeCell>::tree;
  using SphTree<ndim,ParticleType,TreeCell>::ghosttree;
#ifdef MPI_PARALLEL
  using SphTree<ndim,ParticleType,TreeCell>::mpighosttree;
  using SphTree<ndim,ParticleType,TreeCell>::Nmpi;
  using SphTree<ndim,ParticleType,TreeCell>::prunedtree;
#endif


  //-----------------------------------------------------------------------------------------------
  SM2012SphKDTree(int, int, FLOAT, FLOAT, FLOAT, string, string,
                  DomainBox<ndim> *, SmoothingKernel<ndim> *, CodeTiming *);

};



//=================================================================================================
//  Class SM2012SphOctTree
/// \brief   Class containing kd-tree for computing grad-h SPH force loops.
/// \details kd-tree data structure used for efficient neighbour searching
///          and computation of gravitational forces for grad-h SPH.
/// \author  D. A. Hubber
/// \date    17/09/2014
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
class SM2012SphOctTree: public SM2012SphTree<ndim,ParticleType,TreeCell>
{
 public:

  using SphTree<ndim,ParticleType,TreeCell>::tree;
  using SphTree<ndim,ParticleType,TreeCell>::ghosttree;
#ifdef MPI_PARALLEL
  using SphTree<ndim,ParticleType,TreeCell>::mpighosttree;
  using SphTree<ndim,ParticleType,TreeCell>::Nmpi;
  using SphTree<ndim,ParticleType,TreeCell>::prunedtree;
#endif


  //-----------------------------------------------------------------------------------------------
  SM2012SphOctTree(int, int, FLOAT, FLOAT, FLOAT, string, string,
                   DomainBox<ndim> *, SmoothingKernel<ndim> *, CodeTiming *);

};



//=================================================================================================
//  Class GodunovSphTree
/// \brief   Class containing kd-tree for computing grad-h SPH force loops.
/// \details kd-tree data structure used for efficient neighbour searching
///          and computation of gravitational forces for grad-h SPH.
/// \author  D. A. Hubber
/// \date    08/01/2014
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
class GodunovSphTree: public SphTree<ndim,ParticleType,TreeCell>
{
 public:

  using SphTree<ndim,ParticleType,TreeCell>::activelistbuf;
  using SphTree<ndim,ParticleType,TreeCell>::activepartbuf;
  using SphTree<ndim,ParticleType,TreeCell>::allocated_buffer;
  using SphTree<ndim,ParticleType,TreeCell>::box;
  using SphTree<ndim,ParticleType,TreeCell>::cellbuf;
  using SphTree<ndim,ParticleType,TreeCell>::gravity_mac;
  using SphTree<ndim,ParticleType,TreeCell>::kernp;
  using SphTree<ndim,ParticleType,TreeCell>::kernrange;
  using SphTree<ndim,ParticleType,TreeCell>::kernrangesqd;
  using SphTree<ndim,ParticleType,TreeCell>::levelneibbuf;
  using SphTree<ndim,ParticleType,TreeCell>::multipole;
  using SphTree<ndim,ParticleType,TreeCell>::neibcheck;
  using SphTree<ndim,ParticleType,TreeCell>::neibpartbuf;
  using SphTree<ndim,ParticleType,TreeCell>::Ngravcellmaxbuf;
  using SphTree<ndim,ParticleType,TreeCell>::Nleafmax;
  using SphTree<ndim,ParticleType,TreeCell>::Nneibmaxbuf;
  using SphTree<ndim,ParticleType,TreeCell>::Ntot;
  using SphTree<ndim,ParticleType,TreeCell>::Ntotmax;
  using SphTree<ndim,ParticleType,TreeCell>::Ntotmaxold;
  using SphTree<ndim,ParticleType,TreeCell>::Ntotold;
  using SphTree<ndim,ParticleType,TreeCell>::timing;
  using SphTree<ndim,ParticleType,TreeCell>::tree;
  using SphTree<ndim,ParticleType,TreeCell>::ghosttree;
#ifdef MPI_PARALLEL
  using SphTree<ndim,ParticleType,TreeCell>::mpighosttree;
  using SphTree<ndim,ParticleType,TreeCell>::Nmpi;
  using SphTree<ndim,ParticleType,TreeCell>::prunedtree;
#endif


  //-----------------------------------------------------------------------------------------------
  GodunovSphTree(int, int, FLOAT, FLOAT, FLOAT, string, string,
                 DomainBox<ndim> *, SmoothingKernel<ndim> *, CodeTiming *);
  ~GodunovSphTree();


  //-----------------------------------------------------------------------------------------------
  void UpdateAllSphProperties(int, int, SphParticle<ndim> *, Sph<ndim> *, Nbody<ndim> *) {};
  void UpdateAllSphForces(int, int, SphParticle<ndim> *, Sph<ndim> *, Nbody<ndim> *) {};
  void UpdateAllSphHydroForces(int, int, SphParticle<ndim> *, Sph<ndim> *, Nbody<ndim> *) {};
  void UpdateAllSphGravForces(int, int, SphParticle<ndim> *, Sph<ndim> *, Nbody<ndim> *) {};
  void UpdateAllStarGasForces(int, int, SphParticle<ndim> *, Sph<ndim> *, Nbody<ndim> *) {};
  void UpdateAllSphDudt(int, int, SphParticle<ndim> *, Sph<ndim> *) {};
  void UpdateAllSphDerivatives(int, int, SphParticle<ndim> *, Sph<ndim> *) {};

};



//=================================================================================================
//  Class GodunovSphKDTree
/// \brief   Class containing kd-tree for computing grad-h SPH force loops.
/// \details kd-tree data structure used for efficient neighbour searching
///          and computation of gravitational forces for grad-h SPH.
/// \author  D. A. Hubber
/// \date    17/09/2014
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
class GodunovSphKDTree: public GodunovSphTree<ndim,ParticleType,TreeCell>
{
 public:

  using SphTree<ndim,ParticleType,TreeCell>::tree;
  using SphTree<ndim,ParticleType,TreeCell>::ghosttree;
#ifdef MPI_PARALLEL
  using SphTree<ndim,ParticleType,TreeCell>::mpighosttree;
  using SphTree<ndim,ParticleType,TreeCell>::Nmpi;
  using SphTree<ndim,ParticleType,TreeCell>::prunedtree;
#endif


  //-----------------------------------------------------------------------------------------------
  GodunovSphKDTree(int, int, FLOAT, FLOAT, FLOAT, string, string,
                  DomainBox<ndim> *, SmoothingKernel<ndim> *, CodeTiming *);

};



//=================================================================================================
//  Class GodunovSphOctTree
/// \brief   Class containing kd-tree for computing grad-h SPH force loops.
/// \details kd-tree data structure used for efficient neighbour searching
///          and computation of gravitational forces for grad-h SPH.
/// \author  D. A. Hubber
/// \date    17/09/2014
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
class GodunovSphOctTree: public GodunovSphTree<ndim,ParticleType,TreeCell>
{
 public:

  using SphTree<ndim,ParticleType,TreeCell>::tree;
  using SphTree<ndim,ParticleType,TreeCell>::ghosttree;
#ifdef MPI_PARALLEL
  using SphTree<ndim,ParticleType,TreeCell>::mpighosttree;
  using SphTree<ndim,ParticleType,TreeCell>::Nmpi;
  using SphTree<ndim,ParticleType,TreeCell>::prunedtree;
#endif


  //-----------------------------------------------------------------------------------------------
  GodunovSphOctTree(int, int, FLOAT, FLOAT, FLOAT, string, string,
                   DomainBox<ndim> *, SmoothingKernel<ndim> *, CodeTiming *);

};
#endif
