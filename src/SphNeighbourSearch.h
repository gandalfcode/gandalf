//=============================================================================
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
//=============================================================================


#ifndef _SPH_NEIGHBOUR_SEARCH_H_
#define _SPH_NEIGHBOUR_SEARCH_H_

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
#include "Sinks.h"
#include "DomainBox.h"
#include "Parameters.h"
#include "KDTree.h"
using namespace std;



//=============================================================================
//  Class SphNeighbourSearch
/// \brief   SphNeighbourSearch class definition.  
/// \details Contains routines for creating the SPH neighbour search data
///          structure, and for computing local neighbour lists and calling 
///          SPH functions (e.g. computing h, SPH forces, etc..).
/// \author  D. A. Hubber, G. Rosotti
/// \date    03/04/2013
//=============================================================================
template <int ndim>
class SphNeighbourSearch
{
 public:

  SphNeighbourSearch();
  ~SphNeighbourSearch();

  virtual void BuildTree(bool, int, int, int, int, int, 
                         SphParticle<ndim> *, Sph<ndim> *, FLOAT) = 0;
  virtual void UpdateAllSphProperties(int, int, SphParticle<ndim> *, 
                                     Sph<ndim> *, Nbody<ndim> *) = 0;
  virtual void UpdateAllSphForces(int, int, SphParticle<ndim> *, 
                                  Sph<ndim> *, Nbody<ndim> *) = 0;
  virtual void UpdateAllSphHydroForces(int, int, SphParticle<ndim> *, 
                                       Sph<ndim> *, Nbody<ndim> *) = 0;
  virtual void UpdateAllSphGravForces(int, int, SphParticle<ndim> *, 
                                      Sph<ndim> *, Nbody<ndim> *) = 0;
  virtual void UpdateAllSphDudt(int, int, SphParticle<ndim> *, 
                                Sph<ndim> *) = 0;
  virtual void UpdateAllSphDerivatives(int, int, SphParticle<ndim> *, 
                                       Sph<ndim> *) = 0;
  virtual void UpdateActiveParticleCounters(SphParticle<ndim> *, 
                                            Sph<ndim> *) = 0;
  virtual void UpdateAllStarGasForces(int, int, SphParticle<ndim> *, 
                                      Sph<ndim> *, Nbody<ndim> *) = 0;
#if defined(VERIFY_ALL)
  void CheckValidNeighbourList(int, int, SphParticle<ndim> *, Sph<ndim> *, 
                               int, int, int *, string);
#endif

  bool neibcheck;                   ///< Flag to verify neighbour lists
  CodeTiming *timing;               ///< Pointer to code timing object
  DomainBox<ndim> *box;             ///< Pointer to simulation bounding box
  SphKernel<ndim> *kernp;           ///< Pointer to SPH kernel object
  FLOAT kernfac;
  FLOAT kernrange;

};

#if defined MPI_PARALLEL
//Forward declare MpiNode to break circular dependency
template <int ndim>
class MpiNode;
#endif



//=============================================================================
//  Class BruteForceSearch
/// \brief   ..
/// \details Class for computing SPH neighbour lists using brute force only 
///          (i.e. direct summation over all particles).
/// \author  D. A. Hubber, G. Rosotti
/// \date    03/04/2013
§//=============================================================================
template <int ndim, template<int> class ParticleType>
class BruteForceSearch: public SphNeighbourSearch<ndim>
{
  using SphNeighbourSearch<ndim>::neibcheck;
  using SphNeighbourSearch<ndim>::timing;
  using SphNeighbourSearch<ndim>::kernp;
  using SphNeighbourSearch<ndim>::kernfac;
  using SphNeighbourSearch<ndim>::kernrange;
  
 public:

  BruteForceSearch();
  ~BruteForceSearch();

  void BuildTree(bool, int, int, int, int, int, 
                 SphParticle<ndim> *, Sph<ndim> *, FLOAT);
  void UpdateAllSphProperties(int, int, SphParticle<ndim> *, 
                              Sph<ndim> *, Nbody<ndim> *);
  void UpdateAllSphForces(int, int, SphParticle<ndim> *, 
                          Sph<ndim> *, Nbody<ndim> *);
  void UpdateAllSphHydroForces(int, int, SphParticle<ndim> *, 
                               Sph<ndim> *, Nbody<ndim> *);
  void UpdateAllSphGravForces(int, int, SphParticle<ndim> *, 
                              Sph<ndim> *, Nbody<ndim> *);
  void UpdateAllSphDudt(int, int, SphParticle<ndim> *, Sph<ndim> *);
  void UpdateAllSphDerivatives(int, int, SphParticle<ndim> *, Sph<ndim> *);
  void UpdateActiveParticleCounters(SphParticle<ndim> *, Sph<ndim> *);
  void UpdateAllStarGasForces(int, int, SphParticle<ndim> *, 
                              Sph<ndim> *, Nbody<ndim> *);
#if defined MPI_PARALLEL
  void FindGhostParticlesToExport(Sph<ndim>* sph, std::vector<std::vector<SphParticle<ndim>* > >&,
      const std::vector<int>&, MpiNode<ndim>*);
  void FindParticlesToTransfer(Sph<ndim>* sph, std::vector<std::vector<int> >& particles_to_export,
      std::vector<int>& all_particles_to_export, const std::vector<int>& potential_nodes, MpiNode<ndim>* mpinodes);
#endif
};



//=============================================================================
//  Class SphTree
/// \brief   Class containing binary tree
/// \details Binary tree data structure used for efficient neighbour searching 
///          and computation of gravitational forces
/// \author  D. A. Hubber
/// \date    08/01/2014
//=============================================================================
template <int ndim, template<int> class ParticleType>
class SphTree: public SphNeighbourSearch<ndim>
{
 public:

  using SphNeighbourSearch<ndim>::neibcheck;
  using SphNeighbourSearch<ndim>::box;
  using SphNeighbourSearch<ndim>::timing;
  using SphNeighbourSearch<ndim>::kernp;
  using SphNeighbourSearch<ndim>::kernfac;
  using SphNeighbourSearch<ndim>::kernrange;


  SphTree(int, FLOAT, FLOAT, FLOAT, string, string);
  ~SphTree();


  //---------------------------------------------------------------------------
  void BuildTree(bool, int, int, int, int, int, 
                 SphParticle<ndim> *, Sph<ndim> *, FLOAT);
  void UpdateAllSphProperties(int, int, SphParticle<ndim> *, 
                 Sph<ndim> *, Nbody<ndim> *);
  void UpdateAllSphForces(int, int, SphParticle<ndim> *, 
                          Sph<ndim> *, Nbody<ndim> *);
  void UpdateAllSphHydroForces(int, int, SphParticle<ndim> *, 
                               Sph<ndim> *, Nbody<ndim> *);
  void UpdateAllSphGravForces(int, int, SphParticle<ndim> *, 
                              Sph<ndim> *, Nbody<ndim> *);
  void UpdateAllSphDudt(int, int, SphParticle<ndim> *, Sph<ndim> *);
  void UpdateAllSphDerivatives(int, int, SphParticle<ndim> *, Sph<ndim> *);
  void UpdateActiveParticleCounters(SphParticle<ndim> *, Sph<ndim> *);
  void UpdateAllStarGasForces(int, int, SphParticle<ndim> *, 
                              Sph<ndim> *, Nbody<ndim> *);


  // Additional functions for binary tree neighbour search
  //---------------------------------------------------------------------------
  void AllocateMemory(Sph<ndim> *);
  void DeallocateMemory(void);


  // Additional variables for grid
  //---------------------------------------------------------------------------
  string gravity_mac;               ///< Multipole-acceptance criteria for tree
  string multipole;                 ///< Multipole-order for cell gravity
  int Ncell;                        ///< Current no. of grid cells
  int Ncellmax;                     ///< Max. allowed no. of grid cells
  int Nleafmax;                     ///< Max. number of particles per leaf cell
  int Nlistmax;                     ///< Max. length of neighbour list
  int Ntot;                         ///< No. of current points in list
  int Ntotold;                      ///< Prev. no. of particles
  int Ntotmax;                      ///< Max. no. of points in list
  int Ntotmaxold;                   ///< Old value of Ntotmax
  FLOAT macerror;                   ///< Error tolerance for gravity tree-MAC
  FLOAT hmax;                       ///< Store hmax in the tree
  FLOAT theta;                      ///< Geometric opening angle
  FLOAT thetamaxsqd;                ///< Geometric opening angle squared
  FLOAT invthetamaxsqd;             ///< 1 / thetamaxsqd
  KDTree<ndim,ParticleType> *tree;  ///< Pointer to tree

  bool allocated_buffer;            ///< ..
  int Nthreads;                     ///< ..
  int *Nneibmaxbuf;                 ///< ..
  int *Ndirectmaxbuf;               ///< ..
  int *Ngravcellmaxbuf;             ///< ..
  int **activelistbuf;              ///< ..
  int **levelneibbuf;               ///< ..
  ParticleType<ndim> **neibpartbuf;   // Local copy of neighbouring ptcls
  ParticleType<ndim> **activepartbuf; // Local copy of SPH particle  

};


#endif
