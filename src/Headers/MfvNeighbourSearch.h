//=================================================================================================
//  MfvNeighbourSearch.h
//  Header file containing class definitions for all Meshless Finite-Volume neighbour searching
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


#ifndef _MFV_NEIGHBOUR_SEARCH_H_
#define _MFV_NEIGHBOUR_SEARCH_H_


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
#include "NeighbourSearch.h"
#include "SmoothingKernel.h"
#include "Particle.h"
#include "MeshlessFV.h"
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




//=================================================================================================
//  Class MeshlessFVNeighbourSearch
/// \brief   MeshlessFVNeighbourSearch class definition.
/// \details MeshlessFVNeighbourSearch class definition.
/// \author  D. A. Hubber, G. Rosotti
/// \date    03/04/2013
//=================================================================================================
template <int ndim>
class MeshlessFVNeighbourSearch : public virtual NeighbourSearch<ndim>
{
#if defined MPI_PARALLEL
protected:
  vector<int> ids_active_particles;
#endif
 public:


  //-----------------------------------------------------------------------------------------------
  MeshlessFVNeighbourSearch(FLOAT kernrangeaux, DomainBox<ndim> *boxaux,
                            SmoothingKernel<ndim> *kernaux, CodeTiming *timingaux) :
    NeighbourSearch<ndim>(kernrangeaux, boxaux, kernaux, timingaux) {};
  virtual ~MeshlessFVNeighbourSearch() {};


  //-----------------------------------------------------------------------------------------------
  virtual void UpdateAllProperties(int, int, MeshlessFVParticle<ndim> *, MeshlessFV<ndim> *,
                                   Nbody<ndim> *, DomainBox<ndim> &) = 0;
  virtual void UpdateGradientMatrices(int, int, MeshlessFVParticle<ndim> *, MeshlessFV<ndim> *,
                                      Nbody<ndim> *, DomainBox<ndim> &) = 0;
  virtual void UpdateGodunovFluxes(int, int, FLOAT, MeshlessFVParticle<ndim> *,
                                   MeshlessFV<ndim> *, Nbody<ndim> *, DomainBox<ndim> &) = 0;
  virtual void UpdateAllGravForces(int, int, MeshlessFVParticle<ndim> *, MeshlessFV<ndim> *,
                                   Nbody<ndim> *, DomainBox<ndim> &, Ewald<ndim> *) = 0;

};



//=================================================================================================
//  Class MeshlessFVBruteForceSearch
/// \brief   MeshlessFVBruteForceSearch class definition.
/// \details MeshlessFVBruteForceSearch class definition.
/// \author  D. A. Hubber
/// \date    21/04/2015
//=================================================================================================
template <int ndim, template<int> class ParticleType>
class MeshlessFVBruteForce : public MeshlessFVNeighbourSearch<ndim>, public BruteForceSearch<ndim,ParticleType>
{
 public:

  using NeighbourSearch<ndim>::neibcheck;
  using NeighbourSearch<ndim>::timing;
  using NeighbourSearch<ndim>::kernp;
  using NeighbourSearch<ndim>::kernfac;
  using NeighbourSearch<ndim>::kernrange;
  using NeighbourSearch<ndim>::kernrangesqd;


  //-----------------------------------------------------------------------------------------------
  MeshlessFVBruteForce(FLOAT, DomainBox<ndim> *, SmoothingKernel<ndim> *, CodeTiming *);
  virtual ~MeshlessFVBruteForce();


  //-----------------------------------------------------------------------------------------------
  virtual void UpdateAllProperties(int, int, MeshlessFVParticle<ndim> *, MeshlessFV<ndim> *,
                                   Nbody<ndim> *, DomainBox<ndim> &);
  virtual void UpdateGradientMatrices(int, int, MeshlessFVParticle<ndim> *, MeshlessFV<ndim> *,
                                      Nbody<ndim> *, DomainBox<ndim> &);
  virtual void UpdateGodunovFluxes(int, int, FLOAT, MeshlessFVParticle<ndim> *,
                                   MeshlessFV<ndim> *, Nbody<ndim> *, DomainBox<ndim> &);
  virtual void UpdateAllGravForces(int, int, MeshlessFVParticle<ndim> *, MeshlessFV<ndim> *,
                                   Nbody<ndim> *, DomainBox<ndim> &, Ewald<ndim> *);

};



//=================================================================================================
//  Class MeshlessFVTree
/// \brief   MeshlessFVTree class definition.
/// \details MeshlessFVTree class definition.
/// \author  D. A. Hubber
/// \date    21/04/2015
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
class MeshlessFVTree: public MeshlessFVNeighbourSearch<ndim>, public HydroTree<ndim,ParticleType,TreeCell>
{
#if defined MPI_PARALLEL
  vector<vector<int> > ids_sent_particles;
protected:
  using NeighbourSearch<ndim>::ids_active_particles;
  using NeighbourSearch<ndim>::N_imported_part_per_proc;
#endif
 public:

  using HydroTree<ndim,ParticleType,TreeCell>::activelistbuf;
  using HydroTree<ndim,ParticleType,TreeCell>::activepartbuf;
  using HydroTree<ndim,ParticleType,TreeCell>::allocated_buffer;
  using HydroTree<ndim,ParticleType,TreeCell>::box;
  using HydroTree<ndim,ParticleType,TreeCell>::cellbuf;
  using HydroTree<ndim,ParticleType,TreeCell>::gravity_mac;
  using HydroTree<ndim,ParticleType,TreeCell>::kernp;
  using HydroTree<ndim,ParticleType,TreeCell>::kernrange;
  using HydroTree<ndim,ParticleType,TreeCell>::kernrangesqd;
  using HydroTree<ndim,ParticleType,TreeCell>::levelneibbuf;
  using HydroTree<ndim,ParticleType,TreeCell>::multipole;
  using HydroTree<ndim,ParticleType,TreeCell>::neibcheck;
  using HydroTree<ndim,ParticleType,TreeCell>::neibpartbuf;
  using HydroTree<ndim,ParticleType,TreeCell>::Ngravcellmaxbuf;
  using HydroTree<ndim,ParticleType,TreeCell>::Nleafmax;
  using HydroTree<ndim,ParticleType,TreeCell>::Nneibmaxbuf;
  using HydroTree<ndim,ParticleType,TreeCell>::Ntot;
  using HydroTree<ndim,ParticleType,TreeCell>::Ntotmax;
  using HydroTree<ndim,ParticleType,TreeCell>::Ntotmaxold;
  using HydroTree<ndim,ParticleType,TreeCell>::Ntotold;
  using HydroTree<ndim,ParticleType,TreeCell>::timing;
  using HydroTree<ndim,ParticleType,TreeCell>::tree;
  using HydroTree<ndim,ParticleType,TreeCell>::ghosttree;
#ifdef MPI_PARALLEL
  using HydroTree<ndim,ParticleType,TreeCell>::mpighosttree;
  using HydroTree<ndim,ParticleType,TreeCell>::Nmpi;
  using HydroTree<ndim,ParticleType,TreeCell>::prunedtree;
  using HydroTree<ndim,ParticleType,TreeCell>::sendprunedtree;
#endif


  //-----------------------------------------------------------------------------------------------
  MeshlessFVTree(int _Nleafmax, int _Nmpi, int _pruning_level_min, int _pruning_level_max,
                 FLOAT _thetamaxsqd, FLOAT _kernrange, FLOAT _macerror,
                 string _gravity_mac, string _multipole, DomainBox<ndim> *_box,
                 SmoothingKernel<ndim> *_kern, CodeTiming *_timing); //:
    /*NeighbourSearch<ndim>(_kernrange, _box, _kern, _timing),
    MeshlessFVNeighbourSearch<ndim>(_kernrange, _box, _kern, _timing),
    HydroTree<ndim,ParticleType,TreeCell>(_Nleafmax, _Nmpi, _thetamaxsqd, _kernrange, _macerror,
                                          _gravity_mac, _multipole, _box, _kern, _timing) {};*/
  virtual ~MeshlessFVTree(); //{};


  //-----------------------------------------------------------------------------------------------
  virtual void UpdateAllProperties(int, int, MeshlessFVParticle<ndim> *, MeshlessFV<ndim> *,
                                   Nbody<ndim> *, DomainBox<ndim> &);
  virtual void UpdateGradientMatrices(int, int, MeshlessFVParticle<ndim> *, MeshlessFV<ndim> *,
                                      Nbody<ndim> *, DomainBox<ndim> &);
  virtual void UpdateGodunovFluxes(int, int, FLOAT, MeshlessFVParticle<ndim> *,
                                   MeshlessFV<ndim> *, Nbody<ndim> *, DomainBox<ndim> &);
  virtual void UpdateAllGravForces(int, int, MeshlessFVParticle<ndim> *, MeshlessFV<ndim> *,
                                   Nbody<ndim> *, DomainBox<ndim> &, Ewald<ndim> *);

};



//=================================================================================================
//  Class MeshlessFVKDTree
/// \brief   MeshlessFVKDTree class definition.
/// \details MeshlessFVKDTree class definition.
/// \author  D. A. Hubber
/// \date    17/09/2014
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
class MeshlessFVKDTree: public MeshlessFVTree<ndim,ParticleType,TreeCell>
{
 public:

  using MeshlessFVTree<ndim,ParticleType,TreeCell>::tree;
  using MeshlessFVTree<ndim,ParticleType,TreeCell>::ghosttree;
#ifdef MPI_PARALLEL
  using MeshlessFVTree<ndim,ParticleType,TreeCell>::mpighosttree;
  using MeshlessFVTree<ndim,ParticleType,TreeCell>::Nmpi;
  using MeshlessFVTree<ndim,ParticleType,TreeCell>::prunedtree;
  using MeshlessFVTree<ndim,ParticleType,TreeCell>::sendprunedtree;
#endif


  //-----------------------------------------------------------------------------------------------
  MeshlessFVKDTree(int, int, int, int, FLOAT, FLOAT, FLOAT, string, string,
                   DomainBox<ndim> *, SmoothingKernel<ndim> *, CodeTiming *);

};



//=================================================================================================
//  Class MeshlessFVOctTree
/// \brief   MeshlessFVOctTree class definition.
/// \details MeshlessFVOctTree class definition.
/// \author  D. A. Hubber
/// \date    17/09/2014
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
class MeshlessFVOctTree: public MeshlessFVTree<ndim,ParticleType,TreeCell>
{
 public:

  using MeshlessFVTree<ndim,ParticleType,TreeCell>::tree;
  using MeshlessFVTree<ndim,ParticleType,TreeCell>::ghosttree;
#ifdef MPI_PARALLEL
  using MeshlessFVTree<ndim,ParticleType,TreeCell>::mpighosttree;
  using MeshlessFVTree<ndim,ParticleType,TreeCell>::Nmpi;
  using MeshlessFVTree<ndim,ParticleType,TreeCell>::prunedtree;
  using MeshlessFVTree<ndim,ParticleType,TreeCell>::sendprunedtree;
#endif


  //-----------------------------------------------------------------------------------------------
  MeshlessFVOctTree(int, int, int, int, FLOAT, FLOAT, FLOAT, string, string,
                    DomainBox<ndim> *, SmoothingKernel<ndim> *, CodeTiming *);

};
#endif
