//=================================================================================================
//  SphNeighbourSearch.h
//  Header file containing class definitions for all SPH neighbour searching algorithms.
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
#include "NeighbourSearch.h"
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
class SphNeighbourSearch : public virtual NeighbourSearch<ndim>
{
#if defined MPI_PARALLEL
protected:
  vector<int> ids_active_particles;
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
  SphNeighbourSearch(FLOAT kernrangeaux, DomainBox<ndim> *boxaux,
                     SmoothingKernel<ndim> *kernaux, CodeTiming *timingaux) :
    NeighbourSearch<ndim>(kernrangeaux, boxaux, kernaux, timingaux) {};
  virtual ~SphNeighbourSearch() {};


  //-----------------------------------------------------------------------------------------------
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
  //virtual void UpdateAllStarGasForces(int, int, SphParticle<ndim> *,
  //                                    Sph<ndim> *, Nbody<ndim> *) = 0;

};



//=================================================================================================
//  Class SphBruteForceSearch
/// \brief   Class for computing SPH neighbour lists using brute force only.
/// \details Class for computing SPH neighbour lists using brute force only
///          (i.e. direct summation over all particles).
/// \author  D. A. Hubber, G. Rosotti
/// \date    03/04/2013
//=================================================================================================
template <int ndim, template<int> class ParticleType>
class SphBruteForceSearch : public SphNeighbourSearch<ndim>, public BruteForceSearch<ndim,ParticleType>
{
 public:

  using NeighbourSearch<ndim>::neibcheck;
  using NeighbourSearch<ndim>::timing;
  using NeighbourSearch<ndim>::kernp;
  using NeighbourSearch<ndim>::kernfac;
  using NeighbourSearch<ndim>::kernrange;
  using NeighbourSearch<ndim>::kernrangesqd;


  //-----------------------------------------------------------------------------------------------
  SphBruteForceSearch(FLOAT, DomainBox<ndim> *, SmoothingKernel<ndim> *, CodeTiming *);
  virtual ~SphBruteForceSearch();


  //-----------------------------------------------------------------------------------------------
  void UpdateAllSphProperties(int, int, SphParticle<ndim> *, Sph<ndim> *, Nbody<ndim> *);
  void UpdateAllSphForces(int, int, SphParticle<ndim> *, Sph<ndim> *, Nbody<ndim> *);
  void UpdateAllSphHydroForces(int, int, SphParticle<ndim> *, Sph<ndim> *, Nbody<ndim> *);
  void UpdateAllSphGravForces(int, int, SphParticle<ndim> *, Sph<ndim> *, Nbody<ndim> *);
  void UpdateAllSphPeriodicForces(int, int, SphParticle<ndim> *, Sph<ndim> *,
                                  Nbody<ndim> *, DomainBox<ndim> &, Ewald<ndim> *);
  void UpdateAllSphPeriodicGravForces(int, int, SphParticle<ndim> *, Sph<ndim> *,
                                      Nbody<ndim> *, DomainBox<ndim> &, Ewald<ndim> *);
  //void UpdateAllStarGasForces(int, int, SphParticle<ndim> *, Sph<ndim> *, Nbody<ndim> *);

#ifdef MPI_PARALLEL
  using NeighbourSearch<ndim>::ids_active_particles;
  void UpdateGravityExportList(int, int, int, SphParticle<ndim> *, Sph<ndim> *, Nbody<ndim> *);
  void UpdateHydroExportList(int, int, int, SphParticle<ndim> *,  Sph<ndim> *, Nbody<ndim> *);
#endif

};



//=================================================================================================
//  Class GradhSphBruteForce
/// \brief   Class for computing neighbour lists using brute force only for grad-h SPH.
/// \details Class for computing SPH neighbour lists using brute force only
///          (i.e. direct summation over all particles) for grad-h SPH.
/// \author  D. A. Hubber, G. Rosotti
/// \date    12/05/2014
//=================================================================================================
template <int ndim, template<int> class ParticleType>
class GradhSphBruteForce: public SphBruteForceSearch<ndim,ParticleType>
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
  virtual ~GradhSphBruteForce();


  //-----------------------------------------------------------------------------------------------
  void UpdateAllSphProperties(int, int, SphParticle<ndim> *, Sph<ndim> *, Nbody<ndim> *);

};



//=================================================================================================
//  Class SM2012SphBruteForce
/// \brief   Class for computing neighbour lists using brute force only for SM2012 SPH.
/// \details Class for computing neighbour lists using brute force only
///          (i.e. direct summation over all particles) for SM2012 SPH.
/// \author  D. A. Hubber, G. Rosotti
/// \date    12/05/2014
//=================================================================================================
template <int ndim, template<int> class ParticleType>
class SM2012SphBruteForce: public SphBruteForceSearch<ndim,ParticleType>
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
  virtual ~SM2012SphBruteForce();


  //-----------------------------------------------------------------------------------------------
  void UpdateAllSphProperties(int, int, SphParticle<ndim> *, Sph<ndim> *, Nbody<ndim> *);

};



//=================================================================================================
//  Class SphTree
/// \brief   Class containing tree for efficient SPH neighbour searching and gravity calculations.
/// \details Class containing tree for efficient SPH neighbour searching and gravity calculations.
/// \author  D. A. Hubber
/// \date    08/01/2014
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
class SphTree: public SphNeighbourSearch<ndim>, public HydroTree<ndim,ParticleType,TreeCell>
{
#if defined MPI_PARALLEL
  vector<vector<int> > ids_sent_particles;
protected:
  using NeighbourSearch<ndim>::ids_active_particles;
  vector<int> N_imported_part_per_proc;
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
  using HydroTree<ndim,ParticleType,TreeCell>::ghostprunedtree;
  using HydroTree<ndim,ParticleType,TreeCell>::mpighosttree;
  using HydroTree<ndim,ParticleType,TreeCell>::Nghostpruned;
  using HydroTree<ndim,ParticleType,TreeCell>::Nghostprunedmax;
  using HydroTree<ndim,ParticleType,TreeCell>::Nmpi;
  using HydroTree<ndim,ParticleType,TreeCell>::prunedtree;
#endif


  //-----------------------------------------------------------------------------------------------
  SphTree(int _Nleafmax, int _Nmpi, FLOAT _thetamaxsqd, FLOAT _kernrange, FLOAT _macerror,
          string _gravity_mac, string _multipole, DomainBox<ndim> *_box,
          SmoothingKernel<ndim> *_kern, CodeTiming *_timing) :
    NeighbourSearch<ndim>(_kernrange, _box, _kern, _timing),
    SphNeighbourSearch<ndim>(_kernrange, _box, _kern, _timing),
    HydroTree<ndim,ParticleType,TreeCell>(_Nleafmax, _Nmpi, _thetamaxsqd, _kernrange, _macerror,
                                          _gravity_mac, _multipole, _box, _kern, _timing) {};
  virtual ~SphTree() {};


  //-----------------------------------------------------------------------------------------------
  void UpdateAllSphProperties(int, int, SphParticle<ndim> *, Sph<ndim> *, Nbody<ndim> *) {};
  void UpdateAllSphForces(int, int, SphParticle<ndim> *, Sph<ndim> *, Nbody<ndim> *) {};
  void UpdateAllSphHydroForces(int, int, SphParticle<ndim> *, Sph<ndim> *, Nbody<ndim> *) {};
  void UpdateAllSphGravForces(int, int, SphParticle<ndim> *, Sph<ndim> *, Nbody<ndim> *) {};
  void UpdateAllSphPeriodicForces(int, int, SphParticle<ndim> *, Sph<ndim> *,
                                  Nbody<ndim> *, DomainBox<ndim> &, Ewald<ndim> *) {};
  void UpdateAllSphPeriodicGravForces(int, int, SphParticle<ndim> *, Sph<ndim> *,
                                      Nbody<ndim> *, DomainBox<ndim> &, Ewald<ndim> *) {};
  //void UpdateAllStarGasForces(int, int, SphParticle<ndim> *, Sph<ndim> *, Nbody<ndim> *) {};

};



//=================================================================================================
//  Class GradhSphTree
/// \brief   Class containing tree for computing grad-h SPH summation and force loops.
/// \details Class containing tree for computing grad-h SPH summation and force loops.
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
  virtual ~GradhSphTree();


  //-----------------------------------------------------------------------------------------------
  void UpdateAllSphProperties(int, int, SphParticle<ndim> *, Sph<ndim> *, Nbody<ndim> *);
  void UpdateAllSphForces(int, int, SphParticle<ndim> *, Sph<ndim> *, Nbody<ndim> *);
  void UpdateAllSphHydroForces(int, int, SphParticle<ndim> *, Sph<ndim> *, Nbody<ndim> *);
  void UpdateAllSphGravForces(int, int, SphParticle<ndim> *, Sph<ndim> *, Nbody<ndim> *);
  //void UpdateAllStarGasForces(int, int, SphParticle<ndim> *, Sph<ndim> *, Nbody<ndim> *);
  void UpdateAllSphPeriodicForces(int, int, SphParticle<ndim> *, Sph<ndim> *,
                                  Nbody<ndim> *, DomainBox<ndim> &, Ewald<ndim> *);
  void UpdateAllSphPeriodicGravForces(int, int, SphParticle<ndim> *, Sph<ndim> *,
                                      Nbody<ndim> *, DomainBox<ndim> &, Ewald<ndim> *);

};



//=================================================================================================
//  Class GradhSphKDTree
/// \brief   Grad-h SPH neighbour searching class using the KD-tree.
/// \details Grad-h SPH neighbour searching class using the KD-tree.
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
/// \details Class containing octal tree for computing grad-h SPH force loops.
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
/// \details Class containing tree for computing Saitoh & Makino (2012) SPH force loops.
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


  //-----------------------------------------------------------------------------------------------
  void UpdateAllSphProperties(int, int, SphParticle<ndim> *, Sph<ndim> *, Nbody<ndim> *) {};
  void UpdateAllSphForces(int, int, SphParticle<ndim> *, Sph<ndim> *, Nbody<ndim> *) {};
  void UpdateAllSphHydroForces(int, int, SphParticle<ndim> *, Sph<ndim> *, Nbody<ndim> *) {};
  void UpdateAllSphGravForces(int, int, SphParticle<ndim> *, Sph<ndim> *, Nbody<ndim> *) {};
  //void UpdateAllStarGasForces(int, int, SphParticle<ndim> *, Sph<ndim> *, Nbody<ndim> *) {};

};



//=================================================================================================
//  Class SM2012SphKDTree
/// \brief   Class containing kd-tree for computing SM2012 SPH force loops.
/// \details Class containing kd-tree for computing SM2012 SPH force loops.
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
/// \brief   Class containing octal-tree for computing SM2012 SPH summation and force loops.
/// \details Class containing octal-tree for computing SM2012 SPH summation and force loops.
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
#endif
