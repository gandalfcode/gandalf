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
#include "GhostNeighbours.hpp"
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
#include "BruteForceTree.h"
#include "NeighbourManager.h"
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
 public:

  //-----------------------------------------------------------------------------------------------
  virtual ~SphNeighbourSearch() {};


  //-----------------------------------------------------------------------------------------------
  virtual void UpdateAllProperties(Hydrodynamics<ndim>* hydro, Nbody<ndim>* nbody,
                                   DomainBox<ndim> &simbox) {
    UpdateAllSphProperties(static_cast<Sph<ndim>*>(hydro), nbody, simbox) ;
  }
  virtual void UpdateAllSphProperties(Sph<ndim> *, Nbody<ndim> *,DomainBox<ndim> &) = 0;
  virtual void UpdateAllSphHydroForces(Sph<ndim> *, Nbody<ndim> *, DomainBox<ndim> &) =0;
  virtual void UpdateAllSphForces(Sph<ndim> *, Nbody<ndim> *,
                                  DomainBox<ndim> &, Ewald<ndim> *) = 0;

};


//=================================================================================================
//  Class SphTree
/// \brief   Class containing tree for efficient SPH neighbour searching and gravity calculations.
/// \details Class containing tree for efficient SPH neighbour searching and gravity calculations.
/// \author  D. A. Hubber
/// \date    08/01/2014
//=================================================================================================
template <int ndim, template<int> class ParticleType>
class SphTree : public SphNeighbourSearch<ndim>, public HydroTree<ndim,ParticleType>
{
#if defined MPI_PARALLEL
protected:
#endif
 public:

  using HydroTree<ndim,ParticleType>::activelistbuf;
  using HydroTree<ndim,ParticleType>::activepartbuf;
  using HydroTree<ndim,ParticleType>::allocated_buffer;
  using HydroTree<ndim,ParticleType>::box;
  using HydroTree<ndim,ParticleType>::gravity_mac;
  using HydroTree<ndim,ParticleType>::kernp;
  using HydroTree<ndim,ParticleType>::kernrange;
  using HydroTree<ndim,ParticleType>::kernrangesqd;
  using HydroTree<ndim,ParticleType>::levelneibbuf;
  using HydroTree<ndim,ParticleType>::multipole;
  using HydroTree<ndim,ParticleType>::neibcheck;
  //using HydroTree<ndim,ParticleType>::neibpartbuf;
  using HydroTree<ndim,ParticleType>::Ngravcellmaxbuf;
  using HydroTree<ndim,ParticleType>::Nleafmax;
  using HydroTree<ndim,ParticleType>::Nneibmaxbuf;
  using HydroTree<ndim,ParticleType>::Nthreads;
  using HydroTree<ndim,ParticleType>::Ntot;
  using HydroTree<ndim,ParticleType>::Ntotmax;
  using HydroTree<ndim,ParticleType>::Ntotmaxold;
  using HydroTree<ndim,ParticleType>::Ntotold;
  using HydroTree<ndim,ParticleType>::timing;
  using HydroTree<ndim,ParticleType>::tree;
  using HydroTree<ndim,ParticleType>::ghosttree;
#ifdef MPI_PARALLEL
  using HydroTree<ndim,ParticleType>::mpighosttree;
  using HydroTree<ndim,ParticleType>::Nmpi;
  using HydroTree<ndim,ParticleType>::prunedtree;
  using HydroTree<ndim,ParticleType>::sendprunedtree;
#endif


  //-----------------------------------------------------------------------------------------------
  SphTree(string tree_t, int _Nleafmax, int _Nmpi, int _pruning_level_min, int _pruning_level_max,
          FLOAT _thetamaxsqd, FLOAT _kernrange, FLOAT _macerror,
          string _gravity_mac, multipole_method _multipole,
          DomainBox<ndim> *_box, SmoothingKernel<ndim> *_kern, CodeTiming *_timing,
          ParticleTypeRegister& types) :
    HydroTree<ndim,ParticleType>(tree_t, _Nleafmax, _Nmpi, _pruning_level_min, _pruning_level_max,
                                          _thetamaxsqd, _kernrange, _macerror, _gravity_mac,
                                          _multipole, _box, _kern, _timing, types) {} ;
  virtual ~SphTree() {};


  //-----------------------------------------------------------------------------------------------
  virtual void UpdateAllSphProperties(Sph<ndim> *, Nbody<ndim> *, DomainBox<ndim> &) = 0;
  virtual void UpdateAllSphHydroForces(Sph<ndim> *, Nbody<ndim> *, DomainBox<ndim> &) =0;
  virtual void UpdateAllSphForces(Sph<ndim> *, Nbody<ndim> *,
                                  DomainBox<ndim> &, Ewald<ndim> *) = 0;

};



//=================================================================================================
//  Class GradhSphTree
/// \brief   Class containing tree for computing grad-h SPH summation and force loops.
/// \details Class containing tree for computing grad-h SPH summation and force loops.
/// \author  D. A. Hubber
/// \date    08/01/2014
//=================================================================================================
template <int ndim, template<int> class ParticleType>
class GradhSphTree : public SphTree<ndim,ParticleType>
{
private:
	  typedef typename ParticleType<ndim>::HydroForcesParticle HydroParticle;
	  typedef NeighbourManager<ndim, HydroParticle > NeighbourManagerHydro;
	  vector<NeighbourManagerHydro> neibmanagerbufhydro;

	  typedef typename ParticleType<ndim>::DensityParticle DensityParticle;
	  typedef NeighbourManager<ndim, DensityParticle > NeighbourManagerDensity;
	  vector<NeighbourManagerDensity> neibmanagerbufdens;
 public:

  using SphTree<ndim,ParticleType>::activelistbuf;
  using SphTree<ndim,ParticleType>::activepartbuf;
  using SphTree<ndim,ParticleType>::allocated_buffer;
  using SphTree<ndim,ParticleType>::box;
  using SphTree<ndim,ParticleType>::gravity_mac;
  using SphTree<ndim,ParticleType>::kernp;
  using SphTree<ndim,ParticleType>::kernrange;
  using SphTree<ndim,ParticleType>::kernrangesqd;
  using SphTree<ndim,ParticleType>::levelneibbuf;
  using SphTree<ndim,ParticleType>::multipole;
  using SphTree<ndim,ParticleType>::neibcheck;
  //using SphTree<ndim,ParticleType>::neibpartbuf;
  using SphTree<ndim,ParticleType>::Ngravcellmaxbuf;
  using SphTree<ndim,ParticleType>::Nleafmax;
  using SphTree<ndim,ParticleType>::Nneibmaxbuf;
  using SphTree<ndim,ParticleType>::Nthreads;
  using SphTree<ndim,ParticleType>::Ntot;
  using SphTree<ndim,ParticleType>::Ntotmax;
  using SphTree<ndim,ParticleType>::Ntotmaxold;
  using SphTree<ndim,ParticleType>::Ntotold;
  using SphTree<ndim,ParticleType>::timing;
  using SphTree<ndim,ParticleType>::tree;
  using SphTree<ndim,ParticleType>::ghosttree;
  using HydroTree<ndim,ParticleType>::neiblistbuf;
  using HydroTree<ndim,ParticleType>::ptypebuf;
#ifdef MPI_PARALLEL
  using SphTree<ndim,ParticleType>::mpighosttree;
  using SphTree<ndim,ParticleType>::Nmpi;
  using SphTree<ndim,ParticleType>::prunedtree;
  using SphTree<ndim,ParticleType>::sendprunedtree;
#endif


  //-----------------------------------------------------------------------------------------------
  GradhSphTree(string, int, int, int, int, FLOAT, FLOAT, FLOAT, string, multipole_method,
               DomainBox<ndim> *, SmoothingKernel<ndim> *, CodeTiming *,
               ParticleTypeRegister& types);
  virtual ~GradhSphTree();


  //-----------------------------------------------------------------------------------------------
  virtual void UpdateAllSphProperties(Sph<ndim> *, Nbody<ndim> *, DomainBox<ndim> &);
  virtual void UpdateAllSphHydroForces(Sph<ndim> *, Nbody<ndim> *, DomainBox<ndim> &);
  virtual void UpdateAllSphForces(Sph<ndim> *, Nbody<ndim> *,
                                  DomainBox<ndim> &, Ewald<ndim> *);

};

//=================================================================================================
//  Class SM2012SphTree
/// \brief   Class containing tree for computing Saitoh & Makino (2012) SPH force loops.
/// \details Class containing tree for computing Saitoh & Makino (2012) SPH force loops.
/// \author  D. A. Hubber
/// \date    08/01/2014
//=================================================================================================
template <int ndim, template<int> class ParticleType>
class SM2012SphTree: public SphTree<ndim,ParticleType>
{
 public:

  using SphTree<ndim,ParticleType>::activelistbuf;
  using SphTree<ndim,ParticleType>::activepartbuf;
  using SphTree<ndim,ParticleType>::allocated_buffer;
  using SphTree<ndim,ParticleType>::box;
  using SphTree<ndim,ParticleType>::gravity_mac;
  using SphTree<ndim,ParticleType>::kernp;
  using SphTree<ndim,ParticleType>::kernrange;
  using SphTree<ndim,ParticleType>::kernrangesqd;
  using SphTree<ndim,ParticleType>::levelneibbuf;
  using SphTree<ndim,ParticleType>::multipole;
  using SphTree<ndim,ParticleType>::neibcheck;
  //using SphTree<ndim,ParticleType>::neibpartbuf;
  using SphTree<ndim,ParticleType>::Ngravcellmaxbuf;
  using SphTree<ndim,ParticleType>::Nleafmax;
  using SphTree<ndim,ParticleType>::Nneibmaxbuf;
  using SphTree<ndim,ParticleType>::Ntot;
  using SphTree<ndim,ParticleType>::Ntotmax;
  using SphTree<ndim,ParticleType>::Ntotmaxold;
  using SphTree<ndim,ParticleType>::Ntotold;
  using SphTree<ndim,ParticleType>::timing;
  using SphTree<ndim,ParticleType>::tree;
  using SphTree<ndim,ParticleType>::ghosttree;
#ifdef MPI_PARALLEL
  using SphTree<ndim,ParticleType>::mpighosttree;
  using SphTree<ndim,ParticleType>::Nmpi;
  using SphTree<ndim,ParticleType>::prunedtree;
  using SphTree<ndim,ParticleType>::sendprunedtree;
#endif


  //-----------------------------------------------------------------------------------------------
  SM2012SphTree(string, int, int, int, int, FLOAT, FLOAT, FLOAT, string, multipole_method,
                DomainBox<ndim> *, SmoothingKernel<ndim> *, CodeTiming *,
                ParticleTypeRegister&);


  //-----------------------------------------------------------------------------------------------
  virtual void UpdateAllSphProperties(Sph<ndim> *, Nbody<ndim> *, DomainBox<ndim> &){};
  virtual void UpdateAllSphHydroForces(Sph<ndim> *, Nbody<ndim> *, DomainBox<ndim> &){};
  virtual void UpdateAllSphForces(Sph<ndim> *, Nbody<ndim> *,
                                  DomainBox<ndim> &, Ewald<ndim> *){};

};

#endif
