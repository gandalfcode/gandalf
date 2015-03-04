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
#if defined MPI_PARALLEL
#include "MpiExport.h"
#include "MpiNode.h"
#endif
using namespace std;




//=================================================================================================
//  Class NeighbourSearch
/// \brief   SphNeighbourSearch class definition.
/// \details Class for creating the SPH neighbour search data structure, and for computing local
///          neighbour lists and calling SPH functions (e.g. computing h, SPH forces, etc..).
/// \author  D. A. Hubber, G. Rosotti
/// \date    03/04/2013
//=================================================================================================
template <int ndim>
class MeshlessFVNeighbourSearch
{
#if defined MPI_PARALLEL
protected:
  vector<int> ids_active_particles;
#endif
 public:

  //-----------------------------------------------------------------------------------------------
  MeshlessFVNeighbourSearch(FLOAT kernrangeaux, DomainBox<ndim> *boxaux,
                            SmoothingKernel<ndim> *kernaux, CodeTiming *timingaux) :
   kernrange(kernrangeaux),
   kernrangesqd(kernrangeaux*kernrangeaux),
   box(boxaux),
   kernp(kernaux),
   timing(timingaux) {};

  ~MeshlessFVNeighbourSearch() {};


  //-----------------------------------------------------------------------------------------------
  virtual void BuildTree(bool, int, int, int, int, int,
                         MeshlessFVParticle<ndim> *, MeshlessFV<ndim> *, FLOAT) = 0;
  virtual void BuildGhostTree(bool, int, int, int, int, int,
                              MeshlessFVParticle<ndim> *, MeshlessFV<ndim> *, FLOAT) = 0;
  virtual int GetGatherNeighbourList(FLOAT *, FLOAT, MeshlessFVParticle<ndim> *, int, int, int *) = 0;
  virtual void UpdateAllProperties(int, int, MeshlessFVParticle<ndim> *, MeshlessFV<ndim> *, Nbody<ndim> *) = 0;
  virtual void UpdateGradientMatrices(int, int, MeshlessFVParticle<ndim> *, MeshlessFV<ndim> *, Nbody<ndim> *) = 0;
  virtual void UpdateGodunovFluxes(int, int, MeshlessFVParticle<ndim> *, MeshlessFV<ndim> *, Nbody<ndim> *) = 0;


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
//  Class MeshlessFVBruteForceSearch
/// \brief   ..
/// \details Class for computing SPH neighbour lists using brute force only
///          (i.e. direct summation over all particles).
/// \author  D. A. Hubber, G. Rosotti
/// \date    03/04/2013
//=================================================================================================
template <int ndim, template<int> class ParticleType>
class MeshlessFVBruteForce: public MeshlessFVNeighbourSearch<ndim>
{
 public:

  using MeshlessFVNeighbourSearch<ndim>::neibcheck;
  using MeshlessFVNeighbourSearch<ndim>::timing;
  using MeshlessFVNeighbourSearch<ndim>::kernp;
  using MeshlessFVNeighbourSearch<ndim>::kernfac;
  using MeshlessFVNeighbourSearch<ndim>::kernrange;
  using MeshlessFVNeighbourSearch<ndim>::kernrangesqd;


  //-----------------------------------------------------------------------------------------------
  MeshlessFVBruteForce(FLOAT, DomainBox<ndim> *, SmoothingKernel<ndim> *, CodeTiming *);
  ~MeshlessFVBruteForce();


  //-----------------------------------------------------------------------------------------------
  void BuildTree(bool, int, int, int, int, int, MeshlessFVParticle<ndim> *, MeshlessFV<ndim> *, FLOAT) {};
  void BuildGhostTree(bool, int, int, int, int, int, MeshlessFVParticle<ndim> *, MeshlessFV<ndim> *, FLOAT) {};
  int GetGatherNeighbourList(FLOAT *, FLOAT, MeshlessFVParticle<ndim> *, int, int, int *) {return 0;};
  void UpdateAllProperties(int, int, MeshlessFVParticle<ndim> *, MeshlessFV<ndim> *, Nbody<ndim> *);
  void UpdateGradientMatrices(int, int, MeshlessFVParticle<ndim> *, MeshlessFV<ndim> *, Nbody<ndim> *);
  void UpdateGodunovFluxes(int, int, MeshlessFVParticle<ndim> *, MeshlessFV<ndim> *, Nbody<ndim> *);

};
#endif
