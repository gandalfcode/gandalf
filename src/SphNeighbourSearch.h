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
#include "Nbody.h"
#include "SphKernel.h"
#include "SphParticle.h"
#include "Sph.h"
#include "DomainBox.h"
#include "Parameters.h"
#include "KDTree.h"
using namespace std;

#if defined MPI_PARALLEL
#include "MpiExport.h"
#endif


// Forward declaration of Sph to break circular dependency
template <int ndim>
class Sph;

// Forward declare MpiNode to break circular dependency
#if defined MPI_PARALLEL
template <int ndim>
class MpiNode;
#endif



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

  //---------------------------------------------------------------------------
  SphNeighbourSearch(FLOAT, DomainBox<ndim> *, 
                     SphKernel<ndim> *, CodeTiming *);
  ~SphNeighbourSearch();


  //---------------------------------------------------------------------------
  virtual void BuildTree(bool, int, int, int, int, int, 
                         SphParticle<ndim> *, Sph<ndim> *, FLOAT) = 0;
  virtual int GetGatherNeighbourList(FLOAT *, FLOAT, SphParticle<ndim> *, 
                                     int, int, int *) = 0;
  virtual void UpdateAllSphProperties(int, int, SphParticle<ndim> *, 
                                     Sph<ndim> *, Nbody<ndim> *) = 0;
  virtual void UpdateAllSphForces(int, int, SphParticle<ndim> *, 
                                  Sph<ndim> *, Nbody<ndim> *) = 0;
  virtual void UpdateAllSphHydroForces(int, int, SphParticle<ndim> *, 
                                       Sph<ndim> *, Nbody<ndim> *) = 0;
  virtual void UpdateAllSphGravForces(int, int, SphParticle<ndim> *, 
                                      Sph<ndim> *, Nbody<ndim> *) = 0;
  virtual void UpdateActiveParticleCounters(SphParticle<ndim> *, 
                                            Sph<ndim> *) = 0;
  virtual void UpdateAllStarGasForces(int, int, SphParticle<ndim> *, 
                                      Sph<ndim> *, Nbody<ndim> *) = 0;
  virtual void UpdateAllSphDudt(int, int, SphParticle<ndim> *, Sph<ndim> *) = 0;
  virtual void UpdateAllSphDerivatives(int, int, SphParticle<ndim> *, Sph<ndim> *) = 0;


  //---------------------------------------------------------------------------
  const FLOAT kernrange;            ///< Kernel extent (in units of h)
  const FLOAT kernrangesqd;         ///< Kernel extent (squared)


  //---------------------------------------------------------------------------
  CodeTiming *timing;               ///< Pointer to code timing object
  DomainBox<ndim> *box;             ///< Pointer to simulation bounding box
  SphKernel<ndim> *kernp;           ///< Pointer to SPH kernel object

  bool neibcheck;                   ///< Flag to verify neighbour lists
  FLOAT kernfac;                    ///< ..

};



//=============================================================================
//  Class BruteForceSearch
/// \brief   ..
/// \details Class for computing SPH neighbour lists using brute force only 
///          (i.e. direct summation over all particles).
/// \author  D. A. Hubber, G. Rosotti
/// \date    03/04/2013
//=============================================================================
template <int ndim, template<int> class ParticleType>
class BruteForceSearch: public SphNeighbourSearch<ndim>
{
#if defined MPI_PARALLEL
protected:
  vector<ExportParticleInfo<ndim> > particles_to_export;
  vector<int> ids_active_particles;
#endif
 public:

  using SphNeighbourSearch<ndim>::neibcheck;
  using SphNeighbourSearch<ndim>::timing;
  using SphNeighbourSearch<ndim>::kernp;
  using SphNeighbourSearch<ndim>::kernfac;
  using SphNeighbourSearch<ndim>::kernrange;
  using SphNeighbourSearch<ndim>::kernrangesqd;
  

  //---------------------------------------------------------------------------
  BruteForceSearch(FLOAT, DomainBox<ndim> *, SphKernel<ndim> *, CodeTiming *);
  ~BruteForceSearch();


  //---------------------------------------------------------------------------
  void BuildTree(bool, int, int, int, int, int, 
                 SphParticle<ndim> *, Sph<ndim> *, FLOAT);
  int GetGatherNeighbourList(FLOAT *, FLOAT, SphParticle<ndim> *, 
                             int, int, int *);
  void UpdateAllSphProperties(int, int, SphParticle<ndim> *, 
                              Sph<ndim> *, Nbody<ndim> *);
  void UpdateAllSphForces(int, int, SphParticle<ndim> *, 
                          Sph<ndim> *, Nbody<ndim> *);
  void UpdateAllSphHydroForces(int, int, SphParticle<ndim> *, 
                               Sph<ndim> *, Nbody<ndim> *);
  void UpdateAllSphGravForces(int, int, SphParticle<ndim> *, 
                              Sph<ndim> *, Nbody<ndim> *);
  void UpdateActiveParticleCounters(SphParticle<ndim> *, Sph<ndim> *);
  void UpdateAllStarGasForces(int, int, SphParticle<ndim> *, 
                              Sph<ndim> *, Nbody<ndim> *);
  void UpdateAllSphDudt(int, int, SphParticle<ndim> *, Sph<ndim> *);
  void UpdateAllSphDerivatives(int, int, SphParticle<ndim> *, Sph<ndim> *);

#if defined MPI_PARALLEL
  void FindGhostParticlesToExport(Sph<ndim>* sph, std::vector<std::vector<ParticleType<ndim>* > >&,
      const std::vector<int>&, MpiNode<ndim>*);
  void FindParticlesToTransfer(Sph<ndim>* sph, std::vector<std::vector<int> >& particles_to_export,
      std::vector<int>& all_particles_to_export, const std::vector<int>& potential_nodes, MpiNode<ndim>* mpinodes);
  vector<ExportParticleInfo<ndim> >& GetExportInfo(int Nproc, Sph<ndim>* sph);
  void UnpackExported (vector<ExportParticleInfo<ndim> >& arrays, vector<int>& N_received_particles_from_proc,
      Sph<ndim>* sph, int rank);
  void GetBackExportInfo(vector<ExportBackParticleInfo<ndim> >& arrays, vector<int>& N_exported_particles_from_proc,
      Sph<ndim>* sph, int rank);
  void UnpackReturnedExportInfo(vector<ExportBackParticleInfo<ndim> >& received_information, vector<int>& recv_displs,
      Sph<ndim>* sph, int rank);
#endif
};



//=============================================================================
//  Class GradhSphBruteForce
/// \brief   ..
/// \details Class for computing SPH neighbour lists using brute force only 
///          (i.e. direct summation over all particles).
/// \author  D. A. Hubber, G. Rosotti
/// \date    12/05/2014
//=============================================================================
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
  

  //---------------------------------------------------------------------------
  GradhSphBruteForce(FLOAT, DomainBox<ndim> *, 
                     SphKernel<ndim> *, CodeTiming *);
  ~GradhSphBruteForce();


  //---------------------------------------------------------------------------
  void UpdateAllSphProperties(int, int, SphParticle<ndim> *, 
                              Sph<ndim> *, Nbody<ndim> *);
  void UpdateAllSphForces(int, int, SphParticle<ndim> *, 
                          Sph<ndim> *, Nbody<ndim> *);
  void UpdateAllSphHydroForces(int, int, SphParticle<ndim> *, 
                               Sph<ndim> *, Nbody<ndim> *);
  void UpdateAllSphGravForces(int, int, SphParticle<ndim> *, 
                              Sph<ndim> *, Nbody<ndim> *);

};



//=============================================================================
//  Class SM2012SphBruteForce
/// \brief   ..
/// \details Class for computing SPH neighbour lists using brute force only 
///          (i.e. direct summation over all particles).
/// \author  D. A. Hubber, G. Rosotti
/// \date    12/05/2014
//=============================================================================
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
  

  //---------------------------------------------------------------------------
  SM2012SphBruteForce(FLOAT, DomainBox<ndim> *, 
                     SphKernel<ndim> *, CodeTiming *);
  ~SM2012SphBruteForce();


  //---------------------------------------------------------------------------
  void UpdateAllSphProperties(int, int, SphParticle<ndim> *, 
                              Sph<ndim> *, Nbody<ndim> *);
  void UpdateAllSphForces(int, int, SphParticle<ndim> *, 
                          Sph<ndim> *, Nbody<ndim> *);
  void UpdateAllSphHydroForces(int, int, SphParticle<ndim> *, 
                               Sph<ndim> *, Nbody<ndim> *);
  void UpdateAllSphGravForces(int, int, SphParticle<ndim> *, 
                              Sph<ndim> *, Nbody<ndim> *);

};



//=============================================================================
//  Class GodunovSphBruteForce
/// \brief   ..
/// \details Class for computing SPH neighbour lists using brute force only 
///          (i.e. direct summation over all particles).
/// \author  D. A. Hubber, G. Rosotti
/// \date    12/05/2014
//=============================================================================
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
  

  //---------------------------------------------------------------------------
  GodunovSphBruteForce(FLOAT, DomainBox<ndim> *, 
                       SphKernel<ndim> *, CodeTiming *);
  ~GodunovSphBruteForce();


  //---------------------------------------------------------------------------
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
  using SphNeighbourSearch<ndim>::kernrangesqd;


  //---------------------------------------------------------------------------
  SphTree(int, FLOAT, FLOAT, FLOAT, string, string,
          DomainBox<ndim> *, SphKernel<ndim> *, CodeTiming *);
  ~SphTree();


  //---------------------------------------------------------------------------
  void BuildTree(bool, int, int, int, int, int, 
                 SphParticle<ndim> *, Sph<ndim> *, FLOAT);
  int GetGatherNeighbourList(FLOAT *, FLOAT, SphParticle<ndim> *, 
                             int, int, int *);
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
#if defined(VERIFY_ALL)
  void CheckValidNeighbourList(int, int, int, int *, 
			       ParticleType<ndim> *, string);
#endif


  // Additional functions for binary tree neighbour search
  //---------------------------------------------------------------------------
  void AllocateMemory(Sph<ndim> *);
  void DeallocateMemory(void);


  // Additional variables for grid
  //---------------------------------------------------------------------------
  const int Nleafmax;               ///< Max. number of particles per leaf cell
  const FLOAT thetamaxsqd;          ///< Geometric opening angle squared
  const FLOAT invthetamaxsqd;       ///< 1 / thetamaxsqd
  const FLOAT macerror;             ///< Error tolerance for gravity tree-MAC
  const string gravity_mac;         ///< Multipole-acceptance criteria for tree
  const string multipole;           ///< Multipole-order for cell gravity

  // Class variables
  //---------------------------------------------------------------------------
  int Ncell;                        ///< Current no. of grid cells
  int Ncellmax;                     ///< Max. allowed no. of grid cells
  int Nlistmax;                     ///< Max. length of neighbour list
  int Ntot;                         ///< No. of current points in list
  int Ntotold;                      ///< Prev. no. of particles
  int Ntotmax;                      ///< Max. no. of points in list
  int Ntotmaxold;                   ///< Old value of Ntotmax
  FLOAT hmax;                       ///< Store hmax in the tree
  FLOAT theta;                      ///< Geometric opening angle
  KDTree<ndim,ParticleType> *tree;  ///< Pointer to tree
  KDTree<ndim,ParticleType> *ghosttree;  ///< Pointer to tree containing ghosts

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



//=============================================================================
//  Class GradhSphTree
/// \brief   Class containing kd-tree for computing grad-h SPH force loops.
/// \details kd-tree data structure used for efficient neighbour searching 
///          and computation of gravitational forces for grad-h SPH.
/// \author  D. A. Hubber
/// \date    08/01/2014
//=============================================================================
template <int ndim, template<int> class ParticleType>
class GradhSphTree: public SphTree<ndim,ParticleType>
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
  using SphTree<ndim,ParticleType>::neibpartbuf;
  using SphTree<ndim,ParticleType>::Ndirectmaxbuf;
  using SphTree<ndim,ParticleType>::Ngravcellmaxbuf;
  using SphTree<ndim,ParticleType>::Nleafmax;
  using SphTree<ndim,ParticleType>::Nneibmaxbuf;
  using SphTree<ndim,ParticleType>::Ntot;
  using SphTree<ndim,ParticleType>::Ntotmax;
  using SphTree<ndim,ParticleType>::Ntotmaxold;
  using SphTree<ndim,ParticleType>::Ntotold;
  using SphTree<ndim,ParticleType>::timing;
  using SphTree<ndim,ParticleType>::tree;


  //---------------------------------------------------------------------------
  GradhSphTree(int, FLOAT, FLOAT, FLOAT, string, string, 
               DomainBox<ndim> *, SphKernel<ndim> *, CodeTiming *);
  ~GradhSphTree();


  //---------------------------------------------------------------------------
  void UpdateAllSphProperties(int, int, SphParticle<ndim> *, 
                 Sph<ndim> *, Nbody<ndim> *);
  void UpdateAllSphForces(int, int, SphParticle<ndim> *, 
                          Sph<ndim> *, Nbody<ndim> *);
  void UpdateAllSphHydroForces(int, int, SphParticle<ndim> *, 
                               Sph<ndim> *, Nbody<ndim> *);
  void UpdateAllSphGravForces(int, int, SphParticle<ndim> *, 
                              Sph<ndim> *, Nbody<ndim> *);

};



//=============================================================================
//  Class SM2012SphTree
/// \brief   Class containing kd-tree for computing grad-h SPH force loops.
/// \details kd-tree data structure used for efficient neighbour searching 
///          and computation of gravitational forces for grad-h SPH.
/// \author  D. A. Hubber
/// \date    08/01/2014
//=============================================================================
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
  using SphTree<ndim,ParticleType>::neibpartbuf;
  using SphTree<ndim,ParticleType>::Ndirectmaxbuf;
  using SphTree<ndim,ParticleType>::Ngravcellmaxbuf;
  using SphTree<ndim,ParticleType>::Nleafmax;
  using SphTree<ndim,ParticleType>::Nneibmaxbuf;
  using SphTree<ndim,ParticleType>::Ntot;
  using SphTree<ndim,ParticleType>::Ntotmax;
  using SphTree<ndim,ParticleType>::Ntotmaxold;
  using SphTree<ndim,ParticleType>::Ntotold;
  using SphTree<ndim,ParticleType>::timing;
  using SphTree<ndim,ParticleType>::tree;


  //---------------------------------------------------------------------------
  SM2012SphTree(int, FLOAT, FLOAT, FLOAT, string, string, 
                DomainBox<ndim> *, SphKernel<ndim> *, CodeTiming *);
  ~SM2012SphTree();


  //---------------------------------------------------------------------------
  void UpdateAllSphProperties(int, int, SphParticle<ndim> *, 
                 Sph<ndim> *, Nbody<ndim> *);
  void UpdateAllSphForces(int, int, SphParticle<ndim> *, 
                          Sph<ndim> *, Nbody<ndim> *);
  void UpdateAllSphHydroForces(int, int, SphParticle<ndim> *, 
                               Sph<ndim> *, Nbody<ndim> *);
  void UpdateAllSphGravForces(int, int, SphParticle<ndim> *, 
                              Sph<ndim> *, Nbody<ndim> *);

};



//=============================================================================
//  Class GodunovSphTree
/// \brief   Class containing kd-tree for computing grad-h SPH force loops.
/// \details kd-tree data structure used for efficient neighbour searching 
///          and computation of gravitational forces for grad-h SPH.
/// \author  D. A. Hubber
/// \date    08/01/2014
//=============================================================================
template <int ndim, template<int> class ParticleType>
class GodunovSphTree: public SphTree<ndim,ParticleType>
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
  using SphTree<ndim,ParticleType>::neibpartbuf;
  using SphTree<ndim,ParticleType>::Ndirectmaxbuf;
  using SphTree<ndim,ParticleType>::Ngravcellmaxbuf;
  using SphTree<ndim,ParticleType>::Nleafmax;
  using SphTree<ndim,ParticleType>::Nneibmaxbuf;
  using SphTree<ndim,ParticleType>::Ntot;
  using SphTree<ndim,ParticleType>::Ntotmax;
  using SphTree<ndim,ParticleType>::Ntotmaxold;
  using SphTree<ndim,ParticleType>::Ntotold;
  using SphTree<ndim,ParticleType>::timing;
  using SphTree<ndim,ParticleType>::tree;


  //---------------------------------------------------------------------------
  GodunovSphTree(int, FLOAT, FLOAT, FLOAT, string, string, 
               DomainBox<ndim> *, SphKernel<ndim> *, CodeTiming *);
  ~GodunovSphTree();


  //---------------------------------------------------------------------------
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

};
#endif
