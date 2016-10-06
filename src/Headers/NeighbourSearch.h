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
  virtual double GetMaximumSmoothingLength() = 0 ;
#ifdef MPI_PARALLEL
  virtual void BuildPrunedTree(const int, const int, const DomainBox<ndim> &,
                               const MpiNode<ndim> *, Particle<ndim> *) = 0;
  virtual void BuildMpiGhostTree(const bool, const int, const int, const int, const int, const int,
                                 const FLOAT, Particle<ndim> *, Hydrodynamics<ndim> *) = 0;
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
  virtual void UpdateGravityExportList(int, int, int, Particle<ndim> *, Hydrodynamics<ndim> *,
                                       Nbody<ndim> *, const DomainBox<ndim> &) = 0;
  virtual void UpdateHydroExportList(int, int, int, Particle<ndim> *, Hydrodynamics<ndim> *,
                                     Nbody<ndim> *, const DomainBox<ndim> &) = 0;
  virtual void UnpackReturnedExportInfo(vector<char >& received_information,
                                        Hydrodynamics<ndim>* hydro,
                                        const int rank, const int iproc) = 0;
  virtual void FindParticlesToTransfer(Hydrodynamics<ndim> *, vector<vector<int> >& ,
                                       vector<int> &, const vector<int> &, MpiNode<ndim> *) = 0;
  virtual void ResetCountersExportInfo(Hydrodynamics<ndim>* hydro) = 0;
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
//  Class HydroTree
/// \brief   Class containing tree for efficient neighbour searching and gravity calculations.
/// \details Class containing tree for efficient neighbour searching and gravity calculations.
/// \author  D. A. Hubber
/// \date    08/01/2014
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
class HydroTree : public virtual NeighbourSearch<ndim>
{
	void ReallocateMemory(void);
#if defined MPI_PARALLEL
  vector<vector<int> > ids_sent_particles;
  vector<vector<int> > ids_sent_cells;
  vector<int> N_imported_part_per_proc;
  vector<int> N_imported_cells_per_proc;
protected:
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
  virtual double GetMaximumSmoothingLength() ;
#ifdef MPI_PARALLEL
  virtual void BuildPrunedTree(const int, const int, const DomainBox<ndim> &,
                               const MpiNode<ndim> *, Particle<ndim> *);
  virtual void BuildMpiGhostTree(const bool, const int, const int, const int, const int, const int,
                                 const FLOAT, Particle<ndim> *, Hydrodynamics<ndim> *);
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

    const int size_header = 2*sizeof(int);

    int size = size_particles + size_cells + size_header;

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
  virtual void UpdateGravityExportList(int, int, int, Particle<ndim> *, Hydrodynamics<ndim> *,
                                       Nbody<ndim> *, const DomainBox<ndim> &);
  virtual void UpdateHydroExportList(int, int, int, Particle<ndim> *, Hydrodynamics<ndim> *,
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
    tree->InitialiseCellWorkCounters() ;
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
  MultipoleMoment<ndim> **cellbuf;                 ///< Buffers of tree-cell copies

  Tree<ndim,ParticleType,TreeCell> *tree;          ///< Pointer to main (local) tree
  Tree<ndim,ParticleType,TreeCell> *ghosttree;     ///< Pointer to tree containing ghosts
                                                   ///< on local domain

#ifdef MPI_PARALLEL
  int Nprunedcellmax;                              ///< Max. number of cells in pruned tree
  int *Npartexport;                                ///< No. of ptcls to be exported (per MPI node)

  vector<vector<int>> cellexportlist;              ///< List of cell ids
  Tree<ndim,ParticleType,TreeCell> *mpighosttree;  ///< Pointer to tree containing
                                                   ///< ghosts from other MPI procs.
  Tree<ndim,ParticleType,TreeCell> **prunedtree;   ///< 'Pruned' tree for MPI nodes.
  Tree<ndim,ParticleType,TreeCell> **sendprunedtree;  ///< 'Pruned' tree for MPI nodes.
#endif

};


//=================================================================================================
// Gravity force functions
///
/// These functions compute the force on particles due to distant tree cells using one of a
/// variety of approximations:
///   Monopole or Quadropole forces from each cell summed per particle
///   Fast Monopole, per cell, Taylor expanded to the location of each particle.
//=================================================================================================

//=================================================================================================
//  ComputeCellMonopoleForces
/// Compute the force on particle 'parti' due to all cells obtained in the
/// gravity tree walk.  Uses only monopole moments (i.e. COM) of the cell.
//=================================================================================================
template<int ndim>
void ComputeCellMonopoleForces(FLOAT &gpot, ///< [inout] Grav. potential
FLOAT agrav[ndim],                   ///< [inout] Acceleration array
FLOAT rp[ndim],                      ///< [in] Position of point
int Ngravcell,                       ///< [in] No. of tree cells in list
MultipoleMoment<ndim> *gravcell)     ///< [in] List of tree cell ids
{
  FLOAT dr[ndim];                      // Relative position vector

  // Loop over all neighbouring particles in list
  //-----------------------------------------------------------------------------------------------
  for (int cc=0; cc<Ngravcell; cc++) {
    MultipoleMoment<ndim>& cell = gravcell[cc];

    FLOAT mc = cell.m;
    for (int k=0; k<ndim; k++) dr[k] = cell.r[k] - rp[k];
    FLOAT drsqd    = DotProduct(dr,dr,ndim) + small_number;
    FLOAT invdrsqd = (FLOAT) 1.0/drsqd;
    FLOAT invdrmag = sqrt(invdrsqd);
    FLOAT invdr3   = invdrsqd*invdrmag;

    gpot += mc*invdrmag;
    for (int k=0; k<ndim; k++) agrav[k] += mc*dr[k]*invdr3;

  }
  //-----------------------------------------------------------------------------------------------

  return;
}


//=================================================================================================
//  ComputeQuadrupole
/// Compute the quadropole part of the force, depending on dimension.
//=================================================================================================
inline void ComputeQuadropole(const MultipoleMoment<1>& cell, const FLOAT rp[1],
                              FLOAT agrav[1], FLOAT& gpot) {

  FLOAT dr[1] ;
  for (int k=0; k<1; k++) dr[k] = cell.r[k] - rp[k];
  FLOAT drsqd    = DotProduct(dr,dr,1) + small_number;
  FLOAT invdrsqd = (FLOAT) 1.0/drsqd;
  FLOAT invdrmag = sqrt(invdrsqd);
  FLOAT invdr5 = invdrsqd*invdrsqd*invdrmag ;

  // First add monopole term for acceleration
  for (int k=0; k<1; k++) agrav[k] += cell.m*dr[k]*invdrsqd*invdrmag;

  FLOAT qscalar = cell.q[0]*dr[0]*dr[0];
  FLOAT qfactor = 2.5*qscalar*invdr5*invdrsqd;

  agrav[0] += (cell.q[0]*dr[0])*invdr5 - qfactor*dr[0];
  gpot += cell.m*invdrmag + 0.5*qscalar*invdr5;
}
inline void ComputeQuadropole(const MultipoleMoment<2>& cell, const FLOAT rp[2],
                              FLOAT agrav[2], FLOAT& gpot) {
  FLOAT dr[2] ;
  for (int k=0; k<2; k++) dr[k] = cell.r[k] - rp[k];
  FLOAT drsqd    = DotProduct(dr,dr,2) + small_number;
  FLOAT invdrsqd = (FLOAT) 1.0/drsqd;
  FLOAT invdrmag = sqrt(invdrsqd);
  FLOAT invdr5 = invdrsqd*invdrsqd*invdrmag ;

  // First add monopole term for acceleration
  for (int k=0; k<2; k++) agrav[k] += cell.m*dr[k]*invdrsqd*invdrmag;

  FLOAT qscalar = cell.q[0]*dr[0]*dr[0] + cell.q[2]*dr[1]*dr[1] +
    2.0*cell.q[1]*dr[0]*dr[1];
  FLOAT qfactor = 2.5*qscalar*invdr5*invdrsqd;

  agrav[0] += (cell.q[0]*dr[0] + cell.q[1]*dr[1])*invdr5 - qfactor*dr[0];
  agrav[1] += (cell.q[1]*dr[0] + cell.q[2]*dr[1])*invdr5 - qfactor*dr[1];
  gpot += cell.m*invdrmag + 0.5*qscalar*invdr5;
}
inline void ComputeQuadropole(const MultipoleMoment<3>& cell, const FLOAT rp[3],
                              FLOAT agrav[3], FLOAT& gpot) {
  FLOAT dr[3] ;
  for (int k=0; k<3; k++) dr[k] = cell.r[k] - rp[k];
  FLOAT drsqd    = DotProduct(dr,dr,3) + small_number;
  FLOAT invdrsqd = (FLOAT) 1.0/drsqd;
  FLOAT invdrmag = sqrt(invdrsqd);
  FLOAT invdr5 = invdrsqd*invdrsqd*invdrmag ;

  // First add monopole term for acceleration
  for (int k=0; k<3; k++) agrav[k] += cell.m*dr[k]*invdrsqd*invdrmag;

  FLOAT qscalar = cell.q[0]*dr[0]*dr[0] + cell.q[2]*dr[1]*dr[1] -
          (cell.q[0] + cell.q[2])*dr[2]*dr[2] +
           2.0*(cell.q[1]*dr[0]*dr[1] + cell.q[3]*dr[0]*dr[2] + cell.q[4]*dr[1]*dr[2]);

  FLOAT qfactor = 2.5*qscalar*invdr5*invdrsqd;

  agrav[0] +=
      (cell.q[0]*dr[0] + cell.q[1]*dr[1] + cell.q[3]*dr[2])*invdr5 - qfactor*dr[0];
  agrav[1] +=
      (cell.q[1]*dr[0] + cell.q[2]*dr[1] + cell.q[4]*dr[2])*invdr5 - qfactor*dr[1];
  agrav[2] +=
      (cell.q[3]*dr[0] + cell.q[4]*dr[1] - (cell.q[0] + cell.q[2])*dr[2])*invdr5 - qfactor*dr[2];

  gpot += cell.m*invdrmag + 0.5*qscalar*invdr5;
}


//=================================================================================================
//  ComputeCellQuadrupoleForces
/// Compute the force on particle 'parti' due to all cells obtained in the
/// gravity tree walk including the quadrupole moment correction term.
//=================================================================================================
template <int ndim>
void ComputeCellQuadrupoleForces
 (FLOAT &gpot,                         ///< [inout] Grav. potential
  FLOAT agrav[ndim],                   ///< [inout] Acceleration array
  FLOAT rp[ndim],                      ///< [in] Position of point
  int Ngravcell,                       ///< [in] No. of tree cells in list
  MultipoleMoment<ndim> *gravcell)     ///< [in] List of tree cell ids
{

  // Loop over all neighbouring particles in list
  //-----------------------------------------------------------------------------------------------
  for (int cc=0; cc<Ngravcell; cc++) {
    ComputeQuadropole(gravcell[cc], rp, agrav, gpot) ;
  }
  //-----------------------------------------------------------------------------------------------


  return;
}

//=================================================================================================
//  class FastMultipoleForces
/// \brief Class for computing the gravitational forces using the fast monopole expansion
//=================================================================================================
template<int ndim>
class FastMultipoleForces {
public:

  FastMultipoleForces(const FLOAT r[ndim])
  {
    for (int k=0; k<ndim; k++) rc[k]   = r[k];
    for (int k=0; k<ndim; k++) ac[k]   = 0;
    for (int k=0; k<ndim; k++) dphi[k] = 0;
    for (int k=0; k<(ndim*(ndim+1))/2; k++) q[k] = 0;
    pot = 0 ;
  }

  inline void AddMonopoleContribution(const MultipoleMoment<ndim>& cell) ;
  inline void AddQuadrupoleContribution(const MultipoleMoment<ndim>& cell) ;

  inline void ApplyForcesTaylor(const FLOAT r[ndim], FLOAT agrav[ndim], FLOAT& gpot) const ;

private:
  FLOAT rc[ndim] ;
  FLOAT ac[ndim]  ;
  FLOAT dphi[ndim] ;
  FLOAT q[(ndim*(ndim+1))/2];
  FLOAT pot ;
};

//=================================================================================================
//  AddCellContribution
/// Compute the fast monopole cell-cell terms depending on dimension.
//=================================================================================================
template<>
inline void FastMultipoleForces<1>::AddMonopoleContribution(const MultipoleMoment<1>& cell) {
  FLOAT dr[1] ;
  for (int k=0; k<1; k++) dr[k] = cell.r[k] - rc[k];
  FLOAT invdrmag = sqrt((FLOAT) 1.0/DotProduct(dr,dr,1));
  FLOAT invdrsqd = invdrmag*invdrmag;
  FLOAT invdr3   = invdrsqd*invdrmag;

  FLOAT mc = cell.m;
  pot += mc*invdrmag;

  mc *= invdr3;
  for (int k=0; k<1; k++) ac[k] += mc*dr[k];
  for (int k=0; k<1; k++) dphi[k] += mc*dr[k];
  q[0] += mc*(3.0*dr[0]*dr[0]*invdrsqd - 1);
}
template<>
inline void FastMultipoleForces<2>::AddMonopoleContribution(const MultipoleMoment<2>& cell) {
  FLOAT dr[2] ;
  for (int k=0; k<2; k++) dr[k] = cell.r[k] - rc[k];
  FLOAT invdrmag = sqrt((FLOAT) 1.0/DotProduct(dr,dr,2));
  FLOAT invdrsqd = invdrmag*invdrmag;
  FLOAT invdr3   = invdrsqd*invdrmag;

  FLOAT mc = cell.m;
  pot += mc*invdrmag;

  mc *= invdr3;
  for (int k=0; k<2; k++) ac[k] += mc*dr[k];
  for (int k=0; k<2; k++) dphi[k] += mc*dr[k];
  q[0] += mc*(3.0*dr[0]*dr[0]*invdrsqd - 1);
  q[1] += mc*(3.0*dr[0]*dr[1]*invdrsqd);
  q[2] += mc*(3.0*dr[1]*dr[1]*invdrsqd - 1);
}
template<>
inline void FastMultipoleForces<3>::AddMonopoleContribution(const MultipoleMoment<3>& cell)  {
  FLOAT dr[3] ;
  for (int k=0; k<3; k++) dr[k] = cell.r[k] - rc[k];
  FLOAT invdrmag = sqrt((FLOAT) 1.0/DotProduct(dr,dr,3));
  FLOAT invdrsqd = invdrmag*invdrmag;
  FLOAT invdr3   = invdrsqd*invdrmag;

  FLOAT mc = cell.m;
  pot += mc*invdrmag;

  mc *= invdr3;
  for (int k=0; k<3; k++) ac[k] += mc*dr[k];
  for (int k=0; k<3; k++) dphi[k] += mc*dr[k];
  q[0] += mc*(3.0*dr[0]*dr[0]*invdrsqd - 1);
  q[1] += mc*(3.0*dr[0]*dr[1]*invdrsqd);
  q[2] += mc*(3.0*dr[1]*dr[1]*invdrsqd - 1);
  q[3] += mc*(3.0*dr[2]*dr[0]*invdrsqd);
  q[4] += mc*(3.0*dr[2]*dr[1]*invdrsqd);
  q[5] += mc*(3.0*dr[2]*dr[2]*invdrsqd - 1);
}

template<>
inline void FastMultipoleForces<1>::AddQuadrupoleContribution(const MultipoleMoment<1>& cell) {
  AddMonopoleContribution(cell);

  FLOAT dr[1] ;
  for (int k=0; k<1; k++) dr[k] = cell.r[k] - rc[k];
  FLOAT drsqd    = DotProduct(dr,dr,1) + small_number;
  FLOAT invdrsqd = (FLOAT) 1.0/drsqd;
  FLOAT invdrmag = sqrt(invdrsqd);
  FLOAT invdr5 = invdrsqd*invdrsqd*invdrmag ;

  FLOAT qscalar = cell.q[0]*dr[0]*dr[0] ;
  FLOAT qfactor = 2.5*qscalar*invdr5*invdrsqd;

  FLOAT qx[1];
  qx[0] = (cell.q[0]*dr[0])*invdr5;

  pot += 0.5*qscalar*invdr5;

  for (int k=0; k<1; k++) ac[k]   += qx[k] - qfactor*dr[k];
  for (int k=0; k<1; k++) dphi[k] += qx[k] - qfactor*dr[k];
  for (int k=0; k<1; k++) qx[k] *= 5*invdrsqd;

  q[0] += - qfactor*(7*dr[0]*dr[0]*invdrsqd - 1);

  q[0] += qx[0]*dr[0] + qx[0]*dr[0] - cell.q[0]*invdr5;
}
template<>
inline void FastMultipoleForces<2>::AddQuadrupoleContribution(const MultipoleMoment<2>& cell) {
  AddMonopoleContribution(cell);

  FLOAT dr[2] ;
  for (int k=0; k<2; k++) dr[k] = cell.r[k] - rc[k];
  FLOAT drsqd    = DotProduct(dr,dr,2) + small_number;
  FLOAT invdrsqd = (FLOAT) 1.0/drsqd;
  FLOAT invdrmag = sqrt(invdrsqd);
  FLOAT invdr5 = invdrsqd*invdrsqd*invdrmag ;

  FLOAT qscalar = cell.q[0]*dr[0]*dr[0] + cell.q[2]*dr[1]*dr[1] +
      2.0*cell.q[1]*dr[0]*dr[1];

  FLOAT qfactor = 2.5*qscalar*invdr5*invdrsqd;

  FLOAT qx[2];
  qx[0] = (cell.q[0]*dr[0] + cell.q[1]*dr[1])*invdr5;
  qx[1] = (cell.q[1]*dr[0] + cell.q[2]*dr[1])*invdr5;

  pot += 0.5*qscalar*invdr5;

  for (int k=0; k<2; k++) ac[k]   += qx[k] - qfactor*dr[k];
  for (int k=0; k<2; k++) dphi[k] += qx[k] - qfactor*dr[k];
  for (int k=0; k<2; k++) qx[k] *= 5*invdrsqd;

  q[0] += - qfactor*(7*dr[0]*dr[0]*invdrsqd - 1);
  q[1] += - qfactor*(7*dr[0]*dr[1]*invdrsqd);
  q[2] += - qfactor*(7*dr[1]*dr[1]*invdrsqd - 1);

  q[0] += qx[0]*dr[0] + qx[0]*dr[0] - cell.q[0]*invdr5;
  q[1] += qx[0]*dr[1] + qx[1]*dr[0] - cell.q[1]*invdr5;
  q[2] += qx[1]*dr[1] + qx[1]*dr[1] - cell.q[2]*invdr5;
}
template<>
inline void FastMultipoleForces<3>::AddQuadrupoleContribution(const MultipoleMoment<3>& cell) {
  AddMonopoleContribution(cell);

  FLOAT dr[3] ;
  for (int k=0; k<3; k++) dr[k] = cell.r[k] - rc[k];
  FLOAT drsqd    = DotProduct(dr,dr,3) + small_number;
  FLOAT invdrsqd = (FLOAT) 1.0/drsqd;
  FLOAT invdrmag = sqrt(invdrsqd);
  FLOAT invdr5 = invdrsqd*invdrsqd*invdrmag ;

  FLOAT qscalar = cell.q[0]*dr[0]*dr[0] + cell.q[2]*dr[1]*dr[1] -
      (cell.q[0] + cell.q[2])*dr[2]*dr[2] +
      2.0*(cell.q[1]*dr[0]*dr[1] + cell.q[3]*dr[0]*dr[2] + cell.q[4]*dr[1]*dr[2]);

  FLOAT qfactor = 2.5*qscalar*invdr5*invdrsqd;

  FLOAT qx[3];
  qx[0] = (cell.q[0]*dr[0] + cell.q[1]*dr[1] + cell.q[3]*dr[2])*invdr5;
  qx[1] = (cell.q[1]*dr[0] + cell.q[2]*dr[1] + cell.q[4]*dr[2])*invdr5;
  qx[2] = (cell.q[3]*dr[0] + cell.q[4]*dr[1] -(cell.q[0]+cell.q[2])*dr[2])*invdr5;

  pot += 0.5*qscalar*invdr5;

  for (int k=0; k<3; k++) ac[k]   += qx[k] - qfactor*dr[k];
  for (int k=0; k<3; k++) dphi[k] += qx[k] - qfactor*dr[k];
  for (int k=0; k<3; k++) qx[k] *= 5*invdrsqd;

  q[0] += - qfactor*(7*dr[0]*dr[0]*invdrsqd - 1);
  q[1] += - qfactor*(7*dr[0]*dr[1]*invdrsqd);
  q[2] += - qfactor*(7*dr[1]*dr[1]*invdrsqd - 1);
  q[3] += - qfactor*(7*dr[0]*dr[2]*invdrsqd);
  q[4] += - qfactor*(7*dr[1]*dr[2]*invdrsqd);
  q[5] += - qfactor*(7*dr[2]*dr[2]*invdrsqd - 1);

  q[0] += qx[0]*dr[0] + qx[0]*dr[0] - cell.q[0]*invdr5;
  q[1] += qx[0]*dr[1] + qx[1]*dr[0] - cell.q[1]*invdr5;
  q[2] += qx[1]*dr[1] + qx[1]*dr[1] - cell.q[2]*invdr5;
  q[3] += qx[0]*dr[2] + qx[2]*dr[0] - cell.q[3]*invdr5;
  q[4] += qx[1]*dr[2] + qx[2]*dr[1] - cell.q[4]*invdr5;
  q[5] += qx[2]*dr[2] + qx[2]*dr[2] + (cell.q[0] + cell.q[2])*invdr5;
}
//=================================================================================================
//  ApplyMonopoleForces
/// Apply the fast monopole forces at a given location.
//=================================================================================================
template<>
inline void FastMultipoleForces<1>::ApplyForcesTaylor(const FLOAT r[1], FLOAT agrav[1], FLOAT& gpot) const {
  FLOAT dr[1];
  for (int k=0; k<1; k++) dr[k] = r[k] - rc[k];

  agrav[0] += ac[0] + q[0]*dr[0];
  gpot += pot + dphi[0]*dr[0];
}
template<>
inline void FastMultipoleForces<2>::ApplyForcesTaylor(const FLOAT r[2], FLOAT agrav[2], FLOAT& gpot) const {
  FLOAT dr[2];
  for (int k=0; k<2; k++) dr[k] = r[k] - rc[k];

  agrav[0] += ac[0] + q[0]*dr[0] + q[1]*dr[1];
  agrav[1] += ac[1] + q[1]*dr[0] + q[2]*dr[1];
  gpot += pot + dphi[0]*dr[0] + dphi[1]*dr[1];
}
template<>
inline void FastMultipoleForces<3>::ApplyForcesTaylor(const FLOAT r[3], FLOAT agrav[3], FLOAT& gpot) const {
  FLOAT dr[3];
  for (int k=0; k<3; k++) dr[k] = r[k] - rc[k];

  agrav[0] += ac[0] + q[0]*dr[0] + q[1]*dr[1] + q[3]*dr[2];
  agrav[1] += ac[1] + q[1]*dr[0] + q[2]*dr[1] + q[4]*dr[2];
  agrav[2] += ac[2] + q[3]*dr[0] + q[4]*dr[1] + q[5]*dr[2];
  gpot += pot + dphi[0]*dr[0] + dphi[1]*dr[1] + dphi[2]*dr[2];
}


//=================================================================================================
//  ComputeFastMonopoleForces
/// Compute the force on particle 'parti' due to all cells obtained in the
/// gravity tree walk.  Uses only monopole moments (i.e. COM) of the cell.
//=================================================================================================
template <int ndim, template<int> class ParticleType>
void ComputeFastMonopoleForces
 (int Nactive,                         ///< [in] No. of active particles
  int Ngravcell,                       ///< [in] No. of tree cells in list
  MultipoleMoment<ndim> *gravcell,     ///< [in] List of tree cell ids
  TreeCellBase<ndim> &cell,            ///< [in] Current cell pointer
  ParticleType<ndim> *activepart)      ///< [inout] Active Hydrodynamics particle array
{

  FastMultipoleForces<ndim> monopole(cell.r) ;

  //-----------------------------------------------------------------------------------------------
  for (int cc=0; cc<Ngravcell; cc++) {
#ifndef MPI_PARALLEL
    assert(cell.id != gravcell[cc].id);
#endif
    monopole.AddMonopoleContribution(gravcell[cc]);
  }

  for (int j=0; j<Nactive; j++)
    monopole.ApplyForcesTaylor(activepart[j].r, activepart[j].agrav, activepart[j].gpot) ;
  //-----------------------------------------------------------------------------------------------

  return;
}

//=================================================================================================
//  ComputeFastQuadrupoleForces
/// Compute the force on particle 'parti' due to all cells obtained in the
/// gravity tree walk.  Uses quadrupole moments (i.e. COM) of the cell.
//=================================================================================================
template <int ndim, template<int> class ParticleType>
void ComputeFastQuadrupoleForces
 (int Nactive,                         ///< [in] No. of active particles
  int Ngravcell,                       ///< [in] No. of tree cells in list
  MultipoleMoment<ndim> *gravcell,     ///< [in] List of tree cell ids
  TreeCellBase<ndim> &cell,            ///< [in] Current cell pointer
  ParticleType<ndim> *activepart)      ///< [inout] Active Hydrodynamics particle array
{

  FastMultipoleForces<ndim> monopole(cell.r) ;

  //-----------------------------------------------------------------------------------------------
  for (int cc=0; cc<Ngravcell; cc++) {
#ifndef MPI_PARALLEL
    assert(cell.id != gravcell[cc].id);
#endif
    monopole.AddQuadrupoleContribution(gravcell[cc]);
  }

  for (int j=0; j<Nactive; j++)
    monopole.ApplyForcesTaylor(activepart[j].r, activepart[j].agrav, activepart[j].gpot) ;
  //-----------------------------------------------------------------------------------------------

  return;
}




//=================================================================================================
// Tree constructor factory templates
///
/// These are simple functions to construct the correct tree based upon the given arguments.
///
/// A thin proxy struct is used to ensure that the template specialisation can be done cleanly and
/// correctly since partial template specialisation is needed. The struct template below should be
/// specialized as required.
//=================================================================================================

//=================================================================================================
// struct __construct_tree_impl
// The implementation structs for the tree constructor factory functions.
//=================================================================================================
template<int ndim, template<int> class ParticleType, template<int> class TreeCell>
struct __construct_tree_impl
{
	typedef Tree<ndim, ParticleType, TreeCell> return_type ;

	static return_type*  construct(int Nleafmaxaux, FLOAT thetamaxsqdaux,
								   FLOAT kernrangeaux, FLOAT macerroraux,
								   string gravity_mac_aux, string multipole_aux,
								   const DomainBox<ndim>& domain,
								   const ParticleTypeRegister& reg)
	{
	  string message = "Tree cell type for GradhSphTree not recognised." ;
	  ExceptionHandler::getIstance().raise(message);
	  return NULL ;
	}
};


//=================================================================================================
//  new_tree
/// Construct a single tree based on the template types.
/// This function is used to construct at KD, Oct or BruteForce as required from the TreeCell type
/// DO NOT EVER SPECIALIZE THIS TEMPLATE. BAD THINGS ARE GUARANTEED TO HAPPEN. DON'T BLAME ME.
///   USE THE __construct_tree_impl STRUCTURE INSTEAD.
//=================================================================================================
template<int ndim, template<int> class ParticleType, template<int> class TreeCell>
Tree<ndim, ParticleType, TreeCell>* new_tree(int Nleafmax, FLOAT thetamaxsqd,
		   	   	   	   	 	 	 	 	 	 FLOAT kernrange, FLOAT macerror,
		   	   	   	   	 	 	 	 	 	 string gravity_mac, string multipole,
		   	   	   	   	 	 	 	 	 	 const DomainBox<ndim>& domain,
		   	   	   	   	 	 	 	 	 	 const ParticleTypeRegister& reg)
{
  return __construct_tree_impl<ndim, ParticleType, TreeCell>::construct
		  (Nleafmax, thetamaxsqd, kernrange, macerror,
		  gravity_mac,  multipole,domain, reg) ;
}
//=================================================================================================
//  new_tree_array
/// Construct an array of pointer to tree based on the template types.
/// This function is used to construct at KD, Oct or BruteForce as required from the TreeCell type.
/// DO NOT EVER SPECIALIZE THIS TEMPLATE. BAD THINGS ARE GUARANTEED TO HAPPEN. DON'T BLAME ME.
///   USE THE __construct_tree_impl STRUCTURE INSTEAD.
//=================================================================================================
template<int ndim, template<int> class ParticleType, template<int> class TreeCell>
Tree<ndim, ParticleType, TreeCell>** new_tree_array(int NumTrees)
{
  typename __construct_tree_impl<ndim, ParticleType, TreeCell>::return_type** derived =
      __construct_tree_impl<ndim, ParticleType, TreeCell>::construct_array(NumTrees) ;

  return reinterpret_cast<Tree<ndim,ParticleType,TreeCell>**>(derived) ;
}


//=================================================================================================
// struct __construct_tree_impl
// KDTree specialisation
//=================================================================================================
template<int ndim, template<int> class ParticleType>
struct  __construct_tree_impl<ndim, ParticleType, KDTreeCell> {
  typedef KDTree<ndim,ParticleType, KDTreeCell> return_type ;

  static return_type*  construct(int Nleafmax, FLOAT thetamaxsqd,
		                         FLOAT kernrange, FLOAT macerror,
	 	 	 	 	 	         string gravity_mac, string multipole,
	 	 	 	 	 	         const DomainBox<ndim>& domain,
	 	 	 	 	 	         const ParticleTypeRegister& reg)
  {
	return new return_type(Nleafmax, thetamaxsqd, kernrange, macerror,
			 	 	 	   gravity_mac,  multipole,domain, reg) ;
  }
} ;
//=================================================================================================
// struct __construct_tree_impl
// OctTree specialisation
//=================================================================================================
template<int ndim, template<int> class ParticleType>
struct  __construct_tree_impl<ndim, ParticleType, OctTreeCell> {
  typedef OctTree<ndim,ParticleType, OctTreeCell> return_type ;

  static return_type*  construct(int Nleafmax, FLOAT thetamaxsqd,
		                         FLOAT kernrange, FLOAT macerror,
	 	 	 	 	 	         string gravity_mac, string multipole,
	 	 	 	 	 	         const DomainBox<ndim>& domain,
	 	 	 	 	 	         const ParticleTypeRegister& reg)
  {
	return new return_type(Nleafmax, thetamaxsqd, kernrange, macerror,
			 	 	 	   gravity_mac,  multipole,domain, reg) ;
  }
} ;
//=================================================================================================
// struct __construct_tree_impl
// OctTree specialisation with TreeRay Cells
//=================================================================================================
template<int ndim, template<int> class ParticleType>
struct  __construct_tree_impl<ndim, ParticleType, TreeRayCell> {
  typedef OctTree<ndim,ParticleType, TreeRayCell> return_type ;

  static return_type*  construct(int Nleafmax, FLOAT thetamaxsqd,
		                         FLOAT kernrange, FLOAT macerror,
	 	 	 	 	 	         string gravity_mac, string multipole,
	 	 	 	 	 	         const DomainBox<ndim>& domain,
	 	 	 	 	 	         const ParticleTypeRegister& reg)
  {
	return new return_type(Nleafmax, thetamaxsqd, kernrange, macerror,
			 	 	 	   gravity_mac,  multipole,domain, reg) ;
  }
} ;
//=================================================================================================
// struct __construct_tree_impl
// Brute Force Tree specialisation
//=================================================================================================
template<int ndim, template<int> class ParticleType>
struct  __construct_tree_impl<ndim, ParticleType, BruteForceTreeCell> {
  typedef BruteForceTree<ndim,ParticleType, BruteForceTreeCell> return_type ;

  static return_type*  construct(int Nleafmax, FLOAT thetamaxsqd,
		                         FLOAT kernrange, FLOAT macerror,
	 	 	 	 	 	         string gravity_mac, string multipole,
	 	 	 	 	 	         const DomainBox<ndim>& domain,
	 	 	 	 	 	         const ParticleTypeRegister& reg)
  {
	return new return_type(Nleafmax, thetamaxsqd, kernrange, macerror,
			 	 	 	   gravity_mac,  multipole,domain, reg) ;
  }
} ;



#endif
