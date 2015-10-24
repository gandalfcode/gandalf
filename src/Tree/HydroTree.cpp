//=================================================================================================
//  HydroTree.cpp
//  Contains all functions for managing the tree for hydrodynamical particles.
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


#include <cstdlib>
#include <cassert>
#include <iostream>
#include <numeric>
#include <string>
#include <math.h>
#include <vector>
#include "Precision.h"
#include "Exception.h"
#include "Hydrodynamics.h"
#include "NeighbourSearch.h"
#include "Parameters.h"
#include "InlineFuncs.h"
#include "Particle.h"
#include "Debug.h"
#if defined _OPENMP
#include <omp.h>
#endif
using namespace std;




//=================================================================================================
//  HydroTree::HydroTree
/// HydroTree constructor.  Initialises various variables.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
HydroTree<ndim,ParticleType,TreeCell>::HydroTree
 (int _Nleafmax, int _Nmpi, int _pruning_level_min, int _pruning_level_max, FLOAT _thetamaxsqd,
  FLOAT _kernrange, FLOAT _macerror, string _gravity_mac, string _multipole,
  DomainBox<ndim>* _box, SmoothingKernel<ndim>* _kern, CodeTiming* _timing):
  NeighbourSearch<ndim>(_kernrange, _box, _kern, _timing),
  Nleafmax(_Nleafmax),
  Nmpi(_Nmpi),
  pruning_level_min(_pruning_level_min),
  pruning_level_max(_pruning_level_max),
  thetamaxsqd(_thetamaxsqd),
  invthetamaxsqd((FLOAT) 1.0/_thetamaxsqd),
  macerror(_macerror),
  gravity_mac(_gravity_mac),
  multipole(_multipole)
{
  allocated_buffer = false;
  neibcheck        = true;
  Ntot             = 0;
  Ntotmax          = 0;
  Ntotmaxold       = 0;
#if defined _OPENMP
  Nthreads         = omp_get_max_threads();
#else
  Nthreads         = 1;
#endif
#ifdef MPI_PARALLEL
  Ncellexport = new int[Nmpi];
  Npartexport = new int[Nmpi];
  cellexportlist = new TreeCell<ndim>**[Nmpi];
  for (int j=0; j<Nmpi; j++) cellexportlist[j] = NULL;
  ids_sent_particles.resize(Nmpi);
#endif
}



//=================================================================================================
//  HydroTree::~HydroTree
/// HydroTree destructor.  Deallocates tree memory upon object destruction.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
HydroTree<ndim,ParticleType,TreeCell>::~HydroTree()
{
  if (tree->allocated_tree) {
    DeallocateMemory();
    tree->DeallocateTreeMemory();
  }
}



//=================================================================================================
//  HydroTree::AllocateMemory
/// Allocate memory for tree as requested.  If more memory is required
/// than currently allocated, tree is deallocated and reallocated here.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void HydroTree<ndim,ParticleType,TreeCell>::AllocateMemory
 (const int Ngather)                   ///< [in] Average no. of gather neighbours
{
  debug2("[HydroTree::AllocateMemory]");

  if (!allocated_buffer) {

    Nneibmaxbuf     = new int[Nthreads];
    Ngravcellmaxbuf = new int[Nthreads];
    activelistbuf   = new int*[Nthreads];
    levelneibbuf    = new int*[Nthreads];
    activepartbuf   = new ParticleType<ndim>*[Nthreads];
    neibpartbuf     = new ParticleType<ndim>*[Nthreads];
    cellbuf         = new TreeCell<ndim>*[Nthreads];

    for (int ithread=0; ithread<Nthreads; ithread++) {
      Nneibmaxbuf[ithread]     = max(1, 4*Ngather);
      Ngravcellmaxbuf[ithread] = max(1, 4*Ngather);
      activelistbuf[ithread]   = new int[Nleafmax];
      levelneibbuf[ithread]    = new int[Ntotmax];
      activepartbuf[ithread]   = new ParticleType<ndim>[Nleafmax];
      neibpartbuf[ithread]     = new ParticleType<ndim>[Nneibmaxbuf[ithread]];
      cellbuf[ithread]         = new TreeCell<ndim>[Ngravcellmaxbuf[ithread]];
    }
    allocated_buffer = true;

  }

  return;
}



//=================================================================================================
//  HydroTree::DeallocateTreeMemory
/// Deallocates all binary tree memory
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void HydroTree<ndim,ParticleType,TreeCell>::DeallocateMemory(void)
{
  int ithread;                         // Thread id number

  debug2("[HydroTree::DeallocateTreeMemory]");

  if (allocated_buffer) {

    for (ithread=0; ithread<Nthreads; ithread++) {
      delete[] cellbuf[ithread];
      delete[] levelneibbuf[ithread];
      delete[] neibpartbuf[ithread];
      delete[] activepartbuf[ithread];
      delete[] activelistbuf[ithread];
    }
    delete[] cellbuf;
    delete[] levelneibbuf;
    delete[] neibpartbuf;
    delete[] activepartbuf;
    delete[] activelistbuf;
    delete[] Ngravcellmaxbuf;
    delete[] Nneibmaxbuf;

  }

  return;
}



//=================================================================================================
//  HydroTree::BuildTree
/// Main routine to control how the tree is built, re-stocked and interpolated during each step.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void HydroTree<ndim,ParticleType,TreeCell>::BuildTree
 (const bool rebuild_tree,             ///< [in] Flag to rebuild tree
  const int n,                         ///< [in] Integer time
  const int ntreebuildstep,            ///< [in] Tree build frequency
  const int ntreestockstep,            ///< [in] Tree stocking frequency
  const int Npart,                     ///< [in] No. of particles
  const int Npartmax,                  ///< [in] Max. no. of particles
  const FLOAT timestep,                ///< [in] Smallest physical timestep
  Particle<ndim> *part_gen,            ///< [inout] Particle data array
  Hydrodynamics<ndim> *hydro)          ///< [inout] Pointer to Hydrodynamics object
{
  ParticleType<ndim> *partdata = static_cast<ParticleType<ndim>* > (part_gen);

  debug2("[HydroTree::BuildTree]");
  timing->StartTimingSection("BUILD_TREE");

  // Activate nested parallelism for tree building routines
#ifdef _OPENMP
  omp_set_nested(1);
#endif


  // For tree rebuild steps
  //-----------------------------------------------------------------------------------------------
  if (n%ntreebuildstep == 0 || rebuild_tree) {

    // Delete any dead particles from main Hydrodynamics arrays before we re-build tree
    hydro->DeleteDeadParticles();

    Ntotold    = Ntot;
    Ntot       = hydro->Ntot;
    Ntotmaxold = Ntotmax;
    Ntotmax    = max(Ntotmax, Ntot);
    Ntotmax    = max(Ntotmax, hydro->Nhydromax);
    assert(Ntotmax >= Ntot);

    tree->Ntot       = hydro->Nhydro;
    tree->Ntotmaxold = tree->Ntotmax;
    tree->Ntotmax    = max(tree->Ntotmax, tree->Ntot);
    tree->Ntotmax    = max(tree->Ntotmax, hydro->Nhydromax);
    tree->BuildTree(0, hydro->Nhydro-1, Npart, Npartmax, timestep, partdata);

    AllocateMemory(hydro->Ngather);
#ifdef MPI_PARALLEL
    if (Ntotmax > Ntotmaxold) {
      for (int i=Nmpi-1; i>=0; i--) delete[] cellexportlist[i];
      for (int i=0; i<Nmpi; i++) cellexportlist[i] = new TreeCell<ndim>*[tree->gmax];
      assert(tree->gmax > 0);
    }
#endif

  }

  // Else stock the tree
  //-----------------------------------------------------------------------------------------------
  else if (n%ntreestockstep == 0) {

    tree->StockTree(tree->celldata[0],partdata);

  }

  // Otherwise simply extrapolate tree cell properties
  //-----------------------------------------------------------------------------------------------
  else {

    tree->ExtrapolateCellProperties(timestep);

  }
  //-----------------------------------------------------------------------------------------------

#ifdef _OPENMP
  omp_set_nested(0);
#endif

  timing->EndTimingSection("BUILD_TREE");

  return;
}



//=================================================================================================
//  HydroTree::BuildGhostTree
/// Main routine to control how the tree is built, re-stocked and interpolated
/// during each timestep.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void HydroTree<ndim,ParticleType,TreeCell>::BuildGhostTree
 (const bool rebuild_tree,             ///< [in] Flag to rebuild tree
  const int n,                         ///< [in] Integer time
  const int ntreebuildstep,            ///< [in] Tree build frequency
  const int ntreestockstep,            ///< [in] Tree stocking frequency
  const int Npart,                     ///< [in] No. of particles
  const int Npartmax,                  ///< [in] Max. no. of particles
  const FLOAT timestep,                ///< [in] Smallest physical timestep
  Particle<ndim> *part_gen,            ///< [inout] Particle data array
  Hydrodynamics<ndim> *hydro)          ///< [inout] Pointer to Hydrodynamics object
{
  ParticleType<ndim> *partdata = static_cast<ParticleType<ndim>* > (part_gen);

  // If no periodic ghosts exist, do not build tree
  if (hydro->NPeriodicGhost == 0) return;

  debug2("[HydroTree::BuildGhostTree]");
  timing->StartTimingSection("BUILD_GHOST_TREE");

  // Activate nested parallelism for tree building routines
#ifdef _OPENMP
  omp_set_nested(1);
#endif


  // For tree rebuild steps
  //-----------------------------------------------------------------------------------------------
  if (n%ntreebuildstep == 0 || rebuild_tree) {

    ghosttree->Ntot       = hydro->NPeriodicGhost;
    ghosttree->Ntotmaxold = ghosttree->Ntotmax;
    ghosttree->Ntotmax    = max(ghosttree->Ntotmax, ghosttree->Ntot);
    ghosttree->Ntotmax    = max(ghosttree->Ntotmax, hydro->Nhydromax);
    ghosttree->BuildTree(hydro->Nhydro, hydro->Nhydro + hydro->NPeriodicGhost - 1,
                         ghosttree->Ntot, ghosttree->Ntotmax, timestep, partdata);

  }

  // Else stock the tree
  //-----------------------------------------------------------------------------------------------
  else if (n%ntreestockstep == 0) {

    ghosttree->StockTree(ghosttree->celldata[0], partdata);

  }

  // Otherwise simply extrapolate tree cell properties
  //-----------------------------------------------------------------------------------------------
  else {

    ghosttree->ExtrapolateCellProperties(timestep);

  }
  //-----------------------------------------------------------------------------------------------

#ifdef _OPENMP
  omp_set_nested(0);
#endif

  timing->EndTimingSection("BUILD_GHOST_TREE");

  return;
}



//=================================================================================================
//  HydroTree::GetGatherNeighbourList
/// Compute the gather neighbour list at the point 'rp' of all particles within a search radius
/// of 'rsearch'.  Searches through real and ghost neighbour trees.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
int HydroTree<ndim,ParticleType,TreeCell>::GetGatherNeighbourList
 (FLOAT rp[ndim],                      ///< [in] Position vector
  FLOAT rsearch,                       ///< [in] Gather search radius
  Particle<ndim> *part_gen,            ///< [in] Pointer to Hydrodynamics particle array
  int Nhydro,                          ///< [in] No. of hydro particles
  int Nneibmax,                        ///< [in] Max. no. of neighbours
  int *neiblist)                       ///< [out] List of neighbouring particles
{
  int Nneib = 0;                       // No. of (non-dead) neighbours
  ParticleType<ndim>* partdata = static_cast<ParticleType<ndim>* > (part_gen);

  debug2("[HydroTree::GetGatherNeighbourList]");

  Nneib = tree->ComputeGatherNeighbourList(partdata, rp, rsearch, Nneibmax, Nneib, neiblist);
  Nneib = ghosttree->ComputeGatherNeighbourList(partdata, rp, rsearch, Nneibmax, Nneib, neiblist);
#ifdef MPI_PARALLEL
  Nneib = mpighosttree->ComputeGatherNeighbourList(partdata, rp, rsearch, Nneibmax, Nneib, neiblist);
#endif

  return Nneib;
}



//=================================================================================================
//  HydroTree::UpdateActiveParticleCounters
/// Loop through all leaf cells in the tree and update all active particle counters.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void HydroTree<ndim,ParticleType,TreeCell>::UpdateActiveParticleCounters
 (Particle<ndim> *partdata_gen,        ///< [inout] Pointer to hydrodynamics particles array
  Hydrodynamics<ndim> *hydro)          ///< [in] Pointer to hydrodynamics object
{
  ParticleType<ndim>* partdata = static_cast<ParticleType<ndim>* > (partdata_gen);
  tree->UpdateActiveParticleCounters(partdata);
}



//=================================================================================================
//  HydroTree::SearchBoundaryGhostParticles
/// Search domain to create any required ghost particles near any boundaries.
/// Currently only searches to create periodic or mirror ghost particles.
//=================================================================================================
template <int ndim, template <int> class ParticleType, template<int> class TreeCell>
void HydroTree<ndim,ParticleType,TreeCell>::SearchBoundaryGhostParticles
 (FLOAT tghost,                                ///< [in] Ghost particle 'lifetime'
  DomainBox<ndim> &simbox,                     ///< [in] Simulation box structure
  Hydrodynamics<ndim> *hydro)                  ///< [inout] Hydrodynamics object pointer
{
  int c;                                       // Cell counter
  int i;                                       // Particle counter
  TreeCell<ndim> *cellptr;                     // Pointer to tree cell
  const FLOAT grange = ghost_range*kernrange;  // Range of ghost particles (in terms of h)

  // Set all relevant particle counters
  hydro->Nghost         = 0;
  hydro->NPeriodicGhost = 0;
  hydro->Nmpighost      = 0;
  hydro->Nghostmax      = hydro->Nhydromax - hydro->Nhydro;
  hydro->Ntot           = hydro->Nhydro;


  // If all boundaries are open, immediately return to main loop
  if (simbox.boundary_lhs[0] == openBoundary && simbox.boundary_rhs[0] == openBoundary &&
      simbox.boundary_lhs[1] == openBoundary && simbox.boundary_rhs[1] == openBoundary &&
      simbox.boundary_lhs[2] == openBoundary && simbox.boundary_rhs[2] == openBoundary) return;


  debug2("[HydroTree::SearchBoundaryGhostParticles]");


  // Create ghost particles in x-dimension
  //===============================================================================================
  if ((simbox.boundary_lhs[0] == openBoundary &&
       simbox.boundary_rhs[0] == openBoundary) == 0) {

    // Start from root-cell
    c = 0;

    //---------------------------------------------------------------------------------------------
    while (c < tree->Ncell) {
      cellptr = &(tree->celldata[c]);

      // If x-bounding box overlaps edge of x-domain, open cell
      //-------------------------------------------------------------------------------------------
      if (cellptr->bbmin[0] + min((FLOAT) 0.0,cellptr->v[0]*tghost) <
          simbox.boxmin[0] + grange*cellptr->hmax ||
          cellptr->bbmax[0] + max((FLOAT) 0.0,cellptr->v[0]*tghost) >
          simbox.boxmax[0] - grange*cellptr->hmax) {

        // If not a leaf-cell, then open cell to first child cell
        if (cellptr->level != tree->ltot)
          c++;

        else if (cellptr->N == 0)
          c = cellptr->cnext;

        // If leaf-cell, check through particles in turn to find ghosts
        else if (cellptr->level == tree->ltot) {
          i = cellptr->ifirst;
          while (i != -1) {
            hydro->CheckXBoundaryGhostParticle(i,tghost,simbox);
            if (i == cellptr->ilast) break;
            i = tree->inext[i];
          };
          c = cellptr->cnext;
        }
      }

      // If not in range, then open next cell
      //-------------------------------------------------------------------------------------------
      else
        c = cellptr->cnext;

    }
    //---------------------------------------------------------------------------------------------

    hydro->Ntot = hydro->Nhydro + hydro->Nghost;
  }


  // Create ghost particles in y-dimension
  //===============================================================================================
  if (ndim >= 2 && (simbox.boundary_lhs[1] == openBoundary &&
                    simbox.boundary_rhs[1] == openBoundary) == 0) {

    // Start from root-cell
    c = 0;

    //---------------------------------------------------------------------------------------------
    while (c < tree->Ncell) {
      cellptr = &(tree->celldata[c]);

      // If x-bounding box overlaps edge of x-domain, open cell
      //-------------------------------------------------------------------------------------------
      if (cellptr->bbmin[1] + min((FLOAT) 0.0,cellptr->v[1]*tghost) <
          simbox.boxmin[1] + grange*cellptr->hmax ||
          cellptr->bbmax[1] + max((FLOAT) 0.0,cellptr->v[1]*tghost) >
          simbox.boxmax[1] - grange*cellptr->hmax) {

        // If not a leaf-cell, then open cell to first child cell
        if (cellptr->level != tree->ltot)
          c++;

        else if (cellptr->N == 0)
          c = cellptr->cnext;

        // If leaf-cell, check through particles in turn to find ghosts
        else if (cellptr->level == tree->ltot) {
          i = cellptr->ifirst;
          while (i != -1) {
            hydro->CheckYBoundaryGhostParticle(i,tghost,simbox);
            if (i == cellptr->ilast) break;
            i = tree->inext[i];
          };
          c = cellptr->cnext;
        }
      }

      // If not in range, then open next cell
      //-------------------------------------------------------------------------------------------
      else
        c = cellptr->cnext;

    }
    //---------------------------------------------------------------------------------------------


    // Check x-ghosts (which are not part of tree) by direct-sum
    for (i=hydro->Nhydro; i<hydro->Ntot; i++) hydro->CheckYBoundaryGhostParticle(i,tghost,simbox);

    hydro->Ntot = hydro->Nhydro + hydro->Nghost;
  }


  // Create ghost particles in z-dimension
  //===============================================================================================
  if (ndim == 3 && (simbox.boundary_lhs[2] == openBoundary &&
      simbox.boundary_rhs[2] == openBoundary) == 0) {

    // Start from root-cell
    c = 0;

    //---------------------------------------------------------------------------------------------
    while (c < tree->Ncell) {
      cellptr = &(tree->celldata[c]);

      // If x-bounding box overlaps edge of x-domain, open cell
      //-------------------------------------------------------------------------------------------
      if (cellptr->bbmin[2] + min((FLOAT) 0.0,cellptr->v[2]*tghost) <
          simbox.boxmin[2] + grange*cellptr->hmax ||
          cellptr->bbmax[2] + max((FLOAT) 0.0,cellptr->v[2]*tghost) >
          simbox.boxmax[2] - grange*cellptr->hmax) {

        // If not a leaf-cell, then open cell to first child cell
        if (cellptr->level != tree->ltot)
          c++;

        else if (cellptr->N == 0)
          c = cellptr->cnext;

        // If leaf-cell, check through particles in turn to find ghosts
        else if (cellptr->level == tree->ltot) {
          i = cellptr->ifirst;
          while (i != -1) {
            hydro->CheckZBoundaryGhostParticle(i,tghost,simbox);
            if (i == cellptr->ilast) break;
            i = tree->inext[i];
          };
          c = cellptr->cnext;
        }
      }

      // If not in range, then open next cell
      //-------------------------------------------------------------------------------------------
      else
        c = cellptr->cnext;

    }
    //---------------------------------------------------------------------------------------------


    // Check x- and y-ghosts (which are not part of tree) by direct-sum
    for (i=hydro->Nhydro; i<hydro->Ntot; i++) hydro->CheckZBoundaryGhostParticle(i,tghost,simbox);

    hydro->Ntot = hydro->Nhydro + hydro->Nghost;
  }


  // Quit here if we've run out of memory for ghosts
  if (hydro->Ntot > hydro->Nhydromax) {
    string message="Not enough memory for ghost particles";
    ExceptionHandler::getIstance().raise(message);
  }

  hydro->NPeriodicGhost = hydro->Nghost;

  return;
}



//=================================================================================================
//  HydroTree::ComputeCellMonopoleForces
/// Compute the force on particle 'parti' due to all cells obtained in the
/// gravity tree walk.  Uses only monopole moments (i.e. COM) of the cell.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void HydroTree<ndim,ParticleType,TreeCell>::ComputeCellMonopoleForces
 (FLOAT &gpot,                         ///< [inout] Grav. potential
  FLOAT agrav[ndim],                   ///< [inout] Acceleration array
  FLOAT rp[ndim],                      ///< [in] Position of point
  int Ngravcell,                       ///< [in] No. of tree cells in list
  TreeCell<ndim> *gravcell)            ///< [in] List of tree cell ids
{
  int cc;                              // Aux. cell counter
  int k;                               // Dimension counter
  FLOAT dr[ndim];                      // Relative position vector
  FLOAT drsqd;                         // Distance squared
  FLOAT invdrmag;                      // 1 / distance
  FLOAT invdrsqd;                      // 1 / drsqd
  FLOAT invdr3;                        // 1 / dist^3
  FLOAT mc;                            // Mass of cell
  TreeCell<ndim> *cellptr;             // Pointer to gravity tree cell

  // Loop over all neighbouring particles in list
  //-----------------------------------------------------------------------------------------------
  for (cc=0; cc<Ngravcell; cc++) {
    cellptr = &(gravcell[cc]);

    mc = cellptr->m;
    for (k=0; k<ndim; k++) dr[k] = cellptr->r[k] - rp[k];
    drsqd    = DotProduct(dr,dr,ndim) + small_number;
    invdrsqd = (FLOAT) 1.0/drsqd;
    invdrmag = sqrt(invdrsqd);
    invdr3   = invdrsqd*invdrmag;

    gpot += mc*invdrmag;
    for (k=0; k<ndim; k++) agrav[k] += mc*dr[k]*invdr3;

  }
  //-----------------------------------------------------------------------------------------------

  return;
}



//=================================================================================================
//  HydroTree::ComputeCellQuadrupoleForces
/// Compute the force on particle 'parti' due to all cells obtained in the
/// gravity tree walk including the quadrupole moment correction term.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void HydroTree<ndim,ParticleType,TreeCell>::ComputeCellQuadrupoleForces
 (FLOAT &gpot,                         ///< [inout] Grav. potential
  FLOAT agrav[ndim],                   ///< [inout] Acceleration array
  FLOAT rp[ndim],                      ///< [in] Position of point
  int Ngravcell,                       ///< [in] No. of tree cells in list
  TreeCell<ndim> *gravcell)            ///< [in] List of tree cell ids
{
  int cc;                              // Aux. cell counter
  int k;                               // Dimension counter
  FLOAT dr[ndim];                      // Relative position vector
  FLOAT drsqd;                         // Distance squared
  FLOAT invdrsqd;                      // 1 / drsqd
  FLOAT invdrmag;                      // 1 / distance
  FLOAT invdr5;                        // 1 / distance^5
  FLOAT qfactor;                       // Constant factor for optimisation
  FLOAT qscalar;                       // Quadrupole moment scalar quantity
  TreeCell<ndim> *cellptr;             // Pointer to gravity tree cell


  // Loop over all neighbouring particles in list
  //-----------------------------------------------------------------------------------------------
  for (cc=0; cc<Ngravcell; cc++) {
    cellptr = &(gravcell[cc]);

    for (k=0; k<ndim; k++) dr[k] = cellptr->r[k] - rp[k];
    drsqd    = DotProduct(dr,dr,ndim) + small_number;
    invdrsqd = (FLOAT) 1.0/drsqd;
    invdrmag = sqrt(invdrsqd);
    invdr5   = invdrsqd*invdrsqd*invdrmag;

    // First add monopole term for acceleration
    for (k=0; k<ndim; k++) agrav[k] += cellptr->m*dr[k]*invdrsqd*invdrmag;

    // Now add quadrupole moment terms depending on dimensionality
    if (ndim == 3) {
      qscalar = cellptr->q[0]*dr[0]*dr[0] + cellptr->q[2]*dr[1]*dr[1] -
        (cellptr->q[0] + cellptr->q[2])*dr[2]*dr[2] +
         2.0*(cellptr->q[1]*dr[0]*dr[1] + cellptr->q[3]*dr[0]*dr[2] + cellptr->q[4]*dr[1]*dr[2]);
      qfactor = 2.5*qscalar*invdr5*invdrsqd;
      agrav[0] +=
        (cellptr->q[0]*dr[0] + cellptr->q[1]*dr[1] + cellptr->q[3]*dr[2])*invdr5 - qfactor*dr[0];
      agrav[1] +=
        (cellptr->q[1]*dr[0] + cellptr->q[2]*dr[1] + cellptr->q[4]*dr[2])*invdr5 - qfactor*dr[1];
      agrav[2] += (cellptr->q[3]*dr[0] + cellptr->q[4]*dr[1] -
        (cellptr->q[0] + cellptr->q[2])*dr[2])*invdr5 - qfactor*dr[2];
      gpot += cellptr->m*invdrmag + 0.5*qscalar*invdr5;
    }
    else if (ndim == 2) {
      qscalar = cellptr->q[0]*dr[0]*dr[0] + cellptr->q[2]*dr[1]*dr[1] +
        2.0*cellptr->q[1]*dr[0]*dr[1];
      qfactor = 2.5*qscalar*invdr5*invdrsqd;
      agrav[0] += (cellptr->q[0]*dr[0] + cellptr->q[1]*dr[1])*invdr5 - qfactor*dr[0];
      agrav[1] += (cellptr->q[1]*dr[0] + cellptr->q[2]*dr[1])*invdr5 - qfactor*dr[1];
      gpot += cellptr->m*invdrmag + 0.5*qscalar*invdr5;
    }

  }
  //-----------------------------------------------------------------------------------------------


  return;
}



//=================================================================================================
//  HydroTree::ComputeFastMonopoleForces
/// Compute the force on particle 'parti' due to all cells obtained in the
/// gravity tree walk.  Uses only monopole moments (i.e. COM) of the cell.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void HydroTree<ndim,ParticleType,TreeCell>::ComputeFastMonopoleForces
 (int Nactive,                         ///< [in] No. of active particles
  int Ngravcell,                       ///< [in] No. of tree cells in list
  TreeCell<ndim> *gravcell,            ///< [in] List of tree cell ids
  TreeCell<ndim> &cell,                ///< [in] Current cell pointer
  ParticleType<ndim> *activepart)      ///< [inout] Active Hydrodynamics particle array
{
  int cc;                              // Aux. cell counter
  int j;                               // Counter over active particles in cell
  int k;                               // Dimension counter
  FLOAT ac[ndim];                      // Acceleration at cell centre
  FLOAT dr[ndim];                      // Relative position vector
  FLOAT drsqd;                         // Distance squared
  FLOAT invdrmag;                      // 1 / distance
  FLOAT invdrsqd;                      // 1 / drsqd
  FLOAT invdr3;                        // 1 / dist^3
  FLOAT mc;                            // Mass of cell
  FLOAT q[6];                          // Local copy of quadrupole moment
  FLOAT dphi[3];                       // Potential gradient (same as accel?)
  FLOAT cellpot;                       // Potential at cell centre
  FLOAT rc[ndim];                      // Position of cell centre

  for (k=0; k<ndim; k++) rc[k] = cell.r[k];
  for (k=0; k<ndim; k++) ac[k] = (FLOAT) 0.0;
  for (k=0; k<ndim; k++) dphi[k] = (FLOAT) 0.0;
  for (k=0; k<6; k++) q[k] = (FLOAT) 0.0;
  cellpot = (FLOAT) 0.0;


  //-----------------------------------------------------------------------------------------------
  if (ndim == 3) {

    for (cc=0; cc<Ngravcell; cc++) {
#ifndef MPI_PARALLEL
      assert(cell.id != gravcell[cc].id);
#endif
      mc = gravcell[cc].m;
      for (k=0; k<ndim; k++) dr[k] = gravcell[cc].r[k] - rc[k];
      drsqd    = DotProduct(dr,dr,ndim);
      invdrsqd = (FLOAT) 1.0/drsqd;
      invdrmag = sqrt(invdrsqd);
      invdr3   = invdrsqd*invdrmag;
      cellpot += mc*invdrmag;
      for (k=0; k<ndim; k++) ac[k] += mc*dr[k]*invdr3;
      for (k=0; k<ndim; k++) dphi[k] += mc*dr[k]*invdr3;
      q[0] += mc*(3.0*dr[0]*dr[0]*invdr3*invdrsqd - invdrsqd*invdrmag);
      q[1] += mc*(3.0*dr[0]*dr[1]*invdr3*invdrsqd);
      q[2] += mc*(3.0*dr[1]*dr[1]*invdr3*invdrsqd - invdrsqd*invdrmag);
      q[3] += mc*(3.0*dr[2]*dr[0]*invdr3*invdrsqd);
      q[4] += mc*(3.0*dr[2]*dr[1]*invdr3*invdrsqd);
      q[5] += mc*(3.0*dr[2]*dr[2]*invdr3*invdrsqd - invdrsqd*invdrmag);
    }

    for (j=0; j<Nactive; j++) {
      for (k=0; k<ndim; k++) dr[k] = activepart[j].r[k] - rc[k];
      activepart[j].agrav[0] += ac[0] + q[0]*dr[0] + q[1]*dr[1] + q[3]*dr[2];
      activepart[j].agrav[1] += ac[1] + q[1]*dr[0] + q[2]*dr[1] + q[4]*dr[2];
      activepart[j].agrav[2] += ac[2] + q[3]*dr[0] + q[4]*dr[1] + q[5]*dr[2];
      activepart[j].gpot += cellpot + dphi[0]*dr[0] + dphi[1]*dr[1] + dphi[2]*dr[2];
    }

  }
  //-----------------------------------------------------------------------------------------------

  return;
}



//=================================================================================================
//  HydroTree::UpdateAllStarGasForces
/// Calculate the gravitational acceleration on all star particles due to
/// all gas particles via the tree.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void HydroTree<ndim,ParticleType,TreeCell>::UpdateAllStarGasForces
 (int Nhydro,                          ///< [in] No. of SPH particles
  int Ntot,                            ///< [in] No. of SPH + ghost particles
  Particle<ndim> *part_gen,            ///< [inout] Pointer to SPH ptcl array
  Hydrodynamics<ndim> *hydro,          ///< [in] Pointer to SPH object
  Nbody<ndim> *nbody)                  ///< [in] Pointer to N-body object
{
  int Nactive;                         // No. of active particles in cell
  int *activelist;                     // List of active particle ids
  NbodyParticle<ndim> *star;           // Pointer to star particle
  ParticleType<ndim>* partdata = static_cast<ParticleType<ndim>* > (part_gen);


  debug2("[GradhSphTree::UpdateAllStarGasForces]");
  //timing->StartTimingSection("STAR_GAS_GRAV_FORCES");

  // Make list of all active stars
  Nactive = 0;
  activelist = new int[nbody->Nstar];
  for (int i=0; i<nbody->Nstar; i++) {
    if (nbody->nbodydata[i]->active) activelist[Nactive++] = i;
  }


  // Set-up all OMP threads
  //===============================================================================================
#pragma omp parallel default(none) private(star)\
  shared(activelist,hydro,Nactive,Ntot,nbody,partdata,cout)
  {
#if defined _OPENMP
    const int ithread = omp_get_thread_num();
#else
    const int ithread = 0;
#endif
    int i;                                       // Particle id
    int j;                                       // Aux. particle counter
    int okflag;                                  // Flag if h-rho iteration is valid
    int Ndirect;                                 // No. of direct-sum gravity particles
    int Ngravcell;                               // No. of gravity cells
    int Nneib;                                   // No. of neighbours
    int Nneibmax = Ntot; //Nneibmaxbuf[ithread];
    int Ngravcellmax = Ngravcellmaxbuf[ithread]; // ..
    FLOAT macfactor;                             // Gravity MAC factor
    int* neiblist = new int[Nneibmax];           // ..
    int* directlist = new int[Nneibmax];         // ..
    TreeCell<ndim>* gravcell = new TreeCell<ndim>[Ngravcellmax];   // ..


    // Loop over all active cells
    //=============================================================================================
#pragma omp for schedule(dynamic)
    for (j=0; j<Nactive; j++) {
      i = activelist[j];
      star = nbody->nbodydata[i];

      // Compute average/maximum term for computing gravity MAC
      if (gravity_mac == "eigenmac") macfactor = pow((FLOAT) 1.0/star->gpot,twothirds);
      else macfactor = (FLOAT) 0.0;

      // Compute neighbour list for cell depending on physics options
      okflag = tree->ComputeStarGravityInteractionList
       (star, macfactor, Nneibmax, Nneibmax, Ngravcellmax, Nneib,
        Ndirect, Ngravcell, neiblist, directlist, gravcell, partdata);

      // If there are too many neighbours, reallocate the arrays and recompute the neighbour lists.
      while (okflag == -1) {
        delete[] gravcell;
        Ngravcellmax = 2*Ngravcellmax;
        gravcell = new TreeCell<ndim>[Ngravcellmax];
        okflag = tree->ComputeStarGravityInteractionList
         (star, macfactor, Nneibmax, Nneibmax, Ngravcellmax, Nneib,
          Ndirect, Ngravcell, neiblist, directlist, gravcell, partdata);
      };

      // Compute contributions to star force from nearby SPH particles
      nbody->CalculateDirectHydroForces(star, Nneib, Ndirect, neiblist, directlist, hydro);

      // Compute gravitational force due to distant cells
      if (multipole == "monopole" || multipole == "fast_monopole") {
        this->ComputeCellMonopoleForces(star->gpot, star->a, star->r, Ngravcell, gravcell);
      }
      else if (multipole == "quadrupole") {
        this->ComputeCellQuadrupoleForces(star->gpot, star->a, star->r, Ngravcell, gravcell);
      }


    }
    //=============================================================================================


    // Free-up local memory for OpenMP thread
    delete[] gravcell;
    delete[] directlist;
    delete[] neiblist;

  }
  //===============================================================================================

  delete[] activelist;

  //timing->EndTimingSection("STAR_GAS_GRAV_FORCES");

  return;
}



#ifdef MPI_PARALLEL
//=================================================================================================
//  HydroTree::UpdateGravityExportList
/// Compute gravitational forces due to other MPI node trees (using pruned trees).
/// If the other domains are too close (so the pruned trees are not adequate), then flag
/// cell to be exported to that MPI node for a full local tree-walk.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void HydroTree<ndim,ParticleType,TreeCell>::UpdateGravityExportList
 (int rank,                            ///< [in] MPI rank
  int Nhydro,                          ///< [in] No. of hydro particles
  int Ntot,                            ///< [in] No. of hydro + ghost particles
  Particle<ndim> *part_gen,            ///< [inout] Pointer to Hydrodynamics ptcl array
  Hydrodynamics<ndim> *hydro,          ///< [in] Pointer to Hydrodynamics object
  Nbody<ndim> *nbody,                  ///< [in] Pointer to N-body object
  const DomainBox<ndim> &simbox)       ///< [in] Simulation domain box
{
  int cactive;                         // No. of active cells
  TreeCell<ndim> **celllist;           // List of pointers to binary tree cells
  ParticleType<ndim>* partdata = static_cast<ParticleType<ndim>* > (part_gen);

  debug2("[GradhHydroTree::UpdateGravityExportForces]");
  timing->StartTimingSection("HYDRO_DISTANT_FORCES");


  // Find list of all cells that contain active particles
  celllist = new TreeCell<ndim>*[2*tree->gtot];
  cactive = tree->ComputeActiveCellPointers(celllist);

  // Reset all export lists
  for (int j=0; j<Nmpi; j++) {
    Ncellexport[j] = 0;
    Npartexport[j] = 0;
  }


  // Set-up all OMP threads
  //===============================================================================================
#pragma omp parallel default(none) shared(celllist,cactive,cout,hydro,partdata,rand)
  {
#if defined _OPENMP
    const int ithread = omp_get_thread_num();
#else
    const int ithread = 0;
#endif
    int cc;                                    // Aux. cell counter
    int i;                                     // Particle id
    int index_cell;                            // ..
    int j;                                     // Aux. particle counter
    int k;                                     // Dimension counter
    int Nactive;                               // No. of active particles in current cell
    int Ngravcell=0;                           // No. of gravity cells
    int Ngravcellmax = Nprunedcellmax;         // Max. size of gravity cell pointer array
    int Ngravcelltemp;                         // Aux. gravity cell counter
    FLOAT macfactor;                           // Gravity MAC factor for cell
    int *activelist                = activelistbuf[ithread];
    ParticleType<ndim> *activepart = activepartbuf[ithread];
    TreeCell<ndim> *gravcelllist   = new TreeCell<ndim>[Ngravcellmax];


    // Loop over all active cells
    //=============================================================================================
#pragma omp for schedule(guided)
    for (cc=0; cc<cactive; cc++) {
      TreeCell<ndim>* cellptr = celllist[cc];
      TreeCell<ndim>& cell = *cellptr;
      macfactor = (FLOAT) 0.0;
      Ngravcell = 0;

      // Find list of active particles in current cell
      Nactive = tree->ComputeActiveParticleList(cell, partdata, activelist);

      // Make local copies of active particles
      for (j=0; j<Nactive; j++) activepart[j] = partdata[activelist[j]];

      // Compute average/maximum term for computing gravity MAC
      if (gravity_mac == "eigenmac") {
        for (j=0; j<Nactive; j++) {
          macfactor = max(macfactor, pow((FLOAT) 1.0/activepart[j].gpot, twothirds));
        }
      }

      // Zero/initialise all summation variables for active particles
      for (j=0; j<Nactive; j++) {
        activepart[j].gpot = (FLOAT) 0.0;
        for (k=0; k<ndim; k++) activepart[j].a[k]     = (FLOAT) 0.0;
        for (k=0; k<ndim; k++) activepart[j].agrav[k] = (FLOAT) 0.0;
      }


      // Loop over all distant pruned trees and compute list of cells.
      // If pruned tree is too close, record cell id for exporting
      //-------------------------------------------------------------------------------------------
      for (j=0; j<Nmpi; j++) {
        if (j == rank) continue;

        Ngravcelltemp = prunedtree[j]->ComputeDistantGravityInteractionList
          (cellptr, simbox, macfactor, Ngravcellmax, Ngravcell, gravcelllist);

        // If pruned tree is too close to be used (flagged by -1), then record cell id
        // for exporting those particles to other MPI processes
        if (Ngravcelltemp == -1) {
#pragma omp critical
          {
            index_cell = Ncellexport[j]++;
          }
          cellexportlist[j][index_cell] = cellptr;
#pragma omp atomic
          Npartexport[j] += Nactive;
        }
        else {
          Ngravcell = Ngravcelltemp;
          assert(Ngravcell <= Ngravcellmax);
        }

      }
      //-------------------------------------------------------------------------------------------


      // Compute gravitational force due to distant cells
      //-------------------------------------------------------------------------------------------
      for (j=0; j<Nactive; j++) {

        if (multipole == "monopole") {
          this->ComputeCellMonopoleForces(activepart[j].gpot, activepart[j].agrav,
                                          activepart[j].r, Ngravcell, gravcelllist);
        }
        else if (multipole == "quadrupole") {
          this->ComputeCellQuadrupoleForces(activepart[j].gpot, activepart[j].agrav,
                                            activepart[j].r, Ngravcell, gravcelllist);
        }

      }
      //-------------------------------------------------------------------------------------------


      // Compute 'fast' multipole terms here
      if (multipole == "fast_monopole") {
        this->ComputeFastMonopoleForces(Nactive, Ngravcell, gravcelllist, cell, activepart);
      }

      // Add all active particles contributions to main array
      for (j=0; j<Nactive; j++) {
        i = activelist[j];
        for (k=0; k<ndim; k++) partdata[i].a[k]     += activepart[j].agrav[k];
        for (k=0; k<ndim; k++) partdata[i].agrav[k] += activepart[j].agrav[k];
        partdata[i].gpot += activepart[j].gpot;
      }

    }
    //=============================================================================================

    // Free-up local memory for OpenMP thread
    delete[] gravcelllist;

  }
  //===============================================================================================

  delete[] celllist;

  timing->EndTimingSection("HYDRO_DISTANT_FORCES");

  return;
}



//=================================================================================================
//  HydroTree::UpdateHydroExportList
/// Check if any local particles need to be exported to other MPI nodes by comparing the
/// smoothing length box overlaps.  If so, then flag for exporting.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void HydroTree<ndim,ParticleType,TreeCell>::UpdateHydroExportList
 (int rank,                            ///< [in] MPI rank
  int Nhydro,                          ///< [in] No. of hydro particles
  int Ntot,                            ///< [in] No. of hydro + ghost particles
  Particle<ndim> *part_gen,            ///< [inout] Pointer to Hydrodynamics ptcl array
  Hydrodynamics<ndim> *hydro,          ///< [in] Pointer to Hydrodynamics object
  Nbody<ndim> *nbody,                  ///< [in] Pointer to N-body object
  const DomainBox<ndim> &simbox)       ///< [in] Simulation domain box
{
  int cactive;                         // No. of active cells
  TreeCell<ndim> **celllist;           // List of pointers to binary tree cells
  ParticleType<ndim>* partdata = static_cast<ParticleType<ndim>* > (part_gen);

  debug2("[HydroTree::UpdateHydroExportList]");
  timing->StartTimingSection("MPI_HYDRO_EXPORT");


  // Find list of all cells that contain active particles
  celllist = new TreeCell<ndim>*[2*tree->gtot];
  cactive = tree->ComputeActiveCellPointers(celllist);

  // Reset all export lists
  for (int j=0; j<Nmpi; j++) {
    Ncellexport[j] = 0;
    Npartexport[j] = 0;
  }


  // Set-up all OMP threads
  //===============================================================================================
#pragma omp parallel default(none) shared(celllist,cactive,cout,rank,partdata)
  {
#if defined _OPENMP
    const int ithread = omp_get_thread_num();
#else
    const int ithread = 0;
#endif
    bool overlapflag;                            // Flag if cells overlap
    int cc;                                      // Aux. cell counter
    int j;                                       // Aux. particle counter
    int *activelist = activelistbuf[ithread];    // List of active particle ids


    // Loop over all active cells
    //=============================================================================================
#pragma omp for schedule(guided)
    for (cc=0; cc<cactive; cc++) {
      TreeCell<ndim> *cellptr = celllist[cc];
      TreeCell<ndim>& cell = *cellptr;

      // Loop over all distant pruned trees and compute list of cells.
      // If pruned tree is too close, record cell id for exporting
      //-------------------------------------------------------------------------------------------
      for (j=0; j<Nmpi; j++) {
        if (j == rank) continue;

        overlapflag = prunedtree[j]->ComputeHydroTreeCellOverlap(cellptr, simbox);

        // If pruned tree is too close (flagged by -1), then record cell id
        // for exporting to other MPI processes
        if (overlapflag) {
          int index_cell;

#pragma omp critical
          {
            index_cell = Ncellexport[j]++;
          }
          cellexportlist[j][index_cell] = cellptr;
          const int Nactive = tree->ComputeActiveParticleList(cell, partdata, activelist);

#pragma omp atomic
          Npartexport[j] += Nactive;
          assert(Ncellexport[j] <= tree->gmax);
        }

      }
      //-------------------------------------------------------------------------------------------

    }
    //=============================================================================================

#ifdef VERIFY_ALL
    PrintArray("Ncellexport : ",Nmpi,Ncellexport);
    PrintArray("Npartexport : ",Nmpi,Npartexport);
#endif


  }
  //===============================================================================================

  delete[] celllist;


#ifdef MPI_PARALLEL
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  timing->EndTimingSection("MPI_HYDRO_EXPORT");

  return;
}



//=================================================================================================
//  HydroTree::BuildPrunedTree
/// Constructs a pruned version of the local tree ready to be exported to other MPI processes.
/// Copies all levels up to and including 'pruning_level'.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void HydroTree<ndim,ParticleType,TreeCell>::BuildPrunedTree
 (const int rank,                      ///< [in] Rank of local MPI node
  const int Nhydromax,                 ///< [in] Max. no. of hydro particles
  const DomainBox<ndim> &simbox,       ///< [in] Simulation domain box object
  const MpiNode<ndim> *mpinode,        ///< [in] Pointer to MPI node array
  Particle<ndim> *hydro_gen)           ///< [inout] Pointer to Hydrodynamics ptcl array
{
  bool localNode;                              // Is this pruned tree for the local node?
  int c;                                       // Cell counter
  int i;                                       // Particle counter
  Tree<ndim,ParticleType,TreeCell> *treeptr;   // Pointer to tree object in question

  debug2("[HydroTree::BuildPrunedTree]");
  timing->StartTimingSection("BUILD_PRUNED_TREE");

  Nprunedcellmax = 0;
#if defined(OUTPUT_ALL)
  cout << "Levels : " << pruning_level_min << "    " << tree->ltot << endl;
#endif

  // Update all work counters in the tree for load-balancing purposes
  tree->UpdateWorkCounters(tree->celldata[0]);


  // Set level at which tree will be pruned (for all trees)
  //-----------------------------------------------------------------------------------------------
  for (i=0; i<Nmpi; i++) {

    // If creating pruned tree for local node, then set pointer to correct array.
    // Otherwise, set pointer so send buffer is filled with pruned tree
    if (i == rank) {
      localNode = true;
      treeptr = prunedtree[i];
    }
    else {
      localNode = false;
      treeptr = sendprunedtree[i];
    }

    treeptr->ltot     = pruning_level_max;
    treeptr->Ntot     = 0;
    treeptr->gmax     = 0;
    treeptr->Ncellmax = max(treeptr->Ncellmax, treeptr->GetMaxCellNumber(pruning_level_min));
    treeptr->AllocateTreeMemory();

    treeptr->Ncell = tree->CreatePrunedTreeForMpiNode
      (mpinode[i], simbox, (FLOAT) 0.0, localNode, pruning_level_min, pruning_level_max,
       treeptr->Ncellmax, treeptr->celldata);


    // If insufficient memory was allocated, then re-allocate larger array and repeat.
    //---------------------------------------------------------------------------------------------
    while (treeptr->Ncell == -1) {

      treeptr->Ncellmax = max(2*treeptr->Ncellmax, treeptr->GetMaxCellNumber(pruning_level_min));
      treeptr->AllocateTreeMemory();
      treeptr->Ncell = tree->CreatePrunedTreeForMpiNode
        (mpinode[i], simbox, (FLOAT) 0.0, localNode, pruning_level_min, pruning_level_max,
         treeptr->Ncellmax, treeptr->celldata);

    }
    //---------------------------------------------------------------------------------------------

    assert(treeptr->Ncellmax > 0);
    assert(treeptr->Ncell > 0 && treeptr->Ncell <= treeptr->Ncellmax);

    // Allocate (or reallocate if needed) all tree memory
    Nprunedcellmax += treeptr->Ncell;

  }
  //-----------------------------------------------------------------------------------------------

  timing->EndTimingSection("BUILD_PRUNED_TREE");
  MPI_Barrier(MPI_COMM_WORLD);

  return;
}



//=================================================================================================
//  HydroTree::BuildMpiGhostTree
/// Main routine to control how the tree is built, re-stocked and interpolated
/// during each timestep.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void HydroTree<ndim,ParticleType,TreeCell>::BuildMpiGhostTree
 (const bool rebuild_tree,             ///< Flag to rebuild tree
  const int n,                         ///< Integer time
  const int ntreebuildstep,            ///< Tree build frequency
  const int ntreestockstep,            ///< Tree stocking frequency
  const int Npart,                     ///< No. of particles
  const int Npartmax,                  ///< Max. no. of particles
  const FLOAT timestep,                ///< Smallest physical timestep
  Particle<ndim> *part_gen,            ///< Particle data array
  Hydrodynamics<ndim> *hydro)          ///< Pointer to Hydrodynamics object
{
  ParticleType<ndim> *partdata = static_cast<ParticleType<ndim>* > (part_gen);


  // If no MPI ghosts exist, do not build tree
  //if (hydroNmpighost == 0) return;

  debug2("[HydroTree::BuildMpiGhostTree]");
  timing->StartTimingSection("BUILD_MPIGHOST_TREE");

  // Activate nested parallelism for tree building routines
#ifdef _OPENMP
  omp_set_nested(1);
#endif

#ifdef OUTPUT_ALL
  cout << "BUILDING TREE WITH " << hydro->Nmpighost << " MPI GHOSTS!!" << endl;
#endif

  // For tree rebuild steps
  //-----------------------------------------------------------------------------------------------
  if (n%ntreebuildstep == 0 || rebuild_tree) {

    mpighosttree->Ntot       = hydro->Nmpighost;
    mpighosttree->Ntotmaxold = mpighosttree->Ntotmax;
    mpighosttree->Ntotmax    = max(mpighosttree->Ntotmax,mpighosttree->Ntot);
    mpighosttree->Ntotmax    = max(mpighosttree->Ntotmax,hydro->Nhydromax);
    mpighosttree->BuildTree(hydro->Nhydro + hydro->NPeriodicGhost,
                            hydro->Nhydro + hydro->NPeriodicGhost + hydro->Nmpighost - 1,
                            mpighosttree->Ntot, mpighosttree->Ntotmax, timestep, partdata);

  }

  // Else stock the tree
  //-----------------------------------------------------------------------------------------------
  else if (n%ntreestockstep == 0) {

    mpighosttree->StockTree(mpighosttree->celldata[0], partdata);

  }

  // Otherwise simply extrapolate tree cell properties
  //-----------------------------------------------------------------------------------------------
  else {

    mpighosttree->ExtrapolateCellProperties(timestep);

  }
  //-----------------------------------------------------------------------------------------------

#ifdef _OPENMP
  omp_set_nested(0);
#endif

  timing->EndTimingSection("BUILD_MPIGHOST_TREE");


  return;
}



//=================================================================================================
//  HydroTree::SearchMpiGhostParticles
/// Search through local domain for any MPI ghost particles that should be sent to the given
/// domain by checking the predicted position over a time tghost.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
int HydroTree<ndim,ParticleType,TreeCell>::SearchMpiGhostParticles
 (const FLOAT tghost,                  ///< [in] Expected ghost life-time
  const Box<ndim> &mpibox,             ///< [in] Bounding box of MPI domain
  Hydrodynamics<ndim> *hydro,          ///< [in] Pointer to Hydrodynamics object
  vector<int> &export_list)            ///< [out] List of particle ids
{
  int c = 0;                           // Cell counter
  int i;                               // Particle counter
  int k;                               // Dimension counter
  int Nexport = 0;                     // No. of MPI ghosts to export
  FLOAT scattermin[ndim];              // Minimum 'scatter' box size due to particle motion
  FLOAT scattermax[ndim];              // Maximum 'scatter' box size due to particle motion
  TreeCell<ndim> *cellptr;             // Pointer to tree cell
  const FLOAT grange = 2.0*ghost_range*kernrange;


  // Start from root-cell of tree and walk all cells
  //-----------------------------------------------------------------------------------------------
  while (c < tree->Ncell) {
    cellptr = &(tree->celldata[c]);

    // Construct maximum cell bounding box depending on particle velocities
    for (k=0; k<ndim; k++) {
      scattermin[k] = cellptr->bbmin[k] +
        min((FLOAT) 0.0, cellptr->v[k]*tghost) - grange*cellptr->hmax;
      scattermax[k] = cellptr->bbmax[k] +
        max((FLOAT) 0.0, cellptr->v[k]*tghost) + grange*cellptr->hmax;
    }


    // If maximum cell scatter box overlaps MPI domain, open cell
    //---------------------------------------------------------------------------------------------
    if (BoxOverlap(ndim, scattermin, scattermax, mpibox.boxmin, mpibox.boxmax)) {

      // If not a leaf-cell, then open cell to first child cell
      if (cellptr->level != tree->ltot) {
        c++;
      }

      else if (cellptr->N == 0) {
        c = cellptr->cnext;
      }

      // If leaf-cell, check through particles in turn to find ghosts and
      // add to list to be exported
      else if (cellptr->level == tree->ltot) {
        i = cellptr->ifirst;
        while (i != -1) {
          export_list.push_back(i);
          Nexport++;
          if (i == cellptr->ilast) break;
          i = tree->inext[i];
        };
        c = cellptr->cnext;
      }
    }

    // If not in range, then open next cell
    //---------------------------------------------------------------------------------------------
    else {
      c = cellptr->cnext;
    }

  }
  //-----------------------------------------------------------------------------------------------



  // Start from root-cell of tree and walk all cells
  //-----------------------------------------------------------------------------------------------
  c = 0;
  while (c < ghosttree->Ncell) {
    cellptr = &(ghosttree->celldata[c]);

    // Construct maximum cell bounding box depending on particle velocities
    for (k=0; k<ndim; k++) {
      scattermin[k] = cellptr->bbmin[k] +
        min((FLOAT) 0.0, cellptr->v[k]*tghost) - grange*cellptr->hmax;
      scattermax[k] = cellptr->bbmax[k] +
        max((FLOAT) 0.0, cellptr->v[k]*tghost) + grange*cellptr->hmax;
    }


    // If maximum cell scatter box overlaps MPI domain, open cell
    //---------------------------------------------------------------------------------------------
    if (BoxOverlap(ndim, scattermin, scattermax, mpibox.boxmin, mpibox.boxmax)) {

      // If not a leaf-cell, then open cell to first child cell
      if (cellptr->level != ghosttree->ltot) {
        c++;
      }

      else if (cellptr->N == 0) {
        c = cellptr->cnext;
      }

      // If leaf-cell, check through particles in turn to find ghosts and
      // add to list to be exported
      else if (cellptr->level == ghosttree->ltot) {
        i = cellptr->ifirst;
        while (i != -1) {
          export_list.push_back(i);
          Nexport++;
          if (i == cellptr->ilast) break;
          i = ghosttree->inext[i];
        };
        c = cellptr->cnext;
      }
    }

    // If not in range, then open next cell
    //---------------------------------------------------------------------------------------------
    else {
      c = cellptr->cnext;
    }

  }
  //-----------------------------------------------------------------------------------------------


  return Nexport;
}



//=================================================================================================
//  HydroTree::FindMpiTransferParticles
/// ..
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void HydroTree<ndim,ParticleType,TreeCell>::FindMpiTransferParticles
 (Hydrodynamics<ndim> *hydro,                ///< [in] Pointer to Hydrodynamics class
  vector<vector<int> >& particles_to_export, ///< [inout] Vector that for each
                                             ///< node gives the list of particles to export
  vector<int>& all_particles_to_export,      ///< [inout] Vector containing all the particles
                                             ///<         that will be exported by this processor
  const vector<int>& potential_nodes,        ///< [in] Vector containing the potential nodes we
                                             ///<      might be sending particles to
  MpiNode<ndim>* mpinodes)                   ///< [in] Array of other mpi nodes
{
  int c;                                     // Cell counter
  int i;                                     // Particle counter
  int jnode;                                 // Aux. node counter
  int inode;                                 // MPI node id
  TreeCell<ndim> *cellptr;                   // Pointer to cell
  ParticleType<ndim> *partdata = static_cast<ParticleType<ndim>* > (hydro->GetParticleArray());


  // Loop over potential domains and walk the tree for each bounding box
  //-----------------------------------------------------------------------------------------------
  for (jnode=0; jnode<potential_nodes.size(); jnode++) {

    inode = potential_nodes[jnode];
    Box<ndim>& nodebox = mpinodes[inode].domain;

    // Start from root-cell
    c = 0;

    //---------------------------------------------------------------------------------------------
    while (c < tree->Ncell) {
      cellptr = &(tree->celldata[c]);

      // If maximum cell scatter box overlaps MPI domain, open cell
      //-------------------------------------------------------------------------------------------
      if (BoxOverlap(ndim, cellptr->bbmin, cellptr->bbmax, nodebox.boxmin, nodebox.boxmax)) {

        // If not a leaf-cell, then open cell to first child cell
        if (cellptr->level != tree->ltot) {
          c++;
        }

        else if (cellptr->N == 0) {
          c = cellptr->cnext;
        }

        // If leaf-cell, check through particles in turn to find ghosts and
        // add to list to be exported
        else if (cellptr->level == tree->ltot) {
          i = cellptr->ifirst;
          while (i != -1) {
            if (ParticleInBox(partdata[i], mpinodes[inode].domain)) {
              particles_to_export[inode].push_back(i);
              all_particles_to_export.push_back(i);
            }
            if (i == cellptr->ilast) break;
            i = tree->inext[i];
          };
          c = cellptr->cnext;
        }
      }

      // If not in range, then open next cell
      //-------------------------------------------------------------------------------------------
      else {
        c = cellptr->cnext;
      }

    }
    //---------------------------------------------------------------------------------------------

  }
  //-----------------------------------------------------------------------------------------------

  return;
}



//=================================================================================================
//  HydroTree::FindLoadBalancingDivision
/// Find the predicted cell-cell dividing point that approximately balances the CPU work load
/// in order to achieve load balancing amongst all MPI nodes.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
FLOAT HydroTree<ndim,ParticleType,TreeCell>::FindLoadBalancingDivision
 (int k_divide,                        ///< Dimension of cell division
  FLOAT r_old,                         ///< Old position of cell division
  FLOAT boxmin[ndim],                  ///< Minimum extent of parent MPI tree cell
  FLOAT boxmax[ndim])                  ///< Maximum extent of parent MPI tree cell
{
  int i;                               // MPI node counter
  int k;                               // Dimension counter
  FLOAT r_divide = r_old;              // Cell division location
  FLOAT r_max = boxmax[k_divide];      // Max. for bisection iteration of division
  FLOAT r_min = boxmin[k_divide];      // Min. for bisection iteration of division
  FLOAT workleft;                      // Work computed on LHS of division
  FLOAT workright;                     // Work computed on RHS of division
  FLOAT workfrac;                      // Fraction of work on LHS
  FLOAT worktol = (FLOAT) 0.001;        // Work balance tolerance for iteration
  FLOAT boxleftmin[ndim];              // Min. extent of left-box division
  FLOAT boxleftmax[ndim];              // Max. extent of left-box division
  FLOAT boxrightmin[ndim];             // Min. extent of right-box division
  FLOAT boxrightmax[ndim];             // Max. extent of right-box division

  // Set box extents from MPI tree node extent
  for (k=0; k<ndim; k++) {
    boxleftmin[k] = boxmin[k];
    boxleftmax[k] = boxmax[k];
    boxrightmin[k] = boxmin[k];
    boxrightmax[k] = boxmax[k];
  }


  // Find the work-balance position through bisection iteration
  //-----------------------------------------------------------------------------------------------
  do {
    boxleftmax[k_divide] = r_divide;
    boxrightmin[k_divide] = r_divide;
    workleft = (FLOAT) 0.0;
    workright = (FLOAT) 0.0;

    // Compute work included in left-hand side from pruned trees of all MPI domains
    for (i=0; i<Nmpi; i++) {
      workleft += prunedtree[i]->ComputeWorkInBox(boxleftmin, boxleftmax);
    }

    // Compute work included in right-hand side from pruned trees of all MPI domains
    for (i=0; i<Nmpi; i++) {
      workright += prunedtree[i]->ComputeWorkInBox(boxrightmin, boxrightmax);
    }

#ifdef OUTPUT_ALL
    cout << "workleft : " << workleft << "     workright : " << workright
         << "       workfrac : " << workleft / (workleft + workright)
         << "    rold : " << r_old << "     r_new : " << r_divide << endl;
#endif

    // If fraction of work on either side of division is too inbalanced, calculate new position
    // of division and perform another iteration.  Otherwise exit iteration loop.
    workfrac = workleft / (workleft + workright);
    if (workfrac < 0.5 - worktol) r_min = r_divide;
    else if (workfrac > 0.5 + worktol) r_max = r_divide;
    else break;

    r_divide = 0.5*(r_min + r_max);

  } while (fabs(workfrac - 0.5) > worktol);
  //-----------------------------------------------------------------------------------------------

  return r_divide;
}



//=================================================================================================
//  HydroTree::FindParticlesToTransfer
/// Compute on behalf of the MpiControl class the particles that are outside
/// the domain after a load balancing step and need to be transferred to other nodes
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void HydroTree<ndim,ParticleType,TreeCell>::FindParticlesToTransfer
 (Hydrodynamics<ndim> *hydro,                ///< [in] Pointer to sph class
  vector<vector<int> >& id_export_buffers,   ///< [inout] List of ids to export for each node
  vector<int>& all_ids_export_buffer,        ///< [inout] List of all ids to export from proc
  const vector<int>& potential_nodes,        ///< [in] Potential nodes we might send particles to
  MpiNode<ndim>* mpinodes)                   ///< [in] Array of other mpi nodes
{
  ParticleType<ndim> *partdata = static_cast<ParticleType<ndim>* > (hydro->GetParticleArray());

  // Loop over particles and prepare the ones to export
  for (int i=0; i<hydro->Nhydro; i++) {
    ParticleType<ndim>& part = partdata[i];

    // Loop over potential domains and see if we need to transfer this particle to them
    for (int inode=0; inode<potential_nodes.size(); inode++) {
      int node_number = potential_nodes[inode];

      if (ParticleInBox(part, mpinodes[node_number].domain)) {
        id_export_buffers[node_number].push_back(i);
        all_ids_export_buffer.push_back(i);

        // The particle can belong only to one domain, so we can break from this loop
        break;
      }
    }
  }

  return;
}



//=================================================================================================
//  HydroTree::GetExportInfo
/// Get the array with the information that needs to be exported to the given processor
/// (NB: iproc is ignored at the moment, as we must export all ptcls to all other processors)
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
int HydroTree<ndim,ParticleType,TreeCell>::GetExportInfo
 (int iproc,                           ///< [in] No. of processor we want to send data to
  Hydrodynamics<ndim> *hydro,          ///< [in] Pointer to GetParticleArray object
  vector<char >& send_buffer,          ///< [inout] Vector where the ptcls to export will be stored
  MpiNode<ndim>& mpinode,              ///< ..
  int rank,                            ///< ..
  int Nmpi)                            ///< [in] Array with information for the other mpi nodes
{
  int cactive = Ncellexport[iproc];
  int Nactive = Npartexport[iproc];
  int activelist[Nleafmax];
  int exported_particles    = 0;
  const int size_header     = 2*sizeof(int);
  const int size_particles  = Nactive*sizeof(ParticleType<ndim>);
  const int size_cells      = cactive*sizeof(TreeCell<ndim>);
  const int old_size        = send_buffer.size();
  int offset                = size_header + old_size;
  TreeCell<ndim>** celllist = cellexportlist[iproc];
  ParticleType<ndim>* partdata = static_cast<ParticleType<ndim>* > (hydro->GetParticleArray());

  assert(tree->Nimportedcell == 0);

  // Work out size of the information we are sending
  // Header consists of number of particles and number of cells
  send_buffer.resize(size_particles + size_cells + size_header + old_size);

  // Write the header
  copy(&send_buffer[old_size], &Nactive);
  copy(&send_buffer[old_size+sizeof(int)], &cactive);

  // Clear the array needed for bookkeeping (which active particles we sent to which processor)
  vector<int>& ids_active_particles = ids_sent_particles[iproc];
  ids_active_particles.clear();
  ids_active_particles.reserve(Nactive);


  // Loop over all cells to be exported and include all cell and particle data
  //-----------------------------------------------------------------------------------------------
  for (int cc=0; cc<cactive; cc++) {
    copy(&send_buffer[offset], celllist[cc]);
    const int Nactive_cell = tree->ComputeActiveParticleList(*(celllist[cc]), partdata, activelist);

    // Update the ifirst and ilast pointers in the cell
    TreeCell<ndim>* exported_cell = reinterpret_cast<TreeCell<ndim>*> (&send_buffer[offset]);
    exported_cell->ifirst = exported_particles;
    exported_cell->ilast  = exported_particles + Nactive_cell - 1;
    offset += sizeof(TreeCell<ndim>);

    // Copy active particles
    for (int jpart=0; jpart<Nactive_cell; jpart++) {
      ids_active_particles.push_back(activelist[jpart]);
      copy(&send_buffer[offset], &partdata[activelist[jpart]]);
      offset += sizeof(ParticleType<ndim>);
    }
    exported_particles += Nactive_cell;
  }
  //-----------------------------------------------------------------------------------------------

  assert(exported_particles == Nactive);
  assert(offset == send_buffer.size());

  return size_particles + size_cells + size_header;
}



//=================================================================================================
//  HydroTree::UnpackExported
/// Unpack the information exported from the other processors, contaning the particles
/// that were exported and
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void HydroTree<ndim,ParticleType,TreeCell>::UnpackExported
 (vector<char >& received_array,
  vector<int>& Nbytes_exported_from_proc,
  Hydrodynamics<ndim> *hydro)
{
  int offset = 0;
  ParticleType<ndim>* partdata = static_cast<ParticleType<ndim>* > (hydro->GetParticleArray());

  assert(hydro->NImportedParticles == 0);
  tree->Nimportedcell = 0;
  tree->Ncelltot=tree->Ncell;

  N_imported_part_per_proc.resize(Nbytes_exported_from_proc.size());


  //-----------------------------------------------------------------------------------------------
  for (int iproc = 0; iproc<Nbytes_exported_from_proc.size(); iproc++) {

    int N_received_bytes = Nbytes_exported_from_proc[iproc];
    int N_received_cells;
    int N_received_particles;

    if (N_received_bytes == 0) {
      N_imported_part_per_proc[iproc] = 0;
      continue;
    }

    copy(&N_received_particles, &received_array[offset]);
    N_imported_part_per_proc[iproc] = N_received_particles;
    copy(&N_received_cells, &received_array[offset + sizeof(int)]);

    // Ensure there is enough memory
    if (hydro->Ntot + N_received_particles > hydro->Nhydromax) {
      ExceptionHandler::getIstance().raise("Error while receiving imported particles: not enough memory!");
    }
    if (tree->Ncelltot + N_received_cells > tree->Ncellmax) {
      ExceptionHandler::getIstance().raise("Error while receiving imported cells: not enough memory!");
    }


    // Copy received particles inside main arrays and received cells inside the tree array
    // Also update the linked list
    int particle_index = hydro->Ntot;
    offset += 2*sizeof(int);

    //---------------------------------------------------------------------------------------------
    for (int icell=0; icell<N_received_cells; icell++) {
      TreeCell<ndim>& dest_cell = tree->celldata[icell + tree->Ncelltot];
      copy(&dest_cell, &received_array[offset]);
      offset += sizeof(TreeCell<ndim>);

      // Offset the ifirst and ilast pointers
      dest_cell.ifirst += hydro->Ntot;
      dest_cell.ilast += hydro->Ntot;

      // Now copy the received particles inside the hydro particle main arrays
      for (int iparticle=0; iparticle<dest_cell.Nactive; iparticle++) {
        copy(&partdata[particle_index], &received_array[offset]);
        partdata[particle_index].gpot=0;
        for (int k=0; k<ndim; k++) {
          partdata[particle_index].a[k]=0;
          partdata[particle_index].agrav[k]=0;
        }
        tree->inext[particle_index] = particle_index + 1;
        particle_index++;
        offset += sizeof(ParticleType<ndim>);
      }

    }
    //---------------------------------------------------------------------------------------------

    // Update the hydro counters
    hydro->Ntot += N_received_particles;
    hydro->NImportedParticles += N_received_particles;

    // Update the tree counters
    tree->Nimportedcell += N_received_cells;
    tree->Ncelltot      += N_received_cells;
    tree->Ntot           = hydro->Ntot;

  }
  //-----------------------------------------------------------------------------------------------


  assert (offset == std::accumulate(Nbytes_exported_from_proc.begin(),
                                    Nbytes_exported_from_proc.end(), 0));

  return;
}


//=================================================================================================
//  HydroTree::GetBackExportInfo
/// Return the data to transmit back to the other processors (particle acceleration etc.)
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void HydroTree<ndim,ParticleType,TreeCell>::GetBackExportInfo
 (vector<char >& send_buffer,              ///< [inout] These arrays will be overwritten with the information to send
  vector<int>& Nbytes_exported_from_proc,  ///< ..
  vector<int>& Nbytes_to_each_proc,        ///< ..
  Hydrodynamics<ndim> *hydro,              ///< [in] Pointer to the GetParticleArray object
  int rank)                                ///< ..
{
  int removed_particles = 0;               // ..
  int InitialNImportedParticles = hydro->NImportedParticles;
  ParticleType<ndim>* partdata = static_cast<ParticleType<ndim>* > (hydro->GetParticleArray());

  // Loop over the processors, removing particles as we go
  send_buffer.resize(hydro->NImportedParticles * sizeof(ParticleType<ndim>));

  //-----------------------------------------------------------------------------------------------
  for (int iproc=0; iproc<N_imported_part_per_proc.size(); iproc++ ) {

    // Copy the particles inside the send buffer
    const int N_received_particles = N_imported_part_per_proc[iproc];
    const int start_index = hydro->Nhydro + hydro->Nghost + removed_particles;
    int j = 0;

    for (int i=start_index; i<start_index + N_received_particles; i++) {
      copy (&send_buffer[(removed_particles + j)*sizeof(ParticleType<ndim>)], &partdata[i]);
      j++;
    }
    assert(j == N_received_particles);
    removed_particles += j;

    // Decrease the particle counter
    hydro->Ntot -= N_received_particles;
    hydro->NImportedParticles -= N_received_particles;

    // Update the information about how much data we are sending
    Nbytes_exported_from_proc[iproc] = N_received_particles*sizeof(ParticleType<ndim>);

    // Update the information with how much data we are receiving
    Nbytes_to_each_proc[iproc] = ids_sent_particles[iproc].size()*sizeof(ParticleType<ndim>);

  }
  //-----------------------------------------------------------------------------------------------

  tree->Ncelltot      = tree->Ncell;
  tree->Nimportedcell = 0;
  tree->Ntot          -= InitialNImportedParticles;

  assert(hydro->NImportedParticles == 0);
  assert(hydro->Ntot == hydro->Nhydro + hydro->Nghost);
  assert(send_buffer.size() == removed_particles*sizeof(ParticleType<ndim>));

  return;
}



//=================================================================================================
//  HydroTree::UnpackReturnedExportInfo
/// Unpack the data that was returned by the other processors, summing the accelerations to the particles
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void HydroTree<ndim,ParticleType,TreeCell>::UnpackReturnedExportInfo
 (vector<char >& received_information,   ///< ..
  vector<int>& recv_displs,              ///< ..
  Hydrodynamics<ndim> *hydro,            ///< ..
  int rank)                              ///< ..
{
  ParticleType<ndim>* partdata = static_cast<ParticleType<ndim>* > (hydro->GetParticleArray() );


  // For each processor, sum up the received quantities for each particle
  //-----------------------------------------------------------------------------------------------
  for (int iproc=0; iproc<recv_displs.size(); iproc++ ) {

    // (IS THIS CORRECT WHEN NOT ALL PROCS SEND PARTICLES?)
    if (rank == iproc) continue;

    const vector<int>& ids_active_particles = ids_sent_particles[iproc];

    for (int j=0; j<ids_active_particles.size(); j++) {
      const int i = ids_active_particles[j];

      ParticleType<ndim>* received_particle = reinterpret_cast<ParticleType<ndim>*>
        (&received_information[j*sizeof(ParticleType<ndim>) + recv_displs[iproc] ]);

      assert(partdata[i].iorig == received_particle->iorig);

      for (int k=0; k<ndim; k++) {
        partdata[i].a[k]     += received_particle->a[k];
        partdata[i].agrav[k] += received_particle->agrav[k];
      }
      partdata[i].gpot  += received_particle->gpot;
      partdata[i].dudt  += received_particle->dudt;
      partdata[i].div_v += received_particle->div_v;
      partdata[i].levelneib = max(partdata[i].levelneib, received_particle->levelneib);
    }

  }
  //-----------------------------------------------------------------------------------------------


  return;
}



//=================================================================================================
//  HydroTree::CommunicatePrunedTrees
/// Communicate the pruned trees from each MPI process to the others.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void HydroTree<ndim,ParticleType,TreeCell>::CommunicatePrunedTrees
 (vector<int>& my_matches,             ///< [in] Vector containing communication schedule
  int rank)                            ///< [in] MPI rank
{
  bool send_turn;                      // Flag if this node is sending (true) or receiving (false)

  MPI_Barrier(MPI_COMM_WORLD);

  //-----------------------------------------------------------------------------------------------
  for (int iturn=0; iturn<my_matches.size(); iturn++) {

    // See who we need to communicate with
    int inode = my_matches[iturn];
#if defined(OUTPUT_ALL)
    cout << "inode : " << rank << "   " << inode << "    " << iturn << endl;
#endif

    // Decide if sending or receiving first
    if (rank < inode) {
      send_turn = true;
    }
    else {
      send_turn = false;
    }

    // Do the actual communication
    for (int i=0; i<2; i++) {

      if (send_turn) {
        Tree<ndim,ParticleType,TreeCell>* treeptr = sendprunedtree[inode];
        MPI_Send(&(treeptr->Ncell), 1, MPI_INT, inode, 3, MPI_COMM_WORLD);
        MPI_Send(treeptr->celldata, treeptr->Ncell*sizeof(TreeCell<ndim>),
                 MPI_CHAR, inode, 3, MPI_COMM_WORLD);
        send_turn = false;
      }
      else {
        Tree<ndim,ParticleType,TreeCell>* treeptr = prunedtree[inode];
        MPI_Status status;
        MPI_Recv(&(treeptr->Ncell), 1, MPI_INT, inode, 3, MPI_COMM_WORLD, &status);
        treeptr->Ncellmax = max(treeptr->Ncellmax, treeptr->Ncell);
        treeptr->AllocateTreeMemory();
        MPI_Recv(treeptr->celldata, treeptr->Ncell*sizeof(TreeCell<ndim>),
                 MPI_CHAR, inode, 3, MPI_COMM_WORLD, &status);
        send_turn = true;

#ifdef VERIFY_ALL
        int cnext;
#if defined(OUTPUT_ALL)
        cout << "Node : " << inode << "   Ncell received : " << treeptr->Ncell << endl;
#endif

        for (int c=0; c<treeptr->Ncell; c++) {
          if (treeptr->celldata[c].copen != -1) {
            cnext = treeptr->celldata[c].copen;
            assert(cnext == c + 1);
          }
          else {
            cnext = treeptr->celldata[c].cnext;
            assert(cnext > c);
          }
          assert(treeptr->celldata[c].cnext > 0);
          assert(treeptr->celldata[c].copen < treeptr->celldata[c].cnext);
        }
#endif

      }
    }

  }
  //----------------------------------------------------------------------------------------------

  MPI_Barrier(MPI_COMM_WORLD);
#if defined(OUTPUT_ALL)
  int k;
  FLOAT rcom[ndim],mcom = 0.0;
  for (k=0; k<ndim; k++) rcom[k] = 0.0;
  for (int i=0; i<Nmpi; i++) {
    assert(prunedtree[i]->Ncell > 0);
    assert(prunedtree[i]->Ncellmax > 0);
    mcom += prunedtree[i]->celldata[0].m;
    for (k=0; k<ndim; k++) rcom[k] += prunedtree[i]->celldata[0].m*prunedtree[i]->celldata[0].r[k];
    //if (i == rank) continue;
    cout << "Writing pruned tree " << i << " for process " << rank << endl;
    cout << "Ncell : " << prunedtree[i]->Ncell << "   " << prunedtree[i]->Ncellmax << endl;
    cout << "r : " << prunedtree[i]->celldata[0].r[0]
         << "   box : " << prunedtree[i]->celldata[0].bbmin[0]
         << "   " << prunedtree[i]->celldata[0].bbmax[0]
         << "   worktot : " << prunedtree[i]->celldata[0].worktot << endl;
  }
  cout << "mcom : " << mcom << endl;
  for (k=0; k<ndim; k++) cout << "rcom[" << k << "] : " << rcom[k] << endl;
#endif

  return;
}



//=================================================================================================
//  HydroTree::InitialiseCellWorkCounters
/// Initialise all CPU work counters (used for MPI load balancing) to zero for next balance.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void HydroTree<ndim,ParticleType,TreeCell>::InitialiseCellWorkCounters(void)
{
  assert(tree->Ncellmax > 0);
  for (int c=0; c<tree->Ncellmax; c++) tree->celldata[c].worktot = (FLOAT) 0.0;
  return;
}
#endif



#if defined(VERIFY_ALL)
//=================================================================================================
//  NeighbourSearch::CheckValidNeighbourList
/// Checks that the neighbour list generated by the grid is valid in that it
/// (i) does include all true neighbours, and
/// (ii) all true neigbours are only included once and once only.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void HydroTree<ndim,ParticleType,TreeCell>::CheckValidNeighbourList
 (int i,                               ///< [in] Particle i.d.
  int Ntot,                            ///< [in] Total no. of particles
  int Nneib,                           ///< [in] No. of potential neighbours
  int *neiblist,                       ///< [in] List of potential neighbour i.d.s
  ParticleType<ndim> *partdata,        ///< [in] Array of particle data
  string neibtype)                     ///< [in] Neighbour search type
{
  bool invalid_flag = false;           // Flag if neighbour list is invalid
  int count;                           // Valid neighbour counter
  int j;                               // Neighbour particle counter
  int k;                               // Dimension counter
  int Ntrueneib = 0;                   // No. of 'true' neighbours
  int *trueneiblist;                   // List of true neighbour ids
  FLOAT drsqd;                         // Distance squared
  FLOAT dr[ndim];                      // Relative position vector

  // Allocate array to store local copy of potential neighbour ids
  trueneiblist = new int[Ntot];


  // First, create list of 'true' neighbours by looping over all particles
  if (neibtype == "gather") {
    for (j=0; j<Ntot; j++) {
      for (k=0; k<ndim; k++) dr[k] = partdata[j].r[k] - partdata[i].r[k];
      drsqd = DotProduct(dr,dr,ndim);
      if (drsqd <= kernrangesqd*partdata[i].h*partdata[i].h) trueneiblist[Ntrueneib++] = j;
    }
  }
  else if (neibtype == "all") {
    for (j=0; j<Ntot; j++) {
      for (k=0; k<ndim; k++) dr[k] = partdata[j].r[k] - partdata[i].r[k];
      drsqd = DotProduct(dr,dr,ndim);
      if (drsqd < kernrangesqd*partdata[i].h*partdata[i].h ||
          drsqd < kernrangesqd*partdata[j].h*partdata[j].h) trueneiblist[Ntrueneib++] = j;
    }
  }


  // Now compare each given neighbour with true neighbour list for validation
  for (j=0; j<Ntrueneib; j++) {
    count = 0;
    for (k=0; k<Nneib; k++) {
      if (neiblist[k] == trueneiblist[j]) {
        count++;
      }
    }

    // If the true neighbour is not in the list, or included multiple times,
    // then output to screen and terminate program
    if (count != 1) {
      for (k=0; k<ndim; k++) dr[k] = partdata[trueneiblist[j]].r[k] - partdata[i].r[k];
      drsqd = DotProduct(dr,dr,ndim);
      cout << "Could not find neighbour " << j << "   " << trueneiblist[j] << "     " << i
           << "      " << sqrt(drsqd)/kernrange/partdata[i].h << "     "
           << sqrt(drsqd)/kernrange/partdata[trueneiblist[j]].h << "    "
           << partdata[trueneiblist[j]].r[0] << "   type : "
           << partdata[trueneiblist[j]].itype << endl;
      invalid_flag = true;
    }

  }


  // If the true neighbour is not in the list, or included multiple times,
  // then output to screen and terminate program
  if (invalid_flag) {
    cout << "Problem with neighbour lists : " << i << "  " << j << "   "
         << count << "   " << partdata[i].r[0] << "   " << partdata[i].h << endl;
    cout << "Nneib : " << Nneib << "   Ntrueneib : " << Ntrueneib
         << "    neibtype : " << neibtype << endl;
    InsertionSort(Nneib,neiblist);
    PrintArray("neiblist     : ",Nneib,neiblist);
    PrintArray("trueneiblist : ",Ntrueneib,trueneiblist);
    string message = "Problem with neighbour lists in tree search";
    ExceptionHandler::getIstance().raise(message);
  }


  delete[] trueneiblist;

  return;
}
#endif


template class HydroTree<1,SphParticle,KDTreeCell>;
template class HydroTree<2,SphParticle,KDTreeCell>;
template class HydroTree<3,SphParticle,KDTreeCell>;
template class HydroTree<1,SphParticle,OctTreeCell>;
template class HydroTree<2,SphParticle,OctTreeCell>;
template class HydroTree<3,SphParticle,OctTreeCell>;

template class HydroTree<1,GradhSphParticle,KDTreeCell>;
template class HydroTree<2,GradhSphParticle,KDTreeCell>;
template class HydroTree<3,GradhSphParticle,KDTreeCell>;
template class HydroTree<1,GradhSphParticle,OctTreeCell>;
template class HydroTree<2,GradhSphParticle,OctTreeCell>;
template class HydroTree<3,GradhSphParticle,OctTreeCell>;

template class HydroTree<1,SM2012SphParticle,KDTreeCell>;
template class HydroTree<2,SM2012SphParticle,KDTreeCell>;
template class HydroTree<3,SM2012SphParticle,KDTreeCell>;
template class HydroTree<1,SM2012SphParticle,OctTreeCell>;
template class HydroTree<2,SM2012SphParticle,OctTreeCell>;
template class HydroTree<3,SM2012SphParticle,OctTreeCell>;

template class HydroTree<1,MeshlessFVParticle,KDTreeCell>;
template class HydroTree<2,MeshlessFVParticle,KDTreeCell>;
template class HydroTree<3,MeshlessFVParticle,KDTreeCell>;
template class HydroTree<1,MeshlessFVParticle,OctTreeCell>;
template class HydroTree<2,MeshlessFVParticle,OctTreeCell>;
template class HydroTree<3,MeshlessFVParticle,OctTreeCell>;

template class HydroTree<1,GradhSphParticle,TreeRayCell>;
template class HydroTree<2,GradhSphParticle,TreeRayCell>;
template class HydroTree<3,GradhSphParticle,TreeRayCell>;
