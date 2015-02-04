//=================================================================================================
//  SphTree.cpp
//  Contains all functions for building, stocking and walking for the
//  binary KD tree for SPH particles.
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
#include "SphNeighbourSearch.h"
#include "Sph.h"
#include "Parameters.h"
#include "InlineFuncs.h"
#include "Particle.h"
#include "Debug.h"
#if defined _OPENMP
#include <omp.h>
#endif
using namespace std;




//=================================================================================================
//  SphTree::SphTree
/// SphTree constructor.  Initialises various variables.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
SphTree<ndim,ParticleType,TreeCell>::SphTree
 (int Nleafmaxaux,
  int Nmpiaux,
  FLOAT thetamaxsqdaux,
  FLOAT kernrangeaux,
  FLOAT macerroraux,
  string gravity_mac_aux,
  string multipole_aux,
  DomainBox<ndim> *boxaux,
  SmoothingKernel<ndim> *kernaux,
  CodeTiming *timingaux):
  SphNeighbourSearch<ndim>(kernrangeaux,boxaux,kernaux,timingaux),
  Nleafmax(Nleafmaxaux),
  Nmpi(Nmpiaux),
  thetamaxsqd(thetamaxsqdaux),
  invthetamaxsqd((FLOAT) 1.0/thetamaxsqdaux),
  gravity_mac(gravity_mac_aux),
  macerror(macerroraux),
  multipole(multipole_aux)
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
//  SphTree::~SphTree
/// SphTree destructor.  Deallocates tree memory upon object destruction.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
SphTree<ndim,ParticleType,TreeCell>::~SphTree()
{
  if (tree->allocated_tree) {
    DeallocateMemory();
    tree->DeallocateTreeMemory();
  }
}



//=================================================================================================
//  SphTree::AllocateMemory
/// Allocate memory for tree as requested.  If more memory is required
/// than currently allocated, tree is deallocated and reallocated here.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void SphTree<ndim,ParticleType,TreeCell>::AllocateMemory
 (Sph<ndim> *sph)                      ///< Pointer to SPH object
{
  int ithread;                         // Thread id number

  debug2("[SphTree::AllocateMemory]");

  if (!allocated_buffer) {

    Nneibmaxbuf     = new int[Nthreads];
    Ngravcellmaxbuf = new int[Nthreads];
    levelneibbuf    = new int*[Nthreads];
    activelistbuf   = new int*[Nthreads];
    activepartbuf   = new ParticleType<ndim>*[Nthreads];
    neibpartbuf     = new ParticleType<ndim>*[Nthreads];
    cellbuf         = new TreeCell<ndim>*[Nthreads];

    for (ithread=0; ithread<Nthreads; ithread++) {
      Nneibmaxbuf[ithread]     = max(1,4*sph->Ngather);
      Ngravcellmaxbuf[ithread] = max(1,4*sph->Ngather);
      levelneibbuf[ithread]    = new int[Ntotmax];
      activelistbuf[ithread]   = new int[Nleafmax];
      activepartbuf[ithread]   = new ParticleType<ndim>[Nleafmax];
      neibpartbuf[ithread]     = new ParticleType<ndim>[Nneibmaxbuf[ithread]];
      cellbuf[ithread]         = new TreeCell<ndim>[Ngravcellmaxbuf[ithread]];
    }
    allocated_buffer = true;

  }

  return;
}



//=================================================================================================
//  SphTree::DeallocateTreeMemory
/// Deallocates all binary tree memory
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void SphTree<ndim,ParticleType,TreeCell>::DeallocateMemory(void)
{
  int ithread;                         // Thread id number

  debug2("[SphTree::DeallocateTreeMemory]");

  if (allocated_buffer) {

    for (ithread=0; ithread<Nthreads; ithread++) {
      delete[] cellbuf[ithread];
      delete[] neibpartbuf[ithread];
      delete[] activepartbuf[ithread];
      delete[] activelistbuf[ithread];
      delete[] levelneibbuf[ithread];
    }
    delete[] neibpartbuf;
    delete[] activepartbuf;
    delete[] activelistbuf;
    delete[] levelneibbuf;
    delete[] Ngravcellmaxbuf;
    delete[] Nneibmaxbuf;

  }

  return;
}



//=================================================================================================
//  SphTree::BuildTree
/// Main routine to control how the tree is built, re-stocked and interpolated during each step.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void SphTree<ndim,ParticleType,TreeCell>::BuildTree
 (bool rebuild_tree,                   ///< [in] Flag to rebuild tree
  int n,                               ///< [in] Integer time
  int ntreebuildstep,                  ///< [in] Tree build frequency
  int ntreestockstep,                  ///< [in] Tree stocking frequency
  int Npart,                           ///< [in] No. of particles
  int Npartmax,                        ///< [in] Max. no. of particles
  SphParticle<ndim> *sph_gen,          ///< [inout] Particle data array
  Sph<ndim> *sph,                      ///< [inout] Pointer to SPH object
  FLOAT timestep)                      ///< [in] Smallest physical timestep
{
  ParticleType<ndim> *sphdata = static_cast<ParticleType<ndim>* > (sph_gen);

  debug2("[SphTree::BuildTree]");
  timing->StartTimingSection("BUILD_TREE",2);

  // Activate nested parallelism for tree building routines
#ifdef _OPENMP
  omp_set_nested(1);
#endif


  // For tree rebuild steps
  //-----------------------------------------------------------------------------------------------
  if (n%ntreebuildstep == 0 || rebuild_tree) {

    // Delete any dead particles from main SPH arrays before we re-build tree
    sph->DeleteDeadParticles();

    Ntotold    = Ntot;
    Ntot       = sph->Ntot;
    Ntotmaxold = Ntotmax;
    Ntotmax    = max(Ntotmax,Ntot);
    Ntotmax    = max(Ntotmax,sph->Nhydromax);
    assert(Ntotmax >= Ntot);

    tree->Ntot       = sph->Nhydro;
    tree->Ntotmaxold = tree->Ntotmax;
    tree->Ntotmax    = max(tree->Ntotmax,tree->Ntot);
    tree->Ntotmax    = max(tree->Ntotmax,sph->Nhydromax);
    tree->BuildTree(0, sph->Nhydro-1, Npart, Npartmax, sphdata, timestep);

    AllocateMemory(sph);
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

    tree->StockTree(tree->celldata[0],sphdata);

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
//  SphTree::BuildGhostTree
/// Main routine to control how the tree is built, re-stocked and interpolated
/// during each timestep.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void SphTree<ndim,ParticleType,TreeCell>::BuildGhostTree
 (bool rebuild_tree,                   ///< Flag to rebuild tree
  int n,                               ///< Integer time
  int ntreebuildstep,                  ///< Tree build frequency
  int ntreestockstep,                  ///< Tree stocking frequency
  int Npart,                           ///< No. of particles
  int Npartmax,                        ///< Max. no. of particles
  SphParticle<ndim> *sph_gen,          ///< Particle data array
  Sph<ndim> *sph,                      ///< Pointer to SPH object
  FLOAT timestep)                      ///< Smallest physical timestep
{
  ParticleType<ndim> *sphdata = static_cast<ParticleType<ndim>* > (sph_gen);

  // If no periodic ghosts exist, do not build tree
  //if (sph->NPeriodicGhost == 0) return;

  debug2("[SphTree::BuildGhostTree]");
  timing->StartTimingSection("BUILD_GHOST_TREE",2);

  // Activate nested parallelism for tree building routines
#ifdef _OPENMP
  omp_set_nested(1);
#endif


  // For tree rebuild steps
  //-----------------------------------------------------------------------------------------------
  if (n%ntreebuildstep == 0 || rebuild_tree) {

    ghosttree->Ntot       = sph->NPeriodicGhost;
    ghosttree->Ntotmaxold = ghosttree->Ntotmax;
    ghosttree->Ntotmax    = max(ghosttree->Ntotmax,ghosttree->Ntot);
    ghosttree->Ntotmax    = max(ghosttree->Ntotmax,sph->Nhydromax);
    ghosttree->BuildTree(sph->Nhydro, sph->Nhydro + sph->NPeriodicGhost - 1,
                         ghosttree->Ntot, ghosttree->Ntotmax, sphdata, timestep);

  }

  // Else stock the tree
  //-----------------------------------------------------------------------------------------------
  else if (n%ntreestockstep == 0) {

    ghosttree->StockTree(ghosttree->celldata[0],sphdata);

  }

  // Otherwise simply extrapolate tree cell properties
  //-----------------------------------------------------------------------------------------------
  else {

    //ExtrapolateCellProperties(celldata[0],timestep);
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
//  SphTree::GetGatherNeighbourList
/// ..
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
int SphTree<ndim,ParticleType,TreeCell>::GetGatherNeighbourList
 (FLOAT rp[ndim],                      ///< Position vector
  FLOAT rsearch,                       ///< Gather search radius
  SphParticle<ndim> *sph_gen,          ///< Pointer to SPH particle array
  int Nhydro,                            ///< No. of SPH particles
  int Nneibmax,                        ///< Max. no. of neighbours
  int *neiblist)                       ///< List of neighbouring particles
{
  int Nneib = 0;                       // No. of (non-dead) neighbours
  ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sph_gen);

  debug2("[SphTree::GetGatherNeighbourList]");

  Nneib = tree->ComputeGatherNeighbourList(sphdata,rp,rsearch,Nneibmax,neiblist);

  return Nneib;
}



//=================================================================================================
//  SphTree::UpdateActiveParticleCounters
/// Loop through all leaf cells in the tree and update all active particle counters.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void SphTree<ndim,ParticleType,TreeCell>::UpdateActiveParticleCounters
 (SphParticle<ndim> * sphdata_gen,     ///< ..
  Sph<ndim> *sph)                      ///< ..
{
  ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sphdata_gen);
  tree->UpdateActiveParticleCounters(sphdata);
}



//=================================================================================================
//  SphTree::SearchBoundaryGhostParticles
/// Search domain to create any required ghost particles near any boundaries.
/// Currently only searches to create periodic or mirror ghost particles.
//=================================================================================================
template <int ndim, template <int> class ParticleType, template<int> class TreeCell>
void SphTree<ndim,ParticleType,TreeCell>::SearchBoundaryGhostParticles
 (FLOAT tghost,                        ///< Ghost particle 'lifetime'
  DomainBox<ndim> simbox,              ///< Simulation box structure
  Sph<ndim> *sph)                      ///< Sph object pointer
{
  int c;                                       // ..
  int i;                                       // Particle counter
  const FLOAT grange = ghost_range*kernrange;  // ..
  TreeCell<ndim> *cellptr;                     // ..

  // Set all relevant particle counters
  sph->Nghost         = 0;
  sph->NPeriodicGhost = 0;
  sph->Nmpighost      = 0;
  sph->Nghostmax      = sph->Nhydromax - sph->Nhydro;
  sph->Ntot           = sph->Nhydro;


  // If all boundaries are open, immediately return to main loop
  if (simbox.x_boundary_lhs == openBoundary && simbox.x_boundary_rhs == openBoundary &&
      simbox.y_boundary_lhs == openBoundary && simbox.y_boundary_rhs == openBoundary &&
      simbox.z_boundary_lhs == openBoundary && simbox.z_boundary_rhs == openBoundary) return;


  debug2("[SphTree::SearchBoundaryGhostParticles]");


  // Create ghost particles in x-dimension
  //===============================================================================================
  if ((simbox.x_boundary_lhs == openBoundary &&
       simbox.x_boundary_rhs == openBoundary) == 0) {

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
            sph->CheckXBoundaryGhostParticle(i,tghost,simbox);
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

    sph->Ntot = sph->Nhydro + sph->Nghost;
  }


  // Create ghost particles in y-dimension
  //===============================================================================================
  if (ndim >= 2 && (simbox.y_boundary_lhs == openBoundary &&
                    simbox.y_boundary_rhs == openBoundary) == 0) {

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
            sph->CheckYBoundaryGhostParticle(i,tghost,simbox);
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
    for (i=sph->Nhydro; i<sph->Ntot; i++) sph->CheckYBoundaryGhostParticle(i,tghost,simbox);

    sph->Ntot = sph->Nhydro + sph->Nghost;
  }


  // Create ghost particles in z-dimension
  //===============================================================================================
  if (ndim == 3 && (simbox.z_boundary_lhs == openBoundary &&
      simbox.z_boundary_rhs == openBoundary) == 0) {

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
            sph->CheckZBoundaryGhostParticle(i,tghost,simbox);
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
    for (i=sph->Nhydro; i<sph->Ntot; i++) sph->CheckZBoundaryGhostParticle(i,tghost,simbox);

    sph->Ntot = sph->Nhydro + sph->Nghost;
  }


  // Quit here if we've run out of memory for ghosts
  if (sph->Ntot > sph->Nhydromax) {
    string message="Not enough memory for ghost particles";
    ExceptionHandler::getIstance().raise(message);
  }

  sph->NPeriodicGhost = sph->Nghost;

  return;
}



//=================================================================================================
//  SphTree::ComputeCellMonopoleForces
/// Compute the force on particle 'parti' due to all cells obtained in the
/// gravity tree walk.  Uses only monopole moments (i.e. COM) of the cell.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void SphTree<ndim,ParticleType,TreeCell>::ComputeCellMonopoleForces
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
//  SphTree::ComputeCellQuadrupoleForces
/// Compute the force on particle 'parti' due to all cells obtained in the
/// gravity tree walk including the quadrupole moment correction term.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void SphTree<ndim,ParticleType,TreeCell>::ComputeCellQuadrupoleForces
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
//  SphTree::ComputeFastMonopoleForces
/// Compute the force on particle 'parti' due to all cells obtained in the
/// gravity tree walk.  Uses only monopole moments (i.e. COM) of the cell.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void SphTree<ndim,ParticleType,TreeCell>::ComputeFastMonopoleForces
 (int Nactive,                         ///< [in] No. of active particles
  int Ngravcell,                       ///< [in] No. of tree cells in list
  TreeCell<ndim> *gravcell,            ///< [in] List of tree cell ids
  TreeCell<ndim> &cell,                ///< [in] Current cell pointer
  ParticleType<ndim> *activepart)      ///< [inout] Active SPH particle array
{
  int cc;                              // Aux. cell counter
  int j;                               // ..
  int k;                               // Dimension counter
  FLOAT ac[ndim];                      // ..
  FLOAT dr[ndim];                      // Relative position vector
  FLOAT drsqd;                         // Distance squared
  FLOAT invdrmag;                      // 1 / distance
  FLOAT invdrsqd;                      // 1 / drsqd
  FLOAT invdr3;                        // 1 / dist^3
  FLOAT mc;                            // Mass of cell
  FLOAT q[6];                          // Local copy of quadrupole moment
  FLOAT dphi[3];                       // ..
  FLOAT cellpot;                       // ..
  FLOAT rc[ndim];                      // ..

  for (k=0; k<ndim; k++) rc[k] = cell.r[k];
  for (k=0; k<ndim; k++) ac[k] = 0.0;
  for (k=0; k<ndim; k++) dphi[k] = 0.0;
  for (k=0; k<6; k++) q[k] = 0;
  cellpot = 0.0;


  //-----------------------------------------------------------------------------------------------
  if (ndim == 3) {

    for (cc=0; cc<Ngravcell; cc++) {
#ifndef MPI_PARALLEL
      assert(cell.id != gravcell[cc].id);
#endif
      mc = gravcell[cc].m;
      for (k=0; k<ndim; k++) dr[k] = gravcell[cc].r[k] - rc[k];
      drsqd    = DotProduct(dr,dr,ndim);
      invdrsqd = 1.0/drsqd;
      invdrmag = sqrt(invdrsqd);
      invdr3   = invdrsqd*invdrmag;
      cellpot  += mc*invdrmag;
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



#ifdef MPI_PARALLEL
//=================================================================================================
//  SphTree::UpdateGravityExportList
/// Compute all local 'gather' properties of currently active particles, and then compute each
/// particle's contribution to its (active) neighbour hydro forces.  Optimises the algorithm by
/// using grid-cells to construct local neighbour lists for all particles inside the cell.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void SphTree<ndim,ParticleType,TreeCell>::UpdateGravityExportList
 (int rank,                            ///< [in] MPI rank
  int Nhydro,                            ///< [in] No. of SPH particles
  int Ntot,                            ///< [in] No. of SPH + ghost particles
  SphParticle<ndim> *sph_gen,          ///< [inout] Pointer to SPH ptcl array
  Sph<ndim> *sph,                      ///< [in] Pointer to SPH object
  Nbody<ndim> *nbody)                  ///< [in] Pointer to N-body object
{
  int cactive;                         // No. of active cells
  int cc;                              // Aux. cell counter
  int i;                               // Particle id
  int ithread;                         // OpenMP thread i.d.
  int j;                               // Aux. particle counter
  int jj;                              // Aux. particle counter
  int k;                               // Dimension counter
  int okflag;                          // Flag if h-rho iteration is valid
  int Nactive;                         // No. of active particles in current cell
  int Ngravcell=0;                     // No. of gravity cells
  int Ngravcellmax;                    // Max. size of gravity cell pointer array
  int Ngravcelltemp;                   // Aux. gravity cell counter
  FLOAT macfactor;                     // Gravity MAC factor for cell

  int *activelist;                     // List of active particles
  TreeCell<ndim> *cellptr;           // Pointer to binary tree cell
  TreeCell<ndim> **celllist;         // List of pointers to binary tree cells
  TreeCell<ndim> **gravcelllist;     // List of pointers to grav. cells
  ParticleType<ndim> *activepart;      // Local copies of active particles
  ParticleType<ndim> *neibpart;        // Local copies of neighbouring particles
  ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sph_gen);

  debug2("[GradhSphTree::UpdateDistantSphForces]");
  timing->StartTimingSection("SPH_DISTANT_FORCES",2);


  // Find list of all cells that contain active particles
  celllist = new TreeCell<ndim>*[2*tree->gtot];
  cactive = tree->ComputeActiveCellList(celllist);

  // Reset all export lists
  for (j=0; j<Nmpi; j++) {
    Ncellexport[j] = 0;
    Npartexport[j] = 0;
  }


  // Set-up all OMP threads
  //===============================================================================================
#pragma omp parallel default(none) shared(celllist,cactive,sph,sphdata,cout) \
  private(activepart,activelist,cc,cellptr,directlist,draux,drsqd,gravcelllist,hrangesqdi,i)\
  private(interactlist,ithread,j,jj,k,levelneib,macfactor,neiblist,neibpart,Nactive,Ndirect)\
  private(Ndirectaux,Ndirectmax,Ngravcell,Ngravcellmax,Ninteract,Nneib,Nneibmax,okflag,rp)
  {
#if defined _OPENMP
    ithread = omp_get_thread_num();
#else
    ithread = 0;
#endif

    Ngravcellmax = Nprunedcellmax;
    activelist = activelistbuf[ithread];
    activepart = activepartbuf[ithread];
    gravcelllist = new TreeCell<ndim>*[Ngravcellmax];


    // Loop over all active cells
    //=============================================================================================
#pragma omp for schedule(guided)
    for (cc=0; cc<cactive; cc++) {
      cellptr = celllist[cc];
      macfactor = 0.0;
      Ngravcell = 0;

      // Find list of active particles in current cell
      Nactive = tree->ComputeActiveParticleList(cellptr,sphdata,activelist);

      // Make local copies of active particles
      for (j=0; j<Nactive; j++) activepart[j] = sphdata[activelist[j]];

      // Compute average/maximum term for computing gravity MAC
      if (gravity_mac == "eigenmac") {
        for (j=0; j<Nactive; j++)
          macfactor = max(macfactor,pow(1.0/activepart[j].gpot,twothirds));
      }

      // Zero/initialise all summation variables for active particles
      for (j=0; j<Nactive; j++) {
        activepart[j].gpot = activepart[j].m*activepart[j].invh*sph->kernp->wpot(0.0);
      }


      // Loop over all distant pruned trees and compute list of cells.
      // If pruned tree is too close, record cell id for exporting
      //-------------------------------------------------------------------------------------------
      for (j=0; j<Nmpi; j++) {
        if (j == rank) continue;

        Ngravcelltemp = prunedtree[j]->ComputeDistantGravityInteractionList
          (cellptr,macfactor,Ngravcellmax,Ngravcell,gravcelllist);

        // If pruned tree is too close (flagged by -1), then record cell id
        // for exporting to other MPI processes
        if (Ngravcelltemp == -1) {
          cellexportlist[j][Ncellexport[j]++] = cellptr;
          Npartexport[j] += Nactive;
        }
        else Ngravcell = Ngravcelltemp;

      }
      //-------------------------------------------------------------------------------------------


      // Loop over all active particles in the cell
      //-------------------------------------------------------------------------------------------
      for (j=0; j<Nactive; j++) {
        i = activelist[j];

        // Compute gravitational force due to distant cells
        if (multipole == "monopole") {
          tree->ComputeCellMonopoleForces(activepart[j].gpot,activepart[j].agrav,
                                          activepart[j].r,Ngravcell,gravcelllist);
        }
        else if (multipole == "quadrupole") {
          tree->ComputeCellQuadrupoleForces(activepart[j].gpot,activepart[j].agrav,
                                            activepart[j].r,Ngravcell,gravcelllist);
        }

      }
      //-------------------------------------------------------------------------------------------


      // Compute 'fast' multipole terms here
      if (multipole == "fast_monopole") {
        tree->ComputeFastMonopoleForces(Nactive,Ngravcell,gravcelllist,cellptr,activepart);
      }

      // Add all active particles contributions to main array
      for (j=0; j<Nactive; j++) {
        i = activelist[j];
        for (k=0; k<ndim; k++) sphdata[i].a[k] = activepart[j].a[k];
        for (k=0; k<ndim; k++) sphdata[i].agrav[k] = activepart[j].agrav[k];
        for (k=0; k<ndim; k++) sphdata[i].a[k] += sphdata[i].agrav[k];
        sphdata[i].gpot = activepart[j].gpot;
      }

    }

    //=============================================================================================


    // Free-up local memory for OpenMP thread
    delete[] gravcelllist;

  }
  //===============================================================================================

  delete[] celllist;

  timing->EndTimingSection("SPH_DISTANT_FORCES");

  return;
}



//=================================================================================================
//  SphTree::UpdateHydroExportList
/// Compute all local 'gather' properties of currently active particles, and then compute each
/// particle's contribution to its (active) neighbour neighbour hydro forces.  Optimises the
/// algorithm by using grid-cells to construct local neighbour lists for all ptcls inside the cell.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void SphTree<ndim,ParticleType,TreeCell>::UpdateHydroExportList
 (int rank,                            ///< [in] MPI rank
  int Nhydro,                            ///< [in] No. of SPH particles
  int Ntot,                            ///< [in] No. of SPH + ghost particles
  SphParticle<ndim> *sph_gen,          ///< [inout] Pointer to SPH ptcl array
  Sph<ndim> *sph,                      ///< [in] Pointer to SPH object
  Nbody<ndim> *nbody)                  ///< [in] Pointer to N-body object
{
  bool overlapflag;                    // Flag if cells overlap
  int cactive;                         // No. of active cells
  int cc;                              // Aux. cell counter
  int i;                               // Particle id
  int ithread;                         // OpenMP thread i.d.
  int j;                               // Aux. particle counter
  int jj;                              // Aux. particle counter
  TreeCell<ndim> *cellptr;           // Pointer to binary tree cell
  TreeCell<ndim> **celllist;         // List of pointers to binary tree cells
  int *activelist;                     // List of active particles


  debug2("[SphTree::UpdateHydroExportList]");
  timing->StartTimingSection("MPI_HYDRO_EXPORT",2);


  // Find list of all cells that contain active particles
  celllist = new TreeCell<ndim>*[2*tree->gtot];
  cactive = tree->ComputeActiveCellList(celllist);

  ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sph_gen);

  // Reset all export lists
  for (j=0; j<Nmpi; j++) {
    Ncellexport[j] = 0;
    Npartexport[j] = 0;
  }


  // Set-up all OMP threads
  //===============================================================================================
#pragma omp parallel default(none) shared(celllist,cactive,cout) \
  private(cc,cellptr,i,ithread,j,jj,overlapflag)
  {
#if defined _OPENMP
    ithread = omp_get_thread_num();
#else
    ithread = 0;
#endif
    activelist = activelistbuf[ithread];

    // Loop over all active cells
    //=============================================================================================
#pragma omp for schedule(guided)
    for (cc=0; cc<cactive; cc++) {
      cellptr = celllist[cc];

      // Loop over all distant pruned trees and compute list of cells.
      // If pruned tree is too close, record cell id for exporting
      //-------------------------------------------------------------------------------------------
      for (j=0; j<Nmpi; j++) {
        if (j == rank) continue;

        overlapflag = prunedtree[j]->ComputeHydroTreeCellOverlap(cellptr);

        // If pruned tree is too close (flagged by -1), then record cell id
        // for exporting to other MPI processes
        if (overlapflag) {
          cellexportlist[j][Ncellexport[j]++] = cellptr;
          const int Nactive = tree->ComputeActiveParticleList(cellptr,sphdata,activelist);
          Npartexport[j] += Nactive;
          assert(Ncellexport[j] <= tree->gmax);
          cout << "Found overlap : " << cc << "   " << j << "   " << Nactive << "    " << cellptr->N
               << "    " << Npartexport[j] << "    " << Ncellexport[j] << endl;
        }

      }
      //-------------------------------------------------------------------------------------------

    }
    //=============================================================================================

    PrintArray("Ncellexport : ",Nmpi,Ncellexport);
    PrintArray("Npartexport : ",Nmpi,Npartexport);


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
//  SphTree::BuildPrunedTree
/// Constructs a pruned version of the local tree ready to be exported to other MPI processes.
/// Copies all levels up to and including 'pruning_level'.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void SphTree<ndim,ParticleType,TreeCell>::BuildPrunedTree
 (int pruning_level,                   ///< ..
  int rank)                            ///< ..
{
  int c;                               // Cell counter
  int cnew;                            // New i.d. of copied cell in pruned tree
  int cnext;                           // i.d. of next cell in pruned tree
  int c1;                              // i.d. of first child cell
  int c2;                              // i.d. of second child cell
  int i;                               // Particle counter
  int k;                               // Dimension counter
  int l;                               // ..

  debug2("[SphTree::BuildPrunedTree]");
  timing->StartTimingSection("BUILD_PRUNED_TREE",2);

  cnew = 0;
  Nprunedcellmax = 0;
  assert(pruning_level < tree->ltot);


  // Set level at which tree will be pruned (for all trees)
  //-----------------------------------------------------------------------------------------------
  for (i=0; i<Nmpi; i++) {
    prunedtree[i]->ltot_old = prunedtree[i]->ltot;
    prunedtree[i]->ltot     = pruning_level;
    prunedtree[i]->gmax     = pow(2,pruning_level);
    prunedtree[i]->Ncellmax = 2*prunedtree[i]->gmax - 1;
    prunedtree[i]->Ncell    = 2*prunedtree[i]->gmax - 1;
    prunedtree[i]->Ntotmax  = prunedtree[i]->gmax;
    prunedtree[i]->Ntot     = prunedtree[i]->gmax;

    // Allocate (or reallocate if needed) all tree memory
    prunedtree[i]->AllocateTreeMemory();
    Nprunedcellmax += prunedtree[i]->Ncellmax;

    // Sanity-check to ensure pruned tree is not deeper than real local tree
    assert(prunedtree[i]->gmax == pow(2,pruning_level));
    assert(prunedtree[i]->ltot == pruning_level);
    if (pruning_level > tree->ltot) {
      string message = "Invalid pruned tree; pruning_level > ltot";
      ExceptionHandler::getIstance().raise(message);
    }

  }
  //-----------------------------------------------------------------------------------------------


  // Now walk through main tree cell-by-cell and copy all important data to pruned tree cells
  //-----------------------------------------------------------------------------------------------
  for (c=0; c<tree->Ncell; c++) {

    // If cell is on a lower level, skip over
    if (tree->celldata[c].level > pruning_level) continue;

    // Otherwise, record all data from cell, except for cell links which
    // are maintained to ensure a valid tree
    c1 = prunedtree[rank]->celldata[cnew].c1;
    c2 = prunedtree[rank]->celldata[cnew].c2;
    cnext = prunedtree[rank]->celldata[cnew].cnext;
    l = prunedtree[rank]->celldata[cnew].level;

    prunedtree[rank]->celldata[cnew] = tree->celldata[c];
    prunedtree[rank]->celldata[cnew].c1 = c1;
    prunedtree[rank]->celldata[cnew].c2 = c2;
    prunedtree[rank]->celldata[cnew].cnext = cnext;

    if (c < pruning_level) assert (l == cnew);
    assert(cnext >= 0);
    assert(prunedtree[rank]->celldata[cnew].level <= pruning_level);
    if (prunedtree[rank]->celldata[cnew].level < pruning_level) assert(c1 == cnew + 1);
    if (c1 == -1) assert(prunedtree[rank]->celldata[cnew].level == pruning_level);

    cnew++;

  }
  //-----------------------------------------------------------------------------------------------

  //cout << "Pruned tree size : " << prunedtree[rank]->Ncell << "    "
  //     << prunedtree[rank]->ltot << "    " << pruning_level << endl;
  /*for (c=0; c<prunedtree[rank]->Ncell; c++) {
    cout << "bb[" << c << "] : " << prunedtree[rank]->celldata[c].bbmin[0]
         << "    " << prunedtree[rank]->celldata[c].bbmax[0]
         << "    " << prunedtree[rank]->celldata[c].hboxmin[0]
         << "    " << prunedtree[rank]->celldata[c].hboxmax[0]
         << "    N : " << prunedtree[rank]->celldata[c].N << endl;
  }*/

  //cin >> c;
  MPI_Barrier(MPI_COMM_WORLD);

  return;
}



//=================================================================================================
//  SphTree::BuildMpiGhostTree
/// Main routine to control how the tree is built, re-stocked and interpolated
/// during each timestep.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void SphTree<ndim,ParticleType,TreeCell>::BuildMpiGhostTree
 (bool rebuild_tree,                   ///< Flag to rebuild tree
  int n,                               ///< Integer time
  int ntreebuildstep,                  ///< Tree build frequency
  int ntreestockstep,                  ///< Tree stocking frequency
  int Npart,                           ///< No. of particles
  int Npartmax,                        ///< Max. no. of particles
  SphParticle<ndim> *sph_gen,          ///< Particle data array
  Sph<ndim> *sph,                      ///< Pointer to SPH object
  FLOAT timestep)                      ///< Smallest physical timestep
{
  int i;                               // Particle counter
  int k;                               // Dimension counter
  ParticleType<ndim> *sphdata = static_cast<ParticleType<ndim>* > (sph_gen);


  // If no MPI ghosts exist, do not build tree
  //if (sph->Nmpighost == 0) return;

  debug2("[SphTree::BuildGhostTree]");
  timing->StartTimingSection("BUILD_MPIGHOST_TREE",2);

  // Activate nested parallelism for tree building routines
#ifdef _OPENMP
  omp_set_nested(1);
#endif

  cout << "BUILDING TREE WITH " << sph->Nmpighost << " MPI GHOSTS!!" << endl;

  // For tree rebuild steps
  //-----------------------------------------------------------------------------------------------
  if (n%ntreebuildstep == 0 || rebuild_tree) {

    mpighosttree->Ntot       = sph->Nmpighost;
    mpighosttree->Ntotmaxold = mpighosttree->Ntotmax;
    mpighosttree->Ntotmax    = max(mpighosttree->Ntotmax,mpighosttree->Ntot);
    mpighosttree->Ntotmax    = max(mpighosttree->Ntotmax,sph->Nhydromax);
    mpighosttree->BuildTree(sph->Nhydro + sph->NPeriodicGhost,
                            sph->Nhydro + sph->NPeriodicGhost +sph->Nmpighost - 1,
                            mpighosttree->Ntot, mpighosttree->Ntotmax, sphdata, timestep);

  }

  // Else stock the tree
  //-----------------------------------------------------------------------------------------------
  else if (n%ntreestockstep == 0) {

    mpighosttree->StockTree(mpighosttree->celldata[0],sphdata);

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
//  SphTree::SearchMpiGhostParticles
/// ...
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
int SphTree<ndim,ParticleType,TreeCell>::SearchMpiGhostParticles
 (const FLOAT tghost,                  ///< [in] Expected ghost life-time
  const Box<ndim> &mpibox,             ///< [in] Bounding box of MPI domain
  Sph<ndim> *sph,                      ///< [in] Pointer to SPH object
  vector<int> &export_list)            ///< [out] List of particle ids
{
  int c;                               // Cell counter
  int i;                               // ..
  int k;                               // ..
  int Nexport = 0;                     // No. of MPI ghosts to export
  FLOAT scattermin[ndim];              // ..
  FLOAT scattermax[ndim];              // ..
  TreeCell<ndim> *cellptr;              // ..
  const FLOAT grange = ghost_range*kernrange;


  // Start from root-cell
  c = 0;

  //-----------------------------------------------------------------------------------------------
  while (c < tree->Ncell) {
    cellptr = &(tree->celldata[c]);

    // Construct maximum cell bounding box depending on particle velocities
    for (k=0; k<ndim; k++) {
      scattermin[k] = cellptr->bbmin[k] + min(0.0,cellptr->v[k]*tghost) - grange*cellptr->hmax;
      scattermax[k] = cellptr->bbmax[k] + max(0.0,cellptr->v[k]*tghost) + grange*cellptr->hmax;
    }


    // If maximum cell scatter box overlaps MPI domain, open cell
    //---------------------------------------------------------------------------------------------
    if (BoxOverlap(ndim,scattermin,scattermax,mpibox.boxmin,mpibox.boxmax)) {

      // If not a leaf-cell, then open cell to first child cell
      if (cellptr->level != tree->ltot)
        c++;

      else if (cellptr->N == 0)
        c = cellptr->cnext;

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
    else
      c = cellptr->cnext;

  }
  //-----------------------------------------------------------------------------------------------


  // Start from root-cell of ghost-tree
  c = 0;

  //-----------------------------------------------------------------------------------------------
  while (c < ghosttree->Ncell) {
    cellptr = &(ghosttree->celldata[c]);

    // Construct maximum cell bounding box depending on particle velocities
    for (k=0; k<ndim; k++) {
      scattermin[k] = cellptr->bbmin[k] + min(0.0,cellptr->v[k]*tghost) - grange*cellptr->hmax;
      scattermax[k] = cellptr->bbmax[k] + max(0.0,cellptr->v[k]*tghost) + grange*cellptr->hmax;
    }

    // If maximum cell scatter box overlaps MPI domain, open cell
    //---------------------------------------------------------------------------------------------
    if (BoxOverlap(ndim,scattermin,scattermax,mpibox.boxmin,mpibox.boxmax)) {

      // If not a leaf-cell, then open cell to first child cell
      if (cellptr->level != ghosttree->ltot)
        c++;

      else if (cellptr->N == 0)
        c = cellptr->cnext;

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
    else
      c = cellptr->cnext;

  }
  //-----------------------------------------------------------------------------------------------


  return Nexport;
}



//=================================================================================================
//  SphTree::SearchHydroExportParticles
/// ...
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
int SphTree<ndim,ParticleType,TreeCell>::SearchHydroExportParticles
 (const Box<ndim> &mpibox,                 ///< [in] Bounding box of MPI domain
  Sph<ndim> *sph,                          ///< [in] Pointer to SPH object
  vector<TreeCell<ndim> *> &cell_list)   ///< [out] List of particle ids
{
  int c;                                   // Cell counter
  int i;                                   // ..
  int k;                                   // ..
  int Nexport = 0;                         // No. of MPI ghosts to export
  FLOAT scattermin[ndim];                  // ..
  FLOAT scattermax[ndim];                  // ..
  TreeCell<ndim> *cellptr;                  // ..
  const FLOAT grange = ghost_range*kernrange;
  ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sph->GetSphParticleArray());


  // Start from root-cell
  c = 0;

  //-----------------------------------------------------------------------------------------------
  while (c < tree->Ncell) {
    cellptr = &(tree->celldata[c]);

    // Construct maximum cell bounding box depending on particle velocities
    for (k=0; k<ndim; k++) {
      scattermin[k] = cellptr->bbmin[k] - grange*cellptr->hmax;
      scattermax[k] = cellptr->bbmax[k] + grange*cellptr->hmax;
    }


    // If maximum cell scatter box overlaps MPI domain, open cell
    //---------------------------------------------------------------------------------------------
    if (BoxOverlap(ndim,scattermin,scattermax,mpibox.boxmin,mpibox.boxmax)) {

      // If not a leaf-cell, then open cell to first child cell
      if (cellptr->level != tree->ltot)
        c++;

      else if (cellptr->N == 0)
        c = cellptr->cnext;

      // If leaf-cell and active, add the cell to the list of cells being exported
      else if (cellptr->level == tree->ltot) {
        if (cellptr->Nactive>0) {
          Nexport += cellptr->Nactive;
          cell_list.push_back(cellptr);
        }
        c = cellptr->cnext;
      }
    }

    // If not in range, then open next cell
    //---------------------------------------------------------------------------------------------
    else
      c = cellptr->cnext;

  }
  //-----------------------------------------------------------------------------------------------

  return Nexport;

}



//=================================================================================================
//  SphTree::FindMpiTransferParticles
/// ..
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void SphTree<ndim,ParticleType,TreeCell>::FindMpiTransferParticles
 (Sph<ndim>* sph,                            ///< [in] Pointer to sph class
  vector<vector<int> >& particles_to_export, ///< [inout] Vector that for each
                                             ///< node gives the list of particles to export
  vector<int>& all_particles_to_export,      ///< [inout] Vector containing all the particles
                                             ///<         that will be exported by this processor
  const vector<int>& potential_nodes,        ///< [in] Vector containing the potential nodes we
                                             ///<      might be sending particles to
  MpiNode<ndim>* mpinodes)                   ///< [in] Array of other mpi nodes
{
  int c;                                     // ..
  int i;                                     // ..
  int inode;                                 // ..
  int node_number;                           // ..
  TreeCell<ndim> *cellptr;                   // ..
  ParticleType<ndim> *sphdata = static_cast<ParticleType<ndim>* > (sph->GetSphParticleArray());


  // Loop over potential domains and walk the tree for each bounding box
  //-----------------------------------------------------------------------------------------------
  for (inode=0; inode<potential_nodes.size(); inode++) {

    node_number = potential_nodes[inode];
    Box<ndim>& nodebox = mpinodes[node_number].domain;

    // Start from root-cell
    c = 0;

    //---------------------------------------------------------------------------------------------
    while (c < Ncell) {
      cellptr = &(tree->celldata[c]);

      // If maximum cell scatter box overlaps MPI domain, open cell
      //-------------------------------------------------------------------------------------------
      if (BoxOverlap(ndim,cellptr->bbmin,cellptr->bbmax,nodebox.boxmin,nodebox.boxmax)) {

        // If not a leaf-cell, then open cell to first child cell
        if (cellptr->level != tree->ltot)
          c++;

        else if (cellptr->N == 0)
          c = cellptr->cnext;

        // If leaf-cell, check through particles in turn to find ghosts and
        // add to list to be exported
        else if (cellptr->level == tree->ltot) {
          i = cellptr->ifirst;
          while (i != -1) {
            if (ParticleInBox(sphdata[i],mpinodes[node_number].domain)) {
              particles_to_export[node_number].push_back(i);
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
      else
        c = cellptr->cnext;

    }
    //---------------------------------------------------------------------------------------------

  }
  //-----------------------------------------------------------------------------------------------

  return;
}



//=================================================================================================
//  SphTree::FindLoadBalancingDivision
/// Get the array with the information that needs to be exported to the given processor (NB: Nproc
/// is ignored at the moment, as we always need to export all particles to the other processors).
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
FLOAT SphTree<ndim,ParticleType,TreeCell>::FindLoadBalancingDivision
 (int k_divide,                        ///< Dimension of cell division
  FLOAT r_old,                         ///< Old position of cell division
  Box<ndim> &box)                      ///< Parent box to divide
{
  int i;
  int k;
  FLOAT r_divide = r_old;              // Cell division location
  FLOAT r_max = box.boxmax[k_divide];  // Max. for bisection iteration of division
  FLOAT r_min = box.boxmin[k_divide];  // Min. for bisection iteration of division
  FLOAT workleft;                      // Work computed on LHS of division
  FLOAT workright;                     // Work computed on RHS of division
  FLOAT workfrac;                      // Fraction of work on LHS
  FLOAT worktol = 0.001;               // Work balance tolerance for iteration
  Box<ndim> boxleft = box;             // LHS box
  Box<ndim> boxright = box;            // RHS box


  // Find the work-balance position through bisection iteration
  //-----------------------------------------------------------------------------------------------
  do {
    workleft = 0.0;
    workright = 0.0;
    boxleft.boxmax[k_divide] = r_divide;
    boxright.boxmin[k_divide] = r_divide;


    // Compute work included in left-hand side
    for (i=0; i<Nmpi; i++) {
      workleft += 0.0;
    }


    // Compute work included in right-hand side
    for (i=0; i<Nmpi; i++) {
      workright += 0.0;
    }


    // If fraction of work on either side of division is too inbalanced,
    // calculate new position of division and perform another iteration.
    // Otherwise exit iteration loop.
    workfrac = workleft / (workleft + workright);
    if (workfrac < 0.5 - worktol) r_min = r_divide;
    else if (workfrac > 0.5 + worktol) r_max = r_divide;
    else break;

    r_divide = 0.5*(r_min + r_max);

  } while (fabs(workfrac - 0.5) < worktol);
  //-----------------------------------------------------------------------------------------------

  return r_divide;
}



//=================================================================================================
//  SphTree::GetExportInfo
/// Get the array with the information that needs to be exported to the given processor
/// (NB: Nproc is ignored at the moment, as we must export all ptcls to all other processors)
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
int SphTree<ndim,ParticleType,TreeCell>::GetExportInfo
 (int Nproc,                           ///< [in] No. of processor we want to send data to
  Sph<ndim>* sph,                      ///< [in] Pointer to sph object
  vector<char >& send_buffer,          ///< [inout] Vector where the ptcls to export will be stored
  MpiNode<ndim>& mpinode,
  int rank,
  int Nmpi)                            ///< [in] Array with information for the other mpi nodes
{
  ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sph->GetSphParticleArray() );
  const bool first_proc = (Nproc==0) || (rank==0 && Nproc==1);
  const bool hydro_only = !sph->self_gravity && sph->hydro_forces;
  int Nactive=0, cactive;

  assert(tree->Nimportedcell==0);

  // Get active cells and their number (so that we know how much memory to allocate)
//  vector<TreeCell<ndim>*> celllist;
//  celllist.reserve(tree->gtot);
//  if (hydro_only) {
//    Nactive += SearchHydroExportParticles(mpinode.domain,sph,celllist);
//    cactive = celllist.size();
//  }
//  else {
//
//    for (int i=0; i<sph->Nhydro; i++) {
//      if (sphdata[i].active)
//        Nactive++;
//    }
//
//    celllist.resize(tree->gtot);
//    cactive = tree->ComputeActiveCellList(&celllist[0]);
//  }

  TreeCell<ndim>** celllist = cellexportlist[Nproc];
  cactive = Ncellexport[Nproc];
  Nactive = Npartexport[Nproc];


  // Work out size of the information we are sending
  //Header consists of number of particles and number of cells
  const int size_header = 2*sizeof(int);
  const int size_particles = Nactive*sizeof(ParticleType<ndim>);
  const int size_cells = cactive*sizeof(TreeCell<ndim>);
  const int old_size = send_buffer.size();
  send_buffer.resize(size_particles+size_cells+size_header+old_size);

  //Write the header
  copy(&send_buffer[old_size],&Nactive);
  copy(&send_buffer[old_size+sizeof(int)],&cactive);

  //Clear the array needed for bookkeeping (which active particles we sent to which processor)
  vector<int>& ids_active_particles = ids_sent_particles[Nproc];
  ids_active_particles.clear(); ids_active_particles.reserve(Nactive);

  //Copy cells to export inside array
  int offset = size_header+old_size;
  int activelist[Nleafmax];
  int exported_particles = 0;
  for (int i=0; i<cactive; i++) {
    copy(&send_buffer[offset],celllist[i]);
    const int Nactive_cell = tree->ComputeActiveParticleList(celllist[i],sphdata,activelist);
    // Update the ifirst and ilast pointers in the cell
    TreeCell<ndim>* exported_cell = reinterpret_cast<TreeCell<ndim>*> (&send_buffer[offset]);
    exported_cell->ifirst = exported_particles;
    exported_cell->ilast = exported_particles+Nactive_cell-1;
    offset += sizeof(TreeCell<ndim>);
    // Copy active particles
    for (int iparticle=0; iparticle<Nactive_cell; iparticle++) {
      ids_active_particles.push_back(activelist[iparticle]);
      copy(&send_buffer[offset],&sphdata[activelist[iparticle]]);
      offset += sizeof(ParticleType<ndim>);
    }
    exported_particles += Nactive_cell;
  }
  assert(exported_particles == Nactive);
  assert(offset == send_buffer.size());

  return size_particles+size_cells+size_header;

}



//=================================================================================================
//  SphTree::UnpackExported
/// Unpack the information exported from the other processors, contaning the particles
/// that were exported and
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void SphTree<ndim,ParticleType,TreeCell>::UnpackExported
 (vector<char >& received_array,
  vector<int>& Nbytes_from_proc,
  Sph<ndim>* sph)
{
  int offset = 0;

  assert(sph->NImportedParticles==0);
  tree->Nimportedcell=0;
  tree->Ncelltot=tree->Ncell;

  N_imported_part_per_proc.resize(Nbytes_from_proc.size());

  ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sph->GetSphParticleArray() );

  //-----------------------------------------------------------------------------------------------
  for (int Nproc = 0; Nproc<Nbytes_from_proc.size(); Nproc++) {

    int N_received_bytes = Nbytes_from_proc[Nproc];
    int N_received_particles; int N_received_cells;

    if (N_received_bytes == 0) {
      N_imported_part_per_proc[Nproc]=0;
      continue;
    }

    copy(&N_received_particles,&received_array[offset]);
    N_imported_part_per_proc[Nproc]=N_received_particles;
    copy(&N_received_cells,&received_array[offset+sizeof(int)]);

    //Ensure there is enough memory
    if (sph->Ntot + N_received_particles > sph->Nhydromax) {
      ExceptionHandler::getIstance().raise("Error while receiving imported particles: not enough memory!");
    }
    if (tree->Ncelltot + N_received_cells > tree->Ncellmax) {
      ExceptionHandler::getIstance().raise("Error while receiving imported cells: not enough memory!");
    }


    //Copy received particles inside SPH main arrays and received cells inside the tree array
    //Also update the linked list
    int particle_index = sph->Ntot;
    offset += 2*sizeof(int);
    for (int icell=0; icell<N_received_cells; icell++) {
      TreeCell<ndim>& dest_cell = tree->celldata[icell+tree->Ncelltot];
      copy(&dest_cell,&received_array[offset]);
      offset += sizeof(TreeCell<ndim>);
      // Offset the ifirst and ilast pointers
      dest_cell.ifirst += sph->Ntot;
      dest_cell.ilast += sph->Ntot;
      // Now copy the received particles inside the SPH main arrays
      for (int iparticle=0;iparticle<dest_cell.Nactive;iparticle++) {
        copy( &sphdata[particle_index] , &received_array[offset]);
        tree->inext[particle_index] = particle_index + 1;
        particle_index++; offset+=sizeof(ParticleType<ndim>);
      }

    }
    //---------------------------------------------------------------------------------------------

    //Update the SPH counters
    sph->Ntot += N_received_particles;
    sph->NImportedParticles += N_received_particles;

    //Update the tree counters
    tree->Nimportedcell += N_received_cells;
    tree->Ncelltot      += N_received_cells;
    tree->Ntot          = sph->Ntot;

  }

  assert (offset == std::accumulate(Nbytes_from_proc.begin(), Nbytes_from_proc.end(),0));


}


//=================================================================================================
//  SphTree::GetBackExportInfo
/// Return the data to transmit back to the other processors (particle acceleration etc.)
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void SphTree<ndim,ParticleType,TreeCell>::GetBackExportInfo
 (vector<char >& send_buffer,              ///< [inout] These arrays will be overwritten with the information to send
  vector<int>& Nbytes_from_proc,  ///< ..
  vector<int>& Nbytes_to_proc,        ///< ..
  Sph<ndim>* sph,                          ///< [in] Pointer to the SPH object
  int rank)                                ///< ..
{
  int InitialNImportedParticles = sph->NImportedParticles;

  //loop over the processors, removing particles as we go
  int removed_particles=0;
  send_buffer.resize(sph->NImportedParticles * sizeof(ParticleType<ndim>));

  //-----------------------------------------------------------------------------------------------
  for (int Nproc=0 ; Nproc < N_imported_part_per_proc.size(); Nproc++ ) {

    const int N_received_particles = N_imported_part_per_proc[Nproc];

//    //Copy the accelerations and gravitational potential of the particles
//    ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sph->GetSphParticleArray() );
//    int j=0;
//    for (int i=sph->Nhydro - N_received_particles; i<sph->Nhydro; i++) {
//      for (int k=0; k<ndim; k++)
//        send_buffer[removed_particles+j].a[k] = sphdata[i].a[k];
//      send_buffer[removed_particles+j].gpot = sphdata[i].gpot;
//      j++;
//    }

    //Copy the particles inside the send buffer
    ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sph->GetSphParticleArray() );
    int j=0;
    const int start_index = sph->Nhydro + sph->Nghost + removed_particles;
    for (int i=start_index; i<start_index + N_received_particles; i++) {
      copy (&send_buffer[(removed_particles+j)*sizeof(ParticleType<ndim>)],&sphdata[i]);
      j++;
    }

    assert(j==N_received_particles);

    removed_particles += j;

    //Decrease the particle counter
    sph->Ntot -= N_received_particles;
    sph->NImportedParticles -= N_received_particles;

    //Update the information about how much data we are sending
    Nbytes_from_proc[Nproc] = N_received_particles*sizeof(ParticleType<ndim>);

    //Update the information with how much data we are receiving
    Nbytes_to_proc[Nproc] = ids_sent_particles[Nproc].size()*sizeof(ParticleType<ndim>);

  }
  //-----------------------------------------------------------------------------------------------

  tree->Ncelltot      = tree->Ncell;
  tree->Nimportedcell = 0;
  tree->Ntot          -= InitialNImportedParticles;

  assert(sph->NImportedParticles == 0);
  assert(sph->Ntot == sph->Nhydro + sph->Nghost);
  assert(send_buffer.size() == removed_particles*sizeof(ParticleType<ndim>));

}



//=================================================================================================
//  SphTree::UnpackReturnedExportInfo
/// Unpack the data that was returned by the other processors, summing the accelerations to the particles
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void SphTree<ndim,ParticleType,TreeCell>::UnpackReturnedExportInfo
 (vector<char >& received_information,
  vector<int>& recv_displs,
  Sph<ndim>* sph,
  int rank)
{

  ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sph->GetSphParticleArray() );

  //For each processor, sum up the received quantities for each particle
  for(int Nproc=0; Nproc<recv_displs.size(); Nproc++ ) {

    if (rank==Nproc)
      continue;

    const vector<int>& ids_active_particles = ids_sent_particles[Nproc];

    for (int i=0; i< ids_active_particles.size(); i++ ) {
    const int j = ids_active_particles[i];

      ParticleType<ndim>* received_particle = reinterpret_cast<ParticleType<ndim>*>
      (&received_information[i * sizeof(ParticleType<ndim>) + recv_displs[Nproc] ]);

      assert(sphdata[j].iorig == received_particle->iorig);

      for (int k=0; k<ndim; k++) {
        sphdata[j].a[k] += received_particle->a[k];
        sphdata[j].agrav[k] += received_particle->agrav[k];
      }
      sphdata[j].gpot += received_particle->gpot;
      sphdata[j].gpe += received_particle->gpe;
      sphdata[j].dudt += received_particle->dudt;
      sphdata[j].div_v += received_particle->div_v;
      sphdata[j].levelneib = max(sphdata[j].levelneib, received_particle->levelneib);

    }

  }

}



//=================================================================================================
//  SphTree::CommunicatePrunedTrees
/// ..
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void SphTree<ndim,ParticleType,TreeCell>::CommunicatePrunedTrees
 (vector<int>& my_matches,
  int rank)
{

  MPI_Barrier(MPI_COMM_WORLD);

  //-----------------------------------------------------------------------------------------------
  for (int iturn=0; iturn<my_matches.size(); iturn++) {

    //See who we need to communicate with
    int inode = my_matches[iturn];
    cout << "inode : " << rank << "   " << inode << "    " << iturn << endl;

    //Decide if sending or receiving first
    bool send_turn;
    if (rank<inode) {
      send_turn = true;
    }
    else {
      send_turn = false;
    }

    //Do the actual communication
    for (int i=0; i<2; i++) {
      if (send_turn) {
        Tree<ndim, ParticleType, TreeCell>* treeptr = prunedtree[rank];
        MPI_Send(treeptr->celldata,treeptr->Ncell*sizeof(TreeCell<ndim>),
                 MPI_CHAR,inode,3,MPI_COMM_WORLD);
        send_turn=false;
      }
      else {
        Tree<ndim, ParticleType, TreeCell>* treeptr = prunedtree[inode];
        MPI_Status status;
        MPI_Recv(treeptr->celldata,treeptr->Ncell*sizeof(TreeCell<ndim>),
                 MPI_CHAR,inode,3,MPI_COMM_WORLD,&status);
        send_turn=true;
      }
    }

  }
  //----------------------------------------------------------------------------------------------

  MPI_Barrier(MPI_COMM_WORLD);
  for (int i=0; i<Nmpi; i++) {
    cout << "Writing pruned tree " << i << " for process " << rank << endl;
    cout << "Ncell : " << prunedtree[i]->Ncell << endl;
    cout << "r : " << prunedtree[i]->celldata[0].r[0]
         << "   box : " << prunedtree[i]->celldata[0].bbmin[0]
         << "   " << prunedtree[i]->celldata[0].bbmax[0] << endl;
  }

}
#endif



#if defined(VERIFY_ALL)
//=================================================================================================
//  SphNeighbourSearch::CheckValidNeighbourList
/// Checks that the neighbour list generated by the grid is valid in that it
/// (i) does include all true neighbours, and
/// (ii) all true neigbours are only included once and once only.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void SphTree<ndim,ParticleType,TreeCell>::CheckValidNeighbourList
(int i,                             ///< [in] Particle i.d.
 int Ntot,                          ///< [in] Total no. of particles
 int Nneib,                         ///< [in] No. of potential neighbours
 int *neiblist,                     ///< [in] List of potential neighbour i.d.s
 ParticleType<ndim> *partdata,      ///< [in] Array of particle data
 string neibtype)                   ///< [in] Neighbour search type
{
  bool invalid_flag = false;        // ..
  int count;                        // Valid neighbour counter
  int j;                            // Neighbour particle counter
  int k;                            // Dimension counter
  int Ntrueneib = 0;                // No. of 'true' neighbours
  int *trueneiblist;                // List of true neighbour ids
  FLOAT drsqd;                      // Distance squared
  FLOAT dr[ndim];                   // Relative position vector

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
           << partdata[trueneiblist[j]].r[0] << endl;
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
    string message = "Problem with neighbour lists in grid search";
    ExceptionHandler::getIstance().raise(message);
  }


  delete[] trueneiblist;

  return;
}
#endif


template class SphTree<1,SphParticle,KDTreeCell>;
template class SphTree<2,SphParticle,KDTreeCell>;
template class SphTree<3,SphParticle,KDTreeCell>;
template class SphTree<1,SphParticle,OctTreeCell>;
template class SphTree<2,SphParticle,OctTreeCell>;
template class SphTree<3,SphParticle,OctTreeCell>;

template class SphTree<1,GradhSphParticle,KDTreeCell>;
template class SphTree<2,GradhSphParticle,KDTreeCell>;
template class SphTree<3,GradhSphParticle,KDTreeCell>;
template class SphTree<1,GradhSphParticle,OctTreeCell>;
template class SphTree<2,GradhSphParticle,OctTreeCell>;
template class SphTree<3,GradhSphParticle,OctTreeCell>;

template class SphTree<1,SM2012SphParticle,KDTreeCell>;
template class SphTree<2,SM2012SphParticle,KDTreeCell>;
template class SphTree<3,SM2012SphParticle,KDTreeCell>;
template class SphTree<1,SM2012SphParticle,OctTreeCell>;
template class SphTree<2,SM2012SphParticle,OctTreeCell>;
template class SphTree<3,SM2012SphParticle,OctTreeCell>;

template class SphTree<1,GodunovSphParticle,KDTreeCell>;
template class SphTree<2,GodunovSphParticle,KDTreeCell>;
template class SphTree<3,GodunovSphParticle,KDTreeCell>;
template class SphTree<1,GodunovSphParticle,OctTreeCell>;
template class SphTree<2,GodunovSphParticle,OctTreeCell>;
template class SphTree<3,GodunovSphParticle,OctTreeCell>;
