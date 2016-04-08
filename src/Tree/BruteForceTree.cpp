//=================================================================================================
//  BruteForceTree.cpp
//  Contains all functions for building, stocking and walking for the
//  BruteForceTree for neighbour finding.
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
#include <string>
#include <math.h>
#include "Precision.h"
#include "Exception.h"
#include "DomainBox.h"
#include "Parameters.h"
#include "InlineFuncs.h"
#include "Particle.h"
#include "Sph.h"
#include "BruteForceTree.h"
#include "Debug.h"
#if defined _OPENMP
#include <omp.h>
#endif
using namespace std;




//=================================================================================================
//  BruteForceTree::BruteForceTree
/// BruteForceTree constructor.  Initialises various variables.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
BruteForceTree<ndim,ParticleType,TreeCell>::BruteForceTree(int Nleafmaxaux, FLOAT thetamaxsqdaux,
                                           	   	   	   	   FLOAT kernrangeaux, FLOAT macerroraux,
                                           	   	   	   	   string gravity_mac_aux, string multipole_aux,
                                           	   	   	   	   const DomainBox<ndim>& domain,
                                           	   	   	   	   const ParticleTypeRegister& reg):
  Tree<ndim,ParticleType,TreeCell>(Nleafmaxaux, thetamaxsqdaux, kernrangeaux,
                                   macerroraux, gravity_mac_aux, multipole_aux, domain, reg)
{
  allocated_tree = false;
  gmax           = 0;
  gtot           = 0;
  ifirst         = -1;
  ilast          = -1;
  lmax           = 0;
  ltot           = 0;
  ltot_old       = -1;
  Ncell          = 0;
  Ncellmax       = 0;
  Ncellmaxold    = 0;
  Ntot           = 0;
  Ntotmax        = 0;
  Ntotmaxold     = 0;
  Ntotold        = -1;
  hmax           = 0.0;
#if defined _OPENMP
  Nthreads       = omp_get_max_threads();
#else
  Nthreads       = 1;
#endif
#if defined MPI_PARALLEL
  Ncelltot       = 0;
  Nimportedcell  = 0;
#endif
}

//=================================================================================================
//  BruteForceTree::~BruteForceTree
/// BruteForceTree destructor.  Deallocates tree memory upon object destruction.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
BruteForceTree<ndim,ParticleType,TreeCell>::~BruteForceTree()
{
  if (allocated_tree) DeallocateTreeMemory();
}



//=================================================================================================
//  BruteForceTree::AllocateTreeMemory
/// Allocate memory for BruteForce-tree as requested.  If more memory is required
/// than currently allocated, tree is deallocated and reallocated here.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void BruteForceTree<ndim,ParticleType,TreeCell>::AllocateTreeMemory(void)
{
  debug2("[BruteForceTree::AllocateTreeMemory]");

  if (!allocated_tree || Ntotmax > Ntotmaxold || Ntot > Ntotmax || Ncellmax > Ncellmaxold) {
    if (allocated_tree) DeallocateTreeMemory();
    Ntotmax     = max(Ntotmax, Ntot);
    Ntotmaxold  = Ntotmax;
    Ncellmaxold = Ncellmax;

    g2c      = new int[gmax];
    ids      = new int[Ntotmax];
    inext    = new int[Ntotmax];
    celldata = new struct TreeCell<ndim>[Ncellmax];

    allocated_tree = true;
  }


  return;
}



//=================================================================================================
//  BruteForceTree::DeallocateTreeMemory
/// Deallocates all BruteForce-tree memory
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void BruteForceTree<ndim,ParticleType,TreeCell>::DeallocateTreeMemory(void)
{
  debug2("[BruteForceTree::DeallocateTreeMemory]");

  if (allocated_tree) {
    delete[] celldata;
    delete[] inext;
    delete[] ids;
    delete[] g2c;
    allocated_tree = false;
  }

  return;
}


//=================================================================================================
//  BruteForceTree::BuildTree
/// Call all routines to build/re-build the BruteForce-tree on the local node.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void BruteForceTree<ndim,ParticleType,TreeCell>::BuildTree
 (const int _ifirst,                   ///< i.d. of first particle
  const int _ilast,                    ///< i.d. of last particle
  const int Npart,                     ///< No. of particles
  const int Npartmax,                  ///< Max. no. of particles
  const FLOAT timestep,                ///< Smallest physical timestep
  ParticleType<ndim> *partdata)        ///< Particle data array
{
  int i;                               // Particle counter
  int k;                               // Dimension counter
  FLOAT bbmin[ndim];                   // Minimum extent of local bounding box
  FLOAT bbmax[ndim];                   // Maximum extent of local bounding box

  debug2("[BruteForceTree::BuildTree]");
  //timing->StartTimingSection("BUILD_TREE");

  // Set tree size and allocate memory: Only one cell.
  gtot = Ncellmax = 1 ;
  AllocateTreeMemory();

  // Create bounding box of SPH particles
  for (k=0; k<ndim; k++) bbmin[k] = big_number;
  for (k=0; k<ndim; k++) bbmax[k] = -big_number;
  for (i=0; i<Ntot; i++) {
    for (k=0; k<ndim; k++) {
      bbmax[k] = max(bbmax[k], partdata[i].r[k] + kernrange*partdata[i].h);
      bbmin[k] = min(bbmin[k], partdata[i].r[k] - kernrange*partdata[i].h);
    }
  }

  Ncell = 1 ;

  // Set properties for the cell
  ifirst = _ifirst;
  ilast  = _ilast;
  celldata[0].N      = Npart;
  celldata[0].ifirst = ifirst;
  celldata[0].ilast  = ilast;
  celldata[0].cnext  = Ncellmax;
  for (k=0; k<ndim; k++) celldata[0].bbmin[k] = bbmin[k];
  for (k=0; k<ndim; k++) celldata[0].bbmax[k] = bbmax[k];
  for (k=0; k<ndim; k++) celldata[0].cexit[0][k] = -1;
  for (k=0; k<ndim; k++) celldata[0].cexit[1][k] = -1;
  for (i=ifirst; i<ilast; i++) inext[i] = i+1;
  inext[ilast] = -1;

  celldata[0].copen  = -1;
  celldata[0].cnext  = Ncellmax;
  celldata[0].id     = 0;
  celldata[0].level  = 0;
  ltot = 0 ;


  // If number of particles remains unchanged, use old id list
  if (Ntot > 0) {
    if (Ntot != Ntotold) {
      for (i=ifirst; i<=ilast; i++) ids[i] = i;
    }

#if defined(VERIFY_ALL)
    ValidateTree(partdata);
#endif
  }

  return;
}

//=================================================================================================
//  BruteForceTree::StockTree
/// Stock cell in BruteForce-tree.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void BruteForceTree<ndim,ParticleType,TreeCell>::StockTree
 (TreeCell<ndim> &cell,                ///< Reference to current tree cell
  ParticleType<ndim> *partdata)        ///< Particle data array
{
  int i;                               // Particle counter
  int iaux;                            // Aux. particle i.d. variable
  int k;                               // Dimension counter
  FLOAT dr[ndim];                      // Relative position vector
  FLOAT drsqd;                         // Distance squared
  FLOAT mi;                            // Mass of particle i
  FLOAT p = (FLOAT) 0.0;               // ..
  FLOAT lambda = (FLOAT) 0.0;          // ..


  // Zero all summation variables for all cells
  cell.Nactive  = 0;
  cell.N        = 0;
  cell.m        = (FLOAT) 0.0;
  cell.hmax     = (FLOAT) 0.0;
  cell.rmax     = (FLOAT) 0.0;
  cell.dhmaxdt  = (FLOAT) 0.0;
  cell.drmaxdt  = (FLOAT) 0.0;
  cell.mac      = (FLOAT) 0.0;
  cell.cdistsqd = big_number;
  for (k=0; k<5; k++) cell.q[k]          = (FLOAT) 0.0;
  for (k=0; k<ndim; k++) cell.r[k]       = (FLOAT) 0.0;
  for (k=0; k<ndim; k++) cell.v[k]       = (FLOAT) 0.0;
  for (k=0; k<ndim; k++) cell.rcell[k]   = (FLOAT) 0.0;
  for (k=0; k<ndim; k++) cell.bbmin[k]   = big_number;
  for (k=0; k<ndim; k++) cell.bbmax[k]   = -big_number;
  for (k=0; k<ndim; k++) cell.hboxmin[k] = big_number;
  for (k=0; k<ndim; k++) cell.hboxmax[k] = -big_number;
  for (k=0; k<ndim; k++) cell.vboxmin[k] = big_number;
  for (k=0; k<ndim; k++) cell.vboxmax[k] = -big_number;


  // If this is a leaf cell, sum over all particles
  //-----------------------------------------------------------------------------------------------
  assert(cell.level == ltot) ;
  {
    // First, check if any particles have been accreted and remove them
    // from the linked list.  If cell no longer contains any live particles,
    // then set N = 0 to ensure cell is not included in future tree-walks.
    i = cell.ifirst;
    cell.ifirst = -1;
    iaux = -1;
    while (i != -1) {
      if (!partdata[i].flags.is_dead()) {
        if (iaux == -1) cell.ifirst = i;
        else inext[iaux] = i;
        iaux = i;
      }
      if (i == cell.ilast) break;
      i = inext[i];
    };
    cell.ilast = iaux;


    // Loop over all particles in cell summing their contributions
    i = cell.ifirst;
    while (i != -1) {
      if (!partdata[i].flags.is_dead()) {
        cell.N++;
        if (partdata[i].active) cell.Nactive++;
        cell.hmax = max(cell.hmax,partdata[i].h);
        if (gravmask[partdata[i].ptype]) {
          cell.m += partdata[i].m;
          for (k=0; k<ndim; k++) cell.r[k] += partdata[i].m*partdata[i].r[k];
          for (k=0; k<ndim; k++) cell.v[k] += partdata[i].m*partdata[i].v[k];
        }
        for (k=0; k<ndim; k++) {
         if (partdata[i].r[k] < cell.bbmin[k]) cell.bbmin[k] = partdata[i].r[k];
         if (partdata[i].r[k] > cell.bbmax[k]) cell.bbmax[k] = partdata[i].r[k];
         if (partdata[i].r[k] - kernrange*partdata[i].h < cell.hboxmin[k])
         cell.hboxmin[k] = partdata[i].r[k] - kernrange*partdata[i].h;
         if (partdata[i].r[k] + kernrange*partdata[i].h > cell.hboxmax[k])
          cell.hboxmax[k] = partdata[i].r[k] + kernrange*partdata[i].h;
        }
      }
      if (i == cell.ilast) break;
      i = inext[i];
    };

    // Normalise all cell values
    if (cell.N > 0) {
      for (k=0; k<ndim; k++) cell.r[k] /= cell.m;
      for (k=0; k<ndim; k++) cell.v[k] /= cell.m;
      for (k=0; k<ndim; k++) cell.rcell[k] = (FLOAT) 0.5*(cell.bbmin[k] + cell.bbmax[k]);
      for (k=0; k<ndim; k++) dr[k] = (FLOAT) 0.5*(cell.bbmax[k] - cell.bbmin[k]);
      cell.cdistsqd = max(DotProduct(dr,dr,ndim),cell.hmax*cell.hmax)/thetamaxsqd;
      cell.rmax = sqrt(DotProduct(dr,dr,ndim));
    }

    // Compute quadrupole moment terms if selected
    if (multipole == "quadrupole") {
      i = cell.ifirst;

      while (i != -1) {
        if (!partdata[i].flags.is_dead() && gravmask[partdata[i].ptype]) {
          mi = partdata[i].m;
          for (k=0; k<ndim; k++) dr[k] = partdata[i].r[k] - cell.r[k];
          drsqd = DotProduct(dr,dr,ndim);
          if (ndim == 3) {
            cell.q[0] += mi*((FLOAT) 3.0*dr[0]*dr[0] - drsqd);
            cell.q[1] += mi*(FLOAT) 3.0*dr[0]*dr[1];
            cell.q[2] += mi*((FLOAT) 3.0*dr[1]*dr[1] - drsqd);
            cell.q[3] += mi*(FLOAT) 3.0*dr[2]*dr[0];
            cell.q[4] += mi*(FLOAT) 3.0*dr[2]*dr[1];
          }
          else if (ndim == 2) {
            cell.q[0] += mi*((FLOAT) 3.0*dr[0]*dr[0] - drsqd);
            cell.q[1] += mi*(FLOAT) 3.0*dr[0]*dr[1];
            cell.q[2] += mi*((FLOAT) 3.0*dr[1]*dr[1] - drsqd);
          }
        }
        if (i == cell.ilast) break;
        i = inext[i];
      }
    }

  }

  // Calculate eigenvalue MAC criteria
  if (gravity_mac == "eigenmac") {
    if (ndim == 3)
      p = cell.q[0]*cell.q[2] - (cell.q[0] + cell.q[2])*(cell.q[0] + cell.q[2]) -
        cell.q[1]*cell.q[1] - cell.q[3]*cell.q[3] - cell.q[4]*cell.q[4];
    if (p >= (FLOAT) 0.0) cell.mac = (FLOAT) 0.0;
    else {
      lambda = (FLOAT) 2.0*sqrt(-p/(FLOAT) 3.0);
      cell.mac = pow((FLOAT) 0.5*lambda/macerror,(FLOAT) 0.66666666666666);
    }
  }
  else {
    cell.mac = (FLOAT) 0.0;
  }

  return;
}


//=================================================================================================
//  BruteForceTree::UpdateActiveParticleCounters
/// Loop through all leaf cells in BruteForce-tree and update all active particle counters.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void BruteForceTree<ndim,ParticleType,TreeCell>::UpdateActiveParticleCounters
 (ParticleType<ndim> *partdata)        ///< [in] Main particle array
{
  int i;                               // SPH particle index
  int ilast;                           // Last particle in linked list

  debug2("[BruteForceTree::UpdateActiveParticleCounters]");


  celldata[0].Nactive = 0;
  i = celldata[0].ifirst;
  ilast = celldata[0].ilast;

  // Else walk through linked list to obtain list and number of active ptcls.
  while (i != -1) {
    if (i < Ntot && partdata[i].active && !partdata[i].flags.is_dead()) celldata[0].Nactive++;
    if (i == ilast) break;
    i = inext[i];
  };

  return;
}

//=================================================================================================
//  BruteForceTree::UpdateHmaxValues
/// Calculate the physical properties (e.g. total mass, centre-of-mass,
/// opening-distance, etc..) of all cells in the tree.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void BruteForceTree<ndim,ParticleType,TreeCell>::UpdateHmaxValues
 (TreeCell<ndim> &cell,                ///< BruteForce-tree cell
  ParticleType<ndim> *partdata)        ///< SPH particle data array
{
  int i;                               // Particle counter
  int k;                               // Dimension counter


  // Zero all summation variables for all cells
  cell.hmax = (FLOAT) 0.0;
  for (k=0; k<ndim; k++) cell.hboxmin[k] = big_number;
  for (k=0; k<ndim; k++) cell.hboxmax[k] = -big_number;


  // If this is a leaf cell, sum over all particles
  //-----------------------------------------------------------------------------------------------
  assert(cell.level == ltot) ;
  {
    i = cell.ifirst;

    // Loop over all particles in cell summing their contributions
    while (i != -1) {
      cell.hmax = max(cell.hmax,partdata[i].h);
      for (k=0; k<ndim; k++) {
        if (partdata[i].r[k] - kernrange*partdata[i].h < cell.hboxmin[k]) {
          cell.hboxmin[k] = partdata[i].r[k] - kernrange*partdata[i].h;
        }
        if (partdata[i].r[k] + kernrange*partdata[i].h > cell.hboxmax[k]) {
          cell.hboxmax[k] = partdata[i].r[k] + kernrange*partdata[i].h;
        }
      }
      if (i == cell.ilast) break;
      i = inext[i];
    };
  }
  return;
}





template class BruteForceTree<1,Particle,BruteForceTreeCell>;
template class BruteForceTree<2,Particle,BruteForceTreeCell>;
template class BruteForceTree<3,Particle,BruteForceTreeCell>;
template class BruteForceTree<1,SphParticle,BruteForceTreeCell>;
template class BruteForceTree<2,SphParticle,BruteForceTreeCell>;
template class BruteForceTree<3,SphParticle,BruteForceTreeCell>;
template class BruteForceTree<1,GradhSphParticle,BruteForceTreeCell>;
template class BruteForceTree<2,GradhSphParticle,BruteForceTreeCell>;
template class BruteForceTree<3,GradhSphParticle,BruteForceTreeCell>;
template class BruteForceTree<1,SM2012SphParticle,BruteForceTreeCell>;
template class BruteForceTree<2,SM2012SphParticle,BruteForceTreeCell>;
template class BruteForceTree<3,SM2012SphParticle,BruteForceTreeCell>;
template class BruteForceTree<1,MeshlessFVParticle,BruteForceTreeCell>;
template class BruteForceTree<2,MeshlessFVParticle,BruteForceTreeCell>;
template class BruteForceTree<3,MeshlessFVParticle,BruteForceTreeCell>;



