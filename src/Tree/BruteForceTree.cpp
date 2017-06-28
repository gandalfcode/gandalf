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
#include <sstream>
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
                                           	   	   	   	   string gravity_mac_aux,
                                           	   	   	   	   multipole_method multipole_aux,
                                           	   	   	   	   const DomainBox<ndim>& domain,
                                           	   	   	   	   const ParticleTypeRegister& reg,
														   const bool IAmPruned):
  Tree<ndim,ParticleType,TreeCell>(Nleafmaxaux, thetamaxsqdaux, kernrangeaux,
                                   macerroraux, gravity_mac_aux, multipole_aux, domain, reg,
                                   IAmPruned)
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
  Ntot           = 0;
  Ntotmax        = 0;
  hmax           = 0.0;
#if defined _OPENMP
  Nthreads       = omp_get_max_threads();
#else
  Nthreads       = 1;
#endif
#if defined MPI_PARALLEL
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
void BruteForceTree<ndim,ParticleType,TreeCell>::AllocateTreeMemory(int Nparticles, int Ncells, bool force_realloc)
{
  debug2("[BruteForceTree::AllocateTreeMemory]");

  if (!allocated_tree || Nparticles > Ntotmax || Ncells>Ncellmax || force_realloc) {
    if (allocated_tree) DeallocateTreeMemory();

    Nparticles     = max(Nparticles, Ntotmax);
    Ncells		   = max(Ncells, Ncellmax);

    g2c      = new int[gmax];
    celldata = new struct TreeCell<ndim>[Ncells];

    allocated_tree = true;

    Ntotmax = Nparticles;
    Ncellmax = Ncells;
  }

  assert(Ntotmax>=Ntot);
  assert(Ncell<=Ncellmax);


  return;
}

//=================================================================================================
//  BruteForceTree::ReallocateMemory
/// Reallocate memory for KD-tree (when we need to grow the tree) as requested. Preserves the existing
/// information in the tree, differently from the previous function
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void BruteForceTree<ndim,ParticleType,TreeCell>::ReallocateMemory(int Nparticles, int Ncells)
{
  debug2("[BruteForceTree::ReallocateMemory]");

  if (!allocated_tree) {
	  ExceptionHandler::getIstance().raise("This function should not be called if the tree has not been allocated yet!");
  }


  if (Nparticles > Ntotmax ) {
	  Ntotmax = Ntot;
  }

  if (Ncells > Ncellmax) {


    TreeCell<ndim>* celldataold = celldata;

    celldata = new struct TreeCell<ndim>[Ncells];

    std::copy(celldataold,celldataold+Ncellmax,celldata);

    delete[] celldataold;

	Ncellmax=Ncells;


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
  Particle<ndim> *part_gen)            ///< Particle data array
{
  int i;                               // Particle counter
  int k;                               // Dimension counter
  FLOAT bbmin[ndim];                   // Minimum extent of local bounding box
  FLOAT bbmax[ndim];                   // Maximum extent of local bounding box

  debug2("[BruteForceTree::BuildTree]");
  //timing->StartTimingSection("BUILD_TREE");
  ParticleType<ndim>* partdata = reinterpret_cast<ParticleType<ndim>*>(part_gen);


  // Set tree size and allocate memory: Only one cell.
  bool force_realloc=false;
  const int gmax_old = gmax;
  gmax = Npartmax + 1 ;
  if (gmax > gmax_old) force_realloc=true;
  gtot = Ncell = Ntot + 1;
  AllocateTreeMemory(Npartmax,Npartmax+1,force_realloc);

  for (k=0; k<ndim; k++) bbmin[k] = big_number;
  for (k=0; k<ndim; k++) bbmax[k] = -big_number;

  if (Npart > 0) {
    ifirst = _ifirst;
    ilast  = _ifirst + Npart - 1;
    for (i=ifirst; i<=ilast; i++) {
      for (k=0; k<ndim; k++) {
        bbmax[k] = max(bbmax[k], partdata[i].r[k] + kernrange*partdata[i].h);
        bbmin[k] = min(bbmin[k], partdata[i].r[k] - kernrange*partdata[i].h);
      }
    }
  }
  else {
    ifirst = -1;
    ilast  = -1;
  }


  // Set properties for the cell
  celldata[0].N      = Ntot;
  celldata[0].ifirst = ifirst;
  celldata[0].ilast  = ilast;
  celldata[0].cnext  = Ncell;
  celldata[0].copen  = Ntot > 0 ? 1 : -1;
  celldata[0].id     = 0;
  celldata[0].level  = 0;
  celldata[0].hmax = 0;
  celldata[0].parent = -1;
  for (k=0; k<ndim; k++) celldata[0].bb.min[k] = bbmin[k];
  for (k=0; k<ndim; k++) celldata[0].bb.max[k] = bbmax[k];

  // Now do the leaf cells
  if (Npart > 0) {
    i = ifirst;
    for (int c = 1; c < Ncell; c++) {
      celldata[c].N      = 1 ;
      celldata[c].ifirst = i;
      celldata[c].ilast  = i;
      celldata[c].cnext  = c+1;
      celldata[c].copen  = -1;
      celldata[c].id     = c;
      celldata[c].level  = 1;
      celldata[c].parent = 0;
      for (k=0; k<ndim; k++) celldata[c].bb.min[k] = partdata[i].r[k] - kernrange*partdata[i].h ;
      for (k=0; k<ndim; k++) celldata[c].bb.max[k] = partdata[i].r[k] + kernrange*partdata[i].h ;
#ifdef MPI_PARALLEL
      celldata[c].worktot = 0.0;
#endif
      g2c[c-1] = c;
      i++ ;
    }
    ltot = 1;
    if (Ntot > 0) StockTree(celldata[0], partdata, true);
    else ltot = 0;
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
  ParticleType<ndim> *partdata,		   ///< Pointer to particle array
  bool stock_leaf)					   ///< Whether to stock also leaf cells
{
  int c = cell.copen ;

  // Stock children first
  if (c != -1) {
    for (; c < cell.cnext; c++)
      if (stock_leaf) StockTreeProperties(celldata[c], partdata, stock_leaf) ;
  }

  StockTreeProperties(cell, partdata, stock_leaf) ;
  return ;
}
//=================================================================================================
//  BruteForceTree::StockTreeProperties
/// Stock cell in BruteForce-tree.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void BruteForceTree<ndim,ParticleType,TreeCell>::StockTreeProperties
(TreeCell<ndim> &cell,                ///< Reference to current tree cell
 ParticleType<ndim> *partdata,		  ///< Particle data array
 bool stock_leaf)                     ///< Whether to stock also leaf cells
 {
  int i;                               // Particle counter
  int k;                               // Dimension counter
  FLOAT dr[ndim];                      // Relative position vector
  FLOAT drsqd;                         // Distance squared
  FLOAT mi;                            // Mass of particle i
  FLOAT p = (FLOAT) 0.0;               // ..
  FLOAT lambda = (FLOAT) 0.0;          // ..


  const bool need_quadrupole_moments =
      multipole == quadrupole || multipole == fast_quadrupole || gravity_mac == eigenmac ;

  if (stock_leaf || cell.copen != -1) {
    // Zero all summation variables for all cells
    cell.Nactive  = 0;
    cell.N        = 0;
    cell.maxsound = 0.0f;
    cell.m        = (FLOAT) 0.0;
    cell.hmax     = (FLOAT) 0.0;
    cell.rmax     = (FLOAT) 0.0;
    cell.mac      = (FLOAT) 0.0;
    cell.cdistsqd = big_number;
    if (gravity_mac == gadget2)
      cell.amin = big_number ;
    else if (gravity_mac == eigenmac)
      cell.macfactor = 0 ;
    for (k=0; k<5; k++) cell.q[k]          = (FLOAT) 0.0;
    for (k=0; k<ndim; k++) cell.r[k]       = (FLOAT) 0.0;
    for (k=0; k<ndim; k++) cell.bb.min[k]   = big_number;
    for (k=0; k<ndim; k++) cell.bb.max[k]   = -big_number;
    for (k=0; k<ndim; k++) cell.hbox.min[k] = big_number;
    for (k=0; k<ndim; k++) cell.hbox.max[k] = -big_number;
    for (k=0; k<ndim; k++) cell.vbox.min[k] = big_number;
    for (k=0; k<ndim; k++) cell.vbox.max[k] = -big_number;
  }

  if (cell.copen == -1) {
    if (not stock_leaf) return ;

    // Loop over all particles in cell summing their contributions
    i = cell.ifirst;
    int ilast = cell.ilast ;
    for(; i <= ilast; ++i) {
      if (!partdata[i].flags.is_dead()) {
        cell.N++;
        if (partdata[i].flags.check(active)) cell.Nactive++;
        cell.hmax = max(cell.hmax,partdata[i].h);
        cell.maxsound = max(cell.maxsound, partdata[i].sound);
        if (gravmask[partdata[i].ptype]) {
          cell.m += partdata[i].m;
          for (k=0; k<ndim; k++) cell.r[k] += partdata[i].m*partdata[i].r[k];
        }
        for (k=0; k<ndim; k++) {
          if (partdata[i].r[k] < cell.bb.min[k]) cell.bb.min[k] = partdata[i].r[k];
          if (partdata[i].r[k] > cell.bb.max[k]) cell.bb.max[k] = partdata[i].r[k];
          if (partdata[i].r[k] - kernrange*partdata[i].h < cell.hbox.min[k])
            cell.hbox.min[k] = partdata[i].r[k] - kernrange*partdata[i].h;
          if 	(partdata[i].r[k] + kernrange*partdata[i].h > cell.hbox.max[k])
            cell.hbox.max[k] = partdata[i].r[k] + kernrange*partdata[i].h;
          if (partdata[i].v[k] > cell.vbox.max[k]) cell.vbox.max[k] = partdata[i].v[k];
          if (partdata[i].v[k] < cell.vbox.min[k]) cell.vbox.min[k] = partdata[i].v[k];
        }
        if (gravity_mac == gadget2)
          cell.amin = min(cell.amin,
              sqrt(DotProduct(partdata[i].atree,partdata[i].atree,ndim)));
        else if (gravity_mac == eigenmac)
          cell.macfactor = max(cell.macfactor,pow(partdata[i].gpot,-twothirds));
      }
    }

    // Normalise all cell values
    if (cell.m > 0) {
      for (k=0; k<ndim; k++) cell.r[k] /= cell.m;
      for (k=0; k<ndim; k++) dr[k] = (FLOAT) 0.5*(cell.bb.max[k] - cell.bb.min[k]);
      cell.cdistsqd = max(DotProduct(dr,dr,ndim),cell.hmax*cell.hmax)/thetamaxsqd;
      cell.rmax = sqrt(DotProduct(dr,dr,ndim));
    }


    // Compute quadrupole moment terms if selected
    if (need_quadrupole_moments) {
      i = cell.ifirst;

      int ilast = cell.ilast ;
      for(; i <= ilast; ++i) {
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
      }
    }

    // Calculate eigenvalue MAC criteria
    if (gravity_mac == eigenmac) {
      if (ndim == 3) {
        p = cell.q[0]*cell.q[2] - (cell.q[0] + cell.q[2])*(cell.q[0] + cell.q[2]) -
            cell.q[1]*cell.q[1] - cell.q[3]*cell.q[3] - cell.q[4]*cell.q[4];
        if (p >= (FLOAT) 0.0) {
          lambda = 0;
        } else {
          lambda = (FLOAT) 2.0*sqrt(-p/(FLOAT) 3.0);
        }
      } else if (ndim == 2) {
        p = (cell.q[0]-cell.q[2])*(cell.q[0]-cell.q[2]) + 4*cell.q[1]*cell.q[1];
        lambda = (FLOAT) 0.5*max(cell.q[0] + cell.q[2] + sqrt(p), (FLOAT) 0.0);
      } else {
        lambda = fabs(cell.q[0]) ;
      }

      cell.mac = pow((FLOAT) 0.5*lambda/macerror,(FLOAT) 0.66666666666666);

    }
    else {
      cell.mac = (FLOAT) 0.0;
    }
  }
  // For non-leaf cells, sum over all child cells
  //-------------------------------------------------------------------------------------------
  else if (cell.copen != -1) {

    // Set limits for children
    int cc   = cell.copen;
    int cend = cell.cnext;

    assert(cend == Ncell);

#ifdef MPI_PARALLEL
    cell.worktot = 0 ;
#endif

    while (cc != cend) {
      TreeCell<ndim> &child = celldata[cc];

      if (child.N > 0) {
        for (k=0; k<ndim; k++) cell.bb.min[k]   = min(child.bb.min[k], cell.bb.min[k]);
        for (k=0; k<ndim; k++) cell.bb.max[k]   = max(child.bb.max[k], cell.bb.max[k]);
        for (k=0; k<ndim; k++) cell.hbox.min[k] = min(child.hbox.min[k], cell.hbox.min[k]);
        for (k=0; k<ndim; k++) cell.hbox.max[k] = max(child.hbox.max[k], cell.hbox.max[k]);
        for (k=0; k<ndim; k++) cell.vbox.min[k] = min(child.vbox.min[k],cell.vbox.min[k]);
        for (k=0; k<ndim; k++) cell.vbox.max[k] = max(child.vbox.max[k],cell.vbox.max[k]);
        for (k=0; k<ndim; k++) cell.r[k] += child.m*child.r[k];
        cell.hmax = max(child.hmax, cell.hmax);
        cell.maxsound = max(cell.maxsound, child.maxsound);
        if (gravity_mac == gadget2)
          cell.amin = min(cell.amin, child.amin);
        else if (gravity_mac == eigenmac)
          cell.macfactor = max(cell.macfactor, child.macfactor) ;
        cell.m += child.m;
        cell.N += child.N;
      }

#ifdef MPI_PARALLEL
    cell.worktot += child.worktot ;
#endif

      cc = child.cnext;
    };

    if (cell.m > 0.0) {
      for (k=0; k<ndim; k++) cell.r[k] /= cell.m;
    }
    for (k=0; k<ndim; k++) dr[k] = (FLOAT) 0.5*(cell.bb.max[k] - cell.bb.min[k]);
    cell.cdistsqd = max(DotProduct(dr, dr, ndim),cell.hmax*cell.hmax)/thetamaxsqd;

    cell.rmax = sqrt(DotProduct(dr, dr, ndim));

    // Set limits for children
    cc   = cell.copen;
    cend = cell.cnext;

    while (cc != cend) {
      TreeCell<ndim> &child = celldata[cc];

      // Now add individual quadrupole moment terms
      if (need_quadrupole_moments && child.N > 0) {
        mi = child.m;
        for (k=0; k<ndim; k++) dr[k] = child.r[k] - cell.r[k];
        drsqd = DotProduct(dr,dr,ndim);
        if (ndim == 3) {
          for (k=0; k<5; k++) cell.q[k] += child.q[k] ;
          cell.q[0] += mi*((FLOAT) 3.0*dr[0]*dr[0] - drsqd);
          cell.q[1] += mi*(FLOAT) 3.0*dr[0]*dr[1];
          cell.q[2] += mi*((FLOAT) 3.0*dr[1]*dr[1] - drsqd);
          cell.q[3] += mi*(FLOAT) 3.0*dr[2]*dr[0];
          cell.q[4] += mi*(FLOAT) 3.0*dr[2]*dr[1];
        }
        else if (ndim == 2) {
          for (k=0; k<3; k++) cell.q[k] += child.q[k] ;
          cell.q[0] += mi*((FLOAT) 3.0*dr[0]*dr[0] - drsqd);
          cell.q[1] += mi*(FLOAT) 3.0*dr[0]*dr[1];
          cell.q[2] += mi*((FLOAT) 3.0*dr[1]*dr[1] - drsqd);
        }
        else if (ndim == 1) {
          cell.q[0] += child.q[0] ;
          cell.q[0] += mi*((FLOAT) 3.0*dr[0]*dr[0] - drsqd);
        }
      }

      cc = child.cnext;
    }

  }
  //-------------------------------------------------------------------------------------------

  // Calculate eigenvalue MAC criteria
  if (gravity_mac == eigenmac) {
    if (ndim == 3) {
      p = cell.q[0]*cell.q[2] - (cell.q[0] + cell.q[2])*(cell.q[0] + cell.q[2]) -
          cell.q[1]*cell.q[1] - cell.q[3]*cell.q[3] - cell.q[4]*cell.q[4];
      if (p >= (FLOAT) 0.0) {
        lambda = 0;
      } else {
        lambda = (FLOAT) 2.0*sqrt(-p/(FLOAT) 3.0);
      }
    } else if (ndim == 2) {
      p = (cell.q[0]-cell.q[2])*(cell.q[0]-cell.q[2]) + 4*cell.q[1]*cell.q[1];
      lambda = (FLOAT) 0.5*max(cell.q[0] + cell.q[2] + sqrt(p), (FLOAT) 0.0);
    } else {
      lambda = fabs(cell.q[0]) ;
    }

    cell.mac = pow((FLOAT) 0.5*lambda/macerror,(FLOAT) 0.66666666666666);

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
 (Particle<ndim> *part_gen)            ///< [in] Main particle array
{
  int i;                               // SPH particle index
  int c;                               // Cell index
  int ilast;                           // Last particle in linked list

  debug2("[BruteForceTree::UpdateActiveParticleCounters]");

  ParticleType<ndim>* partdata = reinterpret_cast<ParticleType<ndim>*>(part_gen);

#pragma omp parallel for default(none) private(c,i,ilast) shared(partdata)
  for (c=0; c<Ncell; c++) {
    celldata[c].Nactive = 0;

    if (celldata[c].level != ltot) continue;
    i = celldata[c].ifirst;
    ilast = celldata[c].ilast;

    // Else walk through linked list to obtain list and number of active ptcls.
    for(; i <= ilast; ++i) {
      if (i < Ntot && partdata[i].flags.check(active) && !partdata[i].flags.is_dead())
        celldata[c].Nactive++;
    };

  }

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
  ParticleType<ndim> *partdata,        ///< SPH particle data array
  bool stock_leaf)                     ///< Whether to stock the leaf cells.
  {
  UpdateHmaxValuesCell(cell, partdata, stock_leaf) ;

  int c = cell.copen ;
  if (c == -1) c = cell.cnext ;
  for (; c < Ncell; c++)
	UpdateHmaxValuesCell(celldata[c], partdata, stock_leaf) ;
}
//=================================================================================================
//  BruteForceTree::UpdateHmaxValuesCell
/// Calculate the physical properties (e.g. total mass, centre-of-mass,
/// opening-distance, etc..) of all cells in the tree.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void BruteForceTree<ndim,ParticleType,TreeCell>::UpdateHmaxValuesCell
 (TreeCell<ndim> &cell,                ///< BruteForce-tree cell
  ParticleType<ndim> *partdata,        ///< SPH particle data array
  bool stock_leaf)                     ///< Whether to stock the leaf cells.
{
  int i;                               // Particle counter
  int k;                               // Dimension counter

  if (cell.copen == -1 && not stock_leaf) return ;

  // Zero all summation variables for all cells
  cell.hmax = (FLOAT) 0.0;
  for (k=0; k<ndim; k++) cell.hbox.min[k] = big_number;
  for (k=0; k<ndim; k++) cell.hbox.max[k] = -big_number;


  // If this is a leaf cell, sum over all particles
  //-----------------------------------------------------------------------------------------------
  i = cell.ifirst;

  // Loop over all particles in cell summing their contributions
  int ilast = cell.ilast ;
  for(; i <= ilast; ++i) {
    cell.hmax = max(cell.hmax,partdata[i].h);
    for (k=0; k<ndim; k++) {
      if (partdata[i].r[k] - kernrange*partdata[i].h < cell.hbox.min[k]) {
    	cell.hbox.min[k] = partdata[i].r[k] - kernrange*partdata[i].h;
      }
      if (partdata[i].r[k] + kernrange*partdata[i].h > cell.hbox.max[k]) {
    	cell.hbox.max[k] = partdata[i].r[k] + kernrange*partdata[i].h;
      }
    }
  }

  return;
}



#ifdef MPI_PARALLEL
//=================================================================================================
//  BruteForceTree::UpdateWorkCounters
/// Calculate the physical properties (e.g. total mass, centre-of-mass,
/// opening-distance, etc..) of all cells in the tree.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void BruteForceTree<ndim,ParticleType,TreeCell>::UpdateWorkCounters()
{

  double worktot = 0 ;

  for (int c=1; c <Ncell; c++)
    worktot += celldata[c].worktot ;

  celldata[0].worktot = worktot ;

  return;
}
#endif




template class BruteForceTree<1,GradhSphParticle,BruteForceTreeCell>;
template class BruteForceTree<2,GradhSphParticle,BruteForceTreeCell>;
template class BruteForceTree<3,GradhSphParticle,BruteForceTreeCell>;
template class BruteForceTree<1,SM2012SphParticle,BruteForceTreeCell>;
template class BruteForceTree<2,SM2012SphParticle,BruteForceTreeCell>;
template class BruteForceTree<3,SM2012SphParticle,BruteForceTreeCell>;
template class BruteForceTree<1,MeshlessFVParticle,BruteForceTreeCell>;
template class BruteForceTree<2,MeshlessFVParticle,BruteForceTreeCell>;
template class BruteForceTree<3,MeshlessFVParticle,BruteForceTreeCell>;
