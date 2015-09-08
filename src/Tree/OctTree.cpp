//=================================================================================================
//  OctTree.cpp
//  Contains all functions for building, stocking and walking for the
//  octal-spatial tree for SPH particles.
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
#include "Debug.h"
#include "Exception.h"
#include "Parameters.h"
#include "InlineFuncs.h"
#include "Particle.h"
#include "OctTree.h"
#if defined _OPENMP
#include <omp.h>
#endif
using namespace std;



//=================================================================================================
//  OctTree::OctTree
/// OctTree constructor.  Initialises various variables.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
OctTree<ndim,ParticleType,TreeCell>::OctTree(int Nleafmaxaux, FLOAT thetamaxsqdaux,
                                             FLOAT kernrangeaux, FLOAT macerroraux,
                                             string gravity_mac_aux, string multipole_aux) :
  Tree<ndim,ParticleType,TreeCell>(Nleafmaxaux, thetamaxsqdaux, kernrangeaux,
                                   macerroraux, gravity_mac_aux, multipole_aux)
{
  allocated_tree = false;
  ifirst         = -1;
  ilast          = -1;
  lmax           = 40;
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
//  OctTree::~OctTree
/// OctTree destructor.  Deallocates tree memory upon object destruction.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
OctTree<ndim,ParticleType,TreeCell>::~OctTree()
{
  if (allocated_tree) DeallocateTreeMemory();
}



//=================================================================================================
//  OctTree::AllocateTreeMemory
/// Allocate memory for octal tree as requested.  If more memory is required
/// than currently allocated, tree is deallocated and reallocated here.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void OctTree<ndim,ParticleType,TreeCell>::AllocateTreeMemory(void)
{
  debug2("[OctTree::AllocateTreeMemory]");

  assert(Ntot >= 0);
  assert(lmax >= 0);

  if (!allocated_tree || Ntotmax > Ntotmaxold || Ncell > Ncellmax) {
    if (allocated_tree) DeallocateTreeMemory();
    Ntotmax     = max(Ntotmax, Ntot);
    Ntotmaxold  = Ntotmax;
    Ncellmax    = max((int) ((FLOAT) 2.0*(FLOAT) Ncellmax), 4*Ntotmax);
    Ncellmaxold = Ncellmax;
    gtot        = Ntotmax;

    firstCell = new int[lmax];
    lastCell  = new int[lmax];
    ids       = new int[Ntotmax];
    inext     = new int[Ntotmax];
    celldata  = new struct TreeCell<ndim>[Ncellmax];

    allocated_tree = true;
  }

  return;
}



//=================================================================================================
//  OctTree::DeallocateTreeMemory
/// Deallocates all octal tree memory
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void OctTree<ndim,ParticleType,TreeCell>::DeallocateTreeMemory(void)
{
  debug2("[OctTree::DeallocateTreeMemory]");

  if (allocated_tree) {
    delete[] celldata;
    delete[] inext;
    delete[] ids;
    delete[] lastCell;
    delete[] firstCell;
    allocated_tree = false;
  }

  return;
}



//=================================================================================================
//  OctTree::BuildTree
/// Call all routines to build/re-build the octal tree on the local node.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void OctTree<ndim,ParticleType,TreeCell>::BuildTree
 (const int _ifirst,                   ///< [in] i.d. of first particle
  const int _ilast,                    ///< [in] i.d. of last particle
  const int Npart,                     ///< [in] No. of particles
  const int Npartmax,                  ///< [in] Max. no. of particles
  const FLOAT timestep,                ///< [in] Smallest physical timestep
  ParticleType<ndim> *partdata)        ///< [in] Particle data array
{
  bool allDone = false;                // Are all cell divisions completed?
  int c;                               // Cell counter
  int cc;                              // Child cell counter
  int ckid;                            // i.d. fo child cell
  int cnew;                            // i.d. of newly created cell
  int ilast;                           // i.d. of last particle in cell
  int i;                               // Particle counter
  int k;                               // Dimension counter
  int kk;                              // ..
  int Nlist;                           // ..
  int *celllist;                       // List of cells to be processed
  FLOAT cellSize = (FLOAT) 0.0;        // Size of new cell (from centre to edge)

  debug2("[OctTree::BuildTree]");
  //timing->StartTimingSection("BUILD_OCT_TREE");

  // Allocate (or reallocate if needed) all tree memory
  AllocateTreeMemory();


  // Set properties for root cell before constructing tree
  ifirst = _ifirst;
  ilast  = _ilast;
  Ncell  = 0;
  ltot   = 0;
  for (k=0; k<ndim; k++) celldata[0].cexit[0][k] = -1;
  for (k=0; k<ndim; k++) celldata[0].cexit[1][k] = -1;
  for (i=ifirst; i<=ilast; i++) inext[i] = i+1;
  for (c=0; c<Ncellmax; c++) {
    celldata[c].N      = 0;
    celldata[c].ifirst = -1;
    celldata[c].ilast  = -1;
    celldata[c].copen  = -1;
    celldata[c].cnext  = Ncellmax;
    celldata[c].id     = c;
  }

  // Return now if tree contains no particles
  if (Ntot == 0) return;

  celldata[0].N      = Ntot;
  celldata[0].ifirst = _ifirst;
  celldata[0].ilast  = _ilast;
  celldata[0].level  = 0;
  celldata[0].copen  = -1;


  // Compute the bounding box of all particles in root cell and the root cell size
  for (k=0; k<ndim; k++) celldata[0].bbmin[k] = +big_number;
  for (k=0; k<ndim; k++) celldata[0].bbmax[k] = -big_number;
  for (i=_ifirst; i<=_ilast; i++) {
    for (k=0; k<ndim; k++) celldata[0].bbmin[k] = min(celldata[0].bbmin[k], partdata[i].r[k]);
    for (k=0; k<ndim; k++) celldata[0].bbmax[k] = max(celldata[0].bbmax[k], partdata[i].r[k]);
  }
  for (k=0; k<ndim; k++) {
    //celldata[0].rcell[k] = 0.5*(celldata[0].bbmin[k] + celldata[0].bbmax[k]);
    celldata[0].r[k] = (FLOAT) 0.5*(celldata[0].bbmin[k] + celldata[0].bbmax[k]);
    //cellSize = max(cellSize, celldata[0].bbmax[k] - celldata[0].rcell[k]);
    cellSize = max(cellSize, celldata[0].bbmax[k] - celldata[0].r[k]);
  }
  rootCellSize = (FLOAT) 2.0*cellSize;


  // Build tree if there are any particles
  //-----------------------------------------------------------------------------------------------
  if (Ntot > 0) {

    celllist     = new int[Npartmax];
    ltot         = 0;
    Ncell        = 1;
    Nlist        = 1;
    celllist[0]  = 0;
    firstCell[0] = 0;
    lastCell[0]  = 0;


    // Recursively divide tree to lower levels until all leaf cells contain maximum no. of ptcls.
    //---------------------------------------------------------------------------------------------
    while (!allDone) {
      cellSize *= (FLOAT) 0.5;


      // Loop over all unfinished cells to find new child cell occupancy
      //-------------------------------------------------------------------------------------------
      for (cc=0; cc<Nlist; cc++) {
        c = celllist[cc];
        TreeCell<ndim> &cell = celldata[c];

        // Create 8 child cells (N.B. different to SEREN; creates child cells even if empty
        // in order to fill all volume at all levels)
        //-----------------------------------------------------------------------------------------
        for (k=0; k<Noctchild; k++) {
          cnew = Ncell + cc*Noctchild + k;
          assert(cnew < Ncellmax);

          celldata[cnew].level  = ltot + 1;
          for (kk=0; kk<ndim; kk++) celldata[cnew].bbmin[kk] = cell.bbmin[kk];
          for (kk=0; kk<ndim; kk++) celldata[cnew].bbmax[kk] = cell.bbmax[kk];

          // Assign positions of child cells
          if (k == 0 || k == 2 || k == 4 || k == 6) {
            celldata[cnew].r[0] = cell.r[0] - cellSize;
            celldata[cnew].bbmax[0] = cell.r[0];
          }
          else {
            celldata[cnew].r[0] = cell.r[0] + cellSize;
            celldata[cnew].bbmin[0] = cell.r[0];
          }
          if (ndim > 1) {
            if (k == 0 || k == 1 || k == 4 || k == 5) {
              celldata[cnew].r[1] = cell.r[1] - cellSize;
              celldata[cnew].bbmax[1] = cell.r[1];
            }
            else {
              celldata[cnew].r[1] = cell.r[1] + cellSize;
              celldata[cnew].bbmin[1] = cell.r[1];
            }
          }
          if (ndim == 3) {
            if (k < 4) {
              celldata[cnew].r[2] = cell.r[2] - cellSize;
              celldata[cnew].bbmax[2] = cell.r[2];
            }
            else {
              celldata[cnew].r[2] = cell.r[2] + cellSize;
              celldata[cnew].bbmin[2] = cell.r[2];
            }
          }

        }
        //-----------------------------------------------------------------------------------------


        i = cell.ifirst;
        ilast = cell.ilast;

        // Walk through linked list of all particles to find new child cells
        //-----------------------------------------------------------------------------------------
        while (i != -1) {
          ckid = 0;

          // Find child cell i.d. (depending on dimensionality)
          if (partdata[i].r[0] > cell.r[0]) ckid += 1;
          if (ndim > 1) {
            if (partdata[i].r[1] > cell.r[1]) ckid += 2;
          }
          if (ndim == 3) {
            if (partdata[i].r[2] > cell.r[2]) ckid += 4;
          }
          assert(ckid >= 0 && ckid < Noctchild);

          cnew = Ncell + cc*Noctchild + ckid;

          // Walk through linked list of all particles to find new child cells
          if (celldata[cnew].ifirst == -1) {
            celldata[cnew].ifirst = i;
          }
          else {
            inext[celldata[cnew].ilast] = i;
          }
          celldata[cnew].ilast = i;
          celldata[cnew].N++;
          if (i == ilast) break;
          i = inext[i];
        };
        //-----------------------------------------------------------------------------------------


        // Set up linked lists from parent to children cells (avoiding empty cells)
        cell.copen = Ncell + cc*Noctchild;
        for (k=0; k<Noctchild; k++) {
          cnew = Ncell + cc*Noctchild + k;
          //if (cell.copen == -1 && celldata[cnew].N > 0) cell.copen = cnew;
          celldata[cnew].cnext = cnew + 1;
        }
        celldata[Ncell + cc*Noctchild + (Noctchild - 1)].cnext = cell.cnext;

      }
      //-------------------------------------------------------------------------------------------


      // Record first and last cells in newly created level
      Ncell += Noctchild*Nlist;
      ltot++;
      firstCell[ltot] = lastCell[ltot-1] + 1;
      lastCell[ltot] = Ncell - 1;
      assert(Ncell <= Ncellmax);

      // Check if all new cells are children cells
      Nlist = 0;
      allDone = true;
      for (c=firstCell[ltot]; c<=lastCell[ltot]; c++) {
        if (celldata[c].N > Nleafmax) {
          allDone = false;
          celllist[Nlist++] = c;
        }
      }

      // Check we have not reached maximum level
      if (ltot >= lmax) {
        cout << "Reached maximum Oct-tree level.  Exitting program" << endl;
        exit(0);
      }

    }
    //---------------------------------------------------------------------------------------------

    delete[] celllist;

  }
  //-----------------------------------------------------------------------------------------------


  StockTree(celldata[0],partdata);
#if defined(VERIFY_ALL)
  ValidateTree(partdata);
#endif

  return;
}



//=================================================================================================
//  OctTree::StockTree
/// Stock given tree cell in KD-tree.  If cell is not a leaf-cell, recursively
/// calls itself for its two child cells.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void OctTree<ndim,ParticleType,TreeCell>::StockTree
 (TreeCell<ndim> &rootcell,            ///< Reference to cell to be stocked
  ParticleType<ndim> *partdata)        ///< SPH particle data array
{
  int c,cc;                            // Cell counters
  int cend;                            // Last particle in cell
  int i;                               // Particle counter
  int iaux;                            // Aux. particle i.d. variable
  int k;                               // Dimension counter
  int l;                               // Level counter
  FLOAT dr[ndim];                      // Relative position vector
  FLOAT drsqd;                         // Distance squared
  FLOAT mi;                            // Mass of particle i
  FLOAT p = (FLOAT) 0.0;               // ??
  FLOAT lambda = (FLOAT) 0.0;          // Quadrupole moment eigenvalue

  debug2("[OctTree::StockTree]");


  // Loop over all levels in tree starting from lowest up to the top root cell level.
  //===============================================================================================
  for (l=ltot; l>=0; l--) {


    // Loop over all cells on current level
    //---------------------------------------------------------------------------------------------
    for (c=firstCell[l]; c<=lastCell[l]; c++) {
      TreeCell<ndim> &cell = celldata[c];

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
      for (k=0; k<ndim; k++) cell.r[k]       = (FLOAT) 0.0;
      for (k=0; k<ndim; k++) cell.v[k]       = (FLOAT) 0.0;
      for (k=0; k<ndim; k++) cell.rcell[k]   = (FLOAT) 0.0;
      for (k=0; k<ndim; k++) cell.bbmin[k]   = big_number;
      for (k=0; k<ndim; k++) cell.bbmax[k]   = -big_number;
      for (k=0; k<ndim; k++) cell.hboxmin[k] = big_number;
      for (k=0; k<ndim; k++) cell.hboxmax[k] = -big_number;
      for (k=0; k<5; k++) cell.q[k] = (FLOAT) 0.0;


      // If this is a leaf cell, sum over all particles
      //-------------------------------------------------------------------------------------------
      if (cell.copen == -1) {

        // First, check if any particles have been accreted and remove them
        // from the linked list.  If cell no longer contains any live particles,
        // then set N = 0 to ensure cell is not included in future tree-walks.
        i           = cell.ifirst;
        iaux        = -1;
        cell.ifirst = -1;
        cell.N      = 0;
        while (i != -1) {
          if (partdata[i].itype != dead) {
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
          if (partdata[i].itype != dead) {
            cell.N++;
            if (partdata[i].active) cell.Nactive++;
            cell.hmax = max(cell.hmax, partdata[i].h);
            cell.m += partdata[i].m;
            for (k=0; k<ndim; k++) cell.r[k] += partdata[i].m*partdata[i].r[k];
            for (k=0; k<ndim; k++) cell.v[k] += partdata[i].m*partdata[i].v[k];
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
            if (partdata[i].itype != dead) {
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
      // For non-leaf cells, sum over all child cells
      //-------------------------------------------------------------------------------------------
      else {

        // Set limits for children (maximum of 8 but may be less)
        cc   = cell.copen;
        cend = cell.cnext;

        while (cc != cend) {
          TreeCell<ndim> &child = celldata[cc];

          if (child.N > 0) {
            for (k=0; k<ndim; k++) cell.bbmin[k]   = min(child.bbmin[k], cell.bbmin[k]);
            for (k=0; k<ndim; k++) cell.bbmax[k]   = max(child.bbmax[k], cell.bbmax[k]);
            for (k=0; k<ndim; k++) cell.hboxmin[k] = min(child.hboxmin[k], cell.hboxmin[k]);
            for (k=0; k<ndim; k++) cell.hboxmax[k] = max(child.hboxmax[k], cell.hboxmax[k]);
            for (k=0; k<ndim; k++) cell.r[k] += child.m*child.r[k];
            for (k=0; k<ndim; k++) cell.v[k] += child.m*child.v[k];
            cell.hmax = max(child.hmax, cell.hmax);
            cell.m += child.m;
            cell.N += child.N;
          }

          cc = child.cnext;
        };

        if (cell.m > 0.0) {
          for (k=0; k<ndim; k++) cell.r[k] /= cell.m;
          for (k=0; k<ndim; k++) cell.v[k] /= cell.m;
        }
        for (k=0; k<ndim; k++) cell.rcell[k] = (FLOAT) 0.5*(cell.bbmin[k] + cell.bbmax[k]);
        for (k=0; k<ndim; k++) dr[k] = (FLOAT) 0.5*(cell.bbmax[k] - cell.bbmin[k]);
        cell.cdistsqd = max(DotProduct(dr, dr, ndim),cell.hmax*cell.hmax)/thetamaxsqd;
        cell.rmax = sqrt(DotProduct(dr, dr, ndim));
#ifdef MPI_PARALLEL
        //cell.worktot = child1.worktot + child2.worktot;
#endif

        // Set limits for children (maximum of 8 but may be less)
        cc   = cell.copen;
        cend = cell.cnext;

        while (cc != cend) {
          TreeCell<ndim> &child = celldata[cc];

          // Now add individual quadrupole moment terms
          if (multipole == "quadrupole" && child.N > 0) {
            mi = child.m;
            for (k=0; k<ndim; k++) dr[k] = child.r[k] - cell.r[k];
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

          cc = child.cnext;
        };


      }
      //-------------------------------------------------------------------------------------------


      // Calculate eigenvalue MAC criteria
      if (gravity_mac == "eigenmac") {
        if (ndim == 3)
          p = cell.q[0]*cell.q[2] - (cell.q[0] + cell.q[2])*(cell.q[0] + cell.q[2])
            - cell.q[1]*cell.q[1] - cell.q[3]*cell.q[3] - cell.q[4]*cell.q[4];
        if (p >= (FLOAT) 0.0) cell.mac = (FLOAT) 0.0;
        else {
          lambda = (FLOAT) 2.0*sqrt(-onethird*p);
          cell.mac = pow((FLOAT) 0.5*lambda/macerror, twothirds);
        }
      }
      else {
        cell.mac = (FLOAT) 0.0;
      }

    }
    //---------------------------------------------------------------------------------------------

  }
  //===============================================================================================

  //cout << "Finished stocking tree" << endl;

  return;
}



//=================================================================================================
//  OctTree::UpdateHmaxValues
/// Calculate the physical properties (e.g. total mass, centre-of-mass,
/// opening-distance, etc..) of all cells in the tree.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void OctTree<ndim,ParticleType,TreeCell>::UpdateHmaxValues
 (TreeCell<ndim> &rootcell,            ///< KD-tree cell
  ParticleType<ndim> *partdata)        ///< SPH particle data array
{
  int c,cc;                            // Cell counters
  int cfirst,cend;                     // ..
  int i;                               // Particle counter
  int k;                               // ..
  int l;                               // ..


  // Loop over all levels in tree starting from lowest
  //===============================================================================================
  for (l=ltot; l>=0; l--) {


    // Loop over all cells on current level
    //---------------------------------------------------------------------------------------------
    for (c=firstCell[l]; c<=lastCell[l]; c++) {
      TreeCell<ndim> &cell = celldata[c];

      // Zero all summation variables for all cells
      cell.hmax = 0.0;
      for (k=0; k<ndim; k++) cell.hboxmin[k] = big_number;
      for (k=0; k<ndim; k++) cell.hboxmax[k] = -big_number;

      // If this is a leaf cell, sum over all particles
      //-------------------------------------------------------------------------------------------
      if (cell.copen == -1) {
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
      // For non-leaf cells, sum over all child cells
      //-------------------------------------------------------------------------------------------
      else {

        // Set limits for children (maximum of 8 but may be less)
        cfirst = cell.copen;
        cend   = cell.cnext;

        cc = cfirst;
        while (cc != cend) {
          TreeCell<ndim> &child = celldata[cc];

          if (child.N > 0) {
            for (k=0; k<ndim; k++) cell.hboxmin[k] = min(child.hboxmin[k],cell.hboxmin[k]);
            for (k=0; k<ndim; k++) cell.hboxmax[k] = max(child.hboxmax[k],cell.hboxmax[k]);
            cell.hmax = max(child.hmax,cell.hmax);
          }

          cc = child.cnext;
        };

      }
      //-------------------------------------------------------------------------------------------

    }
    //---------------------------------------------------------------------------------------------

  }
  //===============================================================================================


  return;
}


//=================================================================================================
//  OctTree::UpdateActiveParticleCounters
/// Loop through all leaf cells in KD-tree and update all active particle counters.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void OctTree<ndim,ParticleType,TreeCell>::UpdateActiveParticleCounters
 (ParticleType<ndim> *partdata)        ///< ..
{
  int c;                               // Cell counter
  int i;                               // SPH particle index
  int ilast;                           // Last particle in linked list

  debug2("[OctTree::UpdateActiveParticleCounters]");
  //timing->StartTimingSection("TREE_UPDATE_COUNTERS");


  // Loop through all grid cells in turn
  //-----------------------------------------------------------------------------------------------
#pragma omp parallel for default(none) private(c,i,ilast) shared(partdata)
  for (c=0; c<Ncell; c++) {
    celldata[c].Nactive = 0;

    if (celldata[c].copen != -1) continue;
    i = celldata[c].ifirst;
    ilast = celldata[c].ilast;

    // Else walk through linked list to obtain list and number of active ptcls.
    while (i != -1) {
      if (partdata[i].active) celldata[c].Nactive++;
      if (i == ilast) break;
      i = inext[i];
    };

  }
  //-----------------------------------------------------------------------------------------------

  //timing->EndTimingSection("TREE_UPDATE_COUNTERS");

  return;
}



#if defined(VERIFY_ALL)
//=================================================================================================
//  OctTree::ValidateTree
/// Performs various sanity and validation checks on KD-tree structure.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void OctTree<ndim,ParticleType,TreeCell>::ValidateTree
 (ParticleType<ndim> *partdata)        ///< [in] Pointer to SPH class
{
  bool kill_flag = false;              // Flag if tree is invalid to terminate program
  int activecount;                     // Active particles in leaf cell
  int c;                               // Cell counter
  int i;                               // Particle counter
  int k;                               // Dimension counter
  int leafcount;                       // Leaf cell counter
  int Nactivecount=0;                  // Counter for total no. of active ptcls
  int Ncount=0;                        // Total particle counter
  int *ccount;                         // Array for counting cells
  int *pcount;                         // Array for counting particles in tree
  TreeCell<ndim> cell;                 // Local copy of Oct-tree cell

  debug2("[OctTree::ValidateTree]");

  ccount = new int[Ncellmax];
  pcount = new int[Ntotmax];
  for (i=0; i<Ntotmax; i++) pcount[i] = 0;
  for (c=0; c<Ncellmax; c++) ccount[c] = 0;
  Ncount = 0;
  Nactivecount = 0;

  // Count how many times we enter a cell in a full tree walk
  c = 0;
  while (c < Ncell) {
    ccount[c]++;
    if (celldata[c].copen != -1) c = celldata[c].copen;
    else c = celldata[c].cnext;
  }

  // Now check we enter all cells once and once only
  for (c=0; c<Ncell; c++) {
    if (ccount[c] != 1) {
      cout << "Error in cell walk count : " << ccount[c] << endl;
      PrintArray("ccount     : ",Ncell,ccount);
      exit(0);
    }
    if (celldata[c].level < 0 || celldata[c].level > ltot) {
      cout << "Problem with cell levels : " << celldata[c].level << "    " << ltot << endl;
      exit(0);
    }
    if (celldata[c].Nactive > 0 && celldata[c].copen != -1) {
      cout << "Problem with active counters : " << celldata[c].Nactive << "    "
           << celldata[c].N << "    " << celldata[c].level << endl;
      exit(0);
    }
  }

  // Check inext linked list values and ids array are all valid
  for (i=ifirst; i<=ilast; i++) {
    //if (!(ids[i] >= ifirst && ids[i] <= ilast)) {
    //  cout << "Problem with ids array : " << i << "   " << ids[i] << endl;
    //  exit(0);
    //}
    if (!(inext[i] >= -1)) {
      cout << "Problem with inext linked lists : " << i << "   " << inext[i] << endl;
      exit(0);
    }
  }


  // Loop over all cells in tree
  //-----------------------------------------------------------------------------------------------
  for (c=0; c<Ncell; c++) {
    cell = celldata[c];
    activecount = 0;
    leafcount = 0;

    // Check that particles are not in linked lists more than once
    if (cell.copen == -1) {
      i = cell.ifirst;
      while (i != -1) {
        pcount[i]++;
        leafcount++;
        Ncount++;
        if (partdata[i].active) activecount++;
        if (partdata[i].active) Nactivecount++;
        if (partdata[i].h > cell.hmax) {
          cout << "hmax flag error : " << c << "    "
               << partdata[i].h << "   " << cell.hmax << endl;
          exit(0);
        }
        for (k=0; k<ndim; k++) {
          if (partdata[i].r[k] < cell.bbmin[k] || partdata[i].r[k] > cell.bbmax[k]) {
            cout << "Bounding box error : " << c << "   " << i << "    " << k << "    r : "
                 << partdata[i].r[k] << "    " << cell.bbmin[k] << "    " << cell.bbmax[k] << endl;
            exit(0);
          }
        }
        if (i == cell.ilast) break;
        i = inext[i];
      }
      if (leafcount > Nleafmax) {
        cout << "Leaf particle count error : " << leafcount << "   " << Nleafmax << endl;
        exit(0);
      }
      if (activecount > leafcount) {
        cout << "Leaf particle count error : " << leafcount << "   " << Nleafmax << endl;
        exit(0);
      }
    }

  }
  //-----------------------------------------------------------------------------------------------

  // Check particles are included in the tree once and once only
  for (i=ifirst; i<=ilast; i++) {
    if (pcount[i] != 1) {
      cout << "Problem with child cell ptcl counter : " << i << "   " << pcount[i] << endl;
      kill_flag = true;
    }
  }

  // Check all particles accounted for
  if (Ncount != Ntot) {
    cout << "Ncount problem with KD-tree : " << Ncount << "   " << Ntot << endl;
    kill_flag = true;
  }

  // Check active particles don't exceed total number of particles
  if (Nactivecount > Ntot) {
    cout << "Nactivecount problem with KD-tree : " << Nactivecount << "   " << Ntot << endl;
    kill_flag = true;
  }


  delete[] pcount;
  delete[] ccount;

  if (kill_flag) exit(0);

  return;
}
#endif



template class OctTree<1, Particle, OctTreeCell>;
template class OctTree<2, Particle, OctTreeCell>;
template class OctTree<3, Particle, OctTreeCell>;

template class OctTree<1, SphParticle, OctTreeCell>;
template class OctTree<2, SphParticle, OctTreeCell>;
template class OctTree<3, SphParticle, OctTreeCell>;

template class OctTree<1, GradhSphParticle, OctTreeCell>;
template class OctTree<2, GradhSphParticle, OctTreeCell>;
template class OctTree<3, GradhSphParticle, OctTreeCell>;

template class OctTree<1, SM2012SphParticle, OctTreeCell>;
template class OctTree<2, SM2012SphParticle, OctTreeCell>;
template class OctTree<3, SM2012SphParticle, OctTreeCell>;

template class OctTree<1, MeshlessFVParticle, OctTreeCell>;
template class OctTree<2, MeshlessFVParticle, OctTreeCell>;
template class OctTree<3, MeshlessFVParticle, OctTreeCell>;

//template class OctTree<3, GradhSphParticle, OsTreeRayCell>;
//template class OctTree<3, MeshlessFVParticle, OsTreeRayCell>;
