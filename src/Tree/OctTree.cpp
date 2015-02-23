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
#include "Exception.h"
#include "Parameters.h"
#include "InlineFuncs.h"
#include "Particle.h"
#include "Sph.h"
#include "OctTree.h"
#include "Debug.h"
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
  /*gravity_mac(gravity_mac_aux),
  multipole(multipole_aux),
  Nleafmax(Nleafmaxaux),
  invthetamaxsqd(1.0/thetamaxsqdaux),
  kernrange(kernrangeaux),
  macerror(macerroraux),
  theta(sqrt(thetamaxsqdaux)),
  thetamaxsqd(thetamaxsqdaux)*/
{
  allocated_tree = false;
  ltot           = 0;
  lmax           = 40;
  Ncell          = 0;
  Ncellmax       = 0;
  Ntot           = 0;
  Ntotmax        = 0;
  Ntotmaxold     = 0;
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

  if (!allocated_tree || Ntotmax > Ntotmaxold || Ncell > Ncellmax) {
    if (allocated_tree) DeallocateTreeMemory();
    Ntotmax = max(Ntotmax,Ntot);
    Ntotmaxold = Ntotmax;
    Ncellmax = max((int) (1.5*(FLOAT) Ncellmax), 2*Ntotmax);
    gmax = Ntotmax;
    gtot = Ntotmax;

    firstCell = new int[lmax];
    lastCell  = new int[lmax];
    g2c       = new int[gmax];
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
    delete[] g2c;
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
 (int _ifirst,                         ///< [in] i.d. of first particle
  int _ilast,                          ///< [in] i.d. of last particle
  int Npart,                           ///< [in] No. of particles
  int Npartmax,                        ///< [in] Max. no. of particles
  ParticleType<ndim> *partdata,        ///< [in] Particle data array
  FLOAT timestep)                      ///< [in] Smallest physical timestep
{
  bool allDone = false;                // Are all cell divisions completed?
  int c;                               // Cell counter
  int cc;                              // Child cell counter
  int ckid;                            // ..
  int clist[Noctchild];                // ..
  int ilast;                           // ..
  int i;                               // Particle counter
  int k;                               // Dimension counter
  int Nlist;                           // ..
  int Nincell[Noctchild];              // ..
  int Nkids;                           // ..
  int *celllist;                       // ..
  int *whichChild;                     // ..
  FLOAT cellSize = (FLOAT) 0.0;        // Size of new cell (from centre to edge)

  debug2("[OctTree::BuildTree]");
  //timing->StartTimingSection("BUILD_OCT_TREE");

  // Allocate (or reallocate if needed) all tree memory
  AllocateTreeMemory();

  //cout << "Building Oct-tree with " << Npart << " particles.   first/last : "
  //     << _ifirst << "   " << _ilast << endl;

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
    for (k=0; k<Noctchild; k++) celldata[c].childof[k] = -2;
  }

  // Return now if tree contains no particles
  if (Ntot == 0) return;

  celldata[0].N      = Ntot;
  celldata[0].ifirst = ifirst;
  celldata[0].ilast  = ilast;
  celldata[0].level  = 0;

  // Compute the bounding box of all particles in root cell and the root cell size
  for (k=0; k<ndim; k++) celldata[0].bbmin[k] = +big_number;
  for (k=0; k<ndim; k++) celldata[0].bbmax[k] = -big_number;
  for (i=ifirst; i<=ilast; i++) {
    for (k=0; k<ndim; k++) celldata[0].bbmin[k] = min(celldata[0].bbmin[k],partdata[i].r[k]);
    for (k=0; k<ndim; k++) celldata[0].bbmax[k] = max(celldata[0].bbmax[k],partdata[i].r[k]);
  }
  for (k=0; k<ndim; k++) {
    celldata[0].r[k] = 0.5*(celldata[0].bbmin[k] + celldata[0].bbmax[k]);
    cellSize = max(cellSize, celldata[0].bbmax[k] - celldata[0].rcell[k]);
  }

  //cout << "Bounding box : " << celldata[0].bbmin[0] << "    " << celldata[0].bbmax[0] << endl;

  // Build tree if there are any particles
  //-----------------------------------------------------------------------------------------------
  if (Ntot > 0) {

    celllist = new int[Npartmax];
    whichChild = new int[Npartmax];

    ltot         = 0;
    Nlist        = 1;
    Ncell        = 1;
    celllist[0]  = 0;
    firstCell[0] = 0;
    lastCell[0]  = 0;
    for (i=0; i<Npartmax; i++) whichChild[i] = -1;


    // Recursively divide tree to lower levels until all leaf cells contain maximum no. of ptcls.
    //---------------------------------------------------------------------------------------------
    while (!allDone) {
      cellSize *= (FLOAT) 0.5;

      //cout << "Level : " << ltot << "      Ncell : " << Ncell << "     Nlist : " << Nlist << endl;

      // Loop over all unfinished cells to find new child cell occupancy
      //-------------------------------------------------------------------------------------------
      for (cc=0; cc<Nlist; cc++) {
        c = celllist[cc];
        TreeCell<ndim> &cell = celldata[c];

        i = cell.ifirst;
        ilast = cell.ilast;

        //cout << "Investigating cell : " << c << "    r : " << cell.r[0]
        //     << "    N : " << cell.N << endl;

        // Walk through linked list of all particles to find new child cells
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

          whichChild[i] = ckid;
          cell.childof[ckid] = -1;

          if (i == ilast) break;
          i = inext[i];
        };

      }
      //-------------------------------------------------------------------------------------------

      // Find and create all new child cells (must be done in serial)
      //-------------------------------------------------------------------------------------------
      for (cc=0; cc<Nlist; cc++) {
        c = celllist[cc];
        TreeCell<ndim> &cell = celldata[c];

        for (k=0; k<Noctchild; k++) {
          if (cell.childof[k] == -1) {
            assert(Ncell < Ncellmax);
            cell.childof[k] = Ncell++;
          }
        }
      }
      //-------------------------------------------------------------------------------------------


      // Find number of particles in each new child cell
      //-------------------------------------------------------------------------------------------
      for (cc=0; cc<Nlist; cc++) {
        c = celllist[cc];
        TreeCell<ndim> &cell = celldata[c];
        for (k=0; k<Noctchild; k++) Nincell[k] = 0;
        i = cell.ifirst;
        ilast = cell.ilast;

        // Walk through linked list of all particles to find new child cells
        while (i != -1) {
          ckid = cell.childof[whichChild[i]];
          if (celldata[ckid].ifirst == -1) {
            celldata[ckid].ifirst = i;
          }
          else {
            inext[celldata[ckid].ilast] = i;
          }
          celldata[ckid].ilast = i;
          celldata[ckid].N++;
          Nincell[whichChild[i]]++;
          if (i == ilast) break;
          i = inext[i];
        };

        // Create child cell properties
        Nkids = 0;
        for (k=0; k<Noctchild; k++) {
          if (Nincell[k] == 0) continue;
          ckid = cell.childof[k];
          celldata[ckid].level = ltot + 1;

          // If a leaf cell, flag it and store children ids (NOT NEEDED HERE??)


          // Assign positions of child cells
          if (k == 0 || k == 2 || k == 4 || k == 6) {
            celldata[ckid].r[0] = cell.r[0] - cellSize;
          }
          else {
            celldata[ckid].r[0] = cell.r[0] + cellSize;
          }
          if (ndim > 1) {
            if (k == 0 || k == 1 || k == 4 || k == 5) {
              celldata[ckid].r[1] = cell.r[1] - cellSize;
            }
            else {
              celldata[ckid].r[1] = cell.r[1] + cellSize;
            }
          }
          if (ndim == 3) {
            if (k < 4) {
              celldata[ckid].r[2] = cell.r[2] - cellSize;
            }
            else {
              celldata[ckid].r[2] = cell.r[2] + cellSize;
            }
          }

          //cout << "Creating new cell " << ckid << "   " << k << " at " << celldata[ckid].r[0]
          //     << "  with " << Nincell[k] << " particles" << endl;
          /*if (Nincell[k] > 0 && Nincell[k] <= Nleafmax) {
            cout << "Found new leaf cell : " << ckid << "     N : " << Nincell[k]
                 << "      r : " << celldata[ckid].r[0] << "     ifirst : " << celldata[ckid].ifirst << endl;
          }*/
          clist[Nkids++] = ckid;

        }

        // Set up linked lists from parent to children cells
        cell.copen = clist[0];
        for (k=0; k<Nkids-1; k++) {
          ckid = clist[k];
          celldata[ckid].cnext = clist[k+1];
        }
        celldata[clist[Nkids-1]].cnext = cell.cnext;

      }
      //-------------------------------------------------------------------------------------------

      // Record first and last cells in newly created level
      ltot++;
      firstCell[ltot] = lastCell[ltot-1] + 1;
      lastCell[ltot] = Ncell - 1;

      // Check if all new cells are children cells
      Nlist = 0;
      allDone = true;
      for (c=firstCell[ltot]; c<=lastCell[ltot]; c++) {
        if (celldata[c].N > Nleafmax) {
          allDone = false;
          celllist[Nlist++] = c;
        }
      }

      //cout << "Finished creating level : " << ltot << "     allDone : " << allDone
      //     << "   first/lastCell : " << firstCell[ltot] << "    " << lastCell[ltot] << endl;

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

  //cout << "Finished building oct-tree.  Ncell : " << Ncell << "    " << Ncellmax << endl;

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
  int cfirst,cend;                     // ..
  int i;                               // Particle counter
  int iaux;                            // Aux. particle i.d. variable
  int k;                               // Dimension counter
  int l;                               // ..
  FLOAT dr[ndim];                      // Relative position vector
  FLOAT drsqd;                         // Distance squared
  FLOAT mi;                            // Mass of particle i
  FLOAT p = 0.0;                       // ..
  FLOAT lambda = 0.0;                  // ..

  debug2("[OctTree::StockTree]");


  // Loop over all levels in tree starting from lowest
  //===============================================================================================
  for (l=ltot; l>=0; l--) {


    // Loop over all cells on current level
    //---------------------------------------------------------------------------------------------
    for (c=firstCell[l]; c<=lastCell[l]; c++) {
      TreeCell<ndim> &cell = celldata[c];

      // Zero all summation variables for all cells
      cell.Nactive  = 0;
      //cell.N        = 0;
      cell.m        = 0.0;
      cell.hmax     = 0.0;
      cell.rmax     = 0.0;
      cell.dhmaxdt  = 0.0;
      cell.drmaxdt  = 0.0;
      cell.mac      = 0.0;
      cell.cdistsqd = big_number;
      for (k=0; k<ndim; k++) cell.r[k] = 0.0;
      for (k=0; k<ndim; k++) cell.v[k] = 0.0;
      for (k=0; k<ndim; k++) cell.rcell[k] = 0.0;
      for (k=0; k<ndim; k++) cell.bbmin[k] = big_number;
      for (k=0; k<ndim; k++) cell.bbmax[k] = -big_number;
      for (k=0; k<ndim; k++) cell.hboxmin[k] = big_number;
      for (k=0; k<ndim; k++) cell.hboxmax[k] = -big_number;
      for (k=0; k<5; k++) cell.q[k] = 0.0;


      // If this is a leaf cell, sum over all particles
      //-------------------------------------------------------------------------------------------
      if (cell.copen == -1) {

        // First, check if any particles have been accreted and remove them
        // from the linked list.  If cell no longer contains any live particles,
        // then set N = 0 to ensure cell is not included in future tree-walks.
        i           = cell.ifirst;
        cell.ifirst = -1;
        cell.N      = 0;
        iaux        = -1;
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
            cell.hmax = max(cell.hmax,partdata[i].h);
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
          for (k=0; k<ndim; k++) cell.rcell[k] = 0.5*(cell.bbmin[k] + cell.bbmax[k]);
          for (k=0; k<ndim; k++) dr[k] = 0.5*(cell.bbmax[k] - cell.bbmin[k]);
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
                cell.q[0] += mi*(3.0*dr[0]*dr[0] - drsqd);
                cell.q[1] += mi*3.0*dr[0]*dr[1];
                cell.q[2] += mi*(3.0*dr[1]*dr[1] - drsqd);
                cell.q[3] += mi*3.0*dr[2]*dr[0];
                cell.q[4] += mi*3.0*dr[2]*dr[1];
              }
              else if (ndim == 2) {
                cell.q[0] += mi*(3.0*dr[0]*dr[0] - drsqd);
                cell.q[1] += mi*3.0*dr[0]*dr[1];
                cell.q[2] += mi*(3.0*dr[1]*dr[1] - drsqd);
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
        cfirst = cell.copen;
        cend   = cell.cnext;

        cc = cfirst;
        while (cc != cend) {
          TreeCell<ndim> &child = celldata[cc];

          if (child.N > 0) {
            for (k=0; k<ndim; k++) cell.bbmin[k] = min(child.bbmin[k],cell.bbmin[k]);
            for (k=0; k<ndim; k++) cell.bbmax[k] = max(child.bbmax[k],cell.bbmax[k]);
            for (k=0; k<ndim; k++) cell.hboxmin[k] = min(child.hboxmin[k],cell.hboxmin[k]);
            for (k=0; k<ndim; k++) cell.hboxmax[k] = max(child.hboxmax[k],cell.hboxmax[k]);
            for (k=0; k<ndim; k++) cell.r[k] += child.m*child.r[k];
            for (k=0; k<ndim; k++) cell.v[k] += child.m*child.v[k];
            cell.hmax = max(child.hmax,cell.hmax);
            cell.m += child.m;
          }

          cc = child.cnext;
        };

        if (cell.m > 0.0) {
          for (k=0; k<ndim; k++) cell.r[k] /= cell.m;
          for (k=0; k<ndim; k++) cell.v[k] /= cell.m;
        }
        for (k=0; k<ndim; k++) cell.rcell[k] = 0.5*(cell.bbmin[k] + cell.bbmax[k]);
        for (k=0; k<ndim; k++) dr[k] = 0.5*(cell.bbmax[k] - cell.bbmin[k]);
        cell.cdistsqd = max(DotProduct(dr,dr,ndim),cell.hmax*cell.hmax)/thetamaxsqd;
        cell.rmax = sqrt(DotProduct(dr,dr,ndim));
#ifdef MPI_PARALLEL
        //cell.worktot = child1.worktot + child2.worktot;
#endif

        // Set limits for children (maximum of 8 but may be less)
        cfirst = cell.copen;
        cend   = cell.cnext;

        cc = cfirst;
        while (cc != cend) {
          TreeCell<ndim> &child = celldata[cc];

          // Now add individual quadrupole moment terms
          if (multipole == "quadrupole" && child.N > 0) {
            mi = child.m;
            for (k=0; k<ndim; k++) dr[k] = child.r[k] - cell.r[k];
            drsqd = DotProduct(dr,dr,ndim);
            if (ndim == 3) {
              cell.q[0] += mi*(3.0*dr[0]*dr[0] - drsqd);
              cell.q[1] += mi*3.0*dr[0]*dr[1];
              cell.q[2] += mi*(3.0*dr[1]*dr[1] - drsqd);
              cell.q[3] += mi*3.0*dr[2]*dr[0];
              cell.q[4] += mi*3.0*dr[2]*dr[1];
            }
            else if (ndim == 2) {
              cell.q[0] += mi*(3.0*dr[0]*dr[0] - drsqd);
              cell.q[1] += mi*3.0*dr[0]*dr[1];
              cell.q[2] += mi*(3.0*dr[1]*dr[1] - drsqd);
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
        if (p >= 0.0) cell.mac = 0.0;
        else {
          lambda = 2.0*sqrt(-p/3.0);
          cell.mac = pow(0.5*lambda/macerror,0.66666666666666);
        }
      }
      else
        cell.mac = 0.0;

    }
    //---------------------------------------------------------------------------------------------

  }
  //===============================================================================================

  //cout << "Finished stocking tree" << endl;

  return;
}



//=================================================================================================
//  OctTree::ExtrapolateCellProperties
/// Extrapolate important physical properties of all cells in the tree.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void OctTree<ndim,ParticleType,TreeCell>::ExtrapolateCellProperties
 (FLOAT dt)                            ///< Smallest timestep size
{
  int c;                               // Cell counter
  int k;                               // Dimension counter

  debug2("[OctTree::ExtrapolateCellProperties]");


  // ...
  //-----------------------------------------------------------------------------------------------
  for (c=0; c<Ncell; c++) {

    for (k=0; k<ndim; k++) celldata[c].r[k] += celldata[c].v[k]*dt;
    for (k=0; k<ndim; k++) celldata[c].rcell[k] += celldata[c].v[k]*dt;
    for (k=0; k<ndim; k++) celldata[c].bbmin[k] += celldata[c].v[k]*dt;
    for (k=0; k<ndim; k++) celldata[c].bbmax[k] += celldata[c].v[k]*dt;
    for (k=0; k<ndim; k++) celldata[c].hboxmin[k] += celldata[c].v[k]*dt;
    for (k=0; k<ndim; k++) celldata[c].hboxmax[k] += celldata[c].v[k]*dt;
    //celldata[c].rmax += celldata[c].drmaxdt*dt;
    //celldata[c].hmax += celldata[c].dhmaxdt*dt;

  }
  //-----------------------------------------------------------------------------------------------

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
  //---------------------------------------------------------------------------
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
  //---------------------------------------------------------------------------

  //timing->EndTimingSection("TREE_UPDATE_COUNTERS");

  return;
}



//=================================================================================================
//  OctTree::ComputeActiveParticleList
/// Returns the number (Nactive) and list of ids (activelist) of all active
/// SPH particles in the given cell.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
int OctTree<ndim,ParticleType,TreeCell>::ComputeActiveParticleList
 (TreeCell<ndim> &cell,             ///< [in] Pointer to cell
  ParticleType<ndim> *partdata,        ///< [in] Pointer to particle data array
  int *activelist)                     ///< [out] List of active particles in cell
{
  int i = cell.ifirst;             // Local particle id (set to first ptcl id)
  int ilast = cell.ilast;          // i.d. of last particle in cell c
  int Nactive = 0;                     // No. of active particles in cell

  // Walk through linked list to obtain list and number of active ptcls.
  while (i != -1) {
    if (i < Ntot && partdata[i].active && partdata[i].itype != dead) activelist[Nactive++] = i;
    if (i == ilast) break;
    i = inext[i];
  };

  return Nactive;
}



//=================================================================================================
//  OctTree::BoxOverlap
/// Check if two bounding boxes overlap.  If yes, then returns true.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
bool OctTree<ndim,ParticleType,TreeCell>::BoxOverlap
 (const FLOAT box1min[ndim],           ///< Minimum extent of box 1
  const FLOAT box1max[ndim],           ///< Maximum extent of box 1
  const FLOAT box2min[ndim],           ///< Minimum extent of box 2
  const FLOAT box2max[ndim])           ///< Maximum extent of box 2
{
  if (ndim == 1) {
    if (box1min[0] > box2max[0]) return false;
    if (box2min[0] > box1max[0]) return false;
    return true;
  }
  else if (ndim == 2) {
    if (box1min[0] > box2max[0]) return false;
    if (box2min[0] > box1max[0]) return false;
    if (box1min[1] > box2max[1]) return false;
    if (box2min[1] > box1max[1]) return false;
    return true;
  }
  else if (ndim == 3) {
    if (box1min[0] > box2max[0]) return false;
    if (box2min[0] > box1max[0]) return false;
    if (box1min[1] > box2max[1]) return false;
    if (box2min[1] > box1max[1]) return false;
    if (box1min[2] > box2max[2]) return false;
    if (box2min[2] > box1max[2]) return false;
    return true;
  }

}

//=================================================================================================
//  OctTree::ComputeActiveCellList
/// Returns the number of cells containing active particles, 'Nactive', and
/// the i.d. list of cells contains active particles, 'celllist'
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
int OctTree<ndim,ParticleType,TreeCell>::ComputeActiveCellPointers
 (TreeCell<ndim> **celllist)         ///< Cells id array containing active ptcls
{
  int c;                            // Cell counter
  int Nactive = 0;                  // No. of active leaf cells in tree

  for (c=0; c<Ncell; c++) {
    if (celldata[c].Nactive > 0) {
      celllist[Nactive++] = &(celldata[c]);
      if (celllist[Nactive - 1]->copen != -1) {
        cout << "WTF?? : " << c << "    " << celllist[Nactive - 1]->copen << endl;
        exit(0);
      }
    }
  }

#ifdef MPI_PARALLEL
  for (c=Ncell; c<Ncell+Nimportedcell; c++) {
    if (celldata[c].Nactive > 0) celllist[Nactive++] = &(celldata[c]);
  }
#endif

  return Nactive;
}

//=================================================================================================
//  OctTree::ComputeActiveCellList
/// Returns the number of cells containing active particles, 'Nactive', and
/// the i.d. list of cells contains active particles, 'celllist'
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
int OctTree<ndim,ParticleType,TreeCell>::ComputeActiveCellList
 (TreeCell<ndim> *celllist)         ///< Cells id array containing active ptcls
{
  int c;                            // Cell counter
  int Nactive = 0;                  // No. of active leaf cells in tree

  for (c=0; c<Ncell; c++) {
    if (celldata[c].Nactive > 0) {
      celllist[Nactive++] = celldata[c];
      if (celllist[Nactive - 1].copen != -1) {
        cout << "WTF?? : " << c << "    " << celllist[Nactive - 1].copen << endl;
        exit(0);
      }
    }
  }

#ifdef MPI_PARALLEL
  for (c=Ncell; c<Ncell+Nimportedcell; c++) {
    if (celldata[c].Nactive > 0) celllist[Nactive++] = celldata[c];
  }
#endif

  return Nactive;
}



//=================================================================================================
//  OctTree::ComputeGatherNeighbourList
/// Computes and returns number of neighbour, 'Nneib', and the list
/// of neighbour ids, 'neiblist', for all particles inside cell 'c'.
/// Includes all particles in the selected cell, plus all particles
/// contained in adjacent cells (including diagonal cells).
/// Wrapper around the true implementation inside OctTree
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
int OctTree<ndim,ParticleType,TreeCell>::ComputeGatherNeighbourList
 (const ParticleType<ndim> *partdata,  ///< [in] Particle data array
  const FLOAT rp[ndim],                ///< [in] Search position
  const FLOAT rsearch,                 ///< [in] Maximum smoothing length
  const int Nneibmax,                  ///< [in] Max. no. of neighbours
  int *neiblist)                       ///< [out] List of neighbour i.d.s
{
  int cc;                              // Cell counter
  int i;                               // Particle id
  int k;                               // Neighbour counter
  int Nneib = 0;                       // Neighbour counter
  FLOAT dr[ndim];                      // Relative position vector
  FLOAT drsqd;                         // Distance squared
  FLOAT rsearchsqd;                    // Search radius squared

  // Start with root cell and walk through entire tree
  cc = 0;
  rsearchsqd = rsearch*rsearch;

  //===============================================================================================
  while (cc < Ncell) {

    for (k=0; k<ndim; k++) dr[k] = celldata[cc].rcell[k] - rp[k];
    drsqd = DotProduct(dr,dr,ndim);


    // Check if bounding boxes overlap with each other
    //---------------------------------------------------------------------------------------------
    if (drsqd < (rsearch + celldata[cc].rmax)*(rsearch + celldata[cc].rmax)) {

      // If not a leaf-cell, then open cell to first child cell
      if (celldata[cc].copen != -1)
        cc = celldata[cc].copen;

      else if (celldata[cc].N == 0)
        cc = celldata[cc].cnext;

      // If leaf-cell, add particles to list
      else if (celldata[cc].copen == -1 && Nneib + Nleafmax < Nneibmax) {
        i = celldata[cc].ifirst;
        while (i != -1) {
          for (k=0; k<ndim; k++) dr[k] = partdata[i].r[k] - rp[k];
          drsqd = DotProduct(dr,dr,ndim);
          if (drsqd < rsearchsqd && partdata[i].itype != dead) neiblist[Nneib++] = i;
          if (i == celldata[cc].ilast) break;
          i = inext[i];
        };
        cc = celldata[cc].cnext;
      }

      // If leaf-cell, but we've run out of memory, return with error-code (-1)
      else if (celldata[cc].copen == -1 && Nneib + Nleafmax >= Nneibmax)
        return -1;

    }

    // If not in range, then open next cell
    //---------------------------------------------------------------------------------------------
    else
      cc = celldata[cc].cnext;

  };
  //===============================================================================================


  return Nneib;
}



//=================================================================================================
//  OctTree::ComputeGatherNeighbourList
/// Computes and returns number of neighbour, 'Nneib', and the list
/// of neighbour ids, 'neiblist', for all particles inside cell 'c'.
/// Includes all particles in the selected cell, plus all particles
/// contained in adjacent cells (including diagonal cells).
/// Wrapper around the true implementation inside OctTree
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
int OctTree<ndim,ParticleType,TreeCell>::ComputeGatherNeighbourList
 (const TreeCell<ndim> &cell,          ///< [in] Pointer to current cell
  const ParticleType<ndim> *partdata,  ///< [in] Particle data array
  const FLOAT hmax,                    ///< [in] Maximum smoothing length
  const int Nneibmax,                  ///< [in] Max. no. of neighbours
  int &Nneib,                          ///< [inout] No. of neighbours
  int *neiblist)                       ///< [out] List of neighbour i.d.s
{
  int cc = 0;                          // Cell counter
  int i;                               // Particle id
  int j;                               // Aux. particle counter
  int k;                               // Neighbour counter
  int Ntemp = Nneib;                   // Temporary neighbour counter
  FLOAT dr[ndim];                      // Relative position vector
  FLOAT drsqd;                         // Distance squared
  FLOAT gatherboxmin[ndim];            // Minimum gather neighbour box
  FLOAT gatherboxmax[ndim];            // Maximum gather neighbour box
  FLOAT rc[ndim];                      // Position of cell

  // Exit immediately if we have overflowed the neighbour list buffer
  if (Nneib == -1) return -1;

  const FLOAT hrangemaxsqd = pow(cell.rmax + kernrange*hmax,2);
  for (k=0; k<ndim; k++) rc[k] = cell.rcell[k];
  for (k=0; k<ndim; k++) gatherboxmin[k] = cell.bbmin[k] - kernrange*hmax;
  for (k=0; k<ndim; k++) gatherboxmax[k] = cell.bbmax[k] + kernrange*hmax;


  // Start with root cell and walk through entire tree
  //===============================================================================================
  while (cc < Ncell) {

    // Check if bounding boxes overlap with each other
    //---------------------------------------------------------------------------------------------
    if (BoxOverlap(gatherboxmin,gatherboxmax,celldata[cc].bbmin,celldata[cc].bbmax)) {

      // If not a leaf-cell, then open cell to first child cell
      if (celldata[cc].copen != -1)
        cc = celldata[cc].copen;

      else if (celldata[cc].N == 0)
        cc = celldata[cc].cnext;

      // If leaf-cell, add particles to list
      else if (celldata[cc].copen == -1 && Ntemp + Nleafmax < Nneibmax) {
        i = celldata[cc].ifirst;
        while (i != -1) {
          neiblist[Ntemp++] = i;
          if (i == celldata[cc].ilast) break;
          i = inext[i];
        };
        cc = celldata[cc].cnext;
      }

      // If leaf-cell, but we've run out of memory, return with error-code (-1)
      else if (celldata[cc].copen == -1 && Ntemp + Nleafmax >= Nneibmax)
        return -1;

    }

    // If not in range, then open next cell
    //---------------------------------------------------------------------------------------------
    else
      cc = celldata[cc].cnext;

  };
  //===============================================================================================


  // Now, trim the list to remove particles that are definitely not neighbours
  for (j=Nneib; j<Ntemp; j++) {
    i = neiblist[j];
    if (partdata[i].itype == dead) continue;
    for (k=0; k<ndim; k++) dr[k] = partdata[i].r[k] - rc[k];
    drsqd = DotProduct(dr,dr,ndim);
    if (drsqd < hrangemaxsqd) neiblist[Nneib++] = i;
  }

  //cout << "Returning " << Nneib << " neighbours" << endl;

  return Nneib;
}



//=================================================================================================
//  OctTree::ComputeNeighbourList
/// Computes and returns number of neighbour, 'Nneib', and the list
/// of neighbour ids, 'neiblist', for all particles inside cell 'c'.
/// Includes all particles in the selected cell, plus all particles
/// contained in adjacent cells (including diagonal cells).
/// Wrapper around the true implementation inside OctTree.
/// If allocated memory array containing neighbour ids (neiblist) overflows,
/// return with error code (-1) in order to reallocate more memory.
//================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
int OctTree<ndim,ParticleType,TreeCell>::ComputeNeighbourList
 (const TreeCell<ndim> &cell,          ///< [in] Cell pointer
  const ParticleType<ndim> *partdata,  ///< [in] Particle data array
  const int Nneibmax,                  ///< [in] Max. no. of neighbours
  int &Nneib,                          ///< [inout] No. of neighbours
  int *neiblist,                       ///< [out] List of neighbour i.d.s
  ParticleType<ndim> *neibpart)        ///< [out] Array of local copies of neighbours
{
  int cc = 0;                          // Cell counter
  int i;                               // Particle id
  int j;                               // Aux. particle counter
  int k;                               // Neighbour counter
  int Ntemp = Nneib;                   // Temp. no. of neighbouts
  FLOAT dr[ndim];                      // Relative position vector
  FLOAT drsqd;                         // Distance squared
  FLOAT rc[ndim];                      // Position of cell

  // Local (const) copies of cell properties
  const FLOAT hrangemaxsqd = pow(cell.rmax + kernrange*cell.hmax,2);
  const FLOAT rmax = cell.rmax;
  for (k=0; k<ndim; k++) rc[k] = cell.rcell[k];

  // Exit immediately if we have overflowed the neighbour list buffer
  if (Nneib == -1) return -1;


  // Start with root cell and walk through entire tree
  //===============================================================================================
  while (cc < Ncell) {

    for (k=0; k<ndim; k++) dr[k] = celldata[cc].rcell[k] - rc[k];
    drsqd = DotProduct(dr,dr,ndim);

    // Check if bounding boxes overlap with each other
    //---------------------------------------------------------------------------------------------
    if (BoxOverlap(cell.bbmin,cell.bbmax,celldata[cc].hboxmin,celldata[cc].hboxmax) ||
        BoxOverlap(cell.hboxmin,cell.hboxmax,celldata[cc].bbmin,celldata[cc].bbmax)) {

      // If not a leaf-cell, then open cell to first child cell
      if (celldata[cc].copen != -1)
        cc = celldata[cc].copen;

      else if (celldata[cc].N == 0)
        cc = celldata[cc].cnext;

      // If leaf-cell, add particles to list
      else if (celldata[cc].copen == -1 && Ntemp + Nleafmax < Nneibmax) {
        i = celldata[cc].ifirst;
        while (i != -1) {
          //neibpart[Ntemp] = partdata[i];
          neiblist[Ntemp++] = i;
          if (i == celldata[cc].ilast) break;
          i = inext[i];
        };
        cc = celldata[cc].cnext;
      }

      // If leaf-cell, but we've run out of memory, return with error-code (-1)
      else if (celldata[cc].copen  == -1 && Ntemp + Nleafmax >= Nneibmax)
        return -1;

    }

    // If not in range, then open next cell
    //---------------------------------------------------------------------------------------------
    else
      cc = celldata[cc].cnext;

  };
  //===============================================================================================


  // Now, trim the list to remove particles that are definitely not neighbours
  for (j=Nneib; j<Ntemp; j++) {
    i = neiblist[j];
    if (partdata[i].itype == dead) continue;

    for (k=0; k<ndim; k++) dr[k] = partdata[i].r[k] - rc[k];
    drsqd = DotProduct(dr,dr,ndim);
    if (drsqd < hrangemaxsqd || drsqd <
        (rmax + kernrange*partdata[i].h)*(rmax + kernrange*partdata[i].h)) {
      neibpart[Nneib] = partdata[i];
      neiblist[Nneib] = i;
      Nneib++;
    }
  }


  return Nneib;
}



//=================================================================================================
//  OctTree::ComputeGravityInteractionList
/// Computes and returns number of SPH neighbours (Nneib), direct sum particles
/// (Ndirect) and number of cells (Ngravcell), including lists of ids, from
/// the gravity tree walk for active particles inside cell c.
/// Currently defaults to the geometric opening criteria.
/// If any of the interactions list arrays (neiblist,directlist,gravcelllist)
/// overflow, return with error code (-1) to reallocate more memory.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
int OctTree<ndim,ParticleType,TreeCell>::ComputeGravityInteractionList
 (const TreeCell<ndim> &cell,          ///< [in] Pointer to cell
  const ParticleType<ndim> *partdata,  ///< [in] Particle data array
  const FLOAT macfactor,               ///< [in] Gravity MAC particle factor
  const int Nneibmax,                  ///< [in] Max. no. of SPH neighbours
  const int Ngravcellmax,              ///< [in] Max. no. of cell interactions
  int &Nneib,                          ///< [out] Total no. of neighbours
  int &Nhydroneib,                       ///< [out] No. of SPH neighbours
  int &Ndirect,                        ///< [out] No. of direct-sum neighbours
  int &Ngravcell,                      ///< [out] No. of cell interactions
  int *neiblist,                       ///< [out] List of all particle ids
  int *sphneiblist,                    ///< [out] List of SPH neibpart ids
  int *directlist,                     ///< [out] List of direct-sum neibpart ids
  TreeCell<ndim> *gravcell,            ///< [out] Array of local copies of tree cells
  ParticleType<ndim> *neibpart)        ///< [out] Array of local copies of neighbour particles
{
  int cc = 0;                          // Cell counter
  int i;                               // Particle id
  int j;                               // Aux. particle counter
  int k;                               // Neighbour counter
  int Nhydroneibtemp = 0;                // Aux. counter
  FLOAT dr[ndim];                      // Relative position vector
  FLOAT drsqd;                         // Distance squared
  FLOAT rc[ndim];                      // Position of cell

  // Make local copies of important cell properties
  const FLOAT hrangemaxsqd = pow(cell.rmax + kernrange*cell.hmax,2);
  const FLOAT rmax = cell.rmax;
  for (k=0; k<ndim; k++) rc[k] = cell.rcell[k];

  // Start with root cell and walk through entire tree
  Nneib     = 0;
  Nhydroneib  = 0;
  Ndirect   = 0;
  Ngravcell = 0;


  // Walk through all cells in tree to determine particle and cell interaction lists
  //===============================================================================================
  while (cc < Ncell) {

    for (k=0; k<ndim; k++) dr[k] = celldata[cc].rcell[k] - rc[k];
    drsqd = DotProduct(dr,dr,ndim);


    // Check if bounding boxes overlap with each other
    //---------------------------------------------------------------------------------------------
    if (BoxOverlap(cell.bbmin,cell.bbmax,celldata[cc].hboxmin,celldata[cc].hboxmax) ||
        BoxOverlap(cell.hboxmin,cell.hboxmax,celldata[cc].bbmin,celldata[cc].bbmax)) {

      // If not a leaf-cell, then open cell to first child cell
      if (celldata[cc].copen != -1)
        cc = celldata[cc].copen;

      else if (celldata[cc].N == 0)
        cc = celldata[cc].cnext;

      // If leaf-cell, add particles to list
      else if (celldata[cc].copen == -1 && Nneib + Nleafmax <= Nneibmax) {
        i = celldata[cc].ifirst;
        while (i != -1) {
          sphneiblist[Nhydroneib++] = Nneib;
          neiblist[Nneib] = i;
          neibpart[Nneib++] = partdata[i];
          if (i == celldata[cc].ilast) break;
          i = inext[i];
        };
        cc = celldata[cc].cnext;
      }

      // If leaf-cell, but we've run out of memory, return with error-code (-1)
      else if (celldata[cc].copen == -1 && Nneib + Nleafmax > Nneibmax)
        return -1;

    }

    // Check if cell is far enough away to use the COM approximation
    //---------------------------------------------------------------------------------------------
    else if (drsqd > celldata[cc].cdistsqd && drsqd > celldata[cc].mac*macfactor &&
              celldata[cc].N > 0) {

      // If cell is a leaf-cell with only one particle, more efficient to
      // compute the gravitational contribution from the particle than the cell
      if (celldata[cc].copen == -1 && celldata[cc].N == 1 && Nneib + 1 <= Nneibmax) {
        i = celldata[cc].ifirst;
        directlist[Ndirect++] = Nneib;
        neiblist[Nneib] = i;
        neibpart[Nneib++] = partdata[i];
      }
      else if (Ngravcell < Ngravcellmax) {
        gravcell[Ngravcell++] = celldata[cc];
      }
      else {
        return -2;
      }
      cc = celldata[cc].cnext;

    }

    // If cell is too close, open cell to interogate children cells.
    // If cell is too close and a leaf cell, then add particles to direct list.
    //---------------------------------------------------------------------------------------------
    else if (!(drsqd > celldata[cc].cdistsqd && drsqd > celldata[cc].mac*macfactor) &&
              celldata[cc].N > 0) {

      // If not a leaf-cell, then open cell to first child cell
      if (celldata[cc].copen != -1)
        cc = celldata[cc].copen;

      // If leaf-cell, add particles to list
      else if (celldata[cc].copen == -1 && Nneib + Nleafmax <= Nneibmax) {
        i = celldata[cc].ifirst;
        while (i != -1) {
          directlist[Ndirect++] = Nneib;
          neiblist[Nneib] = i;
          neibpart[Nneib++] = partdata[i];
          if (i == celldata[cc].ilast) break;
          i = inext[i];
        };
        cc = celldata[cc].cnext;
      }

      // If leaf-cell, but we've run out of memory, return with error-code (-1)
      else if (celldata[cc].copen == -1 && Nneib + Nleafmax > Nneibmax)
        return -3;

    }

    // If not in range, then open next cell
    //---------------------------------------------------------------------------------------------
    else
      cc = celldata[cc].cnext;

  };
  //===============================================================================================


  // Now, trim the list to remove particles that are definitely not SPH neighbours.
  // If not an SPH neighbour, then add to direct gravity sum list.
  for (j=Nhydroneibtemp; j<Nhydroneib; j++) {
    i = sphneiblist[j];
    if (neibpart[i].itype == dead) continue;
    for (k=0; k<ndim; k++) dr[k] = neibpart[i].r[k] - rc[k];
    drsqd = DotProduct(dr,dr,ndim);
    if (drsqd < hrangemaxsqd || drsqd <
        (rmax + kernrange*neibpart[i].h)*(rmax + kernrange*neibpart[i].h)) {
      sphneiblist[Nhydroneibtemp++] = i;
    }
    else if (Ndirect < Nneibmax) {
      directlist[Ndirect++] = i;
    }
    else {
      return -3;
    }
  }
  Nhydroneib = Nhydroneibtemp;

  return 1;
}



//=================================================================================================
//  OctTree::ComputeStarGravityInteractionList
/// Computes and returns number of SPH neighbours (Nneib), direct sum particles
/// (Ndirect) and number of cells (Ngravcell), including lists of ids, from
/// the gravity tree walk for active particles inside cell c.
/// Currently defaults to the geometric opening criteria.
/// If any of the interactions list arrays (neiblist,directlist,gravcelllist)
/// overflow, return with error code (-1) to reallocate more memory.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
int OctTree<ndim,ParticleType,TreeCell>::ComputeStarGravityInteractionList
 (const NbodyParticle<ndim> *star,     ///< [in] Nbody particle
  const FLOAT macfactor,               ///< [in] Gravity MAC factor
  const int Nneibmax,                  ///< [in] Max. no. of SPH neighbours
  const int Ndirectmax,                ///< [in] Max. no. of direct-sum neighbours
  const int Ngravcellmax,              ///< [in] Max. no. of cell interactions
  int &Nneib,                          ///< [out] No. of SPH neighbours
  int &Ndirect,                        ///< [out] No. of direct-sum neighbours
  int &Ngravcell,                      ///< [out] No. of cell interactions
  int *neiblist,                       ///< [out] List of SPH neighbour ids
  int *directlist,                     ///< [out] List of direct-sum neighbour ids
  TreeCell<ndim> *gravcell,            ///< [out] List of cell ids
  ParticleType<ndim> *partdata)        ///< [in] Particle data array
{
  int cc;                              // Cell counter
  int i;                               // Particle id
  int k;                               // Neighbour counter
  FLOAT dr[ndim];                      // Relative position vector
  FLOAT drsqd;                         // Distance squared
  FLOAT hrangemax;                     // Maximum kernel extent
  FLOAT rs[ndim];                      // Position of star

  // Make local copies of important cell properties
  for (k=0; k<ndim; k++) rs[k] = star->r[k];
  hrangemax = kernrange*star->h;

  // Start with root cell and walk through entire tree
  cc = 0;
  Nneib = 0;
  Ndirect = 0;
  Ngravcell = 0;


  // Walk through all cells in tree to determine particle and cell
  // interaction lists
  //===============================================================================================
  while (cc < Ncell) {

    for (k=0; k<ndim; k++) dr[k] = celldata[cc].rcell[k] - rs[k];
    drsqd = DotProduct(dr,dr,ndim);


    // Check if cells contain SPH neighbours
    //---------------------------------------------------------------------------------------------
    if (drsqd < pow(0.5*hrangemax + celldata[cc].rmax + 0.5*kernrange*celldata[cc].hmax,2)) {

      // If not a leaf-cell, then open cell to first child cell
      if (celldata[cc].copen != -1)
        cc = celldata[cc].copen;

      // If leaf-cell, add particles to list
      else if (celldata[cc].copen == -1 && Nneib + Nleafmax <= Nneibmax) {
        i = celldata[cc].ifirst;
        while (i != -1) {
          if (partdata[i].itype != dead) neiblist[Nneib++] = i;
          if (i == celldata[cc].ilast) break;
          i = inext[i];
        };
        cc = celldata[cc].cnext;
      }

      // If leaf-cell, but we've run out of memory, return with error-code (-1)
      else if (celldata[cc].copen == -1 && Nneib + Nleafmax > Nneibmax)
        return -1;

    }

    // Check if cell is far enough away to use the COM approximation
    //---------------------------------------------------------------------------------------------
    else if (drsqd > celldata[cc].cdistsqd && celldata[cc].N > 0) {

      // If cell is a leaf-cell with only one particle, more efficient to
      // compute the gravitational contribution from the particle than the cell
      if (celldata[cc].copen == -1 && celldata[cc].N == 1 && Ndirect + Nneib < Ndirectmax) {
        if (partdata[celldata[cc].ifirst].itype != dead)
        directlist[Ndirect++] = celldata[cc].ifirst;
      }
      else if (Ngravcell < Ngravcellmax && celldata[cc].N > 0)
        gravcell[Ngravcell++] = celldata[cc];
      else
        return -1;
      cc = celldata[cc].cnext;

    }

    // If cell is too close, open cell to interogate children cells.
    // If cell is too close and a leaf cell, then add particles to direct list.
    //---------------------------------------------------------------------------------------------
    else if (drsqd <= celldata[cc].cdistsqd && celldata[cc].N > 0) {

      // If not a leaf-cell, then open cell to first child cell
      if (celldata[cc].copen != -1)
        cc = celldata[cc].copen;

      // If leaf-cell, add particles to list
      else if (celldata[cc].copen == -1 && Ndirect + Nleafmax <= Ndirectmax) {
        i = celldata[cc].ifirst;
        while (i != -1) {
          if (partdata[i].itype != dead) directlist[Ndirect++] = i;
          if (i == celldata[cc].ilast) break;
          i = inext[i];
        };
        cc = celldata[cc].cnext;
      }

      // If leaf-cell, but we've run out of memory, return with error-code (-1)
      else if (celldata[cc].copen == -1 && Ndirect + Nleafmax > Ndirectmax)
        return -1;

    }

    // If not in range, then open next cell
    //---------------------------------------------------------------------------------------------
    else
      cc = celldata[cc].cnext;

  };
  //===============================================================================================


  return 1;
}



#ifdef MPI_PARALLEL
//=================================================================================================
//  OctTree::ComputeDistantGravityInteractionList
/// Computes and returns number of SPH neighbours (Nneib), direct sum particles
/// (Ndirect) and number of cells (Ngravcell), including lists of ids, from
/// the gravity tree walk for active particles inside cell c.
/// Currently defaults to the geometric opening criteria.
/// If any of the interactions list arrays (neiblist,directlist,gravcelllist)
/// overflow, return with error code (-1) to reallocate more memory.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
int OctTree<ndim,ParticleType,TreeCell>::ComputeDistantGravityInteractionList
 (const TreeCell<ndim> *cellptr,       ///< [in] Pointer to cell
  const FLOAT macfactor,               ///< [in] Gravity MAC particle factor
  const int Ngravcellmax,              ///< [in] Max. no. of cell interactions
  int Ngravcell,                       ///< [in] Current no. of cells in array
  TreeCell<ndim> *gravcelllist)       ///< [out] Array of cells
{
  int cc;                              // Cell counter
  int i;                               // Particle id
  int j;                               // Aux. particle counter
  int k;                               // Neighbour counter
  int Ngravcelltemp = Ngravcell;       // ..
  FLOAT dr[ndim];                      // Relative position vector
  FLOAT drsqd;                         // Distance squared
  FLOAT rc[ndim];                      // Position of cell
  FLOAT hrangemax;                     // Maximum kernel extent
  FLOAT rmax;                          // Radius of sphere containing particles

  // Make local copies of important cell properties
  for (k=0; k<ndim; k++) rc[k] = cellptr->rcell[k];
  hrangemax = cellptr->rmax + kernrange*cellptr->hmax;
  rmax = cellptr->rmax;

  // Start with root cell and walk through entire tree
  cc = 0;


  // Walk through all cells in tree to determine particle and cell
  // interaction lists
  //===========================================================================
  while (cc < Ncell) {

    for (k=0; k<ndim; k++) dr[k] = celldata[cc].rcell[k] - rc[k];
    drsqd = DotProduct(dr,dr,ndim);


    // Check if bounding boxes overlap with each other
    //-------------------------------------------------------------------------
    if (BoxOverlap(cellptr->bbmin,cellptr->bbmax,
                   celldata[cc].hboxmin,celldata[cc].hboxmax) ||
        BoxOverlap(cellptr->hboxmin,cellptr->hboxmax,
                   celldata[cc].bbmin,celldata[cc].bbmax)) {

      // If not a leaf-cell, then open cell to first child cell
      if (celldata[cc].level != ltot)
        cc++;

      else if (celldata[cc].N == 0)
        cc = celldata[cc].cnext;

      // If leaf-cell, add particles to list
      else if (celldata[cc].level == ltot) {
        return -1;
      }

    }

    // Check if cell is far enough away to use the COM approximation
    //-------------------------------------------------------------------------
    else if (drsqd > celldata[cc].cdistsqd && drsqd > celldata[cc].mac*macfactor
	     && celldata[cc].N > 0) {

      gravcelllist[Ngravcelltemp++] = celldata[cc];
      cc = celldata[cc].cnext;

    }

    // If cell is too close, open cell to interogate children cells.
    // If cell is too close and a leaf cell, then add particles to direct list.
    //-------------------------------------------------------------------------
    else if (!(drsqd > celldata[cc].cdistsqd &&
	       drsqd > celldata[cc].mac*macfactor) && celldata[cc].N > 0) {

      // If not a leaf-cell, then open cell to first child cell
      if (celldata[cc].level != ltot)
         cc++;

      // If leaf-cell, add particles to list
      else
       	return -1;

    }

    // If not in range, then open next cell
    //-------------------------------------------------------------------------
    else
      cc = celldata[cc].cnext;

  };
  //===========================================================================

  return Ngravcelltemp;
}



//=================================================================================================
//  OctTree::ComputeHydroTreeCellOverlap
/// ...
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
bool OctTree<ndim,ParticleType,TreeCell>::ComputeHydroTreeCellOverlap
 (const TreeCell<ndim> *cellptr)       ///< [in] Pointer to cell
{
  int cc = 0;                          // Cell counter

  // Walk through all cells in tree to determine particle and cell interaction lists
  //===============================================================================================
  while (cc < Ncell) {

    // Check if bounding boxes overlap with each other
    //---------------------------------------------------------------------------------------------
    if (BoxOverlap(cellptr->bbmin,cellptr->bbmax,celldata[cc].hboxmin,celldata[cc].hboxmax) ||
        BoxOverlap(cellptr->hboxmin,cellptr->hboxmax,celldata[cc].bbmin,celldata[cc].bbmax)) {

      // If not a leaf-cell, then open cell to first child cell
      if (celldata[cc].level != ltot) cc++;

      // If cell contains no particle (dead leaf?) then move to next cell
      else if (celldata[cc].N == 0) cc = celldata[cc].cnext;

      // If cell is overlapping with any leaf, then flag overlap on return
      else if (celldata[cc].level == ltot) return true;

    }

    // If not in range, then open next cell
    //---------------------------------------------------------------------------------------------
    else
      cc = celldata[cc].cnext;

  };
  //===============================================================================================

  // If we've walked the entire tree wihout any leaf overlaps, flag no overlap
  return false;
}
#endif



#if defined(VERIFY_ALL)
//=================================================================================================
//  OctTree::ValidateTree
/// Performs various sanity and validation checks on KD-tree structure.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void OctTree<ndim,ParticleType,TreeCell>::ValidateTree
 (ParticleType<ndim> *partdata)        ///< Pointer to SPH class
{
  bool kill_flag = false;              // ..
  int activecount;                     // Active particles in leaf cell
  int c;                               // Cell counter
  int i;                               // Particle counter
  int k;                               // ..
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



template class OctTree<1,Particle,OctTreeCell>;
template class OctTree<2,Particle,OctTreeCell>;
template class OctTree<3,Particle,OctTreeCell>;
template class OctTree<1,SphParticle,OctTreeCell>;
template class OctTree<2,SphParticle,OctTreeCell>;
template class OctTree<3,SphParticle,OctTreeCell>;
template class OctTree<1,GradhSphParticle,OctTreeCell>;
template class OctTree<2,GradhSphParticle,OctTreeCell>;
template class OctTree<3,GradhSphParticle,OctTreeCell>;
template class OctTree<1,SM2012SphParticle,OctTreeCell>;
template class OctTree<2,SM2012SphParticle,OctTreeCell>;
template class OctTree<3,SM2012SphParticle,OctTreeCell>;
template class OctTree<1,GodunovSphParticle,OctTreeCell>;
template class OctTree<2,GodunovSphParticle,OctTreeCell>;
template class OctTree<3,GodunovSphParticle,OctTreeCell>;
