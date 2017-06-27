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


#include <algorithm>
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
                                             string gravity_mac_aux, multipole_method multipole_aux,
                                             const DomainBox<ndim>& domain,
                                    		 const ParticleTypeRegister& reg,
											 const bool IAmPruned):
  Tree<ndim,ParticleType,TreeCell>(Nleafmaxaux, thetamaxsqdaux, kernrangeaux,
                                   macerroraux, gravity_mac_aux, multipole_aux, domain, reg, IAmPruned)
{
  allocated_tree = false;
  ifirst         = -1;
  ilast          = -1;
  lmax           = 40;
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
void OctTree<ndim,ParticleType,TreeCell>::AllocateTreeMemory(int Nparticles, int Ncells, bool force_realloc)
{
  debug2("[OctTree::AllocateTreeMemory]");

  assert(Ntot >= 0);
  assert(lmax >= 0);

  if (!allocated_tree || Nparticles > Ntotmax || Ncells > Ncellmax || force_realloc) {
    if (allocated_tree) DeallocateTreeMemory();

    Ncells = max(Ncells,Ncellmax);
    Nparticles = max(Nparticles,Ntotmax);
    Ncells    = max((int) ((FLOAT) 2.0*(FLOAT) Ncells), 4*Nparticles);
    gtot        = Ntotmax;

    firstCell = new int[lmax];
    lastCell  = new int[lmax];
    celldata  = new struct TreeCell<ndim>[Ncells];

    allocated_tree = true;

    Ntotmax = Nparticles;
    Ncellmax = Ncells;

  }

  return;
}

//=================================================================================================
//  OctTree::ReallocateMemory
/// Reallocate memory for OctTree (when we need to grow the tree) as requested. Preserves the existing
/// information in the tree, differently from the previous function
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void OctTree<ndim,ParticleType,TreeCell>::ReallocateMemory(int Nparticles, int Ncells)
{
  debug2("[OctTree::ReallocateMemory]");

  if (!allocated_tree) {
	  ExceptionHandler::getIstance().raise("This function should not be called if the tree has not been allocated yet!");
  }


  if (Nparticles > Ntotmax ) {

	  Ntotmax = Nparticles;

  }

  if (Ncells > Ncellmax) {


    TreeCell<ndim>* celldataold = celldata;

    celldata = new struct TreeCell<ndim>[Ncells];

    std::copy(celldataold,celldataold+Ncellmax,celldata);

    delete[] celldataold;

    Ncellmax = Ncells;

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
    delete[] lastCell;
    delete[] firstCell;
    allocated_tree = false;
  }

  return;
}


// Predicate for comparing particles and sorting particles
template<class ParticleType>
class ParticleSorter {

public:
  ParticleSorter(int jdir, FLOAT r)
    : _j(jdir), _r(r) { } ;

  bool operator()(const ParticleType& p) const {
    return p.r[_j] < _r ;
  }
private:
  int _j ;
  FLOAT _r ;
};



//=================================================================================================
//  OctTree::DivideCell
/// Sort the particles into the two sub-cells in the given direction.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void OctTree<ndim,ParticleType,TreeCell>::DivideCell
(int c,                                      ///< [in] Cell to Divide
 int first,                                  ///< [in] Id of first particle
 int last,                                   ///< [in[ ID of last particle
 int cchild,                                 ///< [in] Location of children in the cell array
 ParticleType<ndim>* partdata,               ///< [in] Particle data array
 int jdir)                                   ///< [in] Direction to split
 {
  // First Sort the particles across the range.
  ParticleSorter<ParticleType<ndim> > pred(jdir, celldata[c].rcentre[jdir]) ;
  ParticleType<ndim>* join = std::partition(partdata+first, partdata+last+1, pred);
  int njoin = std::distance(partdata+first, join) ;


  //  Work out where to place the children
  int child1 = cchild;
  int child2 = cchild + (1 << jdir);

  // Split the sub-ranges, or set the particle ranges
  if (jdir > 0) {
    DivideCell(c, first, first+njoin-1, child1, partdata, jdir-1) ;
    DivideCell(c, first+njoin, last,    child2, partdata, jdir-1) ;
  } else {

    if (njoin > 0) {
      celldata[child1].ifirst = first ;
      celldata[child1].ilast  = first + njoin - 1;
      celldata[child1].N      = njoin ;
    }
    else {
      celldata[child1].ifirst = -1 ;
      celldata[child1].ilast  = -1 ;
      celldata[child1].N      = 0;
    }

    if (first + njoin <= last) {
      celldata[child2].ifirst = first + njoin ;
      celldata[child2].ilast  = last ;
      celldata[child2].N      = last+1 - (first+njoin);
    }
    else {
      celldata[child1].ifirst = -1;
      celldata[child1].ilast  = -1;
      celldata[child1].N      = 0;
    }
  }

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
  Particle<ndim> *part_gen)            ///< [in] Particle data array
{
  bool allDone = false;                // Are all cell divisions completed?
  int c;                               // Cell counter
  int cc;                              // Child cell counter
  int cnew;                            // i.d. of newly created cell
  int ilast;                           // i.d. of last particle in cell
  int i;                               // Particle counter
  int k;                               // Dimension counter
  int kk;                              // ..
  int Nlist;                           // ..
  int *celllist;                       // List of cells to be processed
  FLOAT cellSize = (FLOAT) 0.0;        // Size of new cell (from centre to edge)
  FLOAT bbmin[ndim];                   // Minimum extent of local bounding box
  FLOAT bbmax[ndim];                   // Maximum extent of local bounding box
  ParticleType<ndim>* partdata = reinterpret_cast<ParticleType<ndim>*>(part_gen);

  debug2("[OctTree::BuildTree]");
  //timing->StartTimingSection("BUILD_OCT_TREE");

  // Allocate (or reallocate if needed) all tree memory
  AllocateTreeMemory(Npartmax,0,false);

  for (c=0; c<Ncellmax; c++) {
    celldata[c].N      = 0;
    celldata[c].ifirst = -1;
    celldata[c].ilast  = -1;
    celldata[c].copen  = -1;
    celldata[c].cnext  = Ncellmax;
    celldata[c].id     = c;
  }

  // Create bounding box of SPH particles
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


  // Set properties for root cell before constructing tree
  Ncell  = 1;
  ltot   = 0;
  celldata[0].N      = _ilast - _ifirst + 1;
  celldata[0].ifirst = _ifirst;
  celldata[0].ilast  = _ilast;
  celldata[0].level  = 0;
  celldata[0].copen  = -1;
  celldata[0].hmax = 0;
  for (k=0; k<ndim; k++) celldata[0].bb.max[k] = bbmax[k];
  for (k=0; k<ndim; k++) celldata[0].bb.min[k] = bbmin[k];
  celldata[0].parent = -1;
  for (k=0; k<ndim; k++) {
    celldata[0].rcentre[k] = (FLOAT) 0.5*(celldata[0].bb.min[k] + celldata[0].bb.max[k]);
    cellSize = max(cellSize, celldata[0].bb.max[k] - celldata[0].rcentre[k]);
  }
  rootCellSize = (FLOAT) 2.0*cellSize;


  // Build tree if there are any particles
  //-----------------------------------------------------------------------------------------------
  if (Ntot > 0) {

    celllist     = new int[Npartmax];
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
          celldata[cnew].parent = cc;
          for (kk=0; kk<ndim; kk++) celldata[cnew].bb.min[kk] = cell.bb.min[kk];
          for (kk=0; kk<ndim; kk++) celldata[cnew].bb.max[kk] = cell.bb.max[kk];

          // Assign positions of child cells
          if (k == 0 || k == 2 || k == 4 || k == 6) {
            celldata[cnew].rcentre[0] = cell.rcentre[0] - cellSize;
            celldata[cnew].bb.max[0] = cell.rcentre[0];
          }
          else {
            celldata[cnew].rcentre[0] = cell.rcentre[0] + cellSize;
            celldata[cnew].bb.min[0] = cell.rcentre[0];
          }
          if (ndim > 1) {
            if (k == 0 || k == 1 || k == 4 || k == 5) {
              celldata[cnew].rcentre[1] = cell.rcentre[1] - cellSize;
              celldata[cnew].bb.max[1] = cell.rcentre[1];
            }
            else {
              celldata[cnew].rcentre[1] = cell.rcentre[1] + cellSize;
              celldata[cnew].bb.min[1] = cell.rcentre[1];
            }
          }
          if (ndim == 3) {
            if (k < 4) {
              celldata[cnew].rcentre[2] = cell.rcentre[2] - cellSize;
              celldata[cnew].bb.max[2] = cell.rcentre[2];
            }
            else {
              celldata[cnew].rcentre[2] = cell.rcentre[2] + cellSize;
              celldata[cnew].bb.min[2] = cell.rcentre[2];
            }
          }

        }
        //-----------------------------------------------------------------------------------------

        // Divide the particles across the cells
        DivideCell(c, cell.ifirst, cell.ilast, Ncell + cc*Noctchild, partdata, ndim-1);
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
        ExceptionHandler::getIstance().raise("Error : reached maximum oct-tree level");
      }

    }
    //---------------------------------------------------------------------------------------------

    delete[] celllist;

  }
  //-----------------------------------------------------------------------------------------------


  StockTree(celldata[0],partdata,true);
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
  ParticleType<ndim> *partdata,        ///< SPH particle data array
  bool stock_leaf)					   ///< Whether or not to stock leaf cells
{
  int c,cc;                            // Cell counters
  int cend;                            // Last particle in cell
  int i;                               // Particle counter
  int k;                               // Dimension counter
  int l;                               // Level counter
  FLOAT dr[ndim];                      // Relative position vector
  FLOAT drsqd;                         // Distance squared
  FLOAT mi;                            // Mass of particle i
  FLOAT p = (FLOAT) 0.0;               // ??
  FLOAT lambda = (FLOAT) 0.0;          // Quadrupole moment eigenvalue

  debug2("[OctTree::StockTree]");

  const bool need_quadrupole_moments =
      multipole == quadrupole || multipole == fast_quadrupole || gravity_mac == eigenmac ;

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
      for (k=0; k<ndim; k++) cell.r[k]        = (FLOAT) 0.0;
      for (k=0; k<ndim; k++) cell.bb.min[k]   = big_number;
      for (k=0; k<ndim; k++) cell.bb.max[k]   = -big_number;
      for (k=0; k<ndim; k++) cell.hbox.min[k] = big_number;
      for (k=0; k<ndim; k++) cell.hbox.max[k] = -big_number;
      for (k=0; k<5; k++) cell.q[k] = (FLOAT) 0.0;


      // If this is a leaf cell, sum over all particles
      //-------------------------------------------------------------------------------------------
      if (cell.copen == -1 && stock_leaf) {

        // Loop over all particles in cell summing their contributions
        for (i = cell.ifirst; i <= cell.ilast; ++i) {
          if (i == -1) break ;

          if (!partdata[i].flags.is_dead()) {
            cell.N++;
            if (partdata[i].flags.check(active)) cell.Nactive++;
            cell.hmax = max(cell.hmax, partdata[i].h);
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
              if (partdata[i].r[k] + kernrange*partdata[i].h > cell.hbox.max[k])
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
        };

        // Normalise all cell values
        if (cell.N > 0) {
          for (k=0; k<ndim; k++) cell.r[k] /= cell.m;
          for (k=0; k<ndim; k++) dr[k] = (FLOAT) 0.5*(cell.bb.max[k] - cell.bb.min[k]);
          cell.cdistsqd = max(DotProduct(dr,dr,ndim),cell.hmax*cell.hmax)/thetamaxsqd;
          cell.rmax = sqrt(DotProduct(dr,dr,ndim));
        }

        // Compute quadrupole moment terms if selected
        if (need_quadrupole_moments) {
          for (i = cell.ifirst; i <= cell.ilast; ++i) {
            if (i == -1) break ;

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
              else if (ndim == 1) {
                cell.q[0] += mi*((FLOAT) 3.0*dr[0]*dr[0] - drsqd);
              }
            }
          }
        }

      }
      // For non-leaf cells, sum over all child cells
      //-------------------------------------------------------------------------------------------
      else if (cell.copen != -1) {

        // Set limits for children (maximum of 8 but may be less)
        cc   = cell.copen;
        cend = cell.cnext;

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

          cc = child.cnext;
        };

        if (cell.m > 0.0) {
          for (k=0; k<ndim; k++) cell.r[k] /= cell.m;
        }
        //for (k=0; k<ndim; k++) cell.rcell[k] = (FLOAT) 0.5*(cell.bb.min[k] + cell.bb.max[k]);
        for (k=0; k<ndim; k++) dr[k] = (FLOAT) 0.5*(cell.bb.max[k] - cell.bb.min[k]);
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
        };


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
  ParticleType<ndim> *partdata,        ///< SPH particle data array
  bool stock_leaf)                     ///< Whether to stock the leaf
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
      for (k=0; k<ndim; k++) cell.hbox.min[k] = big_number;
      for (k=0; k<ndim; k++) cell.hbox.max[k] = -big_number;

      // If this is a leaf cell, sum over all particles
      //-------------------------------------------------------------------------------------------
      if (cell.copen == -1) {
        if (stock_leaf && cell.ifirst != -1) {
          for (i = cell.ifirst; i <= cell.ilast; ++i) {
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
            for (k=0; k<ndim; k++) cell.hbox.min[k] = min(child.hbox.min[k],cell.hbox.min[k]);
            for (k=0; k<ndim; k++) cell.hbox.max[k] = max(child.hbox.max[k],cell.hbox.max[k]);
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
 (Particle<ndim> *part_gen)            ///< ..
{
  int c;                               // Cell counter
  int i;                               // SPH particle index
  int ilast;                           // Last particle in linked list

  debug2("[OctTree::UpdateActiveParticleCounters]");
  //timing->StartTimingSection("TREE_UPDATE_COUNTERS");
  ParticleType<ndim>* partdata = reinterpret_cast<ParticleType<ndim>*>(part_gen);


  // Loop through all grid cells in turn
  //-----------------------------------------------------------------------------------------------
#pragma omp parallel for default(none) private(c,i,ilast) shared(partdata)
  for (c=0; c<Ncell; c++) {
    celldata[c].Nactive = 0;

    if (celldata[c].copen != -1) continue;
    i = celldata[c].ifirst;
    ilast = celldata[c].ilast;

    // Else walk through linked list to obtain list and number of active ptcls.
    for (i = ifirst; i <= ilast; ++i) {
      if (partdata[i].flags.check(active)) celldata[c].Nactive++;
    };

  }
  //-----------------------------------------------------------------------------------------------

  //timing->EndTimingSection("TREE_UPDATE_COUNTERS");

  return;
}

#ifdef MPI_PARALLEL
//=================================================================================================
//  OctTree::UpdateWorkCounters
/// Calculate the physical properties (e.g. total mass, centre-of-mass,
/// opening-distance, etc..) of all cells in the tree.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void OctTree<ndim,ParticleType,TreeCell>::UpdateWorkCounters
 (TreeCell<ndim> &rootcell)            ///< KD-tree cell
{
  int c,cc;                            // Cell counters
  int cfirst,cend;                     // ..
  int l;                               // ..


  // Loop over all levels in tree starting from lowest
  //===============================================================================================
  for (l=ltot; l>=0; l--) {


    // Loop over all cells on current level
    //---------------------------------------------------------------------------------------------
    for (c=firstCell[l]; c<=lastCell[l]; c++) {
      TreeCell<ndim> &cell = celldata[c];

      // For non-leaf cells, sum over all child cells
      //-------------------------------------------------------------------------------------------
      if (cell.copen != -1) {

        // Set limits for children (maximum of 8 but may be less)
        cfirst = cell.copen;
        cend   = cell.cnext;


        cc = cfirst;
        double worktot = 0 ;
        while (cc != cend) {
          TreeCell<ndim> &child = celldata[cc];

          worktot += child.worktot ;
          cc = child.cnext;
        };
        cell.worktot = worktot ;
      }
      //-------------------------------------------------------------------------------------------

    }
    //---------------------------------------------------------------------------------------------

  }
  //===============================================================================================


  return;
}
#endif


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
      ExceptionHandler::getIstance().raise("Error in cell walk count in OctTree");
    }
    if (celldata[c].level < 0 || celldata[c].level > ltot) {
      cout << "Problem with cell levels : " << celldata[c].level << "    " << ltot << endl;
      ExceptionHandler::getIstance().raise("Error with cell levels in OctTree");
    }
    if (celldata[c].Nactive > 0 && celldata[c].copen != -1) {
      cout << "Problem with active counters : " << celldata[c].Nactive << "    "
           << celldata[c].N << "    " << celldata[c].level << endl;
      ExceptionHandler::getIstance().raise("Error with active counters in OctTree");
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
      if (cell.N == 0) continue ;
      for (i = cell.ifirst; i <= cell.ilast; ++i) {
        pcount[i]++;
        leafcount++;
        Ncount++;
        if (partdata[i].flags.check(active)) activecount++;
        if (partdata[i].flags.check(active)) Nactivecount++;
        if (partdata[i].h > cell.hmax) {
          cout << "hmax flag error : " << c << "    "
               << partdata[i].h << "   " << cell.hmax << endl;
          ExceptionHandler::getIstance().raise("hmax flag error in OctTree");
        }
        for (k=0; k<ndim; k++) {
          if (partdata[i].r[k] < cell.bb.min[k] || partdata[i].r[k] > cell.bb.max[k]) {
            cout << "Bounding box error : " << c << "   " << i << "    " << k << "    r : "
                 << partdata[i].r[k] << "    " << cell.bb.min[k] << "    " << cell.bb.max[k] << endl;
            ExceptionHandler::getIstance().raise("Bounding box error in OctTree");
          }
        }
      }
      if (leafcount > Nleafmax) {
        cout << "Leaf particle count error : " << leafcount << "   " << Nleafmax << endl;
        ExceptionHandler::getIstance().raise("Leaf particle counter error in OctTree");
      }
      if (activecount > leafcount) {
        cout << "Leaf particle count error : " << leafcount << "   " << Nleafmax << endl;
        ExceptionHandler::getIstance().raise("Leaf particle count error in OctTree");
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
    cout << "Ncount problem with Oct-tree : " << Ncount << "   " << Ntot << endl;
    kill_flag = true;
  }

  // Check active particles don't exceed total number of particles
  if (Nactivecount > Ntot) {
    cout << "Nactivecount problem with Oct-tree : " << Nactivecount << "   " << Ntot << endl;
    kill_flag = true;
  }


  delete[] pcount;
  delete[] ccount;

  if (kill_flag) {
    ExceptionHandler::getIstance().raise("kill_flag in OctTree");
  }

  return;
}
#endif





template class OctTree<1, GradhSphParticle, OctTreeCell>;
template class OctTree<2, GradhSphParticle, OctTreeCell>;
template class OctTree<3, GradhSphParticle, OctTreeCell>;

template class OctTree<1, SM2012SphParticle, OctTreeCell>;
template class OctTree<2, SM2012SphParticle, OctTreeCell>;
template class OctTree<3, SM2012SphParticle, OctTreeCell>;

template class OctTree<1, MeshlessFVParticle, OctTreeCell>;
template class OctTree<2, MeshlessFVParticle, OctTreeCell>;
template class OctTree<3, MeshlessFVParticle, OctTreeCell>;

template class OctTree<1, GradhSphParticle, TreeRayCell>;
template class OctTree<2, GradhSphParticle, TreeRayCell>;
template class OctTree<3, GradhSphParticle, TreeRayCell>;

template class OctTree<1, SM2012SphParticle, TreeRayCell>;
template class OctTree<2, SM2012SphParticle, TreeRayCell>;
template class OctTree<3, SM2012SphParticle, TreeRayCell>;

template class OctTree<1, MeshlessFVParticle, TreeRayCell>;
template class OctTree<2, MeshlessFVParticle, TreeRayCell>;
template class OctTree<3, MeshlessFVParticle, TreeRayCell>;
