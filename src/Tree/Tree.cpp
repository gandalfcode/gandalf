//=================================================================================================
//  Tree.cpp
//  Contains all common functions for walking the tree.
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
#include "Hydrodynamics.h"
#include "InlineFuncs.h"
#include "MirrorNeighbours.hpp"
#include "Parameters.h"
#include "Particle.h"
#include "Tree.h"
#include "KDTree.h"
#include "OctTree.h"
#include "Debug.h"
#if defined _OPENMP
#include <omp.h>
#endif
using namespace std;



//=================================================================================================
//  Tree::ComputeActiveParticleList
/// Returns the number (Nactive) and list of ids (activelist) of all active
/// SPH particles in the given cell.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
int Tree<ndim,ParticleType,TreeCell>::ComputeActiveParticleList
 (TreeCellBase<ndim> &cell,                ///< [in] Pointer to cell
  Particle<ndim> *part_gen,        ///< [in] Pointer to particle data array
  int *activelist)                     ///< [out] List of active particles in cell
{
  ParticleType<ndim>* partdata = reinterpret_cast<ParticleType<ndim>* >(part_gen) ;
  const int ilast = cell.ilast;        // i.d. of last particle in cell c
  int i = cell.ifirst;                 // Local particle id (set to first ptcl id)
  int Nactive = 0;                     // No. of active particles in cell

  // Walk through linked list to obtain list and number of active ptcls.
  while (i != -1) {
    if (i < Ntot && partdata[i].active && partdata[i].itype != dead) activelist[Nactive++] = i;
    if (i == ilast) break;
    i = inext[i];
    assert(i < Ntot);
  };

  assert(Nactive <= Nleafmax);
  return Nactive;
}



//=================================================================================================
//  Tree::ComputeActiveCellList
/// Returns the number of cells containing active particles, 'Nactive', and
/// the i.d. list of cells contains active particles, 'celllist'
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
int Tree<ndim,ParticleType,TreeCell>::ComputeActiveCellList
 (TreeCell<ndim> *celllist)            ///< Array containing copies of cells with active ptcls
{
  int c;                               // Cell counter
  int Nactive = 0;                     // No. of active leaf cells in tree

  for (c=0; c<Ncell; c++) {
    if (celldata[c].N <= Nleafmax && celldata[c].copen == -1 && celldata[c].Nactive > 0) {
      celllist[Nactive++] = celldata[c];
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
//  Tree::ComputeActiveCellList
/// Returns the number of cells containing active particles, 'Nactive', and
/// the i.d. list of cells contains active particles, 'celllist'.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
int Tree<ndim,ParticleType,TreeCell>::ComputeActiveCellPointers
 (TreeCellBase<ndim> **celllist)           ///< Array of pointers to cells that will be filled in
{
  int c;                               // Cell counter
  int Nactive = 0;                     // No. of active leaf cells in tree
  assert(celllist != NULL);

  for (c=0; c<Ncell; c++) {
    if (celldata[c].N <= Nleafmax && celldata[c].copen == -1 && celldata[c].Nactive > 0) {
      celllist[Nactive++] = &(celldata[c]);
    }
  }

#ifdef MPI_PARALLEL
  for (c=Ncell; c<Ncell+Nimportedcell; c++) {
    if (celldata[c].Nactive > 0) celllist[Nactive++] = &(celldata[c]);
  }
#endif

  return Nactive;
}


/*
 * Now using template available in InlineFuncs
//=================================================================================================
//  Tree::BoxOverlap
/// Check if two bounding boxes overlap.  If yes, then returns true.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
bool Tree<ndim,ParticleType,TreeCell>::BoxOverlap
 (const FLOAT box1min[ndim],           ///< [in] Minimum extent of box 1
  const FLOAT box1max[ndim],           ///< [in] Maximum extent of box 1
  const FLOAT box2min[ndim],           ///< [in] Minimum extent of box 2
  const FLOAT box2max[ndim])           ///< [in] Maximum extent of box 2
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
*/


//=================================================================================================
//  Tree::ExtrapolateCellProperties
/// Extrapolate important physical properties of all cells in the tree.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void Tree<ndim,ParticleType,TreeCell>::ExtrapolateCellProperties
 (const FLOAT dt)                      ///< [in] Smallest timestep size
{
  int c;                               // Cell counter
  int k;                               // Dimension counter

  debug2("[Tree::ExtrapolateCellProperties]");


  // Loop over all cells and extrapolate all properties
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
//  Tree::ComputeGatherNeighbourList
/// Computes and returns number of neighbour, 'Nneib', and the list of neighbour ids,
/// 'neiblist', for all particles inside cell 'c'.  Includes all particles in the selected cell,
/// plus all particles contained in adjacent cells (including diagonal cells).
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
int Tree<ndim,ParticleType,TreeCell>::ComputeGatherNeighbourList
 (const Particle<ndim> *part_gen,      ///< [in] Particle data array
  const FLOAT rp[ndim],                ///< [in] Search position
  const FLOAT rsearch,                 ///< [in] Maximum smoothing length
  const int Nneibmax,                  ///< [in] Max. no. of neighbours
  int &Nneib,                          ///< [inout] No. of neighbours
  int *neiblist)                       ///< [out] List of neighbour i.d.s
{
  const ParticleType<ndim>* partdata = reinterpret_cast<const ParticleType<ndim>* >(part_gen) ;

  int cc = 0;                          // Cell counter
  int i;                               // Particle id
  int k;                               // Neighbour counter
  FLOAT dr[ndim];                      // Relative position vector
  FLOAT drsqd;                         // Distance squared
  const FLOAT rsearchsqd = rsearch*rsearch;  // Search radius squared
  assert(partdata != NULL);
  assert(neiblist != NULL);


  //===============================================================================================
  while (cc < Ncell) {

    for (k=0; k<ndim; k++) dr[k] = celldata[cc].rcell[k] - rp[k];
    drsqd = DotProduct(dr, dr, ndim);


    // Check if bounding boxes overlap with each other
    //---------------------------------------------------------------------------------------------
    if (drsqd < (rsearch + celldata[cc].rmax)*(rsearch + celldata[cc].rmax)) {

      // If not a leaf-cell, then open cell to first child cell
      if (celldata[cc].copen != -1) {
        cc = celldata[cc].copen;
      }

      // Ensure no empty cells are included in tree walk (due to particle removal)
      else if (celldata[cc].N == 0) {
        cc = celldata[cc].cnext;
      }

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
      else if (celldata[cc].copen == -1 && Nneib + Nleafmax >= Nneibmax) {
        return -1;
      }

    }

    // If not in range, then open next cell
    //---------------------------------------------------------------------------------------------
    else {
      cc = celldata[cc].cnext;
    }

  };
  //===============================================================================================

  assert(Nneib <= Nneibmax);
  return Nneib;
}



//=================================================================================================
//  Tree::ComputeGatherNeighbourList
/// Computes and returns number of neighbour, 'Nneib', and the list of neighbour ids,
/// 'neiblist', for all particles inside cell 'c'.  Includes all particles in the selected cell,
/// plus all particles contained in adjacent cells (including diagonal cells).
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
int Tree<ndim,ParticleType,TreeCell>::ComputeGatherNeighbourList
 (const TreeCellBase<ndim> &cell,      ///< [in] Pointer to current cell
  const Particle<ndim> *part_gen,      ///< [in] Particle data array
  const FLOAT hmax,                    ///< [in] Maximum smoothing length
  const int Nneibmax,                  ///< [in] Max. no. of neighbours
  int &Nneib,                          ///< [inout] No. of neighbours
  int *neiblist)                       ///< [out] List of neighbour i.d.s
{
  const ParticleType<ndim>* partdata = reinterpret_cast<const ParticleType<ndim>* >(part_gen) ;

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
  const FLOAT hrangemaxsqd = pow(cell.rmax + kernrange*hmax,2);
  assert(neiblist != NULL);
  assert(partdata != NULL);

  // Exit immediately if we have overflowed the neighbour list buffer
  if (Nneib == -1) return -1;

  for (k=0; k<ndim; k++) rc[k] = cell.rcell[k];
  for (k=0; k<ndim; k++) gatherboxmin[k] = cell.bbmin[k] - kernrange*hmax;
  for (k=0; k<ndim; k++) gatherboxmax[k] = cell.bbmax[k] + kernrange*hmax;


  //===============================================================================================
  while (cc < Ncell) {

    // Check if bounding boxes overlap with each other
    //---------------------------------------------------------------------------------------------
    if (BoxOverlap(gatherboxmin,gatherboxmax,celldata[cc].bbmin,celldata[cc].bbmax)) {

      // If not a leaf-cell, then open cell to first child cell
      if (celldata[cc].copen != -1) {
        cc = celldata[cc].copen;
      }

      // Ignore any empty cells
      else if (celldata[cc].N == 0) {
        cc = celldata[cc].cnext;
      }

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
      else if (celldata[cc].copen == -1 && Ntemp + Nleafmax >= Nneibmax) {
        return -1;
      }

    }

    // If not in range, then open next cell
    //---------------------------------------------------------------------------------------------
    else {
      cc = celldata[cc].cnext;
    }

  };
  //===============================================================================================


  // Now, trim the list to remove particles that are definitely not neighbours
  assert(Ntemp <= Nneibmax);
  for (j=Nneib; j<Ntemp; j++) {
    i = neiblist[j];
    if (partdata[i].itype == dead) continue;
    for (k=0; k<ndim; k++) dr[k] = partdata[i].r[k] - rc[k];
    drsqd = DotProduct(dr, dr, ndim);
    //cout << "Checking neighbour : " << j << "   " << Nneib << "   " << drsqd << "   " << hrangemaxsqd << endl;
    if (drsqd < hrangemaxsqd) neiblist[Nneib++] = i;
  }

  assert(Nneib <= Nneibmax);
  return Nneib;
}



//=================================================================================================
//  Tree::ComputeNeighbourList
/// Computes and returns number of neighbour, 'Nneib', and the list of neighbour ids,
/// 'neiblist', for all particles inside cell 'c'.  Includes all particles in the selected cell,
/// plus all particles contained in adjacent cells (including diagonal cells).
/// If allocated memory array containing neighbour ids (neiblist) overflows,
/// return with error code (-1) in order to reallocate more memory.
//================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
int Tree<ndim,ParticleType,TreeCell>::ComputeNeighbourList
 (const TreeCellBase<ndim> &cell,      ///< [in] Cell pointer
  const Particle<ndim> *part_gen,      ///< [in] Particle data array
  const int Nneibmax,                  ///< [in] Max. no. of neighbours
  int &Nneib,                          ///< [inout] No. of neighbours
  int *neiblist,                       ///< [out] List of neighbour i.d.s
  Particle<ndim> *neib_out)            ///< [out] Array of local copies of neighbours
{
  const ParticleType<ndim>* partdata = reinterpret_cast<const ParticleType<ndim>* >(part_gen) ;
  ParticleType<ndim>* neibpart = reinterpret_cast<ParticleType<ndim>* >(neib_out) ;

  int cc = 0;                          // Cell counter
  int i;                               // Particle id
  int j;                               // Aux. particle counter
  int k;                               // Neighbour counter
  int Ntemp = Nneib;                   // Temp. no. of neighbouts
  FLOAT dr[ndim];                      // Relative position vector
  FLOAT drsqd;                         // Distance squared
  FLOAT rc[ndim];                      // Position of cell
  const FLOAT hrangemaxsqd = pow(cell.rmax + kernrange*cell.hmax,2);
  const FLOAT rmax = cell.rmax;
  assert(neibpart != NULL);
  assert(partdata != NULL);

  for (k=0; k<ndim; k++) rc[k] = cell.rcell[k];

  // Exit immediately if we have overflowed the neighbour list buffer
  if (Nneib == -1) return -1;


  //===============================================================================================
  while (cc < Ncell) {

    for (k=0; k<ndim; k++) dr[k] = celldata[cc].rcell[k] - rc[k];
    drsqd = DotProduct(dr,dr,ndim);

    // Check if bounding boxes overlap with each other
    //---------------------------------------------------------------------------------------------
    if (BoxOverlap(cell.bbmin, cell.bbmax, celldata[cc].hboxmin, celldata[cc].hboxmax) ||
        BoxOverlap(cell.hboxmin, cell.hboxmax, celldata[cc].bbmin, celldata[cc].bbmax)) {

      // If not a leaf-cell, then open cell to first child cell
      if (celldata[cc].copen != -1) {
        cc = celldata[cc].copen;
      }

      // Ignore empty cells
      else if (celldata[cc].N == 0) {
        cc = celldata[cc].cnext;
      }

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
      else if (celldata[cc].copen == -1 && Ntemp + Nleafmax >= Nneibmax) {
        return -1;
      }

    }

    // If not in range, then open next cell
    //---------------------------------------------------------------------------------------------
    else {
      cc = celldata[cc].cnext;
    }

  };
  //===============================================================================================


  // Now, trim the list to remove particles that are definitely not neighbours
  assert(Ntemp <= Nneibmax);
  for (j=Nneib; j<Ntemp; j++) {
    assert(j < Nneibmax);
    assert(neiblist[j] >= 0);
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

  assert(Nneib <= Nneibmax);
  return Nneib;
}


//=================================================================================================
//  Tree::ComputePeriodicGravityInteractionList
/// Computes and returns number of SPH neighbours (Nneib), including lists of ids, from the
/// tree walk for all active particles inside cell c.  If the interaction list array overflows,
/// returns with error code (-1) to reallocate more memory.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
int Tree<ndim,ParticleType,TreeCell>::ComputePeriodicNeighbourList
 (const TreeCellBase<ndim> &cell,      ///< [in] Pointer to cell
  const Particle<ndim> *part_gen,      ///< [in] Particle data array
  const DomainBox<ndim> &simbox,       ///< [in] Simulation domain box object
  const int Nneibmax,                  ///< [in] Max. no. of SPH neighbours
  int &Nneib,                          ///< [out] Total no. of neighbours
  int *neiblist,                       ///< [out] List of all particle ids
  Particle<ndim> *neib_out)            ///< [out] Array of local copies of neighbour particles
{
  const ParticleType<ndim>* partdata = reinterpret_cast<const ParticleType<ndim>* >(part_gen) ;
  ParticleType<ndim>* neibpart = reinterpret_cast<ParticleType<ndim>* >(neib_out) ;
  int cc = 0;                          // Cell counter
  int i;                               // Particle id
  int j;                               // Aux. particle counter
  int k;                               // Neighbour counter
  int Ntemp = 0;                       // Aux. counter
  FLOAT dr[ndim];                      // Relative position vector
  //FLOAT dr_corr[ndim];                 // Periodic correction vector
  FLOAT drsqd;                         // Distance squared
  FLOAT rc[ndim];                      // Position of cell
  const FLOAT hrangemaxsqd = pow(cell.rmax + kernrange*cell.hmax,2);
  const FLOAT rmax = cell.rmax;
  assert(neibpart != NULL);
  assert(partdata != NULL);

  for (k=0; k<ndim; k++) {
	  rc[k] = cell.rcell[k];
	  //dr_corr[k] = 0 ;
  }

  MirrorNeighbourFinder<ndim> NgbFinder(simbox) ;
  NgbFinder.SetTargetCell(cell) ;

  int MaxMirrors = NgbFinder.MaxNumMirrors() ;

  std::vector<FLOAT> _r_ghost(ndim *pow(3,ndim)) ;
  std::vector<int>   _sign(ndim *pow(3,ndim)) ;

  FLOAT *r_ghost = &(_r_ghost[0]) ;
  int   *sign    = &(_sign[0]) ;

  // Start with root cell and walk through entire tree
  Nneib = 0;
  Ntemp = 0;


  // Walk through all cells in tree to determine particle and cell interaction lists
  //===============================================================================================
  while (cc < Ncell) {

    // Calculate closest periodic replica of cell
    for (k=0; k<ndim; k++) dr[k] = celldata[cc].rcell[k] - rc[k];
    NgbFinder.NearestPeriodicVector(dr);
    drsqd = DotProduct(dr, dr, ndim);

    // Check if bounding spheres overlap with each other (for potential SPH neibs)
    //---------------------------------------------------------------------------------------------
    if (drsqd <= pow(celldata[cc].rmax + cell.rmax + kernrange*cell.hmax,2) ||
        drsqd <= pow(cell.rmax + celldata[cc].rmax + kernrange*celldata[cc].hmax,2)) {

      // If not a leaf-cell, then open cell to first child cell
      if (celldata[cc].copen != -1) {
        cc = celldata[cc].copen;
      }

      // Ignore empty cells
      else if (celldata[cc].N == 0) {
        cc = celldata[cc].cnext;
      }

      // If leaf-cell, add particles to list
      else if (celldata[cc].copen == -1){
    	i = celldata[cc].ifirst;
        while (i != -1) {
          if (Ntemp + MaxMirrors >= Nneibmax){ // Check that we have enough memory
            return -1 ;
          } else {
            int NumMirrors = NgbFinder.ConstructGhostsScatterGather(partdata[i], neibpart + Ntemp) ;

            for (int n(0); n < NumMirrors; ++n){
              for (k=0; k<ndim; k++) dr[k] = neibpart[Ntemp].r[k] - rc[k];

              drsqd = DotProduct(dr, dr, ndim);
              FLOAT h2 = rmax + kernrange*neibpart[Ntemp].h ;
              if (drsqd < hrangemaxsqd || drsqd < h2*h2){
                neiblist[Ntemp] = i ;
                Ntemp++ ;
              } else if (NumMirrors > 1){
            	neibpart[Ntemp] = neibpart[Ntemp+NumMirrors-1] ;
            	NumMirrors-- ;
            	}
              } // Loop over mirrors

            }
          if (i == celldata[cc].ilast) break;
          i = inext[i];
    	  }
        cc = celldata[cc].cnext;
      }

      // If leaf-cell, but we've run out of memory, return with error-code (-1)
      else if (celldata[cc].copen == -1 && Ntemp + Nleafmax >= Nneibmax) {
        return -1;
      }

    }


    // If not in range, then open next cell
    //---------------------------------------------------------------------------------------------
    else {
      cc = celldata[cc].cnext;
    }

  };
  //===============================================================================================

  assert(Ntemp <= Nneibmax);
  // assert(Nneib <= Nneibmax);
  return Ntemp;
}



//=================================================================================================
//  Tree::ComputeGravityInteractionList
/// Computes and returns number of SPH neighbours (Nneib), direct sum particles
/// (Ndirect) and number of cells (Ngravcell), including lists of ids, from
/// the gravity tree walk for active particles inside cell c.
/// Currently defaults to the geometric opening criteria.
/// If any of the interactions list arrays (neiblist,directlist,gravcelllist)
/// overflow, return with error code (-1) to reallocate more memory.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
int Tree<ndim,ParticleType,TreeCell>::ComputeGravityInteractionList
 (const TreeCell<ndim> &cell,          ///< [in] Pointer to cell
  const Particle<ndim> *part_gen,      ///< [in] Particle data array
  const FLOAT macfactor,               ///< [in] Gravity MAC particle factor
  const int Nneibmax,                  ///< [in] Max. no. of SPH neighbours
  const int Ngravcellmax,              ///< [in] Max. no. of cell interactions
  int &Nneib,                          ///< [out] Total no. of neighbours
  int &Nhydroneib,                     ///< [out] No. of SPH neighbours
  int &Ndirect,                        ///< [out] No. of direct-sum neighbours
  int &Ngravcell,                      ///< [out] No. of cell interactions
  int *neiblist,                       ///< [out] List of all particle ids
  int *hydroneiblist,                  ///< [out] List of SPH neibpart ids
  int *directlist,                     ///< [out] List of direct-sum neibpart ids
  TreeCell<ndim> *gravcell,            ///< [out] Array of local copies of tree cells
  Particle<ndim> *neib_out)            ///< [out] Array of local copies of neighbour particles
{
  const ParticleType<ndim>* partdata = reinterpret_cast<const ParticleType<ndim>* >(part_gen) ;
  ParticleType<ndim>* neibpart = reinterpret_cast<ParticleType<ndim>* >(neib_out) ;
  int cc = 0;                          // Cell counter
  int i;                               // Particle id
  int j;                               // Aux. particle counter
  int k;                               // Neighbour counter
  int Nhydroneibtemp = 0;              // Aux. counter
  FLOAT dr[ndim];                      // Relative position vector
  FLOAT drsqd;                         // Distance squared
  FLOAT rc[ndim];                      // Position of cell
  const FLOAT hrangemaxsqd = pow(cell.rmax + kernrange*cell.hmax, 2);
  const FLOAT rmax = cell.rmax;
  assert(gravcell != NULL);
  assert(neibpart != NULL);
  assert(partdata != NULL);

  // Start with root cell and walk through entire tree
  Nneib      = 0;
  Nhydroneib = 0;
  Ndirect    = 0;
  Ngravcell  = 0;


  // Walk through all cells in tree to determine particle and cell interaction lists
  //===============================================================================================
  while (cc < Ncell) {

    for (k=0; k<ndim; k++) dr[k] = celldata[cc].rcell[k] - rc[k];
    drsqd = DotProduct(dr, dr, ndim);


    // Check if bounding boxes overlap with each other
    //---------------------------------------------------------------------------------------------
    if (BoxOverlap(cell.bbmin, cell.bbmax, celldata[cc].hboxmin, celldata[cc].hboxmax) ||
        BoxOverlap(cell.hboxmin, cell.hboxmax, celldata[cc].bbmin, celldata[cc].bbmax)) {

      // If not a leaf-cell, then open cell to first child cell
      if (celldata[cc].copen != -1) {
        cc = celldata[cc].copen;
      }

      // Ignore empty cells
      else if (celldata[cc].N == 0) {
        cc = celldata[cc].cnext;
      }

      // If leaf-cell, add particles to list
      else if (celldata[cc].copen == -1 && Nneib + Nleafmax <= Nneibmax) {
        i = celldata[cc].ifirst;
        while (i != -1) {
          hydroneiblist[Nhydroneib++] = Nneib;
          neiblist[Nneib] = i;
          neibpart[Nneib++] = partdata[i];
          if (i == celldata[cc].ilast) break;
          i = inext[i];
        };
        cc = celldata[cc].cnext;
      }

      // If leaf-cell, but we've run out of memory, return with error-code (-1)
      else if (celldata[cc].copen == -1 && Nneib + Nleafmax > Nneibmax) {
        return -1;
      }

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
      if (celldata[cc].copen != -1) {
         cc = celldata[cc].copen;
      }

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
      else if (celldata[cc].copen == -1 && Nneib + Nleafmax > Nneibmax) {
        return -3;
      }

    }

    // If not in range, then open next cell
    //---------------------------------------------------------------------------------------------
    else {
      cc = celldata[cc].cnext;
    }

  };
  //===============================================================================================


  // Now, trim the list to remove particles that are definitely not SPH neighbours.
  // If not an SPH neighbour, then add to direct gravity sum list.
  for (j=Nhydroneibtemp; j<Nhydroneib; j++) {
    assert(j < Nneibmax);
    i = hydroneiblist[j];
    if (neibpart[i].itype == dead) continue;
    for (k=0; k<ndim; k++) dr[k] = neibpart[i].r[k] - rc[k];
    drsqd = DotProduct(dr,dr,ndim);
    if (drsqd < hrangemaxsqd || drsqd <
        (rmax + kernrange*neibpart[i].h)*(rmax + kernrange*neibpart[i].h)) {
      hydroneiblist[Nhydroneibtemp++] = i;
    }
    else if (Ndirect < Nneibmax) {
      directlist[Ndirect++] = i;
    }
    else {
     return -3;
    }
  }
  Nhydroneib = Nhydroneibtemp;

  assert(Nhydroneib <= Nneibmax);
  assert(Nneib <= Nneibmax);
  assert(Ndirect <= Nneibmax);
  assert(Ngravcell <= Ngravcellmax);
  assert(VerifyUniqueIds(Nneib, Ntot, neiblist));

  return 1;
}



//=================================================================================================
//  Tree::ComputePeriodicGravityInteractionList
/// Computes and returns number of SPH neighbours (Nneib), direct sum particles (Ndirect) and
/// number of cells (Ngravcell), including lists of ids, from the gravity tree walk for active
/// particles inside cell c.  Currently defaults to the geometric opening criteria.
/// If any of the interactions list arrays (neiblist,directlist,gravcelllist)
/// overflow, return with error code (-1) to reallocate more memory.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
int Tree<ndim,ParticleType,TreeCell>::ComputePeriodicGravityInteractionList
 (const TreeCell<ndim> &cell,          ///< [in] Pointer to cell
  const Particle<ndim> *part_gen,      ///< [in] Particle data array
  const DomainBox<ndim> &simbox,       ///< [in] Simulation domain box object
  const FLOAT macfactor,               ///< [in] Gravity MAC particle factor
  const int Nneibmax,                  ///< [in] Max. no. of SPH neighbours
  const int Ngravcellmax,              ///< [in] Max. no. of cell interactions
  int &Nneib,                          ///< [out] Total no. of neighbours
  int &Nhydroneib,                     ///< [out] No. of SPH neighbours
  int &Ndirect,                        ///< [out] No. of direct-sum neighbours
  int &Ngravcell,                      ///< [out] No. of cell interactions
  int *neiblist,                       ///< [out] List of all particle ids
  int *hydroneiblist,                  ///< [out] List of SPH neibpart ids
  int *directlist,                     ///< [out] List of direct-sum neibpart ids
  TreeCell<ndim> *gravcell,            ///< [out] Array of local copies of tree cells
  Particle<ndim> *neib_out)            ///< [out] Array of local copies of neighbour particles
{
  const ParticleType<ndim>* partdata = reinterpret_cast<const ParticleType<ndim>* >(part_gen) ;
  ParticleType<ndim>* neibpart = reinterpret_cast<ParticleType<ndim>* >(neib_out) ;
  int cc = 0;                          // Cell counter
  int i;                               // Particle id
  int j;                               // Aux. particle counter
  int k;                               // Neighbour counter
  int Nhydroneibtemp = 0;              // Aux. counter
  FLOAT dr[ndim];                      // Relative position vector
  FLOAT dr_corr[ndim];                 // Periodic correction vector
  FLOAT drsqd;                         // Distance squared
  FLOAT rc[ndim];                      // Position of cell
  assert(directlist != NULL);
  assert(gravcell != NULL);
  assert(hydroneiblist != NULL);
  assert(neibpart != NULL);
  assert(partdata != NULL);

  MirrorNeighbourFinder<ndim> NgbFinder(simbox) ;
  NgbFinder.SetTargetCell(cell) ;

  assert(NgbFinder.MaxNumMirrors() == 1) ;

  FLOAT r_ghost[ndim] ;
  int   sign[ndim] ;

  for (k=0; k<ndim; k++) rc[k] = cell.rcell[k];


  // Make local copies of important cell properties
  const FLOAT hrangemaxsqd = pow(cell.rmax + kernrange*cell.hmax,2);
  const FLOAT rmax = cell.rmax;
  for (k=0; k<ndim; k++) rc[k] = cell.rcell[k];

  // Start with root cell and walk through entire tree
  Nneib      = 0;
  Nhydroneib = 0;
  Ndirect    = 0;
  Ngravcell  = 0;


  // Walk through all cells in tree to determine particle and cell interaction lists
  //===============================================================================================
  while (cc < Ncell) {

    // Calculate closest periodic replica of cell
    for (k=0; k<ndim; k++) dr[k] = celldata[cc].rcell[k] - rc[k];
    NgbFinder.NearestPeriodicVector(dr);
    drsqd = DotProduct(dr, dr, ndim);

    // Check if bounding spheres overlap with each other (for potential SPH neibs)
    //---------------------------------------------------------------------------------------------
    if (drsqd <= pow(celldata[cc].rmax + cell.rmax + kernrange*cell.hmax,2) ||
        drsqd <= pow(cell.rmax + celldata[cc].rmax + kernrange*celldata[cc].hmax,2)) {

      // If not a leaf-cell, then open cell to first child cell
      if (celldata[cc].copen != -1) {
        cc = celldata[cc].copen;
      }

      // Ignore empty cells
      else if (celldata[cc].N == 0) {
        cc = celldata[cc].cnext;
      }

      // If leaf-cell, add particles to list
      else if (celldata[cc].copen == -1 && Nneib + Nleafmax <= Nneibmax) {

        i = celldata[cc].ifirst;
        while (i != -1) {
          hydroneiblist[Nhydroneib++] = Nneib;
          neiblist[Nneib] = i;

          NgbFinder.ConstructGhostsScatterGather(partdata[i], neibpart + Nneib) ;

          Nneib++;
          if (i == celldata[cc].ilast) break;
          i = inext[i];
        };
        cc = celldata[cc].cnext;
      }

      // If leaf-cell, but we've run out of memory, return with error-code (-1)
      else if (celldata[cc].copen == -1 && Nneib + Nleafmax > Nneibmax) {
        return -1;
      }

    }

    // Check if cell is far enough away to use the COM approximation
    //---------------------------------------------------------------------------------------------
    else if (drsqd > celldata[cc].cdistsqd && drsqd > celldata[cc].mac*macfactor &&
             celldata[cc].N > 0) {

      // If cell is a leaf-cell with only one particle, more efficient to
      // compute the gravitational contribution from the particle than the cell
      if (celldata[cc].copen == -1 && celldata[cc].N == 1 && Ndirect + Nleafmax < Nneibmax) {
        i = celldata[cc].ifirst;
        directlist[Ndirect++] = Nneib;
        neiblist[Nneib] = i;
        neibpart[Nneib] = partdata[i];
        for (k=0; k<ndim; k++) dr[k] = neibpart[Nneib].r[k] - rc[k];
        NgbFinder.NearestPeriodicVector(dr);
        for (k=0; k<ndim; k++) neibpart[Nneib].r[k] = rc[k] + dr[k] ;
        Nneib++;
      }
      else if (Ngravcell < Ngravcellmax) {
        gravcell[Ngravcell] = celldata[cc];
        for (k=0; k<ndim; k++) dr[k] = celldata[cc].rcell[k] - rc[k] ;
        NgbFinder.PeriodicDistanceCorrection(dr, dr_corr);
        for (k=0; k<ndim; k++) gravcell[Ngravcell].r[k] += dr_corr[k] ;
        for (k=0; k<ndim; k++) gravcell[Ngravcell].rcell[k] += dr_corr[k];
        Ngravcell++;
      }
      else {
        return -2;
      }
      cc = celldata[cc].cnext;

    }

    // If cell is too close, open cell to interogate children cells.
    // If cell is too close and a leaf cell, then add particles to direct list.
    //---------------------------------------------------------------------------------------------
    else if (!(drsqd > celldata[cc].cdistsqd && drsqd > celldata[cc].mac*macfactor)
             && celldata[cc].N > 0) {

      // If not a leaf-cell, then open cell to first child cell
      if (celldata[cc].copen != -1) {
        cc = celldata[cc].copen;
      }

      // If leaf-cell, add particles to list
      else if (celldata[cc].copen == -1 && Nneib + Nleafmax <= Nneibmax) {
        i = celldata[cc].ifirst;
        while (i != -1) {
          directlist[Ndirect++] = Nneib;
          neiblist[Nneib] = i;
          neibpart[Nneib] = partdata[i];
          for (k=0; k<ndim; k++) dr[k] = neibpart[Nneib].r[k] - rc[k];
          NgbFinder.NearestPeriodicVector(dr);
          for (k=0; k<ndim; k++) neibpart[Nneib].r[k] = rc[k] + dr[k] ;

          Nneib++;
          if (i == celldata[cc].ilast) break;
          i = inext[i];
        };
        cc = celldata[cc].cnext;
      }

      // If leaf-cell, but we've run out of memory, return with error-code (-1)
      else if (celldata[cc].copen == -1 && Nneib + Nleafmax > Nneibmax) {
        return -3;
      }

    }

    // If not in range, then open next cell
    //---------------------------------------------------------------------------------------------
    else {
      cc = celldata[cc].cnext;
    }

  };
  //===============================================================================================


  // Now, trim the list to remove particles that are definitely not SPH neighbours.
  // If not an SPH neighbour, then add to direct gravity sum list.
  for (j=Nhydroneibtemp; j<Nhydroneib; j++) {
    i = hydroneiblist[j];
    if (neibpart[i].itype == dead) continue;
    for (k=0; k<ndim; k++) dr[k] = neibpart[i].r[k] - rc[k];
    drsqd = DotProduct(dr, dr, ndim);
    if (drsqd < hrangemaxsqd ||
        drsqd < (rmax + kernrange*neibpart[i].h)*(rmax + kernrange*neibpart[i].h)) {
      hydroneiblist[Nhydroneibtemp++] = i;
    }
    else if (Ndirect + Nhydroneibtemp < Nneibmax) {
      directlist[Ndirect++] = i;
    }
    else {
      return -3;
    }
  }
  Nhydroneib = Nhydroneibtemp;

  assert(Nhydroneib <= Nneibmax);
  assert(Nneib <= Nneibmax);
  assert(Ndirect <= Nneibmax);
  assert(Ngravcell <= Ngravcellmax);
  assert(VerifyUniqueIds(Nneib, Ntot, neiblist));

  return 1;
}



//=================================================================================================
//  Tree::ComputeStarGravityInteractionList
/// Computes and returns number of SPH neighbours (Nneib), direct sum particles (Ndirect) and
/// number of cells (Ngravcell), including lists of ids, from the gravity tree walk for active
/// particles inside cell c.  Currently defaults to the geometric opening criteria.
/// If any of the interactions list arrays (neiblist,directlist,gravcelllist)
/// overflow, return with error code (-1) to reallocate more memory.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
int Tree<ndim,ParticleType,TreeCell>::ComputeStarGravityInteractionList
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
  Particle<ndim> *part_gen)            ///< [in] Particle data array
{
  ParticleType<ndim>* partdata = reinterpret_cast<ParticleType<ndim>* >(part_gen) ;
  int cc = 0;                          // Cell counter
  int i;                               // Particle id
  int k;                               // Neighbour counter
  FLOAT dr[ndim];                      // Relative position vector
  FLOAT drsqd;                         // Distance squared
  FLOAT hrangemax;                     // Maximum kernel extent
  FLOAT rs[ndim];                      // Position of star
  assert(directlist != NULL);
  assert(gravcell != NULL);
  assert(neiblist != NULL);
  assert(partdata != NULL);


  // Make local copies of important cell properties
  for (k=0; k<ndim; k++) rs[k] = star->r[k];
  hrangemax = kernrange*star->h;

  Nneib = 0;
  Ndirect = 0;
  Ngravcell = 0;


  // Walk through all cells in tree to determine particle and cell interaction lists
  //===============================================================================================
  while (cc < Ncell) {

    for (k=0; k<ndim; k++) dr[k] = celldata[cc].rcell[k] - rs[k];
    drsqd = DotProduct(dr, dr, ndim);


    // Check if cells contain neighbours
    //---------------------------------------------------------------------------------------------
    if (drsqd < pow((FLOAT) 0.5*hrangemax +
        celldata[cc].rmax + (FLOAT) 0.5*kernrange*celldata[cc].hmax,2)) {

      // If not a leaf-cell, then open cell to first child cell
      if (celldata[cc].copen != -1) {
        cc = celldata[cc].copen;
      }

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
      else if (celldata[cc].copen == -1 && Nneib + Nleafmax > Nneibmax) {
        return -1;
      }

    }

    // Check if cell is far enough away to use the COM approximation
    //---------------------------------------------------------------------------------------------
    else if (drsqd > celldata[cc].cdistsqd && celldata[cc].N > 0) {

      // If cell is a leaf-cell with only one particle, more efficient to
      // compute the gravitational contribution from the particle than the cell
      if (celldata[cc].copen == -1 && celldata[cc].N == 1 && Ndirect + Nneib < Ndirectmax) {
        if (partdata[celldata[cc].ifirst].itype != dead) {
          directlist[Ndirect++] = celldata[cc].ifirst;
        }
      }
      else if (Ngravcell < Ngravcellmax && celldata[cc].N > 0) {
        gravcell[Ngravcell++] = celldata[cc];
      }
      else {
        return -1;
      }
      cc = celldata[cc].cnext;

    }

    // If cell is too close, open cell to interogate children cells.
    // If cell is too close and a leaf cell, then add particles to direct list.
    //---------------------------------------------------------------------------------------------
    else if (drsqd <= celldata[cc].cdistsqd && celldata[cc].N > 0) {

      // If not a leaf-cell, then open cell to first child cell
      if (celldata[cc].copen != -1) {
        cc = celldata[cc].copen;
      }

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
      else if (celldata[cc].copen == -1 && Ndirect + Nleafmax > Ndirectmax) {
        return -1;
      }

    }

    // If not in range, then open next cell
    //---------------------------------------------------------------------------------------------
    else {
      cc = celldata[cc].cnext;
    }

  };
  //===============================================================================================

  assert(Nneib <= Nneibmax);
  assert(Ndirect <= Nneibmax);
  assert(Ngravcell <= Ngravcellmax);
  assert(VerifyUniqueIds(Nneib, Ntot, neiblist));

  return 1;
}



//=================================================================================================
//  Tree::FindLeafCell
/// Finds and returns the leaf cell i.d. containing the given point, rp.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
int Tree<ndim,ParticleType,TreeCell>::FindLeafCell
 (const FLOAT rp[ndim])
{
  int cc = 0;                          // Cell counter (start at root cell)


  //===============================================================================================
  while (cc < Ncell) {

    // Check if bounding boxes overlap with each other
    //---------------------------------------------------------------------------------------------
    if (BoxOverlap(rp, rp, celldata[cc].bbmin, celldata[cc].bbmax)) {

      // If not a leaf-cell, then open cell to first child cell
      if (celldata[cc].copen != -1) {
        cc = celldata[cc].copen;
      }

      // Ignore any empty cells
      else if (celldata[cc].N == 0) {
        cc = celldata[cc].cnext;
      }

      // If leaf-cell, add particles to list
      else if (celldata[cc].copen == -1) {
        return cc;
      }

    }

    // If not in range, then open next cell
    //---------------------------------------------------------------------------------------------
    else {
      cc = celldata[cc].cnext;
    }

  };
  //===============================================================================================

  // If we've somehow not found a leaf cell (e.g. if particle is outside tree domain), then
  // return error code, -1
  return -1;
}



#ifdef MPI_PARALLEL
//=================================================================================================
//  Tree::CreatePrunedTreeForMpiNode
/// Walks through the main local tree and constructs a pruned tree for the given MPI node.
/// If the MPI node is actually the local node, then simply build the pruned tree to the minimum
/// pruning depth.  Otherwise, creates a pruned tree designed to allow gravity to be walked on the
/// given node without requiring transfer of particles (unless very close to the node boundary).
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
int Tree<ndim,ParticleType,TreeCell>::CreatePrunedTreeForMpiNode
 (const MpiNode<ndim> &mpinode,        ///< [in] MPI node to create pruned tree for
  const DomainBox<ndim> &simbox,       ///< [in] Simulation domain box object
  const FLOAT macfactor,               ///< [in] Gravity MAC particle factor
  const bool localNode,                ///< [in] Flag if MPI node is actually the local node
  const int pruning_level_min,         ///< [in] Minimum pruning level
  const int pruning_level_max,         ///< [in] Maximum pruning level
  const int Nprunedcellmax,            ///< [in] Max. no. of cell interactions
  TreeCell<ndim> *prunedcells)         ///< [out] List of cell pointers in pruned tree
{
  int c;                               // Cell counter
  int cnext;                           // id of next cell in tree
  int i;                               // Particle id
  int j;                               // Aux. particle counter
  int k;                               // Neighbour counter
  int Nprunedcell = 0;                 // No. of cells in newly created pruned tree
  int *newCellIds;                     // New cell ids in pruned tree
  FLOAT dr[ndim];                      // Relative position vector
  FLOAT dr_corr[ndim];                 // Periodic correction vector
  FLOAT drsqd;                         // Distance squared
  FLOAT rnode[ndim];                   // Position of cell
  FLOAT rmin[ndim];                    // Minimum extent of MPI node
  FLOAT rmax[ndim];                    // Maximum extent of MPI node
  FLOAT rsize[ndim];                   // Size of MPI node

  newCellIds = new int[Ncellmax + 1];
  for (c=0; c<Ncellmax+1; c++) newCellIds[c] = -1;
  //newCellIds[Ncell] = Ncell;
  //newCellIds[Ncellmax] = Ncellmax;
  assert(Nprunedcellmax <= Ncellmax);

  for (k=0; k<ndim; k++) rmin[k] = mpinode.hbox.boxmin[k];
  for (k=0; k<ndim; k++) rmax[k] = mpinode.hbox.boxmax[k];
  for (k=0; k<ndim; k++) rnode[k] = (FLOAT) 0.5*(mpinode.hbox.boxmin[k] + mpinode.hbox.boxmax[k]);
  for (k=0; k<ndim; k++) rsize[k] = mpinode.hbox.boxmax[k] - rnode[k];

  /*cout << "NODE      : " << rmin[0] << "    " << rnode[0] << "    " << rmax[0] << "    size : " << rsize[0] << endl;
  cout << "FULL TREE : " << celldata[0].bbmin[0] << "    " << celldata[0].bbmax[0] << "    "
       << celldata[0].rmax << "   " << celldata[0].hmax << endl;
  cout << "Main tree? : " << c << "     Ncell : " << Ncell << "   " << Ncellmax << endl;*/

  c = 0;

  // Walk through all cells in tree to determine particle and cell interaction lists
  //===============================================================================================
  while (c < Ncell) {

    // Return with error message if we've run out of memory for the pruned tree
    if (Nprunedcell == Nprunedcellmax) {
      delete[] newCellIds;
      return -1;
    }

    // Add cell to pruned tree
    newCellIds[c] = Nprunedcell;
    prunedcells[Nprunedcell] = celldata[c];
    assert(Nprunedcell <= Nprunedcellmax);


    // Calculate closest periodic replica of cell
    for (k=0; k<ndim; k++) dr[k] = celldata[c].rcell[k] - rnode[k];
    NearestPeriodicVector(simbox, dr, dr_corr);

    // Find vector to nearest edge of the node box
    for (k=0; k<ndim; k++) {
      if (dr[k] > rsize[k]) dr[k] = dr[k] - rsize[k];
      else if (dr[k] < -rsize[k]) dr[k] = dr[k] + rsize[k];
      else dr[k] = 0.0;
    }
    drsqd = DotProduct(dr, dr, ndim);


    // Special case if creating pruned tree for local node to prevent pruned tree from being
    // deeper than the minimum pruning level
    //---------------------------------------------------------------------------------------------
    if (localNode && celldata[c].level == pruning_level_min) {
      prunedcells[Nprunedcell].copen = -1;
      c = celldata[c].cnext;
    }

    // If tree has reached maximum pruning level, or cell is a leaf cell, then record cell
    // and then move to next cell on same or lower level
    //---------------------------------------------------------------------------------------------
    else if (celldata[c].level == pruning_level_max || celldata[c].copen == -1) {
      prunedcells[Nprunedcell].copen = -1;
      c = celldata[c].cnext;
    }

    // If tree has not reached minimum pruning level, then enforce opening to lower level
    //---------------------------------------------------------------------------------------------
    else if (celldata[c].level < pruning_level_min && celldata[c].copen != -1) {
      c = celldata[c].copen;
    }

    // Check if bounding spheres overlap with each other (for potential SPH neibs)
    //---------------------------------------------------------------------------------------------
    else if (drsqd <= pow(celldata[c].rmax + kernrange*celldata[c].hmax,2) ||
             (!(drsqd > celldata[c].cdistsqd && drsqd > celldata[c].mac*macfactor) &&
              celldata[c].N > 0)) {
      c = celldata[c].copen;
    }

    // Otherwise then open next cell
    //---------------------------------------------------------------------------------------------
    else {
      prunedcells[Nprunedcell].copen = -1;
      c = celldata[c].cnext;
    }

    Nprunedcell++;

  };
  //===============================================================================================


  newCellIds[Ncell] = Nprunedcell;
  newCellIds[Ncellmax] = Nprunedcell;


  // Change all cell pointers to new ids in pruned tree
  for (c=0; c<Nprunedcell; c++) {
    if (prunedcells[c].copen != -1) prunedcells[c].copen = newCellIds[prunedcells[c].copen];
    prunedcells[c].cnext = newCellIds[prunedcells[c].cnext];
    if (prunedcells[c].cnext <= 0 || prunedcells[c].copen >= prunedcells[c].cnext) {
      cout << "Problem with new pointers : " << c << "    " << Nprunedcell << "   "
           << "    " << Nprunedcellmax << "    copen : " << prunedcells[c].copen
           << "    cnext : " << prunedcells[c].cnext << endl;
      exit(0);
    }
    assert(prunedcells[c].cnext > 0);
    assert(prunedcells[c].copen < prunedcells[c].cnext);
  }

  // If selected, verify that pruned tree pointers are correctly set-up
  assert(Nprunedcell <= Nprunedcellmax);
  for (c=0; c<Nprunedcell; c++) {
    if (prunedcells[c].copen != -1) cnext = prunedcells[c].copen;
    else cnext = prunedcells[c].cnext;
    assert(cnext == c + 1 || cnext == Ncell || cnext == Ncellmax);
    assert(prunedcells[c].level <= pruning_level_max);
    assert(prunedcells[c].cnext > 0);
    assert(prunedcells[c].copen < prunedcells[c].cnext);
  }


  delete[] newCellIds;

  return Nprunedcell;
}



//=================================================================================================
//  Tree::ComputeDistantGravityInteractionList
/// Compute the list of cells to compute distant gravitational forces for particles in the given
/// cell.  If the cell is too close, then return -1 to signal that the cell must be exported.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
int Tree<ndim,ParticleType,TreeCell>::ComputeDistantGravityInteractionList
 (const TreeCell<ndim> *cellptr,       ///< [in] Pointer to cell
  const DomainBox<ndim> &simbox,       ///< [in] Simulation domain box object
  const FLOAT macfactor,               ///< [in] Gravity MAC particle factor
  const int Ngravcellmax,              ///< [in] Max. no. of cell interactions
  int Ngravcell,                       ///< [in] Current no. of cells in array
  TreeCell<ndim> *gravcelllist)        ///< [out] Array of cells
{
  int cc = 0;                          // Cell counter
  int k;                               // Dimension counter
  int Ngravcelltemp = Ngravcell;       // Aux. cell counter
  FLOAT dr[ndim];                      // Relative position vector
  FLOAT dr_corr[ndim];                 // Correction displacement for periodic boundaries
  FLOAT drsqd;                         // Distance squared
  FLOAT hrangemax;                     // Maximum kernel extent
  FLOAT rc[ndim];                      // Position of cell
  FLOAT rmax;                          // Radius of sphere containing particles

  // Make local copies of important cell properties
  for (k=0; k<ndim; k++) rc[k] = cellptr->rcell[k];
  hrangemax = cellptr->rmax + kernrange*cellptr->hmax;
  rmax = cellptr->rmax;


  // Walk through all cells in tree to determine particle and cell interaction lists
  //===============================================================================================
  while (cc < Ncell) {

    // Calculate closest periodic replica of cell
    for (k=0; k<ndim; k++) dr[k] = celldata[cc].rcell[k] - rc[k];
    NearestPeriodicVector(simbox, dr, dr_corr);
    drsqd = DotProduct(dr, dr, ndim);

    // Check if bounding spheres overlap with each other (for potential SPH neibs)
    //---------------------------------------------------------------------------------------------
    if (drsqd <= pow(celldata[cc].rmax + cellptr->rmax + kernrange*cellptr->hmax,2) ||
        drsqd <= pow(cellptr->rmax + celldata[cc].rmax + kernrange*celldata[cc].hmax,2)) {

      // If not a leaf-cell, then open cell to first child cell
      if (celldata[cc].copen != -1) {
        cc = celldata[cc].copen;
      }

      // Ignore empty cells
      else if (celldata[cc].N == 0) {
        cc = celldata[cc].cnext;
      }

      // If leaf-cell, then pruned tree is invalid so return -1 error code to export cell
      else if (celldata[cc].copen == -1 || celldata[cc].N <= Nleafmax) {
        return -1;
      }

    }

    // Check if cell is far enough away to use the COM approximation.
    // If so, then make copy of cell and move to next cell.
    //---------------------------------------------------------------------------------------------
    else if (drsqd > celldata[cc].cdistsqd && drsqd > celldata[cc].mac*macfactor &&
             celldata[cc].N > 0) {

      gravcelllist[Ngravcelltemp++] = celldata[cc];
      if (Ngravcelltemp >= Ngravcellmax) {
        ExceptionHandler::getIstance().raise("Too many interaction cells "
                                             "in distant gravity interaction list!");
      }
      cc = celldata[cc].cnext;

    }

    // If cell is too close, open cell to interogate children cells.
    // If cell is too close and a leaf cell, then add particles to direct list.
    //---------------------------------------------------------------------------------------------
    else if (!(drsqd > celldata[cc].cdistsqd && drsqd > celldata[cc].mac*macfactor) &&
             celldata[cc].N > 0) {

      // If not a leaf-cell, then open cell to first child cell
      if (celldata[cc].copen != -1) {
        cc = celldata[cc].copen;
      }

      // If leaf-cell, add particles to list
      else {
        return -1;
      }

    }

    // If not in range, then open next cell
    //---------------------------------------------------------------------------------------------
    else {
      cc = celldata[cc].cnext;
    }

    assert(cc != -1);
    assert(cc <= Ncell);

  };
  //===============================================================================================

  return Ngravcelltemp;
}



//=================================================================================================
//  Tree::ComputeHydroTreeCellOverlap
/// Computes if a given cell overlaps with any of the tree's lowest level cells (either leaf
/// cells or lowest level cells in the case of a pruned tree).  Used to decide if a cell should
/// be exported to another MPI node for local hydro forces computations.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
bool Tree<ndim,ParticleType,TreeCell>::ComputeHydroTreeCellOverlap
 (const TreeCell<ndim> *cellptr,       ///< [in] Pointer to cell
  const DomainBox<ndim> &simbox)       ///< [in] Simulation domain box object
{
  int cc = 0;                          // Cell counter
  int k;                               // Neighbour counter
  FLOAT dr[ndim];                      // Relative position vector
  FLOAT dr_corr[ndim];                 // Periodic correction vector
  FLOAT drsqd;                         // Distance squared
  FLOAT rc[ndim];                      // Position of cell

  // Make local copies of important cell properties
  for (k=0; k<ndim; k++) rc[k] = cellptr->rcell[k];


  // Walk through all cells in tree to determine particle and cell interaction lists
  //===============================================================================================
  while (cc < Ncell) {

    // Calculate closest periodic replica of cell
    for (k=0; k<ndim; k++) dr[k] = celldata[cc].rcell[k] - rc[k];
    NearestPeriodicVector(simbox, dr, dr_corr);
    drsqd = DotProduct(dr, dr, ndim);

    // Check if bounding spheres overlap with each other (for potential SPH neibs)
    //---------------------------------------------------------------------------------------------
    if (drsqd <= pow(celldata[cc].rmax + cellptr->rmax + kernrange*cellptr->hmax,2) ||
        drsqd <= pow(cellptr->rmax + celldata[cc].rmax + kernrange*celldata[cc].hmax,2)) {

      // If not a leaf-cell, then open cell to first child cell
      if (celldata[cc].copen != -1) {
        cc = celldata[cc].copen;
      }

      // If cell contains no particle (dead leaf?) then move to next cell
      else if (celldata[cc].N == 0) {
        cc = celldata[cc].cnext;
      }

      // If cell is overlapping with any leaf, then flag overlap on return
      else if (celldata[cc].copen == -1 && celldata[cc].N > 0) {
        return true;
      }

    }

    // If not in range, then open next cell
    //---------------------------------------------------------------------------------------------
    else {
      cc = celldata[cc].cnext;
    }

  };
  //===============================================================================================


  // If we've walked the entire tree wihout any leaf overlaps, return flag for no overlap
  return false;
}



//=================================================================================================
//  Tree::ComputeWorkInBox
/// Compute the total CPU work contained in the given box (representing the potential extent of
/// a new MPI node after load balancing) that overlaps with the tree.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
FLOAT Tree<ndim,ParticleType,TreeCell>::ComputeWorkInBox
 (const FLOAT boxmin[ndim],            ///< [in] Minimum extent of box
  const FLOAT boxmax[ndim])            ///< [in] Maximum extent of box
{
  int c = 0;                           // Cell counter
  FLOAT fracoverlap;                   // Fractional overlap of cell with given box
  FLOAT worktot = (FLOAT) 0.0;         // Total work accumulator


  // Walk through all cells in tree to determine particle and cell interaction lists
  //===============================================================================================
  while (c < Ncell) {

    // Compute what fraction of the cell h-box overlaps with the given MPI node box
    fracoverlap = FractionalBoxOverlap
      (ndim, boxmin, boxmax, celldata[c].hboxmin, celldata[c].hboxmax);

    // If there is zero or full overlap, record the value and move to the next cell
    if (fracoverlap < small_number || fracoverlap > (FLOAT) 0.999999999999999999) {
      worktot += fracoverlap*celldata[c].worktot;
      c = celldata[c].cnext;
    }

    // If there is a parital overlap for a non-leaf cell, then open cell to lower levels
    else if (celldata[c].copen != -1) {
      c = celldata[c].copen;
    }

    // If there is a partial overlap for a leaf-cell, record overlap value
    else if (fracoverlap >= small_number && fracoverlap <= (FLOAT) 1.0 - small_number) {
      worktot += fracoverlap*celldata[c].worktot;
      c = celldata[c].cnext;
    }

    // Code should not technically reach here (unless there's a problem)
    else {
      cout << "Problem with overlap of pruned trees" << endl;
      exit(0);
    }

  };
  //===============================================================================================

  return worktot;
}
#endif



template class Tree<1,Particle,KDTreeCell>;
template class Tree<2,Particle,KDTreeCell>;
template class Tree<3,Particle,KDTreeCell>;
template class Tree<1,SphParticle,KDTreeCell>;
template class Tree<2,SphParticle,KDTreeCell>;
template class Tree<3,SphParticle,KDTreeCell>;
template class Tree<1,GradhSphParticle,KDTreeCell>;
template class Tree<2,GradhSphParticle,KDTreeCell>;
template class Tree<3,GradhSphParticle,KDTreeCell>;
template class Tree<1,SM2012SphParticle,KDTreeCell>;
template class Tree<2,SM2012SphParticle,KDTreeCell>;
template class Tree<3,SM2012SphParticle,KDTreeCell>;
template class Tree<1,MeshlessFVParticle,KDTreeCell>;
template class Tree<2,MeshlessFVParticle,KDTreeCell>;
template class Tree<3,MeshlessFVParticle,KDTreeCell>;


template class Tree<1,Particle,OctTreeCell>;
template class Tree<2,Particle,OctTreeCell>;
template class Tree<3,Particle,OctTreeCell>;
template class Tree<1,SphParticle,OctTreeCell>;
template class Tree<2,SphParticle,OctTreeCell>;
template class Tree<3,SphParticle,OctTreeCell>;
template class Tree<1,GradhSphParticle,OctTreeCell>;
template class Tree<2,GradhSphParticle,OctTreeCell>;
template class Tree<3,GradhSphParticle,OctTreeCell>;
template class Tree<1,SM2012SphParticle,OctTreeCell>;
template class Tree<2,SM2012SphParticle,OctTreeCell>;
template class Tree<3,SM2012SphParticle,OctTreeCell>;
template class Tree<1,MeshlessFVParticle,OctTreeCell>;
template class Tree<2,MeshlessFVParticle,OctTreeCell>;
template class Tree<3,MeshlessFVParticle,OctTreeCell>;

template class Tree<1,Particle,TreeRayCell>;
template class Tree<1,SphParticle,TreeRayCell>;
template class Tree<1,GradhSphParticle,TreeRayCell>;
template class Tree<2,Particle,TreeRayCell>;
template class Tree<2,SphParticle,TreeRayCell>;
template class Tree<2,GradhSphParticle,TreeRayCell>;
template class Tree<3,Particle,TreeRayCell>;
template class Tree<3,SphParticle,TreeRayCell>;
template class Tree<3,GradhSphParticle,TreeRayCell>;
