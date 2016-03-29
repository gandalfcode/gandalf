//=================================================================================================
//  MeshlessFVTree.cpp
//  Contains all functions for building, stocking and walking for the
//  binary KD tree for Meshless-FV particles.
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
#include "MfvNeighbourSearch.h"
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
//  MeshlessFVKDTree::MeshlessFVKDTree
/// MeshlessFVKDTree constructor.  Initialises various variables and creates tree objects.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
MeshlessFVKDTree<ndim,ParticleType,TreeCell>::MeshlessFVKDTree
 (int _Nleafmax, int _Nmpi, int _pruning_level_min, int _pruning_level_max, FLOAT _thetamaxsqd,
  FLOAT _kernrange, FLOAT _macerror, string _gravity_mac, string _multipole,
  DomainBox<ndim>* _box, SmoothingKernel<ndim>* _kern, CodeTiming* _timing):
 NeighbourSearch<ndim>(_kernrange, _box, _kern, _timing),
 MeshlessFVTree<ndim,ParticleType,TreeCell>
  (_Nleafmax, _Nmpi, _pruning_level_min, _pruning_level_max, _thetamaxsqd,
   _kernrange, _macerror, _gravity_mac, _multipole, _box, _kern, _timing)
{
  // Set-up main tree object
  tree = new KDTree<ndim,ParticleType,TreeCell>(_Nleafmax, _thetamaxsqd, _kernrange,
                                                _macerror, _gravity_mac, _multipole, *_box);

  // Set-up ghost-particle tree object
  ghosttree = new KDTree<ndim,ParticleType,TreeCell>(_Nleafmax, _thetamaxsqd, _kernrange,
                                                     _macerror, _gravity_mac, _multipole, *_box);

#ifdef MPI_PARALLEL
  // Set-up ghost-particle tree object
  mpighosttree = new KDTree<ndim,ParticleType,TreeCell>(_Nleafmax, _thetamaxsqd, _kernrange,
                                                        _macerror, _gravity_mac, _multipole, *_box);

  // Set-up multiple pruned trees, one for each MPI process
  KDTree<ndim,ParticleType,TreeCell>** prunedtree_derived = new KDTree<ndim,ParticleType,TreeCell>*[Nmpi];
  prunedtree = (Tree<ndim,ParticleType,TreeCell> **) prunedtree_derived;
  KDTree<ndim,ParticleType,TreeCell>** sendprunedtree_derived = new KDTree<ndim,ParticleType,TreeCell>*[Nmpi];
  sendprunedtree = (Tree<ndim,ParticleType,TreeCell> **) sendprunedtree_derived;

  for (int i=0; i<Nmpi; i++) {
    prunedtree[i] = new KDTree<ndim,ParticleType,TreeCell>
     (_Nleafmax, _thetamaxsqd, _kernrange, _macerror, _gravity_mac, _multipole, *_box);
  }
  for (int i=0; i<Nmpi; i++) {
    sendprunedtree[i] = new KDTree<ndim,ParticleType,TreeCell>
     (_Nleafmax, _thetamaxsqd, _kernrange, _macerror, _gravity_mac, _multipole, *_box);
  }
#endif
}



//=================================================================================================
//  MeshlessFVOctTree::MeshlessFVOctTree
/// SphTree constructor.  Initialises various variables.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
MeshlessFVOctTree<ndim,ParticleType,TreeCell>::MeshlessFVOctTree
 (int _Nleafmax, int _Nmpi, int _pruning_level_min, int _pruning_level_max, FLOAT _thetamaxsqd,
  FLOAT _kernrange, FLOAT _macerror, string _gravity_mac, string _multipole,
  DomainBox<ndim>* _box, SmoothingKernel<ndim>* _kern, CodeTiming* _timing):
 NeighbourSearch<ndim>(_kernrange, _box, _kern, _timing),
 MeshlessFVTree<ndim,ParticleType,TreeCell>
  (_Nleafmax, _Nmpi, _pruning_level_min, _pruning_level_max, _thetamaxsqd,
   _kernrange, _macerror, _gravity_mac, _multipole, _box, _kern, _timing)
{
  // Set-up main tree object
  tree = new OctTree<ndim,ParticleType,TreeCell>(_Nleafmax, _thetamaxsqd, _kernrange,
                                                 _macerror, _gravity_mac, _multipole, *_box);

  // Set-up ghost-particle tree object
  ghosttree = new OctTree<ndim,ParticleType,TreeCell>(_Nleafmax, _thetamaxsqd, _kernrange,
                                                      _macerror, _gravity_mac, _multipole, *_box);

#ifdef MPI_PARALLEL
  // Set-up ghost-particle tree object
  mpighosttree = new OctTree<ndim,ParticleType,TreeCell>(_Nleafmax, _thetamaxsqd, _kernrange,
                                                         _macerror, _gravity_mac, _multipole, *_box);

  // Set-up multiple pruned trees, one for each MPI process
  //*(prunedtree) = *(new OctTree<ndim,ParticleType,TreeCell>*[Nmpi]);
  // Set-up multiple pruned trees, one for each MPI process
  OctTree<ndim,ParticleType,TreeCell>** prunedtree_derived = new OctTree<ndim,ParticleType,TreeCell>*[Nmpi];
  prunedtree = (Tree<ndim,ParticleType,TreeCell> **) prunedtree_derived;
  OctTree<ndim,ParticleType,TreeCell>** sendprunedtree_derived = new OctTree<ndim,ParticleType,TreeCell>*[Nmpi];
  sendprunedtree = (Tree<ndim,ParticleType,TreeCell> **) sendprunedtree_derived;

  for (int j=0; j<Nmpi; j++) {
    prunedtree[j] = new OctTree<ndim,ParticleType,TreeCell>
     (_Nleafmax, _thetamaxsqd, _kernrange, _macerror, _gravity_mac, _multipole, *_box);
  }
  for (int i=0; i<Nmpi; i++) {
    sendprunedtree[i] = new OctTree<ndim,ParticleType,TreeCell>
     (_Nleafmax, _thetamaxsqd, _kernrange, _macerror, _gravity_mac, _multipole, *_box);
  }
#endif
}



//=================================================================================================
//  MeshlessFVTree::MeshlessFVTree
/// MeshlessFVTree constructor.  Initialises various variables.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
MeshlessFVTree<ndim,ParticleType,TreeCell>::MeshlessFVTree
 (int _Nleafmax, int _Nmpi, int _pruning_level_min, int _pruning_level_max, FLOAT _thetamaxsqd,
  FLOAT _kernrange, FLOAT _macerror, string _gravity_mac, string _multipole,
  DomainBox<ndim>* _box, SmoothingKernel<ndim>* _kern, CodeTiming* _timing):
 NeighbourSearch<ndim>(_kernrange, _box, _kern, _timing),
 MeshlessFVNeighbourSearch<ndim>(_kernrange, _box, _kern, _timing),
 HydroTree<ndim,ParticleType,TreeCell>
  (_Nleafmax, _Nmpi, _pruning_level_min, _pruning_level_max, _thetamaxsqd,
   _kernrange, _macerror, _gravity_mac, _multipole, _box, _kern, _timing)
{
}



//=================================================================================================
//  MeshlessFVTree::~MeshlessFVTree
/// MeshlessFVTree destructor.  Deallocates tree memory upon object destruction.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
MeshlessFVTree<ndim,ParticleType,TreeCell>::~MeshlessFVTree()
{
  if (tree->allocated_tree) {
    this->DeallocateMemory();
    tree->DeallocateTreeMemory();
  }
}



//=================================================================================================
//  MeshlessFVTree::UpdateAllSphProperties
/// Update all gather SPH properties (e.g. rho, div_v) for all active particles in domain.
/// Loops over all cells containing active particles, performs a tree walk for all particles in
/// the cell, and then calls SPH class routine to compute properties from neighbours.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void MeshlessFVTree<ndim,ParticleType,TreeCell>::UpdateAllProperties
 (int Nhydro,                              ///< [in] No. of SPH particles
  int Ntot,                                ///< [in] No. of SPH + ghost particles
  MeshlessFVParticle<ndim> *mfvdata,       ///< [inout] Pointer to SPH ptcl array
  MeshlessFV<ndim> *mfv,                   ///< [in] Pointer to SPH object
  Nbody<ndim> *nbody,                      ///< [in] Pointer to N-body object
  DomainBox<ndim> &simbox)                 ///< [in] Simulation domain box
{
  int cactive;                             // No. of active tree cells
  TreeCell<ndim> *celllist;                // List of active tree cells
  //ParticleType<ndim> *partdata = static_cast<ParticleType<ndim>* > (sph_gen);
#ifdef MPI_PARALLEL
  int Nactivetot = 0;                      // Total number of active particles
  double twork = timing->WallClockTime();  // Start time (for load balancing)
#endif

  debug2("[MeshlessFVTree::UpdateAllProperties]");
  timing->StartTimingSection("MFV_PROPERTIES");


  // Find list of all cells that contain active particles
  celllist = new TreeCell<ndim>[tree->gtot];
  cactive = tree->ComputeActiveCellList(celllist);

  // If there are no active cells, return to main loop
  if (cactive == 0) {
    delete[] celllist;
    timing->EndTimingSection("MFV_PROPERTIES");
    return;
  }


  // Set-up all OMP threads
  //===============================================================================================
#pragma omp parallel default(none) shared(cactive,celllist,cout,nbody,mfv,mfvdata,Ntot)
  {
#if defined _OPENMP
    const int ithread = omp_get_thread_num();
#else
    const int ithread = 0;
#endif
    int celldone;                              // Flag if cell is done
    int cc;                                    // Aux. cell counter
    int i;                                     // Particle id
    int j;                                     // Aux. particle counter
    int jj;                                    // Aux. particle counter
    int k;                                     // Dimension counter
    int Nactive;                               // No. of active particles in cell
    int Ngather;                               // No. of gather neighbours
    int Nneib;                                 // No. of neighbours from tree-walk
    int okflag;                                // Flag if particle is done
    FLOAT draux[ndim];                         // Aux. relative position vector var
    FLOAT drsqdaux;                            // Distance squared
    FLOAT hrangesqd;                           // Kernel extent
    FLOAT hmax;                                // Maximum smoothing length
    FLOAT rp[ndim];                            // Local copy of particle position
    //FLOAT *mu,                                 // Mass times specific internal energy arrays
    FLOAT *mu2 = 0;                            // Trimmed array (dummy for grad-h)
    int Nneibmax = Nneibmaxbuf[ithread];       // Local copy of neighbour buffer size
    int* activelist = activelistbuf[ithread];  // Local array of active particle ids
    int* neiblist = new int[Nneibmax];         // Local array of neighbour particle ids
    int* ptype    = new int[Nneibmax];         // Local array of particle types
    FLOAT* gpot   = new FLOAT[Nneibmax];       // Local array of particle potentials
    FLOAT* gpot2  = new FLOAT[Nneibmax];       // Local reduced array of neighbour potentials
    FLOAT* drsqd  = new FLOAT[Nneibmax];       // Local array of distances (squared)
    FLOAT* m      = new FLOAT[Nneibmax];       // Local array of particle masses
    FLOAT* m2     = new FLOAT[Nneibmax];       // Local reduced array of neighbour masses
    FLOAT* r      = new FLOAT[Nneibmax*ndim];                  // Local array of particle positions
    ParticleType<ndim>* activepart = activepartbuf[ithread];   // Local array of active particles


    // Loop over all active cells
    //=============================================================================================
#pragma omp for schedule(guided)
    for (cc=0; cc<cactive; cc++) {
      TreeCell<ndim>& cell = celllist[cc];
      celldone = 1;
      hmax = cell.hmax;


      // If hmax is too small so the neighbour lists are invalid, make hmax
      // larger and then recompute for the current active cell.
      //-------------------------------------------------------------------------------------------
      do {
        hmax = (FLOAT) 1.05*hmax;
        celldone = 1;

        // Find list of active particles in current cell
        Nactive = tree->ComputeActiveParticleList(cell,mfvdata,activelist);
        for (j=0; j<Nactive; j++) activepart[j] = mfvdata[activelist[j]];

        // Compute neighbour list for cell from particles on all trees
        Nneib = 0;
        Nneib = tree->ComputeGatherNeighbourList(cell,mfvdata,hmax,Nneibmax,Nneib,neiblist);
        Nneib = ghosttree->ComputeGatherNeighbourList(cell,mfvdata,hmax,Nneibmax,Nneib,neiblist);
#ifdef MPI_PARALLEL
        Nneib = mpighosttree->ComputeGatherNeighbourList(cell,mfvdata,hmax,Nneibmax,Nneib,neiblist);
#endif

        // If there are too many neighbours so the buffers are filled,
        // reallocate the arrays and recompute the neighbour lists.
        while (Nneib == -1) {
          delete[] r;
          delete[] m2;
          delete[] m;
          delete[] drsqd;
          delete[] gpot2;
          delete[] gpot;
          delete[] ptype;
          delete[] neiblist;
          Nneibmax = 2*Nneibmax;
          neiblist = new int[Nneibmax];
          ptype    = new int[Nneibmax];
          gpot     = new FLOAT[Nneibmax];
          gpot2    = new FLOAT[Nneibmax];
          drsqd    = new FLOAT[Nneibmax];
          m        = new FLOAT[Nneibmax];
          m2       = new FLOAT[Nneibmax];
          r        = new FLOAT[Nneibmax*ndim];
          Nneib = 0;
          Nneib = tree->ComputeGatherNeighbourList(cell,mfvdata,hmax,Nneibmax,Nneib,neiblist);
          Nneib = ghosttree->ComputeGatherNeighbourList(cell,mfvdata,hmax,Nneibmax,Nneib,neiblist);
#ifdef MPI_PARALLEL
          Nneib = mpighosttree->ComputeGatherNeighbourList(cell,mfvdata,hmax,
                                                           Nneibmax,Nneib,neiblist);
#endif
        };


        // Make local copies of important neib information (mass and position)
        for (jj=0; jj<Nneib; jj++) {
          j         = neiblist[jj];
          gpot[jj]  = mfvdata[j].gpot;
          m[jj]     = mfvdata[j].m;
          ptype[jj] = mfvdata[j].ptype;
          for (k=0; k<ndim; k++) r[ndim*jj + k] = mfvdata[j].r[k];
        }


        // Loop over all active particles in the cell
        //-----------------------------------------------------------------------------------------
        for (j=0; j<Nactive; j++) {
          i = activelist[j];
          for (k=0; k<ndim; k++) rp[k] = activepart[j].r[k];

          // Set gather range as current h multiplied by some tolerance factor
          hrangesqd = kernrangesqd*hmax*hmax;
          Ngather = 0;

          // Compute distance (squared) to all
          //---------------------------------------------------------------------------------------
          for (jj=0; jj<Nneib; jj++) {

            // Only include particles of appropriate types in density calculation
            if (!mfv->types[activepart[j].ptype].hmask[ptype[jj]]) continue ;

            for (k=0; k<ndim; k++) draux[k] = r[ndim*jj + k] - rp[k];
            drsqdaux = DotProduct(draux,draux,ndim) + small_number;

            // Record distance squared and masses for all potential gather neighbours
            if (drsqdaux <= hrangesqd) {
              gpot[Ngather]  = gpot[jj];
              drsqd[Ngather] = drsqdaux;
              m2[Ngather]    = m[jj];
              Ngather++;
            }

          }
          //---------------------------------------------------------------------------------------

          // Validate that gather neighbour list is correct
#if defined(VERIFY_ALL)
          if (neibcheck) this->CheckValidNeighbourList(i, Ntot, Nneib, neiblist, mfvdata, "gather");
#endif

          // Compute smoothing length and other gather properties for ptcl i
          okflag = mfv->ComputeH(i, Ngather, hmax, m2, mu2, drsqd, gpot, activepart[j], nbody);

          // If h-computation is invalid, then break from loop and recompute larger neighbour lists
          if (okflag == 0) {
            celldone = 0;
            break;
          }

        }
        //-----------------------------------------------------------------------------------------

      } while (celldone == 0);
      //-------------------------------------------------------------------------------------------

      // Once cell is finished, copy all active particles back to main memory
      for (j=0; j<Nactive; j++) mfvdata[activelist[j]] = activepart[j];


    }
    //=============================================================================================

    // Free-up all memory
    delete[] r;
    delete[] m2;
    delete[] m;
    delete[] drsqd;
    delete[] gpot2;
    delete[] gpot;
    delete[] ptype;
    delete[] neiblist;

  }
  //===============================================================================================

  // Compute time spent in routine and in each cell for load balancing
#ifdef MPI_PARALLEL
  twork = timing->WallClockTime() - twork;
  for (int cc=0; cc<cactive; cc++) Nactivetot += celllist[cc].Nactive;
  for (int cc=0; cc<cactive; cc++) {
    int c = celllist[cc].id;
    tree->celldata[c].worktot += twork*(DOUBLE) tree->celldata[c].Nactive / (DOUBLE) Nactivetot;
  }
  cout << "Time computing smoothing lengths : " << twork << "     Nactivetot : " << Nactivetot << endl;
#endif

  delete[] celllist;

  // Update tree smoothing length values here
  tree->UpdateHmaxValues(tree->celldata[0],mfvdata);

  timing->EndTimingSection("MFV_PROPERTIES");

  return;
}



//=================================================================================================
//  MeshlessFVTree::UpdateGradientMatrices
/// Compute hydro forces for all active SPH particles.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void MeshlessFVTree<ndim,ParticleType,TreeCell>::UpdateGradientMatrices
 (int Nhydro,                              ///< [in] No. of SPH particles
  int Ntot,                                ///< [in] No. of SPH + ghost particles
  MeshlessFVParticle<ndim> *mfvdata,       ///< [inout] Pointer to SPH ptcl array
  MeshlessFV<ndim> *mfv,                   ///< [in] Pointer to SPH object
  Nbody<ndim> *nbody,                      ///< [in] Pointer to N-body object
  DomainBox<ndim> &simbox)                 ///< [in] Simulation domain box
{
  int cactive;                             // No. of active cells
  TreeCell<ndim> *celllist;                // List of active tree cells
#ifdef MPI_PARALLEL
  int Nactivetot = 0;                      // Total number of active particles
  double twork = timing->WallClockTime();  // Start time (for load balancing)
#endif

  debug2("[MeshlessFVTree::UpdateGradientMatrices]");
  timing->StartTimingSection("MFV_UPDATE_GRADIENTS");


  // Find list of all cells that contain active particles
#if defined (MPI_PARALLEL)
  celllist = new TreeCell<ndim>[tree->Ncellmax];
#else
  celllist = new TreeCell<ndim>[tree->gtot];
#endif
  cactive = tree->ComputeActiveCellList(celllist);

  // If there are no active cells, return to main loop
  if (cactive == 0) {
    delete[] celllist;
    timing->EndTimingSection("MFV_UPDATE_GRADIENTS");
    return;
  }

  // Update ghost tree smoothing length values here
  tree->UpdateHmaxValues(tree->celldata[0], mfvdata);
  if (ghosttree->Ntot > 0) ghosttree->UpdateHmaxValues(ghosttree->celldata[0], mfvdata);


  // Set-up all OMP threads
  //===============================================================================================
#pragma omp parallel default(none) shared(cactive,celllist,nbody,mfv,mfvdata,Ntot)
  {
#if defined _OPENMP
    const int ithread = omp_get_thread_num();
#else
    const int ithread = 0;
#endif
    int cc;                                        // Aux. cell counter
    int i;                                         // Particle id
    int j;                                         // Aux. particle counter
    int jj;                                        // Aux. particle counter
    int k;                                         // Dimension counter
    int Nactive;                                   // ..
    int Nneib;                                     // ..
    int Nhydroaux;                                 // ..
    FLOAT draux[ndim];                             // Aux. relative position vector
    FLOAT drsqd;                                   // Distance squared
    FLOAT hrangesqdi;                              // Kernel gather extent
    FLOAT rp[ndim];                                // Local copy of particle position
    Typemask hmask;                                // Neigbour mask for computing gradients
    int Nneibmax    = Nneibmaxbuf[ithread];        // ..
    int* activelist = activelistbuf[ithread];      // ..
    int* levelneib  = levelneibbuf[ithread];       // ..
    int* neiblist   = new int[Nneibmax];           // ..
    int* mfvlist    = new int[Nneibmax];           // ..
    FLOAT* dr       = new FLOAT[Nneibmax*ndim];    // ..
    FLOAT* drmag    = new FLOAT[Nneibmax];         // ..
    FLOAT* invdrmag = new FLOAT[Nneibmax];         // ..
    ParticleType<ndim>* activepart = activepartbuf[ithread];   // ..
    ParticleType<ndim>* neibpart   = neibpartbuf[ithread];     // ..

    for (i=0; i<mfv->Nhydro; i++) levelneib[i] = 0;


    // Loop over all active cells
    //=============================================================================================
#pragma omp for schedule(guided)
    for (cc=0; cc<cactive; cc++) {
      TreeCell<ndim>& cell = celllist[cc];

      // Find list of active particles in current cell
      Nactive = tree->ComputeActiveParticleList(cell, mfvdata, activelist);

      // Make local copies of active particles
      for (j=0; j<Nactive; j++) activepart[j] = mfvdata[activelist[j]];

      // Compute neighbour list for cell from real and periodic ghost particles
      Nneib = 0;
      Nneib = tree->ComputeNeighbourAndGhostList
        (cell, mfvdata, Nneibmax, Nneib, neiblist, neibpart);

      // If there are too many neighbours, reallocate the arrays and
      // recompute the neighbour list.
      while (Nneib == -1) {
        delete[] neibpartbuf[ithread];
        delete[] invdrmag;
        delete[] drmag;
        delete[] dr;
        delete[] mfvlist;
        delete[] neiblist;
        Nneibmax                  = 2*Nneibmax;
        Nneibmaxbuf[ithread]      = Nneibmax;
        Ngravcellmaxbuf[ithread] *= 2;
        neiblist                  = new int[Nneibmax];
        mfvlist                   = new int[Nneibmax];
        dr                        = new FLOAT[Nneibmax*ndim];
        drmag                     = new FLOAT[Nneibmax];
        invdrmag                  = new FLOAT[Nneibmax];
        neibpartbuf[ithread]      = new ParticleType<ndim>[Nneibmax];
        neibpart                  = neibpartbuf[ithread];
        Nneib = 0;
        Nneib = tree->ComputeNeighbourAndGhostList
          (cell, mfvdata, Nneibmax, Nneib, neiblist, neibpart);
      };


      // Loop over all active particles in the cell
      //-------------------------------------------------------------------------------------------
      for (j=0; j<Nactive; j++) {
        i = activelist[j];

        // If particle is NOT a hydro particle (and therefore doesn't need gradients), skip to next
        if (!mfv->types[activepart[j].ptype].hydro_forces) continue;

        // Make local copy of hmask for active particle
        for (k=0; k<Ntypes; k++) hmask[k] = mfv->types[activepart[j].ptype].hmask[k];

        for (k=0; k<ndim; k++) rp[k] = activepart[j].r[k];
        hrangesqdi = activepart[j].hrangesqd;
        Nhydroaux = 0;

        // Validate that gather neighbour list is correct
#if defined(VERIFY_ALL)
        if (neibcheck) this->CheckValidNeighbourList(i, Ntot, Nneib, neiblist, mfvdata, "all");
#endif

        // Compute distances and the inverse between the current particle and all neighbours here,
        // for both gather and inactive scatter neibs.  Only consider particles with j > i to
        // compute pair forces once unless particle j is inactive.
        //-----------------------------------------------------------------------------------------
        for (jj=0; jj<Nneib; jj++) {

          // Only include particles of appropriate types in density calculation
          if (hmask[neibpart[jj].ptype] == false) continue ;

          // Skip current active particle
          if (neiblist[jj] == i) continue;

          for (k=0; k<ndim; k++) draux[k] = neibpart[jj].r[k] - rp[k];
          drsqd = DotProduct(draux,draux,ndim) + small_number;

          // Compute relative position and distance quantities for pair
          if (drsqd <= hrangesqdi || drsqd <= neibpart[jj].hrangesqd) {
            drmag[Nhydroaux] = sqrt(drsqd);
            invdrmag[Nhydroaux] = (FLOAT) 1.0/drmag[Nhydroaux];
            for (k=0; k<ndim; k++) dr[Nhydroaux*ndim + k] = draux[k]*invdrmag[Nhydroaux];
            levelneib[neiblist[jj]] = max(levelneib[neiblist[jj]],activepart[j].level);
            mfvlist[Nhydroaux] = jj;
            Nhydroaux++;
          }

        }
        //-----------------------------------------------------------------------------------------

        // Compute all neighbour contributions to hydro forces
        mfv->ComputePsiFactors(i, Nhydroaux, mfvlist, drmag, invdrmag, dr, activepart[j], neibpart);
        mfv->ComputeGradients(i, Nhydroaux, mfvlist, drmag, invdrmag, dr, activepart[j], neibpart);

      }
      //-------------------------------------------------------------------------------------------


      // Add all active particles contributions to main array
      for (j=0; j<Nactive; j++) {
        i = activelist[j];
        for (k=0; k<ndim; k++) {
          for (int kk=0; kk<ndim; kk++) mfvdata[i].B[k][kk] = activepart[j].B[k][kk];
        }
        for (int var=0; var<ndim+2; var++) {
          for (k=0; k<ndim; k++) mfvdata[i].grad[var][k] = activepart[j].grad[var][k];
          mfvdata[i].Wmin[var] = activepart[j].Wmin[var];
          mfvdata[i].Wmax[var] = activepart[j].Wmax[var];
          mfvdata[i].Wmidmin[var] = activepart[j].Wmidmin[var];
          mfvdata[i].Wmidmax[var] = activepart[j].Wmidmax[var];
        }
        mfvdata[i].vsig_max = activepart[j].vsig_max;
      }


    }
    //=============================================================================================


    // Free-up local memory for OpenMP thread
    delete[] invdrmag;
    delete[] drmag;
    delete[] dr;
    delete[] mfvlist;
    delete[] neiblist;

  }
  //===============================================================================================


  // Compute time spent in routine and in each cell for load balancing
#ifdef MPI_PARALLEL
  twork = timing->WallClockTime() - twork;
  for (int cc=0; cc<cactive; cc++) Nactivetot += celllist[cc].Nactive;
  for (int cc=0; cc<cactive; cc++) {
    int c = celllist[cc].id;
    tree->celldata[c].worktot += twork*(DOUBLE) tree->celldata[c].Nactive / (DOUBLE) Nactivetot;
  }
  cout << "Time computing forces : " << twork << "     Nactivetot : " << Nactivetot << endl;
#endif


  delete[] celllist;

  timing->EndTimingSection("MFV_UPDATE_GRADIENTS");

  return;
}



//=================================================================================================
//  MeshlessFVTree::UpdateGodunovFluxes
/// Compute hydro forces for all active SPH particles.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void MeshlessFVTree<ndim,ParticleType,TreeCell>::UpdateGodunovFluxes
 (const int Nhydro,                        ///< [in] No. of hydro particles
  const int Ntot,                          ///< [in] No. of SPH + ghost particles
  const FLOAT timestep,                    ///< [in] Lowest timestep value
  MeshlessFVParticle<ndim> *mfvdata,       ///< [inout] Pointer to SPH ptcl array
  MeshlessFV<ndim> *mfv,                   ///< [in] Pointer to SPH object
  Nbody<ndim> *nbody,                      ///< [in] Pointer to N-body object
  DomainBox<ndim> &simbox)                 ///< [in] Simulation domain box
{
  int cactive;                             // No. of active cells
  TreeCell<ndim> *celllist;                // List of active tree cells
#ifdef MPI_PARALLEL
  int Nactivetot = 0;                      // Total number of active particles
  double twork = timing->WallClockTime();  // Start time (for load balancing)
#endif

  debug2("[MeshlessFVTree::UpdateGradientMatrices]");
  timing->StartTimingSection("MFV_UPDATE_FLUXES");


  // Find list of all cells that contain active particles
#if defined (MPI_PARALLEL)
  celllist = new TreeCell<ndim>[tree->Ncellmax];
#else
  celllist = new TreeCell<ndim>[tree->gtot];
#endif
  cactive = tree->ComputeActiveCellList(celllist);

  // If there are no active cells, return to main loop
  if (cactive == 0) {
    delete[] celllist;
    timing->EndTimingSection("MFV_UPDATE_FLUXES");
    return;
  }

  // Update ghost tree smoothing length values here
  tree->UpdateHmaxValues(tree->celldata[0], mfvdata);
  //if (ghosttree->Ntot > 0) ghosttree->UpdateHmaxValues(ghosttree->celldata[0], mfvdata);


  // Set-up all OMP threads
  //===============================================================================================
#pragma omp parallel default(none) shared(cactive,celllist,mfv,mfvdata)
  {
#if defined _OPENMP
    const int ithread = omp_get_thread_num();
#else
    const int ithread = 0;
#endif
    int cc;                                        // Aux. cell counter
    int i;                                         // Particle id
    int j;                                         // Aux. particle counter
    int jj;                                        // Aux. particle counter
    int k;                                         // Dimension counter
    int Nactive;                                   // ..
    int Nneib;                                     // ..
    int Nhydroaux;                                 // ..
    FLOAT draux[ndim];                             // Aux. relative position vector
    FLOAT drsqd;                                   // Distance squared
    FLOAT hrangesqdi;                              // Kernel gather extent
    FLOAT rp[ndim];                                // Local copy of particle position
    Typemask hydromask;                            // Mask for computing hydro forces
    int Nneibmax      = Nneibmaxbuf[ithread];      // ..
    int* activelist   = activelistbuf[ithread];    // ..
    int* levelneib    = levelneibbuf[ithread];     // ..
    int* neiblist     = new int[Nneibmax];         // ..
    int* mfvlist      = new int[Nneibmax];         // ..
    FLOAT* dr         = new FLOAT[Nneibmax*ndim];  // ..
    FLOAT* drmag      = new FLOAT[Nneibmax];       // ..
    FLOAT* invdrmag   = new FLOAT[Nneibmax];       // ..
    FLOAT (*fluxBuffer)[ndim+2]    = new FLOAT[Ntot][ndim+2];  // ..
    FLOAT (*rdmdtBuffer)[ndim]     = new FLOAT[Ntot][ndim];    // ..
    ParticleType<ndim>* activepart = activepartbuf[ithread];   // ..
    ParticleType<ndim>* neibpart   = neibpartbuf[ithread];     // ..

    for (i=0; i<mfv->Nhydro; i++) levelneib[i] = 0;
    for (i=0; i<Ntot; i++) {
      for (k=0; k<ndim+2; k++) fluxBuffer[i][k] = (FLOAT) 0.0;
      for (k=0; k<ndim; k++) rdmdtBuffer[i][k] = (FLOAT) 0.0;
    }


    // Loop over all active cells
    //=============================================================================================
#pragma omp for schedule(guided)
    for (cc=0; cc<cactive; cc++) {
      TreeCell<ndim>& cell = celllist[cc];

      // Find list of active particles in current cell
      Nactive = tree->ComputeActiveParticleList(cell,mfvdata,activelist);

      // Make local copies of active particles
      for (j=0; j<Nactive; j++) {
        activepart[j] = mfvdata[activelist[j]];
        for (k=0; k<ndim+2; k++) activepart[j].dQ[k]   = (FLOAT) 0.0;
        for (k=0; k<ndim+2; k++) activepart[j].dQdt[k] = (FLOAT) 0.0;
        for (k=0; k<ndim; k++) activepart[j].rdmdt[k]  = (FLOAT) 0.0;
      }

      // Compute neighbour list for cell from real and periodic ghost particles
      Nneib = 0;
      Nneib = tree->ComputeNeighbourAndGhostList
        (cell, mfvdata, Nneibmax, Nneib, neiblist, neibpart);

      // If there are too many neighbours, reallocate the arrays and
      // recompute the neighbour list.
      while (Nneib == -1) {
        delete[] neibpartbuf[ithread];
        delete[] invdrmag;
        delete[] drmag;
        delete[] dr;
        delete[] mfvlist;
        delete[] neiblist;
        Nneibmax                  = 2*Nneibmax;
        Nneibmaxbuf[ithread]      = Nneibmax;
        Ngravcellmaxbuf[ithread] *= 2;
        neiblist                  = new int[Nneibmax];
        mfvlist                   = new int[Nneibmax];
        dr                        = new FLOAT[Nneibmax*ndim];
        drmag                     = new FLOAT[Nneibmax];
        invdrmag                  = new FLOAT[Nneibmax];
        neibpartbuf[ithread]      = new ParticleType<ndim>[Nneibmax];
        neibpart                  = neibpartbuf[ithread];
        Nneib = 0;
        Nneib = tree->ComputeNeighbourAndGhostList
          (cell, mfvdata, Nneibmax, Nneib, neiblist, neibpart);
      };

      for (j=0; j<Nneibmax; j++) {
        for (k=0; k<ndim+2; k++) neibpart[j].dQ[k]   = (FLOAT) 0.0;
        for (k=0; k<ndim+2; k++) neibpart[j].dQdt[k] = (FLOAT) 0.0;
        for (k=0; k<ndim; k++) neibpart[j].rdmdt[k]  = (FLOAT) 0.0;
      }


      // Loop over all active particles in the cell
      //-------------------------------------------------------------------------------------------
      for (j=0; j<Nactive; j++) {
        i = activelist[j];

        // If particle is not a hydro particle (e.g. cdm), then skip to next active particle
        if (mfv->types[activepart[j].ptype].hydro_forces == false) continue;

        // Make a local copy of the hydro neighbour mask
        for (k=0; k<Ntypes; k++) hydromask[k] = mfv->types[activepart[j].ptype].hydromask[k];

        for (k=0; k<ndim; k++) rp[k] = activepart[j].r[k];
        hrangesqdi = activepart[j].hrangesqd;
        Nhydroaux = 0;

        // Validate that gather neighbour list is correct
#if defined(VERIFY_ALL)
        if (neibcheck) this->CheckValidNeighbourList(i, Ntot, Nneib, neiblist, mfvdata, "all");
#endif

        // Compute distances and the inverse between the current particle and all neighbours here,
        // for both gather and inactive scatter neibs.  Only consider particles with j > i to
        // compute pair forces once unless particle j is inactive.
        //-----------------------------------------------------------------------------------------
        for (jj=0; jj<Nneib; jj++) {

          // Skip if (i) neighbour particle type does not interact hydrodynamically with particle,
          // (ii) neighbour is a dead (e.g. accreted) particle (iii) same i.d. as current active
          // particle, (iv) neighbour is on lower timestep level (i.e. timestep is shorter),
          // or (v) neighbour is on same level as current particle but has larger id. value
          // (to only calculate each pair once).
          if (hydromask[neibpart[jj].ptype] == false || neibpart[jj].itype == dead ||
              neiblist[jj] == i || activepart[j].level < neibpart[jj].level ||
              (neibpart[jj].iorig < i && neibpart[jj].level == activepart[j].level)) continue;

          // Compute relative position and distance quantities for pair
          for (k=0; k<ndim; k++) draux[k] = neibpart[jj].r[k] - rp[k];
          drsqd = DotProduct(draux, draux, ndim) + small_number;

          // Only include gather or scatter neighbours
          if (drsqd < hrangesqdi || drsqd < neibpart[jj].hrangesqd) {
            drmag[Nhydroaux] = sqrt(drsqd);
            invdrmag[Nhydroaux] = (FLOAT) 1.0/drmag[Nhydroaux];
            for (k=0; k<ndim; k++) dr[Nhydroaux*ndim + k] = draux[k]*invdrmag[Nhydroaux];
            levelneib[neiblist[jj]] = max(levelneib[neiblist[jj]], activepart[j].level);
            mfvlist[Nhydroaux] = jj;
            Nhydroaux++;
          }

        }
        //-----------------------------------------------------------------------------------------

        // Compute all neighbour contributions to hydro fluxes
        mfv->ComputeGodunovFlux(i, Nhydroaux, timestep, mfvlist, drmag,
                                invdrmag, dr, activepart[j], neibpart);

      }
      //-------------------------------------------------------------------------------------------


      // Accumulate fluxes for neighbours (only currently works for real and periodic neighbours)
      for (int jj=0; jj<Nneib; jj++) {
        i = neibpart[jj].iorig;
        for (k=0; k<ndim+2; k++) fluxBuffer[i][k] += neibpart[jj].dQ[k];
        for (k=0; k<ndim; k++) rdmdtBuffer[i][k] += neibpart[jj].rdmdt[k];
      }
      /*for (int jj=0; jj<Nneib; jj++) {
        j = neiblist[jj];
        for (k=0; k<ndim+2; k++) fluxBuffer[j][k] += neibpart[jj].dQ[k];
        for (k=0; k<ndim; k++) rdmdtBuffer[j][k] += neibpart[jj].rdmdt[k];
      }*/

      // Add all active particles contributions to main array
      for (j=0; j<Nactive; j++) {
        i = activelist[j];
        for (k=0; k<ndim+2; k++) fluxBuffer[i][k] += activepart[j].dQ[k];
        for (k=0; k<ndim; k++) rdmdtBuffer[i][k] += activepart[j].rdmdt[k];
      }

    }
    //=============================================================================================


    // Add all buffers back to main arrays
#pragma omp barrier
#pragma omp critical
    {
      for (i=0; i<Nhydro; i++) {
        for (k=0; k<ndim+2; k++) mfvdata[i].dQ[k] += fluxBuffer[i][k];
        for (k=0; k<ndim; k++) mfvdata[i].rdmdt[k] += rdmdtBuffer[i][k];
      }
      for (j=Nhydro; j<Nhydro+mfv->NPeriodicGhost; j++) {
        i = mfvdata[j].iorig;
        for (k=0; k<ndim+2; k++) mfvdata[i].dQ[k] += fluxBuffer[j][k];
        for (k=0; k<ndim; k++) mfvdata[i].rdmdt[k] += rdmdtBuffer[j][k];
      }
    }


    // Free-up local memory for OpenMP thread
    delete[] rdmdtBuffer;
    delete[] fluxBuffer;
    delete[] invdrmag;
    delete[] drmag;
    delete[] dr;
    delete[] mfvlist;
    delete[] neiblist;

  }
  //===============================================================================================


  // Compute time spent in routine and in each cell for load balancing
#ifdef MPI_PARALLEL
  twork = timing->WallClockTime() - twork;
  for (int cc=0; cc<cactive; cc++) Nactivetot += celllist[cc].Nactive;
  for (int cc=0; cc<cactive; cc++) {
    int c = celllist[cc].id;
    tree->celldata[c].worktot += twork*(DOUBLE) tree->celldata[c].Nactive / (DOUBLE) Nactivetot;
  }
  cout << "Time computing forces : " << twork << "     Nactivetot : " << Nactivetot << endl;
#endif


  delete[] celllist;

  timing->EndTimingSection("MFV_UPDATE_FLUXES");

  return;
}



//=================================================================================================
//  MeshlessFVTree::UpdateAllGravForces
/// Compute all local 'gather' properties of currently active particles, and
/// then compute each particle's contribution to its (active) neighbour
/// neighbour hydro forces.  Optimises the algorithm by using grid-cells to
/// construct local neighbour lists for all particles  inside the cell.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void MeshlessFVTree<ndim,ParticleType,TreeCell>::UpdateAllGravForces
 (int Nhydro,                          ///< [in] No. of SPH particles
  int Ntot,                            ///< [in] No. of SPH + ghost particles
  MeshlessFVParticle<ndim> *part_gen,  ///< [inout] Pointer to SPH ptcl array
  MeshlessFV<ndim> *mfv,               ///< [in] Pointer to SPH object
  Nbody<ndim> *nbody,                  ///< [in] Pointer to N-body object
  DomainBox<ndim> &simbox,             ///< [in] Simulation domain box
  Ewald<ndim> *ewald)                  ///< [in] Ewald gravity object pointer
{
  int cactive;                         // No. of active cells
  TreeCell<ndim> *celllist;            // List of active tree cells
  ParticleType<ndim>* partdata = static_cast<ParticleType<ndim>* > (part_gen);

  debug2("[MeshlessFVTree::UpdateAllGravForces]");
  timing->StartTimingSection("MFV_GRAV_FORCES");

  // Update ghost tree smoothing length values here
  tree->UpdateHmaxValues(tree->celldata[0], partdata);

  // Find list of all cells that contain active particles
#if defined (MPI_PARALLEL)
  celllist = new TreeCell<ndim>[tree->Ncellmax];
#else
  celllist = new TreeCell<ndim>[tree->gtot];
#endif
  cactive = tree->ComputeActiveCellList(celllist);

  // If there are no active cells, return to main loop
  if (cactive == 0) {
    delete[] celllist;
    timing->EndTimingSection("MFV_GRAV_FORCES");
    return;
  }

  // Set-up all OMP threads
  //===============================================================================================
#pragma omp parallel default(none) shared(celllist,cactive,ewald,mfv,nbody,partdata,simbox,cout)
  {
#if defined _OPENMP
    const int ithread = omp_get_thread_num();
#else
    const int ithread = 0;
#endif
    int cc;                                      // Aux. cell counter
    int i;                                       // Particle id
    int j;                                       // Aux. particle counter
    int jj;                                      // Aux. particle counter
    int k;                                       // Dimension counter
    int okflag;                                  // Flag if h-rho iteration is valid
    int Nactive;                                 // ..
    int Ndirect;                                 // ..
    int Ndirectaux;                              // ..
    int Ngrav;                                   // --
    int Ngravcell;                               // No. of gravity cells
    int Nhydroaux;                               // ..
    int Nhydroneib;                              // ..
    int Nneib;                                   // No. of neighbours
    FLOAT aperiodic[ndim];                       // ..
    FLOAT draux[ndim];                           // Aux. relative position vector
    FLOAT drsqd;                                 // Distance squared
    FLOAT hrangesqdi;                            // Kernel gather extent
    FLOAT macfactor;                             // Gravity MAC factor
    FLOAT potperiodic;                           // ..
    FLOAT rp[ndim];                              // ..
    int Nneibmax     = Nneibmaxbuf[ithread];     // ..
    int Ngravcellmax = Ngravcellmaxbuf[ithread]; // ..
    int *activelist  = activelistbuf[ithread];   // ..
    int *levelneib   = levelneibbuf[ithread];    // ..
    int *neiblist    = new int[Nneibmax];        // ..
    int *mfvlist     = new int[Nneibmax];        // ..
    int *mfvauxlist  = new int[Nneibmax];        // ..
    int *directlist  = new int[Nneibmax];        // ..
    int	*gravlist    = new int[Nneibmax];        // ..
    ParticleType<ndim>* activepart = activepartbuf[ithread];   // ..
    ParticleType<ndim>* neibpart   = neibpartbuf[ithread];     // ..
    TreeCell<ndim>* gravcell       = cellbuf[ithread];         // ..

    // Zero timestep level array
    for (i=0; i<mfv->Nhydro; i++) levelneib[i] = 0;


    // Loop over all active cells
    //=============================================================================================
#pragma omp for schedule(guided)
    for (cc=0; cc<cactive; cc++) {
      TreeCell<ndim>& cell = celllist[cc];
      macfactor = (FLOAT) 0.0;

      // Find list of active particles in current cell
      Nactive = tree->ComputeActiveParticleList(cell, partdata, activelist);

      // Make local copies of active particles
      for (j=0; j<Nactive; j++) activepart[j] = partdata[activelist[j]];

      // Compute average/maximum term for computing gravity MAC
      if (gravity_mac == "eigenmac") {
        for (j=0; j<Nactive; j++)
          macfactor = max(macfactor, pow((FLOAT) 1.0/activepart[j].gpot, twothirds));
      }

      // Zero/initialise all summation variables for active particles
      for (j=0; j<Nactive; j++) {
        activepart[j].levelneib = 0;
        activepart[j].gpot      = activepart[j].m*activepart[j].invh*mfv->kernp->wpot((FLOAT) 0.0);
        for (k=0; k<ndim; k++) activepart[j].a[k] = (FLOAT) 0.0;
        for (k=0; k<ndim; k++) activepart[j].agrav[k] = (FLOAT) 0.0;
      }

      // Compute neighbour list for cell depending on physics options
      okflag = tree->ComputeGravityInteractionAndGhostList
        (cell, partdata, macfactor, Nneibmax, Ngravcellmax, Nneib, Nhydroneib,
         Ndirect, Ngravcell, neiblist, mfvlist, directlist, gravcell, neibpart);

      // If there are too many neighbours, reallocate the arrays and recompute the neighbour lists.
      while (okflag < 0 || Nneib > Nneibmax) {
        delete[] neibpartbuf[ithread];
        delete[] cellbuf[ithread];
        delete[] gravlist;
        delete[] directlist;
        delete[] mfvauxlist;
        delete[] mfvlist;
        delete[] neiblist;
        Nneibmax                 = 2*Nneibmax;
        Ngravcellmax             = 2*Ngravcellmax;
        Nneibmaxbuf[ithread]     = Nneibmax;
        Ngravcellmaxbuf[ithread] = Ngravcellmax;
        neiblist                 = new int[Nneibmax];
        mfvlist                  = new int[Nneibmax];
        mfvauxlist               = new int[Nneibmax];
        directlist               = new int[Nneibmax];
        gravlist                 = new int[Nneibmax];
        neibpartbuf[ithread]     = new ParticleType<ndim>[Nneibmax];
        cellbuf[ithread]         = new TreeCell<ndim>[Ngravcellmax];
        neibpart                 = neibpartbuf[ithread];
        gravcell                 = cellbuf[ithread];
        okflag = tree->ComputeGravityInteractionAndGhostList
          (cell, partdata, macfactor, Nneibmax, Ngravcellmax, Nneib, Nhydroneib,
           Ndirect, Ngravcell, neiblist, mfvlist, directlist, gravcell, neibpart);
      };


      // Loop over all active particles in the cell
      //-------------------------------------------------------------------------------------------
      for (j=0; j<Nactive; j++) {
        i = activelist[j];

        // Only calculate gravity for active particle types that have self-gravity activated
        if (mfv->types[activepart[j].ptype].self_gravity){

          Typemask gravmask;
          for (k=0; k < Ntypes; k++) gravmask[k] = mfv->types[activepart[j].ptype].gravmask[k];

          Nhydroaux = 0;
          Ndirectaux = Ndirect;
          for (k=0; k<ndim; k++) rp[k] = activepart[j].r[k];
          hrangesqdi = activepart[j].hrangesqd;

          //---------------------------------------------------------------------------------------
          for (jj=0; jj<Nhydroneib; jj++) {
            int ii = mfvlist[jj];

            // Skip non-gravitating particles and the current active particle.
            if (gravmask[neibpart[jj].ptype] == false) continue;

            // Compute relative position and distance quantities for pair
            for (k=0; k<ndim; k++) draux[k] = neibpart[ii].r[k] - rp[k];
            drsqd = DotProduct(draux, draux, ndim) + small_number;

            if (drsqd <= small_number) continue;

            // Record if neighbour is direct-sum or and SPH neighbour.
            // If SPH neighbour, also record max. timestep level for neighbour
            if (drsqd > hrangesqdi && drsqd >= neibpart[ii].hrangesqd) {
              directlist[Ndirectaux++] = ii;
            }
            else {
              mfvauxlist[Nhydroaux++] = ii;
              levelneib[neiblist[ii]] = max(levelneib[neiblist[ii]], activepart[j].level);
            }
          }

          //---------------------------------------------------------------------------------------


          // Compute forces between SPH neighbours (hydro and gravity)
          mfv->ComputeSmoothedGravForces(i, Nhydroaux, mfvauxlist, activepart[j], neibpart);

          // Compute direct gravity forces between distant particles
          //mfv->ComputeSmoothedGravForces(i, Ndirectaux, directlist, activepart[j], neibpart);
          mfv->ComputeDirectGravForces(i, Ndirectaux, directlist, activepart[j], neibpart);

          // Compute gravitational force due to distant cells
          if (multipole == "monopole") {
          this->ComputeCellMonopoleForces(activepart[j].gpot, activepart[j].agrav,
                                          activepart[j].r, Ngravcell, gravcell);
          }
          else if (multipole == "quadrupole") {
            this->ComputeCellQuadrupoleForces(activepart[j].gpot, activepart[j].agrav,
                                              activepart[j].r, Ngravcell, gravcell);
          }

          // Add the periodic correction force for SPH and direct-sum neighbours
          if (simbox.PeriodicGravity){
            for (jj=0; jj<Nneib; jj++) {
              for (k=0; k<ndim; k++) draux[k] = neibpart[jj].r[k] - activepart[j].r[k];
              ewald->CalculatePeriodicCorrection(neibpart[jj].m, draux, aperiodic, potperiodic);
              for (k=0; k<ndim; k++) activepart[j].agrav[k] += aperiodic[k];
              activepart[j].gpot += potperiodic;
            }

            // Now add the periodic correction force for all cell COMs
            for (jj=0; jj<Ngravcell; jj++) {
              for (k=0; k<ndim; k++) draux[k] = gravcell[jj].r[k] - activepart[j].r[k];
              ewald->CalculatePeriodicCorrection(gravcell[jj].m, draux, aperiodic, potperiodic);
              for (k=0; k<ndim; k++) activepart[j].agrav[k] += aperiodic[k];
              activepart[j].gpot += potperiodic;
            }
          }
        }
      }
      //-------------------------------------------------------------------------------------------


      // Compute 'fast' multipole terms here
      if (multipole == "fast_monopole") {
        this->ComputeFastMonopoleForces(Nactive, Ngravcell, gravcell, cell, activepart);
      }

      // Compute all star forces for active particles
      for (j=0; j<Nactive; j++) {
        if (activelist[j] < mfv->Nhydro) {
          mfv->ComputeStarGravForces(nbody->Nnbody, nbody->nbodydata, activepart[j]);
        }
      }

      // Add all active particles contributions to main array
      for (j=0; j<Nactive; j++) {
        i = activelist[j];
        for (k=0; k<ndim; k++) partdata[i].a[k]     = activepart[j].a[k];
        for (k=0; k<ndim; k++) partdata[i].agrav[k] = activepart[j].agrav[k];
        for (k=0; k<ndim; k++) partdata[i].a[k]    += partdata[i].agrav[k];
        partdata[i].gpot  = activepart[j].gpot;
        levelneib[i]      = max(levelneib[i], activepart[j].levelneib);
      }

    }
    //=============================================================================================


    // Finally, add all contributions from distant pair-wise forces to arrays
#pragma omp critical
    for (i=0; i<mfv->Nhydro; i++) {
      partdata[i].levelneib = max(partdata[i].levelneib, levelneib[i]);
    }

    // Free-up local memory for OpenMP thread
    delete[] directlist;
    delete[] mfvauxlist;
    delete[] mfvlist;
    delete[] neiblist;

  }
  //===============================================================================================

  delete[] celllist;

  timing->EndTimingSection("MFV_GRAV_FORCES");

  return;
}




template class MeshlessFVTree<1,MeshlessFVParticle,KDTreeCell>;
template class MeshlessFVTree<2,MeshlessFVParticle,KDTreeCell>;
template class MeshlessFVTree<3,MeshlessFVParticle,KDTreeCell>;
template class MeshlessFVKDTree<1,MeshlessFVParticle,KDTreeCell>;
template class MeshlessFVKDTree<2,MeshlessFVParticle,KDTreeCell>;
template class MeshlessFVKDTree<3,MeshlessFVParticle,KDTreeCell>;

template class MeshlessFVTree<1,MeshlessFVParticle,OctTreeCell>;
template class MeshlessFVTree<2,MeshlessFVParticle,OctTreeCell>;
template class MeshlessFVTree<3,MeshlessFVParticle,OctTreeCell>;
template class MeshlessFVOctTree<1,MeshlessFVParticle,OctTreeCell>;
template class MeshlessFVOctTree<2,MeshlessFVParticle,OctTreeCell>;
template class MeshlessFVOctTree<3,MeshlessFVParticle,OctTreeCell>;
