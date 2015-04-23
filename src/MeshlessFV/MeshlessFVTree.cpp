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
 NeighbourSearch<ndim>(kernrangeaux, boxaux, kernaux, timingaux),
 MeshlessFVTree<ndim,ParticleType,TreeCell>(Nleafmaxaux,Nmpiaux,thetamaxsqdaux,kernrangeaux,
                                          macerroraux,gravity_mac_aux,multipole_aux,
                                          boxaux,kernaux,timingaux)
{
  // Set-up main tree object
  tree = new KDTree<ndim,ParticleType,TreeCell>(Nleafmaxaux, thetamaxsqdaux, kernrangeaux,
                                                macerroraux, gravity_mac_aux, multipole_aux);

  // Set-up ghost-particle tree object
  ghosttree = new KDTree<ndim,ParticleType,TreeCell>(Nleafmaxaux, thetamaxsqdaux, kernrangeaux,
                                                     macerroraux, gravity_mac_aux, multipole_aux);

#ifdef MPI_PARALLEL
  // Set-up ghost-particle tree object
  mpighosttree = new KDTree<ndim,ParticleType,TreeCell>(Nleafmaxaux, thetamaxsqdaux, kernrangeaux,
                                                       macerroraux, gravity_mac_aux, multipole_aux);

  // Set-up multiple pruned trees, one for each MPI process
  Nghostpruned = new int[Nmpi];
  KDTree<ndim,ParticleType,TreeCell>** prunedtree_derived = new KDTree<ndim,ParticleType,TreeCell>*[Nmpi];
  prunedtree = (Tree<ndim,ParticleType,TreeCell> **) prunedtree_derived;

  for (int i=0; i<Nmpi; i++) {
    prunedtree[i] = new KDTree<ndim,ParticleType,TreeCell>
      (1, thetamaxsqdaux, kernrangeaux, macerroraux, gravity_mac_aux, multipole_aux);
  }

  KDTree<ndim,ParticleType,TreeCell>*** ghostprunedtree_derived = new KDTree<ndim,ParticleType,TreeCell>**[Nmpi];
  ghostprunedtree = (Tree<ndim,ParticleType,TreeCell> ***) ghostprunedtree_derived;
  for (int i=0; i<Nmpi; i++) {
    KDTree<ndim,ParticleType,TreeCell>** ghostprunedtree_derived2 = new KDTree<ndim,ParticleType,TreeCell>*[Nghostprunedmax];
    ghostprunedtree[i] = (Tree<ndim,ParticleType,TreeCell> **) ghostprunedtree_derived2;
    for (int j=0; j<Nghostprunedmax; j++) {
      ghostprunedtree[i][j] = new KDTree<ndim,ParticleType,TreeCell>
        (1, thetamaxsqdaux, kernrangeaux, macerroraux, gravity_mac_aux, multipole_aux);
    }
  }
#endif
}



//=================================================================================================
//  MeshlessFVOctTree::MeshlessFVOctTree
/// SphTree constructor.  Initialises various variables.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
MeshlessFVOctTree<ndim,ParticleType,TreeCell>::MeshlessFVOctTree
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
 NeighbourSearch<ndim>(kernrangeaux, boxaux, kernaux, timingaux),
 MeshlessFVTree<ndim,ParticleType,TreeCell>(Nleafmaxaux, Nmpiaux, thetamaxsqdaux, kernrangeaux,
                                          macerroraux, gravity_mac_aux, multipole_aux,
                                          boxaux, kernaux, timingaux)
{
  // Set-up main tree object
  tree = new OctTree<ndim,ParticleType,TreeCell>(Nleafmaxaux, thetamaxsqdaux, kernrangeaux,
                                                 macerroraux, gravity_mac_aux, multipole_aux);

  // Set-up ghost-particle tree object
  ghosttree = new OctTree<ndim,ParticleType,TreeCell>(Nleafmaxaux, thetamaxsqdaux, kernrangeaux,
                                                      macerroraux, gravity_mac_aux, multipole_aux);

#ifdef MPI_PARALLEL
  // Set-up ghost-particle tree object
  mpighosttree = new OctTree<ndim,ParticleType,TreeCell>(Nleafmaxaux, thetamaxsqdaux, kernrangeaux,
                                                         macerroraux, gravity_mac_aux, multipole_aux);

  // Set-up multiple pruned trees, one for each MPI process
  *(prunedtree) = *(new OctTree<ndim,ParticleType,TreeCell>*[Nmpi]);
  for (int j=0; j<Nmpi; j++) {
    prunedtree[j] = new OctTree<ndim,ParticleType,TreeCell>
      (Nleafmaxaux, thetamaxsqdaux, kernrangeaux, macerroraux, gravity_mac_aux, multipole_aux);
  }
#endif
}



//=================================================================================================
//  MeshlessFVTree::MeshlessFVTree
/// MeshlessFVTree constructor.  Initialises various variables.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
MeshlessFVTree<ndim,ParticleType,TreeCell>::MeshlessFVTree
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
 NeighbourSearch<ndim>(kernrangeaux, boxaux, kernaux, timingaux),
 MeshlessFVNeighbourSearch<ndim>(kernrangeaux, boxaux, kernaux, timingaux),
 HydroTree<ndim,ParticleType,TreeCell>(Nleafmaxaux, Nmpiaux, thetamaxsqdaux, kernrangeaux,
                                            macerroraux, gravity_mac_aux, multipole_aux,
                                            boxaux, kernaux, timingaux)
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
  Nbody<ndim> *nbody)                      ///< [in] Pointer to N-body object
{
  int cactive;                             // No. of active tree cells
  TreeCell<ndim> *celllist;                // List of active tree cells
  //ParticleType<ndim> *sphdata = static_cast<ParticleType<ndim>* > (sph_gen);
#ifdef MPI_PARALLEL
  int Nactivetot = 0;                      // Total number of active particles
  double twork = timing->WallClockTime();  // Start time (for load balancing)
#endif

  debug2("[MeshlessFVTree::UpdateAllProperties]");
  timing->StartTimingSection("MFV_PROPERTIES");


  // Find list of all cells that contain active particles
  celllist = new TreeCell<ndim>[tree->gtot];
  cactive = tree->ComputeActiveCellList(celllist);


  // Set-up all OMP threads
  //===============================================================================================
#pragma omp parallel default(none) shared(cactive,celllist,cout,nbody,mfv,mfvdata)
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
    FLOAT* gpot   = new FLOAT[Nneibmax];       // Local array of particle potentials
    FLOAT* gpot2  = new FLOAT[Nneibmax];       // Local reduced array of neighbour potentials
    FLOAT* drsqd  = new FLOAT[Nneibmax];       // Local array of distances (squared)
    FLOAT* m      = new FLOAT[Nneibmax];       // Local array of particle masses
    FLOAT* m2     = new FLOAT[Nneibmax];       // Local reduced array of neighbour masses
    FLOAT* r      = new FLOAT[Nneibmax*ndim];                  // Local array of particle positions
    ParticleType<ndim>* activepart = activepartbuf[ithread];   // Local array of active particles
    //ParticleType<ndim>* activepart = new ParticleType<ndim>[Nleafmax];


    // Loop over all active cells
    //=============================================================================================
#pragma omp for schedule(guided)
    for (cc=0; cc<cactive; cc++) {
      TreeCell<ndim>& cell = celllist[cc];
      celldone = 1;
      hmax = cell.hmax;

      // Sanity checks
      //assert(cell.Nactive > 0);

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
          delete[] neiblist;
          Nneibmax = 2*Nneibmax;
          neiblist = new int[Nneibmax];
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
          j        = neiblist[jj];
          gpot[jj] = mfvdata[j].gpot;
          m[jj]    = mfvdata[j].m;
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
          if (neibcheck) this->CheckValidNeighbourList(i,Ntot,Nneib,neiblist,mfvdata,"gather");
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
  Nbody<ndim> *nbody)                      ///< [in] Pointer to N-body object
{
  int cactive;                             // No. of active cells
  TreeCell<ndim> *celllist;                // List of active tree cells
  //ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sph_gen);
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
  if (cactive == 0) return;

  // Update ghost tree smoothing length values here
  //tree->UpdateHmaxValues(tree->celldata[0],mfvdata);
  if (ghosttree->Ntot > 0) ghosttree->UpdateHmaxValues(ghosttree->celldata[0], mfvdata);


  // Set-up all OMP threads
  //===============================================================================================
#pragma omp parallel default(none) shared(cactive,celllist,nbody,mfv,mfvdata)
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
    int Nneibmax      = Nneibmaxbuf[ithread];      // ..
    int* activelist   = activelistbuf[ithread];    // ..
    int* levelneib    = levelneibbuf[ithread];     // ..
    int* neiblist     = new int[Nneibmax];         // ..
    int* mfvlist      = new int[Nneibmax];         // ..
    FLOAT* dr         = new FLOAT[Nneibmax*ndim];  // ..
    FLOAT* drmag      = new FLOAT[Nneibmax];       // ..
    FLOAT* invdrmag   = new FLOAT[Nneibmax];       //..
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
      for (j=0; j<Nactive; j++) {
        activepart[j] = mfvdata[activelist[j]];
      }

      // Compute neighbour list for cell from real and periodic ghost particles
      Nneib = 0;
      Nneib = tree->ComputeNeighbourList(cell, mfvdata, Nneibmax, Nneib, neiblist, neibpart);
      Nneib = ghosttree->ComputeNeighbourList(cell, mfvdata, Nneibmax, Nneib, neiblist, neibpart);


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
        Nneib = tree->ComputeNeighbourList(cell, mfvdata, Nneibmax, Nneib, neiblist, neibpart);
        Nneib = ghosttree->ComputeNeighbourList(cell, mfvdata, Nneibmax, Nneib, neiblist, neibpart);
      };


      // Loop over all active particles in the cell
      //-------------------------------------------------------------------------------------------
      for (j=0; j<Nactive; j++) {
        i = activelist[j];

        for (k=0; k<ndim; k++) rp[k] = activepart[j].r[k];
        hrangesqdi = activepart[j].hrangesqd;
        Nhydroaux = 0;

        // Validate that gather neighbour list is correct
#if defined(VERIFY_ALL)
        if (neibcheck) this->CheckValidNeighbourList(i,Ntot,Nneib,neiblist,mfvdata,"all");
#endif

        // Compute distances and the inverse between the current particle and all neighbours here,
        // for both gather and inactive scatter neibs.  Only consider particles with j > i to
        // compute pair forces once unless particle j is inactive.
        //-----------------------------------------------------------------------------------------
        for (jj=0; jj<Nneib; jj++) {

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
 (const int Nhydro,                    ///< [in] No. of SPH particles
  const int Ntot,                      ///< [in] No. of SPH + ghost particles
  const FLOAT timestep,                ///< [in] ..
  MeshlessFVParticle<ndim> *mfvdata,   ///< [inout] Pointer to SPH ptcl array
  MeshlessFV<ndim> *mfv,               ///< [in] Pointer to SPH object
  Nbody<ndim> *nbody)                  ///< [in] Pointer to N-body object
{
  int cactive;                             // No. of active cells
  TreeCell<ndim> *celllist;                // List of active tree cells
  //ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sph_gen);
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
  if (cactive == 0) return;

  // Update ghost tree smoothing length values here
  //tree->UpdateHmaxValues(tree->celldata[0],mfvdata);
  if (ghosttree->Ntot > 0) ghosttree->UpdateHmaxValues(ghosttree->celldata[0],mfvdata);


  // Set-up all OMP threads
  //===============================================================================================
#pragma omp parallel default(none) shared(cactive,celllist,nbody,mfv,mfvdata)
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
    int Nhydroaux;                                   // ..
    FLOAT draux[ndim];                             // Aux. relative position vector
    FLOAT drsqd;                                   // Distance squared
    FLOAT hrangesqdi;                              // Kernel gather extent
    FLOAT rp[ndim];                                // Local copy of particle position
    int Nneibmax      = Nneibmaxbuf[ithread];      // ..
    int* activelist   = activelistbuf[ithread];    // ..
    int* levelneib    = levelneibbuf[ithread];     // ..
    int* neiblist     = new int[Nneibmax];         // ..
    int* mfvlist      = new int[Nneibmax];         // ..
    FLOAT* dr         = new FLOAT[Nneibmax*ndim];  // ..
    FLOAT* drmag      = new FLOAT[Nneibmax];       // ..
    FLOAT* invdrmag   = new FLOAT[Nneibmax];       // ..
    FLOAT (*fluxBuffer)[ndim+2] = new FLOAT[Ntot][ndim+2];     // ..
    ParticleType<ndim>* activepart = activepartbuf[ithread];   // ..
    ParticleType<ndim>* neibpart   = neibpartbuf[ithread];     // ..

    for (i=0; i<mfv->Nhydro; i++) levelneib[i] = 0;
    for (i=0; i<Ntot; i++) {
      for (k=0; k<ndim+2; k++) fluxBuffer[i][k] = (FLOAT) 0.0;
    }


    // Loop over all active cells
    //=============================================================================================
#pragma omp for schedule(guided)
    for (cc=0; cc<cactive; cc++) {
      TreeCell<ndim>& cell = celllist[cc];

      // Find list of active particles in current cell
      Nactive = tree->ComputeActiveParticleList(cell,mfvdata,activelist);

      // Make local copies of active particles
      for (j=0; j<Nactive; j++) activepart[j] = mfvdata[activelist[j]];


      // Compute neighbour list for cell from real and periodic ghost particles
      Nneib = 0;
      Nneib = tree->ComputeNeighbourList(cell, mfvdata, Nneibmax, Nneib, neiblist, neibpart);
      Nneib = ghosttree->ComputeNeighbourList(cell, mfvdata, Nneibmax, Nneib, neiblist, neibpart);

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
        Nneib = tree->ComputeNeighbourList(cell, mfvdata, Nneibmax, Nneib, neiblist, neibpart);
        Nneib = ghosttree->ComputeNeighbourList(cell, mfvdata, Nneibmax, Nneib, neiblist, neibpart);
      };

      for (j=0; j<Nneibmax; j++) {
        for (k=0; k<ndim+2; k++) neibpart[j].dQdt[k] = (FLOAT) 0.0;
      }


      // Loop over all active particles in the cell
      //-------------------------------------------------------------------------------------------
      for (j=0; j<Nactive; j++) {
        i = activelist[j];

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

          // Skip current active particle
          if (neiblist[jj] == i || (neiblist[jj] > i && neibpart[jj].active)) continue;

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
        mfv->ComputeGodunovFlux(i, Nhydroaux, timestep, mfvlist, drmag, invdrmag, dr,
                                activepart[j], neibpart);

      }
      //-------------------------------------------------------------------------------------------


      // Accumlate fluxes for neighbours
      for (int jj=0; jj<Nneib; jj++) {
        j = neiblist[jj];
        for (k=0; k<ndim+2; k++) fluxBuffer[j][k] += neibpart[jj].dQdt[k];
      }

      // Add all active particles contributions to main array
      for (j=0; j<Nactive; j++) {
        i = activelist[j];
        for (k=0; k<ndim+2; k++) mfvdata[i].dQdt[k] = activepart[j].dQdt[k];
      }

    }
    //=============================================================================================


    // Add all buffers back to main arrays
#pragma omp barrier
#pragma omp critical
    {
      for (i=0; i<Nhydro; i++) {
        for (k=0; k<ndim+2; k++) mfvdata[i].dQdt[k] += fluxBuffer[i][k];
      }
    }


    // Free-up local memory for OpenMP thread
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
