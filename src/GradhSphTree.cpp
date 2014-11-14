//=================================================================================================
//  GradhSphTree.cpp
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
#include <string>
#include <math.h>
#include "Precision.h"
#include "Exception.h"
#include "SphNeighbourSearch.h"
#include "Sph.h"
#include "Parameters.h"
#include "InlineFuncs.h"
#include "SphParticle.h"
#include "Debug.h"
#if defined _OPENMP
#include <omp.h>
#endif
using namespace std;



//=================================================================================================
//  GradhSphKDTree::GradhSphKDTree
/// GradhSphKDTree constructor.  Initialises various variables and creates tree objects.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
GradhSphKDTree<ndim,ParticleType,TreeCell>::GradhSphKDTree
 (int Nleafmaxaux,
  int Nmpiaux,
  FLOAT thetamaxsqdaux,
  FLOAT kernrangeaux,
  FLOAT macerroraux,
  string gravity_mac_aux,
  string multipole_aux,
  DomainBox<ndim> *boxaux,
  SphKernel<ndim> *kernaux,
  CodeTiming *timingaux):
 GradhSphTree<ndim,ParticleType,TreeCell>(Nleafmaxaux,Nmpiaux,thetamaxsqdaux,kernrangeaux,
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
  KDTree<ndim,ParticleType,TreeCell>** prunedtree_derived = new KDTree<ndim,ParticleType,TreeCell>*[Nmpi];
  prunedtree = (Tree<ndim,ParticleType,TreeCell> **) prunedtree_derived;
  for (int j=0; j<Nmpi; j++) {
    prunedtree[j] = new KDTree<ndim,ParticleType,TreeCell>
      (1, thetamaxsqdaux, kernrangeaux, macerroraux, gravity_mac_aux, multipole_aux);
  }
#endif
}



//=================================================================================================
//  GradhSphOctTree::GradhSphOctTree
/// SphTree constructor.  Initialises various variables.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
GradhSphOctTree<ndim,ParticleType,TreeCell>::GradhSphOctTree
 (int Nleafmaxaux,
  int Nmpiaux,
  FLOAT thetamaxsqdaux,
  FLOAT kernrangeaux,
  FLOAT macerroraux,
  string gravity_mac_aux,
  string multipole_aux,
  DomainBox<ndim> *boxaux,
  SphKernel<ndim> *kernaux,
  CodeTiming *timingaux):
 GradhSphTree<ndim,ParticleType,TreeCell>(Nleafmaxaux,Nmpiaux,thetamaxsqdaux,kernrangeaux,
                                          macerroraux,gravity_mac_aux,multipole_aux,
                                          boxaux,kernaux,timingaux)
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
//  GradhSphTree::GradhSphTree
/// GradhSphTree constructor.  Initialises various variables.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
GradhSphTree<ndim,ParticleType,TreeCell>::GradhSphTree
 (int Nleafmaxaux,
  int Nmpiaux,
  FLOAT thetamaxsqdaux,
  FLOAT kernrangeaux,
  FLOAT macerroraux,
  string gravity_mac_aux,
  string multipole_aux,
  DomainBox<ndim> *boxaux,
  SphKernel<ndim> *kernaux,
  CodeTiming *timingaux):
 SphTree<ndim,ParticleType,TreeCell>(Nleafmaxaux,Nmpiaux,thetamaxsqdaux,kernrangeaux,
                                     macerroraux,gravity_mac_aux,multipole_aux,
                                     boxaux,kernaux,timingaux)
{
}



//=================================================================================================
//  GradhSphTree::~GradhSphTree
/// GradhSphTree destructor.  Deallocates tree memory upon object destruction.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
GradhSphTree<ndim,ParticleType,TreeCell>::~GradhSphTree()
{
  if (tree->allocated_tree) {
    this->DeallocateMemory();
    tree->DeallocateTreeMemory();
  }
}



//=================================================================================================
//  GradhSphTree::UpdateAllSphProperties
/// Update all gather SPH properties (e.g. rho, div_v) for all active particles in domain.
/// Loops over all cells containing active particles, performs a tree walk for all particles in
/// the cell, and then calls SPH class routine to compute properties from neighbours.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void GradhSphTree<ndim,ParticleType,TreeCell>::UpdateAllSphProperties
 (int Nsph,                                ///< [in] No. of SPH particles
  int Ntot,                                ///< [in] No. of SPH + ghost particles
  SphParticle<ndim> *sph_gen,              ///< [inout] Pointer to SPH ptcl array
  Sph<ndim> *sph,                          ///< [in] Pointer to SPH object
  Nbody<ndim> *nbody)                      ///< [in] Pointer to N-body object
{
  int celldone;                            // Flag if cell is done
  int okflag;                              // Flag if particle is done
  int cc;                                  // Aux. cell counter
  int cactive;                             // No. of active
  int i;                                   // Particle id
  int ithread;                             // OpenMP thread i.d.
  int j;                                   // Aux. particle counter
  int jj;                                  // Aux. particle counter
  int k;                                   // Dimension counter
  int Ngather;                             // No. of near gather neighbours
  int Nneib;                               // No. of neighbours
  FLOAT draux[ndim];                       // Aux. relative position vector var
  FLOAT drsqdaux;                          // Distance squared
  FLOAT hrangesqd;                         // Kernel extent
  FLOAT hmax;                              // Maximum smoothing length
  FLOAT rp[ndim];                          // Local copy of particle position
  //TreeCell<ndim> *cell;                    // Pointer to binary tree cell
  TreeCell<ndim> **celllist;               // List of binary cell pointers
  ParticleType<ndim> *sphdata = static_cast<ParticleType<ndim>* > (sph_gen);
#ifdef MPI_PARALLEL
  int Nactivetot = 0;                      // Total number of active particles
  double twork = timing->WallClockTime();  // Start time (for load balancing)
#endif

  debug2("[GradhSphTree::UpdateAllSphProperties]");
  timing->StartTimingSection("SPH_PROPERTIES",2);


  // Find list of all cells that contain active particles
  celllist = new TreeCell<ndim>*[tree->gtot];
  cactive = tree->ComputeActiveCellList(celllist);


  // Set-up all OMP threads
  //===============================================================================================
#pragma omp parallel default(none) shared(cactive,celllist,cout,nbody,sph,sphdata)\
  private(cc,celldone,draux,drsqdaux,hmax,hrangesqd,i,j,jj,k,okflag,Ngather,Nneib,rp)
  {
#if defined _OPENMP
    int ithread = omp_get_thread_num();
#else
    int ithread = 0;
#endif
    int Nactive;
    int Nneibmax = Nneibmaxbuf[ithread];
    int* activelist = new int[Nleafmax];
    //int* activelist = activelistbuf[ithread];
    int* neiblist = new int[Nneibmax];
    FLOAT *mu, *mu2;
    FLOAT* gpot   = new FLOAT[Nneibmax];
    FLOAT* gpot2  = new FLOAT[Nneibmax];
    FLOAT* drsqd  = new FLOAT[Nneibmax];
    FLOAT* m      = new FLOAT[Nneibmax];
    FLOAT* m2     = new FLOAT[Nneibmax];
    FLOAT* r      = new FLOAT[Nneibmax*ndim];
    ParticleType<ndim>* activepart = activepartbuf[ithread];
    //ParticleType<ndim>* activepart = new ParticleType<ndim>[Nleafmax];


    // Loop over all active cells
    //=============================================================================================
#pragma omp for schedule(guided)
    for (cc=0; cc<cactive; cc++) {
      TreeCell<ndim> cell = *celllist[cc];
      celldone = 1;
      hmax = cell.hmax;

      // If hmax is too small so the neighbour lists are invalid, make hmax
      // larger and then recompute for the current active cell.
      //-------------------------------------------------------------------------------------------
      do {
        hmax = 1.05*hmax;
        celldone = 1;

        // Find list of active particles in current cell
        Nactive = tree->ComputeActiveParticleList(cell,sphdata,activelist);
        for (j=0; j<Nactive; j++) activepart[j] = sphdata[activelist[j]];

        // Compute neighbour list for cell from particles on all trees
        Nneib = 0;
        Nneib = tree->ComputeGatherNeighbourList(cell,sphdata,hmax,Nneibmax,Nneib,neiblist);
        Nneib = ghosttree->ComputeGatherNeighbourList(cell,sphdata,hmax,Nneibmax,Nneib,neiblist);
#ifdef MPI_PARALLEL
        Nneib = mpighosttree->ComputeGatherNeighbourList(cell,sphdata,hmax,Nneibmax,Nneib,neiblist);
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
          Nneib = tree->ComputeGatherNeighbourList(cell,sphdata,hmax,Nneibmax,Nneib,neiblist);
          Nneib = ghosttree->ComputeGatherNeighbourList(cell,sphdata,hmax,Nneibmax,Nneib,neiblist);
#ifdef MPI_PARALLEL
          Nneib = mpighosttree->ComputeGatherNeighbourList(cell,sphdata,hmax,
                                                           Nneibmax,Nneib,neiblist);
#endif
        };


        // Make local copies of important neib information (mass and position)
        for (jj=0; jj<Nneib; jj++) {
          j        = neiblist[jj];
          gpot[jj] = sphdata[j].gpot;
          m[jj]    = sphdata[j].m;
          for (k=0; k<ndim; k++) r[ndim*jj + k] = (FLOAT) sphdata[j].r[k];
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
            drsqdaux = DotProduct(draux,draux,ndim);

            // Record distance squared for all potential gather neighbours
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
          if (neibcheck) this->CheckValidNeighbourList(i,Ntot,Nneib,neiblist,sphdata,"gather");
#endif

          // Compute smoothing length and other gather properties for ptcl i
          okflag = sph->ComputeH(i,Ngather,hmax,m2,mu2,drsqd,gpot,activepart[j],nbody);

          // If h-computation is invalid, then break from loop and recompute
          // larger neighbour lists
          if (okflag == 0) {
            celldone = 0;
            break;
          }

        }
        //-----------------------------------------------------------------------------------------

      } while (celldone == 0);
      //-------------------------------------------------------------------------------------------

      // Once cell is finished, copy all active particles back to main memory
      for (j=0; j<Nactive; j++) sphdata[activelist[j]] = activepart[j];

#ifdef MPI_PARALLEL
      Nactivetot += Nactive;
#endif

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
    //delete[] activepart;
    delete[] activelist;

  }
  //===============================================================================================

  // Compute time spent in routine and in each cell for load balancing
#ifdef MPI_PARALLEL
  twork = timing->WallClockTime() - twork;
  for (cc=0; cc<cactive; cc++) {
    celllist[cc].worktot += twork*(double) celllist[cc].Nactive / (double) Nactivetot;
  }
#endif

  delete[] celllist;

  // Update tree smoothing length values here
  tree->UpdateHmaxValues(tree->celldata[0],sphdata);

  timing->EndTimingSection("SPH_PROPERTIES");

  return;
}



//=================================================================================================
//  GradhSphTree::UpdateAllSphHydroForces
/// Compute hydro forces for all active SPH particles.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void GradhSphTree<ndim,ParticleType,TreeCell>::UpdateAllSphHydroForces
 (int Nsph,                            ///< [in] No. of SPH particles
  int Ntot,                            ///< [in] No. of SPH + ghost particles
  SphParticle<ndim> *sph_gen,          ///< [inout] Pointer to SPH ptcl array
  Sph<ndim> *sph,                      ///< [in] Pointer to SPH object
  Nbody<ndim> *nbody)                  ///< [in] Pointer to N-body object
{
  int cactive;                         // No. of active cells
  int cc;                              // Aux. cell counter
  int i;                               // Particle id
  int j;                               // Aux. particle counter
  int jj;                              // Aux. particle counter
  int k;                               // Dimension counter
  FLOAT draux[ndim];                   // Aux. relative position vector
  FLOAT drsqd;                         // Distance squared
  FLOAT hrangesqdi;                    // Kernel gather extent
  FLOAT rp[ndim];                      // Local copy of particle position
  TreeCell<ndim> **celllist;           // List of binary tree pointers
  ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sph_gen);

  debug2("[GradhSphTree::UpdateAllSphHydroForces]");
  timing->StartTimingSection("SPH_HYDRO_FORCES",2);


  // Find list of all cells that contain active particles
#if defined (MPI_PARALLEL)
  celllist = new TreeCell<ndim>*[tree->Ncellmax];
#else
  celllist = new TreeCell<ndim>*[tree->gtot];
#endif
  cactive = tree->ComputeActiveCellList(celllist);

  // If there are no active cells, return to main loop
  if (cactive == 0) return;

  // Update ghost tree smoothing length values here
  //tree->UpdateHmaxValues(tree->celldata[0],sphdata);
  if (ghosttree->Ntot > 0) ghosttree->UpdateHmaxValues(ghosttree->celldata[0],sphdata);


  // Set-up all OMP threads
  //===============================================================================================
#pragma omp parallel default(none) shared(cactive,celllist,nbody,sph,sphdata)\
  private(cc,draux,drsqd,hrangesqdi,i,j,jj,k,rp)
  {
#if defined _OPENMP
    int ithread = omp_get_thread_num();
#else
    int ithread = 0;
#endif
    int Nactive;
    int Nneib;
    int Nsphaux;
    int Nneibmax      = Nneibmaxbuf[ithread];
    int* activelist   = activelistbuf[ithread];
    int* levelneib    = levelneibbuf[ithread];
    int* neiblist     = new int[Nneibmax];
    int* sphlist      = new int[Nneibmax];
    FLOAT* dr         = new FLOAT[Nneibmax*ndim];
    FLOAT* drmag      = new FLOAT[Nneibmax];
    FLOAT* invdrmag   = new FLOAT[Nneibmax];
    ParticleType<ndim>* activepart = activepartbuf[ithread];
    ParticleType<ndim>* neibpart   = neibpartbuf[ithread];
    //TreeCell<ndim> *cell;

    for (i=0; i<sph->Nsph; i++) levelneib[i] = 0;


    // Loop over all active cells
    //=============================================================================================
#pragma omp for schedule(guided)
    for (cc=0; cc<cactive; cc++) {
      TreeCell<ndim> cell = *celllist[cc];

      // Find list of active particles in current cell
      Nactive = tree->ComputeActiveParticleList(cell,sphdata,activelist);

      // Make local copies of active particles
      for (j=0; j<Nactive; j++) {
        activepart[j] = sphdata[activelist[j]];
        activepart[j].div_v     = (FLOAT) 0.0;
        activepart[j].dudt      = (FLOAT) 0.0;
        activepart[j].dalphadt  = (FLOAT) 0.0;
        activepart[j].gpot      = (FLOAT) 0.0;
        activepart[j].levelneib = 0;
        for (k=0; k<ndim; k++) activepart[j].a[k] = (FLOAT) 0.0;
        for (k=0; k<ndim; k++) activepart[j].agrav[k] = (FLOAT) 0.0;
      }


      // Compute neighbour list for cell from real and periodic ghost particles
      Nneib = 0;
      Nneib = tree->ComputeNeighbourList(cell,sphdata,Nneibmax,Nneib,neiblist,neibpart);
      Nneib = ghosttree->ComputeNeighbourList(cell,sphdata,Nneibmax,Nneib,neiblist,neibpart);

      // If there are too many neighbours, reallocate the arrays and
      // recompute the neighbour list.
      while (Nneib == -1) {
        delete[] neibpartbuf[ithread];
        delete[] invdrmag;
        delete[] drmag;
        delete[] dr;
        delete[] sphlist;
        delete[] neiblist;
        Nneibmax                  = 2*Nneibmax;
        Nneibmaxbuf[ithread]      = Nneibmax;
        Ngravcellmaxbuf[ithread] *= 2;
        neiblist                  = new int[Nneibmax];
        sphlist                   = new int[Nneibmax];
        dr                        = new FLOAT[Nneibmax*ndim];
        drmag                     = new FLOAT[Nneibmax];
        invdrmag                  = new FLOAT[Nneibmax];
        neibpartbuf[ithread]      = new ParticleType<ndim>[Nneibmax];
        neibpart                  = neibpartbuf[ithread];
        Nneib = 0;
        Nneib = tree->ComputeNeighbourList(cell,sphdata,Nneibmax,Nneib,neiblist,neibpart);
        Nneib = ghosttree->ComputeNeighbourList(cell,sphdata,Nneibmax,Nneib,neiblist,neibpart);
      };

      //for (j=0; j<Nneib; j++) neibpart[j] = sphdata[neiblist[j]];


      // Loop over all active particles in the cell
      //-------------------------------------------------------------------------------------------
      for (j=0; j<Nactive; j++) {
        i = activelist[j];

        for (k=0; k<ndim; k++) rp[k] = activepart[j].r[k];
        hrangesqdi = activepart[j].hrangesqd;
        Nsphaux = 0;

        // Validate that gather neighbour list is correct
#if defined(VERIFY_ALL)
        if (neibcheck) this->CheckValidNeighbourList(i,Ntot,Nneib,neiblist,sphdata,"all");
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
            drmag[Nsphaux] = sqrt(drsqd);
            invdrmag[Nsphaux] = (FLOAT) 1.0/drmag[Nsphaux];
            for (k=0; k<ndim; k++) dr[Nsphaux*ndim + k] = draux[k]*invdrmag[Nsphaux];
            levelneib[neiblist[jj]] = max(levelneib[neiblist[jj]],activepart[j].level);
            sphlist[Nsphaux] = jj;
            Nsphaux++;
          }

        }
        //-----------------------------------------------------------------------------------------

        // Compute all neighbour contributions to hydro forces
        sph->ComputeSphHydroForces(i,Nsphaux,sphlist,drmag,invdrmag,dr,activepart[j],neibpart);

      }
      //-------------------------------------------------------------------------------------------


      // Compute all star forces for active particles
      if (nbody->Nnbody > 0) {
        for (j=0; j<Nactive; j++)
          sph->ComputeStarGravForces(nbody->Nnbody,nbody->nbodydata,activepart[j]);
      }


      // Add all active particles contributions to main array
      for (j=0; j<Nactive; j++) {
        i = activelist[j];
        for (k=0; k<ndim; k++) sphdata[i].a[k]     = activepart[j].a[k];
        for (k=0; k<ndim; k++) sphdata[i].agrav[k] = activepart[j].agrav[k];
        for (k=0; k<ndim; k++) sphdata[i].a[k]     += sphdata[i].agrav[k];
        sphdata[i].gpot     = activepart[j].gpot;
        sphdata[i].dudt     = activepart[j].dudt;
        sphdata[i].dalphadt = activepart[j].dalphadt;
        sphdata[i].div_v    = activepart[j].div_v;
        sphdata[i].active   = false;
        levelneib[i]        = max(levelneib[i],activepart[j].levelneib);
      }

    }
    //=============================================================================================


    // Finally, add all contributions from distant pair-wise forces to arrays
#pragma omp critical
    for (i=0; i<sph->Nsph; i++) sphdata[i].levelneib = max(sphdata[i].levelneib,levelneib[i]);


    // Free-up local memory for OpenMP thread
    delete[] invdrmag;
    delete[] drmag;
    delete[] dr;
    delete[] sphlist;
    delete[] neiblist;

  }
  //===============================================================================================

  delete[] celllist;

  timing->EndTimingSection("SPH_HYDRO_FORCES");

  return;
}



//=================================================================================================
//  GradhSphTree::UpdateAllSphForces
/// Compute all forces on active SPH particles (hydro + gravity)
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void GradhSphTree<ndim,ParticleType,TreeCell>::UpdateAllSphForces
 (int Nsph,                            ///< [in] No. of SPH particles
  int Ntot,                            ///< [in] No. of SPH + ghost particles
  SphParticle<ndim> *sph_gen,          ///< [inout] Pointer to SPH ptcl array
  Sph<ndim> *sph,                      ///< [in] Pointer to SPH object
  Nbody<ndim> *nbody)                  ///< [in] Pointer to N-body object
{
  int cactive;                         // No. of active cells
  int cc;                              // Aux. cell counter
  int i;                               // Particle id
  int j;                               // Aux. particle counter
  int jj;                              // Aux. particle counter
  int k;                               // Dimension counter
  int okflag;                          // Flag if h-rho iteration is valid
  int Ngravcell;                       // No. of gravity cells
  int Nneib;                           // No. of neighbours
  FLOAT macfactor;                     // Gravity MAC factor
  FLOAT draux[ndim];                   // Aux. relative position vector
  FLOAT drsqd;                         // Distance squared
  FLOAT hrangesqdi;                    // Kernel gather extent
  FLOAT rp[ndim];                      // ..
  //TreeCell<ndim> *cell;                // Pointer to kd-tree cell
  TreeCell<ndim> **celllist;           // List of pointers to tree cells
  ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sph_gen);

  debug2("[GradhSphTree::UpdateAllSphForces]");
  timing->StartTimingSection("SPH_ALL_FORCES",2);

  // Update ghost tree smoothing length values here
  if (ghosttree->Ntot > 0) ghosttree->UpdateHmaxValues(ghosttree->celldata[0],sphdata);

  // Find list of all cells that contain active particles
#if defined (MPI_PARALLEL)
  celllist = new TreeCell<ndim>*[tree->Ncellmax];
#else
  celllist = new TreeCell<ndim>*[tree->gtot];
#endif
  cactive = tree->ComputeActiveCellList(celllist);


  // Set-up all OMP threads
  //===============================================================================================
#pragma omp parallel default(none) shared(celllist,cactive,nbody,sph,sphdata,cout)\
  private(cc,draux,drsqd,hrangesqdi,i,j,jj,k,macfactor,Ngravcell,Nneib,okflag,rp)
  {
#if defined _OPENMP
    int ithread = omp_get_thread_num();
#else
    int ithread = 0;
#endif
    int Nactive;
    int Ndirect;
    int Ndirectaux;
    int Nsphaux;
    int Nsphneib;
    int Nneibmax                   = Nneibmaxbuf[ithread];
    int Ngravcellmax               = Ngravcellmaxbuf[ithread];
    int *levelneib                 = levelneibbuf[ithread];
    int *activelist                = activelistbuf[ithread];
    int *neiblist                  = new int[Nneibmax];
    int *sphlist                   = new int[Nneibmax];
    int *sphauxlist                = new int[Nneibmax];
    int *directlist                = new int[Nneibmax];
    ParticleType<ndim>* activepart = activepartbuf[ithread];
    ParticleType<ndim>* neibpart   = neibpartbuf[ithread];
    TreeCell<ndim>* gravcell       = cellbuf[ithread];

    // Zero timestep level array
    for (i=0; i<sph->Nsph; i++) levelneib[i] = 0.0;


    // Loop over all active cells
    //=============================================================================================
#pragma omp for schedule(guided)
    for (cc=0; cc<cactive; cc++) {
      TreeCell<ndim> cell = *celllist[cc];
      macfactor = 0.0;

      // Find list of active particles in current cell
      Nactive = tree->ComputeActiveParticleList(cell,sphdata,activelist);

      // Make local copies of active particles
      for (j=0; j<Nactive; j++) activepart[j] = sphdata[activelist[j]];

      // Compute average/maximum term for computing gravity MAC
      if (gravity_mac == "eigenmac") {
        for (j=0; j<Nactive; j++) macfactor = max(macfactor,pow(1.0/activepart[j].gpot,twothirds));
      }

      // Zero/initialise all summation variables for active particles
      for (j=0; j<Nactive; j++) {
        activepart[j].div_v     = (FLOAT) 0.0;
        activepart[j].dudt      = (FLOAT) 0.0;
        activepart[j].levelneib = 0;
        activepart[j].gpot      = activepart[j].m*activepart[j].invh*sph->kernp->wpot(0.0);
        for (k=0; k<ndim; k++) activepart[j].a[k]     = (FLOAT) 0.0;
        for (k=0; k<ndim; k++) activepart[j].agrav[k] = (FLOAT) 0.0;
        for (k=0; k<ndim; k++) activepart[j].adot[k]  = (FLOAT) 0.0;
      }

      // Compute neighbour list for cell depending on physics options
      okflag = tree->ComputeGravityInteractionList(cell,sphdata,macfactor,Nneibmax,Ngravcellmax,
                                                   Nneib,Nsphneib,Ndirect,Ngravcell,neiblist,
                                                   sphlist,directlist,gravcell,neibpart);

      // If there are too many neighbours, reallocate the arrays and recompute the neighbour lists.
      while (okflag < 0 || Nneib > Nneibmax) {
        delete[] neibpart;
        delete[] gravcell;
        delete[] directlist;
        delete[] sphauxlist;
        delete[] sphlist;
        delete[] neiblist;
        Nneibmax                 = 2*Nneibmax;
        Ngravcellmax             = 2*Ngravcellmax;
        Nneibmaxbuf[ithread]     = Nneibmax;
        Ngravcellmaxbuf[ithread] = Ngravcellmax;
        neiblist                 = new int[Nneibmax];
        sphlist                  = new int[Nneibmax];
        sphauxlist               = new int[Nneibmax];
        directlist               = new int[Nneibmax];
        neibpartbuf[ithread]     = new ParticleType<ndim>[Nneibmax];
        cellbuf[ithread]         = new TreeCell<ndim>[Ngravcellmax];
        neibpart                 = neibpartbuf[ithread];
        gravcell                 = cellbuf[ithread];
        okflag = tree->ComputeGravityInteractionList(cell,sphdata,macfactor,Nneibmax,
                                                     Ngravcellmax,Nneib,Nsphneib,Ndirect,Ngravcell,
                                                     neiblist,sphlist,directlist,gravcell,neibpart);
      };


      // Loop over all active particles in the cell
      //-------------------------------------------------------------------------------------------
      for (j=0; j<Nactive; j++) {
        i = activelist[j];

        Nsphaux = 0;
        Ndirectaux = Ndirect;
        for (k=0; k<ndim; k++) rp[k] = activepart[j].r[k];
        hrangesqdi = activepart[j].hrangesqd;

        //-----------------------------------------------------------------------------------------
        for (jj=0; jj<Nsphneib; jj++) {
          int ii = sphlist[jj];

          // Compute relative position and distance quantities for pair
          for (k=0; k<ndim; k++) draux[k] = neibpart[ii].r[k] - rp[k];
          drsqd = DotProduct(draux,draux,ndim) + small_number;

          // Record if neighbour is direct-sum or and SPH neighbour.
          // If SPH neighbour, also record max. timestep level for neighbour
          if (drsqd > hrangesqdi && drsqd >= neibpart[ii].hrangesqd) {
            directlist[Ndirectaux++] = ii;
          }
          else if (neiblist[ii] != i) {
            sphauxlist[Nsphaux++] = ii;
            levelneib[neiblist[ii]] = max(levelneib[neiblist[ii]],activepart[j].level);
          }
        }
        //-----------------------------------------------------------------------------------------


        // Compute forces between SPH neighbours (hydro and gravity)
        sph->ComputeSphHydroGravForces(i,Nsphaux,sphauxlist,activepart[j],neibpart);

        // Compute direct gravity forces between distant particles
        sph->ComputeDirectGravForces(i,Ndirectaux,directlist,activepart[j],neibpart);

        // Compute gravitational force due to distant cells
        if (multipole == "monopole") {
          tree->ComputeCellMonopoleForces(activepart[j].gpot,activepart[j].agrav,
                                          activepart[j].r,Ngravcell,gravcell);
        }
        else if (multipole == "quadrupole") {
          tree->ComputeCellQuadrupoleForces(activepart[j].gpot,activepart[j].agrav,
                                            activepart[j].r,Ngravcell,gravcell);
        }

      }
      //-------------------------------------------------------------------------------------------


      // Compute 'fast' multipole terms here
      if (multipole == "fast_monopole") {
        tree->ComputeFastMonopoleForces(Nactive,Ngravcell,gravcell,cell,activepart);
      }

      // Compute all star forces for active particles
      for (j=0; j<Nactive; j++) {
        sph->ComputeStarGravForces(nbody->Nnbody,nbody->nbodydata,activepart[j]);
      }

      // Add all active particles contributions to main array
      for (j=0; j<Nactive; j++) {
        i = activelist[j];
        for (k=0; k<ndim; k++) sphdata[i].a[k]     = activepart[j].a[k];
        for (k=0; k<ndim; k++) sphdata[i].agrav[k] = activepart[j].agrav[k];
        for (k=0; k<ndim; k++) sphdata[i].a[k]     += sphdata[i].agrav[k];
        sphdata[i].gpot   = activepart[j].gpot;
        sphdata[i].dudt   = activepart[j].dudt;
        sphdata[i].div_v  = activepart[j].div_v;
        sphdata[i].active = false;
        levelneib[i]      = max(levelneib[i],activepart[j].levelneib);
      }

    }
    //=============================================================================================


    // Finally, add all contributions from distant pair-wise forces to arrays
#pragma omp critical
    for (i=0; i<sph->Nsph; i++) {
      sphdata[i].levelneib = max(sphdata[i].levelneib,levelneib[i]);
    }

    // Free-up local memory for OpenMP thread
    delete[] directlist;
    delete[] sphauxlist;
    delete[] sphlist;
    delete[] neiblist;

  }
  //===============================================================================================

  delete[] celllist;

  timing->EndTimingSection("SPH_ALL_FORCES");

  return;
}



//=================================================================================================
//  GradhSphTree::UpdateAllSphGravForces
/// Compute all local 'gather' properties of currently active particles, and
/// then compute each particle's contribution to its (active) neighbour
/// neighbour hydro forces.  Optimises the algorithm by using grid-cells to
/// construct local neighbour lists for all particles  inside the cell.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void GradhSphTree<ndim,ParticleType,TreeCell>::UpdateAllSphGravForces
 (int Nsph,                            ///< [in] No. of SPH particles
  int Ntot,                            ///< [in] No. of SPH + ghost particles
  SphParticle<ndim> *sph_gen,          ///< [inout] Pointer to SPH ptcl array
  Sph<ndim> *sph,                      ///< [in] Pointer to SPH object
  Nbody<ndim> *nbody)                  ///< [in] Pointer to N-body object
{
  int cactive;                         // No. of active cells
  int cc;                              // Aux. cell counter
  int i;                               // Particle id
  int ithread;                         // OpenMP thread id
  int j;                               // Aux. particle counter
  int jj;                              // Aux. particle counter
  int k;                               // Dimension counter
  int okflag;                          // Flag if h-rho iteration is valid
  int Ngravcell;                       // No. of gravity cells
  int Nneib;                           // No. of neighbours
  int Nneibmax;                        // Max. no. of neighbours
  FLOAT draux[ndim];                   // Aux. relative position vector
  FLOAT drsqd;                         // Distance squared
  FLOAT hrangesqdi;                    // Kernel gather extent
  FLOAT macfactor;                     // Gravity MAC factor
  FLOAT rp[ndim];                      // ..
  //TreeCell<ndim> *cell;                // Pointer to binary tree cell
  TreeCell<ndim> **celllist;           // List of pointers to binary tree cells
  ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sph_gen);

  debug2("[GradhSphTree::UpdateAllSphGravForces]");
  timing->StartTimingSection("SPH_GRAV_FORCES",2);

  // Update ghost tree smoothing length values here
  if (ghosttree->Ntot > 0) ghosttree->UpdateHmaxValues(ghosttree->celldata[0],sphdata);

  // Find list of all cells that contain active particles
#if defined (MPI_PARALLEL)
  celllist = new TreeCell<ndim>*[tree->Ncellmax];
#else
  celllist = new TreeCell<ndim>*[tree->gtot];
#endif
  cactive = tree->ComputeActiveCellList(celllist);


  // Set-up all OMP threads
  //===============================================================================================
#pragma omp parallel default(none) shared(celllist,cactive,nbody,sph,sphdata,cout) \
  private(cc,draux,drsqd,hrangesqdi,i,j,jj,k,macfactor,Ngravcell,Nneib,okflag,rp)
  {
#if defined _OPENMP
    int ithread = omp_get_thread_num();
#else
    int ithread = 0;
#endif
    int Nactive;
    int Ndirect;
    int Ndirectaux;
    int Nsphaux;
    int Nsphneib;
    int Nneibmax                   = Nneibmaxbuf[ithread];
    int Ngravcellmax               = Ngravcellmaxbuf[ithread];
    int* neiblist                  = new int[Nneibmax];
    int* sphlist                   = new int[Nneibmax];
    int *sphauxlist                = new int[Nneibmax];
    int* directlist                = new int[Nneibmax];
    int* levelneib                 = levelneibbuf[ithread];
    int* activelist                = activelistbuf[ithread];
    ParticleType<ndim>* activepart = activepartbuf[ithread];
    ParticleType<ndim>* neibpart   = neibpartbuf[ithread];
    TreeCell<ndim>* gravcell       = cellbuf[ithread];

    // Zero timestep level array
    for (i=0; i<sph->Nsph; i++) levelneib[i] = 0.0;


    // Loop over all active cells
    //=============================================================================================
#pragma omp for schedule(guided)
    for (cc=0; cc<cactive; cc++) {
      TreeCell<ndim> cell = *celllist[cc];
      macfactor = 0.0;

      // Find list of active particles in current cell
      Nactive = tree->ComputeActiveParticleList(cell,sphdata,activelist);

      // Make local copies of active particles
      for (j=0; j<Nactive; j++) activepart[j] = sphdata[activelist[j]];

      // Compute average/maximum term for computing gravity MAC
      if (gravity_mac == "eigenmac") {
        for (j=0; j<Nactive; j++) macfactor = max(macfactor,pow(1.0/activepart[j].gpot,twothirds));
      }

      // Zero/initialise all summation variables for active particles
      for (j=0; j<Nactive; j++) {
        activepart[j].div_v     = (FLOAT) 0.0;
        activepart[j].dudt      = (FLOAT) 0.0;
        activepart[j].levelneib = 0;
        activepart[j].gpot      = activepart[j].m*activepart[j].invh*sph->kernp->wpot(0.0);
        for (k=0; k<ndim; k++) activepart[j].a[k] = (FLOAT) 0.0;
        for (k=0; k<ndim; k++) activepart[j].agrav[k] = (FLOAT) 0.0;
        for (k=0; k<ndim; k++) activepart[j].adot[k] = (FLOAT) 0.0;
      }

      // Compute neighbour list for cell depending on physics options
      okflag = tree->ComputeGravityInteractionList(cell,sphdata,macfactor,Nneibmax,Ngravcellmax,
                                                   Nneib,Nsphneib,Ndirect,Ngravcell,neiblist,
                                                   sphlist,directlist,gravcell,neibpart);

      // If there are too many neighbours, reallocate the arrays and recompute the neighbour lists.
      while (okflag < 0 || Nneib > Nneibmax) {
        delete[] neibpart;
        delete[] gravcell;
        delete[] directlist;
        delete[] sphauxlist;
        delete[] sphlist;
        delete[] neiblist;
        Nneibmax                 = 2*Nneibmax;
        Ngravcellmax             = 2*Ngravcellmax;
        Nneibmaxbuf[ithread]     = Nneibmax;
        Ngravcellmaxbuf[ithread] = Ngravcellmax;
        neiblist                 = new int[Nneibmax];
        sphlist                  = new int[Nneibmax];
        sphauxlist               = new int[Nneibmax];
        directlist               = new int[Nneibmax];
        neibpartbuf[ithread]     = new ParticleType<ndim>[Nneibmax];
        cellbuf[ithread]         = new TreeCell<ndim>[Ngravcellmax];
        neibpart                 = neibpartbuf[ithread];
        gravcell                 = cellbuf[ithread];
        okflag = tree->ComputeGravityInteractionList(cell,sphdata,macfactor,Nneibmax,
                                                     Ngravcellmax,Nneib,Nsphneib,Ndirect,Ngravcell,
                                                     sphlist,neiblist,directlist,gravcell,neibpart);
      };

      // Make local copies of all potential neighbours
      for (j=0; j<Nneib; j++) neibpart[j] = sphdata[neiblist[j]];


      // Loop over all active particles in the cell
      //-------------------------------------------------------------------------------------------
      for (j=0; j<Nactive; j++) {
        i = activelist[j];

        Nsphaux = 0;
        Ndirectaux = Ndirect;
        for (k=0; k<ndim; k++) rp[k] = activepart[j].r[k];
        hrangesqdi = activepart[j].hrangesqd;

        //-----------------------------------------------------------------------------------------
        for (jj=0; jj<Nsphneib; jj++) {
          int ii = sphlist[jj];

          // Compute relative position and distance quantities for pair
          for (k=0; k<ndim; k++) draux[k] = neibpart[ii].r[k] - rp[k];
          drsqd = DotProduct(draux,draux,ndim) + small_number;

          // Record if neighbour is direct-sum or and SPH neighbour.
          // If SPH neighbour, also record max. timestep level for neighbour
          if (drsqd > hrangesqdi && drsqd >= neibpart[ii].hrangesqd) {
            directlist[Ndirectaux++] = ii;
          }
          else if (neiblist[ii] != i) {
            sphauxlist[Nsphaux++] = ii;
            levelneib[neiblist[ii]] = max(levelneib[neiblist[ii]],activepart[j].level);
          }
        }
        //-----------------------------------------------------------------------------------------



        // Compute forces between SPH neighbours (hydro and gravity)
        sph->ComputeSphGravForces(i,Nsphaux,sphauxlist,activepart[j],neibpart);

        // Compute direct gravity forces between distant particles
        sph->ComputeDirectGravForces(i,Ndirectaux,directlist,activepart[j],sphdata);

        // Compute gravitational force due to distant cells
        if (multipole == "monopole") {
          tree->ComputeCellMonopoleForces(activepart[j].gpot,activepart[j].agrav,
                                          activepart[j].r,Ngravcell,gravcell);
        }
        else if (multipole == "quadrupole") {
          tree->ComputeCellQuadrupoleForces(activepart[j].gpot,activepart[j].agrav,
                                            activepart[j].r,Ngravcell,gravcell);
        }

      }
      //-------------------------------------------------------------------------------------------


      // Compute 'fast' multipole terms here
      if (multipole == "fast_monopole") {
        tree->ComputeFastMonopoleForces(Nactive,Ngravcell,gravcell,cell,activepart);
      }

      // Compute all star forces for active particles
      for (j=0; j<Nactive; j++) {
        sph->ComputeStarGravForces(nbody->Nnbody,nbody->nbodydata,activepart[j]);
      }

      // Add all active particles contributions to main array
      for (j=0; j<Nactive; j++) {
        i = activelist[j];
        for (k=0; k<ndim; k++) sphdata[i].a[k]     = activepart[j].a[k];
        for (k=0; k<ndim; k++) sphdata[i].agrav[k] = activepart[j].agrav[k];
        for (k=0; k<ndim; k++) sphdata[i].a[k]     += sphdata[i].agrav[k];
        sphdata[i].gpot   = activepart[j].gpot;
        sphdata[i].dudt   = activepart[j].dudt;
        sphdata[i].div_v  = activepart[j].div_v;
        sphdata[i].active = false;
        levelneib[i]      = max(levelneib[i],activepart[j].levelneib);
      }

    }
    //=============================================================================================


    // Finally, add all contributions from distant pair-wise forces to arrays
#pragma omp critical
    for (i=0; i<sph->Nsph; i++) {
      sphdata[i].levelneib = max(sphdata[i].levelneib,levelneib[i]);
    }

    // Free-up local memory for OpenMP thread
    delete[] directlist;
    delete[] sphauxlist;
    delete[] sphlist;
    delete[] neiblist;

  }
  //===============================================================================================

  delete[] celllist;

  timing->EndTimingSection("SPH_GRAV_FORCES");

  return;
}



//=================================================================================================
//  GradhSphTree::UpdateAllStarGasForces
/// Calculate the gravitational acceleration on all star particles due to
/// all gas particles via the tree.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void GradhSphTree<ndim,ParticleType,TreeCell>::UpdateAllStarGasForces
 (int Nsph,                            ///< [in] No. of SPH particles
  int Ntot,                            ///< [in] No. of SPH + ghost particles
  SphParticle<ndim> *sph_gen,          ///< [inout] Pointer to SPH ptcl array
  Sph<ndim> *sph,                      ///< [in] Pointer to SPH object
  Nbody<ndim> *nbody)                  ///< [in] Pointer to N-body object
{
  int i;                               // Particle id
  int j;                               // Aux. particle counter
  int jj;                              // Aux. particle counter
  int k;                               // Dimension counter
  int okflag;                          // Flag if h-rho iteration is valid
  int Nactive;                         // No. of active particles in cell
  int Ndirect;                         // No. of direct-sum gravity particles
  int Ndirectaux;                      // ..
  int Ngravcell;                       // No. of gravity cells
  int Ngravcellmax;                    // Max. no. of gravity cells
  int Ninteract;                       // No. of interactions with hydro neibs
  int Nneib;                           // No. of neighbours
  int Nneibmax;                        // Max. no. of neighbours
  int *activelist;                     // List of active particle ids
  int *directlist;                     // List of direct sum particle ids
  int *neiblist;                       // List of neighbour ids
  FLOAT macfactor;                     // Gravity MAC factor
  TreeCell<ndim> *gravcelllist;        // List of pointers to grav. cells
  NbodyParticle<ndim> *star;           // Pointer to star particle
  ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sph_gen);


  debug2("[GradhSphTree::UpdateAllStarGasForces]");
  //timing->StartTimingSection("STAR_GAS_GRAV_FORCES",2);

  // Make list of all active stars
  Nactive = 0;
  activelist = new int[nbody->Nstar];
  for (i=0; i<nbody->Nstar; i++) {
    if (nbody->nbodydata[i]->active) activelist[Nactive++] = i;
  }


  // Set-up all OMP threads
  //===============================================================================================
#pragma omp parallel default(none) shared(activelist,Nactive,Ntot,nbody,sph,sphdata,cout) \
  private(directlist,gravcelllist,i,j,jj,k,macfactor,neiblist) \
  private(Ndirect,Ndirectaux,Ngravcell,Ngravcellmax) \
  private(Ninteract,Nneib,Nneibmax,okflag,star)
  {
#if defined _OPENMP
    int ithread = omp_get_thread_num();
#else
    int ithread = 0;
#endif
    Nneibmax = Ntot; //Nneibmaxbuf[ithread];
    Ngravcellmax = Ngravcellmaxbuf[ithread];

    neiblist = new int[Nneibmax];
    directlist = new int[Nneibmax];
    gravcelllist = new TreeCell<ndim>[Ngravcellmax];


    // Loop over all active cells
    //=============================================================================================
#pragma omp for schedule(dynamic)
    for (j=0; j<Nactive; j++) {
      i = activelist[j];
      star = nbody->nbodydata[i];

      // Compute average/maximum term for computing gravity MAC
      if (gravity_mac == "eigenmac") macfactor = pow(1.0/star->gpot,twothirds);
      else macfactor = 0.0;

      // Compute neighbour list for cell depending on physics options
      okflag = tree->ComputeStarGravityInteractionList(star,macfactor,Nneibmax,Nneibmax,
                                                       Ngravcellmax,Nneib,Ndirect,
                                                       Ngravcell,neiblist,directlist,
                                                       gravcelllist,sphdata);

      // If there are too many neighbours, reallocate the arrays and recompute the neighbour lists.
      while (okflag == -1) {
        delete[] gravcelllist;
        Ngravcellmax = 2*Ngravcellmax;
        gravcelllist = new TreeCell<ndim>[Ngravcellmax];
        okflag = tree->ComputeStarGravityInteractionList(star,macfactor,Nneibmax,Nneibmax,
                                                         Ngravcellmax,Nneib,Ndirect,
                                                         Ngravcell,neiblist,directlist,
                                                         gravcelllist,sphdata);
      };

      // Compute contributions to star force from nearby SPH particles
      nbody->CalculateDirectSPHForces(star,Nneib,Ndirect,neiblist,directlist,sph);

      // Compute gravitational force due to distant cells
      if (multipole == "monopole" || multipole == "fast_monopole") {
        tree->ComputeCellMonopoleForces(star->gpot,star->a,star->r,Ngravcell,gravcelllist);
      }
      else if (multipole == "quadrupole") {
        tree->ComputeCellQuadrupoleForces(star->gpot,star->a,star->r,Ngravcell,gravcelllist);
      }


    }
    //=============================================================================================


    // Free-up local memory for OpenMP thread
    delete[] gravcelllist;
    delete[] directlist;
    delete[] neiblist;

  }
  //===============================================================================================

  delete[] activelist;


  //timing->EndTimingSection("STAR_GAS_GRAV_FORCES");

  return;
}



//=================================================================================================
//  GradhSphTree::UpdateAllSphPeriodicForces
/// Compute all forces on active SPH particles (hydro + gravity) for periodic boundary conditions.
//=================================================================================================
template <int ndim, template<int> class ParticleType, template<int> class TreeCell>
void GradhSphTree<ndim,ParticleType,TreeCell>::UpdateAllSphPeriodicForces
 (int Nsph,                            ///< [in] No. of SPH particles
  int Ntot,                            ///< [in] No. of SPH + ghost particles
  SphParticle<ndim> *sph_gen,          ///< [inout] Pointer to SPH ptcl array
  Sph<ndim> *sph,                      ///< [in] Pointer to SPH object
  Nbody<ndim> *nbody,                  ///< [in] Pointer to N-body object
  DomainBox<ndim> &simbox,             ///< [in] Simulation domain box
  Ewald<ndim> *ewald)                  ///< [in] Ewald gravity object pointer
{
  int cactive;                         // No. of active cells
  int cc;                              // Aux. cell counter
  int i;                               // Particle id
  int j;                               // Aux. particle counter
  int jj;                              // Aux. particle counter
  int k;                               // Dimension counter
  int okflag;                          // Flag if h-rho iteration is valid
  int Nactive;                         // No. of active particles in cell
  int Ndirect;                         // No. of direct-sum gravity particles
  int Ndirectaux;                      // Aux. direct-sum neighbour counter
  int Ngravcell;                       // No. of gravity cells
  int Nneib;                           // No. of neighbours
  FLOAT macfactor;                     // Gravity MAC factor
  FLOAT draux[ndim];                   // Aux. relative position vector
  FLOAT drsqd;                         // Distance squared
  FLOAT hrangesqdi;                    // Kernel gather extent
  FLOAT rp[ndim];                      // ..
  TreeCell<ndim> **celllist;           // List of pointers to tree cells
  ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sph_gen);

  debug2("[GradhSphTree::UpdateAllSphPeriodicForces]");
  timing->StartTimingSection("SPH_ALL_PERIODIC_FORCES",2);

  // Update ghost tree smoothing length values here
  if (ghosttree->Ntot > 0) ghosttree->UpdateHmaxValues(ghosttree->celldata[0],sphdata);

  // Find list of all cells that contain active particles
#if defined (MPI_PARALLEL)
  celllist = new TreeCell<ndim>*[tree->Ncellmax];
#else
  celllist = new TreeCell<ndim>*[tree->gtot];
#endif
  cactive = tree->ComputeActiveCellList(celllist);


  // Set-up all OMP threads
  //===============================================================================================
#pragma omp parallel default(none) shared(celllist,cactive,nbody,simbox,sph,sphdata,cout)\
  private(cc,draux,drsqd,hrangesqdi,i,j,jj,k,macfactor)\
  private(Nactive,Ndirect,Ndirectaux,Ngravcell,Nneib,okflag,rp)
  {
#if defined _OPENMP
    int ithread = omp_get_thread_num();
#else
    int ithread = 0;
#endif
    int Nsphaux;
    int Nsphneib;
    int Nneibmax                   = Nneibmaxbuf[ithread];
    int Ngravcellmax               = Ngravcellmaxbuf[ithread];
    int *levelneib                 = levelneibbuf[ithread];
    int *activelist                = activelistbuf[ithread];
    int *neiblist                  = new int[Nneibmax];
    int *sphlist                   = new int[Nneibmax];
    int *sphauxlist                = new int[Nneibmax];
    int *directlist                = new int[Nneibmax];
    ParticleType<ndim>* activepart = activepartbuf[ithread];
    ParticleType<ndim>* neibpart   = neibpartbuf[ithread];
    TreeCell<ndim>* gravcell       = cellbuf[ithread];

    // Zero timestep level array
    for (i=0; i<sph->Nsph; i++) levelneib[i] = 0.0;


    // Loop over all active cells
    //=============================================================================================
#pragma omp for schedule(guided)
    for (cc=0; cc<cactive; cc++) {
      TreeCell<ndim> cell = *celllist[cc];
      macfactor = 0.0;

      // Find list of active particles in current cell
      Nactive = tree->ComputeActiveParticleList(cell,sphdata,activelist);

      // Make local copies of active particles
      for (j=0; j<Nactive; j++) activepart[j] = sphdata[activelist[j]];

      // Compute average/maximum term for computing gravity MAC
      if (gravity_mac == "eigenmac") {
        for (j=0; j<Nactive; j++) macfactor = max(macfactor,pow(1.0/activepart[j].gpot,twothirds));
      }

      // Zero/initialise all summation variables for active particles
      for (j=0; j<Nactive; j++) {
        activepart[j].div_v     = (FLOAT) 0.0;
        activepart[j].dudt      = (FLOAT) 0.0;
        activepart[j].levelneib = 0;
        activepart[j].gpot      = activepart[j].m*activepart[j].invh*sph->kernp->wpot(0.0);
        for (k=0; k<ndim; k++) activepart[j].a[k]     = (FLOAT) 0.0;
        for (k=0; k<ndim; k++) activepart[j].agrav[k] = (FLOAT) 0.0;
        for (k=0; k<ndim; k++) activepart[j].adot[k]  = (FLOAT) 0.0;
      }

      // Compute neighbour list for cell depending on physics options
      okflag = tree->ComputePeriodicGravityInteractionList
        (cell,sphdata,simbox,macfactor,Nneibmax,Ngravcellmax,Nneib,Nsphneib,Ndirect,Ngravcell,
         neiblist,sphlist,directlist,gravcell,neibpart);

      // If there are too many neighbours, reallocate the arrays and recompute the neighbour lists.
      while (okflag < 0 || Nneib > Nneibmax) {
        delete[] neibpart;
        delete[] gravcell;
        delete[] directlist;
        delete[] sphauxlist;
        delete[] sphlist;
        delete[] neiblist;
        Nneibmax                 = 2*Nneibmax;
        Ngravcellmax             = 2*Ngravcellmax;
        Nneibmaxbuf[ithread]     = Nneibmax;
        Ngravcellmaxbuf[ithread] = Ngravcellmax;
        neiblist                 = new int[Nneibmax];
        sphlist                  = new int[Nneibmax];
        sphauxlist               = new int[Nneibmax];
        directlist               = new int[Nneibmax];
        neibpartbuf[ithread]     = new ParticleType<ndim>[Nneibmax];
        cellbuf[ithread]         = new TreeCell<ndim>[Ngravcellmax];
        neibpart                 = neibpartbuf[ithread];
        gravcell                 = cellbuf[ithread];
        okflag = tree->ComputePeriodicGravityInteractionList
          (cell,sphdata,simbox,macfactor,Nneibmax,Ngravcellmax,Nneib,Nsphneib,Ndirect,Ngravcell,
           neiblist,sphlist,directlist,gravcell,neibpart);
      };


      // Loop over all active particles in the cell
      //-------------------------------------------------------------------------------------------
      for (j=0; j<Nactive; j++) {
        i = activelist[j];

        Nsphaux = 0;
        Ndirectaux = Ndirect;
        for (k=0; k<ndim; k++) rp[k] = activepart[j].r[k];
        hrangesqdi = activepart[j].hrangesqd;

        //-----------------------------------------------------------------------------------------
        for (jj=0; jj<Nsphneib; jj++) {
          int ii = sphlist[jj];

          // Compute relative position and distance quantities for pair
          for (k=0; k<ndim; k++) draux[k] = neibpart[ii].r[k] - rp[k];
          drsqd = DotProduct(draux,draux,ndim) + small_number;

          // Record if neighbour is direct-sum or and SPH neighbour.
          // If SPH neighbour, also record max. timestep level for neighbour
          if (drsqd > hrangesqdi && drsqd >= neibpart[ii].hrangesqd) {
            directlist[Ndirectaux++] = ii;
          }
          else if (neiblist[ii] != i) {
            sphauxlist[Nsphaux++] = ii;
            levelneib[neiblist[ii]] = max(levelneib[neiblist[ii]],activepart[j].level);
          }
        }
        //-----------------------------------------------------------------------------------------


        // Compute forces between SPH neighbours (hydro and gravity)
        sph->ComputeSphHydroGravForces(i,Nsphaux,sphauxlist,activepart[j],neibpart);

        // Compute direct gravity forces between distant particles
        sph->ComputeDirectGravForces(i,Ndirectaux,directlist,activepart[j],neibpart);

        // Compute gravitational force due to distant cells
        if (multipole == "monopole") {
          tree->ComputeCellMonopoleForces(activepart[j].gpot,activepart[j].agrav,
                                          activepart[j].r,Ngravcell,gravcell);
        }
        else if (multipole == "quadrupole") {
          tree->ComputeCellQuadrupoleForces(activepart[j].gpot,activepart[j].agrav,
                                            activepart[j].r,Ngravcell,gravcell);
        }

      }
      //-------------------------------------------------------------------------------------------


      // Compute 'fast' multipole terms here
      if (multipole == "fast_monopole") {
        tree->ComputeFastMonopoleForces(Nactive,Ngravcell,gravcell,cell,activepart);
      }

      // Compute all star forces for active particles
      for (j=0; j<Nactive; j++) {
        sph->ComputeStarGravForces(nbody->Nnbody,nbody->nbodydata,activepart[j]);
      }

      // Add all active particles contributions to main array
      for (j=0; j<Nactive; j++) {
        i = activelist[j];
        for (k=0; k<ndim; k++) sphdata[i].a[k]     = activepart[j].a[k];
        for (k=0; k<ndim; k++) sphdata[i].agrav[k] = activepart[j].agrav[k];
        for (k=0; k<ndim; k++) sphdata[i].a[k]     += sphdata[i].agrav[k];
        sphdata[i].gpot   = activepart[j].gpot;
        sphdata[i].dudt   = activepart[j].dudt;
        sphdata[i].div_v  = activepart[j].div_v;
        sphdata[i].active = false;
        levelneib[i]      = max(levelneib[i],activepart[j].levelneib);
      }

    }
    //=============================================================================================


    // Finally, add all contributions from distant pair-wise forces to arrays
#pragma omp critical
    for (i=0; i<sph->Nsph; i++) {
      sphdata[i].levelneib = max(sphdata[i].levelneib,levelneib[i]);
    }

    // Free-up local memory for OpenMP thread
    delete[] directlist;
    delete[] sphauxlist;
    delete[] sphlist;
    delete[] neiblist;

  }
  //===============================================================================================

  delete[] celllist;

  timing->EndTimingSection("SPH_ALL_FORCES");

  return;
}




template class GradhSphTree<1,GradhSphParticle,KDTreeCell>;
template class GradhSphTree<2,GradhSphParticle,KDTreeCell>;
template class GradhSphTree<3,GradhSphParticle,KDTreeCell>;
template class GradhSphKDTree<1,GradhSphParticle,KDTreeCell>;
template class GradhSphKDTree<2,GradhSphParticle,KDTreeCell>;
template class GradhSphKDTree<3,GradhSphParticle,KDTreeCell>;

template class GradhSphTree<1,GradhSphParticle,OctTreeCell>;
template class GradhSphTree<2,GradhSphParticle,OctTreeCell>;
template class GradhSphTree<3,GradhSphParticle,OctTreeCell>;
template class GradhSphOctTree<1,GradhSphParticle,OctTreeCell>;
template class GradhSphOctTree<2,GradhSphParticle,OctTreeCell>;
template class GradhSphOctTree<3,GradhSphParticle,OctTreeCell>;
