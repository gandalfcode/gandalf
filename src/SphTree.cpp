//=============================================================================
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
//=============================================================================


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



//=============================================================================
//  SphTree::SphTree
/// SphTree constructor.  Initialises various variables.
//=============================================================================
template <int ndim, template<int> class ParticleType>
SphTree<ndim,ParticleType>::SphTree(int Nleafmaxaux, FLOAT thetamaxsqdaux, 
                             FLOAT kernrangeaux, FLOAT macerroraux, 
                             string gravity_mac_aux, string multipole_aux):
  SphNeighbourSearch<ndim>()
{
  allocated_buffer = false;
  neibcheck        = true;
  Ntot             = 0;
  Ntotmax          = 0;
  Ntotmaxold       = 0;
  Nleafmax         = Nleafmaxaux;
  kernrange        = kernrangeaux;
  thetamaxsqd      = thetamaxsqdaux;
  invthetamaxsqd   = 1.0/thetamaxsqdaux;
  gravity_mac      = gravity_mac_aux;
  macerror         = macerroraux;
  multipole        = multipole_aux;
#if defined _OPENMP
  Nthreads         = omp_get_max_threads();
#else
  Nthreads         = 1;
#endif

  // Set-up main tree object
  tree = new KDTree<ndim,ParticleType>(Nleafmaxaux, thetamaxsqdaux, 
                                       kernrangeaux, macerroraux, 
                                       gravity_mac_aux, multipole_aux);
}



//=============================================================================
//  SphTree::~SphTree
/// SphTree destructor.  Deallocates tree memory upon object destruction.
//=============================================================================
template <int ndim, template<int> class ParticleType>
SphTree<ndim,ParticleType>::~SphTree()
{
  if (tree->allocated_tree) {
    DeallocateMemory();
    tree->DeallocateTreeMemory();
  }
}



//=============================================================================
//  SphTree::AllocateMemory
/// Allocate memory for binary tree as requested.  If more memory is required 
/// than currently allocated, tree is deallocated and reallocated here.
//=============================================================================
template <int ndim, template<int> class ParticleType>
void SphTree<ndim,ParticleType>::AllocateMemory(Sph<ndim> *sph)
{
  int ithread;                      // Thread id number

  debug2("[SphTree::AllocateMemory]");

  if (!tree->allocated_tree || Ntotmax > Ntotmaxold) {
    if (tree->allocated_tree) {
      DeallocateMemory();
      tree->DeallocateTreeMemory();
    }

    if (!allocated_buffer) {
      Nneibmaxbuf = new int[Nthreads];
      Ndirectmaxbuf = new int[Nthreads];
      Ngravcellmaxbuf = new int[Nthreads];
      levelneibbuf = new int*[Nthreads];
      activelistbuf = new int*[Nthreads];
      activepartbuf = new ParticleType<ndim>*[Nthreads];
      neibpartbuf = new ParticleType<ndim>*[Nthreads];

      for (ithread=0; ithread<Nthreads; ithread++) {
	Nneibmaxbuf[ithread] = max(1,4*sph->Ngather);
	Ndirectmaxbuf[ithread] = max(1,4*sph->Ngather);
	Ngravcellmaxbuf[ithread] = max(1,4*sph->Ngather);
	levelneibbuf[ithread] = new int[Ntotmax];
        activelistbuf[ithread] = new int[Nleafmax];
	activepartbuf[ithread] = new ParticleType<ndim>[Nleafmax];
	neibpartbuf[ithread] = new ParticleType<ndim>[Nneibmaxbuf[ithread]];
      }
      allocated_buffer = true;
    }

  }

  return;
}



//=============================================================================
//  SphTree::DeallocateTreeMemory
/// Deallocates all binary tree memory
//=============================================================================
template <int ndim, template<int> class ParticleType>
void SphTree<ndim,ParticleType>::DeallocateMemory(void)
{
  int ithread;                      // Thread id number

  debug2("[SphTree::DeallocateTreeMemory]");

  if (tree->allocated_tree) {

    for (ithread=0; ithread<Nthreads; ithread++) {
      delete[] neibpartbuf[ithread];
      delete[] activepartbuf[ithread];
      delete[] levelneibbuf[ithread];
    }
    delete[] neibpartbuf;
    delete[] activepartbuf;
    delete[] activelistbuf;
    delete[] levelneibbuf;
    delete[] Ngravcellmaxbuf;
    delete[] Ndirectmaxbuf;
    delete[] Nneibmaxbuf;

    tree->DeallocateTreeMemory();

  }

  return;
}



//=============================================================================
//  SphTree::BuildTree
/// Main routine to control how the tree is built, re-stocked and interpolated 
/// during each timestep.
//=============================================================================
template <int ndim, template<int> class ParticleType>
void SphTree<ndim,ParticleType>::BuildTree
(bool rebuild_tree,                 ///< Flag to rebuild tree
 int n,                             ///< Integer time
 int ntreebuildstep,                ///< Tree build frequency
 int ntreestockstep,                ///< Tree stocking frequency
 int Npart,                         ///< No. of particles
 int Npartmax,                      ///< Max. no. of particles
 SphParticle<ndim> *sph_gen,        ///< Particle data array
 Sph<ndim> *sph,                    ///< Pointer to SPH object
 FLOAT timestep)                    ///< Smallest physical timestep
{
  int i;                            // Particle counter
  int k;                            // Dimension counter
  ParticleType<ndim> *sphdata = static_cast<ParticleType<ndim>* > (sph_gen);

  debug2("[SphTree::BuildTree]");
  timing->StartTimingSection("BUILD_TREE",2);

  // Activate nested parallelism for tree building routines
#ifdef _OPENMP
  omp_set_nested(1);
#endif


  // For tree rebuild steps
  //---------------------------------------------------------------------------
  if (n%ntreebuildstep == 0 || rebuild_tree) {

    // Delete any dead particles from main SPH arrays before we re-build tree
    sph->DeleteDeadParticles();

    Ntotold = Ntot;
    //Nsph = sph->Nsph;
    Ntot = sph->Ntot;
    Ntotmaxold = Ntotmax;
    Ntotmax = max(Ntotmax,Ntot);
    Ntotmax = max(Ntotmax,sph->Nsphmax);
    tree->Ntot = sph->Ntot;
    tree->Ntotmaxold = tree->Ntotmax;
    tree->Ntotmax = max(tree->Ntotmax,tree->Ntot);
    tree->Ntotmax = max(tree->Ntotmax,sph->Nsphmax);
    AllocateMemory(sph);

    tree->BuildTree(Npart,Npartmax,sphdata,timestep);
    
  }

  // Else stock the tree
  //---------------------------------------------------------------------------
  else if (n%ntreestockstep == 0) {

    tree->StockTree(tree->kdcell[0],sphdata);

  }

  // Otherwise simply extrapolate tree cell properties
  //---------------------------------------------------------------------------
  else {

    //ExtrapolateCellProperties(kdcell[0],timestep);
    tree->ExtrapolateCellProperties(timestep);

  }
  //---------------------------------------------------------------------------

#ifdef _OPENMP
  omp_set_nested(0);
#endif

  timing->EndTimingSection("BUILD_TREE");


  return;
}



//=============================================================================
//  SphTree::UpdateAllSphProperties
/// Compute all local 'gather' properties of currently active particles, and 
/// then compute each particle's contribution to its (active) neighbour 
/// neighbour hydro forces.  Optimises the algorithm by using grid-cells to 
/// construct local neighbour lists for all particles  inside the cell.
//=============================================================================
template <int ndim, template<int> class ParticleType>
void SphTree<ndim,ParticleType>::UpdateAllSphProperties
(int Nsph,                          ///< [in] No. of SPH particles
 int Ntot,                          ///< [in] No. of SPH + ghost particles
 SphParticle<ndim> *sph_gen,        ///< [inout] Pointer to SPH ptcl array
 Sph<ndim> *sph,                    ///< [in] Pointer to SPH object
 Nbody<ndim> *nbody)                ///< [in] Pointer to N-body object
{
  int celldone;                    // ..
  int okflag;                      // ..
  int cc;                          // Aux. cell counter
  int cactive;                     // No. of active
  int i;                           // Particle id
  int j;                           // Aux. particle counter
  int jj;                          // Aux. particle counter
  int k;                           // Dimension counter
  int Nactive;                     // No. of active particles in cell
  int Ngather;                     // No. of near gather neighbours
  int Nneib;                       // No. of neighbours
  int Nneibmax;                    // Max. no. of neighbours
  int *activelist;                 // List of active particle ids
  int *neiblist;                   // List of neighbour ids
  FLOAT draux[ndim];               // Aux. relative position vector var
  FLOAT drsqdaux;                  // Distance squared
  FLOAT hrangesqd;                 // Kernel extent
  FLOAT hmax;                      // Maximum smoothing length
  FLOAT rp[ndim];                  // Local copy of particle position
  FLOAT *drsqd;                    // Position vectors to gather neibs
  FLOAT *gpot;                     // Potential for particles
  FLOAT *gpot2;                    // ..
  FLOAT *m;                        // Distances to gather neibs
  FLOAT *m2;                       // ..
  FLOAT *mu;                       // mass*u for gather neibs
  FLOAT *mu2;                      // ..
  FLOAT *r;                        // Positions of neibs
  KDTreeCell<ndim> *cell;      // Pointer to binary tree cell
  KDTreeCell<ndim> **celllist; // List of binary cell pointers
  ParticleType<ndim> *activepart;   // ..
  ParticleType<ndim> *sphdata = static_cast<ParticleType<ndim>* > (sph_gen);

  int Nneibcount = 0;
  int ithread;

  debug2("[SphTree::UpdateAllSphProperties]");
  timing->StartTimingSection("SPH_PROPERTIES",2);

  // Find list of all cells that contain active particles
  celllist = new KDTreeCell<ndim>*[tree->gtot];
  cactive = tree->ComputeActiveCellList(celllist);


  // Set-up all OMP threads
  //===========================================================================
#pragma omp parallel default(none) private(activelist,activepart,cc,cell)\
  private(celldone,draux,drsqd,drsqdaux,gpot,gpot2,hmax,hrangesqd,i,ithread,j)\
  private(jj,k,okflag,m,mu,m2,mu2,Nactive,neiblist,Ngather,Nneib,Nneibmax,r,rp)	\
  shared(cactive,celllist,cout,nbody,sph,sphdata) reduction(+:Nneibcount)
  {
#if defined _OPENMP
    ithread = omp_get_thread_num();
#else
    ithread = 0;
#endif
    Nneibmax = Nneibmaxbuf[ithread];

    //int activelist[Nleafmax];
    //ParticleType<ndim> activepart[Nleafmax];
    activelist = new int[Nleafmax];
    activepart = new ParticleType<ndim>[Nleafmax];
    //activelist = activelistbuf[ithread];
    //activepart = activepartbuf[ithread];

    neiblist = new int[Nneibmax];
    gpot = new FLOAT[Nneibmax];
    gpot2 = new FLOAT[Nneibmax];
    drsqd = new FLOAT[Nneibmax];
    m = new FLOAT[Nneibmax];
    m2 = new FLOAT[Nneibmax];
    mu = new FLOAT[Nneibmax];
    mu2 = new FLOAT[Nneibmax];
    r = new FLOAT[Nneibmax*ndim];

    // Loop over all active cells
    //=========================================================================
#pragma omp for schedule(guided)
    for (cc=0; cc<cactive; cc++) {
      cell = celllist[cc];
      celldone = 1;
      hmax = cell->hmax;


      // If hmax is too small so the neighbour lists are invalid, make hmax
      // larger and then recompute for the current active cell.
      //-----------------------------------------------------------------------
      do {
        hmax = 1.05*hmax;
        celldone = 1;

        // Find list of active particles in current cell
        Nactive = tree->ComputeActiveParticleList(cell,sphdata,activelist);

	for (j=0; j<Nactive; j++)
	  activepart[j] = sphdata[activelist[j]];

        // Compute neighbour list for cell depending on physics options
        Nneib = tree->ComputeGatherNeighbourList(cell,Nneibmax,neiblist,
                                                 hmax,sphdata);

        // If there are too many neighbours, reallocate the arrays and
        // recompute the neighbour lists.
        while (Nneib == -1) {
          delete[] r;
          delete[] mu2;
          delete[] mu;
          delete[] m2;
          delete[] m;
          delete[] drsqd;
          delete[] gpot2;
          delete[] gpot;
          delete[] neiblist;
          Nneibmax = 2*Nneibmax; 
          //Nneibmaxbuf[ithread] = Nneibmax;
          //Ndirectmaxbuf[ithread] *= 2;
          //Ngravcellmaxbuf[ithread] *= 2;
          neiblist = new int[Nneibmax];
          gpot = new FLOAT[Nneibmax];
          gpot2 = new FLOAT[Nneibmax];
          drsqd = new FLOAT[Nneibmax];
          m = new FLOAT[Nneibmax];
          m2 = new FLOAT[Nneibmax];
          mu = new FLOAT[Nneibmax];
          mu2 = new FLOAT[Nneibmax];
          r = new FLOAT[Nneibmax*ndim];
          Nneib = tree->ComputeGatherNeighbourList(cell,Nneibmax,neiblist,
                                                   hmax,sphdata);
        };


        // Make local copies of important neib information (mass and position)
        for (jj=0; jj<Nneib; jj++) {
          j = neiblist[jj];
          //assert(j >= 0 && j < sph->Ntot);
          gpot[jj] = sphdata[j].gpot;
          m[jj] = sphdata[j].m;
          mu[jj] = sphdata[j].m*sphdata[j].u;
          for (k=0; k<ndim; k++) r[ndim*jj + k] = (FLOAT) sphdata[j].r[k];
        }


        // Loop over all active particles in the cell
        //---------------------------------------------------------------------
        for (j=0; j<Nactive; j++) {
	  Nneibcount += Nneib;
          i = activelist[j];
          //assert(i >= 0 && i < sph->Nsph);
          for (k=0; k<ndim; k++) rp[k] = activepart[j].r[k]; //data[i].r[k];

          // Set gather range as current h multiplied by some tolerance factor
          hrangesqd = kernp->kernrangesqd*hmax*hmax;
          Ngather = 0;

          // Compute distance (squared) to all
          //-------------------------------------------------------------------
          for (jj=0; jj<Nneib; jj++) {
            for (k=0; k<ndim; k++) draux[k] = r[ndim*jj + k] - rp[k];
            drsqdaux = DotProduct(draux,draux,ndim);

            // Record distance squared for all potential gather neighbours
            if (drsqdaux <= hrangesqd) {
              gpot[Ngather] = gpot[jj];
              drsqd[Ngather] = drsqdaux;
              m2[Ngather] = m[jj];
              mu2[Ngather] = mu[jj];
              Ngather++;
            }

          }
          //-------------------------------------------------------------------

          // Validate that gather neighbour list is correct
#if defined(VERIFY_ALL)
          if (neibcheck) 
	    this->CheckValidNeighbourList(sph,i,Nneib,neiblist,"gather");
#endif

          // Compute smoothing length and other gather properties for ptcl i
          okflag = sph->ComputeH(i,Ngather,hmax,m2,mu2,
                                 drsqd,gpot,activepart[j],nbody);

          // If h-computation is invalid, then break from loop and recompute
          // larger neighbour lists
          if (okflag == 0) {
            celldone = 0;
            break;
          }

        }
        //---------------------------------------------------------------------

      } while (celldone == 0);
      //-----------------------------------------------------------------------

      for (j=0; j<Nactive; j++) 
	sphdata[activelist[j]] = activepart[j];

    }
    //=========================================================================

    // Free-up all memory
    delete[] r;
    delete[] mu2;
    delete[] mu;
    delete[] m2;
    delete[] m;
    delete[] drsqd;
    delete[] gpot2;
    delete[] gpot;
    delete[] neiblist;
    delete[] activepart;
    delete[] activelist;

  }
  //===========================================================================

  delete[] celllist;

  timing->EndTimingSection("SPH_PROPERTIES");

  return;
}



//=============================================================================
//  SphTree::UpdateAllSphHydroForces
/// Compute all local 'gather' properties of currently active particles, and 
/// then compute each particle's contribution to its (active) neighbour 
/// neighbour hydro forces.  Optimises the algorithm by using grid-cells to 
/// construct local neighbour lists for all particles  inside the cell.
//=============================================================================
template <int ndim, template<int> class ParticleType>
void SphTree<ndim,ParticleType>::UpdateAllSphHydroForces
(int Nsph,                          ///< [in] No. of SPH particles
 int Ntot,                          ///< [in] No. of SPH + ghost particles
 SphParticle<ndim> *sph_gen,        ///< [inout] Pointer to SPH ptcl array
 Sph<ndim> *sph,                    ///< [in] Pointer to SPH object
 Nbody<ndim> *nbody)                ///< [in] Pointer to N-body object
{
  int cactive;                      // No. of active cells
  int cc;                           // Aux. cell counter
  int i;                            // Particle id
  int j;                            // Aux. particle counter
  int jj;                           // Aux. particle counter
  int k;                            // Dimension counter
  int Nactive;                      // No. of active particles in cell
  int Ninteract;                    // No. of near gather neighbours
  int Nneib;                        // No. of neighbours
  int Nneibmax;                     // Max. no. of neighbours
  int *activelist;                  // List of active particle ids
  int *interactlist;                // List of interactng SPH neighbours
  int *neiblist;                    // List of neighbour ids
  FLOAT draux[ndim];                // Aux. relative position vector
  FLOAT drsqd;                      // Distance squared
  FLOAT hrangesqdi;                 // Kernel gather extent
  FLOAT rp[ndim];                   // Local copy of particle position
  FLOAT *dr;                        // Array of relative position vectors
  FLOAT *drmag;                     // Array of neighbour distances
  FLOAT *invdrmag;                  // Array of 1/drmag between particles
  KDTreeCell<ndim> *cell;           // Pointer to binary tree cell
  KDTreeCell<ndim> **celllist;      // List of binary tree pointers
  ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sph_gen);

  int ithread;
  int Nneibcount = 0;
  int *levelneib;
  ParticleType<ndim> *activepart;
  ParticleType<ndim> *neibpart;

  debug2("[SphTree::UpdateAllSphHydroForces]");
  timing->StartTimingSection("SPH_HYDRO_FORCES",2);

  // Update tree smoothing length values here
  tree->UpdateHmaxValues(tree->kdcell[0],sphdata);


  // Find list of all cells that contain active particles
  celllist = new KDTreeCell<ndim>*[tree->gtot];
  cactive = tree->ComputeActiveCellList(celllist);


  // Set-up all OMP threads
  //===========================================================================
#pragma omp parallel default(none) private(activelist,activepart,cc,cell,dr)\
  private(draux,drmag,drsqd,hrangesqdi,i,interactlist,ithread,invdrmag,j,jj,k)\
  private(levelneib,Nactive,neiblist,neibpart,Ninteract,Nneib,Nneibmax,rp) \
  shared(cactive,celllist,sphdata,nbody,sph) reduction(+:Nneibcount)
  {
#if defined _OPENMP
    ithread = omp_get_thread_num();
#else
    ithread = 0;
#endif
    Nneibmax = Nneibmaxbuf[ithread];
    activelist = activelistbuf[ithread];    
    activepart = activepartbuf[ithread];
    neibpart = neibpartbuf[ithread];
    levelneib = levelneibbuf[ithread];
    for (i=0; i<sph->Nsph; i++) levelneib[i] = 0;

    //activelist = new int[Nleafmax];
    neiblist = new int[Nneibmax];
    interactlist = new int[Nneibmax];
    dr = new FLOAT[Nneibmax*ndim];
    drmag = new FLOAT[Nneibmax];
    invdrmag = new FLOAT[Nneibmax];


    // Loop over all active cells
    //=========================================================================
#pragma omp for schedule(guided)
    for (cc=0; cc<cactive; cc++) {
      cell = celllist[cc];

      // Find list of active particles in current cell
      Nactive = tree->ComputeActiveParticleList(cell,sphdata,activelist);

      // Make local copies of active particles
      for (j=0; j<Nactive; j++) {
        //assert(activelist[j] >= 0 && activelist[j] < sph->Nsph);
        activepart[j] = sphdata[activelist[j]];
        activepart[j].div_v = (FLOAT) 0.0;
        activepart[j].dudt = (FLOAT) 0.0;
        activepart[j].levelneib = 0;
        for (k=0; k<ndim; k++) activepart[j].a[k] = (FLOAT) 0.0;
      }

      // Compute neighbour list for cell depending on physics options
      Nneib = tree->ComputeNeighbourList(cell,Nneibmax,neiblist,sphdata);

      // If there are too many neighbours, reallocate the arrays and
      // recompute the neighbour list.
      while (Nneib == -1) {
        delete[] neibpartbuf[ithread];
        delete[] invdrmag;
        delete[] drmag;
        delete[] dr;
        delete[] interactlist;
        delete[] neiblist;
        Nneibmax = 2*Nneibmax;
	Nneibmaxbuf[ithread] = Nneibmax;
	Ndirectmaxbuf[ithread] *= 2;
	Ngravcellmaxbuf[ithread] *= 2;
        neiblist = new int[Nneibmax];
        interactlist = new int[Nneibmax];
        dr = new FLOAT[Nneibmax*ndim];
        drmag = new FLOAT[Nneibmax];
        invdrmag = new FLOAT[Nneibmax];
        neibpartbuf[ithread] = new ParticleType<ndim>[Nneibmax]; 
	neibpart = neibpartbuf[ithread];
	//neibpart = new ParticleType<ndim>[Nneibmax];
        Nneib = tree->ComputeNeighbourList(cell,Nneibmax,
                                           neiblist,sphdata);
      };

      // Make local copies of all potential neighbours
      for (j=0; j<Nneib; j++) {
        //assert(neiblist[j] >= 0 && neiblist[j] < sph->Ntot);
        neibpart[j] = sphdata[neiblist[j]];
      }


      // Loop over all active particles in the cell
      //-----------------------------------------------------------------------
      for (j=0; j<Nactive; j++) {
        i = activelist[j];

        for (k=0; k<ndim; k++) rp[k] = activepart[j].r[k];
        hrangesqdi = activepart[j].hrangesqd;
        Ninteract = 0;

        // Validate that gather neighbour list is correct
#if defined(VERIFY_ALL)
        if (neibcheck) 
	  this->CheckValidNeighbourList(sph,i,Nneib,neiblist,"all");
#endif

        // Compute distances and the inverse between the current particle
        // and all neighbours here, for both gather and inactive scatter neibs.
        // Only consider particles with j > i to compute pair forces once
        // unless particle j is inactive.
        //---------------------------------------------------------------------
        for (jj=0; jj<Nneib; jj++) {

	  // Skip current active particle
	  if (neiblist[jj] == i) continue;

          for (k=0; k<ndim; k++) draux[k] = neibpart[jj].r[k] - rp[k];
          drsqd = DotProduct(draux,draux,ndim) + small_number;

          // Compute relative position and distance quantities for pair
          if (drsqd <= hrangesqdi || drsqd <= neibpart[jj].hrangesqd) {
            drmag[Ninteract] = sqrt(drsqd);
            invdrmag[Ninteract] = (FLOAT) 1.0/drmag[Ninteract];
            for (k=0; k<ndim; k++)
              dr[Ninteract*ndim + k] = draux[k]*invdrmag[Ninteract];
            interactlist[Ninteract] = jj;
            Ninteract++;
	    levelneib[neiblist[jj]] = 
	      max(levelneib[neiblist[jj]],activepart[j].level);
	  }

        }
        //---------------------------------------------------------------------

        // Compute all neighbour contributions to hydro forces
        sph->ComputeSphHydroForces(i,Ninteract,interactlist,
				   drmag,invdrmag,dr,activepart[j],neibpart);

	Nneibcount += Ninteract;

      }
      //-----------------------------------------------------------------------


      // Compute all star forces for active particles
      if (nbody->Nnbody > 0) {
	for (j=0; j<Nactive; j++)
	  sph->ComputeStarGravForces(nbody->Nnbody,nbody->nbodydata,
				     activepart[j]);
      }


      // Add all active particles contributions to main array
      for (j=0; j<Nactive; j++) {
    	i = activelist[j];
        for (k=0; k<ndim; k++) sphdata[i].a[k] = activepart[j].a[k];
        sphdata[i].dudt = activepart[j].dudt;
        sphdata[i].div_v = activepart[j].div_v;
	sphdata[i].active = false;
        levelneib[i] = max(levelneib[i],activepart[j].levelneib);
      }

    }
    //=========================================================================


    // Finally, add all contributions from distant pair-wise forces to arrays
#pragma omp critical
    for (i=0; i<sph->Nsph; i++) {
      sphdata[i].levelneib = 
	max(sphdata[i].levelneib,levelneib[i]);
    }


    // Free-up local memory for OpenMP thread
    //delete[] neibpart;
    delete[] invdrmag;
    delete[] drmag;
    delete[] dr;
    delete[] interactlist;
    delete[] neiblist;
    //delete[] activelist;

  }
  //===========================================================================

  delete[] celllist;

  timing->EndTimingSection("SPH_HYDRO_FORCES");

  return;
}



//=============================================================================
//  SphTree::UpdateAllSphForces
/// Compute all local 'gather' properties of currently active particles, and 
/// then compute each particle's contribution to its (active) neighbour 
/// neighbour hydro forces.  Optimises the algorithm by using grid-cells to 
/// construct local neighbour lists for all particles  inside the cell.
//=============================================================================
template <int ndim, template<int> class ParticleType>
void SphTree<ndim,ParticleType>::UpdateAllSphForces
(int Nsph,                          ///< [in] No. of SPH particles
 int Ntot,                          ///< [in] No. of SPH + ghost particles
 SphParticle<ndim> *sph_gen,        ///< [inout] Pointer to SPH ptcl array
 Sph<ndim> *sph,                    ///< [in] Pointer to SPH object
 Nbody<ndim> *nbody)                ///< [in] Pointer to N-body object
{
  int cactive;                      // No. of active cells
  int cc;                           // Aux. cell counter
  int i;                            // Particle id
  int ithread;                      // OpenMP thread id
  int j;                            // Aux. particle counter
  int jj;                           // Aux. particle counter
  int k;                            // Dimension counter
  int okflag;                       // Flag if h-rho iteration is valid
  int Nactive;                      // No. of active particles in cell
  int Ndirect;                      // No. of direct-sum gravity particles
  int Ndirectaux;                   // ..
  int Ndirectmax;                   // Max. no. of direct sum particles
  int Ngravcell;                    // No. of gravity cells
  int Ngravcellmax;                 // Max. no. of gravity cells
  int Ninteract;                    // No. of interactions with hydro neibs
  int Nneib;                        // No. of neighbours
  int Nneibmax;                     // Max. no. of neighbours
  int *activelist;                  // List of active particle ids
  int *directlist;                  // List of direct sum particle ids
  int *interactlist;                // List of interacting neighbour ids
  int *levelneib;                   // ..
  int *neiblist;                    // List of neighbour ids
  FLOAT macfactor;                  // Gravity MAC factor
  FLOAT draux[ndim];                // Aux. relative position vector
  FLOAT drsqd;                      // Distance squared
  FLOAT hrangesqdi;                 // Kernel gather extent
  FLOAT rp[ndim];                   // .. 
  KDTreeCell<ndim> *cell;       // Pointer to binary tree cell
  KDTreeCell<ndim> **celllist;  // List of pointers to binary tree cells
  KDTreeCell<ndim> **gravcelllist; // List of pointers to grav. cells
  ParticleType<ndim> *activepart;    // ..
  ParticleType<ndim> *neibpart;      // ..
  ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sph_gen);

  int Nactivecount = 0;
  int Ncellcount = 0;
  int Ndirectcount = 0;
  int Nneibcount = 0;


  debug2("[SphTree::UpdateAllSphForces]");
  timing->StartTimingSection("SPH_ALL_FORCES",2);

  // Update tree smoothing length values here
  tree->UpdateHmaxValues(tree->kdcell[0],sphdata);

  // Find list of all cells that contain active particles
  celllist = new KDTreeCell<ndim>*[tree->gtot];
  cactive = tree->ComputeActiveCellList(celllist);


  // Set-up all OMP threads
  //===========================================================================
#pragma omp parallel default(none) shared(celllist,cactive,nbody,sph,cout) \
  private(activepart,activelist,cc,cell,directlist,draux,drsqd,gravcelllist) \
  private(hrangesqdi,i,interactlist,ithread,j,jj,k,levelneib,macfactor)	\
  private(neiblist,neibpart,Nactive,Ndirect,Ndirectaux,Ndirectmax) \
  private(Ngravcell,Ngravcellmax,Ninteract,Nneib,Nneibmax,okflag,rp) \
  reduction(+:Nactivecount,Ncellcount,Ndirectcount,Nneibcount)
  {
#if defined _OPENMP
    ithread = omp_get_thread_num();
#else
    ithread = 0;
#endif
    Nneibmax = Nneibmaxbuf[ithread];
    Ndirectmax = Ndirectmaxbuf[ithread];
    Ngravcellmax = Ngravcellmaxbuf[ithread];

    levelneib = levelneibbuf[ithread];
    activelist = activelistbuf[ithread];
    activepart = activepartbuf[ithread];
    neibpart = neibpartbuf[ithread];
    neiblist = new int[Nneibmax];
    interactlist = new int[Nneibmax];
    directlist = new int[Ndirectmax];
    gravcelllist = new KDTreeCell<ndim>*[Ngravcellmax];

    // Zero timestep level array
    for (i=0; i<sph->Nsph; i++) levelneib[i] = 0.0;


    // Loop over all active cells
    //=========================================================================
#pragma omp for schedule(guided)
    for (cc=0; cc<cactive; cc++) {
      cell = celllist[cc];
      macfactor = 0.0;

      // Find list of active particles in current cell
      Nactive = tree->ComputeActiveParticleList(cell,sphdata,activelist);
      Nactivecount += Nactive;

      // Make local copies of active particles
      for (j=0; j<Nactive; j++) activepart[j] = sphdata[activelist[j]];

      // Compute average/maximum term for computing gravity MAC
      if (gravity_mac == "eigenmac") {
	for (j=0; j<Nactive; j++) 
	  macfactor = max(macfactor,pow(1.0/activepart[j].gpot,twothirds));
      }

      // Zero/initialise all summation variables for active particles
      for (j=0; j<Nactive; j++) {
        activepart[j].div_v = (FLOAT) 0.0;
        activepart[j].dudt = (FLOAT) 0.0;
        activepart[j].levelneib = 0;
        activepart[j].gpot = activepart[j].m*activepart[j].invh*sph->kernp->wpot(0.0);
        for (k=0; k<ndim; k++) activepart[j].a[k] = (FLOAT) 0.0;
        for (k=0; k<ndim; k++) activepart[j].agrav[k] = (FLOAT) 0.0;
      }


      // Compute neighbour list for cell depending on physics options
      okflag = tree->ComputeGravityInteractionList(cell,macfactor,
					           Nneibmax,Ndirectmax,
                                                   Ngravcellmax,Nneib,Ndirect,
                                                   Ngravcell,neiblist,directlist,
                                                   gravcelllist,sphdata);

      // If there are too many neighbours, reallocate the arrays and
      // recompute the neighbour lists.
      while (okflag < 0) {
	delete[] neibpart;  //delete[] neibpartbuf[ithread];
	delete[] gravcelllist;
	delete[] directlist;
	delete[] interactlist;
	delete[] neiblist;
	Nneibmax = 2*Nneibmax;
	Ndirectmax = 2*Ndirectmax;
	Ngravcellmax = 2*Ngravcellmax;
	Nneibmaxbuf[ithread] = Nneibmax;
	Ndirectmaxbuf[ithread] = Ndirectmax;
	Ngravcellmaxbuf[ithread] = Ngravcellmax;
	neiblist = new int[Nneibmax];
	interactlist = new int[Nneibmax];
	directlist = new int[Ndirectmax];
	gravcelllist = new KDTreeCell<ndim>*[Ngravcellmax];
	neibpart = new ParticleType<ndim>[Nneibmax];
	neibpartbuf[ithread] = new ParticleType<ndim>[Nneibmax]; 
	neibpart = neibpartbuf[ithread];
	okflag = tree->ComputeGravityInteractionList(cell,macfactor,
					             Nneibmax,Ndirectmax,
					             Ngravcellmax,Nneib,Ndirect,
					             Ngravcell,neiblist,directlist,
					             gravcelllist,sphdata);
      };


      // Make local copies of all potential neighbours
      for (j=0; j<Nneib; j++) {
        //assert(neiblist[j] >= 0 && neiblist[j] < sph->Ntot);
        neibpart[j] = sphdata[neiblist[j]];
      }

      // Loop over all active particles in the cell
      //-----------------------------------------------------------------------
      for (j=0; j<Nactive; j++) {
        i = activelist[j];

        // Determine SPH neighbour interaction list 
        // (to ensure we don't compute pair-wise forces twice)
        Ninteract = 0;
	Ndirectaux = Ndirect;
        for (k=0; k<ndim; k++) rp[k] = activepart[j].r[k];
        hrangesqdi = activepart[j].hrangesqd;

	//---------------------------------------------------------------------
        for (jj=0; jj<Nneib; jj++) {

          for (k=0; k<ndim; k++) draux[k] = neibpart[jj].r[k] - rp[k];
          drsqd = DotProduct(draux,draux,ndim) + small_number;

          // Compute relative position and distance quantities for pair
          if (drsqd > hrangesqdi && drsqd >= neibpart[jj].hrangesqd)
	    directlist[Ndirectaux++] = neiblist[jj];
	  else if (i != neiblist[jj]) {
            interactlist[Ninteract++] = jj;
	    levelneib[neiblist[jj]] = 
	      max(levelneib[neiblist[jj]],activepart[j].level);
	  }
        }
	//---------------------------------------------------------------------


        // Compute forces between SPH neighbours (hydro and gravity)
        sph->ComputeSphHydroGravForces(i,Ninteract,interactlist,
                                       activepart[j],neibpart);

        // Compute direct gravity forces between distant particles
        sph->ComputeDirectGravForces(i,Ndirectaux,directlist,
                                     activepart[j],sphdata);


        // Compute gravitational force due to distant cells
        if (multipole == "monopole")
          tree->ComputeCellMonopoleForces(activepart[j].gpot,
                                          activepart[j].agrav,
				          activepart[j].r,
                                          Ngravcell,gravcelllist);
        else if (multipole == "quadrupole")
          tree->ComputeCellQuadrupoleForces(activepart[j].gpot,
                                            activepart[j].agrav,
                                            activepart[j].r,
                                            Ngravcell,gravcelllist);

	Nneibcount += Ninteract;
	Ndirectcount += Ndirectaux;
	Ncellcount += Ngravcell;

      }
      //-----------------------------------------------------------------------


      // Compute 'fast' multipole terms here
      if (multipole == "fast_monopole")
	tree->ComputeFastMonopoleForces(Nactive,Ngravcell,gravcelllist,
				        cell,activepart);

      // Compute all star forces for active particles
      for (j=0; j<Nactive; j++)
	sph->ComputeStarGravForces(nbody->Nnbody,nbody->nbodydata,
				   activepart[j]);

      // Add all active particles contributions to main array
      for (j=0; j<Nactive; j++) {
    	i = activelist[j];
        for (k=0; k<ndim; k++) sphdata[i].a[k] = activepart[j].a[k];
        for (k=0; k<ndim; k++) sphdata[i].agrav[k] = activepart[j].agrav[k];
        for (k=0; k<ndim; k++) sphdata[i].a[k] += sphdata[i].agrav[k];
        sphdata[i].gpot = activepart[j].gpot;
        sphdata[i].dudt = activepart[j].dudt;
        sphdata[i].div_v = activepart[j].div_v;
	sphdata[i].active = false;
        levelneib[i] = max(levelneib[i],activepart[j].levelneib);
      }

    }
    //=========================================================================


    // Finally, add all contributions from distant pair-wise forces to arrays
#pragma omp critical
    for (i=0; i<sph->Nsph; i++) {
      sphdata[i].levelneib = max(sphdata[i].levelneib,levelneib[i]);
    }

    // Free-up local memory for OpenMP thread
    //delete[] neibpart;
    delete[] gravcelllist;
    delete[] directlist;
    delete[] interactlist;
    delete[] neiblist;
    //delete[] activepart;
    //delete[] activelist;
    //delete[] levelneib;

  }
  //===========================================================================

  delete[] celllist;

  //cout << "No. of active ptcls : " << Nactivecount << endl;
  //cout << "Average Nneib       : " << Nneibcount/Nactivecount << endl;
  //cout << "Average Ndirect     : " << Ndirectcount/Nactivecount << endl;
  //cout << "Average Ngravcell   : " << Ncellcount/Nactivecount << endl;

  timing->EndTimingSection("SPH_ALL_FORCES");

  return;
}



//=============================================================================
//  SphTree::UpdateAllSphGravForces
/// Compute all local 'gather' properties of currently active particles, and 
/// then compute each particle's contribution to its (active) neighbour 
/// neighbour hydro forces.  Optimises the algorithm by using grid-cells to 
/// construct local neighbour lists for all particles  inside the cell.
//=============================================================================
template <int ndim, template<int> class ParticleType>
void SphTree<ndim,ParticleType>::UpdateAllSphGravForces
(int Nsph,                          ///< [in] No. of SPH particles
 int Ntot,                          ///< [in] No. of SPH + ghost particles
 SphParticle<ndim> *sph_gen,        ///< [inout] Pointer to SPH ptcl array
 Sph<ndim> *sph,                    ///< [in] Pointer to SPH object
 Nbody<ndim> *nbody)                ///< [in] Pointer to N-body object
{
  int cactive;                      // No. of active cells
  int cc;                           // Aux. cell counter
  int i;                            // Particle id
  int ithread;                      // ..
  int j;                            // Aux. particle counter
  int jj;                           // Aux. particle counter
  int k;                            // Dimension counter
  int okflag;                       // Flag if h-rho iteration is valid
  int Nactive;                      // No. of active particles in cell
  int Ndirect;                      // No. of direct-sum gravity particles
  int Ndirectaux;                   // ..
  int Ndirectmax;                   // Max. no. of direct sum particles
  int Ngravcell;                    // No. of gravity cells
  int Ngravcellmax;                 // Max. no. of gravity cells
  int Ninteract;                    // No. of interactions with hydro neibs
  int Nneib;                        // No. of neighbours
  int Nneibmax;                     // Max. no. of neighbours
  int *activelist;                  // List of active particle ids
  int *directlist;                  // List of direct sum particle ids
  int *interactlist;                // List of interacting neighbour ids
  int *levelneib;                   // ..
  int *neiblist;                    // List of neighbour ids
  FLOAT draux[ndim];                // Aux. relative position vector
  FLOAT drsqd;                      // Distance squared
  FLOAT hrangesqdi;                 // Kernel gather extent
  FLOAT macfactor;                  // Gravity MAC factor
  FLOAT rp[ndim];                   // .. 
  KDTreeCell<ndim> *cell;       // Pointer to binary tree cell
  KDTreeCell<ndim> **celllist;  // List of pointers to binary tree cells
  KDTreeCell<ndim> **gravcelllist; // List of pointers to grav. cells
  ParticleType<ndim> *activepart;    // ..
  ParticleType<ndim> *neibpart;      // ..
  ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sph_gen);

  int Nactivecount = 0;
  int Ncellcount = 0;
  int Ndirectcount = 0;
  int Nneibcount = 0;

  debug2("[SphTree::UpdateAllSphGravForces]");
  timing->StartTimingSection("SPH_GRAV_FORCES",2);

  // Update tree smoothing length values here
  tree->UpdateHmaxValues(tree->kdcell[0],sphdata);


  // Find list of all cells that contain active particles
  celllist = new KDTreeCell<ndim>*[tree->gtot];
  cactive = tree->ComputeActiveCellList(celllist);


  // Set-up all OMP threads
  //===========================================================================
#pragma omp parallel default(none) shared(celllist,cactive,nbody,sph,cout) \
  private(activelist,activepart,cc,cell,directlist,draux,drsqd,gravcelllist)\
  private(hrangesqdi,i,interactlist,ithread,j,jj,k,levelneib,macfactor)	\
  private(neiblist,neibpart,Nactive,Ndirect,Ndirectaux,Ndirectmax,Ngravcell) \
  private(Ngravcellmax,Ninteract,Nneib,Nneibmax,okflag,rp)		\
  reduction(+:Nactivecount,Ncellcount,Ndirectcount,Nneibcount)
  {
#if defined _OPENMP
    ithread = omp_get_thread_num();
#else
    ithread = 0;
#endif
    Nneibmax = Nneibmaxbuf[ithread];
    Ndirectmax = Ndirectmaxbuf[ithread];
    Ngravcellmax = Ngravcellmaxbuf[ithread];

    levelneib = levelneibbuf[ithread];
    activelist = activelistbuf[ithread];
    activepart = activepartbuf[ithread];
    neibpart = neibpartbuf[ithread];

    neiblist = new int[Nneibmax];
    interactlist = new int[Nneibmax];
    directlist = new int[Ndirectmax];
    gravcelllist = new KDTreeCell<ndim>*[Ngravcellmax];


    // Zero timestep level array
    for (i=0; i<sph->Nsph; i++) levelneib[i] = 0.0;


    // Loop over all active cells
    //=========================================================================
#pragma omp for schedule(dynamic)
    for (cc=0; cc<cactive; cc++) {
      cell = celllist[cc];
      macfactor = 0.0;

      // Find list of active particles in current cell
      Nactive = tree->ComputeActiveParticleList(cell,sphdata,activelist);
      Nactivecount += Nactive;

      // Make local copies of active particles
      for (j=0; j<Nactive; j++) activepart[j] = sphdata[activelist[j]];

      // Compute average/maximum term for computing gravity MAC
      if (gravity_mac == "eigenmac") {
	for (j=0; j<Nactive; j++) 
	  macfactor = max(macfactor,pow(1.0/activepart[j].gpot,twothirds));
      }

      // Zero/initialise all summation variables for active particles
      for (j=0; j<Nactive; j++) {
        activepart[j].div_v = (FLOAT) 0.0;
        activepart[j].dudt = (FLOAT) 0.0;
        activepart[j].gpot = activepart[j].m*activepart[j].invh*kernp->wpot(0.0);
        for (k=0; k<ndim; k++) activepart[j].a[k] = (FLOAT) 0.0;
        for (k=0; k<ndim; k++) activepart[j].agrav[k] = (FLOAT) 0.0;
      }

      // Compute neighbour list for cell depending on physics options
      okflag = tree->ComputeGravityInteractionList(cell,macfactor,
					           Nneibmax,Ndirectmax,
                                                   Ngravcellmax,Nneib,Ndirect,
                                                   Ngravcell,neiblist,directlist,
                                                   gravcelllist,sphdata);


      // If there are too many neighbours, reallocate the arrays and
      // recompute the neighbour lists.
      while (okflag == -1) {
        delete[] neibpartbuf[ithread];
        delete[] gravcelllist;
        delete[] directlist;
        delete[] interactlist;
        delete[] neiblist;
        Nneibmax = 2*Nneibmax;
        Ndirectmax = 2*Ndirectmax;
        Ngravcellmax = 2*Ngravcellmax;
	Nneibmaxbuf[ithread] = Nneibmax;
	Ndirectmaxbuf[ithread] = Ndirectmax;
	Ngravcellmaxbuf[ithread] = Ngravcellmax;
        neiblist = new int[Nneibmax];
        interactlist = new int[Nneibmax];
        directlist = new int[Ndirectmax];
        gravcelllist = new KDTreeCell<ndim>*[Ngravcellmax];
        neibpartbuf[ithread] = new ParticleType<ndim>[Nneibmax]; 
	neibpart = neibpartbuf[ithread];
        okflag = tree->ComputeGravityInteractionList(cell,macfactor,
					             Nneibmax,Ndirectmax,
                                                     Ngravcellmax,Nneib,Ndirect,
                                                     Ngravcell,neiblist,directlist,
                                                     gravcelllist,sphdata);
      };


      // Make local copies of all potential neighbours
      for (j=0; j<Nneib; j++) {
        //assert(neiblist[j] >= 0 && neiblist[j] < sph->Ntot);
        neibpart[j] = sphdata[neiblist[j]];
      }

      // Loop over all active particles in the cell
      //-----------------------------------------------------------------------
      for (j=0; j<Nactive; j++) {
        i = activelist[j];

        // Determine SPH neighbour interaction list 
        // (to ensure we don't compute pair-wise forces twice)
        Ninteract = 0;
	Ndirectaux = Ndirect;
        for (k=0; k<ndim; k++) rp[k] = activepart[j].r[k];
        hrangesqdi = activepart[j].hrangesqd;

	//---------------------------------------------------------------------
        for (jj=0; jj<Nneib; jj++) {

          for (k=0; k<ndim; k++) draux[k] = neibpart[jj].r[k] - rp[k];
          drsqd = DotProduct(draux,draux,ndim) + small_number;

          // Compute relative position and distance quantities for pair
          if (drsqd > hrangesqdi && drsqd >= neibpart[jj].hrangesqd)
	    directlist[Ndirectaux++] = neiblist[jj];
	  else if (i != neiblist[jj]) {
            interactlist[Ninteract++] = jj;
	    levelneib[neiblist[jj]] = 
	      max(levelneib[neiblist[jj]],activepart[j].level);
	  }
        }
	//---------------------------------------------------------------------


        // Compute forces between SPH neighbours (hydro and gravity)
        sph->ComputeSphGravForces(i,Ninteract,interactlist,
				  activepart[j],neibpart);

        // Compute direct gravity forces between distant particles
        sph->ComputeDirectGravForces(i,Ndirectaux,directlist,
                                     activepart[j],sphdata);

        // Compute gravitational force due to distant cells
        if (multipole == "monopole")
          tree->ComputeCellMonopoleForces(activepart[j].gpot,activepart[j].agrav,
				    activepart[j].r,Ngravcell,gravcelllist);
        else if (multipole == "quadrupole")
          tree->ComputeCellQuadrupoleForces(activepart[j].gpot,activepart[j].agrav,
				      activepart[j].r,Ngravcell,gravcelllist);

	Nneibcount += Ninteract;
	Ndirectcount += Ndirectaux;
	Ncellcount += Ngravcell;

      }
      //-----------------------------------------------------------------------

      // Compute 'fast' multipole terms here
      if (multipole == "fast_monopole")
	tree->ComputeFastMonopoleForces(Nactive,Ngravcell,gravcelllist,
				  cell,activepart);

      // Compute all star forces
      for (j=0; j<Nactive; j++)
        sph->ComputeStarGravForces(nbody->Nnbody,nbody->nbodydata,
                                   activepart[j]);

      // Add all active particles contributions to main array
      for (j=0; j<Nactive; j++) {
    	i = activelist[j];
        for (k=0; k<ndim; k++) sphdata[i].a[k] += activepart[j].a[k];
        for (k=0; k<ndim; k++) sphdata[i].agrav[k] += activepart[j].agrav[k];
	for (k=0; k<ndim; k++) sphdata[i].a[k] += sphdata[i].agrav[k];
        sphdata[i].gpot = activepart[j].gpot;
        sphdata[i].div_v = activepart[j].div_v;
        sphdata[i].levelneib = activepart[j].levelneib;
	sphdata[i].active = false;
      }

    }
    //=========================================================================


    // Finally, add all contributions from distant pair-wise forces to arrays
#pragma omp critical
    for (i=0; i<sph->Nsph; i++) {
      sphdata[i].levelneib = max(sphdata[i].levelneib,levelneib[i]);
    }


    // Free-up local memory for OpenMP thread
    delete[] gravcelllist;
    delete[] directlist;
    delete[] interactlist;
    delete[] neiblist;
    //delete[] activelist;

  }
  //===========================================================================

  delete[] celllist;

  //cout << "No. of active ptcls : " << Nactivecount << endl;
  //cout << "Average Nneib       : " << Nneibcount/Nactivecount << endl;
  //cout << "Average Ndirect     : " << Ndirectcount/Nactivecount << endl;
  //cout << "Average Ngravcell   : " << Ncellcount/Nactivecount << endl;
  timing->EndTimingSection("SPH_GRAV_FORCES");

  return;
}



//=============================================================================
//  SphTree::UpdateAllStarGasForces
/// Calculate the gravitational acceleration on all star particles due to 
/// all gas particles via the tree.
//=============================================================================
template <int ndim, template<int> class ParticleType>
void SphTree<ndim,ParticleType>::UpdateAllStarGasForces
(int Nsph,                          ///< [in] No. of SPH particles
 int Ntot,                          ///< [in] No. of SPH + ghost particles
 SphParticle<ndim> *sph_gen,        ///< [inout] Pointer to SPH ptcl array
 Sph<ndim> *sph,                    ///< [in] Pointer to SPH object
 Nbody<ndim> *nbody)                ///< [in] Pointer to N-body object
{
  int i;                            // Particle id
  int ithread;                      // OpenMP thread id
  int j;                            // Aux. particle counter
  int jj;                           // Aux. particle counter
  int k;                            // Dimension counter
  int okflag;                       // Flag if h-rho iteration is valid
  int Nactive;                      // No. of active particles in cell
  int Ndirect;                      // No. of direct-sum gravity particles
  int Ndirectaux;                   // ..
  int Ndirectmax;                   // Max. no. of direct sum particles
  int Ngravcell;                    // No. of gravity cells
  int Ngravcellmax;                 // Max. no. of gravity cells
  int Ninteract;                    // No. of interactions with hydro neibs
  int Nneib;                        // No. of neighbours
  int Nneibmax;                     // Max. no. of neighbours
  int *activelist;                  // List of active particle ids
  int *directlist;                  // List of direct sum particle ids
  int *neiblist;                    // List of neighbour ids
  FLOAT macfactor;                  // Gravity MAC factor
  KDTreeCell<ndim> **gravcelllist; // List of pointers to grav. cells
  NbodyParticle<ndim> *star;        // Pointer to star particle
  ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sph_gen);

  int Nactivecount = 0;
  int Ncellcount = 0;
  int Ndirectcount = 0;
  int Nneibcount = 0;

  debug2("[SphTree::UpdateAllStarGasForces]");
  timing->StartTimingSection("STAR_GAS_GRAV_FORCES",2);

  // Make list of all active stars
  Nactive = 0;
  activelist = new int[nbody->Nstar];
  for (i=0; i<nbody->Nstar; i++) 
    if (nbody->nbodydata[i]->active) activelist[Nactive++] = i;


  // Set-up all OMP threads
  //===========================================================================
#pragma omp parallel default(none) shared(activelist,Nactive,nbody,sph,cout)\
  private(directlist,gravcelllist,i,ithread,j,jj,k,macfactor,neiblist) \
  private(Ndirect,Ndirectaux,Ndirectmax,Ngravcell,Ngravcellmax) \
  private(Ninteract,Nneib,Nneibmax,okflag,star)	\
  reduction(+:Nactivecount,Ncellcount,Ndirectcount,Nneibcount)
  {
#if defined _OPENMP
    ithread = omp_get_thread_num();
#else
    ithread = 0;
#endif
    Nneibmax = Ntot; //Nneibmaxbuf[ithread];
    Ndirectmax = Ntot; //Ndirectmaxbuf[ithread];
    Ngravcellmax = Ngravcellmaxbuf[ithread];

    neiblist = new int[Nneibmax];
    directlist = new int[Ndirectmax];
    gravcelllist = new KDTreeCell<ndim>*[Ngravcellmax];


    // Loop over all active cells
    //=========================================================================
#pragma omp for schedule(dynamic)
    for (j=0; j<Nactive; j++) {
      i = activelist[j];
      star = nbody->nbodydata[i];

      // Compute average/maximum term for computing gravity MAC
      if (gravity_mac == "eigenmac")
	macfactor = pow(1.0/star->gpot,twothirds);
      else
	macfactor = 0.0;

      // Compute neighbour list for cell depending on physics options
      okflag = tree->ComputeStarGravityInteractionList(star,macfactor,
						 Nneibmax,Ndirectmax,
                                                 Ngravcellmax,Nneib,Ndirect,
                                                 Ngravcell,neiblist,directlist,
                                                 gravcelllist,sphdata);

      // If there are too many neighbours, reallocate the arrays and
      // recompute the neighbour lists.
      while (okflag == -1) {
        delete[] gravcelllist;
        Ngravcellmax = 2*Ngravcellmax;
        gravcelllist = new KDTreeCell<ndim>*[Ngravcellmax];
        okflag = tree->ComputeStarGravityInteractionList(star,macfactor,Nneibmax,
						   Ndirectmax,Ngravcellmax,
						   Nneib,Ndirect,Ngravcell,
						   neiblist,directlist,
	                                           gravcelllist,sphdata);
      };

      if (Nneib > Nneibmax || Ndirect > Ndirectmax) {
	cout << "Nneib problem : " << Nneib << "   " << Nneibmax << "     "
	     << Ndirect << "    " << Ndirectmax << endl;
	exit(0);
      }

      // Compute contributions to star force from nearby SPH particles
      nbody->CalculateDirectSPHForces(star,Nneib,Ndirect,neiblist,
                                      directlist,sphdata);

      // Compute gravitational force due to distant cells
      if (multipole == "monopole" || multipole == "fast_monopole")
	tree->ComputeCellMonopoleForces(star->gpot,star->a,
				        star->r,Ngravcell,gravcelllist);
      else if (multipole == "quadrupole")
	tree->ComputeCellQuadrupoleForces(star->gpot,star->a,
				          star->r,Ngravcell,gravcelllist);

      Nactivecount += 1;
      Nneibcount += Nneib;
      Ndirectcount += Ndirect;
      Ncellcount += Ngravcell;

    }
    //=========================================================================


    // Free-up local memory for OpenMP thread
    delete[] gravcelllist;
    delete[] directlist;
    delete[] neiblist;

  }
  //===========================================================================

  delete[] activelist;


  //cout << "No. of active ptcls : " << Nactivecount << endl;
  //cout << "Average Nneib2      : " << Nneibcount/Nactivecount << endl;
  //cout << "Average Ndirect2    : " << Ndirectcount/Nactivecount << endl;
  //cout << "Average Ngravcell2  : " << Ncellcount/Nactivecount << endl;
  timing->EndTimingSection("STAR_GAS_GRAV_FORCES");

  return;
}



//=============================================================================
//  SphTree::UpdateAllSphDerivatives
/// ..
//=============================================================================
template <int ndim, template<int> class ParticleType>
void SphTree<ndim,ParticleType>::UpdateAllSphDerivatives
(int Nsph,                          ///< [in] No. of SPH particles
 int Ntot,                          ///< [in] No. of SPH + ghost particles
 SphParticle<ndim> *sph_gen,        ///< [inout] Pointer to SPH ptcl array
 Sph<ndim> *sph)                    ///< [in] Pointer to SPH object
{
  return;
}



//=============================================================================
//  SphTree::UpdateAllSphDudt
/// Compute all local 'gather' properties of currently active particles, and 
/// then compute each particle's contribution to its (active) neighbour 
/// neighbour hydro forces.  Optimises the algorithm by using grid-cells to 
/// construct local neighbour lists for all particles  inside the cell.
//=============================================================================
template <int ndim, template<int> class ParticleType>
void SphTree<ndim,ParticleType>::UpdateAllSphDudt
(int Nsph,                          ///< [in] No. of SPH particles
 int Ntot,                          ///< [in] No. of SPH + ghost particles
 SphParticle<ndim> *sph_gen,        ///< [inout] Pointer to SPH ptcl array
 Sph<ndim> *sph)                    ///< [in] Pointer to SPH object
{
  return;
}


//=============================================================================
//  SphTree::UpdateActiveParticleCounters
/// Loop through all leaf cells in the tree and update all active
/// particle counters.
//=============================================================================
template <int ndim, template<int> class ParticleType>
void SphTree<ndim,ParticleType>::UpdateActiveParticleCounters(SphParticle<ndim> * sphdata_gen,
                                            Sph<ndim> * sph) {

  ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sphdata_gen);
  tree->UpdateActiveParticleCounters(sphdata);
}



#if defined(VERIFY_ALL)
//=============================================================================
//  SphTree::ValidateTree
/// Performs various sanity and validation checks on binary tree structure.
//=============================================================================
template <int ndim, template<int> class ParticleType>
void SphTree<ndim,ParticleType>::ValidateTree
(Sph<ndim> *sph)                    ///< Pointer to SPH class
{
  bool overlap_flag = false;        // Flag if cell bounding boxes overlap
  bool hmax_flag = false;           // Flag if ptcls have larger h than hmax
  bool kill_flag = false;
  int activecount;                  // Active particles in leaf cell
  int c;                            // Cell counter
  int cc;                           // Aux. cell counter
  int i;                            // Particle counter
  int j;                            // Aux. particle counter
  int l;                            // Tree level
  int leafcount;                    // Leaf cell counter
  int Nactivecount=0;               // Counter for total no. of active ptcls
  int Ncount=0;                     // Total particle counter
  int *ccount;                      // Array for counting cells
  int *lcount;                      // Array for counting ptcls on each level
  int *pcount;                      // Array for counting particles in tree
  KDTreeCell<ndim> cell;        // Local copy of binary tree cell

  debug2("[SphTree::ValidateTree]");

  ccount = new int[Ncellmax];
  pcount = new int[Ntot];
  lcount = new int[ltot+1];
  for (i=0; i<Ntot; i++) pcount[i] = 0;
  for (c=0; c<Ncellmax; c++) ccount[c] = 0;
  for (l=0; l<ltot; l++) lcount[l] = 0;
  Ncount = 0;
  Nactivecount = 0;

  // Count how many times we enter a cell in a full tree walk
  c = 0;
  while (c < Ncell) {
    ccount[c]++;
    if (tree[c].c1 != -1) c = tree[c].c1;
    else c = tree[c].cnext;
  }

  // Now check we enter all cells once and once only
  for (c=0; c<Ncell; c++) {
    if (ccount[c] != 1) {
      cout << "Error in cell walk count : " << ccount[c] << endl;
      PrintArray("ccount     : ",Ncell,ccount);
      exit(0);
    }
  }

  // Check inext linked list values and ids array are all valid
  for (i=0; i<sph->Ntot; i++) {
    if (!(ids[i] >= 0 && ids[i] < sph->Ntot)) {
      cout << "Problem with ids array : " 
	   << i << "   " << ids[i] << endl;
      exit(0);
    }
    if (!(inext[i] >= -1 && inext[i] < sph->Ntot)) {
      cout << "Problem with inext linked lists : " 
	   << i << "   " << inext[i] << endl;
      exit(0);
    }
  }


  // Loop over all cells in tree
  //---------------------------------------------------------------------------
  for (c=0; c<Ncell; c++) {
    cell = tree[c];
    activecount = 0;
    leafcount = 0;

    // Add particles from current level
    lcount[cell.level] += cell.N;

    // Check that particles are not in linked lists more than once
    if (cell.level == ltot) {
      i = cell.ifirst;
      while (i != -1) {
	pcount[i]++;
	leafcount++;
	Ncount++;
	if (sphdata[i].active) activecount++;
	if (sphdata[i].active) Nactivecount++;
        if (sphdata[i].h > cell.hmax) {
	  cout << "hmax flag error : " << c << "    " 
	       << sphdata[i].h << "   " << cell.hmax << endl;
	  exit(0);
	}
        if (i == cell.ilast) break;
	i = inext[i];
      }
      if (leafcount > Nleafmax) {
	cout << "Leaf particle count error : " 
	     << leafcount << "   " << Nleafmax << endl;
	exit(0);
      }
      if (activecount > leafcount) {
	cout << "Leaf particle count error : " 
	     << leafcount << "   " << Nleafmax << endl;
	exit(0);
      }
    }

    // Check that bounding boxes of cells on each level do not overlap 
    // each other
    for (cc=0; cc<Ncell; cc++) {
      if (c != cc && tree[cc].level == cell.level) {
	if (ndim == 2) {
	  if (cell.bbmin[0] < tree[cc].bbmax[0] &&
	      cell.bbmax[0] > tree[cc].bbmin[0] &&
	      cell.bbmin[1] < tree[cc].bbmax[1] &&
	      cell.bbmax[1] > tree[cc].bbmin[1])
	    overlap_flag = true;
	}
      }
      if (overlap_flag) {
	cout << "Brother/sister cell overlap error!! : " << c << "   " 
             << cc << endl;
	exit(0);
      }
    }
  }
  //---------------------------------------------------------------------------

  // Check particles are included in the tree once and once only
  for (i=0; i<Ntot; i++) {
    if (pcount[i] != 1) {
      cout << "Problem with child cell ptcl counter : " << i << "   " 
	   << pcount[i] << endl;
      kill_flag = true;
    }
  }

  // Check all particles accounted for
  if (Ncount != sph->Ntot) {
    cout << "Ncount problem with binary tree : " 
	 << Ncount << "   " << sph->Ntot << endl;
    kill_flag = true;
  }

  // Check active particles don't exceed total number of particles
  if (Nactivecount > sph->Nsph) {
    cout << "Nactivecount problem with binary tree : " 
	 << Nactivecount << "   " << sph->Nsph << "   " << sph->Ntot << endl;
    kill_flag = true;
  }

  // Check number of particles on all levels is consistent
  for (l=0; l<ltot; l++) {
    if (lcount[l] != sph->Nsph) {
      cout << "Problem with SPH particles on level : " << l << "    " << lcount[l] << "    " << sph->Nsph << endl;
      kill_flag = true;
    }
  }

  delete[] pcount;
  delete[] ccount;

  if (kill_flag) exit(0);

  return;
}
#endif



template class SphTree<1,GradhSphParticle>;
template class SphTree<2,GradhSphParticle>;
template class SphTree<3,GradhSphParticle>;
template class SphTree<1,SM2012SphParticle>;
template class SphTree<2,SM2012SphParticle>;
template class SphTree<3,SM2012SphParticle>;
template class SphTree<1,GodunovSphParticle>;
template class SphTree<2,GodunovSphParticle>;
template class SphTree<3,GodunovSphParticle>;
