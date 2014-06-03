//=============================================================================
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
//  GradhSphTree::GradhSphTree
/// GradhSphTree constructor.  Initialises various variables.
//=============================================================================
template <int ndim, template<int> class ParticleType>
GradhSphTree<ndim,ParticleType>::GradhSphTree
(int Nleafmaxaux,
 FLOAT thetamaxsqdaux,
 FLOAT kernrangeaux,
 FLOAT macerroraux,
 string gravity_mac_aux,
 string multipole_aux,
 DomainBox<ndim> *boxaux,
 SphKernel<ndim> *kernaux,
 CodeTiming *timingaux):
  SphTree<ndim,ParticleType>(Nleafmaxaux,thetamaxsqdaux,kernrangeaux,
                             macerroraux,gravity_mac_aux,multipole_aux,
                             boxaux,kernaux,timingaux)
{
}



//=============================================================================
//  GradhSphTree::~GradhSphTree
/// GradhSphTree destructor.  Deallocates tree memory upon object destruction.
//=============================================================================
template <int ndim, template<int> class ParticleType>
GradhSphTree<ndim,ParticleType>::~GradhSphTree()
{
  if (tree->allocated_tree) {
    this->DeallocateMemory();
    tree->DeallocateTreeMemory();
  }
}



//=============================================================================
//  GradhSphTree::UpdateAllSphProperties
/// Compute all local 'gather' properties of currently active particles, and 
/// then compute each particle's contribution to its (active) neighbour 
/// neighbour hydro forces.  Optimises the algorithm by using grid-cells to 
/// construct local neighbour lists for all particles  inside the cell.
//=============================================================================
template <int ndim, template<int> class ParticleType>
void GradhSphTree<ndim,ParticleType>::UpdateAllSphProperties
(int Nsph,                          ///< [in] No. of SPH particles
 int Ntot,                          ///< [in] No. of SPH + ghost particles
 SphParticle<ndim> *sph_gen,        ///< [inout] Pointer to SPH ptcl array
 Sph<ndim> *sph,                    ///< [in] Pointer to SPH object
 Nbody<ndim> *nbody)                ///< [in] Pointer to N-body object
{
  int celldone;                     // ..
  int okflag;                       // ..
  int cc;                           // Aux. cell counter
  int cactive;                      // No. of active
  int i;                            // Particle id
  int ithread;                      // OpenMP thread i.d.
  int j;                            // Aux. particle counter
  int jj;                           // Aux. particle counter
  int k;                            // Dimension counter
  int Nactive;                      // No. of active particles in cell
  int Ngather;                      // No. of near gather neighbours
  int Nneib;                        // No. of neighbours
  int Nneibmax;                     // Max. no. of neighbours
  int *activelist;                  // List of active particle ids
  int *neiblist;                    // List of neighbour ids
  FLOAT draux[ndim];                // Aux. relative position vector var
  FLOAT drsqdaux;                   // Distance squared
  FLOAT hrangesqd;                  // Kernel extent
  FLOAT hmax;                       // Maximum smoothing length
  FLOAT rp[ndim];                   // Local copy of particle position
  FLOAT *drsqd;                     // Position vectors to gather neibs
  FLOAT *gpot;                      // Potential for particles
  FLOAT *gpot2;                     // ..
  FLOAT *m;                         // Distances to gather neibs
  FLOAT *m2;                        // ..
  FLOAT *mu;                        // mass*u for gather neibs
  FLOAT *mu2;                       // ..
  FLOAT *r;                         // Positions of neibs
  KDTreeCell<ndim> *cell;           // Pointer to binary tree cell
  KDTreeCell<ndim> **celllist;      // List of binary cell pointers
  ParticleType<ndim> *activepart;   // ..
  ParticleType<ndim> *sphdata = static_cast<ParticleType<ndim>* > (sph_gen);

  debug2("[GradhSphTree::UpdateAllSphProperties]");
  timing->StartTimingSection("SPH_PROPERTIES",2);

  // Find list of all cells that contain active particles
  celllist = new KDTreeCell<ndim>*[tree->gtot];
  cactive = tree->ComputeActiveCellList(celllist);


  // Set-up all OMP threads
  //===========================================================================
#pragma omp parallel default(none) private(activelist,activepart,cc,cell)\
  private(celldone,draux,drsqd,drsqdaux,gpot,gpot2,hmax,hrangesqd,i,ithread,j)\
  private(jj,k,okflag,m,mu,m2,mu2,Nactive,neiblist,Ngather,Nneib,Nneibmax,r,rp)	\
  shared(cactive,celllist,cout,nbody,sph,sphdata)
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
          delete[] m;
          delete[] drsqd;
          delete[] gpot2;
          delete[] gpot;
          delete[] neiblist;
          Nneibmax = 2*Nneibmax; 
          neiblist = new int[Nneibmax];
          gpot = new FLOAT[Nneibmax];
          gpot2 = new FLOAT[Nneibmax];
          drsqd = new FLOAT[Nneibmax];
          m = new FLOAT[Nneibmax];
          m2 = new FLOAT[Nneibmax];
          r = new FLOAT[Nneibmax*ndim];
          Nneib = tree->ComputeGatherNeighbourList(cell,Nneibmax,neiblist,
                                                   hmax,sphdata);
        };


        // Make local copies of important neib information (mass and position)
        for (jj=0; jj<Nneib; jj++) {
          j = neiblist[jj];
          gpot[jj] = sphdata[j].gpot;
          m[jj] = sphdata[j].m;
          for (k=0; k<ndim; k++) r[ndim*jj + k] = (FLOAT) sphdata[j].r[k];
        }


        // Loop over all active particles in the cell
        //---------------------------------------------------------------------
        for (j=0; j<Nactive; j++) {
          i = activelist[j];
          for (k=0; k<ndim; k++) rp[k] = activepart[j].r[k];

          // Set gather range as current h multiplied by some tolerance factor
          hrangesqd = kernrangesqd*hmax*hmax;
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
              Ngather++;
            }

          }
          //-------------------------------------------------------------------

          // Validate that gather neighbour list is correct
#if defined(VERIFY_ALL)
          if (neibcheck) 
	    this->CheckValidNeighbourList(i,Ntot,Nneib,neiblist,
                                          sphdata,"gather");
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
//  GradhSphTree::UpdateAllSphHydroForces
/// Compute all local 'gather' properties of currently active particles, and 
/// then compute each particle's contribution to its (active) neighbour 
/// neighbour hydro forces.  Optimises the algorithm by using grid-cells to 
/// construct local neighbour lists for all particles  inside the cell.
//=============================================================================
template <int ndim, template<int> class ParticleType>
void GradhSphTree<ndim,ParticleType>::UpdateAllSphHydroForces
(int Nsph,                          ///< [in] No. of SPH particles
 int Ntot,                          ///< [in] No. of SPH + ghost particles
 SphParticle<ndim> *sph_gen,        ///< [inout] Pointer to SPH ptcl array
 Sph<ndim> *sph,                    ///< [in] Pointer to SPH object
 Nbody<ndim> *nbody)                ///< [in] Pointer to N-body object
{
  int cactive;                      // No. of active cells
  int cc;                           // Aux. cell counter
  int i;                            // Particle id
  int ithread;                      // OpenMP thread i.d.
  int j;                            // Aux. particle counter
  int jj;                           // Aux. particle counter
  int k;                            // Dimension counter
  int Nactive;                      // No. of active particles in cell
  int Ninteract;                    // No. of near gather neighbours
  int Nneib;                        // No. of neighbours
  int Nneibmax;                     // Max. no. of neighbours
  int *activelist;                  // List of active particle ids
  int *interactlist;                // List of interactng SPH neighbours
  int *levelneib;                   // ..
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
  ParticleType<ndim> *activepart;   // ..
  ParticleType<ndim> *neibpart;     // ..
  ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sph_gen);

  debug2("[GradhSphTree::UpdateAllSphHydroForces]");
  timing->StartTimingSection("SPH_HYDRO_FORCES",2);


  // Find list of all cells that contain active particles
  celllist = new KDTreeCell<ndim>*[tree->gtot];
  cactive = tree->ComputeActiveCellList(celllist);

  // If there are no active cells, return to main loop
  if (cactive == 0) return;

  // Update tree smoothing length values here
  tree->UpdateHmaxValues(tree->kdcell[0],sphdata);


  // Set-up all OMP threads
  //===========================================================================
#pragma omp parallel default(none) private(activelist,activepart,cc,cell,dr)\
  private(draux,drmag,drsqd,hrangesqdi,i,interactlist,ithread,invdrmag,j,jj,k)\
  private(levelneib,Nactive,neiblist,neibpart,Ninteract,Nneib,Nneibmax,rp) \
  shared(cactive,celllist,sphdata,nbody,sph)
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
        activepart[j].dalphadt = (FLOAT) 0.0;
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
	  this->CheckValidNeighbourList(i,Ntot,Nneib,neiblist,
					sphdata,"all");
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
        sphdata[i].dalphadt = activepart[j].dalphadt;
        sphdata[i].div_v = activepart[j].div_v;
	sphdata[i].active = false;
        levelneib[i] = max(levelneib[i],activepart[j].levelneib);
      }

    }
    //=========================================================================


    // Finally, add all contributions from distant pair-wise forces to arrays
#pragma omp critical
    for (i=0; i<sph->Nsph; i++)
      sphdata[i].levelneib = max(sphdata[i].levelneib,levelneib[i]);


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
//  GradhSphTree::UpdateAllSphForces
/// Compute all local 'gather' properties of currently active particles, and 
/// then compute each particle's contribution to its (active) neighbour 
/// neighbour hydro forces.  Optimises the algorithm by using grid-cells to 
/// construct local neighbour lists for all particles  inside the cell.
//=============================================================================
template <int ndim, template<int> class ParticleType>
void GradhSphTree<ndim,ParticleType>::UpdateAllSphForces
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
  KDTreeCell<ndim> *cell;           // Pointer to kd-tree cell
  KDTreeCell<ndim> **celllist;      // List of pointers to kd-tree cells
  KDTreeCell<ndim> **gravcelllist;  // List of pointers to grav. cells
  ParticleType<ndim> *activepart;   // ..
  ParticleType<ndim> *neibpart;     // ..
  ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sph_gen);

  debug2("[GradhSphTree::UpdateAllSphForces]");
  timing->StartTimingSection("SPH_ALL_FORCES",2);

  // Update tree smoothing length values here
  tree->UpdateHmaxValues(tree->kdcell[0],sphdata);

  // Find list of all cells that contain active particles
  celllist = new KDTreeCell<ndim>*[tree->gtot];
  cactive = tree->ComputeActiveCellList(celllist);


  // Set-up all OMP threads
  //===========================================================================
#pragma omp parallel default(none) shared(celllist,cactive,nbody,sph,sphdata,cout) \
  private(activepart,activelist,cc,cell,directlist,draux,drsqd,gravcelllist) \
  private(hrangesqdi,i,interactlist,ithread,j,jj,k,levelneib,macfactor)	\
  private(neiblist,neibpart,Nactive,Ndirect,Ndirectaux,Ndirectmax) \
  private(Ngravcell,Ngravcellmax,Ninteract,Nneib,Nneibmax,okflag,rp)
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
	for (k=0; k<ndim; k++) activepart[j].adot[k] = (FLOAT) 0.0;
      }


      // Compute neighbour list for cell depending on physics options
      okflag = tree->ComputeGravityInteractionList(cell,macfactor,
					           Nneibmax,Ndirectmax,
                                                   Ngravcellmax,Nneib,Ndirect,
                                                   Ngravcell,neiblist,directlist,
                                                   gravcelllist,sphdata);

      // If there are too many neighbours, reallocate the arrays and
      // recompute the neighbour lists.
      while (okflag < 0 || Ndirect + Nneib > Ndirectmax) {
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
      for (j=0; j<Nneib; j++) neibpart[j] = sphdata[neiblist[j]];


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

  timing->EndTimingSection("SPH_ALL_FORCES");

  return;
}



//=============================================================================
//  GradhSphTree::UpdateAllSphGravForces
/// Compute all local 'gather' properties of currently active particles, and 
/// then compute each particle's contribution to its (active) neighbour 
/// neighbour hydro forces.  Optimises the algorithm by using grid-cells to 
/// construct local neighbour lists for all particles  inside the cell.
//=============================================================================
template <int ndim, template<int> class ParticleType>
void GradhSphTree<ndim,ParticleType>::UpdateAllSphGravForces
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
  KDTreeCell<ndim> *cell;           // Pointer to binary tree cell
  KDTreeCell<ndim> **celllist;      // List of pointers to binary tree cells
  KDTreeCell<ndim> **gravcelllist;  // List of pointers to grav. cells
  ParticleType<ndim> *activepart;   // ..
  ParticleType<ndim> *neibpart;     // ..
  ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sph_gen);

  debug2("[GradhSphTree::UpdateAllSphGravForces]");
  timing->StartTimingSection("SPH_GRAV_FORCES",2);

  // Update tree smoothing length values here
  tree->UpdateHmaxValues(tree->kdcell[0],sphdata);

  // Find list of all cells that contain active particles
  celllist = new KDTreeCell<ndim>*[tree->gtot];
  cactive = tree->ComputeActiveCellList(celllist);


  // Set-up all OMP threads
  //===========================================================================
#pragma omp parallel default(none) shared(celllist,cactive,nbody,sph,sphdata,cout) \
  private(activepart,activelist,cc,cell,directlist,draux,drsqd,gravcelllist) \
  private(hrangesqdi,i,interactlist,ithread,j,jj,k,levelneib,macfactor)	\
  private(neiblist,neibpart,Nactive,Ndirect,Ndirectaux,Ndirectmax) \
  private(Ngravcell,Ngravcellmax,Ninteract,Nneib,Nneibmax,okflag,rp)
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
	for (k=0; k<ndim; k++) activepart[j].adot[k] = (FLOAT) 0.0;
      }


      // Compute neighbour list for cell depending on physics options
      okflag = tree->ComputeGravityInteractionList(cell,macfactor,
					           Nneibmax,Ndirectmax,
                                                   Ngravcellmax,Nneib,Ndirect,
                                                   Ngravcell,neiblist,directlist,
                                                   gravcelllist,sphdata);

      // If there are too many neighbours, reallocate the arrays and
      // recompute the neighbour lists.
      while (okflag < 0 || Ndirect + Nneib > Ndirectmax) {
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
        sph->ComputeSphGravForces(i,Ninteract,interactlist,
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

  timing->EndTimingSection("SPH_GRAV_FORCES");

  return;
}



template class GradhSphTree<1,GradhSphParticle>;
template class GradhSphTree<2,GradhSphParticle>;
template class GradhSphTree<3,GradhSphParticle>;
