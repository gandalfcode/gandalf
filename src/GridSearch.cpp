//=============================================================================
//  GridSearch.cpp
//  Contains functions for grid neighbour search routines.
//  Creates a uniform grid from particle distribution where the spacing is 
//  the size of the maximum kernel range (i.e. kernrange*h_max) over all ptcls.
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
using namespace std;


static const FLOAT grid_h_tolerance = (FLOAT) 1.1;


//=============================================================================
//  GridSearch::GridSearch
/// GridSearch constructor.  Initialises various variables.
//=============================================================================
template <int ndim>
GridSearch<ndim>::GridSearch()
{
  allocated_grid = false;
  Ncell = 0;
  Ncellmax = 0;
  Noccupymax = 0;
  Ntot = 0;
  Ntotmax = 0;
}



//=============================================================================
//  GridSearch::~GridSearch
/// GridSearch destructor.  Deallocates grid memory upon object destruction.
//=============================================================================
template <int ndim>
GridSearch<ndim>::~GridSearch()
{
  if (allocated_grid) DeallocateGridMemory();
}



//=============================================================================
//  GridSearch::BuildTree
/// Creates a new grid structure each time the neighbour 'tree' needs to be 
/// updated.
//=============================================================================
template <int ndim>
void GridSearch<ndim>::BuildTree(Sph<ndim> *sph, Parameters &simparams)
{
  CreateGrid(sph);
  return;
}



//=============================================================================
//  GridSearch::UpdateTree
/// Creates a new grid structure each time the neighbour 'tree' needs to be 
/// updated.
//=============================================================================
template <int ndim>
void GridSearch<ndim>::UpdateTree(Sph<ndim> *sph, Parameters &simparams)
{
  CreateGrid(sph);
  return;
}



//=============================================================================
//  GridSearch::UpdateActiveParticleCounters
/// ..
//=============================================================================
template <int ndim>
void GridSearch<ndim>::UpdateActiveParticleCounters(Sph<ndim> *sph)
{
  int c;
  int i;
  int ilast;

  // Loop through all grid cells in turn
  //---------------------------------------------------------------------------
  for (c=0; c<Ncell; c++) {
    grid[c].Nactive = 0;

    if (grid[c].Nptcls == 0) continue;
    i = grid[c].ifirst;
    ilast = grid[c].ilast;

    // Else walk through linked list to obtain list and number of active ptcls.
    do {
      if (i < sph->Nsph && sph->sphdata[i].active) grid[c].Nactive++;
      if (i == ilast) break;
      i = inext[i];
    } while (i != -1);

  }
  //---------------------------------------------------------------------------


  return;
}



//=============================================================================
//  GridSearch::UpdateAllSphProperties
/// Compute all local 'gather' properties of currently active particles, and 
/// then compute each particle's contribution to its (active) neighbour 
/// neighbour hydro forces.  Optimises the algorithm by using grid-cells to 
/// construct local neighbour lists for all particles  inside the cell.
//=============================================================================
template <int ndim>
void GridSearch<ndim>::UpdateAllSphProperties
(Sph<ndim> *sph,                    ///< [inout] Pointer to main SPH object
 Nbody<ndim> *nbody)                ///< [in] Pointer to main N-body object
{
  int c;                            // Cell id
  int cactive;                      // No. of active cells
  int cc;                           // Aux. cell counter
  int i;                            // Particle id
  int j;                            // Aux. particle counter
  int jj;                           // Aux. particle counter
  int k;                            // Dimension counter
  int okflag;                       // Flag if h-rho iteration is valid
  int Nactive;                      // No. of active particles in cell
  int Ngather;                      // No. of near gather neighbours
  int Nneib;                        // No. of neighbours
  int Nneibmax;                     // Max. no. of neighbours
  int *activelist;                  // List of active particle ids
  int *celllist;                    // List of active cells
  int *gatherlist;                  // List of nearby neighbour ids
  int *neiblist;                    // List of neighbour ids
  FLOAT draux[ndim];                // Aux. relative position vector var
  FLOAT drsqdaux;                   // Distance squared
  FLOAT hrangesqd;                  // Kernel extent
  FLOAT rp[ndim];                   // Local copy of particle position
  FLOAT *drsqd;                     // Position vectors to gather neibs
  FLOAT *gpot;                      // Potential for particles
  FLOAT *gpot2;                     // ..
  FLOAT *m;                         // Masses of potential neighbours
  FLOAT *m2;                        // Reduced list of masses
  FLOAT *mu;                        // mass*u for potential gather neibs
  FLOAT *mu2;                       // Reduced list of mu
  FLOAT *r;                         // 1/drmag to scatter neibs
  SphParticle<ndim> *data = sph->sphdata;   // Pointer to SPH particle data

  debug2("[GridSearch::UpdateAllSphProperties]");

  // Find list of all cells that contain active particles
  celllist = new int[Ncell];
  cactive = ComputeActiveCellList(celllist);
  Nneibmax = Nlistmax;

  // Set-up all OMP threads
  //===========================================================================
#pragma omp parallel default(none) private(activelist,c,cc,draux,drsqd) \
  private(drsqdaux,gatherlist,hrangesqd,i,j,jj,k,okflag,m,m2,mu,mu2,Nactive) \
  private(neiblist,Ngather,Nneib,r,rp,gpot,gpot2) \
  shared(sph,data,nbody,Nneibmax,cactive,celllist)
  {
    activelist = new int[Noccupymax];
    gatherlist = new int[Nneibmax];
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
#pragma omp for schedule(dynamic)
    for (cc=0; cc<cactive; cc++) {
      c = celllist[cc];

      // Find list of active particles in current cell
      Nactive = ComputeActiveParticleList(c,activelist,sph);

      // Compute neighbour list for cell depending on physics options
      Nneib = ComputeNeighbourList(c,neiblist);

      // Make local copies of important neib information (mass and position)
      for (jj=0; jj<Nneib; jj++) {
        j = neiblist[jj];
        gpot[jj] = data[j].gpot;
        m[jj] = data[j].m;
        mu[jj] = data[j].m*data[j].u;
        for (k=0; k<ndim; k++) r[ndim*jj + k] = (FLOAT) data[j].r[k];
      }

      // Loop over all active particles in the cell
      //-----------------------------------------------------------------------
      for (j=0; j<Nactive; j++) {
        i = activelist[j];
        for (k=0; k<ndim; k++) rp[k] = data[i].r[k];

        // Set gather range as current h multiplied by some tolerance factor
        hrangesqd = grid_h_tolerance*grid_h_tolerance*
          sph->kernp->kernrangesqd*data[i].h*data[i].h;
        Ngather = 0;

        // Compute distance (squared) to all
        //---------------------------------------------------------------------
        for (jj=0; jj<Nneib; jj++) {
          for (k=0; k<ndim; k++) draux[k] = r[ndim*jj + k] - rp[k];
          drsqdaux = DotProduct(draux,draux,ndim);

          // Record distance squared for all potential gather neighbours
          if (drsqdaux < hrangesqd) {
            drsqd[Ngather] = drsqdaux;
            gpot2[Ngather] = gpot[jj];
            m2[Ngather] = m[jj];
            mu2[Ngather] = mu[jj];
            gatherlist[Ngather] = jj;
            Ngather++;
       	  }
	  
        }
        //---------------------------------------------------------------------

        // Validate that gather neighbour list is correct
#if defined(VERIFY_ALL)
        if (neibcheck) CheckValidNeighbourList(sph,i,Ngather,
                                               gatherlist,"gather");
#endif

        // Compute smoothing length and other gather properties for particle i
        okflag = sph->ComputeH(i,Ngather,big_number,m2,mu2,
                               drsqd,gpot,data[i],nbody);

      }
      //-----------------------------------------------------------------------

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
    delete[] gatherlist;
    delete[] activelist;

  }
  //===========================================================================

  delete[] celllist;

  return;
}



//=============================================================================
//  GridSearch::UpdateAllSphHydroForces
/// Compute all local 'gather' properties of currently active particles, and 
/// then compute each particle's contribution to its (active) neighbour 
/// neighbour hydro forces.  Optimises the algorithm by using grid-cells to 
/// construct local neighbour lists for all particles  inside the cell.
//=============================================================================
template <int ndim>
void GridSearch<ndim>::UpdateAllSphHydroForces
(Sph<ndim> *sph)                    ///< Pointer to SPH object
{
  int c;                            // Cell id
  int cactive;                      // No. of active cells
  int cc;                           // Aux. cell counter
  int i;                            // Particle id
  int j;                            // Aux. particle counter
  int jj;                           // Aux. particle counter
  int k;                            // Dimension counter
  int okflag;                       // Flag if h-rho iteration is valid
  int Nactive;                      // No. of active particles in cell
  int Ngrav;                        // No. of direct sum gravity ptcls
  int Ninteract;                    // No. of interaction pairs
  int Nneib;                        // No. of neighbours
  int Nneibmax;                     // Max. no. of neighbours
  int *activelist;                  // List of active particle ids
  int *celllist;                    // List of active cells
  int *interactlist;                // List of interaction pairs to compute
  int *neiblist;                    // List of neighbour ids
  FLOAT draux[ndim];                // Aux. relative position vector var
  FLOAT drsqd;                      // Distance squared
  FLOAT hrangesqdi;                 // Kernel extent for particle i
  //FLOAT hrangesqdj;                 // Kernek extent for particle j
  FLOAT rp[ndim];                   // Local copy of particle position
  FLOAT *dr;                        // Array of relative position vectors
  FLOAT *drmag;                     // Array of neighbour distances
  FLOAT *invdrmag;                  // Array of 1/drmag between particles
  SphParticle<ndim> *activepart;    // Local copyh of active particles
  SphParticle<ndim> *neibpart;      // Local copy of neighbouring ptcls
  SphParticle<ndim> *data = sph->sphdata;   // Pointer to SPH particle data

  debug2("[GridSearch::UpdateAllSphHydroForces]");

  // Find list of all cells that contain active particles
  celllist = new int[Ncell];
  cactive = ComputeActiveCellList(celllist);


  // Set-up all OMP threads
  //===========================================================================
#pragma omp parallel default(none) private(activelist,activepart,c,cc,dr,draux,drmag,drsqd) \
  private(hrangesqdi,hrangesqdj,i,interactlist,invdrmag,j,jj,k,okflag,Nactive) \
  private(neiblist,neibpart,Ninteract,Nneib,Nneibmax,rp) \
  shared(cactive,celllist,data,sph)
  {
    Nneibmax = Nlistmax;
    activelist = new int[Noccupymax];
    neiblist = new int[Nneibmax];
    interactlist = new int[Nneibmax];
    dr = new FLOAT[Nneibmax*ndim];
    drmag = new FLOAT[Nneibmax];
    invdrmag = new FLOAT[Nneibmax];
    activepart = new SphParticle<ndim>[Noccupymax];
    neibpart = new SphParticle<ndim>[Nneibmax];

    // Loop over all active cells
    //=========================================================================
#pragma omp for schedule(dynamic)
    for (cc=0; cc<cactive; cc++) {
      c = celllist[cc];

      // Find list of active particles in current cell
      Nactive = ComputeActiveParticleList(c,activelist,sph);

      // Compute neighbour list for cell depending on physics options
      Nneib = ComputeNeighbourList(c,neiblist);

      // Validate that gather neighbour list is correct
#if defined(VERIFY_ALL)
      if (neibcheck) CheckValidNeighbourList(sph,i,Nneib,neiblist,"gather");
#endif

      // Make local copies of active particles
      for (j=0; j<Nactive; j++) {
        activepart[j] = data[activelist[j]];
      }

      // Make local copies of all potential neighbours
      for (j=0; j<Nneib; j++) {
        neibpart[j] = data[neiblist[j]];
        neibpart[j].div_v = (FLOAT) 0.0;
        neibpart[j].dudt = (FLOAT) 0.0;
        neibpart[j].levelneib = 0;
        for (k=0; k<ndim; k++) neibpart[j].a[k] = (FLOAT) 0.0;
      }

      // Loop over all active particles in the cell
      //-----------------------------------------------------------------------
      for (j=0; j<Nactive; j++) {
        i = activelist[j];
        activepart[j].div_v = (FLOAT) 0.0;
        activepart[j].dudt = (FLOAT) 0.0;
        activepart[j].levelneib = 0;
        for (k=0; k<ndim; k++) activepart[j].a[k] = (FLOAT) 0.0;
        for (k=0; k<ndim; k++) rp[k] = activepart[j].r[k];
        hrangesqdi = activepart[j].hrangesqd;
        //hrangesqdi = sph->kernfacsqd*sph->kernp->kernrangesqd*
	//activepart[j].h*activepart[j].h;
        Ninteract = 0;

        // Compute distances and the inverse between the current particle
        // and all neighbours here, for both gather and inactive scatter neibs.
        // Only consider particles with j > i to compute pair forces once
        // unless particle j is inactive.
        //---------------------------------------------------------------------
        for (jj=0; jj<Nneib; jj++) {

  	      // Skip neighbour if it's not the correct part of an active pair
          if (neiblist[jj] <= i && neibpart[jj].active) continue;

          //hrangesqdj = sph->kernfacsqd*sph->kernp->kernrangesqd*neibpart[jj].h*neibpart[jj].h;
          for (k=0; k<ndim; k++) draux[k] = neibpart[jj].r[k] - rp[k];
          drsqd = DotProduct(draux,draux,ndim);

          // Compute list of particle-neighbour interactions and also
          // compute
          if ((drsqd <= hrangesqdi || drsqd <= neibpart[jj].hrangesqd)) { //&&
              //((neiblist[jj] < i && !neibpart[jj].active) ||
               //neiblist[jj] > i)) {
            drmag[Ninteract] = sqrt(drsqd) + small_number;
            invdrmag[Ninteract] = (FLOAT) 1.0/drmag[Ninteract];
            for (k=0; k<ndim; k++)
              dr[Ninteract*ndim + k] = draux[k]*invdrmag[Ninteract];
            interactlist[Ninteract] = jj;
            Ninteract++;
          }
	  
        }
        //---------------------------------------------------------------------

        // Compute all gather neighbour contributions to hydro forces
        sph->ComputeSphHydroForces(i,Ninteract,interactlist,
				   drmag,invdrmag,dr,activepart[j],neibpart);

      }
      //-----------------------------------------------------------------------


    // Add all particle i contributions to main array
//#pragma omp critical
	{
      for (jj=0; jj<Nactive; jj++) {
        j = activelist[jj];
	    for (k=0; k<ndim; k++) {
#pragma omp atomic
	      data[j].a[k] += activepart[jj].a[k];
	    }
#pragma omp atomic
	    data[j].dudt += activepart[jj].dudt;
#pragma omp atomic
	    data[j].div_v += activepart[jj].div_v;
	    //data[j].levelneib = max(data[i].levelneib,activepart[jj].levelneib);
	  }
	}


      // Now add all active neighbour contributions to the main arrays
//#pragma omp critical
      {
	for (jj=0; jj<Nneib; jj++) {
	  j = neiblist[jj];
	  if (neibpart[jj].active) {
	    for (k=0; k<ndim; k++) {
#pragma omp atomic
	      data[j].a[k] += neibpart[jj].a[k];
	    }
#pragma omp atomic
	    data[j].dudt += neibpart[jj].dudt;
#pragma omp atomic
	    data[j].div_v += neibpart[jj].div_v;
	  }
	  //data[j].levelneib = max(data[j].levelneib,neibpart[jj].levelneib);
	}
      }
      
    }
    //=========================================================================

    // Free-up local memory for OpenMP thread
    delete[] neibpart;
    delete[] activepart;
    delete[] invdrmag;
    delete[] drmag;
    delete[] dr;
    delete[] interactlist;
    delete[] neiblist;
    delete[] activelist;
    
  }
  //===========================================================================

  delete[] celllist;


  // Compute other important SPH quantities after hydro forces are computed
  if (sph->hydro_forces == 1) {
    for (i=0; i<sph->Nsph; i++) {
      if (sph->sphdata[i].active)
    	  sph->ComputePostHydroQuantities(sph->sphdata[i]);
    }
  }

  return;
}



//=============================================================================
//  GridSearch::UpdateAllSphGravityProperties
/// Empty function for now
//=============================================================================
template <int ndim>
void GridSearch<ndim>::UpdateAllSphForces
(Sph<ndim> *sph)                    ///< Pointer to SPH object
{
  return;
}



//=============================================================================
//  GridSearch::UpdateAllSphGravityProperties
/// Empty function for now.
//=============================================================================
template <int ndim>
void GridSearch<ndim>::UpdateAllSphGravForces
(Sph<ndim> *sph)                    ///< Pointer to SPH object
{
  return;
}



//=============================================================================
//  GridSearch::UpdateAllSphDerivatives
/// ..
//=============================================================================
template <int ndim>
void GridSearch<ndim>::UpdateAllSphDerivatives(Sph<ndim> *sph)
{
  int c;                            // Cell id
  int cactive;                      // No. of active cells
  int cc;                           // Aux. cell counter
  int i;                            // Particle id
  int j;                            // Aux. particle counter
  int jj;                           // Aux. particle counter
  int k;                            // Dimension counter
  int okflag;                       // Flag if h-rho iteration is valid
  int Nactive;                      // No. of active particles in cell
  int Ninteract;                    // No. of pair interactions
  int Nneib;                        // No. of neighbours
  int Nneibmax;                     // Max. no. of neighbours
  int *activelist;                  // List of active particle ids
  int *celllist;                    // List of active cells
  int *interactlist;                // List of pair interaction neighbours
  int *neiblist;                    // List of neighbour ids
  FLOAT draux[ndim];                // Aux. relative position vector var
  FLOAT drsqd;                      // Distance squared
  FLOAT hrangesqd;                  // Kernel extent squared
  FLOAT rp[ndim];                   // Local copy of particle position
  FLOAT *dr;                        // Array of relative position vectors
  FLOAT *drmag;                     // Array of neighbour distances
  FLOAT *invdrmag;                  // Array of 1/drmag between particles
  SphParticle<ndim> *neibpart;      // Local copy of neighbouring ptcls
  SphParticle<ndim> parti;                  // Local copy of SPH particle
  SphParticle<ndim> *data = sph->sphdata;   // Pointer to SPH particle data

  debug2("[GridSearch::UpdateAllSphProperties]");

  // Find list of all cells that contain active particles
  celllist = new int[Ncell];
  cactive = ComputeActiveCellList(celllist);
  Nneibmax = Nlistmax;


  // Set-up all OMP threads
  //===========================================================================
#pragma omp parallel default(none) private(activelist,c,cc,dr,draux,drmag,drsqd)\
  private(hrangesqd,i,interactlist,invdrmag,j,jj,k,okflag,Nactive,neiblist)\
  private(neibpart,Ninteract,Nneib,parti,rp) shared(celllist, cactive, Nneibmax, data)\
  shared(sph)
  {
    activelist = new int[Noccupymax];
    neiblist = new int[Nneibmax];
    interactlist = new int[Nneibmax];
    dr = new FLOAT[Nneibmax*ndim];
    drmag = new FLOAT[Nneibmax];
    invdrmag = new FLOAT[Nneibmax];
    neibpart = new SphParticle<ndim>[Nneibmax];

    // Loop over all active cells
    //=========================================================================
#pragma omp for schedule(dynamic)
    for (cc=0; cc<cactive; cc++) {
      c = celllist[cc];

      // Find list of active particles in current cell
      Nactive = ComputeActiveParticleList(c,activelist,sph);

      // Compute neighbour list for cell depending on physics options
      Nneib = ComputeNeighbourList(c,neiblist);
      for (j=0; j<Nneib; j++) neibpart[j] = data[neiblist[j]];

      // Loop over all active particles in the cell
      //-----------------------------------------------------------------------
      for (j=0; j<Nactive; j++) {
        i = activelist[j];
        parti = data[i];

        for (k=0; k<ndim; k++) rp[k] = parti.r[k]; //data[i].r[k];
        hrangesqd = pow(sph->kernp->kernrange*parti.h,2);
        Ninteract = 0;

        // Compute distances and the inverse between the current particle
        // and all neighbours here, for both gather and inactive scatter neibs.
        // Only consider particles with j > i to compute pair forces once
        // unless particle j is inactive.
        //---------------------------------------------------------------------
        for (jj=0; jj<Nneib; jj++) {

          for (k=0; k<ndim; k++) draux[k] = neibpart[jj].r[k] - rp[k];
          drsqd = DotProduct(draux,draux,ndim);

          // Compute list of particle-neighbour interactions
          if (drsqd <= hrangesqd) {
	    interactlist[Ninteract] = jj;
	    drmag[Ninteract] = sqrt(drsqd);
	    invdrmag[Ninteract] = (FLOAT) 1.0/
	      (drmag[Ninteract] + small_number);
	    for (k=0; k<ndim; k++)
	      dr[Ninteract*ndim + k] = draux[k]*invdrmag[Ninteract];
	    Ninteract++;
          }
	  
        }
        //---------------------------------------------------------------------

        // Compute all gather neighbour contributions to hydro forces
        sph->ComputeSphDerivatives(i,Ninteract,interactlist,
				   drmag,invdrmag,dr,parti,neibpart);

	data[i] = parti;

      }
      //-----------------------------------------------------------------------

    }
    //=========================================================================

    // Free-up local memory for OpenMP thread
    delete[] neibpart;
    delete[] invdrmag;
    delete[] drmag;
    delete[] dr;
    delete[] interactlist;
    delete[] neiblist;
    delete[] activelist;
    
  }
  //===========================================================================

  delete[] celllist;

  return;
}




//=============================================================================
//  GridSearch::UpdateAllSphDudt
/// Compute all local 'gather' properties of currently active particles, and 
/// then compute each particle's contribution to its (active) neighbour 
/// neighbour hydro forces.  Optimises the algorithm by using grid-cells to 
/// construct local neighbour lists for all particles  inside the cell.
//=============================================================================
template <int ndim>
void GridSearch<ndim>::UpdateAllSphDudt(Sph<ndim> *sph)
{
  int c;                            // Cell id
  int cactive;                      // No. of active cells
  int cc;                           // Aux. cell counter
  int i;                            // Particle id
  int j;                            // Aux. particle counter
  int jj;                           // Aux. particle counter
  int k;                            // Dimension counter
  int okflag;                       // Flag if h-rho iteration is valid
  int Nactive;                      // No. of active particles in cell
  int Ngrav;                        // No. of direct sum gravity ptcls
  int Ninteract;                    // No. of near gather neighbours
  int Nneib;                        // No. of neighbours
  int Nneibmax;                     // Max. no. of neighbours
  int *activelist;                  // List of active particle ids
  int *celllist;                    // List of active cells
  int *interactlist;                // List of pair interaction ids
  int *neiblist;                    // List of neighbour ids
  FLOAT draux[ndim];                // Aux. relative position vector var
  FLOAT drsqd;                      // Distance squared
  FLOAT hrangesqdi;                 // Kernel extent of particle i squared
  FLOAT hrangesqdj;                 // Kernel extent of particle j squared
  FLOAT rp[ndim];                   // Local copy of particle position
  FLOAT *dr;                        // Array of relative position vectors
  FLOAT *drmag;                     // Array of neighbour distances
  FLOAT *invdrmag;                  // Array of 1/drmag between particles
  SphParticle<ndim> parti;          // Local copy of SPH particle
  SphParticle<ndim> *neibpart;      // Local copy of neighbouring ptcls
  SphParticle<ndim> *data = sph->sphdata;   // Pointer to SPH particle data

  debug2("[GridSearch::UpdateAllSphDudt]");

  // Find list of all cells that contain active particles
  celllist = new int[Ncell];
  cactive = ComputeActiveCellList(celllist);
  Nneibmax = Nlistmax;


  // Set-up all OMP threads
  //===========================================================================
#pragma omp parallel default(none) private(activelist,c,cc,dr,draux,drmag,drsqd)\
  private(hrangesqdi,hrangesqdj,i,interactlist,invdrmag,j,jj,k,okflag,Nactive) \
  private(neiblist,neibpart,Ninteract,Nneib,parti,rp) shared(sph, data, celllist)\
  shared(Nneibmax, cactive)
  {
    activelist = new int[Noccupymax];
    neiblist = new int[Nneibmax];
    interactlist = new int[Nneibmax];
    dr = new FLOAT[Nneibmax*ndim];
    drmag = new FLOAT[Nneibmax];
    invdrmag = new FLOAT[Nneibmax];
    neibpart = new SphParticle<ndim>[Nneibmax];

    // Loop over all active cells
    //=========================================================================
#pragma omp for schedule(dynamic)
    for (cc=0; cc<cactive; cc++) {
      c = celllist[cc];

      // Find list of active particles in current cell
      Nactive = ComputeActiveParticleList(c,activelist,sph);

      // Compute neighbour list for cell depending on physics options
      Nneib = ComputeNeighbourList(c,neiblist);

      // Validate that gather neighbour list is correct
#if defined(VERIFY_ALL)
      if (neibcheck) CheckValidNeighbourList(sph,i,Nneib,neiblist,"gather");
#endif

      // Make local copies of all potential neighbours
      for (j=0; j<Nneib; j++) {
        neibpart[j] = data[neiblist[j]];
        neibpart[j].dudt = (FLOAT) 0.0;
      }

      // Loop over all active particles in the cell
      //-----------------------------------------------------------------------
      for (j=0; j<Nactive; j++) {
        i = activelist[j];
        parti = data[i];
        parti.dudt = (FLOAT) 0.0;

        for (k=0; k<ndim; k++) rp[k] = parti.r[k];
        hrangesqdi = pow(sph->kernfac*sph->kernp->kernrange*parti.h,2);
        Ninteract = 0;

        // Compute distances and the inverse between the current particle
        // and all neighbours here, for both gather and inactive scatter neibs.
        // Only consider particles with j > i to compute pair forces once
        // unless particle j is inactive.
        //---------------------------------------------------------------------
        for (jj=0; jj<Nneib; jj++) {

          hrangesqdj = 
            pow(sph->kernfac*sph->kernp->kernrange*neibpart[jj].h,2);
          for (k=0; k<ndim; k++) draux[k] = neibpart[jj].r[k] - rp[k];
          drsqd = DotProduct(draux,draux,ndim);

          // Compute list of particle-neighbour pair interactions
          if ((drsqd <= hrangesqdi || drsqd <= hrangesqdj) &&
	      ((neiblist[jj] < i && !neibpart[jj].active) ||
	       neiblist[jj] > i)) {
	    interactlist[Ninteract] = jj;
	    drmag[Ninteract] = sqrt(drsqd);
	    invdrmag[Ninteract] = (FLOAT) 1.0/
	      (drmag[Ninteract] + small_number);
	    for (k=0; k<ndim; k++)
	      dr[Ninteract*ndim + k] = draux[k]*invdrmag[Ninteract];
	    Ninteract++;
          }
	  
        }
        //---------------------------------------------------------------------

        // Compute all gather neighbour contributions to hydro forces
        sph->ComputeSphNeibDudt(i,Ninteract,interactlist,
				drmag,invdrmag,dr,parti,neibpart);

        // Add partial sum for particle i to main array
#pragma omp atomic
        data[i].dudt += parti.dudt;
	
      }
      //-----------------------------------------------------------------------

      // Now add all active neighbour contributions to the main arrays
      for (jj=0; jj<Nneib; jj++) {
        if (neibpart[jj].active) {
          j = neiblist[jj];
#pragma omp atomic
          data[j].dudt += neibpart[jj].dudt;
        }
      }
      
    }
    //=========================================================================

    // Free-up local memory for OpenMP thread
    delete[] neibpart;
    delete[] invdrmag;
    delete[] drmag;
    delete[] dr;
    delete[] interactlist;
    delete[] neiblist;
    delete[] activelist;
    
  }
  //===========================================================================

  delete[] celllist;

  return;
}



//=============================================================================
//  GridSearch::AllocateGridMemory
/// Allocate memory for neighbour grid as requested.  If more memory is 
/// required than currently allocated, grid is deallocated and reallocated.
//=============================================================================
template <int ndim>
void GridSearch<ndim>::AllocateGridMemory(int Npart)
{
  debug2("[GridSearch::AllocateGridMemory]");

  Ntot = Npart;

  if (Ntot > Ntotmax || Ncell > Ncellmax) {
    if (allocated_grid) DeallocateGridMemory();
    Ntotmax = 3*Ntot;
    Ncellmax = 3*Ncell;
    inext = new int[Ntotmax];
    grid = new struct GridCell[Ncellmax];
  }

  return;
}



//=============================================================================
//  GridSearch::DeallocateGridMemory
/// Deallocates all neighbour grid memory
//=============================================================================
template <int ndim>
void GridSearch<ndim>::DeallocateGridMemory(void)
{
  debug2("[GridSearch::DeallocateGridMemory]");

  delete[] grid;
  delete[] inext;
  allocated_grid = false;

  return;
}



//=============================================================================
//  GridSearch::CreateGrid
/// Create a regular neighbour grid using all SPH particles contained within 
/// the SPH object.  The grid spacing is equal to the maximum smoothing kernel 
/// range of all particles multiplied by some arbitrary tolerance parameter 
/// (grid_h_tolerance) to allow for some smoothing lengths increasing.
//=============================================================================
template <int ndim>
void GridSearch<ndim>::CreateGrid(Sph<ndim> *sph)
{
  int i;                            // Particle counter
  int k;                            // Dimension counter
  int c;                            // Grid cell counter/id
  FLOAT h_max = 0.0;                // Maximum smoothing length of ptcls

  debug2("[GridSearch::CreateGrid]");
  
  // Compute maximum smoothing length to determine optimum grid spacing
  for (i=0; i<sph->Nsph; i++) h_max = max(h_max,sph->sphdata[i].h);
  dx_grid = grid_h_tolerance*sph->kernfac*sph->kernp->kernrange*h_max;

  // Compute bounding box of all particles
  sph->SphBoundingBox(rmax,rmin,sph->Ntot);

  // Calculate no. of grid cells in each dimension and in total
  Ncell = 1;
  for (k=0; k<ndim; k++) {
    Ngrid[k] = (int)((rmax[k] - rmin[k])/dx_grid) + 1;
    Ncell = Ngrid[k]*Ncell;
  }

  // Allocate memory for grid if not previously done
  AllocateGridMemory(sph->Ntot);
  Nsph = sph->Nsph;

  // Initialise all values in cells
  for (c=0; c<Ncellmax; c++) {
    grid[c].Nactive = 0;
    grid[c].Nptcls = 0;
    grid[c].ifirst = 0;
    grid[c].ilast = 0;
  }
  for (i=0; i<Ntotmax; i++) inext[i] = -1;

  // Now attach all particles to grid cells
  //---------------------------------------------------------------------------
  for (i=0; i<sph->Ntot; i++) {
    c = ComputeParticleGridCell(sph->sphdata[i].r);

    // If cell currently contains no particles, record first particle.
    // Else, add to end of linked list.
    if (c < 0) cout << "rp : " << sph->sphdata[i].r[0] << "   " 
		    << sph->sphdata[i].a[0] << "   " << sph->sphdata[i].dudt 
		    << endl;
    if (grid[c].Nptcls == 0) grid[c].ifirst = i;
    else inext[grid[c].ilast] = i;
    grid[c].ilast = i;
    grid[c].Nptcls++;
    if (i < sph->Nsph && sph->sphdata[i].active) grid[c].Nactive++;

  }
  //---------------------------------------------------------------------------

  // Find maximum occupations of all cells
  Noccupymax = 0;
  for (c=0; c<Ncell; c++) Noccupymax = max(Noccupymax,grid[c].Nptcls);
  Nlistmax = Noccupymax*pow(3,ndim);

#if defined(VERIFY_ALL)
  ValidateGrid();
#endif

  return;
}



//=============================================================================
//  GridSearch::ComputeParticleGridCell
/// Compute and return the grid cell i.d. that contains the position 'rp'.
//=============================================================================
template <int ndim>
int GridSearch<ndim>::ComputeParticleGridCell
(FLOAT *rp)                         ///< [in] Position vector
{
  int k;                            // Dimension counter
  int igrid[ndim];                  // Grid cell coordinate

  for (k=0; k<ndim; k++) {
    igrid[k] = (int) ((rp[k] - rmin[k])/dx_grid);
    if (igrid[k] < 0) igrid[k] = 0;
    if (igrid[k] >= Ngrid[k]) igrid[k] = Ngrid[k] - 1;
  }
  if (ndim == 1) 
    return igrid[0];
  else if (ndim == 2) 
    return igrid[0] + Ngrid[0]*igrid[1];
  else if (ndim == 3) 
    return igrid[0] + Ngrid[0]*igrid[1] + Ngrid[0]*Ngrid[1]*igrid[2];
}



//=============================================================================
//  GridSearch::ComputeCellCoordinate
/// Computes and returns the grid cell coordinate 'igrid' from the 
/// grid cell i.d. 'c'
//=============================================================================
template <int ndim>
void GridSearch<ndim>::ComputeCellCoordinate
(int c,                             ///< [in] Grid-cell id
 int igrid[ndim])                   ///< [out] Grid-cell co-ordinate
{
  if (ndim == 1) 
    igrid[0] = c;
  else if (ndim == 2) {
    igrid[1] = c/Ngrid[0];
    igrid[0] = c%Ngrid[0];
  }
  else if (ndim == 3) {
    igrid[2] = c/Ngrid[0]/Ngrid[1];
    igrid[1] = (c/Ngrid[0])%Ngrid[1];
    igrid[0] = c%Ngrid[0];
  }
#if defined(VERIFY_ALL)
  if ((ndim == 1 && c != igrid[0]) ||
      (ndim == 2 && c != igrid[0] + Ngrid[0]*igrid[1]) || 
      (ndim == 3 && c!= igrid[0] + Ngrid[0]*igrid[1] + 
       Ngrid[0]*Ngrid[1]*igrid[2])) {
    string message= "Problem computing Grid cell coorindate";
    ExceptionHandler::getIstance().raise(message);
  }
#endif
  return;
}



//=============================================================================
//  GridSearch::ComputeActiveCellList
/// Returns the number of cells containing active particles, 'Nactive', and 
/// the i.d. list of cells contains active particles, 'celllist'
//=============================================================================
template <int ndim>
int GridSearch<ndim>::ComputeActiveCellList
(int *celllist)                     ///< [out] List of cell ids
{
  int c;                            // Cell counter
  int Nactive = 0;                  // No. of cells containing active ptcls

  debug2("[GridSearch::ComputeActiveCellList]");

  for (c=0; c<Ncell; c++)
    if (grid[c].Nactive > 0) celllist[Nactive++] = c;

  return Nactive;
}



//=============================================================================
//  GridSearch::ComputeActiveParticleList
/// Returns the number (Nactive) and list of ids (activelist) of all active 
/// SPH particles in the given cell 'c'.
//=============================================================================
template <int ndim>
int GridSearch<ndim>::ComputeActiveParticleList
(int c,                             ///< [in] Cell i.d.
 int *activelist,                   ///< [out] List of active particles in cell
 Sph<ndim> *sph)                    ///< [in] SPH object pointer
{
  int Nactive = 0;                  // No. of active particles in cell
  int i = grid[c].ifirst;           // Particle id (set to first ptcl id)
  int ilast = grid[c].ilast;        // i.d. of last particle in cell c

  // If there are no active particles in this cell, return without walking list
  if (grid[c].Nptcls == 0) return 0;

  // Else walk through linked list to obtain list and number of active ptcls.
  do {
    if (i < sph->Nsph && sph->sphdata[i].active) activelist[Nactive++] = i;
    if (i == ilast) break;
    i = inext[i];
  } while (i != -1);

  return Nactive;
}



//=============================================================================
//  GridSearch::ComputeNeighbourList
/// Computes and returns number of neighbour, 'Nneib', and the list 
/// of neighbour ids, 'neiblist', for all particles inside cell 'c'.
/// Includes all particles in the selected cell, plus all particles 
/// contained in adjacent cells (including diagonal cells).
//=============================================================================
template <int ndim>
int GridSearch<ndim>::ComputeNeighbourList
(int c,                             ///< [in] i.d. of cell
 int *neiblist)                     ///< [out] List of neighbour i.d.s
{
  int i;                            // Particle id
  int ilast;                        // id of last particle in current cell
  int caux,cx,cy,cz;                // Aux. cell counters and coordinates
  int igrid[ndim];                  // Grid cell coordinate
  int gridmin[ndim];                // Minimum neighbour cell coordinate
  int gridmax[ndim];                // Maximum neighbour cell coordinate
  int Nneib = 0;                    // No. of neighbours

  // Compute the location of the cell on the grid using the id
  ComputeCellCoordinate(c,igrid);

  //---------------------------------------------------------------------------
  if (ndim == 1) {
    gridmin[0] = max(0,igrid[0]-1);
    gridmax[0] = min(Ngrid[0]-1,igrid[0]+1);
    
    for (cx=gridmin[0]; cx<=gridmax[0]; cx++) {
      caux = cx;
      if (grid[caux].Nptcls == 0) continue;
      i = grid[caux].ifirst;
      ilast = grid[caux].ilast;
      do {
	neiblist[Nneib++] = i;
	if (i == ilast) break;
	i = inext[i];
      } while (i != -1);
    }
  }
  //---------------------------------------------------------------------------
  else if (ndim == 2) {
    gridmin[0] = max(0,igrid[0]-1);
    gridmax[0] = min(Ngrid[0]-1,igrid[0]+1);
    gridmin[1] = max(0,igrid[1]-1);
    gridmax[1] = min(Ngrid[1]-1,igrid[1]+1);
    
    for (cy=gridmin[1]; cy<=gridmax[1]; cy++) {
      for (cx=gridmin[0]; cx<=gridmax[0]; cx++) {
	caux = cx + cy*Ngrid[0];
	if (grid[caux].Nptcls == 0) continue;
	i = grid[caux].ifirst;
	ilast = grid[caux].ilast;
	do {
	  neiblist[Nneib++] = i;
	  if (i == ilast) break;
	  i = inext[i];
	} while (i != -1);
      }
    }
  }
  //---------------------------------------------------------------------------
  else if (ndim == 3) {
    gridmin[0] = max(0,igrid[0]-1);
    gridmax[0] = min(Ngrid[0]-1,igrid[0]+1);
    gridmin[1] = max(0,igrid[1]-1);
    gridmax[1] = min(Ngrid[1]-1,igrid[1]+1);
    gridmin[2] = max(0,igrid[2]-1);
    gridmax[2] = min(Ngrid[2]-1,igrid[2]+1);
    
    for (cz=gridmin[2]; cz<=gridmax[2]; cz++) {
      for (cy=gridmin[1]; cy<=gridmax[1]; cy++) {
	for (cx=gridmin[0]; cx<=gridmax[0]; cx++) {
	  caux = cx + cy*Ngrid[0] + cz*Ngrid[0]*Ngrid[1];
	  if (grid[caux].Nptcls == 0) continue;
	  i = grid[caux].ifirst;
	  ilast = grid[caux].ilast;
	  do {
	    neiblist[Nneib++] = i;
	    if (i == ilast) break;
	    i = inext[i];
	  } while (i != -1);
	}
      }
    }
  }
  //---------------------------------------------------------------------------

  return Nneib;
}



#if defined(VERIFY_ALL)
//=============================================================================
//  GridSearch::CheckValidNeighbourList
/// Checks that the neighbour list generated by the grid is valid in that it 
/// (i) does include all true neighbours, and 
/// (ii) all true neigbours are only included once and once only.
//=============================================================================
template <int ndim>
void GridSearch<ndim>::CheckValidNeighbourList
(Sph<ndim> *sph,                    ///< [in] SPH object pointer
 int i,                             ///< [in] Particle i.d.
 int Nneib,                         ///< [in] No. of potential neighbours
 int *neiblist,                     ///< [in] List of potential neighbour i.d.s
 string neibtype)                   ///< [in] Neighbour search type
{
  int count;                        // Valid neighbour counter
  int j;                            // Neighbour particle counter
  int k;                            // Dimension counter
  int Ntrueneib = 0;                // No. of 'true' neighbours
  int *trueneiblist;                // List of true neighbour ids
  FLOAT drsqd;                      // Distance squared
  FLOAT dr[ndim];                   // Relative position vector

  // Allocate array to store local copy of potential neighbour ids
  trueneiblist = new int[sph->Ntot];

  // First, create list of 'true' neighbours by looping over all particles
  if (neibtype == "gather") {
    for (j=0; j<sph->Ntot; j++) {
      for (k=0; k<ndim; k++)
	dr[k] = sph->sphdata[j].r[k] - sph->sphdata[i].r[k];
      drsqd = DotProduct(dr,dr,ndim);
      if (drsqd <= 
	  sph->kernp->kernrangesqd*sph->sphdata[i].h*sph->sphdata[i].h)
	trueneiblist[Ntrueneib++] = j;
    }
  }

  // Now compare each given neighbour with true neighbour list for validation
  for (j=0; j<Ntrueneib; j++) {
    count = 0;
    for (k=0; k<Nneib; k++) if (neiblist[k] == trueneiblist[j]) count++;

    // If the true neighbour is not in the list, or included multiple times, 
    // then output to screen and terminate program
    if (count != 1) {
      cout << "Problem with neighbour lists : " << i << "  " << j << "   "
	   << count << "   "
	   << sph->sphdata[i].r[0] << "   " << sph->sphdata[i].h << endl;
      cout << "Nneib : " << Nneib << "   Ntrueneib : " << Ntrueneib << endl;
      PrintArray("neiblist     : ",Nneib,neiblist);
      PrintArray("trueneiblist : ",Ntrueneib,trueneiblist);
      string message = "Problem with neighbour lists in grid search";
      ExceptionHandler::getIstance().raise(message);
    }
  }

  delete[] trueneiblist;

  return;
}



//=============================================================================
//  GridSearch::ValidateGrid
/// Validate that the grid structure is consistent with several sanity checks.
//=============================================================================
template <int ndim>
void GridSearch<ndim>::ValidateGrid(void)
{
  int c;                            // Cell counter
  int i;                            // Particle id
  int ilast;                        // i.d. of last particle in list
  int *gridentry;                   // No. of time ptcl is present in grid

  debug2("[GridSearch::ValidateGrid]");

  gridentry = new int[Ntot];
  for (i=0; i<Ntot; i++) gridentry[i] = 0;

  // Loop over all grid cells and count how many times particles appear 
  // in linked lists
  for (c=0; c<Ncell; c++) {
    if (grid[c].Nptcls == 0) continue;
    i = grid[c].ifirst;
    ilast = grid[c].ilast;
    do {
      gridentry[i]++;
      cout << "COUNTING : " << c << "   " << i << "   " << gridentry[i] 
	   << "    " << grid[c].Nptcls << "   " << ilast << endl;
      if (i == ilast) break;
      i = inext[i];
    } while (i != -1);
  }

  // If particles appear multiple times, or not at all, quit with error msg
  for (i=0; i<Ntot; i++) {
    if (gridentry[i] != 1) {
      cout << "Problem with particle in gridentry : " 
	   << i << "   " << gridentry[i] << endl;
      string message = "Problem with particle in gridentry";
      ExceptionHandler::getIstance().raise(message);
    }
  }

  // Check that linked list termination flag (i.e. inext = -1) is consistent
  // If not, quit with error message to screen
  for (c=0; c<Ncell; c++) {
    if (grid[c].Nptcls != 0 && inext[grid[c].ilast] != -1) {
      cout << "Error in linked list : " << c << "  " << grid[c].ilast << "  "
	   << inext[grid[c].ilast] << endl;
      string message = "Error in grid search linked lists";
      ExceptionHandler::getIstance().raise(message);
    }
  }

  delete[] gridentry;

  return;
}
#endif



template class GridSearch<1>;
template class GridSearch<2>;
template class GridSearch<3>;
