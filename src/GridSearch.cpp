// ============================================================================
// GridSearch.cpp
// Contains functions for grid neighbour search routines.
// Creates a uniform grid from particle distribution where the spacing is 
// the size of the maximum kernel range (i.e. kernrange*h_max) over all ptcls.
// ============================================================================


#include <cstdlib>
#include <iostream>
#include <string>
#include <math.h>
#include "Precision.h"
#include "Dimensions.h"
#include "Exception.h"
#include "SphNeighbourSearch.h"
#include "Sph.h"
#include "Parameters.h"
#include "InlineFuncs.h"
#include "SphParticle.h"
#include "Debug.h"
using namespace std;


static const FLOAT grid_h_tolerance = (FLOAT) 1.3;


// ============================================================================
// GridSearch::GridSearch
// ============================================================================
GridSearch::GridSearch(int ndimaux)
{
#if !defined(FIXED_DIMENSIONS)
  ndim = ndimaux;
#endif
  allocated_grid = false;
  Ncell = 0;
  Ncellmax = 0;
  Noccupymax = 0;
  Ntot = 0;
  Ntotmax = 0;
}



// ============================================================================
// GridSearch::~GridSearch
// ============================================================================
GridSearch::~GridSearch()
{
  if (allocated_grid) DeallocateGridMemory();
}



// ============================================================================
// GridSearch::UpdateTree
// ============================================================================
void GridSearch::UpdateTree(Sph *sph, Parameters &simparams)
{
  CreateGrid(sph);
  return;
}



// ============================================================================
// GridSearch::UpdateAllSphProperties
// Compute all local 'gather' properties of currently active particles, and 
// then compute each particle's contribution to its (active) neighbour 
// neighbour hydro forces.  Optimises the algorithm by using grid-cells to 
// construct local neighbour lists for all particles  inside the cell.
// ============================================================================
void GridSearch::UpdateAllSphProperties(Sph *sph)
{
  int c;                                // Cell id
  int cactive;                          // No. of active cells
  int cc;                               // Aux. cell counter
  int i;                                // Particle id
  int j;                                // Aux. particle counter
  int jj;                               // Aux. particle counter
  int k;                                // Dimension counter
  int okflag;                           // Flag if h-rho iteration is valid
  int Nactive;                          // No. of active particles in cell
  int Ngather;                          // No. of near gather neighbours
  int Nneib;                            // No. of neighbours
  int Nneibmax;                         // Max. no. of neighbours
  int *activelist;                      // List of active particle ids
  int *celllist;                        // List of active cells
  int *neiblist;                        // List of neighbour ids
  int *gatherlist;                      // List of nearby neighbour ids
  FLOAT draux[ndimmax];                 // Aux. relative position vector var
  FLOAT drsqdaux;                       // Distance squared
  FLOAT hrangesqd;                      // Kernel extent
  FLOAT rp[ndimmax];                    // Local copy of particle position
  FLOAT *drsqd;                         // Position vectors to gather neibs
  FLOAT *m;                             // Distances to gather neibs
  FLOAT *m2;                            // ..
  FLOAT *mu;                            // mass*u for gather neibs
  FLOAT *mu2;                           // ..
  FLOAT *r;                             // 1/drmag to scatter neibs
  SphParticle *data = sph->sphdata;     // Pointer to SPH particle data

  debug2("[GridSearch::UpdateAllSphProperties]");

  // Find list of all cells that contain active particles
  celllist = new int[Ncell];
  cactive = ComputeActiveCellList(celllist);
  Nneibmax = Nlistmax;

  // Set-up all OMP threads
  // ==========================================================================
#pragma omp parallel default(shared) private(activelist,c,cc,draux,drsqd,drsqdaux,gatherlist,hrangesqd,i,j,jj,k,okflag,m,m2,mu,mu2,Nactive,neiblist,Ngather,Nneib,r,rp)
  {
    activelist = new int[Noccupymax];
    gatherlist = new int[Nneibmax];
    neiblist = new int[Nneibmax];
    drsqd = new FLOAT[Nneibmax];
    m = new FLOAT[Nneibmax];
    m2 = new FLOAT[Nneibmax];
    mu = new FLOAT[Nneibmax];
    mu2 = new FLOAT[Nneibmax];
    r = new FLOAT[Nneibmax*ndim];

    // Loop over all active cells
    // ========================================================================
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
        m[jj] = data[j].m;
        mu[jj] = data[j].m*data[j].u;
        for (k=0; k<ndim; k++) r[ndim*jj + k] = (FLOAT) data[j].r[k];
      }

      // Loop over all active particles in the cell
      // ----------------------------------------------------------------------
      for (j=0; j<Nactive; j++) {
        i = activelist[j];
        for (k=0; k<ndim; k++) rp[k] = data[i].r[k];

        // Set gather range as current h multiplied by some tolerance factor
        hrangesqd = pow(grid_h_tolerance*sph->kernp->kernrange*data[i].h,2);
        Ngather = 0;

        // Compute distance (squared) to all
        // --------------------------------------------------------------------
        for (jj=0; jj<Nneib; jj++) {
          for (k=0; k<ndim; k++) draux[k] = r[ndim*jj + k] - rp[k];
          drsqdaux = DotProduct(draux,draux,ndim);

          // Record distance squared for all potential gather neighbours
          if (drsqdaux <= hrangesqd) {
            gatherlist[Ngather] = jj;
            drsqd[Ngather] = drsqdaux;
            m2[Ngather] = m[jj];
            mu2[Ngather] = mu[jj];
            Ngather++;
       	  }
	  
        }
        // --------------------------------------------------------------------

        // Compute smoothing length and other gather properties for particle i
        okflag = sph->ComputeH(i,Ngather,m2,mu2,drsqd,data[i]);

      }
      // ----------------------------------------------------------------------

    }
    // ========================================================================

    // Free-up all memory
    delete[] r;
    delete[] mu2;
    delete[] mu;
    delete[] m2;
    delete[] m;
    delete[] drsqd;
    delete[] neiblist;
    delete[] gatherlist;
    delete[] activelist;

  }
  // ==========================================================================

  delete[] celllist;

  return;
}



// ============================================================================
// GridSearch::UpdateAllSphForces
// Compute all local 'gather' properties of currently active particles, and 
// then compute each particle's contribution to its (active) neighbour 
// neighbour hydro forces.  Optimises the algorithm by using grid-cells to 
// construct local neighbour lists for all particles  inside the cell.
// ============================================================================
void GridSearch::UpdateAllSphForces(Sph *sph)
{
  int c;                                  // Cell id
  int cactive;                            // No. of active cells
  int cc;                                 // Aux. cell counter
  int i;                                  // Particle id
  int j;                                  // Aux. particle counter
  int jj;                                 // Aux. particle counter
  int k;                                  // Dimension counter
  int okflag;                             // Flag if h-rho iteration is valid
  int Nactive;                            // No. of active particles in cell
  int Ngrav;                              // No. of direct sum gravity ptcls
  int Ninteract;                          // No. of near gather neighbours
  int Nneib;                              // No. of neighbours
  int Nneibmax;                           // Max. no. of neighbours
  int *activelist;                        // List of active particle ids
  int *celllist;                          // List of active cells
  int *interactlist;                      // ..
  int *neiblist;                          // List of neighbour ids
  FLOAT draux[ndimmax];                   // Aux. relative position vector var
  FLOAT drsqd;                            // Distance squared
  FLOAT hrangesqdi;                       // Kernel extent
  FLOAT hrangesqdj;                       // ..
  FLOAT rp[ndimmax];                      // Local copy of particle position
  FLOAT *dr;                              // Array of relative position vectors
  FLOAT *drmag;                           // Array of neighbour distances
  FLOAT *invdrmag;                        // Array of 1/drmag between particles
  SphParticle *neibpart;                  // Local copy of neighbouring ptcls
  SphParticle parti;                      // Local copy of SPH particle
  SphParticle *data = sph->sphdata;       // Pointer to SPH particle data

  debug2("[GridSearch::UpdateAllSphProperties]");

  // Find list of all cells that contain active particles
  celllist = new int[Ncell];
  cactive = ComputeActiveCellList(celllist);
  Nneibmax = Nlistmax;


  // Set-up all OMP threads
  // ==========================================================================
#pragma omp parallel default(shared) private(activelist,c,cc,dr,draux,drmag,drsqd,hrangesqdi,hrangesqdj,i,interactlist,invdrmag,j,jj,k,okflag,Nactive,neiblist,neibpart,Ninteract,Nneib,parti,rp)
  {
    activelist = new int[Noccupymax];
    neiblist = new int[Nneibmax];
    interactlist = new int[Nneibmax];
    dr = new FLOAT[Nneibmax*ndim];
    drmag = new FLOAT[Nneibmax];
    invdrmag = new FLOAT[Nneibmax];
    neibpart = new SphParticle[Nneibmax];

    // Loop over all active cells
    // ========================================================================
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
        neibpart[j].div_v = (FLOAT) 0.0;
        neibpart[j].dudt = (FLOAT) 0.0;
        for (k=0; k<ndim; k++) neibpart[j].a[k] = (FLOAT) 0.0;
      }

      // Loop over all active particles in the cell
      // ----------------------------------------------------------------------
      for (j=0; j<Nactive; j++) {
        i = activelist[j];
        parti = data[i];
        parti.div_v = (FLOAT) 0.0;
        parti.dudt = (FLOAT) 0.0;
        for (k=0; k<ndim; k++) parti.a[k] = (FLOAT) 0.0;

        for (k=0; k<ndim; k++) rp[k] = parti.r[k]; //data[i].r[k];
        hrangesqdi = pow(sph->kernfac*sph->kernp->kernrange*parti.h,2);
        Ninteract = 0;

        // Compute distances and the inverse between the current particle
        // and all neighbours here, for both gather and inactive scatter neibs.
        // Only consider particles with j > i to compute pair forces once
        // unless particle j is inactive.
        // --------------------------------------------------------------------
        for (jj=0; jj<Nneib; jj++) {

          hrangesqdj = pow(sph->kernfac*sph->kernp->kernrange*neibpart[jj].h,2);
          for (k=0; k<ndim; k++) draux[k] = neibpart[jj].r[k] - rp[k];
          drsqd = DotProduct(draux,draux,ndim);

          // Compute list of particle-neighbour interactions and also
          // compute
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
        // --------------------------------------------------------------------

        // Compute all gather neighbour contributions to hydro forces
        sph->ComputeSphNeibForces(i,Ninteract,interactlist,
				  drmag,invdrmag,dr,parti,neibpart);

        // Add ..
        for (k=0; k<ndim; k++) {
#pragma omp atomic
          data[i].a[k] += parti.a[k];
        }
#pragma omp atomic
        data[i].dudt += parti.dudt;
#pragma omp atomic
        data[i].div_v += parti.div_v;
	
      }
      // ----------------------------------------------------------------------

      // Now add all active neighbour contributions to the main arrays
      for (jj=0; jj<Nneib; jj++) {
        if (neibpart[jj].active) {
          j = neiblist[jj];
          for (k=0; k<ndim; k++) {
#pragma omp atomic
            data[j].a[k] += neibpart[jj].a[k];
          }
#pragma omp atomic
          data[j].dudt += neibpart[jj].dudt;
#pragma omp atomic
          data[j].div_v += neibpart[jj].div_v;
        }
      }
      
    }
    // ========================================================================

    // Free-up local memory for OpenMP thread
    delete[] neibpart;
    delete[] invdrmag;
    delete[] drmag;
    delete[] dr;
    delete[] interactlist;
    delete[] neiblist;
    delete[] activelist;
    
  }
  // ==========================================================================

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



// ============================================================================
// GridSearch::UpdateAllSphDerivatives
// ..
// ============================================================================
void GridSearch::UpdateAllSphDerivatives(Sph *sph)
{
  int c;                                  // Cell id
  int cactive;                            // No. of active cells
  int cc;                                 // Aux. cell counter
  int i;                                  // Particle id
  int j;                                  // Aux. particle counter
  int jj;                                 // Aux. particle counter
  int k;                                  // Dimension counter
  int okflag;                             // Flag if h-rho iteration is valid
  int Nactive;                            // No. of active particles in cell
  int Ninteract;                          // No. of near gather neighbours
  int Nneib;                              // No. of neighbours
  int Nneibmax;                           // Max. no. of neighbours
  int *activelist;                        // List of active particle ids
  int *celllist;                          // List of active cells
  int *interactlist;                      // ..
  int *neiblist;                          // List of neighbour ids
  FLOAT draux[ndimmax];                   // Aux. relative position vector var
  FLOAT drsqd;                            // Distance squared
  FLOAT hrangesqd;                        // Kernel extent
  FLOAT rp[ndimmax];                      // Local copy of particle position
  FLOAT *dr;                              // Array of relative position vectors
  FLOAT *drmag;                           // Array of neighbour distances
  FLOAT *invdrmag;                        // Array of 1/drmag between particles
  SphParticle *neibpart;                  // Local copy of neighbouring ptcls
  SphParticle parti;                      // Local copy of SPH particle
  SphParticle *data = sph->sphdata;       // Pointer to SPH particle data

  debug2("[GridSearch::UpdateAllSphProperties]");

  // Find list of all cells that contain active particles
  celllist = new int[Ncell];
  cactive = ComputeActiveCellList(celllist);
  Nneibmax = Nlistmax;


  // Set-up all OMP threads
  // ==========================================================================
#pragma omp parallel default(shared) private(activelist,c,cc,dr,draux,drmag,drsqd,hrangesqd,i,interactlist,invdrmag,j,jj,k,okflag,Nactive,neiblist,neibpart,Ninteract,Nneib,parti,rp)
  {
    activelist = new int[Noccupymax];
    neiblist = new int[Nneibmax];
    interactlist = new int[Nneibmax];
    dr = new FLOAT[Nneibmax*ndim];
    drmag = new FLOAT[Nneibmax];
    invdrmag = new FLOAT[Nneibmax];
    neibpart = new SphParticle[Nneibmax];

    // Loop over all active cells
    // ========================================================================
#pragma omp for schedule(dynamic)
    for (cc=0; cc<cactive; cc++) {
      c = celllist[cc];

      // Find list of active particles in current cell
      Nactive = ComputeActiveParticleList(c,activelist,sph);

      // Compute neighbour list for cell depending on physics options
      Nneib = ComputeNeighbourList(c,neiblist);
      for (j=0; j<Nneib; j++) neibpart[j] = data[neiblist[j]];

      // Loop over all active particles in the cell
      // ----------------------------------------------------------------------
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
        // --------------------------------------------------------------------
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
        // --------------------------------------------------------------------

        // Compute all gather neighbour contributions to hydro forces
        sph->ComputeSphDerivatives(i,Ninteract,interactlist,
				   drmag,invdrmag,dr,parti,neibpart);

	data[i] = parti;

      }
      // ----------------------------------------------------------------------

    }
    // ========================================================================

    // Free-up local memory for OpenMP thread
    delete[] neibpart;
    delete[] invdrmag;
    delete[] drmag;
    delete[] dr;
    delete[] interactlist;
    delete[] neiblist;
    delete[] activelist;
    
  }
  // ==========================================================================

  delete[] celllist;

  return;
}




// ============================================================================
// GridSearch::UpdateAllSphDudt
// Compute all local 'gather' properties of currently active particles, and 
// then compute each particle's contribution to its (active) neighbour 
// neighbour hydro forces.  Optimises the algorithm by using grid-cells to 
// construct local neighbour lists for all particles  inside the cell.
// ============================================================================
void GridSearch::UpdateAllSphDudt(Sph *sph)
{
  int c;                                  // Cell id
  int cactive;                            // No. of active cells
  int cc;                                 // Aux. cell counter
  int i;                                  // Particle id
  int j;                                  // Aux. particle counter
  int jj;                                 // Aux. particle counter
  int k;                                  // Dimension counter
  int okflag;                             // Flag if h-rho iteration is valid
  int Nactive;                            // No. of active particles in cell
  int Ngrav;                              // No. of direct sum gravity ptcls
  int Ninteract;                          // No. of near gather neighbours
  int Nneib;                              // No. of neighbours
  int Nneibmax;                           // Max. no. of neighbours
  int *activelist;                        // List of active particle ids
  int *celllist;                          // List of active cells
  int *interactlist;                      // ..
  int *neiblist;                          // List of neighbour ids
  FLOAT draux[ndimmax];                   // Aux. relative position vector var
  FLOAT drsqd;                            // Distance squared
  FLOAT hrangesqdi;                       // Kernel extent
  FLOAT hrangesqdj;                       // ..
  FLOAT rp[ndimmax];                      // Local copy of particle position
  FLOAT *dr;                              // Array of relative position vectors
  FLOAT *drmag;                           // Array of neighbour distances
  FLOAT *invdrmag;                        // Array of 1/drmag between particles
  SphParticle *neibpart;                  // Local copy of neighbouring ptcls
  SphParticle parti;                      // Local copy of SPH particle
  SphParticle *data = sph->sphdata;       // Pointer to SPH particle data

  debug2("[GridSearch::UpdateAllSphProperties]");

  // Find list of all cells that contain active particles
  celllist = new int[Ncell];
  cactive = ComputeActiveCellList(celllist);
  Nneibmax = Nlistmax;


  // Set-up all OMP threads
  // ==========================================================================
#pragma omp parallel default(shared) private(activelist,c,cc,dr,draux,drmag,drsqd,hrangesqdi,hrangesqdj,i,interactlist,invdrmag,j,jj,k,okflag,Nactive,neiblist,neibpart,Ninteract,Nneib,parti,rp)
  {
    activelist = new int[Noccupymax];
    neiblist = new int[Nneibmax];
    interactlist = new int[Nneibmax];
    dr = new FLOAT[Nneibmax*ndim];
    drmag = new FLOAT[Nneibmax];
    invdrmag = new FLOAT[Nneibmax];
    neibpart = new SphParticle[Nneibmax];

    // Loop over all active cells
    // ========================================================================
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
      // ----------------------------------------------------------------------
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
        // --------------------------------------------------------------------
        for (jj=0; jj<Nneib; jj++) {

          hrangesqdj = pow(sph->kernfac*sph->kernp->kernrange*neibpart[jj].h,2);
          for (k=0; k<ndim; k++) draux[k] = neibpart[jj].r[k] - rp[k];
          drsqd = DotProduct(draux,draux,ndim);

          // Compute list of particle-neighbour interactions and also
          // compute
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
        // --------------------------------------------------------------------

        // Compute all gather neighbour contributions to hydro forces
        sph->ComputeSphNeibDudt(i,Ninteract,interactlist,
				drmag,invdrmag,dr,parti,neibpart);

        // Add ..
#pragma omp atomic
        data[i].dudt += parti.dudt;
	
      }
      // ----------------------------------------------------------------------

      // Now add all active neighbour contributions to the main arrays
      for (jj=0; jj<Nneib; jj++) {
        if (neibpart[jj].active) {
          j = neiblist[jj];
#pragma omp atomic
          data[j].dudt += neibpart[jj].dudt;
        }
      }
      
    }
    // ========================================================================

    // Free-up local memory for OpenMP thread
    delete[] neibpart;
    delete[] invdrmag;
    delete[] drmag;
    delete[] dr;
    delete[] interactlist;
    delete[] neiblist;
    delete[] activelist;
    
  }
  // ==========================================================================

  delete[] celllist;

  return;
}



// ============================================================================
// GridSearch::UpdateAllSphGravityProperties
// Compute all local 'gather' properties of currently active particles, and 
// then compute each particle's contribution to its (active) neighbour 
// neighbour hydro forces.  Optimises the algorithm by using grid-cells to 
// construct local neighbour lists for all particles  inside the cell.
// ============================================================================
void GridSearch::UpdateAllSphGravityProperties(Sph *sph)
{
  int c;                                  // Cell id
  int cactive;                            // No. of active cells
  int cc;                                 // Aux. cell counter
  int i;                                  // Particle id
  int j;                                  // Aux. particle counter
  int jj;                                 // Aux. particle counter
  int k;                                  // Dimension counter
  int okflag;                             // Flag if h-rho iteration is valid
  int Nactive;                            // No. of active particles in cell
  int Ngrav;                              // No. of direct sum gravity ptcls
  int Ngather;                            // No. of near gather neighbours
  int Nneib;                              // No. of neighbours
  int Nneibmax;                           // Max. no. of neighbours
  int Nscatter;                           // No. of near scatter neighbours
  int *activelist;                        // List of active particle ids
  int *celllist;                          // List of active cells
  int *directgravlist;                    // Direct gravity
  int *neiblist;                          // List of neighbour ids
  int *gatherlist;                        // List of nearby neighbour ids
  int *scatterlist;                       // List of nearby neighbour ids
  FLOAT draux[ndimmax];                   // Aux. relative position vector var
  FLOAT drsqd;                            // Distance squared
  FLOAT hrangesqd;                        // Kernel extent
  FLOAT rp[ndimmax];                      // Local copy of particle position
  FLOAT *dr;                              // Array of relative position vectors
  FLOAT *drmag;                           // Array of neighbour distances
  FLOAT *invdrmag;                        // Array of 1/drmag between particles
  FLOAT *dr2;                             // Position vectors to scatter neibs
  FLOAT *drmag2;                          // Distances to scatter neibs
  FLOAT *invdrmag2;                       // 1/drmag to scatter neibs
  SphParticle *neibpart;                  // Local copy of neighbouring ptcls
  SphParticle *data = sph->sphdata;       // Pointer to SPH particle data

  debug2("[GridSearch::UpdateAllSphProperties]");

  // Find list of all cells that contain active particles
  celllist = new int[Ncell];
  cactive = ComputeActiveCellList(celllist);

  Nneibmax = Nlistmax;
  activelist = new int[Noccupymax];
  neiblist = new int[Nneibmax];
  gatherlist = new int[Nneibmax];
  scatterlist = new int[Nneibmax];
  dr = new FLOAT[Nneibmax*ndim];
  drmag = new FLOAT[Nneibmax];
  invdrmag = new FLOAT[Nneibmax];
  dr2 = new FLOAT[Nneibmax*ndim];
  drmag2 = new FLOAT[Nneibmax];
  invdrmag2 = new FLOAT[Nneibmax];
  neibpart = new SphParticle[Nneibmax];


  // Loop over all active cells
  // ==========================================================================
  for (cc=0; cc<cactive; cc++) {
    c = celllist[cc];

    // Find list of active particles in current cell
    Nactive = ComputeActiveParticleList(c,activelist,sph);

    // Compute neighbour list for cell depending on physics options
    Nneib = ComputeNeighbourList(c,neiblist);

    // Make local copies of all potential neighbours
    for (j=0; j<Nneib; j++) {
      neibpart[j] = data[neiblist[j]];
      neibpart[j].dudt = (FLOAT) 0.0;
      for (k=0; k<ndim; k++) neibpart[j].a[k] = (FLOAT) 0.0;
    }

    // Loop over all active particles in the cell
    // ------------------------------------------------------------------------
    for (j=0; j<Nactive; j++) {
      i = activelist[j];
      for (k=0; k<ndim; k++) rp[k] = data[i].r[k];

      // Set gather range as current h multiplied by some tolerance factor
      hrangesqd = pow(grid_h_tolerance*sph->kernp->kernrange*data[i].h,2);
      Ngather = 0;
      Nscatter = 0;

      // Compute distances and the inverse between the current particle
      // and all neighbours here, for both gather and inactive scatter neibs.
      // ----------------------------------------------------------------------
      for (jj=0; jj<Nneib; jj++) { 
	for (k=0; k<ndim; k++) draux[k] = neibpart[jj].r[k] - rp[k];
	drsqd = DotProduct(draux,draux,ndim);

	// Compute gather list for active particles
	if (drsqd <= hrangesqd) {
	  gatherlist[Ngather] = jj;
	  drmag[Ngather] = sqrt(drsqd);
	  invdrmag[Ngather] = (FLOAT) 1.0/(drmag[Ngather] + small_number);
	  for (k=0; k<ndim; k++)
	    dr[Ngather*ndim + k] = draux[k]*invdrmag[Ngather];
          Ngather++;
	}

	// Compute scatter list for inactive scatter neighbours
	if ((!neibpart[jj].active) && drsqd <= pow(sph->kernp->kernrange*
						   neibpart[jj].h,2)) {
	  scatterlist[Nscatter] = jj;
	  drmag2[Nscatter] = sqrt(drsqd);
	  invdrmag2[Nscatter] = (FLOAT) 1.0/(drmag2[Nscatter] + small_number);
	  for (k=0; k<ndim; k++)
	    dr2[Nscatter*ndim + k] = draux[k]*invdrmag2[Nscatter];
          Nscatter++;
	}
      }
      // ----------------------------------------------------------------------

      // Validate that gather neighbour list is correct
#if defined(VERIFY_ALL)
      if (neibcheck) CheckValidNeighbourList(sph,i,Nneib,neiblist,"gather");
#endif

      // Compute all gather neighbour contributions to hydro forces
      /*sph->ComputeGatherHydroForces(i,data[i],Nneib,Ngather,
				    gatherlist,neibpart,drmag,invdrmag,dr);

      // Compute all scatter neighbour contributions to hydro forces
      sph->ComputeScatterHydroForces(i,data[i],Nneib,Nscatter,
				     scatterlist,neibpart,
				     drmag2,invdrmag2,dr2);*/

    }
    // ------------------------------------------------------------------------

    // Now add all active neighbour contributions to the main arrays
    for (j=0; j<Nneib; j++) {
      if (data[neiblist[j]].active) {
	for (k=0; k<ndim; k++) data[neiblist[j]].a[k] += neibpart[j].a[k];
	data[neiblist[j]].dudt += neibpart[j].dudt;
      }
    }

  }
  // ==========================================================================

  // Free-up all memory
  delete[] neibpart;
  delete[] invdrmag2;
  delete[] drmag2;
  delete[] dr2;
  delete[] invdrmag;
  delete[] drmag;
  delete[] dr;
  delete[] scatterlist;
  delete[] gatherlist;
  delete[] neiblist;
  delete[] activelist;
  delete[] celllist;

  return;
}



// ============================================================================
// GridSearch::AllocateGridMemory
// Allocate memory for neighbour grid as requested.  If more memory is 
// required than currently allocated, grid is deallocated and reallocated here.
// ============================================================================
void GridSearch::AllocateGridMemory(int Npart)
{
  debug2("[GridSearch::AllocateGridMemory]");

  Ntot = Npart;

  if (Ntot > Ntotmax || Ncell > Ncellmax) {
    if (allocated_grid) DeallocateGridMemory();
    Ntotmax = 2*Ntot;
    Ncellmax = 2*Ncell;
    inext = new int[Ntotmax];
    grid = new struct GridCell[Ncellmax];
  }

  return;
}



// ============================================================================
// GridSearch::DeallocateGridMemory
// Deallocates all neighbour grid memory
// ============================================================================
void GridSearch::DeallocateGridMemory(void)
{
  debug2("[GridSearch::DeallocateGridMemory]");

  delete[] grid;
  delete[] inext;
  allocated_grid = false;

  return;
}



// ============================================================================
// GridSearch::CreateGrid
// Create a regular neighbour grid using all SPH particles contained within 
// the SPH object.  The grid spacing is equal to the maximum smoothing kernel 
// range of all particles multiplied by some arbitrary tolerance parameter 
// (grid_h_tolerance) to allow for some smoothing lengths increasing.
// ============================================================================
void GridSearch::CreateGrid(Sph *sph)
{
  int c;                                // Grid cell counter/id
  FLOAT h_max = 0.0;                    // Maximum smoothing length of ptcls

  debug2("[GridSearch::CreateGrid]");
  
  // Compute maximum smoothing length to determine optimum grid spacing
  for (int i=0; i<sph->Nsph; i++) h_max = max(h_max,sph->sphdata[i].h);
  dx_grid = grid_h_tolerance*sph->kernfac*sph->kernp->kernrange*h_max;

  // Compute bounding box of all particles
  sph->SphBoundingBox(rmax,rmin,sph->Ntot);

  // Calculate no. of grid cells in each dimension and in total
  Ncell = 1;
  for (int k=0; k<ndim; k++) {
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
  for (int i=0; i<Ntotmax; i++) inext[i] = -1;

  // Now attach all particles to grid cells
  // --------------------------------------------------------------------------
  for (int i=0; i<sph->Ntot; i++) {
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
  // --------------------------------------------------------------------------

  // Find maximum occupations of all cells
  Noccupymax = 0;
  for (c=0; c<Ncell; c++) Noccupymax = max(Noccupymax,grid[c].Nptcls);
  Nlistmax = Noccupymax*pow(3,ndim);

#if defined(VERIFY_ALL)
  ValidateGrid();
#endif

  return;
}



// ============================================================================
// GridSearch::ComputeParticleGridCell
// Compute and return the grid cell i.d. that contains the position 'rp'.
// ============================================================================
int GridSearch::ComputeParticleGridCell(FLOAT *rp)
{
  int igrid[ndimmax];                   // Grid cell coordinate

  for (int k=0; k<ndim; k++) {
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



// ============================================================================
// GridSearch::ComputeCellCoordinate
// Computes and returns the grid cell coordinate 'igrid' from the 
// grid cell i.d. 'c'
// ============================================================================
void GridSearch::ComputeCellCoordinate(int c, int igrid[ndimmax])
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
  return;
}



// ============================================================================
// GridSearch::ComputeActiveCellList
// Returns the number of cells containing active particles, 'Nactive', and 
// the i.d. list of cells contains active particles, 'celllist'
// ============================================================================
int GridSearch::ComputeActiveCellList(int *celllist)
{
  int Nactive = 0;                      // No. of cells containing active ptcls

  debug2("[GridSearch::ComputeActiveCellList]");

  for (int c=0; c<Ncell; c++)
    if (grid[c].Nactive > 0) celllist[Nactive++] = c;

  return Nactive;
}



// ============================================================================
// GridSearch::ComputeActiveParticleList
// Returns the number (Nactive) and list of ids (activelist) of all active 
// SPH particles in the given cell 'c'.
// ============================================================================
int GridSearch::ComputeActiveParticleList(int c, int *activelist, Sph *sph)
{
  int Nactive = 0;                      // No. of active particles in cell c
  int i = grid[c].ifirst;               // Particle id (set to first ptcl id)
  int ilast = grid[c].ilast;            // i.d. of last particle in cell c

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



// ============================================================================
// GridSearch::ComputeNeighbourList
// Computes and returns number of neighbour, 'Nneib', and the list 
// of neighbour ids, 'neiblist', for all particles inside cell 'c'.
// Includes all particles in the selected cell, plus all particles 
// contained in adjacent cells (including diagonal cells).
// ============================================================================
int GridSearch::ComputeNeighbourList(int c, int *neiblist)
{
  int i;                                // Particle id
  int ilast;                            // id of last particle in current cell
  int caux,cx,cy,cz;                    // Aux. cell counters and coordinates
  int igrid[ndimmax];                   // Grid cell coordinate
  int gridmin[ndimmax];                 // Minimum neighbour cell coordinate
  int gridmax[ndimmax];                 // Maximum neighbour cell coordinate
  int Nneib = 0;                        // No. of neighbours

  // Compute the location of the cell on the grid using the id
  ComputeCellCoordinate(c,igrid);

  // --------------------------------------------------------------------------
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
  // --------------------------------------------------------------------------
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
  // --------------------------------------------------------------------------
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
	  } while (i != 1);
	}
      }
    }
  }
  // --------------------------------------------------------------------------

  return Nneib;
}



#if defined(VERIFY_ALL)
// ============================================================================
// GridSearch::CheckValidNeighbourList
// Checks that the neighbour list generated by the grid is valid in that it 
// (i) does include all true neighbours, and 
// (ii) all true neigbours are only included once and once only.
// ============================================================================
void GridSearch::CheckValidNeighbourList(Sph *sph, int i, int Nneib, 
					 int *neiblist, string neibtype)
{
  int count;                            // Valid neighbour counter
  int j;                                // Neighbour particle counter
  int k;                                // Dimension counter
  int Ntrueneib = 0;                    // No. of 'true' neighbours
  int *trueneiblist;                    // List of true neighbour ids
  FLOAT drsqd;                          // Distance squared
  FLOAT dr[ndimmax];                    // Relative position vector

  // Allocate array to store local copy of potential neighbour ids
  trueneiblist = new int[sph->Ntot];

  // First, create list of 'true' neighbours by looping over all particles
  if (neibtype == "gather") {
    for (j=0; j<sph->Ntot; j++) {
      for (k=0; k<ndimmax; k++)
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



// ============================================================================
// GridSearch::ValidateGrid
// Validate that the grid structure is consistent with several sanity checks.
// ============================================================================
void GridSearch::ValidateGrid(void)
{
  int c;                                // Cell counter
  int i;                                // Particle id
  int *gridentry;                       // No. of time ptcl is present in grid

  debug2("[GridSearch::ValidateGrid]");

  gridentry = new int[Ntot];
  for (i=0; i<Ntot; i++) gridentry[i] = 0;

  // Loop over all grid cells and count how many times particles appear 
  // in linked lists
  for (c=0; c<Ncell; c++) {
    if (grid[c].Nptcls == 0) continue;
    i = grid[c].ifirst;
    do {
      gridentry[i]++;
      if (i == grid[c].ilast) break;
      i = inext[i];
    } while (1);
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
