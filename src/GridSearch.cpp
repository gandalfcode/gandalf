// ============================================================================
// GridSearch.cpp
// Contains functions for grid neighbour search routines.
// Creates a uniform grid from particle distribution where the spacing is 
// the size of the maximum kernel range (i.t. kernrange*h_max) over all ptcls.
// ============================================================================


#include <cstdlib>
#include <iostream>
#include <string>
#include <math.h>
#include "Precision.h"
#include "Dimensions.h"
#include "SphNeighbourSearch.h"
#include "Sph.h"
#include "Parameters.h"
#include "InlineFuncs.h"
#include "SphParticle.h"
#include "Debug.h"
using namespace std;


static FLOAT grid_h_tolerance = 1.2;


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
  int c;                                  // Cell id
  int cactive;                            // No. of active cells
  int cc;                                 // Aux. cell counter
  int i;                                  // Particle id
  int j;                                  // Aux. particle counter
  int jj;                                 // Aux. particle counter
  int k;                                  // Dimension counter
  int okflag;                             // Flag if h-rho iteration is valid
  int *activelist;                        // List of active particle ids
  int *celllist;                          // List of active cells
  int *neiblist;                          // List of neighbour ids
  int *nearlist;                          // List of nearby neighbour ids
  int Nactive;                            // No. of active particles in cell
  int Ngrav;                              // No. of direct sum gravity ptcls
  int Nnear;                              // No. of near neighbours
  int Nneib;                              // No. of neighbours
  int Nneibmax;                           // Max. no. of neighbours
  FLOAT draux[ndimmax];                   // Aux. relative position vector var
  FLOAT drsqd;                            // Distance squared
  FLOAT hrangesqd;                        // Kernel extent
  FLOAT rp[ndimmax];                      // Local copy of particle position
  FLOAT *dr;                              // Array of relative position vectors
  FLOAT *drmag;                           // Array of neighbour distances
  FLOAT *invdrmag;                        // Array of 1/drmag between particles
  SphParticle *neibpart;                  // Local copy of neighbouring ptcls
  SphParticle *data = sph->sphdata;       // Pointer to SPH particle data

  debug2("[GridSearch::UpdateAllSphProperties]");

  // Find list of all cells that contain active particles
  celllist = new int[Ncell];
  cactive = ComputeActiveCellList(celllist);

  Nneibmax = Nlistmax;
  neiblist = new int[Nneibmax];
  nearlist = new int[Nneibmax];
  dr = new FLOAT[Nneibmax*ndim];
  drmag = new FLOAT[Nneibmax];
  invdrmag = new FLOAT[Nneibmax];
  activelist = new int[Noccupymax];
  neibpart = new SphParticle[Nneibmax];


  // Loop over all active cells
  // ==========================================================================
  for (cc=0; cc<cactive; cc++) {
    c = celllist[cc];

    // Find list of active particles
    Nactive = ComputeActiveParticleList(c,activelist,sph);

    // Compute neighbour list for cell depending on physics options
    Nneib = ComputeNeighbourList(c,neiblist);

    // Make local copies of all potential neighbours
    for (j=0; j<Nneib; j++) {
      neibpart[j] = data[neiblist[j]];
      neibpart[j].dudt = 0.0;
      for (k=0; k<ndim; k++) neibpart[j].a[k] = 0.0;
    }

    // Loop over all active particles in the cell
    // ------------------------------------------------------------------------
    for (j=0; j<Nactive; j++) {
      i = activelist[j];
      for (k=0; k<ndim; k++) rp[k] = data[i].r[k];

      // Only compute quantities for definite candidiate neighbours
      hrangesqd = pow(grid_h_tolerance*sph->kernp->kernrange*data[i].h,2);

      // Compute distances and the reciprical between the current particle
      // and all neighbours here
      Nnear = 0;
      for (jj=0; jj<Nneib; jj++) { 
	for (k=0; k<ndim; k++) draux[k] = neibpart[jj].r[k] - rp[k];
	drsqd = DotProduct(draux,draux,ndim);
	if (drsqd <= hrangesqd) {
	  nearlist[Nnear] = jj;
	  drmag[Nnear] = sqrt(drsqd);
	  invdrmag[Nnear] = 1.0/(drmag[Nnear] + small_number);
	  for (k=0; k<ndim; k++) dr[Nnear*ndim + k] = draux[k]*invdrmag[Nnear];
          Nnear++;
	}
      }

      // Compute all SPH gather properties
      okflag = sph->ComputeH(i,data[i],Nneib,Nnear,
			     nearlist,neibpart,drmag,invdrmag,dr);

#if defined(VERIFY_ALL)
      if (neibcheck) CheckValidNeighbourList(sph,i,Nneib,neiblist,"gather");
#endif

      // Compute all current particle contributions to hydro forces
      sph->ComputeHydroForces(i,data[i],Nneib,Nnear,
			      nearlist,neibpart,drmag,invdrmag,dr);

    }
    // ------------------------------------------------------------------------

    // Now add all neighbour contributions to the main arrays
    for (j=0; j<Nneib; j++) {
      for (k=0; k<ndim; k++) data[neiblist[j]].a[k] += neibpart[j].a[k];
      data[neiblist[j]].dudt += neibpart[j].dudt;
    }

  }
  // ==========================================================================

  // Free-up all memory
  delete[] neibpart;
  delete[] activelist;
  delete[] invdrmag;
  delete[] drmag;
  delete[] dr;
  delete[] nearlist;
  delete[] neiblist;
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
  int Nactive;                            // No. of active particles in cell
  int Ngrav;                              // No. of direct sum gravity ptcls
  int Nneib;                              // No. of neighbours
  int Nneibmax;                           // Max. no. of neighbours
  int okflag;                             // Flag if h-rho iteration is valid
  int *activelist;                        // List of active particle ids
  int *celllist;                          // List of active cells
  int *neiblist;                          // List of neighbour ids
  int *gravlist;                          // Direct-sum gravity ids
  FLOAT draux[ndimmax];                   // Aux. relative position vector var
  FLOAT drsqd;                            // Distance squared
  FLOAT hrangesqd;                        // Kernel extent
  FLOAT rp[ndimmax];                      // Local copy of particle position
  FLOAT *dr;                              // Array of relative position vectors
  FLOAT *drmag;                           // Array of neighbour distances
  FLOAT *invdrmag;                        // Array of 1/drmag between particles
  SphParticle *neibpart;                  // Local copy of neighbouring ptcls
  SphParticle *data = sph->sphdata;       // Pointer to SPH particle data

  debug2("[GridSearch::UpdateAllSphProperties]");
  /*
  // Find list of all cells that contain active particles
  celllist = new int[Ncell];
  cactive = ComputeActiveCellList(celllist);

  Nneibmax = Nlistmax;
  neiblist = new int[Nneibmax];
  dr = new FLOAT[Nneibmax*ndim];
  drmag = new FLOAT[Nneibmax];
  invdrmag = new FLOAT[Nneibmax];
  activelist = new int[Noccupymax];
  neibpart = new SphParticle[Nneibmax];
  if (self_gravity == 1) gravlist = new int[sph->Nsph];


  // Loop over all active cells
  // ==========================================================================
  for (cc=0; cc<cactive; cc++) {
    c = celllist[cc];

    // Find list of active particles
    Nactive = ComputeActiveParticleList(c,activelist,sph);

    // Compute neighbour list for cell depending on physics options
    if (self_gravity == 1)
      Nneib = ComputeGravityNeighbourList(c,neiblist,gravlist);
    else
      Nneib = ComputeNeighbourList(c,neiblist);

    for (j=0; j<Nneib; j++) {
      neibpart[j] = data[neiblist[j]];
      neibpart[j].dudt = 0.0;
      for (k=0; k<ndim; k++) neibpart[j].a[k] = 0.0;
    }

    // Loop over all active particles in the cell
    // ------------------------------------------------------------------------
    for (j=0; j<Nactive; j++) {
      i = activelist[j];
      for (k=0; k<ndim; k++) rp[k] = data[i].r[k];

      // Only compute quantities for definite candidiate neighbours
      hrangesqd = pow(grid_h_tolerance*sph->kernp->kernrange*data[i].h,2);

      // Compute distances and the reciprical between the current particle
      // and all neighbours here
      for (jj=0; jj<Nneib; jj++) { 
	for (k=0; k<ndim; k++) draux[k] = neibpart[jj].r[k] - rp[k];
	drsqd = DotProduct(draux,draux,ndim);
	if (drsqd > hrangesqd) 
	  drmag[jj] = big_number;
	else {
	  drmag[jj] = sqrt(drsqd);
	  invdrmag[jj] = 1.0/(drmag[jj] + small_number);
	  for (k=0; k<ndim; k++) dr[jj*ndim + k] = draux[k]*invdrmag[jj];
	}
      }

      // Compute all SPH gather properties
      okflag = sph->ComputeH(i,data[i],Nneib,neibpart,drmag,invdrmag,dr);

#if defined(VERIFY_ALL)
      CheckValidNeighbourList(sph,i,Nneib,neiblist,"gather");
#endif

      // Compute all current particle contributions to hydro forces
      sph->ComputeHydroForces(i,data[i],Nneib,neibpart,drmag,invdrmag,dr);

      // Compute direct-sum gravity for all non-neighbours
      sph->ComputeDirectGravity(i,data[i],Ngrav,gravlist);

    }
    // ------------------------------------------------------------------------

    // Now add all neighbour contributions to the main arrays
    for (j=0; j<Nneib; j++) {
      for (k=0; k<ndim; k++) data[neiblist[j]].a[k] += neibpart[j].a[k];
      data[neiblist[j]].dudt += neibpart[j].dudt;
    }

  }   
  // ==========================================================================

  // Free-up all memory
  if (self_gravity == 1) delete[] gravlist;
  delete[] neibpart;
  delete[] activelist;
  delete[] invdrmag;
  delete[] drmag;
  delete[] dr;
  delete[] neiblist;
  delete[] celllist;


  // Allocate array to store local copy of potential neighbour ids
  // --------------------------------------------------------------------------
  if (sph->self_gravity == 1) {
    Nneib = sph->Ntot;

    // Compute SPH hydro forces for all particles
    for (i=0; i<sph->Nsph; i++)
      sph->ComputeGravForces(i,Nneib,sph->sphdata);
  }
  */
  return;
}



// ============================================================================
// GridSearch::AllocateGridMemory
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
// ============================================================================
void GridSearch::DeallocateGridMemory(void)
{
  debug2("[GridSearch::DeallocateGridMemory]");

  delete[] inext;
  delete[] grid;
  allocated_grid = false;

  return;
}



// ============================================================================
// GridSearch::CreateGrid
// ============================================================================
void GridSearch::CreateGrid(Sph *sph)
{
  int c;
  FLOAT h_max = 0.0;

  debug2("[GridSearch::CreateGrid]");
  
  // Compute maximum smoothing length to determine optimum grid spacing
  for (int i=0; i<sph->Nsph; i++) h_max = max(h_max,sph->sphdata[i].h);
  dx_grid = grid_h_tolerance*sph->kernp->kernrange*h_max;

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
  for (int c=0; c<Ncellmax; c++) {
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
    if (grid[c].Nptcls == 0) grid[c].ifirst = i;
    else inext[grid[c].ilast] = i;
    grid[c].ilast = i;
    grid[c].Nptcls++;
    if (sph->sphdata[i].active) grid[c].Nactive++;

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
  int igrid[ndimmax];

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
  int Nactive = 0;

  debug2("[GridSearch::ComputeActiveCellList]");

  for (int c=0; c<Ncell; c++)
    if (grid[c].Nactive > 0) celllist[Nactive++] = c;

  return Nactive;
}



// ============================================================================
// GridSearch::ComputeActiveParticleList
// ============================================================================
int GridSearch::ComputeActiveParticleList(int c, int *activelist, 
					  Sph *sph)
{
  int Nactive = 0;
  int i = grid[c].ifirst;
  int ilast = grid[c].ilast;

  // If there are no active particles in this cell, return without walking list
  if (grid[c].Nptcls == 0) return 0;

  // Else walk through linked list to obtain list and number of active ptcls.
  do {
    if (sph->sphdata[i].active) activelist[Nactive++] = i;
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
  int i;
  int ilast;
  int caux,cx,cy,cz;
  int igrid[ndimmax];
  int gridmin[ndimmax];
  int gridmax[ndimmax];
  int Nneib = 0;

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
    
    for (cz=gridmin[1]; cz<=gridmax[1]; cz++) {
      for (cy=gridmin[1]; cy<=gridmax[1]; cy++) {
	for (cx=gridmin[0]; cx<=gridmax[0]; cx++) {
	  caux = cx + cy*Ngrid[0] + cz*Ngrid[0]*Ngrid[1];
	  if (grid[caux].Nptcls == 0) continue;
	  i = grid[caux].ifirst;
	  ilast = grid[caux].ilast;
	  do {
	    neiblist[Nneib++] = i;
	    if (i != ilast) i = inext[i];
	  } while (i != ilast);
	}
      }
    }
  }
  // --------------------------------------------------------------------------

  return Nneib;
}



// ============================================================================
// GridSearch::CheckValidNeighbourList
// Checks that the neighbour list generated by the grid is valid in that it 
// (i) does include all true neighbours, and 
// (ii) all true neigbours are only included once and once only.
// ============================================================================
void GridSearch::CheckValidNeighbourList(Sph *sph, int i, int Nneib, 
					 int *neiblist, string neibtype)
{
  int count = 0;
  int j,k;
  int Ntrueneib = 0;
  int *trueneiblist;
  FLOAT drsqd;
  FLOAT dr[ndimmax];

  // Allocate array to store local copy of potential neighbour ids
  trueneiblist = new int[sph->Ntot];

  // First, create list of 'true' neighbours by looping over all particles
  if (neibtype == "gather") {
    for (j=0; j<sph->Ntot; j++) {
      for (k=0; k<ndimmax; k++)
	dr[k] = sph->sphdata[j].r[k] - sph->sphdata[i].r[k];
      drsqd = DotProduct(dr,dr,ndim);
      if (drsqd <= sph->kernp->kernrangesqd*sph->sphdata[i].h*sph->sphdata[i].h)
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
      exit(0);
    }
  }

  delete[] trueneiblist;

  return;
}



// ============================================================================
// GridSearch::ValidateGrid
// Validate that the grid structure is correct with several sanity checks.
// ============================================================================
void GridSearch::ValidateGrid(void)
{
  int c,i;
  int *gridentry;

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
      exit(0);
    }
  }

  // Check that linked list termination flag (i.e. inext = -1) is consistent
  // If not, quit with error message to screen
  for (c=0; c<Ncell; c++) {
    if (grid[c].Nptcls != 0 && inext[grid[c].ilast] != -1) {
      cout << "Error in linked list : " << c << "  " << grid[c].ilast << "  "
	   << inext[grid[c].ilast] << endl;
      exit(0);
    }
  }

  delete[] gridentry;

  return;
}
