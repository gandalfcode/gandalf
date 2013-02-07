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
#include "Dimensions.h"
#include "SphNeighbourSearch.h"
#include "Sph.h"
#include "Parameters.h"
#include "InlineFuncs.h"
#include "SphParticle.h"
#include "Debug.h"
using namespace std;



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
// ============================================================================
void GridSearch::UpdateAllSphProperties(Sph *sph, Parameters &simparams)
{
  int c;
  int i;
  int j;
  int okflag;
  int *celllist;
  int *neiblist;
  int *activelist;
  int Nactive;
  int Nneib;
  int Nneibmax;
  int cactive;

  debug2("[GridSearch::UpdateAllSphProperties]");

  // Find list of all cells that contain active particles
  celllist = new int[Ncell];
  cactive = ComputeActiveCellList(celllist);

  Nneibmax = Nlistmax;
  neiblist = new int[Nneibmax];
  activelist = new int[Noccupymax];

  // Loop over all active cells
  // --------------------------------------------------------------------------
  for (int cc=0; cc<cactive; cc++) {
    c = celllist[cc];

    // Find list of active particles
    Nactive = ComputeActiveParticleList(c,activelist,sph);

    // Compute neighbour list for cell
    Nneib = ComputeNeighbourList(c,neiblist);

    // Creates neighbour list local copy
    SphParticle * neiblistpart = new SphParticle[Nneib];

    // Copies particle from the main array to the new array
    SphParticle * data = sph->sphdata;
    for (j=0; j<Nneib; j++) {
      neiblistpart[j] = data[neiblist[j]];
    }

    // Loop over all active particles in the cell
    for (j=0; j<Nactive; j++) {
      i = activelist[j];
      okflag = sph->ComputeH(i,Nneib,neiblistpart,simparams);
      //CheckValidNeighbourList(sph,i,Nneib,neiblist,"gather");
    }

    delete[] neiblistpart;

    // Compute all other SPH properties
    for (j=0; j<Nactive; j++) {
      i = activelist[j];
      sph->ComputeSphProperties(i,Nneib,neiblist,simparams);
    }

  }   
  // --------------------------------------------------------------------------

  delete[] activelist;
  delete[] neiblist;
  delete[] celllist;

  return;
}



// ============================================================================
// GridSearch::UpdateAllSphForces
// ============================================================================
void GridSearch::UpdateAllSphForces(Sph *sph, Parameters &params)
{
  int c;
  int i;
  int j;
  int okflag;
  int *celllist;
  int *neiblist;
  int *activelist;
  int Nactive;
  int Nneib;
  int Nneibmax;
  int cactive;

  debug2("[GridSearch::UpdateAllSphForces]");

  // Find list of all cells that contain active particles
  celllist = new int[Ncell];
  cactive = ComputeActiveCellList(celllist);

  Nneibmax = Nlistmax;
  neiblist = new int[Nneibmax];
  activelist = new int[Noccupymax];

  // Loop over all active cells
  // --------------------------------------------------------------------------
  if (params.intparams["hydro_forces"] == 1) {
    for (int cc=0; cc<cactive; cc++) {
      c = celllist[cc];
      
      // Find list of active particles
      Nactive = ComputeActiveParticleList(c,activelist,sph);
      
      // Compute neighbour list for cell
      Nneib = ComputeNeighbourList(c,neiblist);


      // Loop over all active particles in the cell
      for (j=0; j<Nactive; j++) {
	i = activelist[j];
	//CheckValidNeighbourList(sph,i,Nneib,neiblist,"all");
	sph->ComputeHydroForces(i,Nneib,neiblist,params);
      }
      
    }   
  }
  // --------------------------------------------------------------------------

  delete[] activelist;
  delete[] neiblist;
  delete[] celllist;


  // Allocate array to store local copy of potential neighbour ids
  // --------------------------------------------------------------------------
  if (params.intparams["self_gravity"] == 1) {
    Nneib = sph->Ntot;
    neiblist = new int[sph->Ntot];
    for (int i=0; i<sph->Ntot; i++) neiblist[i] = i;
    
    // Compute SPH hydro forces for all particles
    for (int i=0; i<sph->Nsph; i++)
      sph->ComputeGravForces(i,Nneib,neiblist);

    // Free up memory from local array
    delete[] neiblist;
  }

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
  float h_max = 0.0;
  static float grid_h_tolerance = 1.6;

  debug2("[GridSearch::CreateGrid]");
  
  // Compute maximum smoothing length to determine optimum grid spacing
  for (int i=0; i<sph->Nsph; i++) h_max = max(h_max,sph->sphdata[i].h);
  dx_grid = grid_h_tolerance*sph->kern->kernrange*h_max;

  //cout << "Grid size : " << dx_grid << "   " << h_max << endl;

  // Compute bounding box of all particles
  sph->SphBoundingBox(rmax,rmin,sph->Ntot);

  // Calculate no. of grid cells in each dimension and in total
  Ncell = 1;
  for (int k=0; k<ndim; k++) {
    Ngrid[k] = (int)((rmax[k] - rmin[k])/dx_grid) + 1;
    Ncell = Ngrid[k]*Ncell;
  }

  //cout << "Ngrid : " << Ncell << endl;

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

  //cout << "Noccupymax : " << Noccupymax << "   " << Nlistmax << endl;
  //ValidateGrid();

  return;
}



// ============================================================================
// GridSearch::ComputeParticleGridCell
// ============================================================================
int GridSearch::ComputeParticleGridCell(float *rp)
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
// GridSearch::ComputeParticleGridCell
// ============================================================================
void GridSearch::ComputeCellCoordinate(int c, int igrid[ndimmax])
{
  if (ndim == 1) 
    igrid[0] = c;
  else if (ndim == 2) {
    igrid[1] = c/Ngrid[0];
    igrid[0] = c%Ngrid[0];
    if (c != igrid[0] + igrid[1]*Ngrid[0]) {
      cout << "Problem with ComputeCellCoordinate : " << c << "   " << 
	igrid[0] + igrid[1]*Ngrid[0] << endl; exit(0);
    }
  }
  else if (ndim == 3) {
    igrid[2] = c/Ngrid[0]/Ngrid[1];
    igrid[1] = (c/Ngrid[0])%Ngrid[1];
    igrid[0] = c%Ngrid[0];
  }
  //if (c != igrid[0]) {
  //  cout << "Problem with ComputeCellCoordinate : " << c << "   " << 
  //    igrid[0] << endl;
  //  exit(0);
  // }
  return;
}



// ============================================================================
// GridSearch::ComputeActiveCellList
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

  if (grid[c].Nptcls == 0) return 0;

  do {
    if (sph->sphdata[i].active) activelist[Nactive++] = i;
    if (i == ilast) break;
    i = inext[i];
  } while (i != -1);

  return Nactive;
}



// ============================================================================
// GridSearch::ComputeNeighbourList
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
// GridSearch::UpdateAllSphProperties
// ============================================================================
void GridSearch::CheckValidNeighbourList(Sph *sph, int i, int Nneib, 
					 int *neiblist, string neibtype)
{
  int count = 0;
  int validflag = true;
  int j,k;
  int Ntrueneib = 0;
  int *trueneiblist;
  float drsqd;
  float dr[ndimmax];

  // Allocate array to store local copy of potential neighbour ids
  trueneiblist = new int[sph->Ntot];

  // First, create 'true' neighbour list using brute-force
  if (neibtype == "gather") {
    for (j=0; j<sph->Ntot; j++) {
      for (k=0; k<ndimmax; k++)
	dr[k] = sph->sphdata[j].r[k] - sph->sphdata[i].r[k];
      drsqd = DotProduct(dr,dr,ndim);
      if (drsqd <= sph->kern->kernrangesqd*sph->sphdata[i].h*sph->sphdata[i].h)
	trueneiblist[Ntrueneib++] = j;
    }
  }
  else if (neibtype == "all") {
    for (j=0; j<sph->Ntot; j++) {
      for (k=0; k<ndimmax; k++)
	dr[k] = sph->sphdata[j].r[k] - sph->sphdata[i].r[k];
      drsqd = DotProduct(dr,dr,ndim);
      if (drsqd <= sph->kern->kernrangesqd*sph->sphdata[i].h*sph->sphdata[i].h 
	  && drsqd <= sph->kern->kernrangesqd*
	  sph->sphdata[j].h*sph->sphdata[j].h)
	trueneiblist[Ntrueneib++] = j;
    }
  }

  // Now check provided neighbour list for validation
  for (j=0; j<Ntrueneib; j++) {
    count = 0;
    validflag = false;
    for (k=0; k<Nneib; k++)
      if (neiblist[k] == trueneiblist[j]) count++;
    if (count != 1) {
      cout << "Problem with neighbour lists : " << i << "  " << j << "   "
	   << count << "   "
	   << sph->sphdata[i].r[0] << "   " << sph->sphdata[i].h << endl;
      cout << "Nneib : " << Nneib << "    Ntrueneib : " << Ntrueneib << endl;
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
// ============================================================================
void GridSearch::ValidateGrid(void)
{
  int c,i;
  int *gridentry;

  debug2("[GridSearch::ValidateGrid]");

  gridentry = new int[Ntot];
  for (i=0; i<Ntot; i++) gridentry[i] = 0;

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
