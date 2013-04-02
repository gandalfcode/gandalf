// ============================================================================
// GhostParticles.cpp
// ============================================================================


#include <cstdlib>
#include <math.h>
#include <map>
#include <string>
#include "Precision.h"
#include "Constants.h"
#include "Debug.h"
#include "Exception.h"
#include "Sph.h"
#include "SphSimulation.h"
#include "SphParticle.h"
using namespace std;


static const FLOAT ghost_range = 1.3;


// ============================================================================
// 
// ============================================================================
template <int ndim>
void SphSimulation<ndim>::CheckBoundaries(void)
{
  int i;
  SphParticle<ndim> *part;

  // x-boundary conditions
  for (i=0; i<sph->Nsph; i++) {
    part = &sph->sphdata[i];

    if (part->r[0] < simbox.boxmin[0])
      if (simbox.x_boundary_lhs == "periodic") part->r[0] += simbox.boxsize[0];
    if (part->r[0] > simbox.boxmax[0])
      if (simbox.x_boundary_rhs == "periodic") part->r[0] -= simbox.boxsize[0];


#if NDIM==2 || NDIM==3
    if (ndim >= 2 && part->r[1] < simbox.boxmin[1])
      if (simbox.y_boundary_lhs == "periodic") part->r[1] += simbox.boxsize[1];
    if (ndim >= 2 && part->r[1] > simbox.boxmax[1])
      if (simbox.y_boundary_rhs == "periodic") part->r[1] -= simbox.boxsize[1];
#endif

#if NDIM==3
    if (ndim == 3 && part->r[2] < simbox.boxmin[2])
      if (simbox.z_boundary_lhs == "periodic") part->r[2] += simbox.boxsize[2];
    if (ndim == 3 && part->r[2] > simbox.boxmax[2])
      if (simbox.z_boundary_rhs == "periodic") part->r[2] -= simbox.boxsize[2];
#endif

  }

  return;
}



// ============================================================================
// SphSimulation::SearchGhostParticles
// ============================================================================
template <int ndim>
void SphSimulation<ndim>::SearchGhostParticles(void)
{
  int i;
  int k;
  SphParticle<ndim>* sphdata = sph->sphdata;
  FLOAT kernrange = sph->kernp->kernrange*sph->kernfac;

  // Set all relevant particle counters
  sph->Nghost    = 0;
  sph->Nghostmax = sph->Nsphmax - sph->Nsph;
  sph->Ntot      = sph->Nsph;

  // If all boundaries are open, immediately return to main loop
  if (simbox.x_boundary_lhs == "open" && simbox.x_boundary_rhs == "open" &&
      simbox.y_boundary_lhs == "open" && simbox.y_boundary_rhs == "open" &&
      simbox.z_boundary_lhs == "open" && simbox.z_boundary_rhs == "open")
    return;

  debug2("[SphSimulation::SearchGhostParticles]");

  // Create ghost particles in x-dimension
  // --------------------------------------------------------------------------
  if ((simbox.x_boundary_lhs == "open" && 
       simbox.x_boundary_rhs == "open") == 0) {
    for (i=0; i<sph->Ntot; i++) {
      if (sphdata[i].r[0] < simbox.boxmin[0] + 
	  ghost_range*kernrange*sphdata[i].h) {
	if (simbox.x_boundary_lhs == "periodic")
	  CreateGhostParticle(i,0,sphdata[i].r[0] + simbox.boxsize[0],
			      sphdata[i].v[0],
			      sphdata[i].r[0] - simbox.boxmin[0]);
	if (simbox.x_boundary_lhs == "mirror")
	  CreateGhostParticle(i,0,2.0*simbox.boxmin[0] - 
			      sphdata[i].r[0],-sphdata[i].v[0],
			      sphdata[i].r[0] - simbox.boxmin[0]);
      }
      if (sphdata[i].r[0] > simbox.boxmax[0] - 
	  ghost_range*kernrange*sphdata[i].h) {
	if (simbox.x_boundary_rhs == "periodic")
	  CreateGhostParticle(i,0,sphdata[i].r[0] - simbox.boxsize[0],
			      sphdata[i].v[0],
			      simbox.boxmax[0] - sphdata[i].r[0]);
	if (simbox.x_boundary_rhs == "mirror")
	  CreateGhostParticle(i,0,2.0*simbox.boxmax[0] - 
			      sphdata[i].r[0],-sphdata[i].v[0],
			      simbox.boxmax[0] - sphdata[i].r[0]);
      }
    }
    sph->Ntot = sph->Nsph + sph->Nghost;
  }


  // Create ghost particles in y-dimension
  // --------------------------------------------------------------------------
#if NDIM==2 || NDIM==3
  if (ndim >= 2 && (simbox.y_boundary_lhs == "open" && 
		    simbox.y_boundary_rhs == "open") == 0) {
    for (i=0; i<sph->Ntot; i++) {
      if (sphdata[i].r[1] < simbox.boxmin[1] + 
	  ghost_range*kernrange*sphdata[i].h) {
	if (simbox.y_boundary_lhs == "periodic")
	  CreateGhostParticle(i,1,sphdata[i].r[1] + simbox.boxsize[1],
			      sphdata[i].v[1],
			      sphdata[i].r[1] - simbox.boxmin[1]);
	if (simbox.y_boundary_lhs == "mirror")
	  CreateGhostParticle(i,1,2.0*simbox.boxmin[1] - 
			      sphdata[i].r[1],-sphdata[i].v[1],
			      sphdata[i].r[1] - simbox.boxmin[1]);
      }
      if (sphdata[i].r[1] > simbox.boxmax[1] - 
	  ghost_range*kernrange*sphdata[i].h) {
	if (simbox.y_boundary_rhs == "periodic")
	  CreateGhostParticle(i,1,sphdata[i].r[1] - simbox.boxsize[1],
			      sphdata[i].v[1],
			      simbox.boxmax[1] - sphdata[i].r[1]);
	if (simbox.y_boundary_rhs == "mirror")
	  CreateGhostParticle(i,1,2.0*simbox.boxmax[1] - 
			      sphdata[i].r[1],-sphdata[i].v[1],
			      simbox.boxmax[1] - sphdata[i].r[1]);
      }
    }
    sph->Ntot = sph->Nsph + sph->Nghost;
  }
#endif


 // Create ghost particles in z-dimension
  // --------------------------------------------------------------------------
#if NDIM==3
  if (ndim == 3 && (simbox.z_boundary_lhs == "open" && 
		    simbox.z_boundary_rhs == "open") == 0) {
    for (i=0; i<sph->Ntot; i++) {
      if (sphdata[i].r[2] < simbox.boxmin[2] + 
	  ghost_range*kernrange*sphdata[i].h) {
	if (simbox.z_boundary_lhs == "periodic")
	  CreateGhostParticle(i,2,sphdata[i].r[2] + simbox.boxsize[2],
			      sphdata[i].v[2],
			      sphdata[i].r[2] - simbox.boxmin[2]);
	if (simbox.z_boundary_lhs == "mirror")
	  CreateGhostParticle(i,2,2.0*simbox.boxmin[2] - 
			      sphdata[i].r[2],-sphdata[i].v[2],
			      sphdata[i].r[2] - simbox.boxmin[2]);
      }
      if (sphdata[i].r[2] > simbox.boxmax[2] - 
	  ghost_range*kernrange*sphdata[i].h) {
	if (simbox.z_boundary_rhs == "periodic")
	  CreateGhostParticle(i,2,sphdata[i].r[2] - simbox.boxsize[2],
			      sphdata[i].v[2],
			      simbox.boxmax[2] - sphdata[i].r[2]);
	if (simbox.z_boundary_rhs == "mirror")
	  CreateGhostParticle(i,2,2.0*simbox.boxmax[2] - 
			      sphdata[i].r[2],-sphdata[i].v[2],
			      simbox.boxmax[2] - sphdata[i].r[2]);
      }
    }
    sph->Ntot = sph->Nsph + sph->Nghost;
  }
#endif

  // Quit here if we've run out of memory for ghosts
  if (sph->Ntot > sph->Nsphmax) {
    string message="Not enough memory for ghost particles";
    ExceptionHandler::getIstance().raise(message);
  }

  return;
}



// ============================================================================
// SphSimulation::CreateGhostParticle
// ============================================================================
template <int ndim>
void SphSimulation<ndim>::CreateGhostParticle(int i, int k,
					FLOAT rk, FLOAT vk, FLOAT bdist)
{
  // Increase ghost counter and check there's enough space in memory
  if (sph->Nghost > sph->Nghostmax) {
    string message= "Not enough memory for new ghost";
    ExceptionHandler::getIstance().raise(message);
  }

  // If there's enough memory, create ghost particle in arrays
  sph->sphdata[sph->Nsph + sph->Nghost] = sph->sphdata[i];
  sph->sphdata[sph->Nsph + sph->Nghost].r[k] = rk;
  sph->sphdata[sph->Nsph + sph->Nghost].v[k] = vk;
  sph->sphdata[sph->Nsph + sph->Nghost].active = false;


  // Record id of original particle for later copying
  if (i >= sph->Nsph)
    sph->sphdata[sph->Nsph + sph->Nghost].iorig = sph->sphdata[i].iorig;
  else
    sph->sphdata[sph->Nsph + sph->Nghost].iorig = i;

  sph->Nghost = sph->Nghost + 1;

  return;
}



// ============================================================================
// SphSimulation::CopySphDataToGhosts
// ============================================================================
template <int ndim>
void SphSimulation<ndim>::CopySphDataToGhosts(void)
{
  int i;
  int iorig;
  int j;
  int k;
  FLOAT rp[ndim];
  FLOAT vp[ndim];

  // --------------------------------------------------------------------------
#pragma omp parallel for default(shared) private(i,iorig,k,rp,vp)
  for (j=0; j<sph->Nghost; j++) {
    i = sph->Nsph + j;
    iorig = sph->sphdata[i].iorig;

    for (k=0; k<ndim; k++) rp[k] = sph->sphdata[i].r[k];
    for (k=0; k<ndim; k++) vp[k] = sph->sphdata[i].v[k];
    
    sph->sphdata[i] = sph->sphdata[iorig];
    sph->sphdata[i].iorig = iorig;
    sph->sphdata[i].active = false;
    for (k=0; k<ndim; k++) sph->sphdata[i].r[k] = rp[k];
    for (k=0; k<ndim; k++) sph->sphdata[i].v[k] = vp[k];
    
  }
  // --------------------------------------------------------------------------

  return;
}



// ============================================================================
// SphSimulation::CopyAccelerationsFromGhosts
// ============================================================================
template <int ndim>
void SphSimulation<ndim>::CopyAccelerationFromGhosts(void)
{
  int i;
  int iorig;
  int j;
  int k;

  // --------------------------------------------------------------------------
#pragma omp parallel for default(shared) private(i,iorig,k)
  for (j=0; j<sph->Nghost; j++) {
    i = sph->Nsph + j;
    iorig = sph->sphdata[i].iorig;
    
    // Only look at active ghosts
    if (!sph->sphdata[iorig].active) continue;
    
    for (k=0; k<ndim; k++) {
#pragma omp atomic
      sph->sphdata[iorig].a[k] += sph->sphdata[i].a[k];
    }
#pragma omp atomic
    sph->sphdata[iorig].dudt += sph->sphdata[i].dudt;
#pragma omp atomic
    sph->sphdata[iorig].div_v += sph->sphdata[i].div_v;
    
  }

  return;
}
