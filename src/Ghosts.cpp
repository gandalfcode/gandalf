//=============================================================================
//  Ghosts.cpp
//  Contains all routines for searching for and creating ghost particles.
//  Also contains routine to correct particle positions/velocities to keep
//  them contained in simulation bounding box.
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
#include <math.h>
#include <map>
#include <string>
#include "Precision.h"
#include "Constants.h"
#include "Debug.h"
#include "Exception.h"
#include "Simulation.h"
#include "SphParticle.h"
#include "Sph.h"
#include "Ghosts.h"
using namespace std;



//=============================================================================
//  Ghosts::CheckBoundaries
/// Check all particles to see if any have crossed the simulation bounding
/// box.  If so, then move the particles to their new location on the other
/// side of the periodic box.
//=============================================================================
template <int ndim>
void PeriodicGhosts<ndim>::CheckBoundaries
(DomainBox<ndim> simbox,
 Sph<ndim> *sph)
{
  // Loop over all particles and check if any lie outside the periodic box.
  // If so, then re-position with periodic wrapping.
  //---------------------------------------------------------------------------
#pragma omp parallel for default(none) shared(simbox,sph)
  for (int i=0; i<sph->Nsph; i++) {
    SphParticle<ndim>& part = sph->GetParticleIPointer(i);


    if (part.r[0] < simbox.boxmin[0])
      if (simbox.x_boundary_lhs == periodicBoundary) {
        part.r[0] += simbox.boxsize[0];
        part.r0[0] += simbox.boxsize[0];
      }
    if (part.r[0] > simbox.boxmax[0])
      if (simbox.x_boundary_rhs == periodicBoundary) {
        part.r[0] -= simbox.boxsize[0];
        part.r0[0] -= simbox.boxsize[0];
      }

    if (ndim >= 2 && part.r[1] < simbox.boxmin[1])
      if (simbox.y_boundary_lhs == periodicBoundary) {
        part.r[1] += simbox.boxsize[1];
        part.r0[1] += simbox.boxsize[1];
      }
    if (ndim >= 2 && part.r[1] > simbox.boxmax[1])
      if (simbox.y_boundary_rhs == periodicBoundary) {
        part.r[1] -= simbox.boxsize[1];
        part.r0[1] -= simbox.boxsize[1];
      }

    if (ndim == 3 && part.r[2] < simbox.boxmin[2])
      if (simbox.z_boundary_lhs == periodicBoundary) {
        part.r[2] += simbox.boxsize[2];
        part.r0[2] += simbox.boxsize[2];
      }
    if (ndim == 3 && part.r[2] > simbox.boxmax[2])
      if (simbox.z_boundary_rhs == periodicBoundary) {
        part.r[2] -= simbox.boxsize[2];
        part.r0[2] -= simbox.boxsize[2];
      }

  }
  //---------------------------------------------------------------------------

  return;
}



//=============================================================================
//  Ghosts::SearchGhostParticles
/// Search domain to create any required ghost particles near any boundaries.
/// Currently only searches to create periodic or mirror ghost particles.
//=============================================================================
template <int ndim, template <int> class ParticleType>
void PeriodicGhostsSpecific<ndim, ParticleType >::SearchGhostParticles
(FLOAT tghost,                      ///< Ghost particle 'lifetime'
 DomainBox<ndim> simbox,            ///< Simulation box structure
 Sph<ndim> *sph)                    ///< Sph object pointer
{
  int i;                                                // Particle counter
  FLOAT kernrange = sph->kernp->kernrange*sph->kernfac; // Kernel extent

  ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sph->GetParticlesArray());

  // Set all relevant particle counters
  sph->Nghost    = 0;
  sph->NPeriodicGhost = 0;
  sph->Nghostmax = sph->Nsphmax - sph->Nsph;
  sph->Ntot      = sph->Nsph;


  // If all boundaries are open, immediately return to main loop
  if (simbox.x_boundary_lhs == openBoundary && simbox.x_boundary_rhs == openBoundary &&
      simbox.y_boundary_lhs == openBoundary && simbox.y_boundary_rhs == openBoundary &&
      simbox.z_boundary_lhs == openBoundary && simbox.z_boundary_rhs == openBoundary)
    return;

  debug2("[SphSimulation::SearchGhostParticles]");

  // Create ghost particles in x-dimension
  //---------------------------------------------------------------------------
  if ((simbox.x_boundary_lhs == openBoundary &&
       simbox.x_boundary_rhs == openBoundary) == 0) {
    for (i=0; i<sph->Ntot; i++) {

      SphParticle<ndim>& part = sph->GetParticleIPointer(i);

      if (part.r[0] + min(0.0,part.v[0]*tghost) <
          simbox.boxmin[0] + ghost_range*kernrange*part.h) {
        if (simbox.x_boundary_lhs == periodicBoundary)
          CreateGhostParticle(i,0,part.r[0] + simbox.boxsize[0],
                              part.v[0],sph,x_lhs_periodic,sphdata);
        if (simbox.x_boundary_lhs == mirrorBoundary)
          CreateGhostParticle(i,0,2.0*simbox.boxmin[0] - part.r[0],
                              -part.v[0],sph,x_lhs_mirror,sphdata);
      }
      if (part.r[0] + max(0.0,part.v[0]*tghost) >
          simbox.boxmax[0] - ghost_range*kernrange*part.h) {
        if (simbox.x_boundary_rhs == periodicBoundary)
          CreateGhostParticle(i,0,part.r[0] - simbox.boxsize[0],
                              part.v[0],sph,x_rhs_periodic,sphdata);
        if (simbox.x_boundary_rhs == mirrorBoundary)
          CreateGhostParticle(i,0,2.0*simbox.boxmax[0] - part.r[0],
                              -part.v[0],sph,x_rhs_mirror,sphdata);
      }
    }
    sph->Ntot = sph->Nsph + sph->Nghost;
  }


  // Create ghost particles in y-dimension
  //---------------------------------------------------------------------------
  if (ndim >= 2 && (simbox.y_boundary_lhs == openBoundary &&
		    simbox.y_boundary_rhs == openBoundary) == 0) {
    for (i=0; i<sph->Ntot; i++) {

      SphParticle<ndim>& part = sph->GetParticleIPointer(i);

      if (part.r[1] + min(0.0,part.v[1]*tghost) <
          simbox.boxmin[1] + ghost_range*kernrange*part.h) {
        if (simbox.y_boundary_lhs == periodicBoundary)
          CreateGhostParticle(i,1,part.r[1] + simbox.boxsize[1],
                              part.v[1],sph,y_lhs_periodic,sphdata);
	    if (simbox.y_boundary_lhs == mirrorBoundary)
          CreateGhostParticle(i,1,2.0*simbox.boxmin[1] - part.r[1],
                              -part.v[1],sph,y_lhs_mirror,sphdata);
      }
      if (part.r[1] + max(0.0,part.v[1]*tghost) >
          simbox.boxmax[1] - ghost_range*kernrange*part.h) {
        if (simbox.y_boundary_rhs == periodicBoundary)
          CreateGhostParticle(i,1,part.r[1] - simbox.boxsize[1],
                              part.v[1],sph,y_rhs_periodic,sphdata);
        if (simbox.y_boundary_rhs == mirrorBoundary)
          CreateGhostParticle(i,1,2.0*simbox.boxmax[1] - part.r[1],
                              -part.v[1],sph,y_rhs_mirror,sphdata);
      }
    }
    sph->Ntot = sph->Nsph + sph->Nghost;
  }


  // Create ghost particles in z-dimension
  //---------------------------------------------------------------------------
  if (ndim == 3 && (simbox.z_boundary_lhs == openBoundary &&
		    simbox.z_boundary_rhs == openBoundary) == 0) {
    for (i=0; i<sph->Ntot; i++) {

      SphParticle<ndim>& part = sph->GetParticleIPointer(i);

      if (part.r[2] + min(0.0,part.v[2]*tghost) <
          simbox.boxmin[2] + ghost_range*kernrange*part.h) {
        if (simbox.z_boundary_lhs == periodicBoundary)
          CreateGhostParticle(i,2,part.r[2] + simbox.boxsize[2],
                              part.v[2],sph,z_lhs_periodic,sphdata);
        if (simbox.z_boundary_lhs == mirrorBoundary)
          CreateGhostParticle(i,2,2.0*simbox.boxmin[2] - part.r[2],
                              -part.v[2],sph,z_lhs_mirror,sphdata);
      }
      if (part.r[2] + max(0.0,part.v[2]*tghost) >
          simbox.boxmax[2] - ghost_range*kernrange*part.h) {
        if (simbox.z_boundary_rhs == periodicBoundary)
          CreateGhostParticle(i,2,part.r[2] - simbox.boxsize[2],
                              part.v[2],sph,z_rhs_periodic,sphdata);
        if (simbox.z_boundary_rhs == mirrorBoundary)
          CreateGhostParticle(i,2,2.0*simbox.boxmax[2] - part.r[2],
                              -part.v[2],sph,z_rhs_mirror,sphdata);
      }
    }
    sph->Ntot = sph->Nsph + sph->Nghost;
  }

  // Quit here if we've run out of memory for ghosts
  if (sph->Ntot > sph->Nsphmax) {
    string message="Not enough memory for ghost particles";
    ExceptionHandler::getIstance().raise(message);
  }

  sph->NPeriodicGhost = sph->Nghost;

  return;
}



//=============================================================================
//  Ghosts::CreateGhostParticle
/// Create a new ghost particle from either
/// (i) a real SPH particle (i < Nsph), or
/// (ii) an existing ghost particle (i >= Nsph).
//=============================================================================
template <int ndim, template <int> class ParticleType>
void PeriodicGhostsSpecific<ndim, ParticleType >::CreateGhostParticle
(int i,                             ///< [in] i.d. of original particle
 int k,                             ///< [in] Boundary dimension for new ghost
 FLOAT rk,                          ///< [in] k-position of original particle
 FLOAT vk,                          ///< [in] k-velocity of original particle
 Sph<ndim> *sph,                    ///< [inout] SPH particle object pointer
 int ghosttype,                     ///< ..
 ParticleType<ndim>* sphdata)       ///< [in] Array with the SPH particles
{
  // Increase ghost counter and check there's enough space in memory
  if (sph->Nghost > sph->Nghostmax) {
    string message= "Not enough memory for new ghost";
    ExceptionHandler::getIstance().raise(message);
  }

  int id_new_ghost = sph->Nsph + sph->Nghost;

  // If there's enough memory, create ghost particle in arrays
  sphdata[id_new_ghost] = sphdata[i];
  sphdata[id_new_ghost].r[k] = rk;
  sphdata[id_new_ghost].v[k] = vk;
  sphdata[id_new_ghost].active = false;
  sphdata[id_new_ghost].itype = ghosttype;

  // Record id of original particle for later copying
  //if (i >= sph->Nsph)
  //  sph->sphdata[sph->Nsph + sph->Nghost].iorig = sphdata[i].iorig;
  //else
  //  sph->sphdata[sph->Nsph + sph->Nghost].iorig = i;
  sphdata[id_new_ghost].iorig = i;

  sph->Nghost = sph->Nghost + 1;

  return;
}



//=============================================================================
//  Ghosts::CopySphDataToGhosts
/// Copy any newly calculated data from original SPH particles to ghosts.
//=============================================================================
template <int ndim, template <int> class ParticleType>
void PeriodicGhostsSpecific<ndim, ParticleType >::CopySphDataToGhosts
(DomainBox<ndim> simbox,
 Sph<ndim> *sph)
{
  int i;                            // Particle id
  int iorig;                        // Original (real) particle id
  int itype;                        // Ghost particle type
  int j;                            // Ghost particle counter
  ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sph->GetParticlesArray());

  debug2("[SphSimulation::CopySphDataToGhosts]");


  //---------------------------------------------------------------------------
//#pragma omp parallel for default(none) private(i,iorig,itype,j) shared(simbox,sph,sphdata)
  for (j=0; j<sph->NPeriodicGhost; j++) {
    i = sph->Nsph + j;
    iorig = sphdata[i].iorig;
    itype = sphdata[i].itype;

    sphdata[i] = sphdata[iorig];
    sphdata[i].iorig = iorig;
    sphdata[i].itype = itype;
    sphdata[i].active = false;

    // Modify ghost position based on ghost type
    if (itype == x_lhs_periodic)
      sphdata[i].r[0] += simbox.boxsize[0];
    else if (itype == x_rhs_periodic)
      sphdata[i].r[0] -= simbox.boxsize[0];
    else if (itype == y_lhs_periodic && ndim > 1)
      sphdata[i].r[1] += simbox.boxsize[1];
    else if (itype == y_rhs_periodic && ndim > 1)
      sphdata[i].r[1] -= simbox.boxsize[1];
    else if (itype == z_lhs_periodic && ndim == 3)
      sphdata[i].r[2] += simbox.boxsize[2];
    else if (itype == z_rhs_periodic && ndim == 3)
      sphdata[i].r[2] -= simbox.boxsize[2];

  }
  //---------------------------------------------------------------------------

  return;
}



//=============================================================================
//  NullGhosts::CheckBoundaries
/// Empty function when no ghost particles are required.
//=============================================================================
template <int ndim>
void NullGhosts<ndim>::CheckBoundaries(DomainBox<ndim> simbox, Sph<ndim> *sph)
{
  return;
}



//=============================================================================
//  NullGhosts::SearchGhostParticles
/// Empty function when no ghost particles are required.
//=============================================================================
template <int ndim>
void NullGhosts<ndim>::SearchGhostParticles
(FLOAT tghost,                      ///< Ghost particle 'lifetime'
 DomainBox<ndim> simbox,            ///< Simulation box structure
 Sph<ndim> *sph)                    ///< Sph object pointer
{

  // Set all relevant particle counters
  sph->Nghost         = 0;
  sph->NPeriodicGhost = 0;
  sph->Nghostmax      = sph->Nsphmax - sph->Nsph;
  sph->Ntot           = sph->Nsph;

  return;
}



//=============================================================================
//  NullGhosts::CopySphDataToGhosts
/// Empty function when no ghost particles are required.
//=============================================================================
template <int ndim>
void NullGhosts<ndim>::CopySphDataToGhosts(DomainBox<ndim> simbox, Sph<ndim> *sph) {
  return;
}



#if defined MPI_PARALLEL
//=============================================================================
//  MpiGhosts::CheckBoundaries
/// ..
//=============================================================================
template <int ndim>
void MPIGhosts<ndim>::CheckBoundaries(DomainBox<ndim> simbox, Sph<ndim> *sph)
{
  return;
}



//=============================================================================
//  MpiGhosts::SearchGhostParticles
/// Handle control to MpiControl to compute particles to send to other nodes
/// and receive from them, then copy received ghost particles inside the main
/// arrays.
//=============================================================================
template <int ndim>
void MPIGhosts<ndim>::SearchGhostParticles
(FLOAT tghost,                      ///< Ghost particle 'lifetime'
 DomainBox<ndim> simbox,            ///< Simulation box structure
 Sph<ndim> *sph)                    ///< Sph object pointer
{
  int i;
  int j;
  SphParticle<ndim>* ghost_array;
  int Nmpighosts = mpicontrol->SendReceiveGhosts(&ghost_array, sph);

  if (sph->Ntot + Nmpighosts > sph->Nsphmax) {
    cout << "Error: not enough memory for MPI ghosts!!! " << Nmpighosts
         << " " << sph->Ntot << " " << sph->Nsphmax<<endl;
    ExceptionHandler::getIstance().raise("");
  }

  SphParticle<ndim>* main_array = sph->sphdata;
  int start_index = sph->Nsph + sph->NPeriodicGhost;

  for (j=0; j<Nmpighosts; j++) {
    i = start_index + j;
    main_array[i] =  ghost_array[j];
    main_array[i].active = false;
  }

  sph->Nghost += Nmpighosts;
  sph->Ntot += Nmpighosts;

  if (sph->Nghost > sph->Nghostmax || sph->Ntot > sph->Nsphmax) {
	cout << "Error: not enough memory for MPI ghosts!!! " << Nmpighosts
             << " " << sph->Ntot << " " << sph->Nsphmax<<endl;
	ExceptionHandler::getIstance().raise("");
  }

}



//=============================================================================
//  MpiGhosts::CopySphDataToGhosts
/// ..
//=============================================================================
template <int ndim>
void MPIGhosts<ndim>::CopySphDataToGhosts
(DomainBox<ndim> simbox, Sph<ndim> *sph)
{
  SphParticle<ndim>* ghost_array;
  SphParticle<ndim>* main_array = sph->sphdata;
  int Nmpighosts = mpicontrol->UpdateGhostParticles(&ghost_array);
  int start_index = sph->Nsph + sph->NPeriodicGhost;

  for (int j=0; j<Nmpighosts; j++) {
    int i = start_index + j;
    main_array[i] =  ghost_array[j];
    main_array[i].active = false;
  }

}
#endif



// Create template class instances of the main SphSimulation object for
// each dimension used (1, 2 and 3)
template class NullGhosts<1>;
template class NullGhosts<2>;
template class NullGhosts<3>;
template class PeriodicGhostsSpecific<1, GradhSphParticle>;
template class PeriodicGhostsSpecific<2, GradhSphParticle>;
template class PeriodicGhostsSpecific<3, GradhSphParticle>;
template class PeriodicGhostsSpecific<1, GodunovSphParticle>;
template class PeriodicGhostsSpecific<2, GodunovSphParticle>;
template class PeriodicGhostsSpecific<3, GodunovSphParticle>;
template class PeriodicGhostsSpecific<1, SM2012SphParticle>;
template class PeriodicGhostsSpecific<2, SM2012SphParticle>;
template class PeriodicGhostsSpecific<3, SM2012SphParticle>;
#ifdef MPI_PARALLEL
template class MPIGhosts<1>;
template class MPIGhosts<2>;
template class MPIGhosts<3>;
#endif
