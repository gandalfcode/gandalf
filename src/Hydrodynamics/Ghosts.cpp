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
#include "Particle.h"
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
  for (int i=0; i<sph->Nhydro; i++) {
    SphParticle<ndim>& part = sph->GetSphParticlePointer(i);


    if (part.r[0] < simbox.boxmin[0])
      if (simbox.boundary_lhs[0] == periodicBoundary) {
        part.r[0] += simbox.boxsize[0];
        part.r0[0] += simbox.boxsize[0];
      }
    if (part.r[0] > simbox.boxmax[0])
      if (simbox.boundary_rhs[0] == periodicBoundary) {
        part.r[0] -= simbox.boxsize[0];
        part.r0[0] -= simbox.boxsize[0];
      }

    if (ndim >= 2 && part.r[1] < simbox.boxmin[1])
      if (simbox.boundary_lhs[1] == periodicBoundary) {
        part.r[1] += simbox.boxsize[1];
        part.r0[1] += simbox.boxsize[1];
      }
    if (ndim >= 2 && part.r[1] > simbox.boxmax[1])
      if (simbox.boundary_rhs[1] == periodicBoundary) {
        part.r[1] -= simbox.boxsize[1];
        part.r0[1] -= simbox.boxsize[1];
      }

    if (ndim == 3 && part.r[2] < simbox.boxmin[2])
      if (simbox.boundary_lhs[2] == periodicBoundary) {
        part.r[2] += simbox.boxsize[2];
        part.r0[2] += simbox.boxsize[2];
      }
    if (ndim == 3 && part.r[2] > simbox.boxmax[2])
      if (simbox.boundary_rhs[2] == periodicBoundary) {
        part.r[2] -= simbox.boxsize[2];
        part.r0[2] -= simbox.boxsize[2];
      }

  }
  //---------------------------------------------------------------------------

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
  ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sph->GetSphParticleArray());

  debug2("[SphSimulation::CopySphDataToGhosts]");


  //---------------------------------------------------------------------------
//#pragma omp parallel for default(none) private(i,iorig,itype,j) shared(simbox,sph,sphdata)
  for (j=0; j<sph->NPeriodicGhost; j++) {
    i = sph->Nhydro + j;
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
template <int ndim, template <int> class ParticleType>
void MPIGhostsSpecific<ndim, ParticleType>::SearchGhostParticles
(FLOAT tghost,                      ///< Ghost particle 'lifetime'
 DomainBox<ndim> simbox,            ///< Simulation box structure
 Sph<ndim> *sph)                    ///< Sph object pointer
{
  int i;
  int j;
  ParticleType<ndim>* ghost_array;
  int Nmpighosts = mpicontrol->SendReceiveGhosts(tghost,sph,&ghost_array);

  if (sph->Ntot + Nmpighosts > sph->Nhydromax) {
    cout << "Error: not enough memory for MPI ghosts!!! " << Nmpighosts
         << " " << sph->Ntot << " " << sph->Nhydromax<<endl;
    ExceptionHandler::getIstance().raise("");
  }

  ParticleType<ndim>* main_array = static_cast<ParticleType<ndim>* > (sph->GetSphParticleArray() );
  int start_index = sph->Nhydro + sph->NPeriodicGhost;

  for (j=0; j<Nmpighosts; j++) {
    i = start_index + j;
    main_array[i] = ghost_array[j];
    main_array[i].active = false;
  }

  sph->Nmpighost = Nmpighosts;
  sph->Nghost += Nmpighosts;
  sph->Ntot += Nmpighosts;

  if (sph->Nghost > sph->Nghostmax || sph->Ntot > sph->Nhydromax) {
	cout << "Error: not enough memory for MPI ghosts!!! " << Nmpighosts
             << " " << sph->Ntot << " " << sph->Nhydromax<<endl;
	ExceptionHandler::getIstance().raise("");
  }

}



//=============================================================================
//  MpiGhosts::CopySphDataToGhosts
/// ..
//=============================================================================
template <int ndim, template<int> class ParticleType >
void MPIGhostsSpecific<ndim, ParticleType>::CopySphDataToGhosts
 (DomainBox<ndim> simbox,              ///< ..
  Sph<ndim> *sph)                      ///< ..
{
  ParticleType<ndim>* ghost_array;
  ParticleType<ndim>* main_array = static_cast<ParticleType<ndim>* > (sph->GetSphParticleArray() );
  int Nmpighosts = mpicontrol->UpdateGhostParticles(&ghost_array);
  int start_index = sph->Nhydro + sph->NPeriodicGhost;

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
//template class PeriodicGhostsSpecific<1, MeshlessFVParticle>;
//template class PeriodicGhostsSpecific<2, MeshlessFVParticle>;
//template class PeriodicGhostsSpecific<3, MeshlessFVParticle>;

#ifdef MPI_PARALLEL
template class MPIGhosts<1>;
template class MPIGhosts<2>;
template class MPIGhosts<3>;
template class MPIGhostsSpecific<1, GradhSphParticle>;
template class MPIGhostsSpecific<2, GradhSphParticle>;
template class MPIGhostsSpecific<3, GradhSphParticle>;
template class MPIGhostsSpecific<1, GodunovSphParticle>;
template class MPIGhostsSpecific<2, GodunovSphParticle>;
template class MPIGhostsSpecific<3, GodunovSphParticle>;
template class MPIGhostsSpecific<1, SM2012SphParticle>;
template class MPIGhostsSpecific<2, SM2012SphParticle>;
template class MPIGhostsSpecific<3, SM2012SphParticle>;
#endif
