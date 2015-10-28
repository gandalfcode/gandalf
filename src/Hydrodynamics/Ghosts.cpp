//=================================================================================================
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
//=================================================================================================


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


// Declare ghost_range constant here (prevents warnings with some compilers)
template <int ndim>
const FLOAT Ghosts<ndim>::ghost_range = 1.6;


//=================================================================================================
//  Ghosts::CheckBoundaries
/// Check all particles to see if any have crossed the simulation bounding box.
/// If so, then move the particles to their new location on the other side of the periodic box.
//================================================================================================
template <int ndim>
void PeriodicGhosts<ndim>::CheckBoundaries
 (DomainBox<ndim> simbox,
  Hydrodynamics<ndim> *hydro)
{
  debug2("[PeriodicGhosts::CheckBoundaries]");

  // Loop over all particles and check if any lie outside the periodic box.
  // If so, then re-position with periodic wrapping.
  //===============================================================================================
#pragma omp parallel for default(none) shared(simbox,hydro)
  for (int i=0; i<hydro->Nhydro; i++) {
    Particle<ndim>& part = hydro->GetParticlePointer(i);

    // --------------------------------------------------------------------------------------------
    for (int k=0; k<ndim; k++) {

      // Check if particle has crossed LHS boundary
      //-------------------------------------------------------------------------------------------
      if (part.r[k] < simbox.boxmin[k]) {

        // Check if periodic boundary
        if (simbox.boundary_lhs[k] == periodicBoundary) {
          part.r[k]  += simbox.boxsize[k];
          part.r0[k] += simbox.boxsize[k];
        }

        // Check if wall or mirror boundary
        if (simbox.boundary_lhs[k] == mirrorBoundary || simbox.boundary_lhs[k] == wallBoundary) {
          part.r[k]  = (FLOAT) 2.0*simbox.boxmin[k] - part.r[k];
          part.r0[k] = (FLOAT) 2.0*simbox.boxmin[k] - part.r0[k];
          part.v[k]  = -part.v[k];
          part.v0[k] = -part.v0[k];
          part.a[k]  = -part.a[k];
          part.a0[k] = -part.a0[k];
        }

      }

      // Check if particle has crossed RHS boundary
      //-------------------------------------------------------------------------------------------
      if (part.r[k] > simbox.boxmax[k]) {

        // Check if periodic boundary
        if (simbox.boundary_rhs[k] == periodicBoundary) {
          part.r[k]  -= simbox.boxsize[k];
          part.r0[k] -= simbox.boxsize[k];
        }

        // Check if wall or mirror boundary
        if (simbox.boundary_rhs[k] == mirrorBoundary || simbox.boundary_rhs[k] == wallBoundary) {
          part.r[k]  = (FLOAT) 2.0*simbox.boxmax[k] - part.r[k];
          part.r0[k] = (FLOAT) 2.0*simbox.boxmax[k] - part.r0[k];
          part.v[k]  = -part.v[k];
          part.v0[k] = -part.v0[k];
          part.a[k]  = -part.a[k];
          part.a0[k] = -part.a0[k];
        }

      }


    }
    //---------------------------------------------------------------------------------------------

  }
  //===============================================================================================

  return;
}



//=================================================================================================
//  Ghosts::CopyHydroDataToGhosts
/// Copy any newly calculated data from original SPH particles to ghosts.
//=================================================================================================
template <int ndim, template <int> class ParticleType>
void PeriodicGhostsSpecific<ndim, ParticleType >::CopyHydroDataToGhosts
(DomainBox<ndim> simbox,
 Hydrodynamics<ndim> *hydro)
{
  int i;                            // Particle id
  int iorig;                        // Original (real) particle id
  int itype;                        // Ghost particle type
  int j;                            // Ghost particle counter
  ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (hydro->GetParticleArray());

  debug2("[SphSimulation::CopyHydroDataToGhosts]");


  //-----------------------------------------------------------------------------------------------
//#pragma omp parallel for default(none) private(i,iorig,itype,j) shared(simbox,sph,sphdata)
  for (j=0; j<hydro->NPeriodicGhost; j++) {
    i = hydro->Nhydro + j;
    iorig = sphdata[i].iorig;
    itype = sphdata[i].itype;

    sphdata[i] = sphdata[iorig];
    sphdata[i].iorig = iorig;
    sphdata[i].itype = itype;
    sphdata[i].active = false;

    // Modify ghost position based on ghost type
    if (itype == x_lhs_periodic) {
      sphdata[i].r[0] += simbox.boxsize[0];
    }
    else if (itype == x_rhs_periodic) {
      sphdata[i].r[0] -= simbox.boxsize[0];
    }
    else if (itype == x_lhs_mirror) {
      sphdata[i].r[0] = 2.0*simbox.boxmin[0] - sphdata[i].r[0];
    }
    else if (itype == x_rhs_mirror) {
      sphdata[i].r[0] = 2.0*simbox.boxmax[0] - sphdata[i].r[0];
    }
    else if (ndim > 1 && itype == y_lhs_periodic) {
      sphdata[i].r[1] += simbox.boxsize[1];
    }
    else if (ndim > 1 && itype == y_rhs_periodic) {
      sphdata[i].r[1] -= simbox.boxsize[1];
    }
    else if (ndim > 1 && itype == y_lhs_mirror) {
      sphdata[i].r[1] = 2.0*simbox.boxmin[1] - sphdata[i].r[1];
    }
    else if (ndim > 1 && itype == y_rhs_mirror) {
      sphdata[i].r[1] = 2.0*simbox.boxmax[1] - sphdata[i].r[1];
    }
    else if (ndim == 3 && itype == z_lhs_periodic) {
      sphdata[i].r[2] += simbox.boxsize[2];
    }
    else if (ndim == 3 && itype == z_rhs_periodic) {
      sphdata[i].r[2] -= simbox.boxsize[2];
    }
    else if (ndim == 3 && itype == z_lhs_mirror) {
      sphdata[i].r[2] = 2.0*simbox.boxmin[2] - sphdata[i].r[2];
    }
    else if (ndim == 3 && itype == z_rhs_mirror) {
      sphdata[i].r[2] = 2.0*simbox.boxmax[2] - sphdata[i].r[2];
    }

  }
  //-----------------------------------------------------------------------------------------------

  return;
}



//=================================================================================================
//  NullGhosts::CheckBoundaries
/// Empty function when no ghost particles are required.
//=================================================================================================
template <int ndim>
void NullGhosts<ndim>::CheckBoundaries(DomainBox<ndim> simbox, Hydrodynamics<ndim> *hydro)
{
  return;
}



//=================================================================================================
//  NullGhosts::CopyHydroDataToGhosts
/// Empty function when no ghost particles are required.
//=================================================================================================
template <int ndim>
void NullGhosts<ndim>::CopyHydroDataToGhosts(DomainBox<ndim> simbox, Hydrodynamics<ndim> *hydro)
{
  return;
}



#if defined MPI_PARALLEL
//=================================================================================================
//  MpiGhosts::CheckBoundaries
/// ..
//=================================================================================================
template <int ndim>
void MpiGhosts<ndim>::CheckBoundaries(DomainBox<ndim> simbox, Hydrodynamics<ndim> *hydro)
{
  return;
}



//=================================================================================================
//  MpiGhosts::SearchGhostParticles
/// Handle control to MpiControl to compute particles to send to other nodes and receive from them,
/// then copy received ghost particles inside the main arrays.
//=================================================================================================
template <int ndim, template <int> class ParticleType>
void MpiGhostsSpecific<ndim, ParticleType>::SearchGhostParticles
 (FLOAT tghost,                        ///< [in] Ghost particle 'lifetime'
  DomainBox<ndim> simbox,              ///< [in] Simulation box structure
  Hydrodynamics<ndim> *hydro)          ///< [inout] Sph object pointer
{
  int i;
  int j;
  ParticleType<ndim>* ghost_array;
  int Nmpighosts = mpicontrol->SendReceiveGhosts(tghost, hydro, &ghost_array);

  if (hydro->Ntot + Nmpighosts > hydro->Nhydromax) {
    cout << "Error: not enough memory for MPI ghosts!!! " << Nmpighosts
         << " " << hydro->Ntot << " " << hydro->Nhydromax<<endl;
    ExceptionHandler::getIstance().raise("");
  }

  ParticleType<ndim>* main_array = static_cast<ParticleType<ndim>* > (hydro->GetParticleArray() );
  int start_index = hydro->Nhydro + hydro->NPeriodicGhost;

  for (j=0; j<Nmpighosts; j++) {
    i = start_index + j;
    main_array[i] = ghost_array[j];
    main_array[i].active = false;
  }

  hydro->Nmpighost = Nmpighosts;
  hydro->Nghost += Nmpighosts;
  hydro->Ntot += Nmpighosts;

  if (hydro->Nghost > hydro->Nghostmax || hydro->Ntot > hydro->Nhydromax) {
    cout << "Error: not enough memory for MPI ghosts!!! " << Nmpighosts
             << " " << hydro->Ntot << " " << hydro->Nhydromax<<endl;
    ExceptionHandler::getIstance().raise("");
  }

}



//=================================================================================================
//  MpiGhosts::CopyHydroDataToGhosts
/// Copy all hydro data from the real particles to the MPI ghost particles
//=================================================================================================
template <int ndim, template<int> class ParticleType >
void MpiGhostsSpecific<ndim, ParticleType>::CopyHydroDataToGhosts
 (DomainBox<ndim> simbox,              ///< [in] Simulation box
  Hydrodynamics<ndim> *hydro)          ///< [inout] Pointer to hydrodynamics object
{
  ParticleType<ndim>* ghost_array;
  ParticleType<ndim>* main_array = static_cast<ParticleType<ndim>* > (hydro->GetParticleArray() );
  int Nmpighosts = mpicontrol->UpdateGhostParticles(&ghost_array);
  int start_index = hydro->Nhydro + hydro->NPeriodicGhost;

  for (int j=0; j<Nmpighosts; j++) {
    int i = start_index + j;
    main_array[i] = ghost_array[j];
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
template class PeriodicGhostsSpecific<1, SM2012SphParticle>;
template class PeriodicGhostsSpecific<2, SM2012SphParticle>;
template class PeriodicGhostsSpecific<3, SM2012SphParticle>;

#ifdef MPI_PARALLEL
template class MpiGhosts<1>;
template class MpiGhosts<2>;
template class MpiGhosts<3>;
template class MpiGhostsSpecific<1, GradhSphParticle>;
template class MpiGhostsSpecific<2, GradhSphParticle>;
template class MpiGhostsSpecific<3, GradhSphParticle>;
template class MpiGhostsSpecific<1, SM2012SphParticle>;
template class MpiGhostsSpecific<2, SM2012SphParticle>;
template class MpiGhostsSpecific<3, SM2012SphParticle>;
#endif
