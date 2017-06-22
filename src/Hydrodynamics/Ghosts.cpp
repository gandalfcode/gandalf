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
  ParticleType<ndim>* sphdata = hydro->template GetParticleArray<ParticleType>();

  debug2("[PeriodicGhosts::CopyHydroDataToGhosts]");


  //-----------------------------------------------------------------------------------------------
//#pragma omp parallel for default(none) private(i,iorig,itype,j) shared(simbox,sph,sphdata)
  for (j=0; j<hydro->NPeriodicGhost; j++) {
    i = hydro->Nhydro + j;
    iorig = sphdata[i].iorig;
    itype = sphdata[i].flags.get();
    assert(itype != none) ;
    assert(sphdata[i].ptype == sphdata[iorig].ptype);

    sphdata[i] = sphdata[iorig];
    sphdata[i].iorig = iorig;
    sphdata[i].flags = type_flag(itype);
    sphdata[i].flags.unset(active);

    // Modify ghost position based on ghost type
    // Ghosts of ghosts refer only to their previous ghosts not the base cell, so
    // only update one direction.
    if (ndim > 2) {
      if (itype & z_periodic_lhs) {
        sphdata[i].r[2] += simbox.size[2];
      }
      else if (itype & z_periodic_rhs) {
    	sphdata[i].r[2] -= simbox.size[2];
      }
      else if (itype & z_mirror_lhs) {
        reflect(sphdata[i], 2, simbox.min[2]);
      }
      else if (itype & z_mirror_rhs) {
        reflect(sphdata[i], 2, simbox.max[2]);
      }

    }
    if (ndim > 1) {
      if (itype & y_periodic_lhs) {
    	sphdata[i].r[1] += simbox.size[1];
      }
      else if (itype & y_periodic_rhs) {
    	sphdata[i].r[1] -= simbox.size[1];
      }
      else if (itype & y_mirror_lhs) {
    	reflect(sphdata[i], 1, simbox.min[1]);
      }
      else if (itype & y_mirror_rhs) {
        reflect(sphdata[i], 1, simbox.max[1]);
      }
    }

    if (itype & x_periodic_lhs) {
      sphdata[i].r[0] += simbox.size[0];
    }
    else if (itype & x_periodic_rhs) {
      sphdata[i].r[0] -= simbox.size[0];
    }
    else if (itype & x_mirror_lhs) {
      reflect(sphdata[i], 0, simbox.min[0]);
    }
    else if (itype & x_mirror_rhs) {
      reflect(sphdata[i], 0, simbox.max[0]);
    }
  }
  //-----------------------------------------------------------------------------------------------

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
	  hydro->AllocateMemory(hydro->Ntot + Nmpighosts);
  }

  ParticleType<ndim>* main_array = hydro->template GetParticleArray<ParticleType>();
  int start_index = hydro->Nhydro + hydro->NPeriodicGhost;

  for (j=0; j<Nmpighosts; j++) {
    i = start_index + j;
    main_array[i] = ghost_array[j];
    main_array[i].flags.unset(active);
  }

  hydro->Nmpighost = Nmpighosts;
  hydro->Nghost += Nmpighosts;
  hydro->Ntot += Nmpighosts;


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
  debug2("[MpiGhostsSpecific::CopyHydroDataToGhosts]");


  ParticleType<ndim>* ghost_array;
  ParticleType<ndim>* main_array = hydro->template GetParticleArray<ParticleType>();
  int Nmpighosts = mpicontrol->UpdateGhostParticles(main_array,&ghost_array);
  int start_index = hydro->Nhydro + hydro->NPeriodicGhost;

  assert(hydro->Nmpighost == Nmpighosts);

  for (int j=0; j<Nmpighosts; j++) {
    int i = start_index + j;

    int jghost = main_array[i].iorig ;

    assert(main_array[i].ptype == ghost_array[jghost].ptype);
    assert(main_array[i].iorig == ghost_array[jghost].iorig);

    main_array[i] = ghost_array[jghost];
    main_array[i].flags.unset(active);
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
template class PeriodicGhostsSpecific<1, MeshlessFVParticle>;
template class PeriodicGhostsSpecific<2, MeshlessFVParticle>;
template class PeriodicGhostsSpecific<3, MeshlessFVParticle>;


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
template class MpiGhostsSpecific<1, MeshlessFVParticle>;
template class MpiGhostsSpecific<2, MeshlessFVParticle>;
template class MpiGhostsSpecific<3, MeshlessFVParticle>;
#endif
