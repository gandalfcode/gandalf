//=================================================================================================
//  Hydrodynamics.cpp
//  Contains important default routines for Hydrodynamics class.
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


#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <math.h>
#include "Precision.h"
#include "Hydrodynamics.h"
#include "SmoothingKernel.h"
#include "Particle.h"
#include "Parameters.h"
#include "EOS.h"
#include "Debug.h"
#include "InlineFuncs.h"
using namespace std;



// Declare invndim constant here (prevents warnings with some compilers)
template <int ndim>
const FLOAT Hydrodynamics<ndim>::invndim = 1.0/ndim;



//=================================================================================================
//  Hydrodynamics::Hydrodynamics
/// Constructor for parent Hydrodynamics class.  Initialises important variables and
/// sets important parameters using initialialisation lists.
//=================================================================================================
template <int ndim>
Hydrodynamics<ndim>::Hydrodynamics(int hydro_forces_aux, int self_gravity_aux, FLOAT h_fac_aux,
                                   string gas_eos_aux, string KernelName, int size_hydro):
  size_hydro_part(size_hydro),
  hydro_forces(hydro_forces_aux),
  self_gravity(self_gravity_aux),
  gas_eos(gas_eos_aux),
  h_fac(h_fac_aux),
  kerntab(TabulatedKernel<ndim>(KernelName))
{
  allocated          = false;
  Nhydro             = 0;
  Nhydromax          = 0;
  NImportedParticles = 0;
  Nmpighost          = 0;
  NPeriodicGhost     = 0;

  for (int itype = 0; itype < Ntypes; ++itype)
	  types[itype] = ParticleType() ;

  // Set flags for gas particle type
  //-----------------------------------------------------------------------------------------------
  types[gas_type].hydro_forces   = true;
  types[gas_type].self_gravity   = true;
  //for (int k=0; k<Ntypes; k++) types[gas].hmask[k] = true;
  //for (int k=0; k<Ntypes; k++) types[gas].hydromask[k] = true;
  types[gas_type].hmask[gas_type]     = true;
  types[gas_type].hmask[cdm_type]     = true;
  types[gas_type].hydromask[gas_type] = true;
  types[gas_type].gravmask[gas_type]  = true;
  types[gas_type].gravmask[cdm_type]  = true;
  types[gas_type].gravmask[dust_type] = true;


  // Set flags for cdm particle type
  //-----------------------------------------------------------------------------------------------
  types[cdm_type].self_gravity        = true;
  types[cdm_type].hmask[gas_type]     = true;
  types[cdm_type].hmask[cdm_type]     = true;
  types[cdm_type].gravmask[gas_type]  = true;
  types[cdm_type].gravmask[cdm_type]  = true;
  types[cdm_type].gravmask[dust_type] = true;


  // Set flags for dust particle type
  //-----------------------------------------------------------------------------------------------
  types[dust_type].self_gravity        = true;
  types[dust_type].hmask[dust_type]    = true;
  types[dust_type].gravmask[gas_type]  = true;
  types[dust_type].gravmask[cdm_type]  = true;
  types[dust_type].gravmask[dust_type] = true;

}



//=================================================================================================
//  Hydrodynamics::ComputeBoundingBox
/// Calculate the bounding box containing all hydro particles.
//=================================================================================================
template <int ndim>
void Hydrodynamics<ndim>::ComputeBoundingBox
 (FLOAT rmax[ndim],                    ///< [out] Maximum extent of bounding box
  FLOAT rmin[ndim],                    ///< [out] Minimum extent of bounding box
  const int Nmax)                      ///< [in] Maximum particle i.d. in loop
{
  int i;                               // Particle counter
  int k;                               // Dimension counter

  debug2("[Hydrodynamics::ComputeBoundingBox]");

  for (k=0; k<ndim; k++) rmin[k] = big_number;
  for (k=0; k<ndim; k++) rmax[k] = -big_number;

  for (i=0; i<Nmax; i++) {
    Particle<ndim>& part = GetParticlePointer(i);
    for (k=0; k<ndim; k++) rmin[k] = min(rmin[k],part.r[k]);
    for (k=0; k<ndim; k++) rmax[k] = max(rmax[k],part.r[k]);
  }

  return;
}



//=================================================================================================
//  Hydrodynamics::CheckXBoundaryGhostParticle
/// Check if we must create a ghost replica of particle i in the x-direction.  Checks how deep
/// the ghost region is likely to be based on the particle's x-velocity.
//=================================================================================================
template <int ndim>
void Hydrodynamics<ndim>::CheckXBoundaryGhostParticle
 (const int i,                         ///< [in] i.d. of particles to check
  const FLOAT tghost,                  ///< [in] Expected lifetime of ghost
  const DomainBox<ndim> &simbox)       ///< [in] Simulation domain box
{
  Particle<ndim>& part = GetParticlePointer(i);

  if (part.r[0] + min((FLOAT) 0.0, part.v[0]*tghost) <
      simbox.boxmin[0] + ghost_range*kernrange*part.h) {
    if (simbox.boundary_lhs[0] == periodicBoundary) {
      CreateBoundaryGhostParticle(i, 0, x_lhs_periodic, part.r[0] + simbox.boxsize[0], part.v[0]);
    }
    if (simbox.boundary_lhs[0] == mirrorBoundary) {
      CreateBoundaryGhostParticle(i, 0, x_lhs_mirror,
                                  (FLOAT) 2.0*simbox.boxmin[0] - part.r[0], -part.v[0]);
    }
  }
  if (part.r[0] + max((FLOAT) 0.0, part.v[0]*tghost) >
      simbox.boxmax[0] - ghost_range*kernrange*part.h) {
    if (simbox.boundary_rhs[0] == periodicBoundary) {
      CreateBoundaryGhostParticle(i, 0, x_rhs_periodic, part.r[0] - simbox.boxsize[0], part.v[0]);
    }
    if (simbox.boundary_rhs[0] == mirrorBoundary) {
      CreateBoundaryGhostParticle(i, 0, x_rhs_mirror,
                                  (FLOAT) 2.0*simbox.boxmax[0] - part.r[0], -part.v[0]);
    }
  }

  return;
}



//=================================================================================================
//  Hydrodynamics::CheckYBoundaryGhostParticle
/// Check if we must create a ghost replica of particle i in the y-direction.  Checks how deep
/// the ghost region is likely to be based on the particle's y-velocity.
//=================================================================================================
template <int ndim>
void Hydrodynamics<ndim>::CheckYBoundaryGhostParticle
 (const int i,                         ///< [in] i.d. of particles to check
  const FLOAT tghost,                  ///< [in] Expected lifetime of ghost
  const DomainBox<ndim> &simbox)       ///< [in] Simulation domain box
{
  if (ndim > 1) {
    Particle<ndim>& part = GetParticlePointer(i);

    if (part.r[1] + min((FLOAT) 0.0, part.v[1]*tghost) <
        simbox.boxmin[1] + ghost_range*kernrange*part.h) {
      if (simbox.boundary_lhs[1] == periodicBoundary) {
        CreateBoundaryGhostParticle(i, 1, y_lhs_periodic, part.r[1] + simbox.boxsize[1], part.v[1]);
      }
      if (simbox.boundary_lhs[1] == mirrorBoundary) {
        CreateBoundaryGhostParticle(i, 1, y_lhs_mirror,
                                    (FLOAT) 2.0*simbox.boxmin[1] - part.r[1], -part.v[1]);
      }
    }
    if (part.r[1] + max((FLOAT) 0.0, part.v[1]*tghost) >
        simbox.boxmax[1] - ghost_range*kernrange*part.h) {
      if (simbox.boundary_rhs[1] == periodicBoundary) {
        CreateBoundaryGhostParticle(i, 1, y_rhs_periodic, part.r[1] - simbox.boxsize[1], part.v[1]);
      }
      if (simbox.boundary_rhs[1] == mirrorBoundary) {
        CreateBoundaryGhostParticle(i, 1, y_rhs_mirror,
                                    (FLOAT) 2.0*simbox.boxmax[1] - part.r[1], -part.v[1]);
      }
    }
  }

  return;
}



//=================================================================================================
//  Hydrodynamics::CheckZBoundaryGhostParticle
/// Check if we must create a ghost replica of particle i in the z-direction.  Checks how deep
/// the ghost region is likely to be based on the particle's z-velocity.
//=================================================================================================
template <int ndim>
void Hydrodynamics<ndim>::CheckZBoundaryGhostParticle
 (const int i,                         ///< [in] i.d. of particles to check
  const FLOAT tghost,                  ///< [in] Expected lifetime of ghost
  const DomainBox<ndim> &simbox)       ///< [in] Simulation domain box
{
  if (ndim == 3) {
    Particle<ndim>& part = GetParticlePointer(i);

    if (part.r[2] + min((FLOAT) 0.0, part.v[2]*tghost) <
        simbox.boxmin[2] + ghost_range*kernrange*part.h) {
      if (simbox.boundary_lhs[2] == periodicBoundary) {
        CreateBoundaryGhostParticle(i, 2, z_lhs_periodic, part.r[2] + simbox.boxsize[2], part.v[2]);
      }
      if (simbox.boundary_lhs[2] == mirrorBoundary) {
        CreateBoundaryGhostParticle(i, 2, z_lhs_mirror,
                                    (FLOAT) 2.0*simbox.boxmin[2] - part.r[2], -part.v[2]);
      }
    }
    if (part.r[2] + max((FLOAT) 0.0, part.v[2]*tghost) >
        simbox.boxmax[2] - ghost_range*kernrange*part.h) {
      if (simbox.boundary_rhs[2] == periodicBoundary) {
        CreateBoundaryGhostParticle(i, 2, z_rhs_periodic, part.r[2] - simbox.boxsize[2], part.v[2]);
      }
      if (simbox.boundary_rhs[2] == mirrorBoundary) {
        CreateBoundaryGhostParticle(i, 2, z_rhs_mirror,
                                    (FLOAT) 2.0*simbox.boxmax[2] - part.r[2], -part.v[2]);
      }
    }
  }

  return;
}



//=================================================================================================
//  Hydrodynamics::CreateBoundaryGhostParticle
/// Create a new ghost particle from either
/// (i) a real Hydrodynamics particle (i < Nhydro), or
/// (ii) an existing ghost particle (i >= Nhydro).
//=================================================================================================
template <int ndim>
void Hydrodynamics<ndim>::CreateBoundaryGhostParticle
 (const int i,                         ///< [in] i.d. of original particle
  const int k,                         ///< [in] Boundary dimension for new ghost
  const int ghosttype,                 ///< [in] Type of ghost particle (periodic, mirror, etc..)
  const FLOAT rk,                      ///< [in] k-position of original particle
  const FLOAT vk)                      ///< [in] k-velocity of original particle
{
  // Increase ghost counter and check there's enough space in memory
  if (Nghost > Nghostmax) {
    cout << "Nghost : " << Nghost << "     Nghostmax : " << Nghostmax << endl;
    string message = "Not enough memory for new ghost";
    ExceptionHandler::getIstance().raise(message);
  }

  int id_new_ghost = Nhydro + Nghost;
  Particle<ndim>& origpart  = GetParticlePointer(i);
  Particle<ndim>& ghostpart = GetParticlePointer(id_new_ghost);

  // If there's enough memory, create ghost particle in arrays
  ghostpart        = origpart;
  ghostpart.r[k]   = rk;
  ghostpart.v[k]   = vk;
  ghostpart.active = false;
  ghostpart.itype  = ghosttype;
  ghostpart.iorig  = i;

  Nghost = Nghost + 1;

  return;
}



template class Hydrodynamics<1>;
template class Hydrodynamics<2>;
template class Hydrodynamics<3>;
