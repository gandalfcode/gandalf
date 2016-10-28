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
                                   string _gas_eos, string KernelName, int size_hydro,
                                   SimUnits &units, Parameters *params):
  size_hydro_part(size_hydro),
  hydro_forces(hydro_forces_aux),
  self_gravity(self_gravity_aux),
  gas_eos(_gas_eos),
  h_fac(h_fac_aux),
  types(params),
  kerntab(TabulatedKernel<ndim>(KernelName))
{
  // Local references to parameter variables for brevity
  map<string, int> &intparams = params->intparams;
  map<string, double> &floatparams = params->floatparams;
  map<string, string> &stringparams = params->stringparams;
  string gas_radiation = stringparams["radiation"];


  // Zero or initialise all common Hydrodynamics variables
  //-----------------------------------------------------------------------------------------------
  allocated          = false;
  Nhydro             = 0;
  Nhydromax          = 0;
  NImportedParticles = 0;
  Nmpighost          = 0;
  NPeriodicGhost     = 0;
  create_sinks       = intparams["create_sinks"];

  // Select and construct equation of state object from given parameters
  //-----------------------------------------------------------------------------------------------
  if ((_gas_eos == "energy_eqn" || _gas_eos == "constant_temp" ||
       _gas_eos == "isothermal" || _gas_eos == "barotropic" ||
       _gas_eos == "barotropic2") && gas_radiation == "ionisation") {
    eos = new IonisingRadiation<ndim>
      (_gas_eos, floatparams["temp0"], floatparams["mu_bar"],
       floatparams["gamma_eos"], floatparams["rho_bary"], &units);
  }
  else if ((_gas_eos == "energy_eqn" || _gas_eos == "constant_temp" ||
            _gas_eos == "isothermal" || _gas_eos == "barotropic" ||
            _gas_eos == "barotropic2") && gas_radiation == "monoionisation") {
    eos = new MCRadiationEOS<ndim>
      (_gas_eos, floatparams["temp0"], floatparams["temp_ion"], floatparams["mu_bar"],
       floatparams["mu_ion"], floatparams["gamma_eos"], floatparams["rho_bary"], &units);
  }
  else if (_gas_eos == "energy_eqn" || _gas_eos == "constant_temp") {
    eos = new Adiabatic<ndim>
      (floatparams["temp0"], floatparams["mu_bar"], floatparams["gamma_eos"]);
  }
  else if (_gas_eos == "isothermal") {
    eos = new Isothermal<ndim>
      (floatparams["temp0"], floatparams["mu_bar"], floatparams["gamma_eos"], &units);
  }
  else if (_gas_eos == "barotropic") {
    eos = new Barotropic<ndim>(floatparams["temp0"], floatparams["mu_bar"],
                               floatparams["gamma_eos"], floatparams["rho_bary"], &units);
  }
  else if (_gas_eos == "barotropic2") {
    eos = new Barotropic2<ndim>(floatparams["temp0"], floatparams["mu_bar"],
                                     floatparams["gamma_eos"], floatparams["rho_bary"], &units);
  }
  else if (_gas_eos == "rad_ws" || _gas_eos == "radws") {
    eos = new Radws<ndim>(floatparams["temp0"], floatparams["mu_bar"], floatparams["gamma_eos"]);
  }
  else {
    string message = "Unrecognised parameter : gas_eos = " + _gas_eos;
    ExceptionHandler::getIstance().raise(message);
  }

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
  const FLOAT r0 = part.r[0];
  const FLOAT v0 = part.v[0];
  const FLOAT h = part.h;

  if (r0 + min((FLOAT) 0.0, v0*tghost) <
      simbox.boxmin[0] + ghost_range*kernrange*h) {
    if (simbox.boundary_lhs[0] == periodicBoundary) {
      CreateBoundaryGhostParticle(i, 0, x_periodic_lhs,r0 + simbox.boxsize[0], v0);
    }
    if (simbox.boundary_lhs[0] == mirrorBoundary) {
      CreateBoundaryGhostParticle(i, 0, x_mirror_lhs,
                                  (FLOAT) 2.0*simbox.boxmin[0] - r0, -v0);
    }
  }
  if (r0 + max((FLOAT) 0.0, v0*tghost) >
      simbox.boxmax[0] - ghost_range*kernrange*h) {
    if (simbox.boundary_rhs[0] == periodicBoundary) {
      CreateBoundaryGhostParticle(i, 0, x_periodic_rhs, r0 - simbox.boxsize[0], v0);
    }
    if (simbox.boundary_rhs[0] == mirrorBoundary) {
      CreateBoundaryGhostParticle(i, 0, x_mirror_rhs,
                                  (FLOAT) 2.0*simbox.boxmax[0] - r0, -v0);
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

	assert(ndim>1);

	Particle<ndim>& part = GetParticlePointer(i);
	const FLOAT r1 = part.r[1];
	const FLOAT v1 = part.v[1];
	const FLOAT h = part.h;


  if (ndim > 1) {

    if (r1 + min((FLOAT) 0.0, v1*tghost) <
        simbox.boxmin[1] + ghost_range*kernrange*h) {
      if (simbox.boundary_lhs[1] == periodicBoundary) {
        CreateBoundaryGhostParticle(i, 1, y_periodic_lhs, r1 + simbox.boxsize[1], v1);
      }
      if (simbox.boundary_lhs[1] == mirrorBoundary) {
        CreateBoundaryGhostParticle(i, 1, y_mirror_lhs,
                                    (FLOAT) 2.0*simbox.boxmin[1] - r1, -v1);
      }
    }
    if (r1 + max((FLOAT) 0.0, v1*tghost) >
        simbox.boxmax[1] - ghost_range*kernrange*h) {
      if (simbox.boundary_rhs[1] == periodicBoundary) {
        CreateBoundaryGhostParticle(i, 1, y_periodic_rhs, r1 - simbox.boxsize[1], v1);
      }
      if (simbox.boundary_rhs[1] == mirrorBoundary) {
        CreateBoundaryGhostParticle(i, 1, y_mirror_rhs,
                                    (FLOAT) 2.0*simbox.boxmax[1] - r1, -r1);
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
	assert(ndim>2);

	Particle<ndim>& part = GetParticlePointer(i);
	const FLOAT r2 = part.r[2];
	const FLOAT v2 = part.v[2];
	const FLOAT h = part.h;

  if (ndim == 3) {

    if (r2 + min((FLOAT) 0.0, v2*tghost) <
        simbox.boxmin[2] + ghost_range*kernrange*h) {
      if (simbox.boundary_lhs[2] == periodicBoundary) {
        CreateBoundaryGhostParticle(i, 2, z_periodic_lhs, r2 + simbox.boxsize[2], v2);
      }
      if (simbox.boundary_lhs[2] == mirrorBoundary) {
        CreateBoundaryGhostParticle(i, 2, z_mirror_lhs,
                                    (FLOAT) 2.0*simbox.boxmin[2] - r2, -v2);
      }
    }
    if (r2 + max((FLOAT) 0.0, v2*tghost) >
        simbox.boxmax[2] - ghost_range*kernrange*h) {
      if (simbox.boundary_rhs[2] == periodicBoundary) {
        CreateBoundaryGhostParticle(i, 2, z_periodic_rhs, r2 - simbox.boxsize[2],v2);
      }
      if (simbox.boundary_rhs[2] == mirrorBoundary) {
        CreateBoundaryGhostParticle(i, 2, z_mirror_rhs,
                                    (FLOAT) 2.0*simbox.boxmax[2] - r2, -v2);
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
  // Reallocate memory if necessary
  if (Nhydro+Nghost >= Nhydromax) {
	  //TODO: should create a layer on top of AllocateMemory, it should not be called directly
	  AllocateMemory(1.1*Nhydromax+1);
  }

  int id_new_ghost = Nhydro + Nghost;
  Particle<ndim>& origpart  = GetParticlePointer(i);
  Particle<ndim>& ghostpart = GetParticlePointer(id_new_ghost);

  // If there's enough memory, create ghost particle in arrays
  ghostpart        = origpart;
  ghostpart.r[k]   = rk;
  ghostpart.v[k]   = vk;
  ghostpart.flags.unset_flag(active);
  ghostpart.flags.set_flag(ghosttype); // Allow ghost to have multiple ghost flags
  ghostpart.iorig  = i;

  Nghost++;

  return;
}



template class Hydrodynamics<1>;
template class Hydrodynamics<2>;
template class Hydrodynamics<3>;
