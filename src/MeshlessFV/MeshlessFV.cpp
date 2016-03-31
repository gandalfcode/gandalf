//=================================================================================================
//  MeshlessFV.cpp
//  Contains all functions for calculating Meshless Finite-Volume Hydrodynamics quantities.
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


#include <assert.h>
#include <cstdio>
#include <cstring>
#include <cstdlib>
//#include <cassert>
#include <iostream>
#include <math.h>
#include "Precision.h"
#include "MeshlessFV.h"
#include "Particle.h"
#include "Parameters.h"
#include "SmoothingKernel.h"
#include "EOS.h"
#include "Debug.h"
#include "Exception.h"
#include "InlineFuncs.h"
using namespace std;



//=================================================================================================
//  MeshlessFV::MeshlessFV
/// MeshlessFV class constructor.  Calls main SPH class constructor and also
/// sets additional kernel-related quantities
//=================================================================================================
template <int ndim>
MeshlessFV<ndim>::MeshlessFV(int _hydro_forces, int _self_gravity, FLOAT _accel_mult,
                             FLOAT _courant_mult, FLOAT _h_fac, FLOAT _h_converge,
                             FLOAT _gamma, string _gas_eos, string KernelName, int size_part,
                             SimUnits &units, Parameters *params):
  FV<ndim>(_hydro_forces, _self_gravity, _accel_mult, _courant_mult, _h_fac,
           _h_converge, _gamma, _gas_eos, KernelName, size_part, units, params),
  accel_mult(_accel_mult),
  courant_mult(_courant_mult),
  h_converge(_h_converge),
  staticParticles(params->intparams["static_particles"])
{
  // Local references to parameter variables for brevity
  map<string, int> &intparams = params->intparams;
  map<string, double> &floatparams = params->floatparams;
  //map<string, string> &stringparams = params->stringparams;

  Nhydromax       = intparams["Nhydromax"];
  create_sinks    = intparams["create_sinks"];
  fixed_sink_mass = intparams["fixed_sink_mass"];
  msink_fixed     = floatparams["m1"];
}



//=================================================================================================
//  MeshlessFV::~MeshlessFV
/// MeshlessFV class destructor
//=================================================================================================
template <int ndim>
MeshlessFV<ndim>::~MeshlessFV()
{
  //DeallocateMemory();
}



//=================================================================================================
//  MeshlessFV::AllocateMemory
/// Allocate main SPH particle array.  Estimates the maximum number of boundary ghost particles
/// assuming a roughly uniform depth of ghosts at each boundary.
//=================================================================================================
template <int ndim>
void MeshlessFV<ndim>::AllocateMemory(int N)
{
  debug2("[MeshlessFV::AllocateMemory]");

  if (N > Nhydromax || !allocated) {
    if (allocated) DeallocateMemory();

    // Set conservative estimate for maximum number of particles, assuming
    // extra space required for (periodic) ghost particles
    if (Nhydromax < N) {
      Nhydromax = 2*(int) powf(powf((FLOAT) N,invndim) + (FLOAT) 16.0*kernp->kernrange,ndim);
    }

    iorder    = new int[Nhydromax];
    hydrodata = new struct MeshlessFVParticle<ndim>[Nhydromax];
    allocated = true;
    hydrodata_unsafe = hydrodata;
  }

  return;
}



//=================================================================================================
//  MeshlessFV::DeallocateMemory
/// Deallocate main array containing SPH particle data.
//=================================================================================================
template <int ndim>
void MeshlessFV<ndim>::DeallocateMemory(void)
{
  debug2("[MeshlessFV::DeallocateMemory]");

  if (allocated) {
    delete[] hydrodata;
    delete[] iorder;
  }
  allocated = false;

  return;
}



//=================================================================================================
//  MeshlessFV::DeleteDeadParticles
/// Delete 'dead' (e.g. accreted) SPH particles from the main arrays.
//=================================================================================================
template <int ndim>
void MeshlessFV<ndim>::DeleteDeadParticles(void)
{
  int i;                               // Particle counter
  int itype;                           // Current particle type
  int Ndead = 0;                       // No. of 'dead' particles
  int ilast = Nhydro;                  // Aux. counter of last free slot

  debug2("[MeshlessFV::DeleteDeadParticles]");


  // Determine new order of particles in arrays.
  // First all live particles and then all dead particles.
  for (i=0; i<Nhydro; i++) {
    itype = hydrodata[i].itype;
    while (itype == dead) {
      Ndead++;
      ilast--;
      if (i < ilast) {
        hydrodata[i] = hydrodata[ilast];
        hydrodata[ilast].itype = dead;
        hydrodata[ilast].m = (FLOAT) 0.0;
      }
      else break;
      itype = hydrodata[i].itype;
    };
    if (i >= ilast - 1) break;
  }

  // Reorder all arrays following with new order, with dead particles at end
  if (Ndead == 0) return;

  // Reduce particle counters once dead particles have been removed
  Nhydro -= Ndead;
  Ntot -= Ndead;
  for (i=0; i<Nhydro; i++) {
    iorder[i] = i;
    assert(hydrodata[i].itype != dead);
  }

  return;
}



//=================================================================================================
//  MeshlessFV::ReorderParticles
/// Reorder particles in main arrays according to order given in iorder array.
//=================================================================================================
template <int ndim>
void MeshlessFV<ndim>::ReorderParticles(void)
{
  int i;                                   // Particle counter
  MeshlessFVParticle<ndim> *hydrodataaux;  // Aux. SPH particle array

  hydrodataaux = new MeshlessFVParticle<ndim>[Nhydro];

  for (i=0; i<Nhydro; i++) hydrodataaux[i] = hydrodata[i];
  for (i=0; i<Nhydro; i++) hydrodata[i] = hydrodataaux[iorder[i]];

  delete[] hydrodataaux;

  return;
}



//=================================================================================================
//  MeshlessFV::ComputeThermalProperties
/// Compute all thermal properties for grad-h SPH method for given particle.
//=================================================================================================
template <int ndim>
void MeshlessFV<ndim>::ComputeThermalProperties
 (MeshlessFVParticle<ndim> &part)          ///< [inout] Particle i data
{
  part.u     = eos->SpecificInternalEnergy(part);
  part.sound = eos->SoundSpeed(part);
  part.press = eos->Pressure(part);

  FLOAT ekin = (FLOAT) 0.0;
  for (int k=0; k<ndim; k++) ekin += part.v[k]*part.v[k];
  part.Qcons[ietot] = part.m*part.u + (FLOAT) 0.5*part.m*ekin;
  //part.Qcons[ietot] = part.press*part.volume/(gamma_eos - (FLOAT) 1.0) +
  //  (FLOAT) 0.5*part.rho*part.volume*ekin;

  assert(part.u > (FLOAT) 0.0);
  assert(part.sound > (FLOAT) 0.0);
  assert(part.press > (FLOAT) 0.0);

  return;
}



//=================================================================================================
//  MeshlessFV<ndim>::Timestep
/// Compute timestep for particle based on Courant and acceleration conditions.
//=================================================================================================
template <int ndim>
FLOAT MeshlessFV<ndim>::Timestep(MeshlessFVParticle<ndim> &part)
{
  const FLOAT dt_cfl = courant_mult*part.h/part.vsig_max;
  const FLOAT dt_grav = accel_mult*
    sqrtf(part.h/sqrt(DotProduct(part.a0, part.a0, ndim) + small_number));

  if (hydro_forces && self_gravity) return min(dt_cfl, dt_grav);
  else if (hydro_forces) return dt_cfl;
  else if (self_gravity) return dt_grav;
  else return big_number;
}



//=================================================================================================
//  MeshlessFV<ndim>::IntegrateParticles
/// Calculate or reset all quantities for all particles that reach the end of their timesteps.
//=================================================================================================
template <int ndim>
void MeshlessFV<ndim>::IntegrateParticles
 (const int n,                         ///< [in] Integer time in block time struct
  const int Npart,                     ///< [in] Number of particles
  const FLOAT t,                       ///< [in] Current simulation time
  const FLOAT timestep,                ///< [in] Base timestep value
  const DomainBox<ndim> &simbox,       ///< [in] Simulation box
  MeshlessFVParticle<ndim> *partdata)  ///< [inout] Pointer to SPH particle array
{
  int dn;                              // Integer time since beginning of step
  int i;                               // Particle counter
  int k;                               // Dimension counter

  debug2("[MeshlessFV::IntegrateParticles]");

  // Integrate all conserved variables to end of timestep
  //-----------------------------------------------------------------------------------------------
  if (!staticParticles) {
    for (i=0; i<Nhydro; i++) {
      MeshlessFVParticle<ndim> &part = partdata[i];
      dn = n - part.nlast;
      part.m = part.Qcons[FV<ndim>::irho] + part.dQ[FV<ndim>::irho];

      //-------------------------------------------------------------------------------------------
      for (k=0; k<ndim; k++) {
        part.r[k] = part.r0[k] + part.v0[k]*timestep*(FLOAT) dn;

        // Check if particle has crossed LHS boundary
        //-----------------------------------------------------------------------------------------
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
        //-----------------------------------------------------------------------------------------
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
        //-----------------------------------------------------------------------------------------

      }
      //-------------------------------------------------------------------------------------------

    }
  }
  //-----------------------------------------------------------------------------------------------


  return;
}



//=================================================================================================
//  MeshlessFV<ndim>::EndTimestep
/// Calculate or reset all quantities for all particles that reach the end of their timesteps.
//=================================================================================================
template <int ndim>
void MeshlessFV<ndim>::EndTimestep
 (const int n,                         ///< [in] Integer time in block time struct
  const int Npart,                     ///< [in] Number of particles
  const FLOAT t,                       ///< [in] Current simulation time
  const FLOAT timestep,                ///< [in] Base timestep value
  MeshlessFVParticle<ndim> *partdata)  ///< [inout] Pointer to SPH particle array
{
  int i;                               // Particle counter

  debug2("[MeshlessFV::EndTimestep]");
  //timing->StartTimingSection("MFV_END_TIMESTEP");


  //-----------------------------------------------------------------------------------------------
#pragma omp parallel for default(none) shared(partdata)
  for (i=0; i<Npart; i++) {

    MeshlessFVParticle<ndim> &part = partdata[i];    // Local reference to particle
    int dn = n - part.nlast;                         // Integer time since beginning of step
    int k;                                           // Dimension counter
    int nstep = part.nstep;                          // Particle (integer) step size

    if (part.itype == dead) continue;


    // If particle is at the end of its timestep
    //---------------------------------------------------------------------------------------------
    if (dn == nstep) {

      // Integrate all conserved quantities to end of the step (adding sums from neighbours)
      for (int var=0; var<nvar; var++) {
        part.Qcons[var] += part.dQ[var];
        part.dQ[var]     = (FLOAT) 0.0;
      }

      // Further update conserved quantities if computing gravitational contributions
      if (self_gravity == 1) {
        for (k=0; k<ndim; k++) part.Qcons[k] += (FLOAT) 0.5*(FLOAT) dn*timestep*
          (part.Qcons0[irho]*part.a0[k] + part.Qcons[irho]*part.a[k]);
        part.Qcons[ietot] += (FLOAT) 0.5*(FLOAT) dn*timestep*
          (part.Qcons0[irho]*DotProduct(part.v0, part.a0, ndim) +
           part.Qcons[irho]*DotProduct(part.v, part.a, ndim) +
           DotProduct(part.a0, part.rdmdt0, ndim) +
           DotProduct(part.a, part.rdmdt, ndim));
        //part.Qcons[ietot] += (FLOAT) 0.5*(FLOAT) dn*timestep*
        //  (part.Qcons0[irho]*DotProduct(part.v0, part.a0, ndim) +
        //   part.Qcons[irho]*DotProduct(part.v, part.a, ndim));
      }

      // Compute primtive values and update all main array quantities
      this->ConvertConservedToPrimitive(part.volume, part.Qcons, part.Wprim);
      this->UpdateArrayVariables(part);

      // Update all values to the beginning of the next step
      part.nlast  = n;
      part.tlast  = t;
      part.active = true;
      for (k=0; k<ndim; k++) part.r0[k]     = part.r[k];
      for (k=0; k<ndim; k++) part.v0[k]     = part.v[k];
      for (k=0; k<ndim; k++) part.a0[k]     = part.a[k];
      for (k=0; k<ndim; k++) part.rdmdt0[k] = part.rdmdt[k];
      for (k=0; k<nvar; k++) part.Qcons0[k] = part.Qcons[k];

    }
    //---------------------------------------------------------------------------------------------
    else {
      part.active = false;
    }
    //---------------------------------------------------------------------------------------------

  }
  //-----------------------------------------------------------------------------------------------

  //timing->EndTimingSection("MFV_END_TIMESTEP");

  return;
}



//=================================================================================================
//  MeshlessFV::UpdatePrimitiveVector
/// Updates the primitive vector from particle quantities.
//=================================================================================================
template <int ndim>
void MeshlessFV<ndim>::UpdatePrimitiveVector(MeshlessFVParticle<ndim> &part)
{
  for (int k=0; k<ndim; k++) part.Wprim[k] = part.v[k];
  part.Wprim[irho] = part.rho;
  part.Wprim[ipress] = part.press;
}



//=================================================================================================
//  MeshlessFV::UpdateArrayVariables
/// Updates all particle quantities based on the primmitive/conserved variables.
//=================================================================================================
template <int ndim>
void MeshlessFV<ndim>::UpdateArrayVariables(MeshlessFVParticle<ndim> &part)
{
  part.m = part.Qcons[irho] + part.dQ[irho];
  part.rho = part.m/part.volume;
  for (int k=0; k<ndim; k++) part.v[k] = (part.Qcons[k] + part.dQ[k])/part.m;

  FLOAT ekin = (FLOAT) 0.0;
  for (int k=0; k<ndim; k++) ekin += part.v[k]*part.v[k];
  part.u = (part.Qcons[ietot] + part.dQ[ietot] - (FLOAT) 0.5*part.m*ekin)/part.m;
  part.u = eos->SpecificInternalEnergy(part);
  part.press = (gamma_eos - (FLOAT) 1.0)*part.rho*part.u;


  assert(part.m > (FLOAT) 0.0);
  assert(part.u > (FLOAT) 0.0);
  assert(part.press > (FLOAT) 0.0);

}



//=================================================================================================
//  MfvMuscl::InitialSmoothingLengthGuess
/// Perform initial guess of smoothing.  In the abscence of more sophisticated techniques, we guess
/// the smoothing length assuming a uniform density medium with the same volume and total mass.
//=================================================================================================
template <int ndim>
void MeshlessFV<ndim>::InitialSmoothingLengthGuess(void)
{
  int i;                           // Particle counter
  FLOAT h_guess;                   // Global guess of smoothing length
  FLOAT volume;                    // Volume of global bounding box
  FLOAT rmin[ndim];                // Min. extent of bounding box
  FLOAT rmax[ndim];                // Max. extent of bounding box

  debug2("[MeshlessFV::InitialSmoothingLengthGuess]");

  // Calculate bounding box containing all SPH particles
  this->ComputeBoundingBox(rmax,rmin,Nhydro);

  // Depending on the dimensionality, calculate the average smoothing
  // length assuming a uniform density distribution filling the bounding box.
  //-----------------------------------------------------------------------------------------------
  if (ndim == 1) {
    Ngather = (int) ((FLOAT) 2.0*kernp->kernrange*h_fac);
    volume = rmax[0] - rmin[0];
    h_guess = (volume*(FLOAT) Ngather)/((FLOAT) 4.0*(FLOAT) Nhydro);
  }
  //-----------------------------------------------------------------------------------------------
  else if (ndim == 2) {
    Ngather = (int) (pi*pow(kernp->kernrange*h_fac,2));
    volume = (rmax[0] - rmin[0])*(rmax[1] - rmin[1]);
    h_guess = sqrtf((volume*(FLOAT) Ngather)/((FLOAT) 4.0*(FLOAT) Nhydro));
  }
  //-----------------------------------------------------------------------------------------------
  else if (ndim == 3) {
    Ngather = (int) ((FLOAT) 4.0*pi*pow(kernp->kernrange*h_fac,3)/(FLOAT) 3.0);
    volume = (rmax[0] - rmin[0])*(rmax[1] - rmin[1])*(rmax[2] - rmin[2]);
    h_guess = powf((3.0*volume*(FLOAT) Ngather)/(32.0*pi*(FLOAT) Nhydro),onethird);
  }
  //-----------------------------------------------------------------------------------------------

  // Set all smoothing lengths equal to average value
  for (i=0; i<Nhydro; i++) {
    MeshlessFVParticle<ndim>& part = GetMeshlessFVParticlePointer(i);
    part.h         = h_guess;
    part.invh      = (FLOAT) 1.0/h_guess;
    part.hrangesqd = kernfacsqd*kernp->kernrangesqd*part.h*part.h;
  }

  return;
}



template class MeshlessFV<1>;
template class MeshlessFV<2>;
template class MeshlessFV<3>;
