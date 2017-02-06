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
  map<string, string> &stringparams = params->stringparams;

  Nhydromax       = intparams["Nhydromax"];
  create_sinks    = intparams["create_sinks"];
  fixed_sink_mass = intparams["fixed_sink_mass"];
  msink_fixed     = floatparams["m1"];

  timestep_limiter = stringparams["time_step_limiter"] ;

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

    MeshlessFVParticle<ndim>* newhydrodata =
        new struct MeshlessFVParticle<ndim>[N];

    // Swap so that hydrodata points to the new memory
    std::swap(newhydrodata, hydrodata) ;
    if (allocated) {
      // Copy back particle data
      std::copy(newhydrodata,newhydrodata+Nhydromax,hydrodata);
      delete[] newhydrodata;
    }


    Nhydromax=N;
    allocated        = true;
    hydrodata_unsafe = hydrodata;
  }
  assert(Nhydromax >= Nhydro);
  assert(hydrodata);


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
    itype = hydrodata[i].flags.get();
    while (itype & dead) {
      Ndead++;
      ilast--;
      if (i < ilast) {
        hydrodata[i] = hydrodata[ilast];
        hydrodata[ilast].flags.set_flag(dead);
        hydrodata[ilast].m = (FLOAT) 0.0;
      }
      else break;
      itype = hydrodata[i].flags.get();
    };
    if (i >= ilast - 1) break;
  }

  // Reorder all arrays following with new order, with dead particles at end
  if (Ndead == 0) return;

  // Reduce hydro particle counters once dead particles have been removed and reset all
  // other particle counters since a ghost and tree rebuild is required.
  this->NPeriodicGhost = 0;
  this->Nmpighost      = 0;
  Nhydro               -= Ndead;
  Ntot                 = Nhydro;

  // Some sanity checking to ensure there are no dead particles remaining
  for (i=0; i<Nhydro; i++) {
    assert(!hydrodata[i].flags.is_dead());
  }

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
  const FLOAT dt_cfl = 2*courant_mult*part.h/part.vsig_max;
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
  debug2("[MeshlessFV::IntegrateParticles]");
  CodeTiming::BlockTimer timer = timing->StartNewTimer("MFV_INTEGRATE_PARTICLES");

  // Integrate all conserved variables to end of timestep
  //-----------------------------------------------------------------------------------------------
  for (int i=0; i<Nhydro; i++) {
    MeshlessFVParticle<ndim> &part = partdata[i];
    if (part.flags.is_dead()) continue;

    const int dn = n - part.nlast;
    const FLOAT dt = timestep*(FLOAT) dn;
    FLOAT Qcons[nvar];

    if (dn == part.nstep) {
      part.flags.set_flag(active);
      for (int k=0; k<nvar; k++) Qcons[k] = part.Qcons0[k] + part.dQ[k];
    }
    else {
      part.flags.unset_flag(active);
      for (int k=0; k<nvar; k++) Qcons[k] = part.Qcons0[k] + part.dQdt[k]*dt;
    }
    for (int k=0; k<ndim; k++) Qcons[k] += part.Qcons0[irho]*part.a0[k]*dt;


    // Some sanity-checking
    assert(isnormal(Qcons[irho]));
    assert(isnormal(Qcons[ipress]));


    // Compute primitive values and update all main array quantities
    this->UpdateArrayVariables(part, Qcons);
    this->ComputeThermalProperties(part);
    this->UpdatePrimitiveVector(part);


    //---------------------------------------------------------------------------------------------
    if (!staticParticles) {
      part.flags.set_flag(update_density);

      //-------------------------------------------------------------------------------------------
      for (int k=0; k<ndim; k++) {
        part.r[k] = part.r0[k] + (FLOAT) 0.5*(part.v0[k] + part.v[k])*dt;


        // Check if particle has crossed LHS boundary
        //-----------------------------------------------------------------------------------------
        if (part.r[k] < simbox.min[k]) {

          // Check if periodic boundary
          if (simbox.boundary_lhs[k] == periodicBoundary) {
            part.r[k]  += simbox.size[k];
            part.r0[k] += simbox.size[k];
          }

          // Check if wall or mirror boundary
          if (simbox.boundary_lhs[k] == mirrorBoundary || simbox.boundary_lhs[k] == wallBoundary) {
            part.r[k]  = (FLOAT) 2.0*simbox.min[k] - part.r[k];
            part.r0[k] = (FLOAT) 2.0*simbox.min[k] - part.r0[k];
            part.v[k]  = -part.v[k];
            part.v0[k] = -part.v0[k];
            part.a[k]  = -part.a[k];
            part.a0[k] = -part.a0[k];
          }
        }

        // Check if particle has crossed RHS boundary
        //-----------------------------------------------------------------------------------------
        if (part.r[k] > simbox.max[k]) {

          // Check if periodic boundary
          if (simbox.boundary_rhs[k] == periodicBoundary) {
            part.r[k]  -= simbox.size[k];
            part.r0[k] -= simbox.size[k];
          }

          // Check if wall or mirror boundaryq
          if (simbox.boundary_rhs[k] == mirrorBoundary || simbox.boundary_rhs[k] == wallBoundary) {
            part.r[k]  = (FLOAT) 2.0*simbox.max[k] - part.r[k];
            part.r0[k] = (FLOAT) 2.0*simbox.max[k] - part.r0[k];
            part.v[k]  = -part.v[k];
            part.v0[k] = -part.v0[k];
            part.a[k]  = -part.a[k];
            part.a0[k] = -part.a0[k];
          }

        }
        //-----------------------------------------------------------------------------------------

      }
    }
    //---------------------------------------------------------------------------------------------

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
  debug2("[MeshlessFV::EndTimestep]");
  CodeTiming::BlockTimer timer = timing->StartNewTimer("MFV_END_TIMESTEP");


  //-----------------------------------------------------------------------------------------------
#pragma omp parallel for default(none) shared(partdata)
  for (int i=0; i<Npart; i++) {
    MeshlessFVParticle<ndim> &part = partdata[i];    // Local reference to particle
    if (part.flags.is_dead()) continue;

    int dn = n - part.nlast;                         // Integer time since beginning of step
    int k;                                           // Dimension counter
    int nstep = part.nstep;                          // Particle (integer) step size


    // If particle is at the end of its timestep
    //---------------------------------------------------------------------------------------------
    if (dn == nstep) {

      // Integrate all conserved quantities to end of the step (adding sums from neighbours)
      FLOAT Qcons[nvar] ;
      for (int var=0; var<nvar; var++) {
        Qcons[var] = part.Qcons0[var] + part.dQ[var];
        part.dQ[var]    = (FLOAT) 0.0;
        part.dQdt[var]  = (FLOAT) 0.0;
      }

      // Further update conserved quantities if computing gravitational/nbody  contributions
      for (k=0; k<ndim; k++) {
    	Qcons[k] += (FLOAT) 0.5*(FLOAT) dn*timestep*
    			(part.Qcons0[irho]*part.a0[k] + Qcons[irho]*part.a[k]);
        part.v[k] = Qcons[k] / Qcons[irho] ;
      }
      Qcons[ietot] += (FLOAT) 0.5*(FLOAT) dn*timestep*
    	(part.Qcons0[irho]*DotProduct(part.v0, part.a0, ndim) +
    	   Qcons[irho]*DotProduct(part.v, part.a, ndim) +
         DotProduct(part.a0, part.rdmdt0, ndim) +
         DotProduct(part.a, part.rdmdt, ndim));

      // Compute primitive values and update all main array quantities
      this->UpdateArrayVariables(part, Qcons);
      this->ComputeThermalProperties(part);
      this->UpdatePrimitiveVector(part) ;

      // Update all values to the beginning of the next step
      part.nlast  = n;
      part.tlast  = t;
      part.flags.set_flag(active);
      for (k=0; k<ndim; k++) part.r0[k]     = part.r[k];
      for (k=0; k<ndim; k++) part.v0[k]     = part.v[k];
      for (k=0; k<ndim; k++) part.a0[k]     = part.a[k];
      for (k=0; k<ndim; k++) part.rdmdt0[k] = part.rdmdt[k];
      for (k=0; k<ndim; k++) part.rdmdt[k] = 0.0;
      for (k=0; k<nvar; k++) part.Qcons0[k] = Qcons[k];
      for (k=0; k<ndim; k++) part.a[k] = 0.0;
      part.gpot=0.0;

      for (k=0; k<ndim; k++) part.rdmdt[k] = (FLOAT) 0.0;

      for (k=0; k<ndim; k++) part.rdmdt[k] = (FLOAT) 0.0;

    }
    //---------------------------------------------------------------------------------------------
    else {
      part.flags.unset_flag(active);
    }
    //---------------------------------------------------------------------------------------------

  }
  //-----------------------------------------------------------------------------------------------

  //timing->EndTimingSection("MFV_END_TIMESTEP");

  return;
}

//=================================================================================================
//  MeshlessFV<ndim>::EndTimestep
/// Calculate or reset all quantities for all particles that reach the end of their timesteps.
//=================================================================================================
template <int ndim>
int MeshlessFV<ndim>::CheckTimesteps
(const int level_diff_max,            ///< [in] Max. allowed SPH neib dt diff
 const int level_step,                ///< [in] Level of base timestep
 const int n,                         ///< [in] Integer time in block time struct
 double timestep,                     ///< [in] Timestep
 int mode_)
 {
  int dn;                              // Integer time since beginning of step
  int level_new;                       // New timestep level
  int nnewstep;                        // New integer timestep
  int activecount = 0;                 // No. of newly active particles
  int i;                               // Particle counter


  MeshlessFVParticle<ndim> *mfvdata = GetMeshlessFVParticleArray() ;

  debug2("[MeshlessFV::CheckTimesteps]");
  CodeTiming::BlockTimer timer = timing->StartNewTimer("MESHLESS_CHECK_TIMESTEPS");

  const int mode = mode_ ;
  if (mode == 1 && (timestep_limiter != "simple")) return 0 ;

  const double tstep = timestep ;
  //-----------------------------------------------------------------------------------------------
  #pragma omp parallel for default(none) private(dn,i,level_new,nnewstep) \
    shared(mfvdata,cout) reduction(+:activecount)
  for (i=0; i<Nhydro; i++) {
    MeshlessFVParticle<ndim>& part = mfvdata[i];
    if (part.flags.is_dead()) continue;

    dn = n - part.nlast;

    // Check if neighbour timesteps are too small.  If so, then reduce timestep if possible
    if (part.levelneib - part.level > level_diff_max) {
      level_new = part.levelneib - level_diff_max;
      nnewstep  = pow(2,level_step - level_new);

      // Force recalculation of fluxes for particles at the end of their step
      if (mode == 0) {
        if(dn == 0) {
          part.nstep = nnewstep ;
          part.level = level_new;
          part.flags.set_flag(active);
          activecount++;
        }
      }
      // Saitoh & Makino type reactive limiting
      //   Note: The current timestep ends at dn+1
      else if (mode == 1) {
        if(dn%nnewstep == 0 && dn != part.nstep) {
          part.nstep = dn;
          part.level = level_new;
          double dt = part.nstep * tstep ;

          // Use current predicted value for dQ
          for (int var=0; var<nvar; var++)
            part.dQ[var] = dt * part.dQdt[var];

          part.flags.set_flag(active);
        }
      }
    }
  }
    //-----------------------------------------------------------------------------------------------

  return activecount;
 }

//=================================================================================================
//  MeshlessFV::UpdateArrayVariables
/// Updates all particle quantities based on the primitive/conserved variables.
//=================================================================================================
template <int ndim>
void MeshlessFV<ndim>::UpdateArrayVariables(MeshlessFVParticle<ndim> &part, FLOAT Qcons[nvar])
{
  // TODO: Check all callers.
  //   This now uses the currently predicted value of Qcons, no need to add dQ.
  part.m = Qcons[irho] ;
  part.rho = part.m*part.ndens;
  for (int k=0; k<ndim; k++) part.v[k] = Qcons[k]/part.m;

  FLOAT ekin = (FLOAT) 0.0;
  for (int k=0; k<ndim; k++) ekin += part.v[k]*part.v[k];
  part.u = (Qcons[ietot] - (FLOAT) 0.5*part.m*ekin)/part.m;
  part.u = eos->SpecificInternalEnergy(part);
  part.press = (gamma_eos - (FLOAT) 1.0)*part.rho*part.u;

  assert(isnormal(part.m));
  assert(isnormal(part.u));
  assert(isnormal(part.press));
  /*assert(part.m > (FLOAT) 0.0);
  assert(part.u > (FLOAT) 0.0);
  assert(part.press > (FLOAT) 0.0);*/

  return;
}



//=================================================================================================
//  MfvMuscl::InitialSmoothingLengthGuess
/// Perform initial guess of smoothing.  In the absence of more sophisticated techniques, we guess
/// the smoothing length assuming a uniform density medium with the same volume and total mass.
//=================================================================================================
template<>
void MeshlessFV<1>::InitialSmoothingLengthGuess(void)
{
  int i;                           // Particle counter
  FLOAT h_guess;                   // Global guess of smoothing length
  FLOAT volume;                    // Volume of global bounding box
  FLOAT rmin[1];                   // Min. extent of bounding box
  FLOAT rmax[1];                   // Max. extent of bounding box

  debug2("[MeshlessFV::InitialSmoothingLengthGuess]");

  // Calculate bounding box containing all SPH particles
  this->ComputeBoundingBox(rmax,rmin,Nhydro);

  // Depending on the dimensionality, calculate the average smoothing
  // length assuming a uniform density distribution filling the bounding box.
  //-----------------------------------------------------------------------------------------------
  Ngather = (int) ((FLOAT) 2.0*kernp->kernrange*h_fac);
  volume = rmax[0] - rmin[0];
  h_guess = (volume*(FLOAT) Ngather)/((FLOAT) 4.0*(FLOAT) Nhydro);
  //-----------------------------------------------------------------------------------------------

  // Set all smoothing lengths equal to average value
  for (i=0; i<Nhydro; i++) {
    MeshlessFVParticle<1>& part = GetMeshlessFVParticlePointer(i);
    part.h         = h_guess;
    part.hrangesqd = kernp->kernrangesqd*part.h*part.h;
  }

  return;
}
template<>
void MeshlessFV<2>::InitialSmoothingLengthGuess(void)
{
  int i;                           // Particle counter
  FLOAT h_guess;                   // Global guess of smoothing length
  FLOAT volume;                    // Volume of global bounding box
  FLOAT rmin[2];                   // Min. extent of bounding box
  FLOAT rmax[2];                   // Max. extent of bounding box

  debug2("[MeshlessFV::InitialSmoothingLengthGuess]");

  // Calculate bounding box containing all SPH particles
  this->ComputeBoundingBox(rmax,rmin,Nhydro);

  // Depending on the dimensionality, calculate the average smoothing
  // length assuming a uniform density distribution filling the bounding box.
  //-----------------------------------------------------------------------------------------------
  Ngather = (int) (pi*pow(kernp->kernrange*h_fac,2));
  volume = (rmax[0] - rmin[0])*(rmax[1] - rmin[1]);
  h_guess = sqrtf((volume*(FLOAT) Ngather)/((FLOAT) 4.0*(FLOAT) Nhydro));
  //-----------------------------------------------------------------------------------------------

  // Set all smoothing lengths equal to average value
  for (i=0; i<Nhydro; i++) {
    MeshlessFVParticle<2>& part = GetMeshlessFVParticlePointer(i);
    part.h         = h_guess;
    part.hrangesqd = kernp->kernrangesqd*part.h*part.h;
  }

  return;
}
template <>
void MeshlessFV<3>::InitialSmoothingLengthGuess(void)
{
  int i;                           // Particle counter
  FLOAT h_guess;                   // Global guess of smoothing length
  FLOAT volume;                    // Volume of global bounding box
  FLOAT rmin[3];                   // Min. extent of bounding box
  FLOAT rmax[3];                   // Max. extent of bounding box

  debug2("[MeshlessFV::InitialSmoothingLengthGuess]");

  // Calculate bounding box containing all SPH particles
  this->ComputeBoundingBox(rmax,rmin,Nhydro);

  // Depending on the dimensionality, calculate the average smoothing
  // length assuming a uniform density distribution filling the bounding box.
  //-----------------------------------------------------------------------------------------------
  Ngather = (int) ((FLOAT) 4.0*pi*pow(kernp->kernrange*h_fac,3)/(FLOAT) 3.0);
  volume = (rmax[0] - rmin[0])*(rmax[1] - rmin[1])*(rmax[2] - rmin[2]);
  h_guess = powf((3.0*volume*(FLOAT) Ngather)/(32.0*pi*(FLOAT) Nhydro),onethird);
  //-----------------------------------------------------------------------------------------------

  // Set all smoothing lengths equal to average value
  for (i=0; i<Nhydro; i++) {
    MeshlessFVParticle<3>& part = GetMeshlessFVParticlePointer(i);
    part.h         = h_guess;
    part.hrangesqd = kernp->kernrangesqd*part.h*part.h;
  }

  return;
}

//=================================================================================================
//  MeshlessFV::ZeroAccelerations
/// Initialise key variables before force calculations
//=================================================================================================
template <int ndim>
void MeshlessFV<ndim>::ZeroAccelerations()
{
  for (int i=0; i<Nhydro; i++) {
    MeshlessFVParticle<ndim>& part = GetMeshlessFVParticlePointer(i);
    if (part.flags.check_flag(active)) {
      for (int k=0; k<ndim; k++) part.a[k] = 0;
      for (int k=0; k<ndim; k++) part.atree[k] = 0;
      part.gpot = 0 ;
    }
  }
}




template class MeshlessFV<1>;
template class MeshlessFV<2>;
template class MeshlessFV<3>;
