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
//  MeshlessFV::ComputeThermalProperties
/// Compute all thermal properties for grad-h SPH method for given particle.
//=================================================================================================
template <int ndim>
void MeshlessFV<ndim>::ComputeThermalProperties
 (MeshlessFVParticle<ndim> &part)          ///< [inout] Particle i data
{
  // Skip non hydro particles
  if (!types[part.ptype].hydro_forces)
    return ;


  part.u     = eos->SpecificInternalEnergy(part);
  part.sound = eos->SoundSpeed(part);
  part.pressure = eos->Pressure(part);

  assert(part.u > (FLOAT) 0.0);
  assert(part.sound > (FLOAT) 0.0);
  assert(part.pressure > (FLOAT) 0.0);

  return;
}




//
//=================================================================================================
//  MeshlessFV::UpdateArrayVariables
/// Updates all particle quantities based on the primitive/conserved variables.
//=================================================================================================
template <int ndim>
void MeshlessFV<ndim>::UpdateArrayVariables(MeshlessFVParticle<ndim> &part, FLOAT Qcons[nvar])
{

  part.m = Qcons[irho] ;
  part.rho = part.m*part.ndens;
  for (int k=0; k<ndim; k++) part.v[k] = Qcons[k]/part.m;
  assert(isnormal(part.m));

  if (types[part.ptype].hydro_forces) {
    FLOAT ekin = (FLOAT) 0.0;
    for (int k=0; k<ndim; k++) ekin += part.v[k]*part.v[k];

    part.u = (Qcons[ietot] - (FLOAT) 0.5*part.m*ekin)/part.m;
  }

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
    if (part.flags.check(active)) {
      for (int k=0; k<ndim; k++) part.a[k] = 0;
      for (int k=0; k<ndim; k++) part.atree[k] = 0;
      part.gpot = 0 ;
      part.gpot_hydro = 0;
    }
  }
}




template class MeshlessFV<1>;
template class MeshlessFV<2>;
template class MeshlessFV<3>;
