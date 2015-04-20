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


#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cassert>
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
MeshlessFV<ndim>::MeshlessFV(int hydro_forces_aux, int self_gravity_aux, FLOAT _accel_mult,
                             FLOAT _courant_mult, FLOAT h_fac_aux, FLOAT h_converge_aux,
                             FLOAT gamma_aux, string gas_eos_aux, string KernelName, int size_part):
  FV<ndim>(hydro_forces_aux, self_gravity_aux, _accel_mult, _courant_mult, h_fac_aux,
           h_converge_aux, gamma_aux, gas_eos_aux, KernelName, size_part),
  accel_mult(_accel_mult),
  courant_mult(_courant_mult),
  h_converge(h_converge_aux)
  //gamma_eos(gamma_aux),
  //gammam1(gamma_aux - 1.0)
  //size_hydro_part(size_part)
{
  /*this->kernp      = &kern;
  this->kernfac    = (FLOAT) 1.0;
  this->kernfacsqd = (FLOAT) 1.0;
  this->kernrange  = this->kernp->kernrange;*/
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
    // extra space required for periodic ghost particles
    if (Nhydromax < N) {
      Nhydromax = 2*(int) powf(powf((FLOAT) N,invndim) + (FLOAT) 8.0*kernp->kernrange,ndim);
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
/// Delete selected SPH particles from the main arrays.
//=================================================================================================
template <int ndim>
void MeshlessFV<ndim>::ReorderParticles(void)
{
  int i;                               // Particle counter
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
 (MeshlessFVParticle<ndim> &part)         ///< [inout] Particle i data
{
  //MeshlessFVParticle<ndim>& part = static_cast<MeshlessFVParticle<ndim> &> (part_gen);

  part.u     = eos->SpecificInternalEnergy(part);
  part.sound = eos->SoundSpeed(part);
  part.press = eos->Pressure(part);

  assert(part.u > 0.0);
  assert(part.sound > 0.0);
  assert(part.press > 0.0);

  return;
}



//=================================================================================================
//  MeshlessFV<ndim>::Timestep
/// ...
//=================================================================================================
template <int ndim>
FLOAT MeshlessFV<ndim>::Timestep(MeshlessFVParticle<ndim> &part)
{
  //cout << "Timestep : " << 0.5*part.h/(part.sound + sqrt(DotProduct(part.v, part.v, ndim)))
  //     << "    " << part.sound << "    " << part.h <<"    "
  //      << part.v[0] << endl;
  //cout << "Timestep : " << courant_mult*part.h/part.vsig_max << "   " << courant_mult << "   " << part.h << "    " << part.vsig_max << endl;
  return courant_mult*part.h/part.vsig_max;
  return courant_mult*part.h/(part.sound + sqrt(DotProduct(part.v, part.v, ndim)));
}



//=================================================================================================
//  MeshlessFV<ndim>::Timestep
/// ...
//=================================================================================================
template <int ndim>
void MeshlessFV<ndim>::IntegrateConservedVariables
 (MeshlessFVParticle<ndim> &part,
  FLOAT timestep)
{
  FLOAT dUdt = part.dQdt[ietot] - DotProduct(part.v, part.dQdt, ndim) +
    0.5*DotProduct(part.v, part.v, ndim)*part.dQdt[irho];
  part.Utot += dUdt*timestep;
  for (int var=0; var<nvar; var++) {
    part.Qcons[var] += part.dQdt[var]*timestep;
  }

  return;
}



//=================================================================================================
//  MeshlessFV::UpdatePrimitiveVector
/// ...
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
/// ...
//=================================================================================================
template <int ndim>
void MeshlessFV<ndim>::UpdateArrayVariables(MeshlessFVParticle<ndim> &part)
{
  part.m = part.Qcons[irho];
  part.rho = part.m/part.volume; //part.Wprim[irho];
  for (int k=0; k<ndim; k++) part.v[k] = part.Qcons[k]/part.m;

  FLOAT ekin = 0.0;
  for (int k=0; k<ndim; k++) ekin += part.v[k]*part.v[k];
  part.u = (part.Qcons[ietot] - 0.5*part.m*ekin)/part.m;
  //part.u = part.U/part.m;
  part.press = (gamma_eos - 1.0)*part.rho*part.u;

  if (part.u < 0.0 || part.m < 0.0) {
    cout << "Mistake? : " << part.Qcons[ietot] << "    " << 0.5*part.m*ekin << "    " << part.m << "    " << part.u << endl;
    cout << "r : " << part.r[0] << "     v : " << part.v[0] << endl;
    cout << "Internal energy : " << part.u << "     " << part.Utot/part.m << endl;
  }

  assert(part.m > 0.0);
  assert(part.u > 0.0);

}



//=================================================================================================
//  MeshlessFV::InitialSmoothingLengthGuess
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

  debug2("[Sph::InitialSmoothingLengthGuess]");

  // Calculate bounding box containing all SPH particles
  this->ComputeBoundingBox(rmax,rmin,Nhydro);

  // Depending on the dimensionality, calculate the average smoothing
  // length assuming a uniform density distribution filling the bounding box.
  //-----------------------------------------------------------------------------------------------
  if (ndim == 1) {
    Ngather = (int) (2.0*kernp->kernrange*h_fac);
    volume = rmax[0] - rmin[0];
    h_guess = (volume*(FLOAT) Ngather)/(4.0*(FLOAT) Nhydro);
  }
  //-----------------------------------------------------------------------------------------------
  else if (ndim == 2) {
    Ngather = (int) (pi*pow(kernp->kernrange*h_fac,2));
    volume = (rmax[0] - rmin[0])*(rmax[1] - rmin[1]);
    h_guess = sqrtf((volume*(FLOAT) Ngather)/(4.0*(FLOAT) Nhydro));
  }
  //-----------------------------------------------------------------------------------------------
  else if (ndim == 3) {
    Ngather = (int) (4.0*pi*pow(kernp->kernrange*h_fac,3)/3.0);
    volume = (rmax[0] - rmin[0])*(rmax[1] - rmin[1])*(rmax[2] - rmin[2]);
    h_guess = powf((3.0*volume*(FLOAT) Ngather)/(32.0*pi*(FLOAT) Nhydro),onethird);
  }
  //-----------------------------------------------------------------------------------------------

  // Set all smoothing lengths equal to average value
  for (i=0; i<Nhydro; i++) {
    MeshlessFVParticle<ndim>& part = GetMeshlessFVParticlePointer(i);
    part.h         = h_guess;
    part.invh      = 1.0/h_guess;
    part.hrangesqd = kernfacsqd*kernp->kernrangesqd*part.h*part.h;
  }

  return;
}



template class MeshlessFV<1>;
template class MeshlessFV<2>;
template class MeshlessFV<3>;
