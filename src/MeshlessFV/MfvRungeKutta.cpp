//=================================================================================================
//  MfvRungeKutta.cpp
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
//  MfvRungeKutta::MfvRungeKutta
/// MfvRungeKutta class constructor.  Calls main SPH class constructor and also
/// sets additional kernel-related quantities
//=================================================================================================
template <int ndim, template<int> class kernelclass, class SlopeLimiter>
MfvRungeKutta<ndim, kernelclass,SlopeLimiter>::MfvRungeKutta
 (int _hydro_forces, int _self_gravity, FLOAT _accel_mult, FLOAT _courant_mult,
  FLOAT _h_fac, FLOAT _h_converge, FLOAT _gamma, string _gas_eos, string KernelName,
  int size_part, SimUnits &units, Parameters *params):
  MfvCommon<ndim,kernelclass,SlopeLimiter>(_hydro_forces, _self_gravity, _accel_mult, _courant_mult, _h_fac,
                              _h_converge, _gamma, _gas_eos, KernelName, size_part, units, params)
{
}



//=================================================================================================
//  MfvRungeKutta::ComputeGodunovFlux
/// ...
//=================================================================================================
template <int ndim, template<int> class kernelclass, class SlopeLimiter>
void MfvRungeKutta<ndim, kernelclass,SlopeLimiter>::ComputeGodunovFlux
 (const int i,                         ///< [in] id of particle
  const int Nneib,                     ///< [in] No. of neins in neibpart array
  const int *neiblist,                 ///< [in] id of gather neibs in neibpart
  const FLOAT timestep,                ///< [in] Minimum timestep size
  MeshlessFVParticle<ndim> &part,      ///< [inout] Particle i data
  MeshlessFVParticle<ndim> *neibpart)  ///< [inout] Neighbour particle data
{
  int j;                               // Neighbour list id
  int jj;                              // Aux. neighbour counter
  int k;                               // Dimension counter
  int var;                             // Particle state vector variable counter
  FLOAT Aij[ndim];                     // Pseudo 'Area' vector
  FLOAT Aunit[ndim];                   // ..
  FLOAT draux[ndim];                   // Position vector of part relative to neighbour
  FLOAT drsqd;                         // Distance squared
  FLOAT psitildai[ndim];               // Normalised gradient psi value for particle i
  FLOAT psitildaj[ndim];               // Normalised gradient psi value for neighbour j
  FLOAT rface[ndim];                   // Position of working face (to compute Godunov fluxes)
  FLOAT vface[ndim];                   // Velocity of working face (to compute Godunov fluxes)
  FLOAT flux[nvar][ndim];              // Flux tensor
  FLOAT Wi[nvar];                      // Primitive vector for LHS of Riemann problem
  FLOAT Wj[nvar];                      // Primitive vector for RHS of Riemann problem
  FLOAT Wdot[nvar];                    // Time derivative of primitive vector
  FLOAT gradW[nvar][ndim];             // Gradient of primitive vector
  FLOAT dW[nvar];                      // Change in primitive quantities
  const FLOAT dt = timestep*(FLOAT) part.nstep;    // Timestep of given particle
  const FLOAT invh_i   = 1.0/part.h;
  const FLOAT volume_i = 1.0/part.ndens;

  // Loop over all potential neighbours in the list
  //-----------------------------------------------------------------------------------------------
  for (jj=0; jj<Nneib; jj++) {
    j = neiblist[jj];

    const FLOAT invh_j   = (FLOAT) 1.0/neibpart[j].h;
    const FLOAT volume_j = (FLOAT) 1/neibpart[j].ndens;

    for (k=0; k<ndim; k++) draux[k] = neibpart[j].r[k] - part.r[k];
    drsqd = DotProduct(draux, draux, ndim);


    // Calculate psitilda values
    for (k=0; k<ndim; k++) {
      psitildai[k] = (FLOAT) 0.0;
      psitildaj[k] = (FLOAT) 0.0;
      for (int kk=0; kk<ndim; kk++) {
        psitildai[k] += neibpart[j].B[k][kk]*draux[kk]*neibpart[j].hfactor*
          kern.w0_s2(drsqd*invh_j*invh_j)*volume_j;
        psitildaj[k] -= part.B[k][kk]*draux[kk]*part.hfactor*
          kern.w0_s2(drsqd*invh_i*invh_i)*volume_i;
      }
      Aij[k] = volume_i*psitildaj[k] - volume_j*psitildai[k];
    }

    FLOAT Amag = sqrt(DotProduct(Aij, Aij, ndim) + small_number);
    for (k=0; k<ndim; k++) Aunit[k] = Aij[k] / Amag;

    // Calculate position and velocity of the face
    if (staticParticles) {
      for (k=0; k<ndim; k++) rface[k] = (FLOAT) 0.5*(part.r[k] + neibpart[j].r[k]);
      for (k=0; k<ndim; k++) vface[k] = (FLOAT) 0.0;
    }
    else {
      for (k=0; k<ndim; k++) rface[k] = (FLOAT) 0.5*(part.r[k] + neibpart[j].r[k]);
      for (k=0; k<ndim; k++) vface[k] = (FLOAT) 0.5*(part.v[k] + neibpart[j].v[k]);
    }

    // Compute slope-limited values for LHS
    for (k=0; k<ndim; k++) draux[k] = rface[k] - part.r[k];
    limiter.ComputeLimitedSlopes(part, neibpart[j], draux, gradW, dW);
    for (var=0; var<nvar; var++) Wi[var] = part.Wprim[var] + dW[var];
    for (k=0; k<ndim; k++) Wi[k] -= vface[k];


    // Compute slope-limited values for RHS
    for (k=0; k<ndim; k++) draux[k] = rface[k] - neibpart[j].r[k];
    limiter.ComputeLimitedSlopes(neibpart[j], part, draux, gradW, dW);
    for (var=0; var<nvar; var++) Wj[var] = neibpart[j].Wprim[var] + dW[var];
    for (k=0; k<ndim; k++) Wj[k] -= vface[k];


    // Pressure and density floors incase of very strong gradients/shocks or pathological cases
    Wi[irho] = max(Wi[irho], small_number);
    Wj[irho] = max(Wj[irho], small_number);
    Wi[ipress] = max(Wi[ipress], small_number);
    Wj[ipress] = max(Wj[ipress], small_number);

    assert(isnormal(Wi[irho]));
    assert(isnormal(Wi[ipress]));
    assert(isnormal(Wj[irho]));
    assert(isnormal(Wj[ipress]));

    // Calculate Godunov flux using the selected Riemann solver
    if (RiemannSolverType == exact) {
      riemannExact.ComputeFluxes(Wj, Wi, Aunit, vface, flux);
    }
    else {
      riemannHLLC.ComputeFluxes(Wj, Wi, Aunit, vface, flux);
    }

    // Finally calculate flux terms for all quantities based on Lanson & Vila gradient operators
    for (var=0; var<nvar; var++) {
      const FLOAT f = DotProduct(flux[var], Aij, ndim);
      part.dQ[var] += f*dt;
      part.dQdt[var] += f;
      neibpart[j].dQ[var] -= f*dt;
      neibpart[j].dQdt[var] -= f;
    }

    // Compute mass-loss moments for gravitational correction terms
    for (k=0; k<ndim; k++) {
      part.rdmdt[k] -= (part.r[k] - neibpart[j].r[k])*DotProduct(flux[irho], Aij, ndim);
      neibpart[j].rdmdt[k] += (part.r[k] - neibpart[j].r[k])*DotProduct(flux[irho], Aij, ndim);
    }
  }
  //-----------------------------------------------------------------------------------------------


  return;
}



//=================================================================================================
//  MeshlessFV<ndim>::IntegrateParticles
/// Calculate or reset all quantities for all particles that reach the end of their timesteps.
//=================================================================================================
template <int ndim, template<int> class kernelclass, class SlopeLimiter>
void MfvRungeKutta<ndim, kernelclass,SlopeLimiter>::IntegrateParticles
 (const int n,                         ///< [in] Integer time in block time struct
  const FLOAT t,                       ///< [in] Current simulation time
  const FLOAT timestep,                ///< [in] Base timestep value
  const DomainBox<ndim> &simbox)       ///< [in] Simulation box
{
  debug2("[MfvRungeKutta:::IntegrateParticles]");
  CodeTiming::BlockTimer timer = timing->StartNewTimer("MFVRK_INTEGRATE_PARTICLES");

  MeshlessFVParticle<ndim>* partdata = GetMeshlessFVParticleArray() ;

  // Integrate all conserved variables to end of timestep
  //-----------------------------------------------------------------------------------------------
  for (int i=0; i<Nhydro; i++) {
    MeshlessFVParticle<ndim> &part = partdata[i];
    if (part.flags.is_dead()) continue;
    const int dn = n - part.nlast;
    const FLOAT dt = timestep*(FLOAT) dn;
    FLOAT Qcons[nvar];

    for (int k=0; k<nvar; k++) Qcons[k] = part.Qcons0[k] + part.dQdt[k]*dt;
    for (int k=0; k<ndim; k++) Qcons[k] += part.Qcons0[irho]*part.a0[k]*dt;

    part.flags.unset_flag(active);
    if (dn == part.nstep) {
      part.flags.set_flag(active);
      // TODO: This should be done at the step mid-point
      for (int k=0; k<ndim; k++) part.rdmdt0[k] = part.rdmdt[k];
      for (int k=0; k<ndim; k++) part.rdmdt[k]  = 0;
    }


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

          // Check if wall or mirror boundary
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
template <int ndim, template<int> class kernelclass, class SlopeLimiter>
void MfvRungeKutta<ndim, kernelclass,SlopeLimiter>::EndTimestep
 (const int n,                         ///< [in] Integer time in block time struct
  const FLOAT t,                       ///< [in] Current simulation time
  const FLOAT timestep)                ///< [in] Base timestep value
{
  debug2("[MfvRungeKutta::EndTimestep]");
  CodeTiming::BlockTimer timer = timing->StartNewTimer("MFVRK_END_TIMESTEP");

  MeshlessFVParticle<ndim>* partdata = GetMeshlessFVParticleArray() ;

  //-----------------------------------------------------------------------------------------------
#pragma omp parallel for default(none) shared(partdata)
  for (int i=0; i<Nhydro; i++) {
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
        // Factor of 0.5 here is because we are averaging the flux from both the start and end of
        // the step
        Qcons[var] = part.Qcons0[var] + 0.5*part.dQ[var];
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
      for (k=0; k<ndim; k++) part.a[k]      = 0;
      for (k=0; k<nvar; k++) part.Qcons0[k] = Qcons[k];
      for (k=0; k<ndim; k++) part.rdmdt0[k] = 0;
      for (k=0; k<ndim; k++) part.rdmdt[k]  = 0;
      part.gpot = 0;
    }
    //---------------------------------------------------------------------------------------------
    else {
      part.flags.unset_flag(active);
    }
    //---------------------------------------------------------------------------------------------

  }
  //-----------------------------------------------------------------------------------------------

  return;
}





template class MfvRungeKutta<1, M4Kernel, NullLimiter<1,MeshlessFVParticle> >;
template class MfvRungeKutta<2, M4Kernel, NullLimiter<2,MeshlessFVParticle> >;
template class MfvRungeKutta<3, M4Kernel, NullLimiter<3,MeshlessFVParticle> >;
template class MfvRungeKutta<1, QuinticKernel, NullLimiter<1,MeshlessFVParticle> >;
template class MfvRungeKutta<2, QuinticKernel, NullLimiter<2,MeshlessFVParticle> >;
template class MfvRungeKutta<3, QuinticKernel, NullLimiter<3,MeshlessFVParticle> >;
template class MfvRungeKutta<1, TabulatedKernel, NullLimiter<1,MeshlessFVParticle> >;
template class MfvRungeKutta<2, TabulatedKernel, NullLimiter<2,MeshlessFVParticle> >;
template class MfvRungeKutta<3, TabulatedKernel, NullLimiter<3,MeshlessFVParticle> >;


template class MfvRungeKutta<1, M4Kernel, ZeroSlopeLimiter<1,MeshlessFVParticle> >;
template class MfvRungeKutta<2, M4Kernel, ZeroSlopeLimiter<2,MeshlessFVParticle> >;
template class MfvRungeKutta<3, M4Kernel, ZeroSlopeLimiter<3,MeshlessFVParticle> >;
template class MfvRungeKutta<1, QuinticKernel, ZeroSlopeLimiter<1,MeshlessFVParticle> >;
template class MfvRungeKutta<2, QuinticKernel, ZeroSlopeLimiter<2,MeshlessFVParticle> >;
template class MfvRungeKutta<3, QuinticKernel, ZeroSlopeLimiter<3,MeshlessFVParticle> >;
template class MfvRungeKutta<1, TabulatedKernel, ZeroSlopeLimiter<1,MeshlessFVParticle> >;
template class MfvRungeKutta<2, TabulatedKernel, ZeroSlopeLimiter<2,MeshlessFVParticle> >;
template class MfvRungeKutta<3, TabulatedKernel, ZeroSlopeLimiter<3,MeshlessFVParticle> >;

template class MfvRungeKutta<1, M4Kernel, TVDScalarLimiter<1,MeshlessFVParticle> >;
template class MfvRungeKutta<2, M4Kernel, TVDScalarLimiter<2,MeshlessFVParticle> >;
template class MfvRungeKutta<3, M4Kernel, TVDScalarLimiter<3,MeshlessFVParticle> >;
template class MfvRungeKutta<1, QuinticKernel, TVDScalarLimiter<1,MeshlessFVParticle> >;
template class MfvRungeKutta<2, QuinticKernel, TVDScalarLimiter<2,MeshlessFVParticle> >;
template class MfvRungeKutta<3, QuinticKernel, TVDScalarLimiter<3,MeshlessFVParticle> >;
template class MfvRungeKutta<1, TabulatedKernel,TVDScalarLimiter<1,MeshlessFVParticle> >;
template class MfvRungeKutta<2, TabulatedKernel, TVDScalarLimiter<2,MeshlessFVParticle> >;
template class MfvRungeKutta<3, TabulatedKernel, TVDScalarLimiter<3,MeshlessFVParticle> >;

template class MfvRungeKutta<1, M4Kernel, ScalarLimiter<1,MeshlessFVParticle> >;
template class MfvRungeKutta<2, M4Kernel, ScalarLimiter<2,MeshlessFVParticle> >;
template class MfvRungeKutta<3, M4Kernel, ScalarLimiter<3,MeshlessFVParticle> >;
template class MfvRungeKutta<1, QuinticKernel, ScalarLimiter<1,MeshlessFVParticle> >;
template class MfvRungeKutta<2, QuinticKernel, ScalarLimiter<2,MeshlessFVParticle> >;
template class MfvRungeKutta<3, QuinticKernel, ScalarLimiter<3,MeshlessFVParticle> >;
template class MfvRungeKutta<1, TabulatedKernel, ScalarLimiter<1,MeshlessFVParticle> >;
template class MfvRungeKutta<2, TabulatedKernel, ScalarLimiter<2,MeshlessFVParticle> >;
template class MfvRungeKutta<3, TabulatedKernel, ScalarLimiter<3,MeshlessFVParticle> >;

template class MfvRungeKutta<1, M4Kernel, Springel2009Limiter<1,MeshlessFVParticle> >;
template class MfvRungeKutta<2, M4Kernel, Springel2009Limiter<2,MeshlessFVParticle> >;
template class MfvRungeKutta<3, M4Kernel, Springel2009Limiter<3,MeshlessFVParticle> >;
template class MfvRungeKutta<1, QuinticKernel, Springel2009Limiter<1,MeshlessFVParticle> >;
template class MfvRungeKutta<2, QuinticKernel, Springel2009Limiter<2,MeshlessFVParticle> >;
template class MfvRungeKutta<3, QuinticKernel, Springel2009Limiter<3,MeshlessFVParticle> >;
template class MfvRungeKutta<1, TabulatedKernel, Springel2009Limiter<1,MeshlessFVParticle> >;
template class MfvRungeKutta<2, TabulatedKernel, Springel2009Limiter<2,MeshlessFVParticle> >;
template class MfvRungeKutta<3, TabulatedKernel, Springel2009Limiter<3,MeshlessFVParticle> >;

template class MfvRungeKutta<1, M4Kernel, GizmoLimiter<1,MeshlessFVParticle> >;
template class MfvRungeKutta<2, M4Kernel, GizmoLimiter<2,MeshlessFVParticle> >;
template class MfvRungeKutta<3, M4Kernel, GizmoLimiter<3,MeshlessFVParticle> >;
template class MfvRungeKutta<1, QuinticKernel, GizmoLimiter<1,MeshlessFVParticle> >;
template class MfvRungeKutta<2, QuinticKernel, GizmoLimiter<2,MeshlessFVParticle> >;
template class MfvRungeKutta<3, QuinticKernel, GizmoLimiter<3,MeshlessFVParticle> >;
template class MfvRungeKutta<1, TabulatedKernel, GizmoLimiter<1,MeshlessFVParticle> >;
template class MfvRungeKutta<2, TabulatedKernel, GizmoLimiter<2,MeshlessFVParticle> >;
template class MfvRungeKutta<3, TabulatedKernel, GizmoLimiter<3,MeshlessFVParticle> >;

