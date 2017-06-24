//=================================================================================================
//  MfvCommon.cpp
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
//  MfvCommon::MfvCommon
/// MfvCommon class constructor.  Calls main SPH class constructor and also
/// sets additional kernel-related quantities
//=================================================================================================
template <int ndim, template<int> class kernelclass, class SlopeLimiter>
MfvCommon<ndim, kernelclass,SlopeLimiter>::MfvCommon
  (int hydro_forces_aux, int self_gravity_aux, FLOAT _accel_mult, FLOAT _courant_mult,
   FLOAT _h_fac, FLOAT h_converge_aux, FLOAT gamma_aux, string gas_eos_aux, string KernelName,
   int size_part, SimUnits &units, Parameters *params):
  MeshlessFV<ndim>(hydro_forces_aux, self_gravity_aux, _accel_mult, _courant_mult, _h_fac,
                   h_converge_aux, gamma_aux, gas_eos_aux, KernelName, size_part, units, params),
  kern(kernelclass<ndim>(KernelName)),
  riemannExact(gamma_aux, params->intparams["zero_mass_flux"]),
  riemannHLLC(gamma_aux, params->intparams["zero_mass_flux"], gas_eos_aux == "isothermal"),
  viscosity(params->floatparams["shear_visc"], params->floatparams["bulk_visc"]),
  need_viscosity(params->floatparams["shear_visc"] || params->floatparams["bulk_visc"])
{
  this->kernp      = &kern;
  this->kernrange  = this->kernp->kernrange;

  // Local references to parameter variables for brevity
  map<string, string> &stringparams = params->stringparams;


  // Riemann solver object
  //-----------------------------------------------------------------------------------------------
  string riemann_solver = stringparams["riemann_solver"];
  if (riemann_solver == "exact") {
    RiemannSolverType = exact ;
  }
  else if (riemann_solver == "hllc") {
   RiemannSolverType = hllc ;
  }
  else {
    string message = "Unrecognised parameter : riemann_solver = " + riemann_solver;
    ExceptionHandler::getIstance().raise(message);
  }

}



//=================================================================================================
//  MfvCommon::~MfvCommon
/// MfvCommon class destructor
//=================================================================================================
template <int ndim, template<int> class kernelclass, class SlopeLimiter>
MfvCommon<ndim, kernelclass,SlopeLimiter>::~MfvCommon()
{ } ;



//=================================================================================================
//  MfvCommon::ComputeH
/// Compute the value of the smoothing length of particle 'i' by iterating the relation :
/// h = h_fac*(m/rho)^(1/ndim).
/// Uses the previous value of h as a starting guess and then uses either a Newton-Rhapson solver,
/// or fixed-point iteration, to converge on the correct value of h.  The maximum tolerance used
/// for deciding whether the iteration has converged is given by the 'h_converge' parameter.
//=================================================================================================
template <int ndim, template<int> class kernelclass, class SlopeLimiter>
int MfvCommon<ndim, kernelclass,SlopeLimiter>::ComputeH
 (MeshlessFVParticle<ndim> &part,        ///< [inout] Particle i data
  FLOAT hmax,                            ///< [in] Maximum allowable smoothing length
  const NeighbourList<DensNeib>& ngbs,   ///< [in] List of neighbours
  Nbody<ndim> *nbody)                    ///< [in] Main N-body object
{
  int j;                               // Neighbour id
  int iteration = 0;                   // h-rho iteration counter
  int iteration_max = 30;              // Max. no of iterations
  FLOAT h_lower_bound = (FLOAT) 0.0;   // Lower bound on h
  FLOAT h_upper_bound = hmax;          // Upper bound on h
  FLOAT invhsqd;                       // (1 / h)^2
  FLOAT ssqd;                          // Kernel parameter squared, (r/h)^2
  FLOAT dr[ndim];

  // If there are sink particles present, check if the particle is inside one.
  // If so, then adjust the iteration bounds and ensure they are valid (i.e. hmax is large enough)
  if (part.flags.check(inside_sink)) {
    h_lower_bound = hmin_sink;
    if (hmax < hmin_sink) return -1;
  }

  // Some basic sanity-checking in case of invalid input into routine
  assert(hmax > (FLOAT) 0.0);
  assert(!part.flags.is_dead());
  assert(part.m > (FLOAT) 0.0);

  FLOAT ndens, rho, volume, invomega, zeta, h, invh, hfactor ;

  int Nneib = ngbs.size();

  h = part.h ;
  // Main smoothing length iteration loop
  //===============================================================================================
  do {

    // Initialise all variables for this value of h
    iteration++;
    ndens    = 0;
    invomega = 0;
    zeta     = 0;
    invh     = 1/h;
    hfactor  = pow(invh,ndim);
    invhsqd  = invh*invh ;

    // Loop over all nearest neighbours in list to calculate density, omega and zeta.
    //---------------------------------------------------------------------------------------------
    for (j=0; j<Nneib; j++) {
      const DensNeib &ngb = ngbs[j];
      for (int k=0; k<ndim; k++) dr[k] = ngb.r[k] - part.r[k];
      ssqd = DotProduct(dr,dr,ndim)*invhsqd;
      ndens    += kern.w0_s2(ssqd);
      invomega += invh*kern.womega_s2(ssqd);
      zeta     += ngb.m*kern.wzeta_s2(ssqd);
    }
    //---------------------------------------------------------------------------------------------

    ndens    *= hfactor;
    invomega *= hfactor;
    zeta     *= invhsqd;
    volume    = 1/ndens;
    rho       = part.m*ndens;

    // If h changes below some fixed tolerance, exit iteration loop
    if (part.rho > (FLOAT) 0.0 && part.h > h_lower_bound &&
        fabs(part.h - h_fac*pow(volume,MeshlessFV<ndim>::invndim)) < h_converge) break;


    // Use fixed-point iteration, i.e. h_new = h_fac*(m/rho_old)^(1/ndim), for now.  If this does
    // not converge in a reasonable number of iterations (iteration_max), then assume something is
    // wrong and switch to a bisection method, which should be guaranteed to converge, albeit much
    // more slowly.  (N.B. will implement Newton-Raphson soon)
    //---------------------------------------------------------------------------------------------
    if (iteration < iteration_max) {
      h = h_fac*pow(volume,MeshlessFV<ndim>::invndim);
    }
    else if (iteration == iteration_max) {
      h = (FLOAT) 0.5*(h_lower_bound + h_upper_bound);
    }
    else if (iteration < 5*iteration_max) {
      if (ndens < small_number || ndens*pow(h,ndim) > pow(h_fac,ndim)) {
        h_upper_bound = h;
      }
      else {
        h_lower_bound = h;
      }
      h = (FLOAT) 0.5*(h_lower_bound + h_upper_bound);
    }
    else {
      cout << "H ITERATION : " << iteration << "    h : " << h
           << "   rho : " << rho << "   h_upper " << h_upper_bound << "    hmax :  " << hmax
           << "   h_lower : " << h_lower_bound << "    " << hfactor << "    m : " << part.m
           << "     " << part.m*hfactor*kern.w0(0.0) << "    " << Nneib << endl;
      string message = "Problem with convergence of h-rho iteration";
      ExceptionHandler::getIstance().raise(message);
    }

    // If the smoothing length is too large for the neighbour list, exit routine and flag neighbour
    // list error in order to generate a larger neighbour list (not properly implemented yet).
    if (h > hmax) return 0;

  } while (h > h_lower_bound && h < h_upper_bound);
  //===============================================================================================


  // Compute other terms once number density and smoothing length are known
  part.ndens     = ndens ;
  part.rho       = rho ;
  part.h         = h_fac*pow(volume, MeshlessFV<ndim>::invndim);
  part.hfactor   = pow(1/part.h, ndim+1);
  part.hrangesqd = kern.kernrangesqd*part.h*part.h;
  part.div_v     = (FLOAT) 0.0;
  part.invomega  = (FLOAT) 1.0 + (FLOAT) MeshlessFV<ndim>::invndim*part.h*part.invomega/part.ndens;
  part.invomega  = (FLOAT) 1.0/part.invomega;
  part.zeta      = -(FLOAT) MeshlessFV<ndim>::invndim*part.m*part.h*part.zeta*part.invomega/part.ndens;


  // Set important thermal variables here
  this->ComputeThermalProperties(part);
  this->UpdatePrimitiveVector(part);


  // If h is invalid (i.e. larger than maximum h), then return error code (0)
  if (part.h <= hmax) return 1;
  else return -1;
}



//=================================================================================================
//  MfvCommon::ComputeDerivatives
/// Compute Psi factors required for computing derivatives needed in Meshess FV equations.
//=================================================================================================
template <int ndim, template<int> class kernelclass, class SlopeLimiter>
void MfvCommon<ndim, kernelclass,SlopeLimiter>::ComputeGradients
(MeshlessFVParticle<ndim>& part,                                ///< [inout] Particle data
 NeighbourList<typename MeshlessFV<ndim>::GradNeib>& neibpart)  ///< [inout] Neighbour data
{
  int k;                                       // Dimension counter
  int var;                                     // Primitive variable counter
  FLOAT draux[ndim], dv[ndim];                 // Relative position / velocity vector
  FLOAT drsqd;                                 // Distance squared
  FLOAT dvdr ;                                 // Delta v , Delta r
  FLOAT E[ndim][ndim];                         // E-matrix for computing normalised B-matrix
  const FLOAT invhsqd = 1/(part.h*part.h);     // Local copy of 1/h^2
  FLOAT grad_tmp[nvar][ndim] ;                 // Workspace for computing gradient

  // Initialise/zero all variables to be updated in this routine
  part.vsig_max = (FLOAT) 0.0;
  for (var=0; var<nvar; var++) {
    for (k=0; k<ndim; k++) {
      grad_tmp[var][k] = (FLOAT) 0.0;
    }
  }
  for (k=0; k<ndim; k++) {
    for (int kk=0; kk<ndim; kk++) {
      E[k][kk] = (FLOAT) 0.0;
      part.B[k][kk] = (FLOAT) 0.0;
    }
  }
  if (create_sinks==1) part.flags.set(potmin);


  // Loop over all potential neighbours in the list
  //-----------------------------------------------------------------------------------------------
  int Nneib = neibpart.size();
  for (int j=0; j<Nneib; j++) {

    for (k=0; k<ndim; k++) draux[k] = neibpart[j].r[k] - part.r[k];
    for (k=0; k<ndim; k++) dv[k] = neibpart[j].v[k] - part.v[k];
    drsqd = DotProduct(draux, draux, ndim);
    dvdr = DotProduct(dv, draux, ndim);

    double w = part.hfactor*kern.w0_s2(drsqd*invhsqd)/part.ndens;

    for (k=0; k<ndim; k++) {
      for (int kk=0; kk<ndim; kk++) {
        E[k][kk] += draux[k]*draux[kk]*w ;
      }
    }

    // Compute the first part of gradient, we do the matrix mult later, when part.B is available.
    for (var=0; var <nvar; var++) {
      for (k=0; k<ndim; k++) {
        grad_tmp[var][k] += draux[k]*(neibpart[j].Wprim[var] - part.Wprim[var])*w ;
      }
    }

    // Calculate maximum signal velocity
    part.vsig_max = max((FLOAT) part.vsig_max, part.sound + neibpart[j].sound -
        min((FLOAT) 0.0, dvdr/(sqrtf(drsqd) + small_number)));
    part.levelneib = max(part.levelneib, neibpart[j].level) ;
    neibpart[j].levelneib = max(neibpart[j].levelneib, part.level);

    // Calculate the minimum neighbour potential (used later to identify new sinks)
    if (create_sinks == 1) {
        if (neibpart[j].gpot > (FLOAT) 1.000000001*part.gpot &&
            drsqd*invhsqd < kern.kernrangesqd) part.flags.unset(potmin);
    }

  }
  //-----------------------------------------------------------------------------------------------


  // Invert the matrix (depending on dimensionality)
  InvertMatrix<ndim>(E,part.B) ;

  // Check the accuracy of the integral gradients (using the square of the condition number)
  double modE(0), modB(0) ;
  for (int i=0; i<ndim; i++)
    for (int j=0; j<ndim; j++){
      modE +=      E[i][j]*     E[i][j];
      modB += part.B[i][j]*part.B[i][j];
    }
  double sqd_condition_number = modE*modB / (ndim*ndim) ;

  if (sqd_condition_number < 1e4) {
    part.flags.unset(bad_gradients);

    // Complete the calculation of the gradients:
    for (var=0; var<nvar; var++) {
      for (k=0; k<ndim; k++) {
        part.grad[var][k] = DotProduct(part.B[k], grad_tmp[var], ndim) ;
      }
    }
  }
  else {
    // Compute SPH (corrected) gradients instead
    part.flags.set(bad_gradients);


    cout << "Bad integral gradient, using SPH gradients:\n"
         << "\tiorig: " << part.iorig << ", ||E||, ||E^{-1}||: " << modE << " " << modB
         << ", N_cond^2: " << sqd_condition_number << "\n";

    for (var=0; var<nvar; var++)
      for (k=0; k<ndim; k++)
        part.grad[var][k] = 0;

    for (int j=0; j<Nneib; j++) {
      for (k=0; k<ndim; k++) draux[k] = neibpart[j].r[k] - part.r[k];
      drsqd = DotProduct(draux, draux, ndim);
      double dr = sqrt(drsqd);

      if (dr == 0) continue ;

      double w = (part.hfactor/part.ndens)*kern.w1(dr/part.h);

      for (var=0; var<nvar; var++)
        for (k=0; k<ndim; k++) {
          part.grad[var][k] -= (neibpart[j].Wprim[var] - part.Wprim[var])*(draux[k]/dr) * w ;
        }
    }
  }

  // Finally apply the slope limiter
  //-----------------------------------------------------------------------------------------------
  limiter.CellLimiter(part, neibpart) ;

  assert(part.vsig_max >= part.sound);

  return;
}



//=================================================================================================
//  MfvCommon::ComputeSmoothedGravForces
/// Compute SPH neighbour force pairs for
/// (i) All neighbour interactions of particle i with i.d. j > i,
/// (ii) Active neighbour interactions of particle j with i.d. j > i
/// (iii) All inactive neighbour interactions of particle i with i.d. j < i.
/// This ensures that particle-particle pair interactions are computed once only for efficiency.
//=================================================================================================
template <int ndim, template<int> class kernelclass, class SlopeLimiter>
void MfvCommon<ndim, kernelclass,SlopeLimiter>::ComputeSmoothedGravForces
(MeshlessFVParticle<ndim>& parti,                               ///< [inout] Particle data
 NeighbourList<typename MeshlessFV<ndim>::GravNeib>& neibpart)  ///< [inout] Neighbour data
{
  int k;                               // Dimension counter
  FLOAT dr[ndim];                      // Relative position vector
  FLOAT drmag;                         // Distance
  FLOAT invdrmag;                      // 1 / distance
  FLOAT gaux;                          // Aux. grav. potential variable
  FLOAT paux;                          // Aux. pressure force variable
  const FLOAT invh_i = 1.0/parti.h;
  const FLOAT invm = 1.0/parti.m;


  // Loop over all potential neighbours in the list
  //-----------------------------------------------------------------------------------------------
  int Nneib = neibpart.size();
  for (int j=0; j<Nneib; j++) {
    assert(!neibpart[j].flags.is_dead());

    const FLOAT mj = neibpart[j].m;
    const FLOAT invh_j = 1/neibpart[j].h;

    for (k=0; k<ndim; k++) dr[k] = neibpart[j].r[k] - parti.r[k];
    drmag = sqrt(DotProduct(dr, dr, ndim) + small_number);
    invdrmag = (FLOAT) 1.0/drmag;
    for (k=0; k<ndim; k++) dr[k] *= invdrmag;

    // Main SPH gravity terms
    //---------------------------------------------------------------------------------------------
    paux = (FLOAT) 0.5*(mj*invh_i*invh_i*kern.wgrav(drmag*invh_i) +
                        mj*invh_j*invh_j*kern.wgrav(drmag*invh_j) +
                        invm*parti.zeta*parti.hfactor*kern.w1(drmag*invh_i) +
                        invm*neibpart[j].zeta*neibpart[j].hfactor*kern.w1(drmag*invh_j));
    gaux = (FLOAT) 0.5*(invh_i*kern.wpot(drmag*invh_i) + invh_j*kern.wpot(drmag*invh_j));

    // Add total hydro contribution to acceleration for particle i
    for (k=0; k<ndim; k++) parti.atree[k] += dr[k]*paux;
    parti.gpot += neibpart[j].m*gaux;

  }
  //-----------------------------------------------------------------------------------------------

  return;
}



//=================================================================================================
//  MfvCommon::ComputeDirectGravForces
/// Compute the contribution to the total gravitational force of particle 'i'
/// due to 'Nneib' neighbouring particles in the list 'neiblist'.
//=================================================================================================
template <int ndim, template<int> class kernelclass, class SlopeLimiter>
void MfvCommon<ndim, kernelclass,SlopeLimiter>::ComputeDirectGravForces
(MeshlessFVParticle<ndim>& parti,                                 ///< [inout] Particle data
 NeighbourList<typename MeshlessFV<ndim>::DirectNeib>& neibdata)  ///< [inout] Neighbour data
{
  int k;                               // Dimension counter
  FLOAT dr[ndim];                      // Relative position vector
  FLOAT drsqd;                         // Distance squared
  FLOAT invdrmag;                      // 1 / distance
  FLOAT invdr3;                        // 1 / dist^3

  // Loop over all neighbouring particles in list
  //-----------------------------------------------------------------------------------------------
  int Nneib = neibdata.size();
  for (int j=0; j<Nneib; j++) {
    assert(!neibdata[j].flags.is_dead());

    for (k=0; k<ndim; k++) dr[k] = neibdata[j].r[k] - parti.r[k];
    drsqd    = DotProduct(dr,dr,ndim) + small_number;
    invdrmag = (FLOAT) 1.0/sqrt(drsqd);
    invdr3   = invdrmag*invdrmag*invdrmag;

    // Add contribution to current particle
    for (k=0; k<ndim; k++) parti.atree[k] += neibdata[j].m*dr[k]*invdr3;
    parti.gpot += neibdata[j].m*invdrmag;

    assert(drsqd >= parti.hrangesqd && drsqd >= neibdata[j].hrangesqd);

  }
  //-----------------------------------------------------------------------------------------------

  return;
}



//=================================================================================================
//  MfvCommon::ComputeStarGravForces
/// Computes contribution of gravitational force and potential due to stars.
//=================================================================================================
template <int ndim, template<int> class kernelclass, class SlopeLimiter>
void MfvCommon<ndim, kernelclass,SlopeLimiter>::ComputeStarGravForces
 (const int N,                         ///< [in] No. of stars
  NbodyParticle<ndim> **nbodydata,     ///< [in] Array of star pointers
  MeshlessFVParticle<ndim> &part)      ///< [inout] SPH particle reference
{
  int j;                               // Star counter
  int k;                               // Dimension counter
  FLOAT dr[ndim];                      // Relative position vector
  FLOAT drmag;                         // Distance
  FLOAT drsqd;                         // Distance squared
  //FLOAT drdt;                          // Rate of change of relative distance
  //FLOAT dv[ndim];                      // Relative velocity vector
  FLOAT invdrmag;                      // 1 / drmag
  FLOAT invhmean;                      // 1 / hmean
  FLOAT ms;                            // Star mass
  FLOAT paux;                          // Aux. force variable
  MeshlessFVParticle<ndim>& parti = static_cast<MeshlessFVParticle<ndim>& > (part);

  // Loop over all stars and add each contribution
  //-----------------------------------------------------------------------------------------------
  for (j=0; j<N; j++) {

    ms = nbodydata[j]->m;

    for (k=0; k<ndim; k++) dr[k] = nbodydata[j]->r[k] - parti.r[k];
    drsqd    = DotProduct(dr,dr,ndim) + small_number;
    drmag    = sqrt(drsqd);
    invdrmag = (FLOAT) 1.0/drmag;
    invhmean = (FLOAT) 2.0/(parti.h + nbodydata[j]->h);
    paux     = ms*invhmean*invhmean*kern.wgrav(drmag*invhmean)*invdrmag;

    // Add total hydro contribution to acceleration for particle i
    for (k=0; k<ndim; k++) parti.atree[k] += paux*dr[k];
    parti.gpot += ms*invhmean*kern.wpot(drmag*invhmean);

    assert(drmag > (FLOAT) 0.0);
    assert(drmag*invhmean > (FLOAT) 0.0);

  }
  //-----------------------------------------------------------------------------------------------

  return;
}



template class MfvCommon<1, M4Kernel, NullLimiter<1> >;
template class MfvCommon<2, M4Kernel, NullLimiter<2> >;
template class MfvCommon<3, M4Kernel, NullLimiter<3> >;
template class MfvCommon<1, QuinticKernel, NullLimiter<1> >;
template class MfvCommon<2, QuinticKernel, NullLimiter<2> >;
template class MfvCommon<3, QuinticKernel, NullLimiter<3> >;
template class MfvCommon<1, TabulatedKernel, NullLimiter<1> >;
template class MfvCommon<2, TabulatedKernel, NullLimiter<2> >;
template class MfvCommon<3, TabulatedKernel, NullLimiter<3> >;


template class MfvCommon<1, M4Kernel, ZeroSlopeLimiter<1> >;
template class MfvCommon<2, M4Kernel, ZeroSlopeLimiter<2> >;
template class MfvCommon<3, M4Kernel, ZeroSlopeLimiter<3> >;
template class MfvCommon<1, QuinticKernel, ZeroSlopeLimiter<1> >;
template class MfvCommon<2, QuinticKernel, ZeroSlopeLimiter<2> >;
template class MfvCommon<3, QuinticKernel, ZeroSlopeLimiter<3> >;
template class MfvCommon<1, TabulatedKernel, ZeroSlopeLimiter<1> >;
template class MfvCommon<2, TabulatedKernel, ZeroSlopeLimiter<2> >;
template class MfvCommon<3, TabulatedKernel, ZeroSlopeLimiter<3> >;

template class MfvCommon<1, M4Kernel, TVDScalarLimiter<1> >;
template class MfvCommon<2, M4Kernel, TVDScalarLimiter<2> >;
template class MfvCommon<3, M4Kernel, TVDScalarLimiter<3> >;
template class MfvCommon<1, QuinticKernel, TVDScalarLimiter<1> >;
template class MfvCommon<2, QuinticKernel, TVDScalarLimiter<2> >;
template class MfvCommon<3, QuinticKernel, TVDScalarLimiter<3> >;
template class MfvCommon<1, TabulatedKernel,TVDScalarLimiter<1> >;
template class MfvCommon<2, TabulatedKernel, TVDScalarLimiter<2> >;
template class MfvCommon<3, TabulatedKernel, TVDScalarLimiter<3> >;

template class MfvCommon<1, M4Kernel, ScalarLimiter<1> >;
template class MfvCommon<2, M4Kernel, ScalarLimiter<2> >;
template class MfvCommon<3, M4Kernel, ScalarLimiter<3> >;
template class MfvCommon<1, QuinticKernel, ScalarLimiter<1> >;
template class MfvCommon<2, QuinticKernel, ScalarLimiter<2> >;
template class MfvCommon<3, QuinticKernel, ScalarLimiter<3> >;
template class MfvCommon<1, TabulatedKernel, ScalarLimiter<1> >;
template class MfvCommon<2, TabulatedKernel, ScalarLimiter<2> >;
template class MfvCommon<3, TabulatedKernel, ScalarLimiter<3> >;

template class MfvCommon<1, M4Kernel, Springel2009Limiter<1> >;
template class MfvCommon<2, M4Kernel, Springel2009Limiter<2> >;
template class MfvCommon<3, M4Kernel, Springel2009Limiter<3> >;
template class MfvCommon<1, QuinticKernel, Springel2009Limiter<1> >;
template class MfvCommon<2, QuinticKernel, Springel2009Limiter<2> >;
template class MfvCommon<3, QuinticKernel, Springel2009Limiter<3> >;
template class MfvCommon<1, TabulatedKernel, Springel2009Limiter<1> >;
template class MfvCommon<2, TabulatedKernel, Springel2009Limiter<2> >;
template class MfvCommon<3, TabulatedKernel, Springel2009Limiter<3> >;

template class MfvCommon<1, M4Kernel, GizmoLimiter<1> >;
template class MfvCommon<2, M4Kernel, GizmoLimiter<2> >;
template class MfvCommon<3, M4Kernel, GizmoLimiter<3> >;
template class MfvCommon<1, QuinticKernel, GizmoLimiter<1> >;
template class MfvCommon<2, QuinticKernel, GizmoLimiter<2> >;
template class MfvCommon<3, QuinticKernel, GizmoLimiter<3> >;
template class MfvCommon<1, TabulatedKernel, GizmoLimiter<1> >;
template class MfvCommon<2, TabulatedKernel, GizmoLimiter<2> >;
template class MfvCommon<3, TabulatedKernel, GizmoLimiter<3> >;
