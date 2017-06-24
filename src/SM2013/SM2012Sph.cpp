//=================================================================================================
//  SM2012Sph.cpp
//  Contains all functions for calculating SPH quantities using the method
//  proposed by Saitoh & Makino (2012) in the conservative formalation
//  of Hopkins (2013).
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
#include "Debug.h"
#include "Precision.h"
#include "Exception.h"
#include "Sph.h"
#include "SmoothingKernel.h"
#include "Particle.h"
#include "Parameters.h"
#include "EOS.h"
#include "InlineFuncs.h"
using namespace std;



//=================================================================================================
//  SM2012Sph::SM2012Sph
/// Saitoh & Makino (2012) SPH object constructor.
//=================================================================================================
template <int ndim, template<int> class kernelclass>
SM2012Sph<ndim, kernelclass >::SM2012Sph(int hydro_forces_aux, int self_gravity_aux,
  FLOAT alpha_visc_aux, FLOAT beta_visc_aux, FLOAT h_fac_aux, FLOAT h_converge_aux,
  aviscenum avisc_aux, acondenum acond_aux, tdaviscenum tdavisc_aux,
  string gas_eos_aux, string KernelName, SimUnits &units, Parameters *params):
  Sph<ndim>(hydro_forces_aux, self_gravity_aux, alpha_visc_aux, beta_visc_aux,
            h_fac_aux, h_converge_aux, avisc_aux, acond_aux, tdavisc_aux,
            gas_eos_aux, KernelName, sizeof(SM2012SphParticle<ndim>), units, params),
  kern(kernelclass<ndim>(KernelName))
{
  this->kernp = &kern;
}



//=================================================================================================
//  SM2012Sph::~SM2012Sph
/// Saitoh & Makino (2012) SPH object destructor.
//=================================================================================================
template <int ndim, template<int> class kernelclass>
SM2012Sph<ndim, kernelclass >::~SM2012Sph()
{
}



//=================================================================================================
//  SM2012Sph::AllocateMemory
/// Allocate main SPH particle array.  Currently sets the maximum memory to
/// be 10 times the numbers of particles to allow space for ghost particles
/// and new particle creation.
//=================================================================================================
template <int ndim, template<int> class kernelclass>
void SM2012Sph<ndim, kernelclass>::AllocateMemory(int N)
{
  debug2("[SM2012Sph::AllocateMemory]");

  if (N > Nhydromax || !allocated) {
    if (allocated) DeallocateMemory();

    // Set conservative estimate for maximum number of particles, assuming
    // extra space required for periodic ghost particles
    if (Nhydromax < N)
      Nhydromax = pow(pow(N,invndim) + 8.0*kernp->kernrange,ndim);
    //Nhydromax = N;

    sphdata = new struct SM2012SphParticle<ndim>[Nhydromax];
    allocated = true;
    hydrodata_unsafe = sphdata;
    sphdata_unsafe = sphdata;
  }

  return;
}



//=================================================================================================
//  SM2012Sph::DeallocateMemory
/// Deallocate main array containing SPH particle data.
//=================================================================================================
template <int ndim, template<int> class kernelclass>
void SM2012Sph<ndim, kernelclass>::DeallocateMemory(void)
{
  debug2("[SM2012Sph::DeallocateMemory]");

  if (allocated) {
    delete[] sphdata;
  }
  allocated = false;

  return;
}


//=================================================================================================
//  SM2012Sph::ComputeH
/// Compute the value of the smoothing length of particle 'i' by iterating
/// the relation : h = h_fac*(m/rho)^(1/ndim).
/// Uses the previous value of h as a starting guess and then uses either
/// a Newton-Rhapson solver, or fixed-point iteration, to converge on the
/// correct value of h.  The maximum tolerance used for deciding whether the
/// iteration has converged is given by the 'h_converge' parameter.
//=================================================================================================
template <int ndim, template<int> class kernelclass>
int SM2012Sph<ndim, kernelclass >::ComputeH
 (SphParticle<ndim> &part,                                ///< [inout] Particle i data
  FLOAT hmax,                                             ///< [in] Maximum smoothing length
  const NeighbourList<DensityParticle> &ngbs,             ///< [in] Neighbour properties
  Nbody<ndim> *nbody)                                     ///< [in] Main N-body object
{
  int j;                            // Neighbour id
  int k;                            // Dimension counter
  int iteration = 0;                // h-rho iteration counter
  int iteration_max = 30;           // Max. no of iterations
  FLOAT dr[ndim];                   // Relative position vector
  FLOAT h_lower_bound = 0.0;        // Lower bound on h
  FLOAT h_upper_bound = hmax;       // Upper bound on h
  FLOAT invhsqd;                    // (1 / h)^2
  FLOAT w;                          // Kernel parameter squared, (r/h)^2

  SM2012SphParticle<ndim>& parti = static_cast<SM2012SphParticle<ndim>& > (part);

  FLOAT invh ;
  int Nneib = ngbs.size();

  // Main smoothing length iteration loop
  //===============================================================================================
  do {

    // Initialise all variables for this value of h
    iteration++;
    invh           = (FLOAT) 1.0/parti.h;
    parti.rho      = (FLOAT) 0.0;
    //parti.invomega = (FLOAT) 0.0;
    parti.q        = (FLOAT) 0.0;
    parti.hfactor  = pow(invh,ndim);
    invhsqd        = invh*invh;

    // Loop over all nearest neighbours in list to calculate
    // density.
    //---------------------------------------------------------------------------------------------
    for (j=0; j<Nneib; j++) {
      const DensityParticle &ngb = ngbs[j];
      for (k=0; k<ndim; k++) dr[k] = ngb.r[k] - parti.r[k];
      w          = kern.w0_s2(invhsqd*DotProduct(dr, dr, ndim));
      parti.rho += ngb.m*w;
      parti.q   += ngb.m*ngb.u*w;
    }
    //---------------------------------------------------------------------------------------------

    parti.rho *= parti.hfactor;
    parti.q *= parti.hfactor;

    FLOAT invrho = 0;
    if (parti.rho > (FLOAT) 0.0) invrho = (FLOAT) 1.0/parti.rho;

    // If h changes below some fixed tolerance, exit iteration loop
    if (parti.rho > (FLOAT) 0.0 && parti.h > h_lower_bound &&
        fabs(parti.h - h_fac*pow(parti.m*invrho,invndim)) < h_converge) break;

    // Use fixed-point iteration, i.e. h_new = h_fac*(m/rho_old)^(1/ndim),
    // for now.  If this does not converge in a reasonable number of
    // iterations (iteration_max), then assume something is wrong and switch
    // to a bisection method, which should be guaranteed to converge,
    // albeit much more slowly.  (N.B. will implement Newton-Raphson soon)
    //---------------------------------------------------------------------------------------------
    if (iteration < iteration_max)
      parti.h = h_fac*pow(parti.m*invrho,Sph<ndim>::invndim);

    else if (iteration == iteration_max)
      parti.h = (FLOAT) 0.5*(h_lower_bound + h_upper_bound);

    else if (iteration < 5*iteration_max) {
      if (parti.rho < small_number || parti.rho*pow(parti.h,ndim) > pow(h_fac,ndim)*parti.m)
        h_upper_bound = parti.h;
      else
        h_lower_bound = parti.h;
      parti.h = (FLOAT) 0.5*(h_lower_bound + h_upper_bound);
    }

    else {
      string message = "Problem with convergence of h-rho iteration";
      ExceptionHandler::getIstance().raise(message);
    }

    // If the smoothing length is too large for the neighbour list, exit
    // routine and flag neighbour list error in order to generate a larger
    // neighbour list (not properly implemented yet).
    if (parti.h > hmax) return 0;

  } while (parti.h > h_lower_bound && parti.h < h_upper_bound);
  //===============================================================================================


  // Normalise all SPH sums correctly
  parti.h         = max(h_fac*powf(parti.m/parti.rho,Sph<ndim>::invndim), h_lower_bound);
  invh            = (FLOAT) 1.0/parti.h;
  parti.hfactor   = pow(invh,ndim+1);
  parti.hrangesqd = kern.kernrangesqd*parti.h*parti.h;
  parti.div_v     = (FLOAT) 0.0;
  parti.dudt      = (FLOAT) 0.0;

  // Set important thermal variables here
  ComputeThermalProperties(parti);

  // Calculate the minimum neighbour potential
  // (used later to identify new sinks)
  if (create_sinks == 1) {
    parti.flags.set(potmin);
    for (j=0; j<Nneib; j++) {
      const DensityParticle &ngb = ngbs[j];
      for (k=0; k<ndim; k++) dr[k] = ngb.r[k] - parti.r[k];
      FLOAT drsqd = DotProduct(dr,dr,ndim);
      if (ngb.gpot > 1.000000001*parti.gpot &&
          drsqd*invhsqd < kern.kernrangesqd) parti.flags.unset(potmin);
    }
  }

  // If there are star particles, compute N-body chi correction term
  //-----------------------------------------------------------------------------------------------
  if (nbody->nbody_softening == 1) {
    for (j=0; j<nbody->Nstar; j++) {
      invhsqd = pow(2.0 / (parti.h + nbody->stardata[j].h),2);
      for (k=0; k<ndim; k++) dr[k] = nbody->stardata[j].r[k] - parti.r[k];
      //ssqd = DotProduct(dr,dr,ndim)*invhsqd;
    }
  }
  else {
    invhsqd = 4.0*invh*invh;
    for (j=0; j<nbody->Nstar; j++) {
      for (k=0; k<ndim; k++) dr[k] = nbody->stardata[j].r[k] - parti.r[k];
     // ssqd = DotProduct(dr,dr,ndim)*invhsqd;
    }
  }


  // If h is invalid (i.e. larger than maximum h), then return error code (0)
  if (parti.h <= hmax) return 1;
  else return -1;
}



//=================================================================================================
//  SM2012Sph::ComputeThermalProperties
/// Compute all thermal properties for grad-h SPH method for given particle.
//=================================================================================================
template <int ndim, template<int> class kernelclass>
void SM2012Sph<ndim, kernelclass>::ComputeThermalProperties
(SphParticle<ndim> &part_gen)        ///< [inout] Particle i data
{
  SM2012SphParticle<ndim>& part = static_cast<SM2012SphParticle<ndim> &> (part_gen);

  part.invq    = (FLOAT) 1.0/part.q;
  part.u       = eos->SpecificInternalEnergy(part);
  part.sound   = eos->SoundSpeed(part);
  part.pfactor = eos->Pressure(part)*part.invq/part.rho;

  return;
}



//=================================================================================================
//  SM2012Sph::ComputeSphHydroForces
/// Compute SPH neighbour force pairs for
/// (i) All neighbour interactions of particle i with id j > i,
/// (ii) Active neighbour interactions of particle j with id j > i
/// (iii) All inactive neighbour interactions of particle i with id j < i.
/// This ensures that all particle-particle pair interactions are only
/// computed once only for efficiency.
//=================================================================================================
template <int ndim, template<int> class kernelclass>
void SM2012Sph<ndim, kernelclass >::ComputeSphHydroForces
(const int i,                       ///< [in] id of particle
 const int Nneib,                   ///< [in] No. of neins in neibpart array
 const int *neiblist,               ///< [in] id of gather neibs in neibpart
 const FLOAT *drmag,                ///< [in] Distances of gather neighbours
 const FLOAT *invdrmag,             ///< [in] Inverse distances of gather neibs
 const FLOAT *dr,                   ///< [in] Position vector of gather neibs
 SphParticle<ndim> &part,           ///< [inout] Particle i data
 SphParticle<ndim> *neibpart_gen)   ///< [inout] Neighbour particle data
{
  int j;                            // Neighbour list id
  int jj;                           // Aux. neighbour counter
  int k;                            // Dimension counter
  FLOAT alpha_mean;                 // Mean articial viscosity alpha value
  FLOAT draux[ndim];                // Relative position vector
  FLOAT dv[ndim];                   // Relative velocity vector
  FLOAT dvdr;                       // Dot product of dv and dr
  FLOAT wkerni;                     // Value of w1 kernel function
  FLOAT wkernj;                     // Value of w1 kernel function
  FLOAT vsignal;                    // Signal velocity
  FLOAT paux;                       // Aux. pressure force variable
  FLOAT uaux;                       // Aux. internal energy variable
  FLOAT winvrho;                    // 0.5*(wkerni + wkernj)*invrhomean

  SM2012SphParticle<ndim>& parti = static_cast<SM2012SphParticle<ndim>& > (part);
  SM2012SphParticle<ndim>* neibpart = static_cast<SM2012SphParticle<ndim>* > (neibpart_gen);


  FLOAT invh_i   = 1/parti.h;
  FLOAT invrho_i = 1/parti.rho;

  // Loop over all potential neighbours in the list
  //-----------------------------------------------------------------------------------------------
  for (jj=0; jj<Nneib; jj++) {
    j = neiblist[jj];

    FLOAT invh_j   = 1/neibpart[j].h;
    FLOAT invrho_j = 1/neibpart[j].rho;

    wkerni = parti.hfactor*kern.w1(drmag[jj]*invh_i);
    wkernj = neibpart[j].hfactor*kern.w1(drmag[jj]*invh_j);

    for (k=0; k<ndim; k++) draux[k] = dr[jj*ndim + k];
    for (k=0; k<ndim; k++) dv[k] = neibpart[j].v[k] - parti.v[k];
    dvdr = DotProduct(dv,draux,ndim);

    // Add contribution to velocity divergence
    parti.div_v -= neibpart[j].m*dvdr*wkerni;
    neibpart[j].div_v -= parti.m*dvdr*wkernj;

    // Main SPH pressure force term
    paux = (FLOAT) 0.5*(eos->gamma - 1.0)*parti.u*neibpart[j].u*
    	        (parti.invq + neibpart[j].invq)*(wkerni + wkernj);


    // Add dissipation terms (for approaching particle pairs)
    //---------------------------------------------------------------------------------------------
    if (dvdr < (FLOAT) 0.0) {

      winvrho = (FLOAT) 0.25*(wkerni + wkernj)*(invrho_i + invrho_j);
      //winvrho = (FLOAT) (wkerni + wkernj)/(parti.rho + neibpart[j].rho);

      // Artificial viscosity term
      if (avisc == mon97) {
        vsignal = parti.sound + neibpart[j].sound - beta_visc*alpha_visc*dvdr;
        paux -= (FLOAT) alpha_visc*vsignal*dvdr*winvrho;
        uaux = (FLOAT) 0.5*alpha_visc*vsignal*dvdr*dvdr*winvrho;
        parti.dudt -= neibpart[j].m*uaux;
        neibpart[j].dudt -= parti.m*uaux;
      }
      else if (avisc == mon97mm97) {
        alpha_mean = 0.5*(parti.alpha + neibpart[j].alpha);
        vsignal = parti.sound + neibpart[j].sound - beta_visc*alpha_mean*dvdr;
        paux -= (FLOAT) alpha_mean*vsignal*dvdr*winvrho;
        uaux = (FLOAT) 0.5*alpha_mean*vsignal*dvdr*dvdr*winvrho;
        parti.dudt -= neibpart[j].m*uaux;
        neibpart[j].dudt -= parti.m*uaux;
      }

      // Artificial conductivity term
      if (acond == wadsley2008) {
        uaux = (FLOAT) 0.5*dvdr*(neibpart[j].u - parti.u)*
	      (invrho_i*wkerni + invrho_j*wkernj);
        parti.dudt += neibpart[j].m*uaux;
        neibpart[j].dudt -= parti.m*uaux;
      }
      else if (acond == price2008) {
    	vsignal = sqrt(fabs(eos->Pressure(parti) -
      		      eos->Pressure(neibpart[j]))*0.5*
      		 (invrho_i + invrho_j));
        parti.dudt += 0.5*neibpart[j].m*vsignal*
          (parti.u - neibpart[j].u)*winvrho;
        neibpart[j].dudt -= 0.5*parti.m*vsignal*
          (parti.u - neibpart[j].u)*winvrho;
      }

    }
    //---------------------------------------------------------------------------------------------


    uaux = 0.5*neibpart[j].m*neibpart[j].u*dvdr*(wkerni + wkernj);

    // Add total hydro contribution to acceleration for particle i
    for (k=0; k<ndim; k++) parti.a[k] += neibpart[j].m*draux[k]*paux;
    parti.dudt += 0.5*neibpart[j].m*neibpart[j].u*dvdr*(wkerni + wkernj)*
      parti.pfactor;

    // If neighbour is also active, add contribution to force here
    for (k=0; k<ndim; k++) neibpart[j].a[k] -= parti.m*draux[k]*paux;
    neibpart[j].dudt += 0.5*parti.m*parti.u*dvdr*(wkerni + wkernj)*
      neibpart[j].pfactor;

  }
  //-----------------------------------------------------------------------------------------------


  return;
}



template class SM2012Sph<1, M4Kernel>;
template class SM2012Sph<1, QuinticKernel>;
template class SM2012Sph<1, GaussianKernel>;
template class SM2012Sph<1, TabulatedKernel>;
template class SM2012Sph<2, M4Kernel>;
template class SM2012Sph<2, QuinticKernel>;
template class SM2012Sph<2, GaussianKernel>;
template class SM2012Sph<2, TabulatedKernel>;
template class SM2012Sph<3, M4Kernel>;
template class SM2012Sph<3, QuinticKernel>;
template class SM2012Sph<3, GaussianKernel>;
template class SM2012Sph<3, TabulatedKernel>;
