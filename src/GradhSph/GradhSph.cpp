//=================================================================================================
//  GradhSph.cpp
//  Contains all functions for calculating conservative 'grad-h' SPH quantities
//  (See Springel & Hernquist (2002) and Price & Monaghan (2007)).
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
#include "Sph.h"
#include "Particle.h"
#include "Parameters.h"
#include "SmoothingKernel.h"
#include "EOS.h"
#include "Debug.h"
#include "Exception.h"
#include "InlineFuncs.h"
using namespace std;



//=================================================================================================
//  GradhSph::GradhSph
/// GradhSph class constructor.  Calls main SPH class constructor and also
/// sets additional kernel-related quantities
//=================================================================================================
template <int ndim, template<int> class kernelclass>
GradhSph<ndim, kernelclass>::GradhSph(int hydro_forces_aux, int self_gravity_aux,
  FLOAT alpha_visc_aux, FLOAT beta_visc_aux, FLOAT h_fac_aux, FLOAT h_converge_aux,
  aviscenum avisc_aux, acondenum acond_aux, tdaviscenum tdavisc_aux,
  string gas_eos_aux, string KernelName, SimUnits &units, Parameters *params):
  Sph<ndim>(hydro_forces_aux, self_gravity_aux, alpha_visc_aux, beta_visc_aux,
            h_fac_aux, h_converge_aux, avisc_aux, acond_aux, tdavisc_aux,
            gas_eos_aux, KernelName, sizeof(GradhSphParticle<ndim>), units, params),
  kern(kernelclass<ndim>(KernelName))
{
  this->kernp      = &kern;
  this->kernfac    = (FLOAT) 1.0;
  this->kernfacsqd = (FLOAT) 1.0;
  this->kernrange  = this->kernp->kernrange;
}



//=================================================================================================
//  GradhSph::~GradhSph
/// GradhSph class destructor
//=================================================================================================
template <int ndim, template<int> class kernelclass>
GradhSph<ndim, kernelclass>::~GradhSph()
{
  DeallocateMemory();
}



//=================================================================================================
//  GradhSph::AllocateMemory
/// Allocate main SPH particle array.  Estimates the maximum number of boundary ghost particles
/// assuming a roughly uniform depth of ghosts at each boundary.
//=================================================================================================
template <int ndim, template<int> class kernelclass>
void GradhSph<ndim, kernelclass>::AllocateMemory(int N)
{
  debug2("[GradhSph::AllocateMemory]");

  if (N > Nhydromax || !allocated) {
    if (allocated) DeallocateMemory();

    // Set conservative estimate for maximum number of particles, assuming
    // extra space required for periodic ghost particles
    if (Nhydromax < N) {
      Nhydromax = 2*(int) powf(powf((FLOAT) N,invndim) + (FLOAT) 16.0*kernp->kernrange,ndim);
    }

    iorder           = new int[Nhydromax];
    sphdata          = new struct GradhSphParticle<ndim>[Nhydromax];
    allocated        = true;
    hydrodata_unsafe = sphdata;
    sphdata_unsafe   = sphdata;
  }

  assert(Nhydromax > Nhydro);
  assert(iorder);
  assert(sphdata);

  return;
}



//=================================================================================================
//  GradhSph::DeallocateMemory
/// Deallocate main array containing SPH particle data.
//=================================================================================================
template <int ndim, template<int> class kernelclass>
void GradhSph<ndim, kernelclass>::DeallocateMemory(void)
{
  debug2("[GradhSph::DeallocateMemory]");

  if (allocated) {
    delete[] sphdata;
    delete[] iorder;
  }
  allocated = false;

  return;
}



//=================================================================================================
//  GradhSph::DeleteDeadParticles
/// Delete 'dead' (e.g. accreted) SPH particles from the main arrays.
//=================================================================================================
template <int ndim, template<int> class kernelclass>
void GradhSph<ndim, kernelclass>::DeleteDeadParticles(void)
{
  int i;                               // Particle counter
  int itype;                           // Current particle type
  int Ndead = 0;                       // No. of 'dead' particles
  int ilast = Nhydro;                  // Aux. counter of last free slot

  debug2("[GradhSph::DeleteDeadParticles]");


  // Determine new order of particles in arrays.
  // First all live particles and then all dead particles.
  for (i=0; i<Nhydro; i++) {
    itype = sphdata[i].itype;
    while (itype == dead) {
      Ndead++;
      ilast--;
      if (i < ilast) {
        sphdata[i] = sphdata[ilast];
        sphdata[ilast].itype = dead;
        sphdata[ilast].m = (FLOAT) 0.0;
      }
      else break;
      itype = sphdata[i].itype;
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
    assert(sphdata[i].itype != dead);
  }

  return;
}



//=================================================================================================
//  GradhSph::ReorderParticles
/// Delete selected SPH particles from the main arrays.
//=================================================================================================
template <int ndim, template<int> class kernelclass>
void GradhSph<ndim, kernelclass>::ReorderParticles(void)
{
  int i;                               // Particle counter
  GradhSphParticle<ndim> *sphdataaux;  // Aux. SPH particle array

  sphdataaux = new GradhSphParticle<ndim>[Nhydro];

  for (i=0; i<Nhydro; i++) sphdataaux[i] = sphdata[i];
  for (i=0; i<Nhydro; i++) sphdata[i] = sphdataaux[iorder[i]];

  delete[] sphdataaux;

  return;
}



//=================================================================================================
//  GradhSph::ComputeH
/// Compute the value of the smoothing length of particle 'i' by iterating the relation :
/// h = h_fac*(m/rho)^(1/ndim).
/// Uses the previous value of h as a starting guess and then uses either a Newton-Rhapson solver,
/// or fixed-point iteration, to converge on the correct value of h.  The maximum tolerance used
/// for deciding whether the iteration has converged is given by the 'h_converge' parameter.
//=================================================================================================
template <int ndim, template<int> class kernelclass>
int GradhSph<ndim, kernelclass>::ComputeH
 (const int i,                         ///< [in] id of particle
  const int Nneib,                     ///< [in] No. of potential neighbours
  const FLOAT hmax,                    ///< [in] Max. h permitted by neib list
  FLOAT *m,                            ///< [in] Array of neib. masses
  FLOAT *mu,                           ///< [in] Array of m*u (not needed here)
  FLOAT *drsqd,                        ///< [in] Array of neib. distances squared
  FLOAT *gpot,                         ///< [in] Array of neib. grav. potentials
  SphParticle<ndim> &part,             ///< [inout] Particle i data
  Nbody<ndim> *nbody)                  ///< [in] Main N-body object
{
  int j;                               // Neighbour id
  int k;                               // Dimension counter
  int iteration = 0;                   // h-rho iteration counter
  int iteration_max = 30;              // Max. no of iterations
  FLOAT dr[ndim];                      // Relative position vector
  FLOAT h_lower_bound = (FLOAT) 0.0;   // Lower bound on h
  FLOAT h_upper_bound = hmax;          // Upper bound on h
  FLOAT invhsqd;                       // (1 / h)^2
  FLOAT ssqd;                          // Kernel parameter squared, (r/h)^2
  GradhSphParticle<ndim>& parti = static_cast<GradhSphParticle<ndim>& > (part);


  // If there are sink particles present, check if the particle is inside one.
  // If so, then adjust the iteration bounds and ensure they are valid (i.e. hmax is large enough)
  if (parti.sinkid != -1) {
    h_lower_bound = hmin_sink;
    //h_lower_bound = nbody->stardata[parti.sinkid].h;  //hmin_sink;
    if (hmax < hmin_sink) return -1;
  }

  // Some basic sanity-checking in case of invalid input into routine
  assert(Nneib > 0);
  assert(hmax > (FLOAT) 0.0);
  assert(parti.itype != dead);
  assert(parti.m > (FLOAT) 0.0);


  // Main smoothing length iteration loop
  //===============================================================================================
  do {

    // Initialise all variables for this value of h
    iteration++;
    parti.invh     = (FLOAT) 1.0/parti.h;
    parti.rho      = (FLOAT) 0.0;
    parti.invomega = (FLOAT) 0.0;
    parti.zeta     = (FLOAT) 0.0;
    parti.hfactor  = pow(parti.invh,ndim);
    invhsqd        = parti.invh*parti.invh;

    // Loop over all nearest neighbours in list to calculate density, omega and zeta.
    //---------------------------------------------------------------------------------------------
    for (j=0; j<Nneib; j++) {
      ssqd           = drsqd[j]*invhsqd;
      parti.rho      += m[j]*kern.w0_s2(ssqd);
      parti.invomega += m[j]*parti.invh*kern.womega_s2(ssqd);
      parti.zeta     += m[j]*kern.wzeta_s2(ssqd);
    }
    //---------------------------------------------------------------------------------------------

    parti.rho      *= parti.hfactor;
    parti.invomega *= parti.hfactor;
    parti.zeta     *= invhsqd;

    // Density must at least equal its self-contribution
    // (failure could indicate neighbour list problem)
    assert(parti.rho >= parti.m*parti.hfactor*kern.w0(0.0));

    if (parti.rho > (FLOAT) 0.0) parti.invrho = (FLOAT) 1.0/parti.rho;


    // If h changes below some fixed tolerance, exit iteration loop
    if (parti.rho > (FLOAT) 0.0 && parti.h > h_lower_bound &&
        fabs(parti.h - h_fac*pow(parti.m*parti.invrho,Sph<ndim>::invndim)) < h_converge) break;

    // Use fixed-point iteration, i.e. h_new = h_fac*(m/rho_old)^(1/ndim), for now.  If this does
    // not converge in a reasonable number of iterations (iteration_max), then assume something is
    // wrong and switch to a bisection method, which should be guaranteed to converge, albeit much
    // more slowly.  (N.B. will implement Newton-Raphson soon)
    //---------------------------------------------------------------------------------------------
    if (iteration < iteration_max) {
      parti.h = h_fac*pow(parti.m*parti.invrho,Sph<ndim>::invndim);
    }
    else if (iteration == iteration_max) {
      parti.h = (FLOAT) 0.5*(h_lower_bound + h_upper_bound);
    }
    else if (iteration < 5*iteration_max) {
      if (parti.rho < small_number || parti.rho*pow(parti.h,ndim) > pow(h_fac,ndim)*parti.m) {
        h_upper_bound = parti.h;
      }
      else {
        h_lower_bound = parti.h;
      }
      parti.h = (FLOAT) 0.5*(h_lower_bound + h_upper_bound);
    }
    else {
      cout << "H ITERATION : " << iteration << "    h : " << parti.h
           << "   rho : " << parti.rho << "   h_upper " << h_upper_bound << "    hmax :  " << hmax
           << "   h_lower : " << h_lower_bound << "    " << parti.hfactor << "    m : " << parti.m
           << "     " << parti.m*parti.hfactor*kern.w0(0.0) << "    " << Nneib << endl;

      string message = "Problem with convergence of h-rho iteration";
      ExceptionHandler::getIstance().raise(message);
    }

    // If the smoothing length is too large for the neighbour list, exit routine and flag neighbour
    // list error in order to generate a larger neighbour list (not properly implemented yet).
    if (parti.h > hmax) return 0;

  } while (parti.h > h_lower_bound && parti.h < h_upper_bound);
  //===============================================================================================


  // Normalise all SPH sums correctly
  parti.h         = max(h_fac*powf(parti.m*parti.invrho, Sph<ndim>::invndim), h_lower_bound);
  parti.invh      = (FLOAT) 1.0/parti.h;
  parti.hfactor   = pow(parti.invh, ndim+1);
  parti.hrangesqd = kernfacsqd*kern.kernrangesqd*parti.h*parti.h;
  parti.div_v     = (FLOAT) 0.0;
  assert(!(isinf(parti.h)) && !(isnan(parti.h)));
  assert(part.h >= h_lower_bound);

  // Calculate the minimum neighbour potential (used later to identify new sinks)
  if (create_sinks == 1) {
    parti.potmin = true;
    for (j=0; j<Nneib; j++) {
      if (gpot[j] > (FLOAT) 1.000000001*parti.gpot &&
          drsqd[j]*invhsqd < kern.kernrangesqd) parti.potmin = false;
    }
  }

  // If there are star particles, compute N-body zeta correction term
  //-----------------------------------------------------------------------------------------------
  if (part.sinkid == -1) {
    if (nbody->nbody_softening == 1) {
      for (j=0; j<nbody->Nstar; j++) {
        invhsqd = pow((FLOAT) 2.0 / (parti.h + nbody->stardata[j].h),2);
        for (k=0; k<ndim; k++) dr[k] = nbody->stardata[j].r[k] - parti.r[k];
        ssqd = DotProduct(dr,dr,ndim)*invhsqd;
        parti.zeta += nbody->stardata[j].m*invhsqd*kern.wzeta_s2(ssqd);
      }
    }
    else {
      invhsqd = (FLOAT) 4.0*parti.invh*parti.invh;
      for (j=0; j<nbody->Nstar; j++) {
        for (k=0; k<ndim; k++) dr[k] = nbody->stardata[j].r[k] - parti.r[k];
        ssqd = DotProduct(dr,dr,ndim)*invhsqd;
        parti.zeta += nbody->stardata[j].m*invhsqd*kern.wzeta_s2(ssqd);
      }
    }
    parti.invomega = (FLOAT) 1.0 + Sph<ndim>::invndim*parti.h*parti.invomega*parti.invrho;
    parti.invomega = (FLOAT) 1.0/parti.invomega;
    parti.zeta = -Sph<ndim>::invndim*parti.h*parti.zeta*parti.invrho*parti.invomega;
  }
  else {
    parti.invomega = (FLOAT) 1.0;
    parti.zeta     = (FLOAT) 0.0;
  }


  // Set important thermal variables here
  ComputeThermalProperties(parti);


  // If h is invalid (i.e. larger than maximum h), then return error code (0)
  if (parti.h <= hmax) return 1;
  else return -1;
}



//=================================================================================================
//  GradhSph::ComputeThermalProperties
/// Compute all thermal properties for grad-h SPH method for given particle.
//=================================================================================================
template <int ndim, template<int> class kernelclass>
void GradhSph<ndim, kernelclass>::ComputeThermalProperties
 (SphParticle<ndim> &part_gen)         ///< [inout] Particle i data
{
  GradhSphParticle<ndim>& part = static_cast<GradhSphParticle<ndim> &> (part_gen);

  part.u       = eos->SpecificInternalEnergy(part);
  part.sound   = eos->SoundSpeed(part);
  //part.press   = eos->Pressure(part);
  part.pfactor = eos->Pressure(part)*part.invrho*part.invrho*part.invomega;

  return;
}



//=================================================================================================
//  GradhSph::ComputeSphHydroForces
/// Compute SPH neighbour force pairs for
/// (i) All neighbour interactions of particle i with i.d. j > i,
/// (ii) Active neighbour interactions of particle j with i.d. j > i
/// (iii) All inactive neighbour interactions of particle i with i.d. j < i.
/// This ensures that all particle-particle pair interactions are only
/// computed once only for efficiency.
//=================================================================================================
template <int ndim, template<int> class kernelclass>
void GradhSph<ndim, kernelclass>::ComputeSphHydroForces
 (const int i,                         ///< [in] id of particle
  const int Nneib,                     ///< [in] No. of neins in neibpart array
  const int *neiblist,                 ///< [in] id of gather neibs in neibpart
  const FLOAT *drmag,                  ///< [in] Distances of gather neighbours
  const FLOAT *invdrmag,               ///< [in] Inverse distances of gather neibs
  const FLOAT *dr,                     ///< [in] Position vector of gather neibs
  SphParticle<ndim> &part,             ///< [inout] Particle i data
  SphParticle<ndim> *neibpart_gen)     ///< [inout] Neighbour particle data
{
  int j;                               // Neighbour list id
  int jj;                              // Aux. neighbour counter
  int k;                               // Dimension counter
  FLOAT alpha_mean;                    // Mean articial viscosity alpha value
  FLOAT draux[ndim];                   // Relative position vector
  FLOAT dvdr;                          // Dot product of dv and dr
  FLOAT wkerni;                        // Value of w1 kernel function for part i
  FLOAT wkernj;                        // Value of w1 kernel function for neighbour j
  FLOAT vsignal;                       // Signal velocity
  FLOAT paux;                          // Aux. pressure force variable
  FLOAT winvrho;                       // 0.5*(wkerni + wkernj)*invrhomean
  GradhSphParticle<ndim>& parti = static_cast<GradhSphParticle<ndim>& > (part);
  GradhSphParticle<ndim>* neibpart = static_cast<GradhSphParticle<ndim>* > (neibpart_gen);

  // Some basic sanity-checking in case of invalid input into routine
  assert(neibpart_gen != NULL);
  assert(parti.itype != dead);

  // Loop over all potential neighbours in the list
  //-----------------------------------------------------------------------------------------------
  for (jj=0; jj<Nneib; jj++) {
    j = neiblist[jj];
    assert(neibpart[j].itype != dead);

    wkerni = parti.hfactor*kern.w1(drmag[jj]*parti.invh);
    wkernj = neibpart[j].hfactor*kern.w1(drmag[jj]*neibpart[j].invh);

    for (k=0; k<ndim; k++) draux[k] = dr[jj*ndim + k];
    dvdr = DotProduct(neibpart[j].v,draux, ndim);
    dvdr -= DotProduct(parti.v,draux,ndim);

    // Add contribution to velocity divergence
    parti.div_v -= neibpart[j].m*dvdr*wkerni;

    // Main SPH pressure force term
    paux = parti.pfactor*wkerni + neibpart[j].pfactor*wkernj;

    // Add dissipation terms (for approaching particle pairs)
    //---------------------------------------------------------------------------------------------
    if (dvdr < (FLOAT) 0.0) {

      winvrho = (FLOAT) 0.25*(wkerni + wkernj)*(parti.invrho + neibpart[j].invrho);
      //winvrho = (FLOAT) (wkerni + wkernj)/(parti.rho + neibpart[j].rho);

      // Artificial viscosity term
      if (avisc == mon97) {
        vsignal    = parti.sound + neibpart[j].sound - beta_visc*alpha_visc*dvdr;
        paux       -= alpha_visc*vsignal*dvdr*winvrho;
        parti.dudt -= neibpart[j].m*alpha_visc*vsignal*dvdr*dvdr*winvrho;
      }
      else if (avisc == mon97mm97) {
        alpha_mean = (FLOAT) 0.5*(parti.alpha + neibpart[j].alpha);
        vsignal    = parti.sound + neibpart[j].sound - beta_visc*alpha_mean*dvdr;
        paux       -= alpha_mean*vsignal*dvdr*winvrho;
        parti.dudt -= neibpart[j].m*alpha_mean*vsignal*dvdr*dvdr*winvrho;
      }

      // Artificial conductivity term
      if (acond == wadsley2008) {
        parti.dudt += neibpart[j].m*dvdr*(neibpart[j].u - parti.u)*
          (parti.invrho*wkerni + neibpart[j].invrho*wkernj);
      }
      else if (acond == price2008) {
        parti.dudt += (FLOAT) 0.5*neibpart[j].m*(parti.u - neibpart[j].u)*
          winvrho*(parti.invrho + neibpart[j].invrho)*
          sqrt(fabs(eos->Pressure(parti) - eos->Pressure(neibpart[j])));
      }

    }
    //---------------------------------------------------------------------------------------------

    // Add total hydro contribution to acceleration for particle i
    for (k=0; k<ndim; k++) parti.a[k] += neibpart[j].m*draux[k]*paux;
    parti.levelneib = max(parti.levelneib,neibpart[j].level);

  }
  //-----------------------------------------------------------------------------------------------

  // Set velocity divergence and compressional heating rate terms
  parti.div_v    *= parti.invrho;
  parti.dudt     *= (FLOAT) 0.5;
  parti.dudt     -= eos->Pressure(parti)*parti.div_v*parti.invrho*parti.invomega;
  parti.dalphadt = (FLOAT) 0.1*parti.sound*(alpha_visc_min - parti.alpha)*parti.invh +
    max(-parti.div_v, (FLOAT) 0.0)*(alpha_visc - parti.alpha);


  return;
}



//=================================================================================================
//  GradhSph::ComputeSphHydroGravForces
/// Compute SPH neighbour force pairs for
/// (i) All neighbour interactions of particle i with i.d. j > i,
/// (ii) Active neighbour interactions of particle j with i.d. j > i
/// (iii) All inactive neighbour interactions of particle i with i.d. j < i.
/// This ensures that all particle-particle pair interactions are only
/// computed once only for efficiency.
//=================================================================================================
template <int ndim, template<int> class kernelclass>
void GradhSph<ndim, kernelclass>::ComputeSphHydroGravForces
 (const int i,                         ///< [in] id of particle
  const int Nneib,                     ///< [in] No. of neins in neibpart array
  int *neiblist,                       ///< [in] id of gather neibs in neibpart
  SphParticle<ndim> &part,             ///< [inout] Particle i data
  SphParticle<ndim> *neibpart_gen)     ///< [inout] Neighbour particle data
{
  int j;                               // Neighbour list id
  int jj;                              // Aux. neighbour counter
  int k;                               // Dimension counter
  FLOAT alpha_mean;                    // Mean articial viscosity alpha value
  FLOAT dr[ndim];                      // Relative position vector
  FLOAT drmag;                         // Distance
  FLOAT dv[ndim];                      // Relative velocity vector
  FLOAT dvdr;                          // Dot product of dv and dr
  FLOAT invdrmag;                      // 1 / distance
  FLOAT wkerni;                        // Value of w1 kernel function
  FLOAT wkernj;                        // Value of w1 kernel function
  FLOAT vsignal;                       // Signal velocity
  FLOAT paux;                          // Aux. pressure force variable
  FLOAT winvrho;                       // 0.5*(wkerni + wkernj)*invrhomean

  GradhSphParticle<ndim>& parti = static_cast<GradhSphParticle<ndim>& > (part);
  GradhSphParticle<ndim>* neibpart = static_cast<GradhSphParticle<ndim>* > (neibpart_gen);
  const FLOAT invhsqdi = parti.invh*parti.invh;


  // Loop over all potential neighbours in the list
  //-----------------------------------------------------------------------------------------------
  for (jj=0; jj<Nneib; jj++) {
    j = neiblist[jj];

    assert(neibpart[j].itype != dead);
    assert(neibpart[j].h > (FLOAT) 0.0);
    assert(neibpart[j].rho > (FLOAT) 0.0);
    assert(neibpart[j].m > (FLOAT) 0.0);

    for (k=0; k<ndim; k++) dr[k] = neibpart[j].r[k] - parti.r[k];
    for (k=0; k<ndim; k++) dv[k] = neibpart[j].v[k] - parti.v[k];
    drmag = sqrt(DotProduct(dr,dr,ndim) + small_number);
    invdrmag = (FLOAT) 1.0/drmag;
    for (k=0; k<ndim; k++) dr[k] *= invdrmag;
    dvdr = DotProduct(dv,dr,ndim);

    wkerni = parti.hfactor*kern.w1(drmag*parti.invh);
    wkernj = neibpart[j].hfactor*kern.w1(drmag*neibpart[j].invh);

    // Main SPH pressure force term
    paux = parti.pfactor*wkerni + neibpart[j].pfactor*wkernj;

    // Add dissipation terms (for approaching particle pairs)
    //---------------------------------------------------------------------------------------------
    if (dvdr < (FLOAT) 0.0) {

      winvrho = (FLOAT) 0.25*(wkerni + wkernj)*(parti.invrho + neibpart[j].invrho);

      // Artificial viscosity term
      if (avisc == mon97) {
        vsignal     = parti.sound + neibpart[j].sound - beta_visc*alpha_visc*dvdr;
        paux       -= alpha_visc*vsignal*dvdr*winvrho;
        parti.dudt -= neibpart[j].m*alpha_visc*vsignal*dvdr*dvdr*winvrho;
      }
      else if (avisc == mon97mm97) {
        alpha_mean  = (FLOAT) 0.5*(parti.alpha + neibpart[j].alpha);
        vsignal     = parti.sound + neibpart[j].sound - beta_visc*alpha_mean*dvdr;
        paux       -= alpha_mean*vsignal*dvdr*winvrho;
        parti.dudt -= neibpart[j].m*alpha_mean*vsignal*dvdr*dvdr*winvrho;
      }

      // Artificial conductivity term
      if (acond == wadsley2008) {
        parti.dudt += neibpart[j].m*dvdr*(neibpart[j].u - parti.u)*
          (parti.invrho*wkerni + neibpart[j].invrho*wkernj);
      }
      else if (acond == price2008) {
        parti.dudt += (FLOAT) 0.5*neibpart[j].m*(parti.u - neibpart[j].u)*
          winvrho*(parti.invrho + neibpart[j].invrho)*
          sqrt(fabs(eos->Pressure(parti) - eos->Pressure(neibpart[j])));
      }

    }
    //---------------------------------------------------------------------------------------------

    // Add total hydro contribution to acceleration for particle i
    for (k=0; k<ndim; k++) parti.a[k] += neibpart[j].m*dr[k]*paux;


    // Compute SPH gravity terms
    //---------------------------------------------------------------------------------------------
    paux = (FLOAT) 0.5*(invhsqdi*kern.wgrav(drmag*parti.invh) + parti.zeta*wkerni +
                        neibpart[j].invh*neibpart[j].invh*kern.wgrav(drmag*neibpart[j].invh) +
                        neibpart[j].zeta*wkernj);
    for (k=0; k<ndim; k++) parti.agrav[k] += neibpart[j].m*dr[k]*paux;
    parti.gpot += (FLOAT) 0.5*neibpart[j].m*(parti.invh*kern.wpot(drmag*parti.invh) +
                                             neibpart[j].invh*kern.wpot(drmag*neibpart[j].invh));

    // Add velocity divergence and neighbour timestep level terms
    parti.div_v -= neibpart[j].m*dvdr*wkerni;
    parti.levelneib = max(parti.levelneib,neibpart[j].level);

  }
  //-----------------------------------------------------------------------------------------------


  // Set velocity divergence and compressional heating rate terms
  parti.div_v   *= parti.invrho;
  parti.dudt    -= eos->Pressure(parti)*parti.div_v*parti.invrho*parti.invomega;
  parti.dalphadt = (FLOAT) 0.1*parti.sound*(alpha_visc_min - parti.alpha)*
    parti.invh + max(parti.div_v,(FLOAT) 0.0)*(alpha_visc - parti.alpha);


  return;
}



//=================================================================================================
//  GradhSph::ComputeSphGravForces
/// Compute SPH neighbour force pairs for
/// (i) All neighbour interactions of particle i with i.d. j > i,
/// (ii) Active neighbour interactions of particle j with i.d. j > i
/// (iii) All inactive neighbour interactions of particle i with i.d. j < i.
/// This ensures that particle-particle pair interactions are computed once only for efficiency.
//=================================================================================================
template <int ndim, template<int> class kernelclass>
void GradhSph<ndim, kernelclass>::ComputeSphGravForces
 (const int i,                         ///< [in] id of particle
  const int Nneib,                     ///< [in] No. of neins in neibpart array
  int *neiblist,                       ///< [in] id of gather neibs in neibpart
  SphParticle<ndim> &part,             ///< [inout] Particle i data
  SphParticle<ndim> *neibpart_gen)     ///< [inout] Neighbour particle data
{
  int j;                               // Neighbour list id
  int jj;                              // Aux. neighbour counter
  int k;                               // Dimension counter
  FLOAT dr[ndim];                      // Relative position vector
  FLOAT drmag;                         // Distance
  //FLOAT dv[ndim];                      // Relative velocity vector
  //FLOAT dvdr;                          // Dot product of dv and dr
  FLOAT invdrmag;                      // 1 / distance
  FLOAT gaux;                          // Aux. grav. potential variable
  FLOAT paux;                          // Aux. pressure force variable
  GradhSphParticle<ndim>& parti = static_cast<GradhSphParticle<ndim>& > (part);
  GradhSphParticle<ndim>* neibpart = static_cast<GradhSphParticle<ndim>* > (neibpart_gen);


  // Loop over all potential neighbours in the list
  //-----------------------------------------------------------------------------------------------
  for (jj=0; jj<Nneib; jj++) {
    j = neiblist[jj];
    assert(neibpart[j].itype != dead);

    for (k=0; k<ndim; k++) dr[k] = neibpart[j].r[k] - parti.r[k];
    //for (k=0; k<ndim; k++) dv[k] = neibpart[j].v[k] - parti.v[k];
    //dvdr = DotProduct(dv,dr,ndim);
    drmag = sqrt(DotProduct(dr,dr,ndim) + small_number);
    invdrmag = (FLOAT) 1.0/drmag;
    for (k=0; k<ndim; k++) dr[k] *= invdrmag;

    // Main SPH gravity terms
    //---------------------------------------------------------------------------------------------
    paux = (FLOAT) 0.5*(parti.invh*parti.invh*kern.wgrav(drmag*parti.invh) +
                        parti.zeta*parti.hfactor*kern.w1(drmag*parti.invh) +
                        neibpart[j].invh*neibpart[j].invh*kern.wgrav(drmag*neibpart[j].invh) +
                        neibpart[j].zeta*neibpart[j].hfactor*kern.w1(drmag*neibpart[j].invh));
    gaux = (FLOAT) 0.5*(parti.invh*kern.wpot(drmag*parti.invh) +
                        neibpart[j].invh*kern.wpot(drmag*neibpart[j].invh));

    // Add total hydro contribution to acceleration for particle i
    for (k=0; k<ndim; k++) parti.agrav[k] += neibpart[j].m*dr[k]*paux;
    parti.gpot += neibpart[j].m*gaux;

  }

  //===============================================================================================

  return;
}



//=================================================================================================
//  GradhSph::ComputeDirectGravForces
/// Compute the contribution to the total gravitational force of particle 'i'
/// due to 'Nneib' neighbouring particles in the list 'neiblist'.
//=================================================================================================
template <int ndim, template<int> class kernelclass>
void GradhSph<ndim, kernelclass>::ComputeDirectGravForces
 (const int i,                         ///< id of particle
  const int Ndirect,                   ///< No. of nearby 'gather' neighbours
  int *directlist,                     ///< id of gather neighbour in neibpart
  SphParticle<ndim> &part,             ///< Particle i data
  SphParticle<ndim> *sph_gen)          ///< Neighbour particle data
{
  int j;                               // Neighbour particle id
  int jj;                              // Aux. neighbour loop counter
  int k;                               // Dimension counter
  FLOAT dr[ndim];                      // Relative position vector
  FLOAT drsqd;                         // Distance squared
  FLOAT invdrmag;                      // 1 / distance
  FLOAT invdr3;                        // 1 / dist^3
  GradhSphParticle<ndim>& parti = static_cast<GradhSphParticle<ndim>& > (part);
  GradhSphParticle<ndim>* sphdata = static_cast<GradhSphParticle<ndim>* > (sph_gen);


  // Loop over all neighbouring particles in list
  //-----------------------------------------------------------------------------------------------
  for (jj=0; jj<Ndirect; jj++) {
    j = directlist[jj];
    assert(sphdata[j].itype != dead);

    for (k=0; k<ndim; k++) dr[k] = sphdata[j].r[k] - parti.r[k];
    drsqd    = DotProduct(dr,dr,ndim) + small_number;
    invdrmag = (FLOAT) 1.0/sqrt(drsqd);
    invdr3   = invdrmag*invdrmag*invdrmag;

    // Add contribution to current particle
    for (k=0; k<ndim; k++) parti.agrav[k] += sphdata[j].m*dr[k]*invdr3;
    parti.gpot += sphdata[j].m*invdrmag;

    // Sanity-checkt to ensure particles are really un-softened direct-sum neighbours
    assert(drsqd >= parti.hrangesqd && drsqd >= sphdata[j].hrangesqd);

  }
  //-----------------------------------------------------------------------------------------------

  return;
}



//=================================================================================================
//  GradhSph::ComputeStarGravForces
/// Computes contribution of gravitational force and potential due to stars.
//=================================================================================================
template <int ndim, template<int> class kernelclass>
void GradhSph<ndim, kernelclass>::ComputeStarGravForces
 (const int N,                         ///< [in] No. of stars
  NbodyParticle<ndim> **nbodydata,     ///< [in] Array of star pointers
  SphParticle<ndim> &part)             ///< [inout] SPH particle reference
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
  GradhSphParticle<ndim>& parti = static_cast<GradhSphParticle<ndim>& > (part);

  // Loop over all stars and add each contribution
  //-----------------------------------------------------------------------------------------------
  for (j=0; j<N; j++) {

    if (fixed_sink_mass) ms = msink_fixed;
    else ms = nbodydata[j]->m;

    for (k=0; k<ndim; k++) dr[k] = nbodydata[j]->r[k] - parti.r[k];
    //for (k=0; k<ndim; k++) dv[k] = nbodydata[j]->v[k] - parti.v[k];
    drsqd    = DotProduct(dr,dr,ndim) + small_number;
    drmag    = sqrt(drsqd);
    invdrmag = (FLOAT) 1.0/drmag;
    invhmean = (FLOAT) 2.0/(parti.h + nbodydata[j]->h);
    //drdt     = DotProduct(dv,dr,ndim)*invdrmag;
    paux     = ms*invhmean*invhmean*kern.wgrav(drmag*invhmean)*invdrmag;

    // Add total hydro contribution to acceleration for particle i
    for (k=0; k<ndim; k++) parti.agrav[k] += paux*dr[k];
    //for (k=0; k<ndim; k++) parti.adot[k] += paux*dv[k] - (FLOAT) 3.0*paux*drdt*invdrmag*dr[k] +
    //  (FLOAT) 2.0*twopi*ms*drdt*kern.w0(drmag*invhmean)*powf(invhmean,ndim)*invdrmag*dr[k];
    parti.gpot += ms*invhmean*kern.wpot(drmag*invhmean);

    assert(invhmean > (FLOAT) 0.0);

  }
  //-----------------------------------------------------------------------------------------------

  return;
}



template class GradhSph<1, M4Kernel>;
template class GradhSph<2, M4Kernel>;
template class GradhSph<3, M4Kernel>;
template class GradhSph<1, QuinticKernel>;
template class GradhSph<2, QuinticKernel>;
template class GradhSph<3, QuinticKernel>;
template class GradhSph<1, GaussianKernel>;
template class GradhSph<2, GaussianKernel>;
template class GradhSph<3, GaussianKernel>;
template class GradhSph<1, TabulatedKernel>;
template class GradhSph<2, TabulatedKernel>;
template class GradhSph<3, TabulatedKernel>;
