//=================================================================================================
//  GodunovSph.cpp
//  Contains all functions for calculating Godunov SPH (Inutsuka 2002) terms.
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
#include "SphKernel.h"
#include "SphParticle.h"
#include "NbodyParticle.h"
#include "Parameters.h"
#include "EOS.h"
#include "RiemannSolver.h"
#include "InlineFuncs.h"
using namespace std;



//=================================================================================================
//  GodunovSph::GodunovSph
/// GodunovSph class contructor
//=================================================================================================
template <int ndim, template<int> class kernelclass>
GodunovSph<ndim, kernelclass >::GodunovSph(int hydro_forces_aux,
	    int self_gravity_aux, FLOAT alpha_visc_aux, FLOAT beta_visc_aux,
	    FLOAT h_fac_aux, FLOAT h_converge_aux, aviscenum avisc_aux,
	    acondenum acond_aux, tdaviscenum tdavisc_aux, string gas_eos_aux,
            string KernelName):
  Sph<ndim>(hydro_forces_aux,self_gravity_aux, alpha_visc_aux, beta_visc_aux,
	    h_fac_aux, h_converge_aux, avisc_aux, acond_aux, tdavisc_aux,
            gas_eos_aux, KernelName, sizeof(GodunovSphParticle<ndim>)),
  kern(kernelclass<ndim>(KernelName))
{
  this->kernp = &kern;
  this->kernfac = sqrttwo;
  this->kernfacsqd = (FLOAT) 2.0;
}



//=================================================================================================
//  GodunovSph::~GodunovSph
/// GodunovSph class destructor
//=================================================================================================
template <int ndim, template<int> class kernelclass>
GodunovSph<ndim, kernelclass >::~GodunovSph()
{
}


//=================================================================================================
//  GodunovSph::AllocateMemory
/// Allocate main SPH particle array.  Estimates the maximum number of boundary ghost particles
/// assuming a roughly uniform depth of ghosts at each boundary.
//=================================================================================================
template <int ndim, template<int> class kernelclass>
void GodunovSph<ndim, kernelclass>::AllocateMemory(int N)
{
  debug2("[GodunovSph::AllocateMemory]");

  if (N > Nsphmax || !allocated) {
    if (allocated) DeallocateMemory();

    // Set conservative estimate for maximum number of particles, assuming
    // extra space required for periodic ghost particles
    if (Nsphmax < N) {
      Nsphmax = 2*(int) powf(powf((FLOAT) N,invndim) + (FLOAT) 6.0*kernp->kernrange,ndim);
    }

    iorder    = new int[Nsphmax];
    rsph      = new FLOAT[ndim*Nsphmax];
    sphdata   = new struct GodunovSphParticle<ndim>[Nsphmax];
    allocated = true;
    sphdata_unsafe = sphdata;
  }

  return;
}



//=================================================================================================
//  GodunovSph::DeallocateMemory
/// Deallocate main array containing SPH particle data.
//=================================================================================================
template <int ndim, template<int> class kernelclass>
void GodunovSph<ndim, kernelclass>::DeallocateMemory(void)
{
  debug2("[GodunovSph::DeallocateMemory]");

  if (allocated) {
    delete[] sphdata;
    delete[] rsph;
    delete[] iorder;
  }
  allocated = false;

  return;
}



//=================================================================================================
//  GodunovSph::DeleteDeadParticles
/// Delete 'dead' (e.g. accreted) SPH particles from the main arrays.
//=================================================================================================
template <int ndim, template<int> class kernelclass>
void GodunovSph<ndim, kernelclass>::DeleteDeadParticles(void)
{
  int i;                               // Particle counter
  int itype;                           // ..
  int Ndead = 0;                       // No. of 'dead' particles
  int ilast = Nsph;                    // Aux. counter of last free slot

  debug2("[GodunovSph::DeleteDeadParticles]");

  // Determine new order of particles in arrays.
  // First all live particles and then all dead particles
  for (i=0; i<Nsph; i++) {
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
  Nsph -= Ndead;
  Ntot -= Ndead;
  for (i=0; i<Nsph; i++) iorder[i] = i;


  return;
}



//=================================================================================================
//  GodunovSph::ReorderParticles
/// Delete selected SPH particles from the main arrays.
//=================================================================================================
template <int ndim, template<int> class kernelclass>
void GodunovSph<ndim, kernelclass>::ReorderParticles(void)
{
  int i;                                   // Particle counter
  GodunovSphParticle<ndim> *sphdataaux;    // Aux. SPH particle array

  sphdataaux = new GodunovSphParticle<ndim>[Nsph];

  for (i=0; i<Nsph; i++) sphdataaux[i] = sphdata[i];
  for (i=0; i<Nsph; i++) sphdata[i] = sphdataaux[iorder[i]];

  delete[] sphdataaux;

  return;
}



//=============================================================================
//  GodunovSph::ComputeH
/// Compute the value of the smoothing length of particle 'i' by iterating
/// the relation : h = h_fac*(m/rho)^(1/ndim).
/// Uses the previous value of h as a starting guess and then uses either
/// a Newton-Rhapson solver, or fixed-point iteration, to converge on the
/// correct value of h.  The maximum tolerance used for deciding whether the
/// iteration has converged is given by the 'h_converge' parameter.
//=============================================================================
template <int ndim, template<int> class kernelclass>
int GodunovSph<ndim, kernelclass >::ComputeH
(const int i,                                 // id of particle
 const int Nneib,                             // No. of potential neighbours
 const FLOAT hmax,                            // Maximum allowed smoothing length
 FLOAT *m,                              // Array of neib. masses
 FLOAT *mu,                             // Array of m*u (not needed here)
 FLOAT *drsqd,                          // Array of neib. distances (squared)
 FLOAT *gpot,                           // Grav. potential
 SphParticle<ndim> &part,              // Particle i data
 Nbody<ndim> *nbody)                    // Pointer to N-body object
{
  int j;                            // Neighbour id
  int iteration = 0;                // h-rho iteration counter
  int iteration_max = 30;           // Max. no of iterations
  FLOAT h_lower_bound = 0.0;        // Lower bound on h
  FLOAT h_upper_bound = hmax;       // Upper bound on h
  FLOAT invhsqd;                    // (1 / h)^2
  FLOAT ssqd;                       // Kernel parameter squared, (r/h)^2

  GodunovSphParticle<ndim>& parti = static_cast<GodunovSphParticle<ndim>& > (part);


  // If there are sink particles present, check if the particle is inside one
  if (parti.sinkid != -1) h_lower_bound = hmin_sink;


  // Main smoothing length iteration loop
  //===========================================================================
  do {

    // Initialise all variables for this value of h
    iteration++;
    parti.invh = (FLOAT) 1.0/parti.h;
    parti.rho = (FLOAT) 0.0;
    parti.hfactor = pow(parti.invh,ndim);
    invhsqd = parti.invh*parti.invh;

    // Loop over all nearest neighbours in list to calculate
    // density, omega and zeta.
    //-------------------------------------------------------------------------
    for (j=0; j<Nneib; j++) {
      ssqd = drsqd[j]*invhsqd;
      parti.rho += m[j]*kern.w0_s2(ssqd);
    }
    //-------------------------------------------------------------------------

    parti.rho *= parti.hfactor;

    if (parti.rho > (FLOAT) 0.0) parti.invrho = (FLOAT) 1.0/parti.rho;

    // If h changes below some fixed tolerance, exit iteration loop
    if (parti.rho > (FLOAT) 0.0 && parti.h > h_lower_bound &&
    		fabs(parti.h - h_fac*pow(parti.m*parti.invrho,
    				Sph<ndim>::invndim)) < h_converge) break;

    // Use fixed-point iteration, i.e. h_new = h_fac*(m/rho_old)^(1/ndim),
    // for now.  If this does not converge in a reasonable number of
    // iterations (iteration_max), then assume something is wrong and switch
    // to a bisection method, which should be guaranteed to converge,
    // albeit much more slowly.  (N.B. will implement Newton-Raphson soon)
    //-------------------------------------------------------------------------
    if (iteration < iteration_max)
      parti.h = h_fac*pow(parti.m*parti.invrho,Sph<ndim>::invndim);

    else if (iteration == iteration_max)
      parti.h = (FLOAT) 0.5*(h_lower_bound + h_upper_bound);

    else if (iteration < 5*iteration_max) {
      if (parti.rho < small_number ||
	  parti.rho*pow(parti.h,ndim) > pow(h_fac,ndim)*parti.m)
        h_upper_bound = parti.h;
      else
        h_lower_bound = parti.h;
      parti.h = (FLOAT) 0.5*(h_lower_bound + h_upper_bound);
    }
    else {
      cout << "H ITERATION : " << iteration << "    " << parti.h << "    "
	   << parti.rho << "    " << h_upper_bound << "     " << hmax << "   "
	   << h_lower_bound << "    " << parti.hfactor << "     "
	   << parti.m*parti.hfactor*kern.w0(0.0) << endl;
      cout << "rp : " << parti.r[0] << "     " << parti.v[0]
           << "    " << parti.a[0] << endl;
      string message = "Problem with convergence of h-rho iteration";
      ExceptionHandler::getIstance().raise(message);
    }

    // If the smoothing length is too large for the neighbour list, exit
    // routine and flag neighbour list error in order to generate a larger
    // neighbour list (not properly implemented yet).
    if (parti.h > hmax) return 0;

  } while (parti.h > h_lower_bound && parti.h < h_upper_bound);
  //===========================================================================


  // Normalise all SPH sums correctly
  parti.h = max(h_fac*pow(parti.m*parti.invrho,Sph<ndim>::invndim),
                h_lower_bound);
  parti.invh = (FLOAT) 1.0/parti.h;
  parti.hfactor = pow(parti.invh,ndim+1);
  parti.hrangesqd = kernfacsqd*kern.kernrangesqd*parti.h*parti.h;
  parti.div_v = (FLOAT) 0.0;

  // Set important thermal variables here
  ComputeThermalProperties(parti);

  // Calculate the minimum neighbour potential
  // (used later to identify new sinks)
  if (create_sinks == 1) {
    parti.potmin = true;
    for (j=0; j<Nneib; j++)
      if (gpot[j] > 1.000000001*parti.gpot &&
	  drsqd[j]*invhsqd < kern.kernrangesqd) parti.potmin = false;
  }

  // If h is invalid (i.e. larger than maximum h), then return error code (0)
  if (parti.h <= hmax) return 1;
  else return -1;
}



//=============================================================================
//  GodunovSph::ComputeThermalProperties
/// Compute all thermal properties for grad-h SPH method for given particle.
//=============================================================================
template <int ndim, template<int> class kernelclass>
void GodunovSph<ndim, kernelclass>::ComputeThermalProperties
(SphParticle<ndim> &part_gen)        ///< [inout] Particle i data
{
  GodunovSphParticle<ndim>& part =
    static_cast<GodunovSphParticle<ndim> &> (part_gen);

  part.u = eos->SpecificInternalEnergy(part);
  part.sound = eos->SoundSpeed(part);
  part.press = eos->Pressure(part);
  part.pfactor = part.press*part.invrho*part.invrho;

  return;
}



//=============================================================================
//  GodunovSph::ComputeSphHydroForces
/// Compute SPH neighbour force pairs for
/// (i) All neighbour interactions of particle i with i.d. j > i,
/// (ii) Active neighbour interactions of particle j with i.d. j > i
/// (iii) All inactive neighbour interactions of particle i with i.d. j < i.
/// This ensures that all particle-particle pair interactions are only
/// computed once only for efficiency.
//=============================================================================
template <int ndim, template<int> class kernelclass>
void GodunovSph<ndim, kernelclass >::ComputeSphHydroForces
 (const int i,                         ///< [in] id of particle
  const int Nneib,                     ///< [in] No. of neins in neibpart array
  const int *neiblist,                 ///< [in] id of gather neibs in neibpart
  const FLOAT *drmag,                  ///< [in] Distances of gather neighbours
  const FLOAT *invdrmag,               ///< [in] Inverse distances of gather neibs
  const FLOAT *dr,                     ///< [in] Position vector of gather neibs
  SphParticle<ndim> &part,             ///< [inout] Particle i data
  SphParticle<ndim> *neibpart_gen)     ///< [inout] Neighbour particle data
{
  int j;                            // Neighbour list id
  int jj;                           // Aux. neighbour counter
  int k;                            // Dimension counter
  FLOAT draux[ndim];                // Relative position vector
  FLOAT dv[ndim];                   // Relative velocity vector
  FLOAT dvdr;                       // Dot product of dv and dr
  FLOAT wkerni;                     // Value of w1 kernel function
  FLOAT wkernj;                     // Value of w1 kernel function
  FLOAT paux;                       // Aux. pressure force variable

  FLOAT Aij;
  FLOAT Bij;
  FLOAT Cij;                        // ..
  FLOAT Dij;                        // ..
  FLOAT Vsqdi;                      // ..
  FLOAT Vsqdj;                      // ..
  FLOAT Vprimei;
  FLOAT Vprimej;
  FLOAT Sij;                        // ..
  FLOAT pl,pr;                      // ..
  FLOAT rhol,rhor;                  // ..
  FLOAT vl,vr;                      // ..
  FLOAT pstar;                      // ..
  FLOAT vstar;                      // ..
  //FLOAT vtemp[ndim];                // ..
  FLOAT hconv = pow(invsqrttwo,ndim+1);

  GodunovSphParticle<ndim>& parti = static_cast<GodunovSphParticle<ndim>& > (part);
  GodunovSphParticle<ndim>* neibpart = static_cast<GodunovSphParticle<ndim>* > (neibpart_gen);

  string interpolation = "linear"; //"linear";

  // Loop over all potential neighbours in the list
  //---------------------------------------------------------------------------
  for (jj=0; jj<Nneib; jj++) {
    j = neiblist[jj];

    wkerni = hconv*parti.hfactor*kern.w1(invsqrttwo*drmag[jj]*parti.invh);
    wkernj = hconv*neibpart[j].hfactor*kern.w1(invsqrttwo*drmag[jj]*neibpart[j].invh);

    for (k=0; k<ndim; k++) draux[k] = dr[jj*ndim + k];
    for (k=0; k<ndim; k++) dv[k] = neibpart[j].v[k] - parti.v[k];
    dvdr = DotProduct(dv,draux,ndim);

    // Interpolate quantites between left and right states
    if (interpolation == "linear") {
      Cij = (neibpart[j].invrho - parti.invrho)*invdrmag[jj];
      Dij = (FLOAT) 0.5*(parti.invrho + neibpart[j].invrho);
      Vsqdi = (FLOAT) 0.25*parti.h*parti.h*Cij*Cij + Dij*Dij;
      Vsqdj = (FLOAT) 0.25*neibpart[j].h*neibpart[j].h*Cij*Cij + Dij*Dij;
      Sij = (FLOAT) 0.25*Cij*Dij*(parti.h*parti.h/Vsqdi +
        neibpart[j].h*neibpart[j].h/Vsqdj);
    }
    else { //if (interpolation == "cubic") {
      Vprimei = -DotProduct(parti.gradrho,draux,ndim)*
        parti.invrho*parti.invrho;
      Vprimej = -DotProduct(neibpart[j].gradrho,draux,ndim)*
        neibpart[j].invrho*neibpart[j].invrho;
      if (Vprimei*Vprimej > 0.0) {
        Aij = -2.0*(neibpart[j].invrho - parti.invrho)*pow(invdrmag[jj],3) +
          (Vprimei + Vprimej)*invdrmag[jj]*invdrmag[jj];
        Bij = 0.5*(Vprimej - Vprimei)*invdrmag[jj];
        Cij = 1.5*(neibpart[j].invrho - parti.invrho)*invdrmag[jj] - 0.25*(Vprimei + Vprimej);
        Dij = 0.5*(neibpart[j].invrho + parti.invrho) - 0.125*(Vprimej - Vprimei)*drmag[jj];
        Vsqdi = 15.0*pow(parti.h,6)*Aij*Aij/64.0 +
          3.0*pow(parti.h,4)*(2.0*Aij*Cij + Bij*Bij)/16.0 +
          0.25*parti.h*parti.h*(2.0*Bij*Dij + Cij*Cij) + Dij*Dij;
        Vsqdj = 15.0*pow(neibpart[j].h,6)*Aij*Aij/64.0 +
          3.0*pow(neibpart[j].h,4)*(2.0*Aij*Cij + Bij*Bij)/16.0 +
          0.25*neibpart[j].h*neibpart[j].h*(2.0*Bij*Dij + Cij*Cij) + Dij*Dij;
        Sij = 0.5*((15.0*pow(parti.h,6)*Aij*Bij/32.0 +
          3.0*pow(parti.h,4)*(Aij*Dij + Bij*Cij)/8.0 + 0.5*parti.h*parti.h*Cij*Dij)/Vsqdi +
          (15.0*pow(neibpart[j].h,6)*Aij*Bij/32.0 +
           3.0*pow(neibpart[j].h,4)*(Aij*Dij + Bij*Cij)/8.0 +
           0.5*neibpart[j].h*neibpart[j].h*Cij*Dij)/Vsqdi);
      }
      else {
        Cij = (neibpart[j].invrho - parti.invrho)*invdrmag[jj];
        Dij = (FLOAT) 0.5*(parti.invrho + neibpart[j].invrho);
        Vsqdi = (FLOAT) 0.25*parti.h*parti.h*Cij*Cij + Dij*Dij;
        Vsqdj = (FLOAT) 0.25*neibpart[j].h*neibpart[j].h*Cij*Cij + Dij*Dij;
        Sij = (FLOAT) 0.25*Cij*Dij*(parti.h*parti.h/Vsqdi + neibpart[j].h*neibpart[j].h/Vsqdj);
      }
    }
    assert(fabs(Sij) < 0.5*drmag[jj]);
    Sij = -Sij;

    // Initialise the LHS and RHS of the Riemann problem
    InitialiseRiemannProblem(parti,neibpart[j],draux,drmag[jj],Sij,dvdr,parti.sound,
                             neibpart[j].sound,pl,pr,rhol,rhor,vl,vr);

    // Now solve Riemann problem and return intermediate state variables
    riemann->SolveRiemannProblem(pl,pr,rhol,rhor,parti.sound,
                                 neibpart[j].sound,vl,vr,pstar,vstar);

    // Main SPH pressure force term
    paux = pstar*(Vsqdi*wkerni + Vsqdj*wkernj);

    // Add total hydro contribution to acceleration for particle i
    for (k=0; k<ndim; k++) parti.a[k] += neibpart[j].m*draux[k]*paux;

    // If neighbour is also active, add contribution to force here
    for (k=0; k<ndim; k++) neibpart[j].a[k] -= parti.m*draux[k]*paux;

  }
  //---------------------------------------------------------------------------


  return;
}



//=============================================================================
//  GodunovSph::InitialiseRiemannProblem
/// ..
//=============================================================================
template <int ndim, template<int> class kernelclass>
void GodunovSph<ndim, kernelclass>::InitialiseRiemannProblem
(GodunovSphParticle<ndim>& partl,
 GodunovSphParticle<ndim>& partr,
 FLOAT draux[ndim],
 FLOAT drmag,
 FLOAT Sij,
 FLOAT dvdr,
 FLOAT soundl,
 FLOAT soundr,
 FLOAT &pl,
 FLOAT &pr,
 FLOAT &rhol,
 FLOAT &rhor,
 FLOAT &vl,
 FLOAT &vr)
{


  //int k;
  //FLOAT f;
  //FLOAT R;
  //FLOAT deltal;
  //FLOAT deltar;
  //FLOAT limiter;
  //FLOAT vec1[ndim];
  //FLOAT vec2[ndim];
  //FLOAT limiterl,limiterr,vlmin,vlmax,vrmin,vrmax;
  //FLOAT R1,R2,R3;
  //FLOAT deltamean;

  // 1st-order approximation for initialising Riemann problem
  pl = partl.press;
  pr = partr.press;
  rhol = partl.rho;
  rhor = partr.rho;
  vl = DotProduct(partl.v,draux,ndim); //(FLOAT) 0.0;
  vr = DotProduct(partr.v,draux,ndim); //dvdr;
  Sij = 0.0;

  return;
  if (riemann_order == 1) return;

  /*
  FLOAT plorig = pl;
  FLOAT prorig = pr;
  FLOAT rholorig = rhol;
  FLOAT rhororig = rhor;
  FLOAT vlorig = vl;
  FLOAT vrorig = vr;

  //cout << "Orig state; p   : " << pl << "   " << pr << endl;
  //cout << "            rho : " << rhol << "   " << rhor << endl;
  //cout << "            v   : " << vl << "   " << vr << endl;
  //cout << "riemann_order : " << riemann_order << endl;


  // For 2nd-order, extrapolate quantities using SPH gradients.
  //---------------------------------------------------------------------------
  if (riemann_order == 2 && slope_limiter == "none") {

    // Slope limiter for pressure gradients
    deltal = DotProduct(partl.gradP,draux,ndim)*
      (Sij + 0.5*drmag - 0.5*partl.sound*partl.dt);
    deltar = DotProduct(partr.gradP,draux,ndim)*
      (-Sij - 0.5*drmag + 0.5*partr.sound*partr.dt);
    pl += deltal;
    pr += deltar;

    // Slope limiter for pressure gradients
    deltal = DotProduct(partl.gradrho,draux,ndim)*
      (Sij + 0.5*drmag - 0.5*partl.sound*partl.dt);
    deltar = DotProduct(partr.gradrho,draux,ndim)*
      (-Sij - 0.5*drmag + 0.5*partr.sound*partr.dt);
    rhol += deltal;
    rhor += deltar;

    // Slope limiter for pressure gradients
    for (k=0; k<ndim; k++) vec1[k] = DotProduct(partl.gradv[k],draux,ndim);
    for (k=0; k<ndim; k++) vec2[k] = DotProduct(partr.gradv[k],draux,ndim);
    deltal = DotProduct(vec1,draux,ndim)*
      (Sij + 0.5*drmag - 0.5*partl.sound*partl.dt);
    deltar = DotProduct(vec2,draux,ndim)*
      (-Sij - 0.5*drmag + 0.5*partr.sound*partr.dt);
    vl += deltal;
    vr += deltar;
  }
  // ..
  //---------------------------------------------------------------------------
  else if (riemann_order == 2 && slope_limiter == "simple") {

    // Slope limiter for pressure gradients
    deltal = DotProduct(partl.gradP,draux,ndim)*
      (Sij + 0.5*drmag - 0.5*partl.sound*partl.dt);
    deltar = DotProduct(partr.gradP,draux,ndim)*
      (-Sij - 0.5*drmag + 0.5*partr.sound*partr.dt);
    deltamean = pr - pl;

    if (deltal*deltamean > 0.0)
      deltal = min(0.5*fabs(deltamean),fabs(deltal))*sgn(deltal);
    else
      deltal = 0.0;

    if (deltar*deltamean < 0.0)
      deltar = min(0.5*fabs(deltamean),fabs(deltar))*sgn(deltar);
    else
      deltar = 0.0;

    pl += deltal;
    pr += deltar;


    // Slope limiter for density gradients
    deltal = DotProduct(partl.gradrho,draux,ndim)*
      (Sij + 0.5*drmag - 0.5*partl.sound*partl.dt);
    deltar = DotProduct(partr.gradrho,draux,ndim)*
      (-Sij - 0.5*drmag + 0.5*partr.sound*partr.dt);
    deltamean = rhor - rhol;

    if (deltal*deltamean > 0.0)
      deltal = min(0.5*fabs(deltamean),fabs(deltal))*sgn(deltal);
    else
      deltal = 0.0;

    if (deltar*deltamean < 0.0)
      deltar = min(0.5*fabs(deltamean),fabs(deltar))*sgn(deltar);
    else
      deltar = 0.0;

    rhol += deltal;
    rhor += deltar;


    // Slope limiter for velocity gradients
    for (k=0; k<ndim; k++) vec1[k] = DotProduct(partl.gradv[k],draux,ndim);
    for (k=0; k<ndim; k++) vec2[k] = DotProduct(partr.gradv[k],draux,ndim);
    deltal = DotProduct(vec1,draux,ndim)*
      (Sij + 0.5*drmag - 0.5*partl.sound*partl.dt);
    deltar = DotProduct(vec2,draux,ndim)*
      (-Sij - 0.5*drmag + 0.5*partr.sound*partr.dt);
    deltamean = vr - vl;

    if (deltal*deltamean > 0.0)
      deltal = min(0.5*fabs(deltamean),fabs(deltal))*sgn(deltal);
    else
      deltal = 0.0;

    if (deltar*deltamean < 0.0)
      deltar = min(0.5*fabs(deltamean),fabs(deltar))*sgn(deltar);
    else
      deltar = 0.0;

    vl += deltal;
    vr += deltar;

  }
  // For 2nd-order, extrapolate quantities using SPH gradients.
  //---------------------------------------------------------------------------
  else if (riemann_order == 2 && slope_limiter == "I02") {

    // Slope limiter for pressure gradients
    deltal = DotProduct(partl.gradP,draux,ndim)*drmag;
    deltar = DotProduct(partr.gradP,draux,ndim)*drmag;
    pl += 0.5*deltal*(1.0 - min(0.5,partl.sound*partl.dt/drmag));
    pr -= 0.5*deltar*(1.0 - min(0.5,partr.sound*partr.dt/drmag));


    // Slope limiter for pressure gradients
    deltal = DotProduct(partl.gradrho,draux,ndim)*drmag;
    deltar = DotProduct(partr.gradrho,draux,ndim)*drmag;
    rhol += 0.5*deltal*(1.0 - min(0.5,partl.sound*partl.dt/drmag));
    rhor -= 0.5*deltar*(1.0 - min(0.5,partr.sound*partr.dt/drmag));


    // Slope limiter for pressure gradients
    for (k=0; k<ndim; k++) vec1[k] = DotProduct(partl.gradv[k],draux,ndim);
    for (k=0; k<ndim; k++) vec2[k] = DotProduct(partr.gradv[k],draux,ndim);
    deltal = DotProduct(vec1,draux,ndim)*drmag;
    deltar = DotProduct(vec2,draux,ndim)*drmag;
    if (deltal*deltar < 0.0 || fabs(vr - vl) > min(soundl,soundr)) limiter = 0.0;
    else limiter = 1.0;
    vl += 0.5*limiter*deltal*(1.0 - min(0.5,partl.sound*partl.dt/drmag));
    vr -= 0.5*limiter*deltar*(1.0 - min(0.5,partr.sound*partr.dt/drmag));

  }
  // For 2nd-order, extrapolate quantities using SPH gradients.
  //---------------------------------------------------------------------------
  else if (riemann_order == 2 && slope_limiter == "mine") {

    // Slope limiter for pressure gradients
    deltal = DotProduct(partl.gradP,draux,ndim)*drmag;
    deltar = DotProduct(partr.gradP,draux,ndim)*drmag;
    R = min(deltal/(deltar + small_number),
                  deltar/(deltal + small_number));
    R = max(R,0.0);
    limiter = 4.0*R/(R + 1.0)/(R + 1.0);

    pl += 0.5*limiter*deltal*(1.0 - min(1.0,partl.sound*partl.dt/drmag));
    pr -= 0.5*limiter*deltar*(1.0 - min(1.0,partr.sound*partr.dt/drmag));

    R1 = R;

    // Slope limiter for pressure gradients
    deltal = DotProduct(partl.gradrho,draux,ndim)*drmag;
    deltar = DotProduct(partr.gradrho,draux,ndim)*drmag;
    R = min(deltal/(deltar + small_number),
                  deltar/(deltal + small_number));
    R = max(R,0.0);
    limiter = 4.0*R/(R + 1.0)/(R + 1.0);

    rhol += 0.5*limiter*deltal*(1.0 - min(1.0,partl.sound*partl.dt/drmag));
    rhor -= 0.5*limiter*deltar*(1.0 - min(1.0,partr.sound*partr.dt/drmag));

    R2 = R;

    // Slope limiter for pressure gradients
    for (k=0; k<ndim; k++) vec1[k] = DotProduct(partl.gradv[k],draux,ndim);
    for (k=0; k<ndim; k++) vec2[k] = DotProduct(partr.gradv[k],draux,ndim);
    deltal = DotProduct(vec1,draux,ndim)*drmag;
    deltar = DotProduct(vec2,draux,ndim)*drmag;
    R = min(deltal/(deltar + small_number),
                  deltar/(deltal + small_number));
    R = max(R,0.0);
    limiter = 4.0*R/(R + 1.0)/(R + 1.0);

    vl += 0.5*limiter*deltal*(1.0 - min(1.0,partl.sound*partl.dt/drmag));
    vr -= 0.5*limiter*deltar*(1.0 - min(1.0,partr.sound*partr.dt/drmag));

    R3 = R;
  }
  // ..
  //---------------------------------------------------------------------------
  else if (riemann_order == 2 && slope_limiter == "vanleer1979a") {

    // Slope limiter for pressure gradients
    deltal = DotProduct(partl.gradP,draux,ndim)*drmag;
    deltar = DotProduct(partr.gradP,draux,ndim)*drmag;
    deltamean = pr - pl;

    if (deltal*deltamean > 0.0)
      limiterl = 2.0*deltal*deltamean/(deltal + deltamean);
    else
      limiterl = 0.0;

    if (deltar*deltamean > 0.0)
      limiterr = 2.0*deltar*deltamean/(deltar + deltamean);
    else
      limiterr = 0.0;

    pl += 0.5*limiterl*(1.0 - min(0.5,partl.sound*partl.dt/drmag));
    pr -= 0.5*limiterr*(1.0 - min(0.5,partr.sound*partr.dt/drmag));


    // Slope limiter for pressure gradients
    deltal = DotProduct(partl.gradrho,draux,ndim)*drmag;
    deltar = DotProduct(partr.gradrho,draux,ndim)*drmag;
    deltamean = rhor - rhol;

    if (deltal*deltamean > 0.0)
      limiterl = 2.0*deltal*deltamean/(deltal + deltamean);
    else
      limiterl = 0.0;

    if (deltar*deltamean > 0.0)
      limiterr = 2.0*deltar*deltamean/(deltar + deltamean);
    else
      limiterr = 0.0;

    rhol += 0.5*limiterl*(1.0 - min(0.5,partl.sound*partl.dt/drmag));
    rhor -= 0.5*limiterr*(1.0 - min(0.5,partr.sound*partr.dt/drmag));


    // Slope limiter for pressure gradients
    for (k=0; k<ndim; k++) vec1[k] = DotProduct(partl.gradv[k],draux,ndim);
    for (k=0; k<ndim; k++) vec2[k] = DotProduct(partr.gradv[k],draux,ndim);
    deltal = DotProduct(vec1,draux,ndim)*drmag;
    deltar = DotProduct(vec2,draux,ndim)*drmag;
    deltamean = vr - vl;

    if (deltal*deltamean > 0.0)
      limiterl = 2.0*deltal*deltamean/(deltal + deltamean);
    else
      limiterl = 0.0;

    if (deltar*deltamean > 0.0)
      limiterr = 2.0*deltar*deltamean/(deltar + deltamean);
    else
      limiterr = 0.0;

    vl += 0.5*limiterl*(1.0 - min(0.5,partl.sound*partl.dt/drmag));
    vr -= 0.5*limiterr*(1.0 - min(0.5,partr.sound*partr.dt/drmag));

  }
  // Van Leer (1979)b limiter
  //---------------------------------------------------------------------------
  else if (riemann_order == 2 && slope_limiter == "vanleer1979b") {

    // Slope limiter for pressure gradients
    deltal = DotProduct(partl.gradP,draux,ndim)*drmag;
    deltar = DotProduct(partr.gradP,draux,ndim)*drmag;

    if (sgn(pr - pl) == sgn(deltal))
      limiterl = min(fabs(pr - pl),fabs(deltal))*sgn(deltal);
    else
      limiterl = 0.0;

    if (sgn(pr - pl) == sgn(deltar))
      limiterr = min(fabs(pr - pl),fabs(deltar))*sgn(deltar);
    else
      limiterr = 0.0;

    pl += 0.5*limiterl*(1.0 - min(1.0,partl.sound*partl.dt/drmag));
    pr -= 0.5*limiterr*(1.0 - min(1.0,partr.sound*partr.dt/drmag));


    // Slope limiter for pressure gradients
    deltal = DotProduct(partl.gradrho,draux,ndim)*drmag;
    deltar = DotProduct(partr.gradrho,draux,ndim)*drmag;

    if (sgn(rhor - rhol) == sgn(deltal))
      limiterl = min(fabs(rhor - rhol),fabs(deltal))*sgn(deltal);
    else
      limiterl = 0.0;

    if (sgn(rhor - rhol) == sgn(deltar))
      limiterr = min(fabs(rhor - rhol),fabs(deltar))*sgn(deltar);
    else
      limiterr = 0.0;

    rhol += 0.5*limiterl*(1.0 - min(1.0,partl.sound*partl.dt/drmag));
    rhor -= 0.5*limiterr*(1.0 - min(1.0,partr.sound*partr.dt/drmag));


    // Slope limiter for pressure gradients
    for (k=0; k<ndim; k++) vec1[k] = DotProduct(partl.gradv[k],draux,ndim);
    for (k=0; k<ndim; k++) vec2[k] = DotProduct(partr.gradv[k],draux,ndim);
    deltal = DotProduct(vec1,draux,ndim)*drmag;
    deltar = DotProduct(vec2,draux,ndim)*drmag;

    if (sgn(vr - vl) == sgn(deltar))
      limiterl = min(fabs(vr - vl),fabs(deltal))*sgn(deltal);
    else
      limiterl = 0.0;

    if (sgn(vr - vl) == sgn(deltar))
      limiterr = min(fabs(vr - vl),fabs(deltar))*sgn(deltar);
    else
      limiterr = 0.0;

    vl += 0.5*limiterl*(1.0 - min(1.0,partl.sound*partl.dt/drmag));
    vr -= 0.5*limiterr*(1.0 - min(1.0,partr.sound*partr.dt/drmag));

  }
  //---------------------------------------------------------------------------
  else if (riemann_order == 2) {
    cout << "Unrecognised slope limiter" << endl;
    exit(0);
  }


  if (slope_limiter != "none") {
    if ((pr - pl)*(prorig - plorig) < -0.00001 || (rhor - rhol)*(rhororig - rholorig) < -0.00001 || (vr - vl)*(vrorig - vlorig) < -0.00001) {
      cout << "Something wrong with initial values?? : " << endl;
      cout << "R    : " << R1 << "   " << R2 << "   " << R3 << endl;
      cout << "pl   : " << pl << "   " << plorig << endl;
      cout << "pr   : " << pr << "   " << prorig << endl;
      cout << "rhol   : " << rhol << "   " << rholorig << endl;
      cout << "rhor   : " << rhor << "   " << rhororig << endl;
      cout << "vl   : " << vl << "   " << vlorig << endl;
      cout << "vr   : " << vr << "   " << vrorig << endl;
      cout << "limiter : " << limiterl << "   " << limiterr << endl;
      cout << "delta : " << deltal << "    " << deltar << endl;
      cout << "extra : " << (pr - pl)*(prorig - plorig) << "   " << (rhor - rhol)*(rhororig - rholorig) << "    " << (vr - vl)*(vrorig - vlorig) << endl;
      deltal = DotProduct(partl.gradrho,draux,ndim)*drmag;
      deltar = DotProduct(partr.gradrho,draux,ndim)*drmag;
      deltamean = rhororig - rholorig;
      cout << "grads : " << deltal << "    " << deltar << "    " << deltamean << endl;
      exit(0);

    }
  }



  //cout << "Fin. state; p   : " << pl << "   " << pr << endl;
  //cout << "            rho : " << rhol << "   " << rhor << endl;
  //cout << "            v   : " << vl << "   " << vr << endl;
  //cout << "       deltav   : " << deltal << "   " << deltar << endl;
  //cout << "         vmin   : " << vlmin << "   " << vrmin << endl;
  //cout << "         vmax   : " << vlmax << "   " << vrmax << endl;
  //cout << "      limiter   : " << limiterl << "   " << limiterr << endl;
  //cout << "          draux : " << draux[0] << endl;
  //cin >> limiterl;

  //vl -= DotProduct(partl.v,draux,ndim); //(FLOAT) 0.0;
  //vr -= DotProduct(partl.v,draux,ndim); //(FLOAT) 0.0;


  // Check for NaNs or invalid values
  if (pl != pl || pr != pr || vl != vl || vr != vr || rhol != rhol ||
      rhor != rhor || pl <= 0.0 || pr <= 0.0 || rhol <= 0.0 || rhor <= 0.0) {
    cout << "Problem initialising Riemann problem" << endl;
    cout << "p   : " << pl << "    " << pr << endl;
    cout << "v   : " << vl << "    " << vr << endl;
    cout << "rho : " << rhol << "    " << rhor << endl;
    cout << "porig   : " << plorig << "   " << prorig << endl;
    cout << "rhoorig : " << rholorig << "   " << rhororig << endl;
    cout << "vorig   : " << vlorig << "   " << vrorig << endl;
    cout << "riemann_order : " << riemann_order << endl;
    cout << "drmag : " << drmag << "   " << partl.sound*partl.dt << "   " << partr.sound*partr.dt << endl;
    cout << "Sij   : " << Sij << endl;
    cout << "u     : " << partl.u << "   " << partr.u << endl;
    cout << "dudt  : " << partl.dudt << "   " << partr.dudt << endl;
    cout << "dPdr  : " << DotProduct(partl.gradP,draux,ndim)*drmag << "    "
	 << DotProduct(partr.gradP,draux,ndim)*drmag << endl;
    cout << "x     : " << partl.r[0] << "    " << partr.r[0] << endl;
    cout << "dr    : " << draux[0] << "    dPdr : " << partl.gradP[0] << "    " << partr.gradP[1] << endl;
    cout << (Sij + 0.5*drmag - 0.5*partl.sound*partl.dt)/drmag << "    "
	 << (-Sij - 0.5*drmag + 0.5*partr.sound*partr.dt)/drmag << endl;

    exit(0);
  }


  */


  return;
}



//=============================================================================
//  GodunovSph::ComputeSphDerivatives
/// ..
//=============================================================================
template <int ndim, template<int> class kernelclass>
void GodunovSph<ndim, kernelclass >::ComputeSphDerivatives
(const int i,                             ///< id of particle
 const int Nneib,                         ///< No. of neighbours in neibpart array
 int *neiblist,                     ///< id of gather neighbour in neibpart
 FLOAT *drmag,                      ///< Distances of gather neighbours
 FLOAT *invdrmag,                   ///< Inverse distances of gather neibs
 FLOAT *dr,                         ///< Position vector of gather neibs
 SphParticle<ndim> &part,          ///< Particle i data
 SphParticle<ndim> *neibpart_gen)       ///< Neighbour particle data
{
  int j;                            // Neighbour list id
  int jj;                           // Aux. neighbour counter
  int k;                            // Dimension counter
  //int kk;                           // ..
  FLOAT draux[ndim];                // Relative position vector
  FLOAT dv[ndim];                   // Relative velocity vector
  FLOAT wkern;                      // Value of w1 kernel function
  FLOAT dvdr;                       // ..

  GodunovSphParticle<ndim>& parti = static_cast<GodunovSphParticle<ndim>& > (part);
  GodunovSphParticle<ndim>* neibpart = static_cast<GodunovSphParticle<ndim>* > (neibpart_gen);

  parti.div_v = (FLOAT) 0.0;
  for (k=0; k<ndim; k++) {
    parti.gradrho[k] = (FLOAT) 0.0;
    //parti.gradP[k] = (FLOAT) 0.0;
    //for (kk=0; kk<ndim; kk++) parti.gradv[k][kk] = (FLOAT) 0.0;
  }


  // Loop over all potential neighbours in the list
  //---------------------------------------------------------------------------
  for (jj=0; jj<Nneib; jj++) {
    j = neiblist[jj];
    wkern = parti.hfactor*kern.w1(drmag[jj]*parti.invh);

    for (k=0; k<ndim; k++) draux[k] = dr[jj*ndim + k];
    for (k=0; k<ndim; k++) dv[k] = neibpart[j].v[k] - parti.v[k];
    dvdr = DotProduct(dv,draux,ndim);

    parti.div_v -= neibpart[j].m*dvdr*wkern;

    // Compute gradients of density, pressure and velocity for slope limiter
    for (k=0; k<ndim; k++) parti.gradrho[k] -= neibpart[j].m*
      (neibpart[j].rho - parti.rho)*wkern*draux[k];
    //for (k=0; k<ndim; k++) parti.gradP[k] -= neibpart[j].m*
    //  (neibpart[j].press - parti.press)*wkern*draux[k];
    //for (k=0; k<ndim; k++)
    //  for (kk=0; kk<ndim; kk++) parti.gradv[kk][k] -= neibpart[j].m*
    //    dv[kk]*wkern*draux[k];
  }
  //---------------------------------------------------------------------------

  // Normalise summations
  parti.div_v *= parti.invrho;
  for (k=0; k<ndim; k++) {
    parti.gradrho[k] *= parti.invrho;
    //parti.gradP[k] *= parti.invrho;
    //for (kk=0; kk<ndim; kk++) parti.gradv[k][kk] *= parti.invrho;
  }

  return;
}



//=============================================================================
//  GodunovSph::ComputeSphNeibDudt
/// Compute SPH neighbour contributions to dudt for
/// (i) All neighbour interactions of particle i with i.d. j > i,
/// (ii) Active neighbour interactions of particle j with i.d. j > i
/// (iii) All inactive neighbour interactions of particle i with i.d. j < i.
/// This ensures that all particle-particle pair interactions are only
/// computed once only for efficiency.
//=============================================================================
template <int ndim, template<int> class kernelclass>
void GodunovSph<ndim, kernelclass >::ComputeSphNeibDudt
(const int i,                             ///< id of particle
 const int Nneib,                         ///< No. of neighbours in neibpart array
 int *neiblist,                     ///< id of gather neighbour in neibpart
 FLOAT *drmag,                      ///< Distances of gather neighbours
 FLOAT *invdrmag,                   ///< Inverse distances of gather neibs
 FLOAT *dr,                         ///< Position vector of gather neibs
 SphParticle<ndim> &part,          ///< Particle i data
 SphParticle<ndim> *neibpart_gen)       ///< Neighbour particle data
{
  int j;                            // Neighbour list id
  int jj;                           // Aux. neighbour counter
  int k;                            // Dimension counter
  //FLOAT da[ndim];                   // ..
  FLOAT draux[ndim];                // Relative position vector
  FLOAT dv[ndim];                   // Relative velocity vector
  FLOAT dvdr;                       // Dot product of dv and dr
  FLOAT wkerni;                     // Value of w1 kernel function
  FLOAT wkernj;                     // Value of w1 kernel function
  FLOAT uaux;                       // Aux. internal energy variable

  FLOAT Aij;
  FLOAT Bij;
  FLOAT Cij;                        // ..
  FLOAT Dij;                        // ..
  FLOAT Vprimei;
  FLOAT Vprimej;
  FLOAT Vsqdi;                      // ..
  FLOAT Vsqdj;                      // ..
  FLOAT Sij;                        // ..
  FLOAT pl,pr;
  FLOAT rhol,rhor;
  FLOAT vl,vr;
  FLOAT pstar;
  FLOAT vstar;
  FLOAT vhalfi;
  FLOAT vhalfj;
  //FLOAT vtemp[ndim];
  //FLOAT gradi, gradj;

  GodunovSphParticle<ndim>& parti = static_cast<GodunovSphParticle<ndim>& > (part);
  GodunovSphParticle<ndim>* neibpart = static_cast<GodunovSphParticle<ndim>* > (neibpart_gen);

  FLOAT hconv = powf(invsqrttwo,ndim+1);
  string interpolation = "linear"; //"linear";


  // Loop over all potential neighbours in the list
  //---------------------------------------------------------------------------
  for (jj=0; jj<Nneib; jj++) {
    j = neiblist[jj];
    wkerni = hconv*parti.hfactor*kern.w1(invsqrttwo*drmag[jj]*parti.invh);
    wkernj = hconv*neibpart[j].hfactor*
      kern.w1(invsqrttwo*drmag[jj]*neibpart[j].invh);

    for (k=0; k<ndim; k++) draux[k] = dr[jj*ndim + k];
    for (k=0; k<ndim; k++) dv[k] = neibpart[j].v[k] - parti.v[k];
    //for (k=0; k<ndim; k++) da[k] = neibpart[j].a[k] - parti.a[k];
    dvdr = DotProduct(dv,draux,ndim);
    //vhalfi = (FLOAT) 0.5*DotProduct(parti.a,draux,ndim)*parti.dt;
    //vhalfj = dvdr + (FLOAT) 0.5*DotProduct(neibpart[j].a,draux,ndim)*neibpart[j].dt;
    vhalfi = DotProduct(parti.v,draux,ndim) + (FLOAT) 0.5*DotProduct(parti.a,draux,ndim)*parti.dt;
    vhalfj = DotProduct(neibpart[j].v,draux,ndim) + (FLOAT) 0.5*DotProduct(neibpart[j].a,draux,ndim)*neibpart[j].dt;

    // Interpolate quantites between left and right states
    if (interpolation == "linear") {
      Cij = (neibpart[j].invrho - parti.invrho)*invdrmag[jj];
      Dij = (FLOAT) 0.5*(parti.invrho + neibpart[j].invrho);
      Vsqdi = (FLOAT) 0.25*parti.h*parti.h*Cij*Cij + Dij*Dij;
      Vsqdj = (FLOAT) 0.25*neibpart[j].h*neibpart[j].h*Cij*Cij + Dij*Dij;
      Sij = (FLOAT) 0.25*Cij*Dij*(parti.h*parti.h/Vsqdi + neibpart[j].h*neibpart[j].h/Vsqdj);
    }
    else { //if (interpolation == "cubic") {
      Vprimei = -DotProduct(parti.gradrho,draux,ndim)*parti.invrho*parti.invrho;
      Vprimej = -DotProduct(neibpart[j].gradrho,draux,ndim)*neibpart[j].invrho*neibpart[j].invrho;
      if (Vprimei*Vprimej >= 0.0) {
        Aij = -2.0*(neibpart[j].invrho - parti.invrho)*pow(invdrmag[jj],3) +
          (Vprimei + Vprimej)*invdrmag[jj]*invdrmag[jj];
        Bij = 0.5*(Vprimej - Vprimei)*invdrmag[jj];
        Cij = 1.5*(neibpart[j].invrho - parti.invrho)*invdrmag[jj] - 0.25*(Vprimei + Vprimej);
        Dij = 0.5*(neibpart[j].invrho + parti.invrho) - 0.125*(Vprimej - Vprimei)*drmag[jj];
        Vsqdi = 15.0*pow(parti.h,6)*Aij*Aij/64.0 +
          3.0*pow(parti.h,4)*(2.0*Aij*Cij + Bij*Bij)/16.0 +
          0.25*parti.h*parti.h*(2.0*Bij*Dij + Cij*Cij) + Dij*Dij;
        Vsqdj = 15.0*pow(neibpart[j].h,6)*Aij*Aij/64.0 +
          3.0*pow(neibpart[j].h,4)*(2.0*Aij*Cij + Bij*Bij)/16.0 +
          0.25*neibpart[j].h*neibpart[j].h*(2.0*Bij*Dij + Cij*Cij) + Dij*Dij;
        Sij = 0.5*((15.0*pow(parti.h,6)*Aij*Bij/32.0 + 3.0*pow(parti.h,4)*(Aij*Dij + Bij*Cij)/8.0 +
          0.5*parti.h*parti.h*Cij*Dij)/Vsqdi + (15.0*pow(neibpart[j].h,6)*Aij*Bij/32.0 +
          3.0*pow(neibpart[j].h,4)*(Aij*Dij + Bij*Cij)/8.0 +
          0.5*neibpart[j].h*neibpart[j].h*Cij*Dij)/Vsqdi);
      }
      else {
        Cij = (neibpart[j].invrho - parti.invrho)*invdrmag[jj];
        Dij = (FLOAT) 0.5*(parti.invrho + neibpart[j].invrho);
        Vsqdi = (FLOAT) 0.25*parti.h*parti.h*Cij*Cij + Dij*Dij;
        Vsqdj = (FLOAT) 0.25*neibpart[j].h*neibpart[j].h*Cij*Cij + Dij*Dij;
        Sij = (FLOAT) 0.25*Cij*Dij*(parti.h*parti.h/Vsqdi + neibpart[j].h*neibpart[j].h/Vsqdj);
      }
    }
    assert(fabs(Sij) < 0.5*drmag[jj]);


    // Initialise the LHS and RHS of the Riemann problem
    InitialiseRiemannProblem(parti,neibpart[j],draux,drmag[jj],Sij,dvdr,
                             parti.sound,neibpart[j].sound,pl,pr,rhol,rhor,vl,vr);

    // Now solve Riemann problem and return intermediate state variables
    riemann->SolveRiemannProblem(pl,pr,rhol,rhor,parti.sound,neibpart[j].sound,vl,vr,pstar,vstar);

    // Main SPH pressure force term
    uaux = pstar*(Vsqdi*wkerni + Vsqdj*wkernj);

    // Add total hydro contribution to acceleration for particle i
    parti.dudt += neibpart[j].m*uaux*(vstar - vhalfi);
    //parti.dudt += 2.0*neibpart[j].m*parti.press*Vsqdi*wkerni*(vstar - vhalfi);

    // If neighbour is also active, add contribution to force here
    neibpart[j].dudt -= parti.m*uaux*(vstar - vhalfj);
    //neibpart[j].dudt -= 2.0*parti.m*neibpart[j].press*Vsqdj*wkernj*(vstar - vhalfj);

    //cout << "dudt : " << vhalfi << "   " << vstar << "    " << vhalfj << endl;

    /*if (i == 0) {
      cout << "Heating rates : " << i << "    r : " << parti.r[0] << "   Sij : " << Sij/drmag[jj]
	   << "   " << uaux << endl;
      cout << "   vstar : " << vstar << "   vhalf : " << vhalfi << "    " << vhalfj << endl;
      cout << "   vel : " << parti.v[0] << "    " << neibpart[j].v[0] << "    " << vl << "    " << vr << endl;
      cout << "  dvdr : " << parti.gradv[0][0] << "    " << neibpart[j].gradv[0][0] << "     " << parti.gradv[0][0]*drmag[jj] << "    " << neibpart[j].gradv[0][0]*drmag[jj] << endl;
      cout << "   dudt : " << neibpart[j].m*uaux*(vstar - vhalfi) << endl;
      }*/

  }
  //---------------------------------------------------------------------------

  //if (i==0) cin >> k;

  return;
}



template class GodunovSph<1, M4Kernel>;
template class GodunovSph<1, QuinticKernel>;
template class GodunovSph<1, GaussianKernel>;
template class GodunovSph<1, TabulatedKernel>;
template class GodunovSph<2, M4Kernel>;
template class GodunovSph<2, QuinticKernel>;
template class GodunovSph<2, GaussianKernel>;
template class GodunovSph<2, TabulatedKernel>;
template class GodunovSph<3, M4Kernel>;
template class GodunovSph<3, QuinticKernel>;
template class GodunovSph<3, GaussianKernel>;
template class GodunovSph<3, TabulatedKernel>;
