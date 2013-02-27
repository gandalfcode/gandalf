// ============================================================================
// GradhSph.cpp
// ============================================================================


#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <math.h>
#include "Precision.h"
#include "Sph.h"
#include "SphKernel.h"
#include "SphParticle.h"
#include "Parameters.h"
#include "EOS.h"
#include "InlineFuncs.h"
using namespace std;



// ============================================================================
// GradhSph::GradhSph
// ============================================================================
template <typename kernelclass>
GradhSph<kernelclass>::GradhSph(int ndimaux, int vdimaux, int bdimaux):
#if !defined(FIXED_DIMENSIONS)
  Sph(ndimaux, vdimaux, bdimaux),
#endif
  kern (kernelclass(ndimaux))
{
  allocated = false;
  Nsph = 0;
  Nsphmax = 0;
  kernp = &kern;
}



// ============================================================================
// GradhSph::~GradhSph
// ============================================================================
template <typename kernelclass>
GradhSph<kernelclass>::~GradhSph()
{
}



// ============================================================================
// GradhSph::ComputeH
// Compute the value of the smoothing length of particle 'i' by iterating  
// the relation : h = h_fac*(m/rho)^(1/ndim).
// Uses the previous value of h as a starting guess and then uses either 
// a Newton-Rhapson solver, or fixed-point iteration, to converge on the 
// correct value of h.  The maximum tolerance used for deciding whether the 
// iteration has converged is given by the 'h_converge' parameter.
// ============================================================================
template <typename kernelclass>
int GradhSph<kernelclass>::ComputeH(int i, SphParticle &parti, int Nneib,
				    int Nnear, int *nearlist, 
				    SphParticle *neibpart, FLOAT *drmag, 
				    FLOAT *invdrmag, FLOAT *dr)
{
  int j;                                      // Neighbour id
  int jj;                                     // Aux. neighbour counter
  int k;                                      // Dimension counter
  int iteration = 0;                          // h-rho iteration counter
  int iteration_max = 30;                     // Max. no of iterations
  FLOAT draux[ndimmax];                       // Relative position vector
  FLOAT dv[ndimmax];                          // Relative velocity vector
  FLOAT h_max = big_number;                   // Max. allowed value of h
  FLOAT h_lower_bound = 0.0;                  // Lower bound on h
  FLOAT h_upper_bound = big_number;           // Upper bound on h
  FLOAT hfactor;                              // (1 / h)^ndim
  FLOAT hrange;                               // Kernel extent
  FLOAT invrho;                               // 1 / rho
  FLOAT skern;                                // Kernel parameter, r/h


  // Main smoothing length iteration loop
  // --------------------------------------------------------------------------
  do {

    // Initialise all variables for this value of h
    iteration++;
    parti.invh = (FLOAT) 1.0/parti.h;
    parti.rho = (FLOAT) 0.0;
    parti.invomega = (FLOAT) 0.0;
    parti.zeta = (FLOAT) 0.0;
    parti.div_v = (FLOAT) 0.0;
    hrange = kern.kernrange*parti.h;
    hfactor = pow(parti.invh,ndim);

    // Loop over all nearest neighbours in list to calculate 
    // density, omega and div_v
    // ------------------------------------------------------------------------
    for (jj=0; jj<Nnear; jj++) {
      j = nearlist[jj];
      
      for (k=0; k<ndim; k++) draux[k] = dr[jj*ndim + k];
      for (k=0; k<ndim; k++) dv[k] = neibpart[j].v[k] - parti.v[k];
      skern = drmag[jj]*parti.invh;
      
      parti.div_v -= neibpart[j].m*DotProduct(dv,draux,ndim)*
	kern.w1(skern)*hfactor*parti.invh;
      parti.rho += neibpart[j].m*hfactor*kern.w0(skern);
      parti.invomega += neibpart[j].m*hfactor*
	parti.invh*kern.womega(skern);
      parti.zeta += neibpart[j].m*parti.invh*parti.invh*kern.wzeta(skern);
    }
    // ------------------------------------------------------------------------

    if (parti.rho > (FLOAT) 0.0) invrho = (FLOAT) 1.0/parti.rho;

    // If h changes below some fixed tolerance, exit iteration loop
    if (parti.rho > (FLOAT) 0.0 && parti.h > h_lower_bound &&
	fabs(parti.h - h_fac*pow(parti.m*invrho,
				 invndim)) < h_converge) break;

    // Use fixed-point iteration for now.  If this does not converge in a 
    // reasonable number of iterations (iteration_max), then assume something 
    // is wrong and switch to a bisection method, which should be guaranteed 
    // to converge, albeit much more slowly.
    // ------------------------------------------------------------------------
    if (iteration < iteration_max)
      parti.h = h_fac*pow(parti.m*invrho,invndim);

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

    // If the smoothing length is too large for the neighbour list, exit 
    // routine and flag neighbour list error in order to generate a larger
    // neighbour list
    if (parti.h > h_max) return 0;
    
  } while (parti.h > h_lower_bound && parti.h < h_upper_bound);
  // --------------------------------------------------------------------------


  // Normalise all SPH sums correctly
  parti.h = h_fac*pow(parti.m*invrho,invndim);
  parti.invh = (FLOAT) 1.0/parti.h;
  parti.invomega = (FLOAT) 1.0 + invndim*parti.h*parti.invomega*invrho;
  parti.invomega = (FLOAT) 1.0/parti.invomega;
  parti.zeta = -invndim*parti.h*parti.zeta*invrho*parti.invomega;
  parti.div_v *= invrho;

  // Set important thermal variables here
  parti.u = eos->SpecificInternalEnergy(parti);
  parti.sound = eos->SoundSpeed(parti);
  
  return 1;
}



// ============================================================================
// GradhSph::ComputeHydroForces
// ..
// ============================================================================
template <typename kernelclass>
void GradhSph<kernelclass>::ComputeHydroForces(int i, SphParticle &parti,
					       int Nneib, int Nnear, 
					       int *nearlist, 
					       SphParticle *neiblist,
					       FLOAT *drmag, FLOAT *invdrmag, 
					       FLOAT *dr)
{
  int j;
  int jj;
  int k;
  FLOAT draux[ndimmax];
  FLOAT dv[ndimmax];
  FLOAT dvdr;
  FLOAT wkern;
  FLOAT vsignal;
  FLOAT paux,uaux;
  FLOAT hfactor = pow(parti.invh,ndim+1);
  FLOAT hrange = kern.kernrange*parti.h;
  FLOAT invrho = (FLOAT) 1.0/parti.rho;
  FLOAT pfactor = eos->Pressure(parti)*invrho*invrho*parti.invomega;

  // Compute contribution to compressional heating rate
  parti.dudt -= eos->Pressure(parti)*parti.div_v*invrho*parti.invomega;

  // Loop over all potential neighbours in the list
  // ==========================================================================
  for (jj=0; jj<Nnear; jj++) {
    j = nearlist[jj];

    for (k=0; k<ndim; k++) draux[k] = dr[jj*ndim + k];
    for (k=0; k<ndim; k++) dv[k] = neiblist[j].v[k] - parti.v[k];
    dvdr = DotProduct(dv,draux,ndim);
    wkern = hfactor*kern.w1(drmag[jj]*parti.invh);

    // Compute hydro forces
    // ------------------------------------------------------------------------
    if (hydro_forces == 1) {
      paux = pfactor*wkern;

      // Add dissipation terms (for approaching particle-pairs)
      if (dvdr < (FLOAT) 0.0) {
	
	// Artificial viscosity term
	if (avisc == "mon97" || avisc == "pf2010") {
	  vsignal = parti.sound - beta_visc*dvdr;
	  paux -= (FLOAT) 0.5*alpha_visc*vsignal*dvdr*invrho*parti.invomega*wkern;
	  parti.dudt -= (FLOAT) 0.25*neiblist[j].m*alpha_visc*vsignal*dvdr*dvdr*
	    invrho*parti.invomega*wkern;
	  neiblist[j].dudt -= (FLOAT) 0.25*parti.m*alpha_visc*vsignal*dvdr*dvdr*
	    invrho*parti.invomega*wkern;
	}
	
	// Artificial conductivity term
	if (acond == "wadsley2008") {
	  parti.dudt += (FLOAT) 0.5*neiblist[j].m*fabs(dvdr)*
	    (parti.u - neiblist[j].u)*wkern*invrho;
	  neiblist[j].dudt -= (FLOAT) 0.5*parti.m*fabs(dvdr)*
	    (parti.u - neiblist[j].u)*wkern*invrho;
	}
      }
    }

    // Compute gravitational contribution
    // ------------------------------------------------------------------------
    if (self_gravity == 1) {
      paux += parti.invh*parti.invh*kern.wgrav(drmag[j]*parti.invh) +
	parti.zeta*wkern;
    }

    // Add total contribution to acceleration for particle i
    for (k=0; k<ndim; k++) parti.a[k] += neiblist[j].m*draux[k]*paux;
    
    // If neighbour is also active, add contribution to force here
    if (neiblist[j].active)
      for (k=0; k<ndim; k++) neiblist[j].a[k] -= parti.m*draux[k]*paux;

  }
  // ==========================================================================

  return;
}



// ============================================================================
// GradhSph::ComputeGravForces
// Compute the contribution to the total gravitational force of particle 'i' 
// due to 'Nneib' neighbouring particles in the list 'neiblist'.
// ============================================================================
template <typename kernelclass>
void GradhSph<kernelclass>::ComputeDirectGravForces(int i, int Nneib, 
						    int *neiblist, 
						    SphParticle &parti,
						    SphParticle *sph)
{
  int j;
  int jj;
  int k;
  FLOAT dr[ndimmax];
  FLOAT drmag;
  FLOAT invdrmag;

  // Loop over all neighbouring particles in list
  // --------------------------------------------------------------------------
  for (jj=0; jj<Nneib; jj++) {
    j = neiblist[jj];

    for (k=0; k<ndim; k++) dr[k] = sph[j].r[k] - parti.r[k];
    drmag = DotProduct(dr,dr,ndim);
    drmag = sqrt(drmag);
    invdrmag = 1.0/(drmag + small_number);

    for (k=0; k<ndim; k++) parti.agrav[k] += sph[j].m*dr[k]*pow(invdrmag,3);
    parti.gpot -= sph[j].m*invdrmag;
  }

  return;
}



template class GradhSph<M4Kernel>;
template class GradhSph<QuinticKernel>;
