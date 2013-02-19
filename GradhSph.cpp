// ============================================================================
// GradhSph.cpp
// ============================================================================


#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <math.h>
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
GradhSph::GradhSph(int ndimaux, int vdimaux, int bdimaux)
{
#if !defined(FIXED_DIMENSIONS)
  ndim = ndimaux;
  vdim = vdimaux;
  bdim = bdimaux;
  invndim = 1.0/(float)ndim;
#endif
  allocated = false;
  Nsph = 0;
  Nsphmax = 0;
}



// ============================================================================
// GradhSph::~GradhSph
// ============================================================================
GradhSph::~GradhSph()
{
}



// ============================================================================
// GradhSph::ComputeH
// Compute the value of the smoothing length of particle 'i' by iterating  
// the relation :  h = h_fac*(m/rho)^(1/ndim).
// Uses the previous value of h as a starting guess and then uses either 
// a Newton-Rhapson solver, or fixed-point iteration, to converge on the 
// correct value of h.  The maximum tolerance used for deciding whether the 
// iteration has converged is given by the 'h_converge' parameter.
// ============================================================================
int GradhSph::ComputeH(int i, SphParticle &parti, 
		       int Nneib, SphParticle *neiblistpart)
{
  int j;                                      // Neighbour id
  int jj;                                     // Aux. neighbour counter
  int k;                                      // Dimension counter
  int iteration = 0;                          // h-rho iteration counter
  int iteration_max = 30;                     // Max. no of iterations
  float h_lower_bound = 0.0;                  // Lower bound on h
  float h_upper_bound = big_number;           // Upper bound on h
  float h_max = big_number;                   // Max. allowed value of h
  float hrangesqd;                            // Kernel extent (squared)
  float dr[ndimmax];                          // Relative position vector
  float dv[ndimmax];                          // Relative velocity vector
  float drmag;                                // Neighbour distance
  float invdrmag;                             // ..

  // Particle local copy
  //SphParticle parti = sphdata[i];

  // Main smoothing length iteration loop
  // --------------------------------------------------------------------------
  do {

    // Initialise all variables for this value of h
    iteration++;
    hrangesqd = kern->kernrangesqd*parti.h*parti.h;
    parti.invh = 1.0/parti.h;
    parti.hfactor = pow(parti.invh,ndim);
    parti.rho = 0.0;
    parti.invomega = 0.0;
    parti.div_v = 0.0;

    // Loop over all neighbours in list
    // ------------------------------------------------------------------------
    for (jj=0; jj<Nneib; jj++) {
//      j = neiblist[jj];

      for (k=0; k<ndim; k++) dr[k] = neiblistpart[jj].r[k] - parti.r[k];
      drmag = DotProduct(dr,dr,ndim);
      
      // Skip particle if not a neighbour
      if (drmag > hrangesqd) continue;

      drmag = sqrt(drmag);
      invdrmag = 1.0/(drmag + small_number);
      for (k=0; k<ndim; k++) dr[k] *= invdrmag;
      for (k=0; k<ndim; k++) dv[k] = neiblistpart[jj].v[k] - parti.v[k];
      
      parti.div_v -= neiblistpart[jj].m*DotProduct(dv,dr,ndim)*
	kern->w1(drmag*parti.invh)*parti.hfactor*parti.invh;
      parti.rho += neiblistpart[jj].m*parti.hfactor*
	kern->w0(drmag*parti.invh);
      parti.invomega += neiblistpart[jj].m*parti.hfactor*
	parti.invh*kern->womega(drmag*parti.invh);

    }
    // ------------------------------------------------------------------------

    if (parti.rho > 0.0) parti.invrho = 1.0/parti.rho;

    // If h changes below some fixed tolerance, exit iteration loop
    if (parti.rho > 0.0 && parti.h > h_lower_bound &&
	fabs(parti.h - h_fac*pow(parti.m*parti.invrho,
				      invndim)) < h_converge) break;

    // Use fixed-point iteration for now.  If this does not converge in a 
    // reasonable number of iterations (iteration_max), then assume something 
    // is wrong and switch to a bisection method, which should be guaranteed 
    // to converge, albeit slowly.
    // ------------------------------------------------------------------------
    if (iteration < iteration_max)
      parti.h = h_fac*pow(parti.m*parti.invrho,invndim);

    else if (iteration == iteration_max)
      parti.h = 0.5*(h_lower_bound + h_upper_bound);

    else if (iteration < 5*iteration_max) {
      if (parti.rho < small_number ||
	  parti.rho*pow(parti.h,ndim) > pow(h_fac,ndim)*parti.m)
	h_upper_bound = parti.h;
      else 
	h_lower_bound = parti.h;
      parti.h = 0.5*(h_lower_bound + h_upper_bound);
    }

    // If the smoothing length is too large for the neighbour list, exit 
    // routine and flag neighbour list error in order to generate a larger
    // neighbour list
    if (parti.h > h_max) return 0;
    

  } while (parti.h > h_lower_bound && parti.h < h_upper_bound);
  // --------------------------------------------------------------------------

  // Normalise all SPH sums correctly
  parti.invrho = 1.0/parti.rho;
  parti.h = h_fac*pow(parti.m*parti.invrho,invndim);
  parti.invh = 1.0/parti.h;
  parti.invomega = 1.0 + invndim*parti.h*parti.invomega*parti.invrho;
  parti.invomega = 1.0/parti.invomega;
  parti.div_v *= parti.invrho;
  parti.hfactor = pow(parti.invh,ndim+1);
  parti.u = eos->SpecificInternalEnergy(parti);
  parti.sound = eos->SoundSpeed(parti);
  parti.press = eos->Pressure(parti);
  parti.pfactor = eos->Pressure(parti)*
    parti.invrho*parti.invrho*parti.invomega;

  // Copies back particle data
  //sphdata[i] = parti;

  return 1;
}



// ============================================================================
// GradhSph::ComputeHydroForces
// ============================================================================
void GradhSph::ComputeHydroForces(int i, SphParticle &parti,
				  int Nneib, SphParticle *neiblist)
{
  int j;
  int jj;
  int k;
  float hrangesqd;
  float dr[ndimmax];
  float dv[ndimmax];
  float dvdr;
  float drmag;
  float invhmean;
  float invrhomean;
  float invdrmag;
  float wmean;
  float vsignal;

  parti.dudt = -parti.press*parti.div_v*
    parti.invrho*parti.invomega;
  hrangesqd = kern->kernrangesqd*parti.h*parti.h;
  if (self_gravity == 1) parti.zeta = 0.0;

  // Loop over all potential neighbours in the list
  // --------------------------------------------------------------------------
  for (jj=0; jj<Nneib; jj++) {
//    j = neiblist[jj];
//    if (i == j) continue;

    // Calculate relative position vector and determine if particles
    // are neighbours or not. If not, skip to next potential neighbour.
    for (k=0; k<ndim; k++) dr[k] = neiblist[jj].r[k] - parti.r[k];
    drmag = DotProduct(dr,dr,ndim);
    if (drmag > hrangesqd && 
	drmag > kern->kernrangesqd*neiblist[jj].h*neiblist[jj].h) continue;

    // If particles are neighbours, continue computing hydro quantities
    drmag = sqrt(drmag);
    invdrmag = 1.0/(drmag + small_number);
    for (k=0; k<ndim; k++) dr[k] *= invdrmag;
    for (k=0; k<ndim; k++) dv[k] = neiblist[jj].v[k] - parti.v[k];
    dvdr = DotProduct(dv,dr,ndim);

    // Compute hydro acceleration
    for (k=0; k<ndim; k++) parti.a[k] += neiblist[jj].m*dr[k]*
      (parti.pfactor*parti.hfactor*kern->w1(drmag*parti.invh) +
       neiblist[jj].pfactor*neiblist[jj].hfactor*kern->w1(drmag*neiblist[jj].invh));

    // Add dissipation terms (for approaching particle-pairs)
    // ------------------------------------------------------------------------
    if (dvdr < 0.0) {
      wmean = 0.5*(kern->w1(drmag*parti.invh)*parti.hfactor +
          kern->w1(drmag*neiblist[jj].invh)*neiblist[jj].hfactor);
      invrhomean = 2.0 / (parti.rho + neiblist[jj].rho);

      // Artificial viscosity term
      if (avisc == "mon97") {
	vsignal = parti.sound + neiblist[jj].sound - beta_visc*dvdr;
	for (k=0; k<ndim; k++) parti.a[k] -= 
	  neiblist[jj].m*alpha_visc*vsignal*dvdr*dr[k]*wmean*invrhomean;
	parti.dudt -= 0.5*neiblist[jj].m*alpha_visc*
	  vsignal*wmean*invrhomean*dvdr*dvdr;
      }

      // Artificial conductivity term
      if (acond == "price2008") {
	vsignal = sqrt(fabs(parti.press - neiblist[jj].press)*invrhomean);
	parti.dudt += neiblist[jj].m*vsignal*
	  (parti.u - neiblist[jj].u)*wmean*invrhomean;
      }
    }
    // ------------------------------------------------------------------------

    // Calculate zeta term for mean-h gravity
    if (self_gravity == 1) {
      invhmean = 2.0/(parti.h + neiblist[jj].h);
      parti.zeta += neiblist[jj].m*invhmean*invhmean*
	kern->wzeta(drmag*invhmean);
    }

  }
  // --------------------------------------------------------------------------

  // Normalise zeta term
  parti.zeta *= -invndim*parti.h*
    parti.invrho*parti.invomega;

  return;
}



// ============================================================================
// GradhSph::ComputeGravForces
// Compute the contribution to the total gravitational force of particle 'i' 
// due to 'Nneib' neighbouring particles in the list 'neiblist'.
// ============================================================================
void GradhSph::ComputeGravForces(int i, int Nneib, SphParticle *neiblist)
{
  int j;
  int jj;
  int k;
  float hrangesqd;
  float dr[ndimmax];
  float drmag;
  float invdrmag;
  float invhmean;
  float kernrange = kern->kernrange;

  // Loop over all neighbouring particles in list
  // --------------------------------------------------------------------------
  for (jj=0; jj<Nneib; jj++) {
//    j = neiblist[jj];
    if (i == j) continue;

    for (k=0; k<ndim; k++) dr[k] = neiblist[jj].r[k] - sphdata[i].r[k];
    drmag = DotProduct(dr,dr,ndim);
    drmag = sqrt(drmag);
    invdrmag = 1.0/(drmag + small_number);

    // Calculate kernel-softened gravity if within mean-h kernel.
    // Otherwise, use point-mass Newton's law of gravitation.
    if (2.0*drmag < kernrange*(sphdata[i].h + neiblist[jj].h)) {
      invhmean = 2.0/(sphdata[i].h + neiblist[jj].h);
      for (k=0; k<ndim; k++) sphdata[i].agrav[k] += neiblist[jj].m*dr[k]*
	invdrmag*invhmean*invhmean*kern->wgrav(drmag*invhmean);
      sphdata[i].gpot -= neiblist[jj].m*invhmean*kern->wpot(drmag*invhmean);
    }
    else {
      for (k=0; k<ndim; k++) sphdata[i].agrav[k] += 
	neiblist[jj].m*dr[k]*pow(invdrmag,3);
      sphdata[i].gpot -= neiblist[jj].m*invdrmag;
    }

    if (drmag*sphdata[i].invh < kernrange) {
      for (k=0; k<ndim; k++) 
	sphdata[i].agrav[k] += 0.5*neiblist[jj].m*dr[k]*invdrmag*
	  sphdata[i].zeta*kern->w1(drmag*sphdata[i].invh)*sphdata[i].hfactor;
    }

    if (drmag*neiblist[jj].invh < kernrange) {
      for (k=0; k<ndim; k++) 
	sphdata[i].agrav[k] += 0.5*neiblist[jj].m*dr[k]*invdrmag*
	  neiblist[jj].zeta*kern->w1(drmag*neiblist[jj].invh)*neiblist[jj].hfactor;
    }

  }

  return;
}



// ============================================================================
// GradhSph::ComputeMeanhZeta
// ============================================================================
void GradhSph::ComputeMeanhZeta(int i, int Nneib, int *neiblist)
{
  int j;
  int jj;
  int k;
  float dr[ndimmax];
  float drmag;
  float invhmean;
  float kernrangesqd = kern->kernrangesqd;

  sphdata[i].zeta = 0.0;

  // Loop over all potential neighbours in the list
  // --------------------------------------------------------------------------
  for (jj=0; jj<Nneib; jj++) {
    j = neiblist[jj];
    if (i == j) continue;

    // Calculate relative position vector and determine if particles
    // are neighbours or not. If not, skip to next potential neighbour.
    for (k=0; k<ndim; k++) dr[k] = sphdata[j].r[k] - sphdata[i].r[k];
    drmag = DotProduct(dr,dr,ndim);
    if (4.0*drmag > kernrangesqd*pow(sphdata[i].h + sphdata[j].h,2)) continue;

    // If particles are neighbours, continue computing hydro quantities
    drmag = sqrt(drmag);
    for (k=0; k<ndim; k++) dr[k] /= (drmag + small_number);
    invhmean = 2.0/(sphdata[i].h + sphdata[j].h);
    sphdata[i].zeta += sphdata[j].m*invhmean*invhmean*
      kern->wzeta(drmag*invhmean);

  }
  // --------------------------------------------------------------------------

  // Normalise zeta term
  sphdata[i].zeta *= -invndim*sphdata[i].h*
    sphdata[i].invrho*sphdata[i].invomega;

  return;
}
