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
int GradhSph::ComputeH(int i, int Nneib, int *neiblist, Parameters &params)
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

  // Local copies of necessary input parameters
  float h_fac = params.floatparams["h_fac"];
  float h_converge = params.floatparams["h_converge"];

  for (int k=0; k<ndimmax; k++) dr[k] = 0.0;
  for (int k=0; k<ndimmax; k++) dv[k] = 0.0;

  // Main smoothing length iteration loop
  // --------------------------------------------------------------------------
  do {

    // Initialise all variables for this value of h
    iteration++;
    hrangesqd = kern->kernrangesqd*sphdata[i].h*sphdata[i].h;
    sphdata[i].invh = 1.0/sphdata[i].h;
    sphdata[i].hfactor = pow(sphdata[i].invh,ndim);
    sphdata[i].rho = 0.0;
    sphdata[i].invomega = 0.0;
    sphdata[i].div_v = 0.0;

    // Loop over all neighbours in list
    // ------------------------------------------------------------------------
    for (jj=0; jj<Nneib; jj++) {
      j = neiblist[jj];

      for (k=0; k<ndim; k++) dr[k] = sphdata[j].r[k] - sphdata[i].r[k];
      drmag = DotProduct(dr,dr,ndim);
      
      // Skip particle if not a neighbour
      if (drmag > hrangesqd) continue;

      drmag = sqrt(drmag);
      for (k=0; k<ndim; k++) dr[k] /= (drmag + small_number);
      for (k=0; k<ndim; k++) dv[k] = sphdata[j].v[k] - sphdata[i].v[k];
      
      sphdata[i].div_v -= sphdata[j].m*DotProduct(dv,dr,ndim)*
	kern->w1(drmag*sphdata[i].invh)*sphdata[i].hfactor*sphdata[i].invh;
      sphdata[i].rho += sphdata[j].m*sphdata[i].hfactor*
	kern->w0(drmag*sphdata[i].invh);
      sphdata[i].invomega += sphdata[j].m*sphdata[i].hfactor*
	sphdata[i].invh*kern->womega(drmag*sphdata[i].invh);

    }
    // ------------------------------------------------------------------------
    
    if (sphdata[i].rho > 0.0) sphdata[i].invrho = 1.0/sphdata[i].rho;

    // If h changes below some fixed tolerance, exit iteration loop
    if (sphdata[i].rho > 0.0 && sphdata[i].h > h_lower_bound && 
	fabs(sphdata[i].h - h_fac*pow(sphdata[i].m*sphdata[i].invrho,
				      invndim)) < h_converge) break;

    // Use fixed-point iteration for now.  If this does not converge in a 
    // reasonable number of iterations (iteration_max), then assume something 
    // is wrong and switch to a bisection method, which should be guaranteed 
    // to converge, albeit slowly.
    // ------------------------------------------------------------------------
    if (iteration < iteration_max)
      sphdata[i].h = h_fac*pow(sphdata[i].m*sphdata[i].invrho,invndim);

    else if (iteration == iteration_max)
      sphdata[i].h = 0.5*(h_lower_bound + h_upper_bound);

    else if (iteration < 5*iteration_max) {
      if (sphdata[i].rho < small_number || 
	  sphdata[i].rho*pow(sphdata[i].h,ndim) > pow(h_fac,ndim)*sphdata[i].m)
	h_upper_bound = sphdata[i].h;
      else 
	h_lower_bound = sphdata[i].h;
      sphdata[i].h = 0.5*(h_lower_bound + h_upper_bound);
    }

    // If the smoothing length is too large for the neighbour list, exit 
    // routine and flag neighbour list error in order to generate a larger
    // neighbour list
    if (sphdata[i].h > h_max) return 0;
    

  } while (sphdata[i].h > h_lower_bound && sphdata[i].h < h_upper_bound);
  // --------------------------------------------------------------------------

  // Normalise all SPH sums correctly
  sphdata[i].invrho = 1.0/sphdata[i].rho;
  sphdata[i].h = h_fac*pow(sphdata[i].m*sphdata[i].invrho,invndim);
  sphdata[i].invh = 1.0/sphdata[i].h;
  sphdata[i].invomega = 1.0 + invndim*sphdata[i].h*
    sphdata[i].invomega*sphdata[i].invrho;
  sphdata[i].invomega = 1.0/sphdata[i].invomega;
  sphdata[i].div_v *= sphdata[i].invrho;

  return 1;
}



// ============================================================================
// GradhSph::ComputeSphProperties
// ============================================================================
void GradhSph::ComputeSphProperties(int i, int Nneib,int *neiblist, Parameters &simparams)
{
  sphdata[i].hfactor = pow(sphdata[i].invh,ndim+1);
  sphdata[i].u = eos->SpecificInternalEnergy(sphdata[i]);
  sphdata[i].sound = eos->SoundSpeed(sphdata[i]);
  sphdata[i].pfactor = eos->Pressure(sphdata[i])*
    sphdata[i].invrho*sphdata[i].invrho*sphdata[i].invomega;

  return;
}



// ============================================================================
// GradhSph::ComputeHydroForces
// ============================================================================
void GradhSph::ComputeHydroForces(int i, int Nneib, 
				  int *neiblist, Parameters &params)
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
  float wmean;
  float vsignal;
  int self_gravity = params.intparams["self_gravity"];
  string avisc = params.stringparams["avisc"];
  string acond = params.stringparams["acond"];

  sphdata[i].dudt = -eos->Pressure(sphdata[i])*sphdata[i].div_v*
    sphdata[i].invrho*sphdata[i].invomega;
  hrangesqd = kern->kernrangesqd*sphdata[i].h*sphdata[i].h;
  if (self_gravity == 1) sphdata[i].zeta = 0.0;

  for (int k=0; k<ndimmax; k++) dr[k] = 0.0;
  for (int k=0; k<ndimmax; k++) dv[k] = 0.0;

  // Loop over all potential neighbours in the list
  // --------------------------------------------------------------------------
  for (jj=0; jj<Nneib; jj++) {
    j = neiblist[jj];
    if (i == j) continue;

    // Calculate relative position vector and determine if particles
    // are neighbours or not. If not, skip to next potential neighbour.
    for (k=0; k<ndim; k++) dr[k] = sphdata[j].r[k] - sphdata[i].r[k];
    drmag = DotProduct(dr,dr,ndim);
    if (drmag > hrangesqd && 
	drmag > kern->kernrangesqd*sphdata[j].h*sphdata[j].h) continue;

    // If particles are neighbours, continue computing hydro quantities
    drmag = sqrt(drmag);
    for (k=0; k<ndim; k++) dr[k] /= (drmag + small_number);
    for (k=0; k<ndim; k++) dv[k] = sphdata[j].v[k] - sphdata[i].v[k];
    dvdr = DotProduct(dv,dr,ndim);

    // Compute hydro acceleration
    for (k=0; k<ndim; k++) sphdata[i].a[k] += sphdata[j].m*dr[k]*
      (sphdata[i].pfactor*sphdata[i].hfactor*kern->w1(drmag*sphdata[i].invh) +
       sphdata[j].pfactor*sphdata[j].hfactor*kern->w1(drmag*sphdata[j].invh));

    // Add dissipation terms (for approaching particle-pairs)
    // ------------------------------------------------------------------------
    if (dvdr < 0.0) {
      wmean = 0.5*(kern->w1(drmag*sphdata[i].invh)*sphdata[i].hfactor +
          kern->w1(drmag*sphdata[j].invh)*sphdata[j].hfactor);
      invrhomean = 2.0 / (sphdata[i].rho + sphdata[j].rho);

      // Artificial viscosity term
      if (avisc == "mon97") {
	vsignal = sphdata[i].sound + sphdata[j].sound - beta_visc*dvdr;
	for (k=0; k<ndim; k++) sphdata[i].a[k] -= 
	  sphdata[j].m*alpha_visc*vsignal*dvdr*dr[k]*wmean*invrhomean;
	sphdata[i].dudt -= 0.5*sphdata[j].m*alpha_visc*
	  vsignal*wmean*invrhomean*dvdr*dvdr;
      }

      // Artificial conductivity term
      if (acond == "price2008") {
	vsignal = sqrt(fabs(eos->Pressure(sphdata[i]) - 
			    eos->Pressure(sphdata[j]))*invrhomean);
	sphdata[i].dudt += sphdata[j].m*vsignal*
	  (sphdata[i].u - sphdata[j].u)*wmean*invrhomean;
      }
    }
    // ------------------------------------------------------------------------

    // Calculate zeta term for mean-h gravity
    if (self_gravity == 1) {
      invhmean = 2.0/(sphdata[i].h + sphdata[j].h);
      sphdata[i].zeta += sphdata[j].m*invhmean*invhmean*
	kern->wzeta(drmag*invhmean);
    }

  }
  // --------------------------------------------------------------------------

  // Normalise zeta term
  sphdata[i].zeta *= -invndim*sphdata[i].h*
    sphdata[i].invrho*sphdata[i].invomega;

  return;
}



// ============================================================================
// GradhSph::ComputeGravForces
// Compute the contribution to the total gravitational force of particle 'i' 
// due to 'Nneib' neighbouring particles in the list 'neiblist'.
// ============================================================================
void GradhSph::ComputeGravForces(int i, int Nneib, int *neiblist)
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
    j = neiblist[jj];
    if (i == j) continue;

    for (k=0; k<ndim; k++) dr[k] = sphdata[j].r[k] - sphdata[i].r[k];
    drmag = DotProduct(dr,dr,ndim);
    drmag = sqrt(drmag);
    invdrmag = 1.0/(drmag + small_number);

    // Calculate kernel-softened gravity if within mean-h kernel.
    // Otherwise, use point-mass Newton's law of gravitation.
    if (2.0*drmag < kernrange*(sphdata[i].h + sphdata[j].h)) {
      invhmean = 2.0/(sphdata[i].h + sphdata[j].h);
      for (k=0; k<ndim; k++) sphdata[i].agrav[k] += sphdata[j].m*dr[k]*
	invdrmag*invhmean*invhmean*kern->wgrav(drmag*invhmean);
      sphdata[i].gpot -= sphdata[j].m*invhmean*kern->wpot(drmag*invhmean);
    }
    else {
      for (k=0; k<ndim; k++) sphdata[i].agrav[k] += 
	sphdata[j].m*dr[k]*pow(invdrmag,3);
      sphdata[i].gpot -= sphdata[j].m*invdrmag;
    }

    if (drmag*sphdata[i].invh < kernrange) {
      for (k=0; k<ndim; k++) 
	sphdata[i].agrav[k] += 0.5*sphdata[j].m*dr[k]*invdrmag*
	  sphdata[i].zeta*kern->w1(drmag*sphdata[i].invh)*sphdata[i].hfactor;
    }

    if (drmag*sphdata[j].invh < kernrange) {
      for (k=0; k<ndim; k++) 
	sphdata[i].agrav[k] += 0.5*sphdata[j].m*dr[k]*invdrmag*
	  sphdata[j].zeta*kern->w1(drmag*sphdata[j].invh)*sphdata[j].hfactor;
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
