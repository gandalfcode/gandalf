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
#endif
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
// Sph::RandomBox
// ============================================================================
void GradhSph::RandomBox(void)
{
  printf("[GradhSph::RandomBox]");

  AllocateMemory(Nsph);

  for (int i=0; i<Nsph; i++) {
    for (int k=0; k<NDIM; k++) {
      sphdata[i].r[k] = (float)(rand()%RAND_MAX)/(float)RAND_MAX;
      sphdata[i].v[k] = (float)(rand()%RAND_MAX)/(float)RAND_MAX;
      sphdata[i].a[k] = 0.0f;
    }
    sphdata[i].m = 1.0f / (float) Nsph;
    sphdata[i].invomega = 0.5f;
  }

  return;
}



// ============================================================================
// GradhSph::ComputeH
// ============================================================================
void GradhSph::ComputeH(int i, int Nneib, int *neiblist, Parameters &simparams)
{
  int j;
  int jj;
  int k;
  int iteration = 0;
  float h_lower_bound = 0.0;
  float h_upper_bound = big_number;
  float hrangesqd;
  float dr[ndimmax];
  float drmag;
  float h_fac = simparams.floatparams["h_fac"];
  float h_converge = simparams.floatparams["h_converge"];

  while (1) {

    hrangesqd = kern->kernrangesqd*sphdata[i].h*sphdata[i].h;
    sphdata[i].invh = 1.0/sphdata[i].h;
    sphdata[i].hfactor = pow(sphdata[i].invh,ndim);
    sphdata[i].rho = 0.0;
    sphdata[i].div_v = 0.0;
    sphdata[i].invomega = 0.0;

    // Loop over all neighbours in list
    for (jj=0; jj<Nneib; jj++) {
      j = neiblist[jj];

      for (k=0; k<ndim; k++) dr[k] = sphdata[j].r[k] - sphdata[i].r[k];
      drmag = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
      
      // Skip particle if not a neighbour
      if (drmag > hrangesqd) continue;

      drmag = sqrt(drmag);
      for (k=0; k<ndim; k++) dr[k] /= (drmag + small_number);
      sphdata[i].rho += sphdata[j].m*sphdata[i].hfactor*
	kern->w0(drmag*sphdata[i].invh);
    }

    if (sphdata[i].rho > 0.0) sphdata[i].invrho = 1.0/sphdata[i].rho;

    cout << "h convergence : " << i << "   " << sphdata[i].h << "   " << 
      sphdata[i].rho << "   " << 
      fabs(sphdata[i].h - h_fac*pow(sphdata[i].m*sphdata[i].invrho,
				    1.0/(float)ndim)) << endl;

    // If h changes below some fixed tolerance, exit iteration loop
    if (sphdata[i].rho > 0.0 && sphdata[i].h > h_lower_bound && 
	fabs(sphdata[i].h - h_fac*pow(sphdata[i].m*sphdata[i].invrho,
				     1.0/(float)ndim)) < h_converge) break;

    // Use fixed-point iteration for now
    sphdata[i].h = h_fac*pow(sphdata[i].m*sphdata[i].invrho,1.0/(float)ndim);

  }

  sphdata[i].invomega = 1.0;

  return;
}



// ============================================================================
// GradhSph::ComputeSphProperties
// ============================================================================
void GradhSph::ComputeSphProperties(int i, int Nneib,int *neiblist, Parameters &simparams)
{
  sphdata[i].pfactor = eos->Pressure(sphdata[i])*
    sphdata[i].invrho*sphdata[i].invrho;
  sphdata[i].hfactor = pow(sphdata[i].invh,ndim+1);
  sphdata[i].sound = eos->SoundSpeed(sphdata[i]);

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
  float invrhomean;
  float wmean;
  float vsignal;

  sphdata[i].dudt = 0.0;
  hrangesqd = kern->kernrangesqd*sphdata[i].h*sphdata[i].h;


  // Loop over all potential neighbours in the list
  for (jj=0; jj<Nneib; jj++) {
    j = neiblist[jj];
    if (i == j) continue;

    // Calculate relative position vector and determine if particles
    // are neighbours or not. If not, skip to next potential neighbour.
    for (k=0; k<ndim; k++) dr[k] = sphdata[j].r[k] - sphdata[i].r[k];
    drmag = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
    if (drmag > hrangesqd && 
	drmag > kern->kernrangesqd*sphdata[j].h*sphdata[j].h) continue;

    // If particles are neighbours, continue computing hydro quantities
    drmag = sqrt(drmag);
    for (k=0; k<ndim; k++) dr[k] /= (drmag + small_number);
    for (k=0; k<ndim; k++) dv[k] = sphdata[j].v[k] - sphdata[i].v[k];
    dvdr = dv[0]*dr[0] + dv[1]*dr[1] * dv[2]*dr[2];

    // Compute hydro acceleration
    for (k=0; k<ndim; k++) sphdata[i].a[k] += sphdata[j].m*dr[k]*
      (sphdata[i].pfactor*sphdata[i].hfactor*kern->w1(drmag*sphdata[i].invh) +
       sphdata[j].pfactor*sphdata[j].hfactor*kern->w1(drmag*sphdata[j].invh));

    // Add dissipation terms
    if (dvdr < 0.) {
      wmean = 0.5*(kern->w1(drmag*sphdata[i].invh)*sphdata[i].hfactor +
          kern->w1(drmag*sphdata[j].invh)*sphdata[j].hfactor);
      invrhomean = 2. / (sphdata[i].rho+sphdata[j].rho);
      vsignal = sphdata[i].sound + sphdata[j].sound - beta_visc*dvdr;
      for (k=0; k<ndim; k++) sphdata[i].a[k] -= sphdata[j].m*alpha_visc*
          vsignal*dvdr*dr[k]*wmean*invrhomean;
    }

  }

  return;
}



// ============================================================================
// GradhSph::ComputeGravForces
// ============================================================================
void GradhSph::ComputeGravForce(int i, int j, float *agrav)
{
  return;
}
