// ============================================================================
// SphLeapfrogKDK.cpp
// ============================================================================


#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <math.h>
#include "Dimensions.h"
#include "Sph.h"
#include "SphKernel.h"
#include "SphIntegration.h"
#include "SphParticle.h"
#include "EOS.h"
#include "Debug.h"
using namespace std;



// ============================================================================
// SphLeapfrogKDK::SphLeapfrogKDK()
// ============================================================================
SphLeapfrogKDK::SphLeapfrogKDK(int ndimaux, int vdimaux, 
			       DOUBLE accel_mult_aux, DOUBLE courant_mult_aux) :
  SphIntegration(ndimaux, vdimaux, accel_mult_aux, courant_mult_aux)
{
}



// ============================================================================
// SphLeapfrogKDK::~SphLeapfrog()
// ============================================================================
SphLeapfrogKDK::~SphLeapfrogKDK()
{
}



// ============================================================================
// SphLeapfrogKDK::AdvanceParticles
// ============================================================================
void SphLeapfrogKDK::AdvanceParticles(int n, int level_step, 
				      int Nsph, SphParticle *sph, FLOAT dt)
{
  int i;
  int k;
  int nstep;

  debug2("[SphLeapfrogKDK::AdvanceParticles]");

  for (i=0; i<Nsph; i++) {
    nstep = pow(2,level_step - sph[i].level);
    for (k=0; k<ndim; k++) sph[i].r[k] = sph[i].r0[k] + sph[i].v0[k]*dt
      + (FLOAT) 0.5*sph[i].a0[k]*dt*dt;
    for (k=0; k<vdim; k++) sph[i].v[k] = sph[i].v0[k] + sph[i].a0[k]*dt;
    if (n%nstep == 0) sph[i].active = true;
  }

  return;
}
 


// ============================================================================
// SphLeapfrogKDK::CorrectionTerms
// ============================================================================
void SphLeapfrogKDK::CorrectionTerms(int n, int level_step, 
				     int Nsph, SphParticle *sph, FLOAT dt)
{
  int i;
  int k;
  int nstep;

  debug2("[SphLeapfrogKDK::CorrectionTerms]");

  for (i=0; i<Nsph; i++) {
    nstep = pow(2,level_step - sph[i].level);
    if (n%nstep == 0)
      for (k=0; k<ndim; k++) sph[i].v[k] += 
	(FLOAT) 0.5*(sph[i].a[k] - sph[i].a0[k])*dt*(FLOAT) nstep;
  }

  return;
}



// ============================================================================
// SphLeapfrogKDK::EndTimestep
// ============================================================================
void SphLeapfrogKDK::EndTimestep(int n, int level_step, 
				 int Nsph, SphParticle *sph)
{
  int i;
  int k;
  int nstep;

  debug2("[SphLeapfrogKDK::EndTimestep]");

  for (i=0; i<Nsph; i++) {
    nstep = pow(2,level_step - sph[i].level);
    if (n%nstep == 0) {
      for (k=0; k<ndim; k++) sph[i].r0[k] = sph[i].r[k];
      for (k=0; k<ndim; k++) sph[i].v0[k] = sph[i].v[k];
      for (k=0; k<ndim; k++) sph[i].a0[k] = sph[i].a[k];
      sph[i].active = false;
    }
  }

  return;
}
