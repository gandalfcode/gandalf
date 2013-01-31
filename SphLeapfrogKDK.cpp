// ============================================================================
// SphLeapfrogKDK.cpp
// ============================================================================


#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <iostream>
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
SphLFKDK::SphLFKDK(double accel_mult_aux, double courant_mult_aux) :
  SphIntegration(accel_mult_aux, courant_mult_aux)
{
}



// ============================================================================
// SphLeapfrogKDK::~SphLeapfrog()
// ============================================================================
SphLFKDK::~SphLFKDK()
{
}



// ============================================================================
// SphLeapfrogKDK::AdvanceParticles
// ============================================================================
void SphLFKDK::AdvanceParticles(int Nsph, SphParticle *sph, double dt)
{
  int i;
  int k;

  for (i=0; i<Nsph; i++) {
    for (k=0; k<ndim; k++) sph[i].r[k] = sph[i].r0[k] + 
      sph[i].v[k]*dt + 0.5*sph[i].a[k]*dt*dt;
    for (k=0; k<ndim; k++) sph[i].v[k] = sph[i].v0[k] + sph[i].a[k]*dt;
  }

  return;
}
 


// ============================================================================
// SphLeapfrogKDK::CorrectionTerms
// ============================================================================
void SphLFKDK::CorrectionTerms(int Nsph, SphParticle *sph, double dt)
{
  int i;
  int k;

  for (i=0; i<Nsph; i++) {
    for (k=0; k<ndim; k++) sph[i].v[k] += 
      0.5*(sph[i].a[k] - sph[i].a0[k])*dt;
  }

  return;
}




// ============================================================================
// SphLeapfrogKDK::EndTimestep
// ============================================================================
void SphLFKDK::EndTimestep(int n, int Nsph, SphParticle *sph)
{
  int i,k;

  debug2("[SphLFKDK::EndTimestep]\n");

  for (i=0; i<Nsph; i++) {
    for (k=0; k<ndim; k++) sph[i].r0[k] = sph[i].r[k];
    for (k=0; k<ndim; k++) sph[i].v0[k] = sph[i].v[k];
    for (k=0; k<ndim; k++) sph[i].a0[k] = sph[i].a[k];
  }

  return;
}


