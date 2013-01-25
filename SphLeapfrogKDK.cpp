// ============================================================================
// SphLFKDK.cpp
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
// ..
// ============================================================================
SphLFKDK::SphLFKDK()
{
}



// ============================================================================
// ..
// ============================================================================
SphLFKDK::~SphLFKDK()
{
}





// ============================================================================
// ..
// ============================================================================
void SphLFKDK::AdvanceParticles(int Nsph, SphParticle *sph, double dt)
{
  int i,k;
  for (i=0; i<Nsph; i++) {
    for (k=0; k<ndim; k++) sph[i].r[k] = sph[i].r0[k] + sph[i].v[k]*dt;
    for (k=0; k<ndim; k++) sph[i].v[k] = sph[i].v0[k] + sph[i].a[k]*dt;
  }
}



// ============================================================================
// ..
// ============================================================================
void SphLFKDK::CorrectionTerms(int Nsph, SphParticle *sph, double dt)
{
  return;
}




// ============================================================================
// ..
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


