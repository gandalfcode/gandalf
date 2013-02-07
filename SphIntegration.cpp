// ============================================================================
// SphIntegration.cpp
// ============================================================================

#include <algorithm>
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
#include "InlineFuncs.h"
#include "Debug.h"
using namespace std;



// ============================================================================
// SphIntegration::SphIntegration
// ============================================================================
SphIntegration::SphIntegration(int ndimaux, int vdimaux, double accel_mult_aux, double courant_mult_aux):
#if !defined(FIXED_DIMENSIONS)
  ndim(ndimaux),
  vdim(vdimaux),
#endif
  accel_mult(accel_mult_aux),
  courant_mult(courant_mult_aux)
{
}



// ============================================================================
// SphIntegration::~SphIntegration
// ============================================================================
SphIntegration::~SphIntegration()
{
}



// ============================================================================
// ..
// ============================================================================
double SphIntegration::Timestep(SphParticle &part, int hydro_forces)
{
  double timestep;
  double amag;

  //Courant condition
  //if (params.intparams["hydro_forces"] == 1)
  if (hydro_forces == 1)
    timestep = courant_mult*part.h/
      (part.sound + part.h*fabs(part.div_v) + small_number_dp);
  else timestep = courant_mult*part.h/
    (part.h*fabs(part.div_v) + small_number_dp);

  //Acceleration condition
  amag = sqrt(DotProduct(part.a,part.a,ndim));
  timestep = min(timestep, accel_mult*sqrt(part.h/(amag + small_number_dp)));

  return timestep;
}
