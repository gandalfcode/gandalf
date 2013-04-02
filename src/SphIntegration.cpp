// ============================================================================
// SphIntegration.cpp
// Contains default functions for SphIntegration class.
// ============================================================================


#include <algorithm>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <math.h>
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
template <int ndim>
SphIntegration<ndim>::SphIntegration(DOUBLE accel_mult_aux, DOUBLE courant_mult_aux):
//#if !defined(FIXED_DIMENSIONS)
//  ndim(ndimaux),
//  vdim(ndim),
//#endif
  accel_mult(accel_mult_aux),
  courant_mult(courant_mult_aux)
{
}



// ============================================================================
// SphIntegration::~SphIntegration
// ============================================================================
template <int ndim>
SphIntegration<ndim>::~SphIntegration()
{
}



// ============================================================================
// SphIntegration::Timestep
// Default timestep size for SPH particles.  Takes the minimum of : 
// (i)  const*h/(sound_speed + h*|div_v|)    (Courant condition)
// (ii) const*sqrt(h/|a|)                    (Acceleration condition)
// ============================================================================
template <int ndim>
DOUBLE SphIntegration<ndim>::Timestep(SphParticle<ndim> &part, int hydro_forces)
{
  DOUBLE timestep;
  DOUBLE amag;

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

template class SphIntegration<1>;
template class SphIntegration<2>;
template class SphIntegration<3>;
