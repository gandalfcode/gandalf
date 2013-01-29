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
SphIntegration::SphIntegration(double accel_mult_aux, double courant_mult_aux):
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
double SphIntegration::Timestep(SphParticle &part)
{
  double timestep;
  double amag;

  //Courant condition
  timestep = courant_mult*part.h/
    (part.sound + part.h*fabs(part.div_v) + small_number_dp);

  //Acceleration condition
  amag = sqrt(DotProduct(part.a,part.a));
  timestep = min(timestep, accel_mult*sqrt(part.h/(amag + small_number_dp)));


  //TODO: implement energy condition (once we will be solving the energy equation)

  return timestep;
}
