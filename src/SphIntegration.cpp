//=============================================================================
//  SphIntegration.cpp
//  Contains default functions for SphIntegration class.
//=============================================================================


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



//=============================================================================
//  SphIntegration::SphIntegration
/// SphIntegration constructor
//=============================================================================
template <int ndim>
SphIntegration<ndim>::SphIntegration
(DOUBLE accel_mult_aux,             ///< Copy of accel timestep multiplier
 DOUBLE courant_mult_aux):          ///< Copy of Courant timestep multipiler
  accel_mult(accel_mult_aux),
  courant_mult(courant_mult_aux)
{
}



//=============================================================================
//  SphIntegration::~SphIntegration
/// SphIntegration destructor
//=============================================================================
template <int ndim>
SphIntegration<ndim>::~SphIntegration()
{
}



//=============================================================================
// SphIntegration::Timestep
/// Default timestep size for SPH particles.  Takes the minimum of : 
/// (i)  const*h/(sound_speed + h*|div_v|)    (Courant condition)
/// (ii) const*sqrt(h/|a|)                    (Acceleration condition)
//=============================================================================
template <int ndim>
DOUBLE SphIntegration<ndim>::Timestep
(SphParticle<ndim> &part,               ///< Reference to SPH particle
 int hydro_forces)                      ///< Hydro forces flag
{
  DOUBLE timestep;                      // Minimum value of particle timesteps
  DOUBLE amag;                          // Magnitude of particle acceleration

  // Courant condition.  If hydro forces are not used, compute the 
  // timescale using only div_v, i.e. the compression timescale.
  if (hydro_forces == 1)
    //timestep = courant_mult*part.h/
      //(part.sound + part.h*fabs(part.div_v) +
       //0.6*(part.sound + 2.0*part.h*fabs(part.div_v))) ;
    timestep = courant_mult*part.h/
    (part.sound + part.h*fabs(part.div_v) + small_number_dp);
  else timestep = courant_mult*part.h/
    (part.h*fabs(part.div_v) + small_number_dp);

  // Acceleration condition
  amag = sqrt(DotProduct(part.a,part.a,ndim));
  timestep = min(timestep, accel_mult*sqrt(part.h/(amag + small_number_dp)));

  //cout << "TIMESTEP?? : " << part.h << "    " << part.sound << "     " << part.div_v << "    " << amag << "     " << timestep << endl;

  return timestep;
}




// Create instances of SphIntegration templates for all dimensions (1,2 and 3)
template class SphIntegration<1>;
template class SphIntegration<2>;
template class SphIntegration<3>;
