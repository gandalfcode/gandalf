//=============================================================================
//  SphIntegration.cpp
//  Contains default functions for SphIntegration class.
//
//  This file is part of GANDALF :
//  Graphical Astrophysics code for N-body Dynamics and Lagrangian Fluids
//  https://github.com/gandalfcode/gandalf
//  Contact : gandalfcode@gmail.com
//
//  Copyright (C) 2013  D. A. Hubber, G Rosotti
//
//  GANDALF is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 2 of the License, or
//  (at your option) any later version.
//
//  GANDALF is distributed in the hope that it will be useful, but
//  WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  General Public License (http://www.gnu.org/licenses) for more details.
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
    //imestep = courant_mult*part.h/
      //(part.sound + part.h*fabs(part.div_v) +
       //0.6*(part.sound + 2.0*part.h*fabs(part.div_v))) ;
    timestep = courant_mult*part.h/
      (part.sound + part.h*fabs(part.div_v) + small_number_dp);
  else timestep = courant_mult*part.h/
    (part.h*fabs(part.div_v) + small_number_dp);

  // Acceleration condition
  amag = sqrt(DotProduct(part.a,part.a,ndim));
  timestep = min(timestep, accel_mult*sqrt(part.h/(amag + small_number_dp)));

  return timestep;
}




// Create instances of SphIntegration templates for all dimensions (1,2 and 3)
template class SphIntegration<1>;
template class SphIntegration<2>;
template class SphIntegration<3>;
