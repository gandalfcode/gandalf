//=================================================================================================
//  SphIntegration.cpp
//  Contains default functions for SphIntegration class.
//
//  This file is part of GANDALF :
//  Graphical Astrophysics code for N-body Dynamics And Lagrangian Fluids
//  https://github.com/gandalfcode/gandalf
//  Contact : gandalfcode@gmail.com
//
//  Copyright (C) 2013  D. A. Hubber, G. Rosotti
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
//=================================================================================================


#include <algorithm>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <math.h>
#include <iostream>

#include "../Headers/Integration.h"
#include "Sph.h"
#include "SmoothingKernel.h"
#include "Particle.h"
#include "EOS.h"
#include "InlineFuncs.h"
#include "Debug.h"
using namespace std;



//=================================================================================================
//  SphIntegration::SphIntegration
/// SphIntegration constructor
//=================================================================================================
template <int ndim>
SphIntegration<ndim>::SphIntegration
(DOUBLE accel_mult_aux,             ///< Copy of accel timestep multiplier
 DOUBLE courant_mult_aux,           ///< Copy of Courant timestep multipiler
 DOUBLE energy_mult_aux,            ///< Copy of Energy timestep multipiler
 eosenum gas_eos_aux,
 tdaviscenum tdavisc_aux) :
  accel_mult(accel_mult_aux),
  courant_mult(courant_mult_aux),
  energy_mult(energy_mult_aux),
  gas_eos(gas_eos_aux),
  tdavisc(tdavisc_aux)
{
}



//=================================================================================================
//  SphIntegration::~SphIntegration
/// SphIntegration destructor
//=================================================================================================
template <int ndim>
SphIntegration<ndim>::~SphIntegration()
{
}



//=================================================================================================
//  SphIntegration::Timestep
/// Default timestep size for SPH particles.  Takes the minimum of :
/// (i)  const*h/(sound_speed + h*|div_v|)    (Courant condition)
/// (ii) const*sqrt(h/|a|)                    (Acceleration condition)
//=================================================================================================
template <int ndim>
DOUBLE SphIntegration<ndim>::Timestep
 (Particle<ndim> &partaux,                         ///< [inout] Reference to SPH particle
  Hydrodynamics<ndim> *hydro)                      ///< [in] Pointer to main SPH object
{
  DOUBLE timestep;                     // Minimum value of particle timesteps
  //DOUBLE adotmag;                      // Magnitude of particle jerk
  DOUBLE amag;                         // Magnitude of particle acceleration


  Sph<ndim>* sph = reinterpret_cast<Sph<ndim>*>(hydro) ;
  SphParticle<ndim>& part= *reinterpret_cast<SphParticle<ndim>*>(&partaux);

  // Courant condition.  If hydro forces are not used, compute the
  // timescale using only div_v, i.e. the compression timescale.
  if (sph->hydro_forces == 1 && sph->avisc == mon97 && part.flags.check(inside_sink)) {
    timestep = courant_mult*part.h/(part.sound + part.h*fabs(part.div_v) + small_number_dp);
  }
  else if (sph->hydro_forces == 1 && sph->avisc == mon97) {
    //timestep = courant_mult*part.h/
      //(part.sound + part.h*fabs(part.div_v) + 0.6*(part.sound + 2.0*part.h*fabs(part.div_v)));
    timestep = courant_mult*part.h/(part.sound + part.h*fabs(part.div_v) + small_number_dp);
  }
  else if (sph->hydro_forces == 1) {
    timestep = courant_mult*part.h/(part.sound + part.h*fabs(part.div_v) + small_number_dp);
  }
  else {
    timestep = courant_mult*part.h/(part.h*fabs(part.div_v) + small_number_dp);
  }


  // Acceleration condition
  amag = sqrt(DotProduct(part.a, part.a, ndim));
  timestep = min(timestep, accel_mult*sqrt(part.h/(amag + small_number_dp)));


  // Explicit energy integration timestep condition
  if (gas_eos == energy_eqn && part.ptype == gas_type) {
    timestep = min(timestep, this->energy_mult*(DOUBLE) (part.u/(fabs(part.dudt) + small_number)));
  }


  if (timestep > 1.0e20 || timestep <= 0.0) {
    cout << "Timestep problem : " << timestep << "   " << amag << "   " << part.h << "    "
         << part.sound << "    " << part.h*fabs(part.div_v) << endl;
    cout << "tcourant : " << courant_mult*part.h/(part.sound + part.h*fabs(part.div_v) + small_number_dp)
         << "    taccel : " << accel_mult*sqrt(part.h/(amag + small_number_dp))
         << "    tenergy : " << energy_mult*(DOUBLE) (part.u/(fabs(part.dudt) + small_number)) << endl;
    ExceptionHandler::getIstance().raise("Error : Hydro timestep too large");
  }

  assert(isfinite(amag) && isfinite(timestep));

  return timestep;
}


// Create instances of SphIntegration templates for all dimensions (1,2 and 3)
template class SphIntegration<1>;
template class SphIntegration<2>;
template class SphIntegration<3>;
