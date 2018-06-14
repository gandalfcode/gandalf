//=================================================================================================
//  RadiativeFB.cpp
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

#include <iostream>
#include "RadiativeFB.h"
#include "Constants.h"
#include "Debug.h"
#include "Hydrodynamics.h"
#include "Parameters.h"
#include "Particle.h"
#include "SimUnits.h"
#include "Nbody.h"
using namespace std;

//=================================================================================================
//  RadiativeFB::RadiativeFB()
/// RadiativeFB class constructor.
//=================================================================================================
template <int ndim>
RadiativeFB<ndim>::RadiativeFB
(SimUnits *simunits,
 Parameters *params)
{
  // Ensure 2 or 3 dimensions are used
  if (ndim < 2) {
    ExceptionHandler::getIstance().raise("Radiative feedback requires at least 2 dimensions!");
  }

  // Get and limit number of central objects
  int Ncentral = params->intparams["disc_heating"];
  if (Ncentral < 0) Ncentral = 0;
  if (Ncentral > 2) Ncentral = 2;

  // Convert to units first
  temp_unit = simunits->temp.outscale * simunits->temp.outSI;
  temp_inf = params->floatparams["temp_ambient"] / temp_unit;
  temp_inf4 = pow(temp_inf, 4.0);

  ambient_heating = params->intparams["ambient_heating"];

  if (params->intparams["disc_heating"]) {
    disc_heating = new DiscHeating<ndim>(simunits, params, Ncentral);
  }
  else {
    disc_heating = NULL;
  }

  if (params->intparams["sink_heating"]) {
    sink_heating = new SinkHeating<ndim>(simunits, params, Ncentral);
  }
  else {
    sink_heating = NULL;
  }

  return;
}

//=================================================================================================
//  RadiativeFB::~RadiativeFB()
/// RadiativeFB class destructor.
//=================================================================================================
template <int ndim>
RadiativeFB<ndim>::~RadiativeFB()
{
  if (disc_heating != NULL) delete disc_heating;
  if (sink_heating != NULL) delete sink_heating;
}

//=================================================================================================
//  RadiativeFB::AmbientTemp()
/// Combines heating from an external ambient source, a disc enforced temperature profile, and
/// heating from sinks.
//=================================================================================================
template <int ndim>
FLOAT RadiativeFB<ndim>::AmbientTemp
(Particle<ndim> &part)
{
  debug2("[RadiativeFB::AmbientTemp]");

  FLOAT temp = 0.0;

  if (ambient_heating) temp += temp_inf4;
  if (disc_heating) temp += disc_heating->AmbientTemp(part, sinks);
  if (sink_heating) temp += sink_heating->AmbientTemp(part, sinks);

  return pow(temp, 0.25);
}

//=================================================================================================
//  DiscHeating::DiscHeating()
/// DiscHeating class constructor.
//=================================================================================================
template <int ndim>
DiscHeating<ndim>::DiscHeating(
 SimUnits *simunits,
 Parameters *params,
 int Ncentral_aux)
{
  Ncentral = Ncentral_aux;

  temp_unit = simunits->temp.outscale * simunits->temp.outSI;
  temp_au = params->floatparams["temp_au"] / temp_unit;
  temp_q = params->floatparams["temp_q"];

  runit = simunits->r.outscale * simunits->r.outSI;
  rsmooth = params->floatparams["r_smooth"] / runit;

  temp_au4 = pow(temp_au, 4.0);
  temp_exp = -2.0 * temp_q;

  rsmooth2 = pow(rsmooth, 2.0);
}

//=================================================================================================
//  DiscHeating::DiscHeating()
/// DiscHeating class destructor.
//=================================================================================================
template <int ndim>
DiscHeating<ndim>::~DiscHeating()
{

}

//=================================================================================================
//  DiscHeating::DiscHeating()
/// Calculates the ambient temperature for particles assuming a fixed temperature profile.
/// Depending on the case, the temperature profile may be imposed due to a single or binary
/// star. This routine should only be used for disc systems where the central object temperature
/// profile is fixed and does not depend on accretion.
//=================================================================================================
template <int ndim>
FLOAT DiscHeating<ndim>::AmbientTemp
(Particle<ndim> &part,
 Sinks<ndim> *sinks) {
  FLOAT temp = 0.0;

  if (Ncentral <= sinks->Nsink) {
    for (int i = 0; i < Ncentral; ++i) {
      SinkParticle<ndim> sink = sinks->sink[i];
      FLOAT dist = Distance(part.r, sink.star->r, 2); // Midplane distance only
      temp = temp_au4 * pow((dist * dist + rsmooth2), temp_exp);
    }
  }
  return temp;
}

//=================================================================================================
//  SinkHeating::SinkHeating()
/// SinkHeating class constructor.
//=================================================================================================
template <int ndim>
SinkHeating<ndim>::SinkHeating
(SimUnits *simunits,
 Parameters *params,
 int Ncentral_aux)
{
  Ncentral = Ncentral_aux;

  // Set variables from params file
  f_acc    = params->floatparams["f_acc"];
  r_star   = params->floatparams["r_star"];
  r_bdwarf = params->floatparams["r_bdwarf"];
  r_planet = params->floatparams["r_planet"];

  // Unit conversion, who doesn't love this?
  DOUBLE num, denom;

  temp_unit = simunits->temp.outscale * simunits->temp.outSI;

  // Calculate Boltzmann constant in code units
  num       = pow(simunits->r.outscale * simunits->r.outSI, 2.0) *
              simunits->t.outscale * simunits->t.outSI;
  denom     = simunits->E.outscale * simunits->E.outSI;
  rad_const = stefboltz * (num * pow(temp_unit, 4.0)) / denom;

  // Calculate L_sun in code units
  denom     = simunits->L.outscale * simunits->L.outSI;
  lsun 	    = L_sun / denom;

  // Calculate M_sun and M_jup in code units
  denom     = simunits->m.outscale * simunits->m.outSI;
  msun    	= m_sun / denom;
  mjup      = m_jup / denom;

  // Calculate R_sun in code units
  denom     = simunits->r.outscale * simunits->r.outSI;
  rsun 	    = r_sun / denom;

  // Convert object sizes to code units
  r_planet *= rsun;
  r_bdwarf *= rsun;
  r_star   *= rsun;
}

//=================================================================================================
//  SinkHeating::SinkHeating()
/// SinkHeating class destructor.
//=================================================================================================
template <int ndim>
SinkHeating<ndim>::~SinkHeating()
{

}

//=================================================================================================
//  SinkHeating::SinkTemperature()
/// Returns the sink temperature which is found via it's luminosity, such that:
/// T = (L / (4 * pi * boltz * R^2))^1/4
//=================================================================================================
template <int ndim>
FLOAT SinkHeating<ndim>::SinkTemperature
(FLOAT L,
 FLOAT r)
{
  return pow(L / (4.0 * pi * rad_const * r * r), 0.25);
}

//=================================================================================================
//  SinkHeating::SinkLuminosity()
/// Returns the sink luminosity which consists of two parts: an intrinsic luminosity which only
/// applies to stars (objects M > 80 Mjup) and an accretion luminosity. Specifically
/// L = f_n * (m / m_sun)^3 * L_sun + f_acc * (G * m * m_dot / r) * (1 - (r / 2 * rsink))
/// f_n is the flag for intrinsic accretion. f_acc is the fraction of gravitational energy
/// converted to heat energy (typically f_acc = 0.75, Offner et al. (2010)). R is the radius of
/// the object at the center of the sink and Rsink is the radius of the sink itself.
//=================================================================================================
template <int ndim>
FLOAT SinkHeating<ndim>::SinkLuminosity
(FLOAT m,
 FLOAT mdot,
 FLOAT rsink,
 FLOAT r,
 int f_n)
{
  return f_n * pow(m / msun, 3.0) * lsun +
         f_acc * ((m * mdot) / r) * (1 - (r / (2.0 * rsink)));
}

//=================================================================================================
//  SinkHeating::AmbientTemp()
/// Calculates the ambient temperature of particles assuming continuous radiative feedback, i.e.
/// when particles are accreted onto an object, the energy is instantly released into the system.
//=================================================================================================
template <int ndim>
FLOAT SinkHeating<ndim>::AmbientTemp
(Particle<ndim> &part,
 Sinks<ndim> *sinks)
{
  FLOAT temp = 0.0;

  for (int i = Ncentral; i < sinks->Nsink; ++i) {
    SinkParticle<ndim> sink = sinks->sink[i];

    FLOAT dist      = Distance(part.r, sink.star->r, ndim);
    FLOAT sink_m    = sink.star->m;
    FLOAT sink_dmdt = sink.dmdt;
    FLOAT sink_r    = sink.radius;

    // Set source radius and intrinsic luminosity flag (f_n) depending on sink mass. If the sink
    // is over the hydrogen burning limit, we want to give it some intrinsic luminosity.
    // Planet masses below 13.0 Mjupiter
    // Brown dwarf masses between 13.0 and 80.0 Mjupiter
    // Stellar masses aobe 80.0 Mjupiter
    FLOAT r_source = r_planet;
    int f_n = 0;
    if (sink_m >= 13.0 * mjup) r_source = r_bdwarf;
    if (sink_m >= 80.0 * mjup) {
      r_source = r_star;
      f_n = 1;
    }

    FLOAT sink_lum = SinkLuminosity(sink_m, sink_dmdt, sink_r, r_source, f_n);
    FLOAT sink_temp = SinkTemperature(sink_lum, r_source);

    temp += 0.25 * pow(r_source / dist, 2.0) * pow(sink_temp, 4.0);
  }
  return temp;
}

template class RadiativeFB<1>;
template class RadiativeFB<2>;
template class RadiativeFB<3>;
template class DiscHeating<1>;
template class DiscHeating<2>;
template class DiscHeating<3>;
template class SinkHeating<1>;
template class SinkHeating<2>;
template class SinkHeating<3>;
