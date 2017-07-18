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
/// RadiativeFB class constructor
//=================================================================================================
template <int ndim>
RadiativeFB<ndim>::RadiativeFB
(SimUnits *simunits,
 Parameters *params)
{
  // Set the radiative feedback regime
  regime = params->stringparams["rad_fb"];

  // Convert to units first
  temp_unit = simunits->temp.outscale * simunits->temp.outcgs;
  temp_inf = params->floatparams["temp_ambient"] / temp_unit;
  temp_au = params->floatparams["temp_au"] / temp_unit;
  temp_q = params->floatparams["temp_q"];

  runit = simunits->r.outscale * simunits->r.outSI;
  runit2 = pow(runit, 2.0);
  rsmooth = params->floatparams["r_smooth"];

  temp_inf4 = pow(temp_inf, 4.0);
  temp_au4 = pow(temp_au, 4.0);
  temp_exp = -2.0 * temp_q;
  rsmooth2 = pow(rsmooth, 2.0);

  return;
}

//=================================================================================================
//  RadiativeFB::~RadiativeFB()
/// RadiativeFB class destructor
//=================================================================================================
template <int ndim>
RadiativeFB<ndim>::~RadiativeFB()
{
  delete sinks;
  delete simunits;
}

//=================================================================================================
//  RadiativeFB::AmbientTemp()
/// Calculates the ambient temperature for particles assuming a fixed temperature profile.
/// Depending on the case, the temperature profile may be imposed due to a single or binary
/// star. This routine should only be used for disc systems where the central object temperature
/// profile is fixed and does not depend on accretion.
//=================================================================================================
template <int ndim>
void RadiativeFB<ndim>::AmbientTemp
(Hydrodynamics<ndim> *hydro,
 Sinks<ndim> *sinks)
{
  debug2("[RadiativeFB::AmbientTemp]");

  if (regime == "hdisc_single") {
    FLOAT sink_x = sinks->sink[0].star->r[0];
    FLOAT sink_y = sinks->sink[0].star->r[1];

#pragma omp parallel for default(none) shared(sink_x, sink_y, hydro)
    for (int i = 0; i < hydro->Nhydro; ++i) {
      Particle<ndim> &part = hydro->GetParticlePointer(i);
      FLOAT part_x = part.r[0];
      FLOAT part_y = part.r[1];
      FLOAT r2 = pow(part_x - sink_x, 2.0) + pow(part_y - sink_y, 2.0);

      // Requires unit conversion, seems to work okay with AU units though
      part.temp_ambient = pow(temp_inf4 + temp_au4 * pow((r2 + rsmooth2), temp_exp), 0.25);
    }
  }
  else if (regime == "hdisc_binary") {
    // TODO (MERCER): Implement heating from central binary system
    // How can we parameterise primary and secondary mass? T(1AU) from
    // luminosity data. q will have to be assumed.
  }

  return;
}

//=================================================================================================
//  ContinuousFB::ContinuousFB()
/// ContinuousFB class constructor
//=================================================================================================
template <int ndim>
ContinuousFB<ndim>::ContinuousFB
(SimUnits *simunits,
 Parameters *params) : RadiativeFB<ndim>(simunits, params)
{
  // Set variables from params file
  type     = params->stringparams["fb_type"];
  f_acc    = params->floatparams["f_acc"];
  r_star   = params->floatparams["r_star"];
  r_bdwarf = params->floatparams["r_bdwarf"];
  r_planet = params->floatparams["r_planet"];

  // Unit conversion, who doesn't love this?
  DOUBLE num, denom;

  // Calculate Boltzmann constant in code units
  num       = pow(simunits->r.outscale * simunits->r.outSI, 2.0) *
              simunits->t.outscale * simunits->t.outSI;
  denom     = simunits->E.outscale * simunits->E.outSI;
  rad_const = stefboltz * (num * pow(temp_unit, 4.0)) / denom;

  // Calculate gravitational constant in code units
  num        = pow(simunits->t.outscale * simunits->t.outSI, 2.0) *
               simunits->m.outscale * simunits->m.outSI;
  denom      = pow(simunits->r.outscale * simunits->r.outSI, 3.0);
  grav_const = G_const * (num / denom);

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

// //=================================================================================================
// //  ContinuousFB::~ContinuousFB()
// /// ContinuousFB class destructor
// //=================================================================================================
template <int ndim>
ContinuousFB<ndim>::~ContinuousFB()
{

}

//=================================================================================================
//  ContinuousFB::AmbientTemp()
/// Calculates the ambient temperature of particles assuming continuous radiative feedback, i.e.
/// when particles are accreted onto an object, the energy is instantly released into the system.
//=================================================================================================
template <int ndim>
void ContinuousFB<ndim>::AmbientTemp
(Hydrodynamics<ndim> *hydro,
 Sinks<ndim> *sinks)
{
  if (regime == "sink_heating") {
#pragma omp parallel for default(none) shared(sinks, hydro)
    for (int i = 0; i < hydro->Nhydro; ++i) {
      Particle<ndim> &part = hydro->GetParticlePointer(i);
      FLOAT part_x = part.r[0];
      FLOAT part_y = part.r[1];
      FLOAT part_z = part.r[2];

      FLOAT sum = 0.0;
      for (int j = 1; j < sinks->Nsink; ++j) {
        FLOAT sink_x    = sinks->sink[j].star->r[0];
        FLOAT sink_y    = sinks->sink[j].star->r[1];
        FLOAT sink_z    = sinks->sink[j].star->r[2];
        FLOAT sink_m    = sinks->sink[j].star->m;
        FLOAT sink_mdot = sinks->sink[j].dmdt;
        FLOAT sink_r    = sinks->sink[j].radius;

        // Set source radius and intrinsic luminosity flag depending on sink mass
        FLOAT r_source = r_planet;
        int f_n = 0;
        if (sink_m >= 13.0 * mjup) r_source = r_bdwarf;
        if (sink_m >= 80.0 * mjup) {
          r_source = r_star;
          f_n = 1;
        }

        FLOAT sink_lum = SinkLuminosity(sink_m, sink_mdot, sink_r, r_source, f_n);
        FLOAT sink_temp = SinkTemperature(sink_lum, r_source);

        FLOAT r_source2 = r_source * r_source;
        FLOAT dist = pow(part_x - sink_x, 2.0) +
                     pow(part_y - sink_y, 2.0) +
                     pow(part_z - sink_z, 2.0);

        sum += 0.25 * (r_source2 / dist) * pow(sink_temp, 4.0);
      }

      part.temp_ambient = pow(temp_inf4 + sum, 0.25);
    }
  }
  else if (regime == "hdisc_single_plus_sink_heating") {
    // TODO (MERCER): Implement heating from single central star and formed sinks
  }
  else if (regime == "hdisc_binary_plus_sink_heating") {
    // TODO (MERCER): Implement heating from binary central system and formed sinks
  }
}

//=================================================================================================
//  ContinuousFB::SinkTemperature()
/// Returns the sink temperature which is found via it's luminosity, such that:
/// T = (L / (4 * pi * boltz * R^2))^1/4
//=================================================================================================
template <int ndim>
FLOAT ContinuousFB<ndim>::SinkTemperature
(FLOAT L,
 FLOAT r)
{
  return pow(L / (4.0 * pi * rad_const * r), 0.25);
}

//=================================================================================================
//  ContinuousFB::SinkLuminosity()
/// Returns the sink luminosity which consists of two parts: an intrinsic luminosity which only
/// applies to stars (objects M > 80 Mjup) and an accretion luminosity. Specifically
/// L = f_n * (m / m_sun)^3 * L_sun + f_acc * (G * m * m_dot / r) * (1 - (r / 2 * rsink))
/// f_n is the flag for intrinsic accretion. f_acc is the fraction of gravitational energy
/// converted to heat energy (typically f_acc = 0.75, Offner et al. (2010)). R is the radius of
/// the object at the center of the sink and Rsink is the radius of the sink itself.
//=================================================================================================
template <int ndim>
FLOAT ContinuousFB<ndim>::SinkLuminosity
(FLOAT m,
 FLOAT mdot,
 FLOAT rsink,
 FLOAT r,
 int f_n)
{
  return f_n * pow(m / msun, 3.0) * lsun +
         f_acc * ((grav_const * m * mdot) / r) * (1 - (r / (2.0 * rsink)));
}

template class RadiativeFB<1>;
template class RadiativeFB<2>;
template class RadiativeFB<3>;
template class ContinuousFB<1>;
template class ContinuousFB<2>;
template class ContinuousFB<3>;
