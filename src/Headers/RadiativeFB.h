//=================================================================================================
//  RadiativeFB.h
//  Class definitions of radiative feedback methods.
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

#ifndef _RADIATIVEFB_H_
#define _RADIATIVEFB_H_

#include "Constants.h"
#include "Particle.h"
#include "Precision.h"
#include "SimUnits.h"
#include "Hydrodynamics.h"
#include "Sinks.h"

template <int ndim> class DiscHeating;
template <int ndim> class SinkHeating;

//=================================================================================================
//  Class RadiativeFB
/// \brief   Handles the radiative feedback object and find the ambient temperature for the
///          defined methods: ambient heating, disc heating and sink heating.
/// \details Radiative feedback is heating provided by the accretion of material from sink
///          particles within the simulation. This class is the main parent to different radiative
///          feedback methods. The temperature contribution is parameterised by an "ambient"
///          particle temperature. This is used within the RadWS EoS scheme to control the heating
///          and cooling of a particle.
/// \author  A. P. Mercer
/// \date    10/07/2017
//=================================================================================================
template <int ndim>
class RadiativeFB
{
public:
  RadiativeFB(SimUnits *, Parameters *);
  ~RadiativeFB();

  FLOAT AmbientTemp(Particle<ndim> &);

  void SetSinks(Sinks<ndim> *sinks_aux) { sinks = sinks_aux; }

  FLOAT temp_unit;
  FLOAT temp_inf;
  FLOAT temp_inf4;

  SimUnits *simunits;
  Sinks<ndim> *sinks;

  int ambient_heating;
  DiscHeating<ndim> *disc_heating;
  SinkHeating<ndim> *sink_heating;
};

//=================================================================================================
//  Class DiscHeating
/// \brief   Handles heating from a central protostellar system which affects the disc. The
///          central system may be single or a binary.
/// \details This class creates a proxy temperature profile across the disc from the central
///          protostellar system.
/// \author  A. P. Mercer
/// \date    28/07/2017
//=================================================================================================
template <int ndim>
class DiscHeating {
public:
  DiscHeating(SimUnits *, Parameters *, int);
  ~DiscHeating();

  FLOAT AmbientTemp(Particle<ndim> &, Sinks<ndim> *);

private:
  int Ncentral;

  FLOAT temp_au;
  FLOAT temp_q;
  FLOAT rsmooth;

  FLOAT temp_unit;
  FLOAT temp_au4;
  FLOAT temp_exp;
  FLOAT runit;
  FLOAT runit2;
  FLOAT rsmooth2;
};

//=================================================================================================
//  Class SinkHeating
/// \brief   Handles heating from sinks within the system due to accretion, and, if they are
///          massive enough, intrinsic luminosity.
/// \details This class provides heating from sinks: accretion heating and heating from hydrogen
///          burning. If DiscHeating is also turned on, the central protostellar system is excluded
///          from this method. Only formed sinks will then provide extra heating.
/// \author  A. P. Mercer
/// \date    28/07/2017
//=================================================================================================
template <int ndim>
class SinkHeating  {
public:
  SinkHeating(SimUnits *, Parameters *, int);
  ~SinkHeating();

  FLOAT AmbientTemp(Particle<ndim> &, Sinks<ndim> *);

  FLOAT SinkTemperature(FLOAT, FLOAT);
  FLOAT SinkLuminosity(FLOAT, FLOAT, FLOAT, FLOAT, int);

  int Ncentral;
  FLOAT temp_unit;
  DOUBLE rad_const;
  FLOAT mjup, msun, lsun, rsun;

  FLOAT r_star;
  FLOAT r_bdwarf;
  FLOAT r_planet;
  FLOAT f_acc;

};

#endif
