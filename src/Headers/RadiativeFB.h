//=================================================================================================
//  RadiativeFB.h
//  Class definitions of main energy equation class plus inherited children
//  classes for various energy integration algorithms.
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
#include "Nbody.h"

//=================================================================================================
//  Class RadiativeFB
/// \brief   Parent class to radiative feedback methods.
/// \details Radiative feedback is heating provided by the accretion of material from sink
///          particles within the simulation. This class is the main parent to different radiative
///          feedback methods.
/// \author  A. P. Mercer
/// \date    10/07/2017
//=================================================================================================
template <int ndim>
class RadiativeFB
{
public:
  RadiativeFB(SimUnits *, Parameters *);
  ~RadiativeFB();

  virtual void AmbientTemp(Hydrodynamics<ndim> *, Sinks<ndim> *);

  string regime;
  FLOAT temp_inf;
  FLOAT temp_au;
  FLOAT temp_q;
  FLOAT rsmooth;

  // Useful values
  FLOAT temp_inf4;
  FLOAT temp_au4;
  FLOAT temp_exp;
  FLOAT temp_unit;
  FLOAT runit;
  FLOAT runit2;
  FLOAT rsmooth2;

  Sinks<ndim> *sinks;
  SimUnits *simunits;
};

template <int ndim>
class ContinuousFB : public RadiativeFB<ndim>
{
  using RadiativeFB<ndim>::regime;
  using RadiativeFB<ndim>::temp_inf4;
  using RadiativeFB<ndim>::temp_unit;

  using RadiativeFB<ndim>::simunits;
  using RadiativeFB<ndim>::sinks;

public:
  ContinuousFB(SimUnits *, Parameters *);
  ~ContinuousFB();

  void AmbientTemp(Hydrodynamics<ndim> *, Sinks<ndim> *);

  string type;

  FLOAT rad_const, grav_const;
  FLOAT mjup, msun, lsun, rsun;

  FLOAT r_star;
  FLOAT r_bdwarf;
  FLOAT r_planet;
  FLOAT f_acc;

private:
  FLOAT SinkTemperature(FLOAT, FLOAT);
  FLOAT SinkLuminosity(FLOAT, FLOAT, FLOAT, FLOAT, int);
};

#endif
