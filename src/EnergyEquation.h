//=============================================================================
//  EnergyEquation.h
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
//=============================================================================


#ifndef _ENERGY_EQUATION_H_
#define _ENERGY_EQUATION_H_


#include "Constants.h"
#include "Precision.h"
#include "Sph.h"
#include "EOS.h"
#include "SphParticle.h"



//=============================================================================
//  EnergyEquation
/// Main energy equation class, with virtual functions that require full 
/// definitions in the children classes.
//=============================================================================
template <int ndim>
class EnergyEquation
{
 public:

  EnergyEquation(DOUBLE);
  ~EnergyEquation();

  virtual void EnergyIntegration(int, int, SphParticle<ndim> *, FLOAT) = 0;
  virtual void EnergyCorrectionTerms(int, int, SphParticle<ndim> *, FLOAT) = 0;
  virtual void EndTimestep(int, int, FLOAT, SphParticle<ndim> *) = 0;
  virtual DOUBLE Timestep(SphParticle<ndim> &) = 0;

  const DOUBLE energy_mult;
  CodeTiming *timing;               ///< Pointer to code timing object

};



//=============================================================================
//  EnergyPEC
/// Class definition for energy equation integration class using a 
/// Predict-Evaluate-Correct (PEC) scheme.
//=============================================================================
template <int ndim>
class EnergyPEC: public EnergyEquation<ndim>
{
 public:

  using EnergyEquation<ndim>::timing;

  EnergyPEC(DOUBLE);
  ~EnergyPEC();

  void EnergyIntegration(int, int, SphParticle<ndim> *, FLOAT);
  void EnergyCorrectionTerms(int, int, SphParticle<ndim> *, FLOAT);
  void EndTimestep(int, int, FLOAT, SphParticle<ndim> *);
  DOUBLE Timestep(SphParticle<ndim> &);

};



//=============================================================================
//  EnergyLeapfrogDKD
/// Class definition for energy equation integration class using a
/// Leapfrog drift-kick-drift scheme
//=============================================================================
template <int ndim>
class EnergyLeapfrogDKD: public EnergyEquation<ndim>
{
 public:

  EnergyLeapfrogDKD(DOUBLE);
  ~EnergyLeapfrogDKD();

  void EnergyIntegration(int, int, SphParticle<ndim> *, FLOAT);
  void EnergyCorrectionTerms(int, int, SphParticle<ndim> *, FLOAT);
  void EndTimestep(int, int, FLOAT, SphParticle<ndim> *);
  DOUBLE Timestep(SphParticle<ndim> &);

};



//=============================================================================
//  EnergyGodunovIntegration
/// ..
//=============================================================================
template <int ndim>
class EnergyGodunovIntegration: public EnergyEquation<ndim>
{
 public:

  EnergyGodunovIntegration(DOUBLE);
  ~EnergyGodunovIntegration();

  void EnergyIntegration(int, int, SphParticle<ndim> *, FLOAT);
  void EnergyCorrectionTerms(int, int, SphParticle<ndim> *, FLOAT);
  void EndTimestep(int, int, FLOAT, SphParticle<ndim> *);
  DOUBLE Timestep(SphParticle<ndim> &);

};
#endif
