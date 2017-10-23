//=================================================================================================
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
//=================================================================================================


#ifndef _ENERGY_EQUATION_H_
#define _ENERGY_EQUATION_H_


#include "Constants.h"
#include "CodeTiming.h"
#include "EOS.h"
#include "Hydrodynamics.h"
#include "Particle.h"
#include "Precision.h"
#include "RadiativeFB.h"
#include "SimUnits.h"


//=================================================================================================
//  EnergyEquation
/// Main energy equation class, with virtual functions that require full
/// definitions in the children classes.
//=================================================================================================
template <int ndim>
class EnergyEquation
{
 public:

  EnergyEquation(DOUBLE);
  virtual ~EnergyEquation();

  virtual void EnergyIntegration(const int,  const FLOAT, const FLOAT, Hydrodynamics<ndim> *) = 0;
  virtual void EnergyCorrectionTerms(const int, const FLOAT, const FLOAT, Hydrodynamics<ndim> *) = 0;
  virtual void EndTimestep(const int, const FLOAT, const FLOAT, Hydrodynamics<ndim> *) = 0;
  virtual DOUBLE Timestep(Particle<ndim> &) = 0;


  const DOUBLE energy_mult;            ///< Explicit integration timestep multiplier
  CodeTiming *timing;                  ///< Pointer to code timing object

};



//=================================================================================================
//  EnergyRadws
/// Energy equation class using Stamatellos et al. (2007) radiation cooling scheme.
//=================================================================================================
template <int ndim, template <int> class ParticleType>
class EnergyRadws : public EnergyEquation<ndim>
{
 public:

  using EnergyEquation<ndim>::timing;

  EnergyRadws(Parameters*, SimUnits *, Radws<ndim> *, RadiativeFB<ndim> *);
  ~EnergyRadws();

  //  void ReadTable();
  void EnergyIntegration(const int,  const FLOAT, const FLOAT, Hydrodynamics<ndim> *);
  void EnergyCorrectionTerms(const int, const FLOAT, const FLOAT, Hydrodynamics<ndim> *) {};
  void EndTimestep(const int, const FLOAT, const FLOAT, Hydrodynamics<ndim> *);
  void EnergyFindEqui(const FLOAT, const FLOAT, const FLOAT, const FLOAT,
                      const FLOAT, const FLOAT, FLOAT &, FLOAT &);
  void EnergyFindEquiTemp(const int, const FLOAT, const FLOAT, const FLOAT,
                          const FLOAT, const FLOAT, FLOAT &);
  DOUBLE Timestep(Particle<ndim> &) {return big_number_dp;}

  FLOAT ebalance(const FLOAT, const FLOAT, const FLOAT, const FLOAT, const FLOAT, const FLOAT);
  FLOAT GetCol2(Particle<ndim> &part);


  //-----------------------------------------------------------------------------------------------

  int ndens, ntemp;
  int lombardi;
  FLOAT fcol2;
  FLOAT rad_const;
  FLOAT temp_ambient0;
  FLOAT temp_min;
  Radws<ndim> *eos;
  OpacityTable<ndim> *table;
  RadiativeFB<ndim> *radfb;
};



//=================================================================================================
//  NullEnergy
/// Null (empty) class when no energy option is selected.
//=================================================================================================
template <int ndim>
class NullEnergy : public EnergyEquation<ndim>
{
 public:

  NullEnergy(DOUBLE dt_mult) : EnergyEquation<ndim>(dt_mult) {};
  ~NullEnergy() {};

  void EnergyIntegration(const int,  const FLOAT, const FLOAT, Hydrodynamics<ndim> *) {};
  void EnergyCorrectionTerms(const int, const FLOAT, const FLOAT, Hydrodynamics<ndim> *) {};
  void EndTimestep(const int, const FLOAT, const FLOAT, Hydrodynamics<ndim> *) {};
  DOUBLE Timestep(Particle<ndim> &) {return big_number_dp;}

};
#endif
