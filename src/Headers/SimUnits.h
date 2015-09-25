//=================================================================================================
//  SimUnits.h
//  Contains definitions for all unit scaling classes.
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


#ifndef _SIM_UNITS_H_
#define _SIM_UNITS_H_


#include <iostream>
#include <map>
#include <string>
#include "Precision.h"
#include "Parameters.h"
using namespace std;


//=================================================================================================
//  Class SimUnit
/// \brief   Main parent class for each individual unit class.
/// \details Main parent class for each individual unit class.
/// \author  D. A. Hubber, G. Rosotti
/// \date    03/04/2013
//=================================================================================================
class SimUnit
{
 public:

  SimUnit();
  virtual ~SimUnit() {};

  virtual DOUBLE SIUnit(string) = 0;
  virtual string LatexLabel(string) = 0;
  DOUBLE OutputScale(string);

  DOUBLE inscale;                          ///< Input scaling factor
  DOUBLE inSI;                             ///< Input SI scaling factor
  DOUBLE outcgs;                           ///< Output cgs scaling factor
  DOUBLE outscale;                         ///< Output scaling factor
  DOUBLE outSI;                            ///< Output SI scaling factor
  string inunit;                           ///< Input unit string
  string outunit;                          ///< Output unit string

};



//=================================================================================================
//  Class LengthUnit
//=================================================================================================
class LengthUnit: public SimUnit
{
 public:
  LengthUnit() : SimUnit() {};
  DOUBLE SIUnit(string);
  string LatexLabel(string);

};



//=================================================================================================
//  Class MassUnit
//=================================================================================================
class MassUnit: public SimUnit
{
 public:
  MassUnit() : SimUnit() {};
  DOUBLE SIUnit(string);
  string LatexLabel(string);

};



//=================================================================================================
//  Class TimeUnit
//=================================================================================================
class TimeUnit: public SimUnit
{
 public:
  TimeUnit() : SimUnit() {};
  DOUBLE SIUnit(string);
  string LatexLabel(string);

};



//=================================================================================================
//  Class VelocityUnit
//=================================================================================================
class VelocityUnit: public SimUnit
{
 public:
  VelocityUnit() : SimUnit() {};
  DOUBLE SIUnit(string);
  string LatexLabel(string);

};



//=================================================================================================
//  Class AccelerationUnit
//=================================================================================================
class AccelerationUnit: public SimUnit
{
 public:
  AccelerationUnit() : SimUnit() {};
  DOUBLE SIUnit(string);
  string LatexLabel(string);

};



//=================================================================================================
//  Class DensityUnit
//=================================================================================================
class DensityUnit: public SimUnit
{
 public:
  DensityUnit() : SimUnit() {};
  DOUBLE SIUnit(string);
  string LatexLabel(string);

};



//=================================================================================================
//  Class ColumnDensityUnit
//=================================================================================================
class ColumnDensityUnit: public SimUnit
{
 public:
  ColumnDensityUnit() : SimUnit() {};
  DOUBLE SIUnit(string);
  string LatexLabel(string);

};



//=================================================================================================
//  Class PressureUnit
//=================================================================================================
class PressureUnit: public SimUnit
{
 public:
  PressureUnit() : SimUnit() {};
  DOUBLE SIUnit(string);
  string LatexLabel(string);

};



//=================================================================================================
//  Class ForceUnit
//=================================================================================================
class ForceUnit: public SimUnit
{
 public:
  ForceUnit() : SimUnit() {};
  DOUBLE SIUnit(string);
  string LatexLabel(string);

};



//=================================================================================================
//  Class EnergyUnit
//=================================================================================================
class EnergyUnit: public SimUnit
{
 public:
  EnergyUnit() : SimUnit() {};
  DOUBLE SIUnit(string);
  string LatexLabel(string);

};



//=================================================================================================
//  Class MomentumUnit
//=================================================================================================
class MomentumUnit: public SimUnit
{
 public:
  MomentumUnit() : SimUnit() {};
  DOUBLE SIUnit(string);
  string LatexLabel(string);

};




//=================================================================================================
//  Class AngularMomentumUnit
//=================================================================================================
class AngularMomentumUnit: public SimUnit
{
 public:
  AngularMomentumUnit() : SimUnit() {};
  DOUBLE SIUnit(string);
  string LatexLabel(string);

};



//=================================================================================================
//  Class AngularVelocityUnit
//=================================================================================================
class AngularVelocityUnit: public SimUnit
{
 public:
  AngularVelocityUnit() : SimUnit() {};
  DOUBLE SIUnit(string);
  string LatexLabel(string);

};



//=================================================================================================
//  Class MassRateUnit
//=================================================================================================
class MassRateUnit: public SimUnit
{
 public:
  MassRateUnit() : SimUnit() {};
  DOUBLE SIUnit(string);
  string LatexLabel(string);

};



//=================================================================================================
//  Class LuminosityUnit
//=================================================================================================
class LuminosityUnit: public SimUnit
{
 public:
  LuminosityUnit() : SimUnit() {};
  DOUBLE SIUnit(string);
  string LatexLabel(string);

};




//=================================================================================================
//  Class OpacityUnit
//=================================================================================================
class OpacityUnit: public SimUnit
{
 public:
  OpacityUnit() : SimUnit() {};
  DOUBLE SIUnit(string);
  string LatexLabel(string);

};



//=================================================================================================
//  Class MagneticFieldUnit
//=================================================================================================
class MagneticFieldUnit: public SimUnit
{
 public:
  MagneticFieldUnit() : SimUnit() {};
  DOUBLE SIUnit(string);
  string LatexLabel(string);

};



//=================================================================================================
//  Class ChargeUnit
//=================================================================================================
class ChargeUnit: public SimUnit
{
 public:
  ChargeUnit() : SimUnit() {};
  DOUBLE SIUnit(string);
  string LatexLabel(string);

};



//=================================================================================================
//  Class CurrentDensityUnit
//=================================================================================================
class CurrentDensityUnit: public SimUnit
{
 public:
  CurrentDensityUnit() : SimUnit() {};
  DOUBLE SIUnit(string);
  string LatexLabel(string);

};



//=================================================================================================
//  Class SpecificEnergyUnit
//=================================================================================================
class SpecificEnergyUnit: public SimUnit
{
 public:
  SpecificEnergyUnit() : SimUnit() {};
  DOUBLE SIUnit(string);
  string LatexLabel(string);

};



//=================================================================================================
//  Class SpecificEnergyRateUnit
//=================================================================================================
class SpecificEnergyRateUnit: public SimUnit
{
 public:
  SpecificEnergyRateUnit() : SimUnit() {};
  DOUBLE SIUnit(string);
  string LatexLabel(string);

};



//=================================================================================================
//  Class TemperatureUnit
//=================================================================================================
class TemperatureUnit: public SimUnit
{
 public:
  TemperatureUnit() : SimUnit() {};
  DOUBLE SIUnit(string);
  string LatexLabel(string);

};



//=================================================================================================
//  Class DimensionlessUnit
//=================================================================================================
class DimensionlessUnit: public SimUnit
{
 public:
  DimensionlessUnit() : SimUnit() {};
  DOUBLE SIUnit(string) { return 1.0; }
  string LatexLabel(string) { return ""; }

};



//=================================================================================================
//  Class SimUnits
/// \brief   Main simulation scaling class containing an instance of each unit.
/// \details Main simulation scaling class containing an instance of each unit.
/// \author  D. A. Hubber, G. Rosotti
/// \date    03/04/2013
//=================================================================================================
class SimUnits
{
 public:

  SimUnits();
  ~SimUnits();

  void SetupUnits(Parameters *);
  void OutputScalingFactors(Parameters *);

  int dimensionless;                ///< Are we using dimensionless units?
  bool ReadInputUnits;              ///< Are input units read from snapshot?


  // Instances of all unit classes
  //-----------------------------------------------------------------------------------------------
  LengthUnit r;                     ///< Length unit
  MassUnit m;                       ///< Mass unit
  TimeUnit t;                       ///< Time unit
  VelocityUnit v;                   ///< Velocity unit
  AccelerationUnit a;               ///< Acceleration unit
  DensityUnit rho;                  ///< Density unit
  ColumnDensityUnit sigma;          ///< Column density unit
  PressureUnit press;               ///< Pressure unit
  ForceUnit f;                      ///< Force unit
  EnergyUnit E;                     ///< Energy unit
  MomentumUnit mom;                 ///< Linear momentum unit
  AngularMomentumUnit angmom;       ///< Angular momentum unit
  AngularVelocityUnit angvel;       ///< Angular velocity unit
  MassRateUnit dmdt;                ///< Mass (accretion) rate unit
  LuminosityUnit L;                 ///< Luminosity unit
  OpacityUnit kappa;                ///< Volume opacity unit
  MagneticFieldUnit B;              ///< Magnetic field unit
  ChargeUnit Q;                     ///< Charge unit
  CurrentDensityUnit Jcur;          ///< Current density unit
  SpecificEnergyUnit u;             ///< Specific internal energy unit
  SpecificEnergyRateUnit dudt;      ///< Rate of change of internal energy unit
  TemperatureUnit temp;             ///< Temperature unit
  DimensionlessUnit nounits;        ///< Dimensionless unit

};
#endif
