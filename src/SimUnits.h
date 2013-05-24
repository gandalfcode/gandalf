//=============================================================================
//  SimUnits.h
//  Contains class definitions for all unit scaling classes.
//=============================================================================


#ifndef _SIM_UNITS_H_
#define _SIM_UNITS_H_


#include <iostream>
#include <map>
#include <string>
#include "Precision.h"
#include "Parameters.h"
using namespace std;


//=============================================================================
//  Class SimUnit
/// \brief   Main parent class for each individual unit class.
/// \details ..
/// \author  D. A. Hubber, G. Rosotti
/// \date    03/04/2013
//=============================================================================
class SimUnit
{
 public:

  SimUnit();

  virtual DOUBLE SIUnit(string) = 0;
  virtual string LatexLabel(string) = 0;
  DOUBLE OutputScale(string);

  DOUBLE inscale;                       ///< Input scaling factor
  DOUBLE inSI;                          ///< Input SI scaling factor
  DOUBLE outcgs;                        ///< Output cgs scaling factor
  DOUBLE outscale;                      ///< Output scaling factor
  DOUBLE outSI;                         ///< Output SI scaling factor
  string inunit;                        ///< Input unit string
  string outunit;                       ///< Output unit string

};



//=============================================================================
//  Class LengthUnit
//=============================================================================
class LengthUnit: public SimUnit
{
 public:
  LengthUnit() : SimUnit() {};
  DOUBLE SIUnit(string);
  string LatexLabel(string);

};



//=============================================================================
//  Class MassUnit
//=============================================================================
class MassUnit: public SimUnit
{
 public:
  MassUnit() : SimUnit() {};
  DOUBLE SIUnit(string);
  string LatexLabel(string);
 
};



//=============================================================================
//  Class TimeUnit
//=============================================================================
class TimeUnit: public SimUnit
{
 public:
  TimeUnit() : SimUnit() {};
  DOUBLE SIUnit(string);
  string LatexLabel(string);

};



//=============================================================================
//  Class VelocityUnit
//=============================================================================
class VelocityUnit: public SimUnit
{
 public:
  VelocityUnit() : SimUnit() {};
  DOUBLE SIUnit(string);
  string LatexLabel(string);

};



//=============================================================================
//  Class AccelerationUnit
//=============================================================================
class AccelerationUnit: public SimUnit
{
 public:
  AccelerationUnit() : SimUnit() {};
  DOUBLE SIUnit(string);
  string LatexLabel(string);

};



//=============================================================================
//  Class DensityUnit
//=============================================================================
class DensityUnit: public SimUnit
{
 public:
  DensityUnit() : SimUnit() {};
  DOUBLE SIUnit(string);
  string LatexLabel(string);

};


/*
//=============================================================================
//  Class ColumnDensityUnit
//=============================================================================
class ColumnDensityUnit: public SimUnit
{
 public:

  DOUBLE SIUnit(string);
  string LatexLabel(string);

};
*/



//=============================================================================
//  Class PressureUnit
//=============================================================================
class PressureUnit: public SimUnit
{
 public:
  PressureUnit() : SimUnit() {};
  DOUBLE SIUnit(string);
  string LatexLabel(string);

};



//=============================================================================
//  Class ForceUnit
//=============================================================================
class ForceUnit: public SimUnit
{
 public:

  DOUBLE SIUnit(string);
  string LatexLabel(string);

};



//=============================================================================
//  Class EnergyUnit
//=============================================================================
class EnergyUnit: public SimUnit
{
 public:
  EnergyUnit() : SimUnit() {};
  DOUBLE SIUnit(string);
  string LatexLabel(string);

};



//=============================================================================
//  Class MomentumUnit
//=============================================================================
class MomentumUnit: public SimUnit
{
 public:
  MomentumUnit() : SimUnit() {};
  DOUBLE SIUnit(string);
  string LatexLabel(string);

};




//=============================================================================
//  Class ColumnDensityUnit
//=============================================================================
class AngularMomentumUnit: public SimUnit
{
 public:
  AngularMomentumUnit() : SimUnit() {};
  DOUBLE SIUnit(string);
  string LatexLabel(string);

};



//=============================================================================
//  Class AngularVelocityUnit
//=============================================================================
class AngularVelocityUnit: public SimUnit
{
 public:
  AngularVelocityUnit() : SimUnit() {};
  DOUBLE SIUnit(string);
  string LatexLabel(string);

};



//=============================================================================
//  Class MassRateUnit
//=============================================================================
class MassRateUnit: public SimUnit
{
 public:
  MassRateUnit() : SimUnit() {};
  DOUBLE SIUnit(string);
  string LatexLabel(string);

};


/*
//=============================================================================
//  Class LuminosityUnit
//=============================================================================
class LuminosityUnit: public SimUnit
{
 public:
  LuminosityUnit() : SimUnit() {};
  DOUBLE SIUnit(string);
  string LatexLabel(string);

};
*/


/*
//=============================================================================
//  Class OpacityUnit
//=============================================================================
class OpacityUnit: public SimUnit
{
 public:

  DOUBLE SIUnit(string);
  string LatexLabel(string);

};



//=============================================================================
//  Class MagneticFieldUnit
//=============================================================================
class MagneticFieldUnit: public SimUnit
{
 public:

  DOUBLE SIUnit(string);
  string LatexLabel(string);

};



//=============================================================================
//  Class ChargeUnit
//=============================================================================
class ChargeUnit: public SimUnit
{
 public:

  DOUBLE SIUnit(string);
  string LatexLabel(string);

};



//=============================================================================
//  Class CurrentDensityUnit
//=============================================================================
class CurrentDensityUnit: public SimUnit
{
 public:

  DOUBLE SIUnit(string);
  string LatexLabel(string);

};
*/


//=============================================================================
//  Class SpecificEnergyUnit
//=============================================================================
class SpecificEnergyUnit: public SimUnit
{
 public:
  SpecificEnergyUnit() : SimUnit() {};
  DOUBLE SIUnit(string);
  string LatexLabel(string);

};



//=============================================================================
//  Class SpecificEnergyRateUnit
//=============================================================================
class SpecificEnergyRateUnit: public SimUnit
{
 public:
  SpecificEnergyRateUnit() : SimUnit() {};
  DOUBLE SIUnit(string);
  string LatexLabel(string);

};



//=============================================================================
//  Class TemperatureUnit
//=============================================================================
class TemperatureUnit: public SimUnit
{
 public:
  TemperatureUnit() : SimUnit() {};
  DOUBLE SIUnit(string);
  string LatexLabel(string);

};



//=============================================================================
//  Class SimUnits
/// \brief   Main simulation scaling class containing an instance of each unit.
/// \details ..
/// \author  D. A. Hubber, G. Rosotti
/// \date    03/04/2013
//=============================================================================
class SimUnits
{
 public:

  SimUnits();
  ~SimUnits();

  void SetupUnits(Parameters *);

  int dimensionless;                ///< Are we using dimensionless units?
  bool ReadInputUnits;              ///< Are input units read from snapshot?


  // Instances of all unit classes
  // --------------------------------------------------------------------------
  LengthUnit r;
  MassUnit m;
  TimeUnit t;
  VelocityUnit v;
  AccelerationUnit a;
  DensityUnit rho;
  //ColumnDensityUnit sigma;
  PressureUnit press;
  ForceUnit f;
  EnergyUnit E;
  MomentumUnit mom;
  AngularMomentumUnit angmom;
  AngularVelocityUnit angvel;
  MassRateUnit dmdt;
  //LuminosityUnit L;
  //OpacityUnit kappa;
  //MagneticFieldUnit B;
  //ChargeUnit Q;
  //CurrentDensityUnit Jcur;
  SpecificEnergyUnit u;
  SpecificEnergyRateUnit dudt;
  TemperatureUnit temp;

};


#endif




