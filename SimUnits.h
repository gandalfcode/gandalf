// ============================================================================
// SimUnits.h
// ============================================================================


#ifndef _SIM_UNITS_H_
#define _SIM_UNITS_H_


#include <iostream>
#include <map>
#include <string>
#include "Parameters.h"
using namespace std;


// ============================================================================
// Class SimUnit
// ============================================================================
class SimUnit
{
 public:

  SimUnit();

  virtual double SIUnit(string) = 0;
  double OutputScale(string);

  double inscale;
  double inSI;
  double outcgs;
  double outscale;
  double outSI;
  string inunit;
  string outunit;

};



// ============================================================================
// Class LengthUnit
// ============================================================================
class LengthUnit: public SimUnit
{
 public:

  double SIUnit(string);

};



// ============================================================================
// Class MassUnit
// ============================================================================
class MassUnit: public SimUnit
{
 public:

  double SIUnit(string);
 
};



// ============================================================================
// Class TimeUnit
// ============================================================================
class TimeUnit: public SimUnit
{
 public:

  double SIUnit(string);

};



// ============================================================================
// Class VelocityUnit
// ============================================================================
class VelocityUnit: public SimUnit
{
 public:

  double SIUnit(string);

};



// ============================================================================
// Class AccelerationUnit
// ============================================================================
class AccelerationUnit: public SimUnit
{
 public:

  double SIUnit(string);

};




// ============================================================================
// Class DensityUnit
// ============================================================================
class DensityUnit: public SimUnit
{
 public:

  double SIUnit(string);

};



// ============================================================================
// Class ColumnDensityUnit
// ============================================================================
class ColumnDensityUnit: public SimUnit
{
 public:

  double SIUnit(string);

};



// ============================================================================
// Class PressureUnit
// ============================================================================
class PressureUnit: public SimUnit
{
 public:

  double SIUnit(string);

};



// ============================================================================
// Class ForceUnit
// ============================================================================
class ForceUnit: public SimUnit
{
 public:

  double SIUnit(string);

};



// ============================================================================
// Class EnergyUnit
// ============================================================================
class EnergyUnit: public SimUnit
{
 public:

  double SIUnit(string);

};



// ============================================================================
// Class MomentumUnit
// ============================================================================
class MomentumUnit: public SimUnit
{
 public:

  double SIUnit(string);

};




// ============================================================================
// Class ColumnDensityUnit
// ============================================================================
class AngularMomentumUnit: public SimUnit
{
 public:

  double SIUnit(string);

};



// ============================================================================
// Class AngularVelocityUnit
// ============================================================================
class AngularVelocityUnit: public SimUnit
{
 public:

  double SIUnit(string);

};



// ============================================================================
// Class MassAccretionRateUnit
// ============================================================================
class MassAccretionRateUnit: public SimUnit
{
 public:

  double SIUnit(string);

};



// ============================================================================
// Class LuminosityUnit
// ============================================================================
class LuminosityUnit: public SimUnit
{
 public:

  double SIUnit(string);

};



// ============================================================================
// Class OpacityUnit
// ============================================================================
class OpacityUnit: public SimUnit
{
 public:

  double SIUnit(string);

};



// ============================================================================
// Class MagneticFieldUnit
// ============================================================================
class MagneticFieldUnit: public SimUnit
{
 public:

  double SIUnit(string);

};



// ============================================================================
// Class ChargeUnit
// ============================================================================
class ChargeUnit: public SimUnit
{
 public:

  double SIUnit(string);

};



// ============================================================================
// Class CurrentDensityUnit
// ============================================================================
class CurrentDensityUnit: public SimUnit
{
 public:

  double SIUnit(string);

};



// ============================================================================
// Class SpecificEnergyUnit
// ============================================================================
class SpecificEnergyUnit: public SimUnit
{
 public:

  double SIUnit(string);

};



// ============================================================================
// Class SpecificEnergyRateUnit
// ============================================================================
class SpecificEnergyRateUnit: public SimUnit
{
 public:

  double SIUnit(string);

};



// ============================================================================
// Class TemperatureUnit
// ============================================================================
class TemperatureUnit: public SimUnit
{
 public:

  double SIUnit(string);

};





// ============================================================================
// Class SimUnits
// ============================================================================
class SimUnits
{
 public:

  SimUnits();
  ~SimUnits();

  void SetupUnits(Parameters &);

  bool ReadInputUnits;

  LengthUnit r;
  MassUnit m;
  TimeUnit t;
  VelocityUnit v;
  AccelerationUnit a;
  DensityUnit rho;
  //ColumnDensityUnit sigma;
  //PressureUnit press;
  //ForceUnit f;
  EnergyUnit E;
  MomentumUnit mom;
  AngularMomentumUnit angmom;
  AngularVelocityUnit angvel;
  //MassAccretionRateUnit dmdt;
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




