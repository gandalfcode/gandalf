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


/*
// ============================================================================
// Class VelocityUnit
// ============================================================================
class VelocityUnit: public SimUnit
{
 public:

  double SIUnit(string);

};
*/


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
  //Velocity v;

};


#endif




