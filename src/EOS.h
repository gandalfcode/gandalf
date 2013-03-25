// ============================================================================
// EOS.h
// Contains (virtual) definitions for equation of state class.  
// Also includes inherited class defintions for various equation of  
// state options.
// ============================================================================


#ifndef _EOS_H_
#define _EOS_H_


#include "Precision.h"
#include "Constants.h"
#include "Dimensions.h"
#include "Parameters.h"
#include "SphParticle.h"



// ============================================================================
// Class EOS
// Main Equation of state class.  Only contains virtual function defintions. 
// All functions must be defined by the inherited EOS classes.
// ============================================================================
class EOS
{
 public:

  virtual FLOAT Pressure(SphParticle &) = 0;
  virtual FLOAT EntropicFunction(SphParticle &) = 0;
  virtual FLOAT SoundSpeed(SphParticle &) = 0;
  virtual FLOAT Temperature(SphParticle &) = 0;
  virtual FLOAT SpecificInternalEnergy(SphParticle &) = 0;

  FLOAT gamma;
  FLOAT gammam1;

};



// ============================================================================
// Class Isothermal
// Isothermal EOS class defintion
// ============================================================================
class Isothermal: public EOS
{
 public:

  Isothermal(FLOAT, FLOAT, FLOAT);
  ~Isothermal();

  FLOAT Pressure(SphParticle &);
  FLOAT EntropicFunction(SphParticle &);
  FLOAT SoundSpeed(SphParticle &);
  FLOAT Temperature(SphParticle &);
  FLOAT SpecificInternalEnergy(SphParticle &);

  FLOAT temp0;
  FLOAT mu_bar;

};



// ============================================================================
// Class Adiabatic
// Adiabatic equation of state class definition
// ============================================================================
class Adiabatic: public EOS
{
 public:

  Adiabatic(FLOAT, FLOAT, FLOAT);
  ~Adiabatic();

  FLOAT Pressure(SphParticle &);
  FLOAT EntropicFunction(SphParticle &);
  FLOAT SoundSpeed(SphParticle &);
  FLOAT Temperature(SphParticle &);
  FLOAT SpecificInternalEnergy(SphParticle &);

  FLOAT temp0;
  FLOAT mu_bar;

};


#endif
