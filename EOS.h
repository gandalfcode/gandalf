// ============================================================================
// EOS.h
// ============================================================================


#ifndef _EOS_H_
#define _EOS_H_


#include "Constants.h"
#include "Dimensions.h"
#include "Parameters.h"
#include "SphParticle.h"



// ============================================================================
// CLASS EOS
// ============================================================================
class EOS
{
 public:

  virtual float Pressure(SphParticle &) = 0;
  virtual float EntropicFunction(SphParticle &) = 0;
  virtual float SoundSpeed(SphParticle &) = 0;
  virtual float Temperature(SphParticle &) = 0;
  virtual float SpecificInternalEnergy(SphParticle &) = 0;

};



// ============================================================================
// CLASS EOS
// ============================================================================
class Isothermal: public EOS
{
 public:

  Isothermal(float,float,float);
  ~Isothermal();

  float Pressure(SphParticle &);
  float EntropicFunction(SphParticle &);
  float SoundSpeed(SphParticle &);
  float Temperature(SphParticle &);
  float SpecificInternalEnergy(SphParticle &);

  float gamma;
  float gammam1;
  float temp0;
  float mu_bar;

};



// ============================================================================
// Class Adiabatic
// ============================================================================
class Adiabatic: public EOS
{
 public:

  Adiabatic(float,float,float);
  ~Adiabatic();

  float Pressure(SphParticle &);
  float EntropicFunction(SphParticle &);
  float SoundSpeed(SphParticle &);
  float Temperature(SphParticle &);
  float SpecificInternalEnergy(SphParticle &);

  float gamma;
  float gammam1;
  float temp0;
  float mu_bar;

};


#endif
