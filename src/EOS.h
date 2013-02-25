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

  virtual FLOAT Pressure(SphParticle &) = 0;
  virtual FLOAT EntropicFunction(SphParticle &) = 0;
  virtual FLOAT SoundSpeed(SphParticle &) = 0;
  virtual FLOAT Temperature(SphParticle &) = 0;
  virtual FLOAT SpecificInternalEnergy(SphParticle &) = 0;

};



// ============================================================================
// CLASS EOS
// ============================================================================
class Isothermal: public EOS
{
 public:

  Isothermal(FLOAT,FLOAT,FLOAT);
  ~Isothermal();

  FLOAT Pressure(SphParticle &);
  FLOAT EntropicFunction(SphParticle &);
  FLOAT SoundSpeed(SphParticle &);
  FLOAT Temperature(SphParticle &);
  FLOAT SpecificInternalEnergy(SphParticle &);

  FLOAT gamma;
  FLOAT gammam1;
  FLOAT temp0;
  FLOAT mu_bar;

};



// ============================================================================
// Class Adiabatic
// ============================================================================
class Adiabatic: public EOS
{
 public:

  Adiabatic(FLOAT,FLOAT,FLOAT);
  ~Adiabatic();

  FLOAT Pressure(SphParticle &);
  FLOAT EntropicFunction(SphParticle &);
  FLOAT SoundSpeed(SphParticle &);
  FLOAT Temperature(SphParticle &);
  FLOAT SpecificInternalEnergy(SphParticle &);

  FLOAT gamma;
  FLOAT gammam1;
  FLOAT temp0;
  FLOAT mu_bar;

};


#endif
