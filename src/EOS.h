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
template <int ndim>
class EOS
{
 public:

  EOS (FLOAT _gamma):
    gamma(_gamma),
    gammam1 (gamma-1){};

  virtual FLOAT Pressure(SphParticle<ndim> &) = 0;
  virtual FLOAT EntropicFunction(SphParticle<ndim> &) = 0;
  virtual FLOAT SoundSpeed(SphParticle<ndim> &) = 0;
  virtual FLOAT Temperature(SphParticle<ndim> &) = 0;
  virtual FLOAT SpecificInternalEnergy(SphParticle<ndim> &) = 0;

  FLOAT gamma;
  FLOAT gammam1;

};



// ============================================================================
// Class Isothermal
// Isothermal EOS class defintion
// ============================================================================
template <int ndim>
class Isothermal: public EOS<ndim>
{

  using EOS<ndim>::gamma;
  using EOS<ndim>::gammam1;

 public:

  Isothermal(FLOAT, FLOAT, FLOAT);
  ~Isothermal();

  FLOAT Pressure(SphParticle<ndim> &);
  FLOAT EntropicFunction(SphParticle<ndim> &);
  FLOAT SoundSpeed(SphParticle<ndim> &);
  FLOAT Temperature(SphParticle<ndim> &);
  FLOAT SpecificInternalEnergy(SphParticle<ndim> &);

  FLOAT temp0;
  FLOAT mu_bar;

};



// ============================================================================
// Class Adiabatic
// Adiabatic equation of state class definition
// ============================================================================
template <int ndim>
class Adiabatic: public EOS<ndim>
{

  using EOS<ndim>::gamma;
  using EOS<ndim>::gammam1;

 public:

  Adiabatic(FLOAT, FLOAT, FLOAT);
  ~Adiabatic();

  FLOAT Pressure(SphParticle<ndim> &);
  FLOAT EntropicFunction(SphParticle<ndim> &);
  FLOAT SoundSpeed(SphParticle<ndim> &);
  FLOAT Temperature(SphParticle<ndim> &);
  FLOAT SpecificInternalEnergy(SphParticle<ndim> &);

  FLOAT temp0;
  FLOAT mu_bar;

};


#endif
