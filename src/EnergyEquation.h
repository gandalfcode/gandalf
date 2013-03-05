// ============================================================================
// EnergyEquation.h
// Class definitions of main energy equation class plus inherited children 
// classes for various energy integration algorithms.
// ============================================================================


#ifndef _ENERGY_EQUATION_H_
#define _ENERGY_EQUATION_H_

#include "Constants.h"
#include "Dimensions.h"
#include "Precision.h"
#include "Sph.h"
#include "EOS.h"
#include "SphParticle.h"



// ============================================================================
// EnergyEquation
// Main energy equation class, with virtual functions that require full 
// definitions in the children classes.
// ============================================================================
class EnergyEquation
{
 public:

  EnergyEquation(DOUBLE);
  ~EnergyEquation();

  virtual void EnergyIntegration(int, int, int, SphParticle *, FLOAT) = 0;
  virtual void EnergyCorrectionTerms(int, int, int, SphParticle *, FLOAT) = 0;
  virtual void EndTimestep(int, int, int, SphParticle *) = 0;
  virtual DOUBLE Timestep(SphParticle &) = 0;

  const DOUBLE energy_mult;

};



// ============================================================================
// EnergyPEC
// Class definition for energy equation integration class using a 
// Predict-Evaluate-Correct (PEC) scheme.
// ============================================================================
class EnergyPEC: public EnergyEquation
{
 public:

  EnergyPEC(DOUBLE);
  ~EnergyPEC();

  void EnergyIntegration(int, int, int, SphParticle *, FLOAT);
  void EnergyCorrectionTerms(int, int, int, SphParticle *, FLOAT);
  void EndTimestep(int, int, int, SphParticle *);
  DOUBLE Timestep(SphParticle &);

};


#endif
