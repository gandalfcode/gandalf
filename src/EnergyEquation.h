// ============================================================================
// EnergyEquation.h
// Class definitions of main energy equation class plus inherited children 
// classes for various energy integration algorithms.
// ============================================================================


#ifndef _ENERGY_EQUATION_H_
#define _ENERGY_EQUATION_H_

#include "Constants.h"
#include "Precision.h"
#include "Sph.h"
#include "EOS.h"
#include "SphParticle.h"



// ============================================================================
// EnergyEquation
// Main energy equation class, with virtual functions that require full 
// definitions in the children classes.
// ============================================================================
template <int ndim>
class EnergyEquation
{
 public:

  EnergyEquation(DOUBLE);
  ~EnergyEquation();

  virtual void EnergyIntegration(int, int, int, SphParticle<ndim> *, FLOAT) = 0;
  virtual void EnergyCorrectionTerms(int, int, int, SphParticle<ndim> *, FLOAT) = 0;
  virtual void EndTimestep(int, int, int, SphParticle<ndim> *) = 0;
  virtual DOUBLE Timestep(SphParticle<ndim> &) = 0;

  const DOUBLE energy_mult;

};



// ============================================================================
// EnergyPEC
// Class definition for energy equation integration class using a 
// Predict-Evaluate-Correct (PEC) scheme.
// ============================================================================
template <int ndim>
class EnergyPEC: public EnergyEquation<ndim>
{
 public:

  EnergyPEC(DOUBLE);
  ~EnergyPEC();

  void EnergyIntegration(int, int, int, SphParticle<ndim> *, FLOAT);
  void EnergyCorrectionTerms(int, int, int, SphParticle<ndim> *, FLOAT);
  void EndTimestep(int, int, int, SphParticle<ndim> *);
  DOUBLE Timestep(SphParticle<ndim> &);

};



// ============================================================================
// EnergyGodunovIntegration
// Class definition for energy equation integration class using a 
// Predict-Evaluate-Correct (PEC) scheme.
// ============================================================================
template <int ndim>
class EnergyGodunovIntegration: public EnergyEquation<ndim>
{
 public:

  EnergyGodunovIntegration(DOUBLE);
  ~EnergyGodunovIntegration();

  void EnergyIntegration(int, int, int, SphParticle<ndim> *, FLOAT);
  void EnergyCorrectionTerms(int, int, int, SphParticle<ndim> *, FLOAT);
  void EndTimestep(int, int, int, SphParticle<ndim> *);
  DOUBLE Timestep(SphParticle<ndim> &);

};



#endif
