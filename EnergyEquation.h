// ============================================================================
// EnergyEquation.h
// ============================================================================


#ifndef _ENERGY_EQUATION_H_
#define _ENERGY_EQUATION_H_

#include "Constants.h"
#include "Dimensions.h"
#include "Sph.h"
#include "EOS.h"
#include "SphParticle.h"



// ============================================================================
// EnergyEquation
// ============================================================================
class EnergyEquation
{
 public:

  EnergyEquation(double);
  ~EnergyEquation();

  virtual void EnergyIntegration(int,SphParticle *,double) = 0;
  virtual void EnergyCorrectionTerms(int,SphParticle *,double) = 0;
  virtual void EndTimestep(int,int,SphParticle *) = 0;
  virtual double Timestep(SphParticle &) = 0;

  const double energy_mult;

};



// ============================================================================
// EnergyPEC
// ============================================================================
class EnergyPEC: public EnergyEquation
{
 public:

  EnergyPEC(double);
  ~EnergyPEC();

  void EnergyIntegration(int,SphParticle *,double);
  void EnergyCorrectionTerms(int,SphParticle *,double);
  void EndTimestep(int,int,SphParticle *);
  double Timestep(SphParticle &);

};


#endif
