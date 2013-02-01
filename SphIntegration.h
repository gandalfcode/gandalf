// ============================================================================
// SphIntegration.h
// ============================================================================


#ifndef _SPH_INTEGRATOR_H_
#define _SPH_INTEGRATOR_H_


#include "Constants.h"
#include "Dimensions.h"
#include "Sph.h"
#include "EOS.h"
#include "Parameters.h"
#include "SphParticle.h"



// ============================================================================
// Class SphIntegration
// ============================================================================
class SphIntegration
{
 public:

  //SphIntegration(double accel_mult_aux, double courant_mult_aux):
  //  accel_mult(accel_mult_aux),courant_mult(courant_mult_aux) {}
  SphIntegration(double, double);
  ~SphIntegration();

  virtual void AdvanceParticles(int,SphParticle *,double) = 0;
  virtual void CorrectionTerms(int,SphParticle *,double) = 0;
  virtual void EndTimestep(int,int,SphParticle *) = 0;

  virtual double Timestep(SphParticle &, Parameters &);
  
  const double courant_mult;
  const double accel_mult;

};



// ============================================================================
// Class SphLeapfrogKDK
// ============================================================================
class SphLFKDK: public SphIntegration
{
 public:

  SphLFKDK(double, double);
  ~SphLFKDK();

  void AdvanceParticles(int,SphParticle *,double);
  void CorrectionTerms(int,SphParticle *,double);
  void EndTimestep(int,int,SphParticle *);


};


#endif
