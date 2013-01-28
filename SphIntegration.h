// ============================================================================
// SphIntegration.h
// ============================================================================


#ifndef _SPH_INTEGRATOR_H_
#define _SPH_INTEGRATOR_H_

#include "Constants.h"
#include "Dimensions.h"
#include "Sph.h"
#include "EOS.h"
#include "SphParticle.h"



// ============================================================================
// CLASS SphIntegration
// ============================================================================
class SphIntegration
{
 public:

  SphIntegration(double, double);
  ~SphIntegration();

  virtual void AdvanceParticles(int,SphParticle *,double) = 0;
  virtual void CorrectionTerms(int,SphParticle *,double) = 0;

  virtual double Timestep(SphParticle &);
  
  virtual void EndTimestep(int,int,SphParticle *) = 0;

  const double courant_mult;
  const double accel_mult;

};



// ============================================================================
// CLASS SphLFKDK
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
