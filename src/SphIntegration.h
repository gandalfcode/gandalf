// ============================================================================
// SphIntegration.h
// ============================================================================


#ifndef _SPH_INTEGRATOR_H_
#define _SPH_INTEGRATOR_H_


#include "Precision.h"
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

  //SphIntegration(DOUBLE accel_mult_aux, DOUBLE courant_mult_aux):
  //  accel_mult(accel_mult_aux),courant_mult(courant_mult_aux) {}
  SphIntegration(int, int, DOUBLE, DOUBLE);
  ~SphIntegration();

  virtual void AdvanceParticles(int,int,int,SphParticle *,DOUBLE) = 0;
  virtual void CorrectionTerms(int,int,int,SphParticle *,DOUBLE) = 0;
  virtual void EndTimestep(int,int,int,SphParticle *) = 0;

  virtual DOUBLE Timestep(SphParticle &, int);
  
  int level_step;
  const DOUBLE courant_mult;
  const DOUBLE accel_mult;
#if !defined(FIXED_DIMENSIONS)
  const int ndim;
  const int vdim;
#endif

};



// ============================================================================
// Class SphLeapfrogKDK
// ============================================================================
class SphLFKDK: public SphIntegration
{
 public:

  SphLFKDK(int, int, DOUBLE, DOUBLE);
  ~SphLFKDK();

  void AdvanceParticles(int,int,int,SphParticle *,DOUBLE);
  void CorrectionTerms(int,int,int,SphParticle *,DOUBLE);
  void EndTimestep(int,int,int,SphParticle *);

};


#endif
