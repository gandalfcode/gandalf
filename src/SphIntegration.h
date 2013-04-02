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
template <int ndim>
class SphIntegration
{
 public:

  //SphIntegration(DOUBLE accel_mult_aux, DOUBLE courant_mult_aux):
  //  accel_mult(accel_mult_aux),courant_mult(courant_mult_aux) {}
  SphIntegration(DOUBLE, DOUBLE);
  ~SphIntegration();

  virtual void AdvanceParticles(int,int,int,SphParticle<ndim> *,FLOAT) = 0;
  virtual void CorrectionTerms(int,int,int,SphParticle<ndim> *,FLOAT) = 0;
  virtual void EndTimestep(int,int,int,SphParticle<ndim> *) = 0;

  virtual DOUBLE Timestep(SphParticle<ndim> &, int);
  
  int level_step;
  const DOUBLE courant_mult;
  const DOUBLE accel_mult;
//#if !defined(FIXED_DIMENSIONS)
//  const int ndim;
  static const int vdim=ndim;
//#endif

};



// ============================================================================
// Class SphLeapfrogKDK
// ============================================================================
template <int ndim>
class SphLeapfrogKDK: public SphIntegration<ndim>
{
 public:

  SphLeapfrogKDK(DOUBLE, DOUBLE);
  ~SphLeapfrogKDK();

  void AdvanceParticles(int,int,int,SphParticle<ndim> *,FLOAT);
  void CorrectionTerms(int,int,int,SphParticle<ndim> *,FLOAT);
  void EndTimestep(int,int,int,SphParticle<ndim> *);

};



// ============================================================================
// Class SphGodunovIntegration
// ============================================================================
template <int ndim>
class SphGodunovIntegration: public SphIntegration<ndim>
{
 public:

  SphGodunovIntegration(DOUBLE, DOUBLE);
  ~SphGodunovIntegration();

  void AdvanceParticles(int,int,int,SphParticle<ndim> *,FLOAT);
  void CorrectionTerms(int,int,int,SphParticle<ndim> *,FLOAT);
  void EndTimestep(int,int,int,SphParticle<ndim> *);
  static const int vdim = ndim;
  DOUBLE Timestep(SphParticle<ndim> &, int);
};


#endif
