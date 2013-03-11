// ============================================================================
// NbodyIntegration.h
// ============================================================================


#ifndef _NBODY_INTEGRATOR_H_
#define _NBODY_INTEGRATOR_H_


#include "Precision.h"
#include "Constants.h"
#include "Dimensions.h"
#include "Parameters.h"
#include "StarParticle.h"



// ============================================================================
// Class NbodyIntegration
// ============================================================================
class NbodyIntegration
{
 public:

  NbodyIntegration(int, int, DOUBLE);
  ~NbodyIntegration();

  virtual void AdvanceParticles(int,int,int,StarParticle *,DOUBLE) = 0;
  virtual void CorrectionTerms(int,int,int,StarParticle *,DOUBLE) = 0;
  virtual void EndTimestep(int,int,int,SphParticle *) = 0;

  virtual DOUBLE Timestep(SphParticle &, int);
  
  int level_step;
  const DOUBLE nbody_mult;
#if !defined(FIXED_DIMENSIONS)
  const int ndim;
  const int vdim;
#endif

};



// ============================================================================
// Class NbodyLeapfrogKDK
// ============================================================================
class NbodyLeapfrogKDK: public NbodyIntegration
{
 public:

  NbodyLeapfrogKDK(int, int, DOUBLE);
  ~NbodyLeapfrogKDK();

  void AdvanceParticles(int,int,int,SphParticle *,FLOAT);
  void CorrectionTerms(int,int,int,SphParticle *,FLOAT);
  void EndTimestep(int,int,int,SphParticle *);

};



// ============================================================================
// Class NbodyHermite4
// ============================================================================
class NbodyHermite4: public NbodyIntegration
{
 public:

  NbodyHermite4(int, int, DOUBLE);
  ~NbodyHermite4();

  void AdvanceParticles(int,int,int,SphParticle *,FLOAT);
  void CorrectionTerms(int,int,int,SphParticle *,FLOAT);
  void EndTimestep(int,int,int,SphParticle *);

};



#endif
