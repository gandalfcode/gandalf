//=============================================================================
//  SphIntegration.h
//=============================================================================


#ifndef _SPH_INTEGRATION_H_
#define _SPH_INTEGRATION_H_


#include "Precision.h"
#include "Constants.h"
#include "Sph.h"
#include "EOS.h"
#include "Parameters.h"
#include "SphParticle.h"



//=============================================================================
//  Class SphIntegration
/// \brief   Main parent Sph Integration class
/// \details ..
/// \author  D. A. Hubber, G. Rosotti
/// \date    03/04/2013
//=============================================================================
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
  static const int vdim=ndim;

};



//=============================================================================
//  Class SphLeapfrogKDK
/// \brief   Leapfrog kick-drift-kick SPH particle integration scheme.
/// \details ..
/// \author  D. A. Hubber, G. Rosotti
/// \date    03/04/2013
//=============================================================================
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



//=============================================================================
//  Class SphGodunovIntegration
/// \brief   ..
/// \details ..
/// \author  D. A. Hubber, G. Rosotti
/// \date    03/04/2013
//=============================================================================
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
