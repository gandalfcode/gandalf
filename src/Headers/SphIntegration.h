//=================================================================================================
//  SphIntegration.h
//  Contains class definitions for all SPH integration schemes.
//
//  This file is part of GANDALF :
//  Graphical Astrophysics code for N-body Dynamics And Lagrangian Fluids
//  https://github.com/gandalfcode/gandalf
//  Contact : gandalfcode@gmail.com
//
//  Copyright (C) 2013  D. A. Hubber, G. Rosotti
//
//  GANDALF is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 2 of the License, or
//  (at your option) any later version.
//
//  GANDALF is distributed in the hope that it will be useful, but
//  WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  General Public License (http://www.gnu.org/licenses) for more details.
//=================================================================================================


#ifndef _SPH_INTEGRATION_H_
#define _SPH_INTEGRATION_H_


#include "Precision.h"
#include "CodeTiming.h"
#include "Constants.h"
#include "DomainBox.h"
#include "Sph.h"
#include "EOS.h"
#include "Parameters.h"
#include "Particle.h"



//=================================================================================================
//  Class SphIntegration
/// \brief   Main parent Sph Integration class
/// \details Main parent Sph Integration class.  All employed child classes
///          (e.g. SphLeapfrogKDK) inherit from this class.
/// \author  D. A. Hubber, G. Rosotti
/// \date    03/04/2013
//=================================================================================================
template <int ndim>
class SphIntegration
{
 public:

  //SphIntegration(DOUBLE accel_mult_aux, DOUBLE courant_mult_aux):
  //  accel_mult(accel_mult_aux),courant_mult(courant_mult_aux) {}
  SphIntegration(DOUBLE, DOUBLE, DOUBLE, eosenum, tdaviscenum);
  ~SphIntegration();

  virtual void AdvanceParticles(const unsigned int, const int, const FLOAT,
                                const FLOAT, SphParticle<ndim> *) = 0;
  virtual void CorrectionTerms(const unsigned int, const int, const FLOAT,
                               const FLOAT, SphParticle<ndim> *) = 0;
  virtual void EndTimestep(const unsigned int, const int, const FLOAT,
                           const FLOAT, SphParticle<ndim> *) = 0;
  virtual int CheckTimesteps(const unsigned int, const unsigned int, const unsigned int,
                             const int, SphParticle<ndim> *) = 0;
  virtual DOUBLE Timestep(SphParticle<ndim> &, Sph<ndim> *);
  virtual void CheckBoundaries(DomainBox<ndim> &, Sph<ndim> *);

  const DOUBLE accel_mult;
  const DOUBLE courant_mult;
  const DOUBLE energy_mult;
  const eosenum gas_eos;
  const tdaviscenum tdavisc;
  static const int vdim=ndim;

  CodeTiming *timing;               ///< Pointer to code timing object

};



//=================================================================================================
//  Class SphLeapfrogKDK
/// \brief   Leapfrog kick-drift-kick SPH particle integration scheme.
/// \details Class definition for leapfrog kick-drift-kick SPH particle
///          integration scheme.  Inherits from main parent SphIntegration
///          class and provides implementations of all virtual functions.
/// \author  D. A. Hubber, G. Rosotti
/// \date    03/04/2013
//=================================================================================================
template <int ndim, template <int> class ParticleType>
class SphLeapfrogKDK: public SphIntegration<ndim>
{
 public:

  using SphIntegration<ndim>::gas_eos;
  using SphIntegration<ndim>::tdavisc;
  using SphIntegration<ndim>::timing;

  SphLeapfrogKDK(DOUBLE, DOUBLE, DOUBLE, eosenum, tdaviscenum);
  ~SphLeapfrogKDK();

  void AdvanceParticles(const unsigned int, const int, const FLOAT,
                        const FLOAT, SphParticle<ndim> *);
  void CorrectionTerms(const unsigned int, const int, const FLOAT, const FLOAT, SphParticle<ndim> *);
  void EndTimestep(const unsigned int, const int, const FLOAT, const FLOAT, SphParticle<ndim> *);
  int CheckTimesteps(const unsigned int, const unsigned int, const unsigned int,
                     const int, SphParticle<ndim> *);

};



//=================================================================================================
//  Class SphLeapfrogDKD
/// \brief   Leapfrog drift-kick-drift SPH particle integration scheme.
/// \details Class definition for leapfrog drift-kick-drift SPH particle
///          integration scheme.  Inherits from main parent SphIntegration
///          class and provides implementations of all virtual functions.
/// \author  D. A. Hubber, G. Rosotti
/// \date    03/04/2013
//=================================================================================================
template <int ndim, template <int> class ParticleType>
class SphLeapfrogDKD: public SphIntegration<ndim>
{
 public:

  using SphIntegration<ndim>::gas_eos;
  using SphIntegration<ndim>::tdavisc;
  using SphIntegration<ndim>::timing;

  SphLeapfrogDKD(DOUBLE, DOUBLE, DOUBLE, eosenum, tdaviscenum);
  ~SphLeapfrogDKD();

  void AdvanceParticles(const unsigned int, const int, const FLOAT, const FLOAT, SphParticle<ndim> *);
  void CorrectionTerms(const unsigned int, const int, const FLOAT, const FLOAT, SphParticle<ndim> *) {};
  void EndTimestep(const unsigned int, const int, const FLOAT, const FLOAT, SphParticle<ndim> *);
  int CheckTimesteps(const unsigned int, const unsigned int, const unsigned int,
                     const int, SphParticle<ndim> *);

};
#endif
