//=============================================================================
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
//=============================================================================


#ifndef _SPH_INTEGRATION_H_
#define _SPH_INTEGRATION_H_


#include "Precision.h"
#include "CodeTiming.h"
#include "Constants.h"
#include "Sph.h"
#include "EOS.h"
#include "Parameters.h"
#include "SphParticle.h"



//=============================================================================
//  Class SphIntegration
/// \brief   Main parent Sph Integration class
/// \details Main parent Sph Integration class.  All employed child classes 
///          (e.g. SphLeapfrogKDK) inherit from this class.
/// \author  D. A. Hubber, G. Rosotti
/// \date    03/04/2013
//=============================================================================
template <int ndim>
class SphIntegration
{
 public:

  //SphIntegration(DOUBLE accel_mult_aux, DOUBLE courant_mult_aux):
  //  accel_mult(accel_mult_aux),courant_mult(courant_mult_aux) {}
  SphIntegration(DOUBLE, DOUBLE, DOUBLE, eosenum, tdaviscenum);
  ~SphIntegration();

  virtual void AdvanceParticles(int, FLOAT, Sph<ndim> *) = 0;
  virtual void CorrectionTerms(int, FLOAT, Sph<ndim> *) = 0;
  virtual void EndTimestep(int, FLOAT, Sph<ndim> *) = 0;
  virtual int CheckTimesteps(int, int, Sph<ndim> *) = 0;
  virtual DOUBLE Timestep(SphParticle<ndim> &, Sph<ndim> *);
  
  const DOUBLE accel_mult;
  const DOUBLE courant_mult;
  const DOUBLE energy_mult;
  const eosenum gas_eos;
  const tdaviscenum tdavisc;
  static const int vdim=ndim;

  CodeTiming *timing;               ///< Pointer to code timing object

};



//=============================================================================
//  Class SphLeapfrogKDK
/// \brief   Leapfrog kick-drift-kick SPH particle integration scheme.
/// \details Class definition for leapfrog kick-drift-kick SPH particle 
///          integration scheme.  Inherits from main parent SphIntegration 
///          class and provides implementations of all virtual functions.
/// \author  D. A. Hubber, G. Rosotti
/// \date    03/04/2013
//=============================================================================
template <int ndim>
class SphLeapfrogKDK: public SphIntegration<ndim>
{
 public:

  using SphIntegration<ndim>::gas_eos;
  using SphIntegration<ndim>::tdavisc;
  using SphIntegration<ndim>::timing;

  SphLeapfrogKDK(DOUBLE, DOUBLE, DOUBLE, eosenum, tdaviscenum);
  ~SphLeapfrogKDK();

  void AdvanceParticles(int, FLOAT, Sph<ndim> *);
  void CorrectionTerms(int, FLOAT, Sph<ndim> *);
  void EndTimestep(int, FLOAT, Sph<ndim> *);
  int CheckTimesteps(int, int, Sph<ndim> *);

};



//=============================================================================
//  Class SphLeapfrogDKD
/// \brief   Leapfrog drift-kick-drift SPH particle integration scheme.
/// \details Class definition for leapfrog drift-kick-drift SPH particle 
///          integration scheme.  Inherits from main parent SphIntegration 
///          class and provides implementations of all virtual functions.
/// \author  D. A. Hubber, G. Rosotti
/// \date    03/04/2013
//=============================================================================
template <int ndim>
class SphLeapfrogDKD: public SphIntegration<ndim>
{
 public:

  using SphIntegration<ndim>::gas_eos;
  using SphIntegration<ndim>::tdavisc;
  using SphIntegration<ndim>::timing;

  SphLeapfrogDKD(DOUBLE, DOUBLE, DOUBLE, eosenum, tdaviscenum);
  ~SphLeapfrogDKD();

  void AdvanceParticles(int, FLOAT, Sph<ndim> *);
  void CorrectionTerms(int, FLOAT, Sph<ndim> *);
  void EndTimestep(int, FLOAT, Sph<ndim> *);
  int CheckTimesteps(int, int, Sph<ndim> *);

};



//=============================================================================
//  Class SphGodunovIntegration
/// \brief   Inutsuka (2002) Godunov SPH conservative SPH integration scheme.
/// \details Class definition for Inutsuka (2002) Godunov SPH conservative 
///          SPH integration scheme.  Algorithm conserves energy to 
///          machine precision (for direct summation and global timesteps).
/// \author  D. A. Hubber, G. Rosotti
/// \date    03/04/2013
//=============================================================================
template <int ndim>
class SphGodunovIntegration: public SphIntegration<ndim>
{
 public:

  using SphIntegration<ndim>::gas_eos;
  using SphIntegration<ndim>::tdavisc;
  using SphIntegration<ndim>::timing;

  SphGodunovIntegration(DOUBLE, DOUBLE, DOUBLE, eosenum, tdaviscenum);
  ~SphGodunovIntegration();

  void AdvanceParticles(int, FLOAT, Sph<ndim> *);
  void CorrectionTerms(int, FLOAT, Sph<ndim> *);
  void EndTimestep(int, FLOAT, Sph<ndim> *);
  int CheckTimesteps(int, int, Sph<ndim> *);
  DOUBLE Timestep(SphParticle<ndim> &, Sph<ndim> *);

  static const int vdim = ndim;

};
#endif
