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


#ifndef _INTEGRATION_H_
#define _INTEGRATION_H_


#include "Precision.h"
#include "CodeTiming.h"
#include "Constants.h"
#include "DomainBox.h"
#include "Sph.h"
#include "EOS.h"
#include "Parameters.h"
#include "Particle.h"



//=================================================================================================
//  Class TimeIntegration
/// \brief   Main parent Sph Integration class
/// \details Main parent Sph Integration class.  All employed child classes
///          (e.g. SphLeapfrogKDK) inherit from this class.
/// \author  D. A. Hubber, G. Rosotti
/// \date    03/04/2013
//=================================================================================================
template <int ndim>
class TimeIntegration
{
 public:

  TimeIntegration() : timing(NULL) { } ;
  virtual ~TimeIntegration(){};

  virtual void SetActiveParticles(const int, Hydrodynamics<ndim> *) = 0 ;
  virtual void AdvanceParticles(const int, const FLOAT, const FLOAT, Hydrodynamics<ndim> *) = 0;
  virtual void CorrectionTerms(const int, const FLOAT, const FLOAT,Hydrodynamics<ndim> *) = 0;
  virtual void EndTimestep(const int, const FLOAT, const FLOAT, Hydrodynamics<ndim> *) = 0;
  virtual int CheckTimesteps(const int, const int, const int, const FLOAT,
                             Hydrodynamics<ndim> *) = 0;
  virtual DOUBLE Timestep(Particle<ndim> &, Hydrodynamics<ndim> *) = 0;
  virtual void CheckBoundaries(DomainBox<ndim> &, Hydrodynamics<ndim> *);

  static const int vdim=ndim;
  CodeTiming *timing;               ///< Pointer to code timing object
};


template <int ndim>
class SphIntegration : public TimeIntegration<ndim>
{

  using TimeIntegration<ndim>::vdim;
 public:
  using TimeIntegration<ndim>::timing;


  //SphIntegration(DOUBLE accel_mult_aux, DOUBLE courant_mult_aux):
  //  accel_mult(accel_mult_aux),courant_mult(courant_mult_aux) {}
  SphIntegration(DOUBLE, DOUBLE, DOUBLE, eosenum, tdaviscenum);
  virtual ~SphIntegration();

  void SetActiveParticles(const int, Hydrodynamics<ndim> *) = 0 ;
  void AdvanceParticles(const int, const FLOAT, const FLOAT, Hydrodynamics<ndim> *) = 0;
  void CorrectionTerms(const int, const FLOAT, const FLOAT,Hydrodynamics<ndim> *) = 0;
  void EndTimestep(const int, const FLOAT, const FLOAT, Hydrodynamics<ndim> *) = 0;
  int CheckTimesteps(const int, const int, const int, const FLOAT, Hydrodynamics<ndim> *) = 0;
  DOUBLE Timestep(Particle<ndim> &, Hydrodynamics<ndim> *);


  const DOUBLE accel_mult ;
  const DOUBLE courant_mult ;
  const DOUBLE energy_mult ;
  const eosenum gas_eos ;
  const tdaviscenum tdavisc;
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
  virtual ~SphLeapfrogKDK();

  void SetActiveParticles(const int, Hydrodynamics<ndim> *);
  void AdvanceParticles(const int, const FLOAT, const FLOAT, Hydrodynamics<ndim> *);
  void CorrectionTerms(const int, const FLOAT, const FLOAT, Hydrodynamics<ndim> *);
  void EndTimestep(const int, const FLOAT, const FLOAT, Hydrodynamics<ndim> *);
  int CheckTimesteps(const int, const int, const int, const FLOAT, Hydrodynamics<ndim> *);

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
  virtual  ~SphLeapfrogDKD();

  void SetActiveParticles(const int, Hydrodynamics<ndim> *);
  void AdvanceParticles(const int, const FLOAT, const FLOAT, Hydrodynamics<ndim> *);
  void CorrectionTerms(const int, const FLOAT, const FLOAT, Hydrodynamics<ndim> *){};
  void EndTimestep(const int, const FLOAT, const FLOAT, Hydrodynamics<ndim> *);
  int CheckTimesteps(const int, const int, const int, const FLOAT, Hydrodynamics<ndim> *);


};

//=================================================================================================
//  Class  MfvIntegration
/// \brief  Time-integration for meshless scheme
/// \details Technically only valid for MfvMuscl scheme.
/// \author  R. Booth
/// \date    27/03/2017
//=================================================================================================
template <int ndim, template <int> class ParticleType>
class MfvIntegration: public TimeIntegration<ndim>
{
public:
  using TimeIntegration<ndim>::timing;


  MfvIntegration(Parameters *simparams)
  : accel_mult(simparams->floatparams["accel_mult"]),
    courant_mult(simparams->floatparams["courant_mult"]),
    visc_mult(simparams->floatparams["visc_mult"]),
    visc_coeff(simparams->floatparams["shear_visc"] + simparams->floatparams["bulk_visc"]),
    staticParticles(simparams->intparams["static_particles"])
  { } ;
  virtual ~MfvIntegration() {} ;

  void SetActiveParticles(const int, Hydrodynamics<ndim> *);
  void AdvanceParticles(const int, const FLOAT, const FLOAT, Hydrodynamics<ndim> *);
  void CorrectionTerms(const int, const FLOAT, const FLOAT, Hydrodynamics<ndim> *){};
  void EndTimestep(const int, const FLOAT, const FLOAT, Hydrodynamics<ndim> *);
  int CheckTimesteps(const int, const int, const int, const FLOAT, Hydrodynamics<ndim> *);
  DOUBLE Timestep(Particle<ndim> &, Hydrodynamics<ndim> *);

 private:
  const DOUBLE accel_mult ;
  const DOUBLE courant_mult ;
  const DOUBLE visc_mult ;
  const DOUBLE visc_coeff;
  bool staticParticles;

};




#endif
