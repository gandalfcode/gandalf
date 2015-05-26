//=================================================================================================
//  EnergyEquation.h
//  Class definitions of main energy equation class plus inherited children
//  classes for various energy integration algorithms.
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


#ifndef _ENERGY_EQUATION_H_
#define _ENERGY_EQUATION_H_


#include "Constants.h"
#include "Precision.h"
#include "Sph.h"
#include "EOS.h"
#include "Particle.h"
#include "SimUnits.h"


//=================================================================================================
//  EnergyEquation
/// Main energy equation class, with virtual functions that require full
/// definitions in the children classes.
//=================================================================================================
template <int ndim>
class EnergyEquation
{
 public:

  EnergyEquation(DOUBLE);
  ~EnergyEquation();

  virtual void EnergyIntegration(const int, const int, const FLOAT, const FLOAT,
                                SphParticle<ndim> *) = 0;
  virtual void EnergyCorrectionTerms(const int, const int, const FLOAT, const FLOAT,
                                     SphParticle<ndim> *) = 0;
  virtual void EndTimestep(const int, const int, const FLOAT, const FLOAT, SphParticle<ndim> *) = 0;
  virtual DOUBLE Timestep(SphParticle<ndim> &) = 0;


  const DOUBLE energy_mult;            ///< ..
  CodeTiming *timing;                  ///< Pointer to code timing object

};



//=================================================================================================
//  EnergyPEC
/// Class definition for energy equation integration class using a
/// Predict-Evaluate-Correct (PEC) scheme.
//=================================================================================================
template <int ndim, template <int> class ParticleType>
class EnergyPEC: public EnergyEquation<ndim>
{
 public:

  using EnergyEquation<ndim>::timing;

  EnergyPEC(DOUBLE);
  ~EnergyPEC();

  void EnergyIntegration(const int, const int, const FLOAT, const FLOAT, SphParticle<ndim> *);
  void EnergyCorrectionTerms(const int, const int, const FLOAT, const FLOAT, SphParticle<ndim> *);
  void EndTimestep(const int, const int, const FLOAT, const FLOAT, SphParticle<ndim> *);
  DOUBLE Timestep(SphParticle<ndim> &);

};




//=================================================================================================
//  EnergyRadws
/// ..
//=================================================================================================
template <int ndim, template <int> class ParticleType>
class EnergyRadws: public EnergyEquation<ndim>
{
 public:

   using EnergyEquation<ndim>::timing;

  EnergyRadws(DOUBLE, string, FLOAT, SimUnits *);
  ~EnergyRadws();

  //  void ReadTable();
  void EnergyIntegration(const int, const int, const FLOAT, const FLOAT, SphParticle<ndim> *);
  void EnergyCorrectionTerms(const int, const int, const FLOAT, const FLOAT, SphParticle<ndim> *) {};
  void EndTimestep(const int, const int, const FLOAT, const FLOAT, SphParticle<ndim> *);
  void EnergyFindEqui(FLOAT , FLOAT , FLOAT , FLOAT , FLOAT, FLOAT &, FLOAT &);
  void EnergyFindEquiTemp(int, FLOAT, FLOAT, FLOAT, FLOAT, FLOAT &);

  int GetIDens(FLOAT);
  int GetITemp(FLOAT);
  void GetKappa(int, int, FLOAT, FLOAT, FLOAT &, FLOAT &, FLOAT &);
  FLOAT GetEnergy(int , int , FLOAT , FLOAT );
  DOUBLE Timestep(SphParticle<ndim> &) {return big_number_dp;}
  FLOAT ebalance(FLOAT, FLOAT, FLOAT, FLOAT, FLOAT, FLOAT);


  int ndens, ntemp;
  FLOAT  *eos_dens, *eos_temp ;
  FLOAT **eos_energy, **eos_mu, **kappa_table, **kappar_table, **kappap_table;
  FLOAT rad_const;
  FLOAT temp_ambient;

};



//=================================================================================================
//  NullEnergy
/// Null (empty) class when no energy option is selected.
//=================================================================================================
template <int ndim>
class NullEnergy: public EnergyEquation<ndim>
{
 public:

  NullEnergy(DOUBLE dt_mult) : EnergyEquation<ndim>(dt_mult) {};
  ~NullEnergy();

  void EnergyIntegration(const int, const int, const FLOAT, const FLOAT, SphParticle<ndim> *) {};
  void EnergyCorrectionTerms(const int, const int, const FLOAT, const FLOAT, SphParticle<ndim> *) {};
  void EndTimestep(const int, const int, const FLOAT, const FLOAT, SphParticle<ndim> *) {};
  DOUBLE Timestep(SphParticle<ndim> &) {return big_number_dp;}

};
#endif
