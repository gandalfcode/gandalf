//=================================================================================================
//  EOS.h
//  Contains (virtual) definitions for equation of state class.
//  Also includes inherited class defintions for various equation of
//  state options.
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


#ifndef _EOS_H_
#define _EOS_H_


#include "Precision.h"
#include "Constants.h"
#include "Nbody.h"
#include "Parameters.h"
#include "Particle.h"
//#include "SphNeighbourSearch.h"
#include "SimUnits.h"


// Forward declaration of Hydrodynamics to break circular dependency
template <int ndim>
class Hydrodynamics;

template <int ndim>
class SphNeighbourSearch;

template <int ndim>
class EOS;


enum eosenum{noeos, isothermal, barotropic, barotropic2, energy_eqn,
             constant_temp, radws, Nhydroeos};



//=================================================================================================
//  Class EOS
/// \brief   Main equation of state
/// \details Main equation of state class.  Only contains virtual function
///          defintions. All functions must be defined by the inherited EOS classes.
/// \author  D. A. Hubber, G. Rosotti
/// \date    03/04/2013
//=================================================================================================
template <int ndim>
class EOS
{
 public:

  EOS(FLOAT _gamma):
    gamma(_gamma),
    gammam1(_gamma - (FLOAT) 1.0),
    gammaMinusOne(_gamma - (FLOAT) 1.0),
    oneMinusGamma((FLOAT) 1.0 - _gamma) {};

  template <class ParticleType> FLOAT Pressure(const ParticleType& part) {return gammam1*part.rho*part.u;} ;
  virtual FLOAT EntropicFunction(Particle<ndim> &) = 0;
  virtual FLOAT SoundSpeed(Particle<ndim> &) = 0;
  virtual FLOAT Temperature(Particle<ndim> &) = 0;
  virtual FLOAT SpecificInternalEnergy(Particle<ndim> &) = 0;

  const FLOAT gamma;
  const FLOAT gammam1;
  const FLOAT gammaMinusOne;
  const FLOAT oneMinusGamma;

};



//=================================================================================================
//  Class Isothermal
/// \brief   Isothermal equation of state
/// \details Isothermal equation of state
/// \author  D. A. Hubber, G. Rosotti
/// \date    03/04/2013
//=================================================================================================
template <int ndim>
class Isothermal: public EOS<ndim>
{
  using EOS<ndim>::gamma;
  using EOS<ndim>::gammam1;
  using EOS<ndim>::gammaMinusOne;
  using EOS<ndim>::oneMinusGamma;

 public:

  Isothermal(FLOAT, FLOAT, FLOAT, SimUnits *);
  ~Isothermal();

  FLOAT EntropicFunction(Particle<ndim> &);
  FLOAT SoundSpeed(Particle<ndim> &);
  FLOAT Temperature(Particle<ndim> &);
  FLOAT SpecificInternalEnergy(Particle<ndim> &);

  const FLOAT temp0;
  const FLOAT mu_bar;

};



//=================================================================================================
//  Class Barotropic
/// \brief   Barotropic equation of state
/// \details Barotropic equation of state
/// \author  D. A. Hubber, G. Rosotti
/// \date    03/04/2013
//=================================================================================================
template <int ndim>
class Barotropic: public EOS<ndim>
{
  using EOS<ndim>::gamma;
  using EOS<ndim>::gammam1;

 public:

  Barotropic(FLOAT, FLOAT, FLOAT, FLOAT, SimUnits *);
  ~Barotropic();

  FLOAT EntropicFunction(Particle<ndim> &);
  FLOAT SoundSpeed(Particle<ndim> &);
  FLOAT Temperature(Particle<ndim> &);
  FLOAT SpecificInternalEnergy(Particle<ndim> &);

  FLOAT temp0;
  FLOAT mu_bar;
  FLOAT rho_bary;
  FLOAT invrho_bary;

};



//=================================================================================================
//  Class Barotropic2
/// \brief   2nd form of barotropic equation of state (piecewise polynomials)
/// \details 2nd form of barotropic equation of state (piecewise polynomials)
/// \author  D. A. Hubber, G. Rosotti
/// \date    08/04/2014
//=================================================================================================
template <int ndim>
class Barotropic2: public EOS<ndim>
{
  using EOS<ndim>::gamma;
  using EOS<ndim>::gammam1;

 public:

  Barotropic2(FLOAT, FLOAT, FLOAT, FLOAT, SimUnits *);
  ~Barotropic2();

  FLOAT EntropicFunction(Particle<ndim> &);
  FLOAT SoundSpeed(Particle<ndim> &);
  FLOAT Temperature(Particle<ndim> &);
  FLOAT SpecificInternalEnergy(Particle<ndim> &);

  FLOAT temp0;
  FLOAT mu_bar;
  FLOAT rho_bary;
  FLOAT invrho_bary;

};



//=================================================================================================
//  Class Adiabatic
/// \brief   Adiabatic equation of state
/// \details Adiabatic equation of state.  Requires integrating the energy
///          equation parallel to the main dynamical quantities.
/// \author  D. A. Hubber, G. Rosotti
/// \date    03/04/2013
//=================================================================================================
template <int ndim>
class Adiabatic: public EOS<ndim>
{
  using EOS<ndim>::gamma;
  using EOS<ndim>::gammam1;

 public:

  Adiabatic(FLOAT, FLOAT, FLOAT);
  ~Adiabatic();

  FLOAT EntropicFunction(Particle<ndim> &);
  FLOAT SoundSpeed(Particle<ndim> &);
  FLOAT Temperature(Particle<ndim> &);
  FLOAT SpecificInternalEnergy(Particle<ndim> &);

  FLOAT temp0;
  const FLOAT mu_bar;

};



//=================================================================================================
//  Class IonisingRadiation
/// \brief   Ionising radiation from stars/sinks including general EOS
/// \details Ionising radiation from stars/sinks including general EOS
/// \author  S. Balfour, D. A. Hubber
/// \date    24/04/2014
//=================================================================================================
template <int ndim>
class IonisingRadiation: public EOS<ndim>
{
  using EOS<ndim>::gamma;
  using EOS<ndim>::gammam1;

 public:

  IonisingRadiation(string, FLOAT, FLOAT, FLOAT, FLOAT, SimUnits *);
  ~IonisingRadiation();

  FLOAT EntropicFunction(Particle<ndim> &);
  FLOAT SoundSpeed(Particle<ndim> &);
  FLOAT Temperature(Particle<ndim> &);
  FLOAT SpecificInternalEnergy(Particle<ndim> &);

  FLOAT temp0;
  FLOAT mu_bar;
  EOS<ndim> *eos;

};



//=================================================================================================
//  Class MCRadiationEOS
/// \brief   Ionising radiation from stars/sinks including general EOS
/// \details Ionising radiation from stars/sinks including general EOS
/// \author  S. Balfour, D. A. Hubber
/// \date    24/04/2014
//=================================================================================================
template <int ndim>
class MCRadiationEOS: public EOS<ndim>
{
  using EOS<ndim>::gamma;
  using EOS<ndim>::gammam1;
  using EOS<ndim>::gammaMinusOne;
  using EOS<ndim>::oneMinusGamma;

 public:

  MCRadiationEOS(string, FLOAT, FLOAT, FLOAT, FLOAT, FLOAT, FLOAT, SimUnits *);
  ~MCRadiationEOS();

  FLOAT EntropicFunction(Particle<ndim> &);
  FLOAT SoundSpeed(Particle<ndim> &);
  FLOAT Temperature(Particle<ndim> &);
  FLOAT SpecificInternalEnergy(Particle<ndim> &);

  FLOAT mu_bar;
  FLOAT mu_ion;
  FLOAT temp0;
  FLOAT temp_ion;
  EOS<ndim> *eos;

};



//=================================================================================================
//  Class Radws
/// \brief   Equation of state for Stamatellos wt al. (2007) radiative cooling scheme
/// \details ...
/// \author  D. A. Hubber, G. Rosotti
/// \date    03/04/2013
//=================================================================================================
template <int ndim>
class Radws : public EOS<ndim>
{
  using EOS<ndim>::gamma;
  using EOS<ndim>::gammam1;

 public:

  Radws(FLOAT, FLOAT, FLOAT);
  ~Radws();

  FLOAT EntropicFunction(Particle<ndim> &);
  FLOAT SoundSpeed(Particle<ndim> &);
  FLOAT Temperature(Particle<ndim> &);
  FLOAT SpecificInternalEnergy(Particle<ndim> &);

  FLOAT temp0;
  FLOAT mu_bar;

};
#endif
