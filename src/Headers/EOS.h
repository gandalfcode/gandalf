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
#include "OpacityTable.h"
#include "Particle.h"
#include "SimUnits.h"


// Forward declaration of Hydrodynamics to break circular dependency
template <int ndim>
class Hydrodynamics;

template <int ndim>
class SphNeighbourSearch;

template <int ndim>
class Nbody;

template <int ndim>
class EOS;

class Paramters;


enum eosenum{noeos, isothermal, locally_isothermal, disc_locally_isothermal, polytropic, barotropic, barotropic2,
             energy_eqn, constant_temp, radws, Nhydroeos};


//=================================================================================================
// struct PrimitiveVariables
/// \brief A good ol' bag o' data wrapper for the primitive hydrodynamic quantities
//=================================================================================================
template<int ndim>
struct PrimitiveVariables {
  PrimitiveVariables()
    : density(0), pressure(0)
  {
    for (int i=0; i < ndim; ++i) velocity[i] = 0;
  }
  FLOAT density;
  FLOAT pressure;
  FLOAT velocity[ndim];
};

//=================================================================================================
// struct ConservedVariables
/// \brief A good ol' bag o' data wrapper for the conserved hydrodynamic quantities
//=================================================================================================
template<int ndim>
struct ConservedVariables {
  ConservedVariables()
    : density(0), energy(0)
  {
    for (int i=0; i < ndim; ++i) momentum[i] = 0;
  }
  FLOAT density;
  FLOAT energy;
  FLOAT momentum[ndim];
};

//=================================================================================================
// struct StateVector
/// \brief  Holds the hydrodynamic state variables
/// \detail The StateVector is designed so to enable Riemann Solvers like the HLLC to compute
//          the fluxes without making any assumptions about the underlying equation of state.
//=================================================================================================
template <int ndim>
struct StateVector {
  StateVector()
    : sound_speed(0)
  { } ;

  PrimitiveVariables<ndim> Wprim;
  ConservedVariables<ndim> Ucons;

  FLOAT sound_speed ;
};

//=================================================================================================
// struct EosParticleProxy
/// \brief  Wraps the particle state data needed by the EOS.
/// \detail This proxy simplifies calculation of the thermal properties between particles, e.g. at
///         the locations of faces in the meshless, without creating of a full particle to pass in.
//=================================================================================================
template<int ndim>
struct EosParticleProxy {
  EosParticleProxy()
  : rho(0), u(0), p(0), ionfrac(0), ionstate(0)
  {
    for (int i=0; i < ndim; ++i) r[i] = 0;
    for (int i=0; i < ndim; ++i) v[i] = 0;
  }
  EosParticleProxy(Particle<ndim>& p)
    : rho(p.rho), u(p.u), p(p.pressure), ionfrac(p.ionfrac), ionstate(p.ionstate)
  {
    for (int i=0; i < ndim; ++i) r[i] = p.r[i];
    for (int i=0; i < ndim; ++i) v[i] = p.v[i];
  }

  FLOAT r[ndim];
  FLOAT v[ndim];
  FLOAT rho;
  FLOAT u;
  FLOAT p;
  FLOAT ionfrac;
  int   ionstate;
};

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
    gammam1(_gamma - (FLOAT) 1.0)
   {};
  virtual ~EOS() {};

  virtual FLOAT Pressure(const EosParticleProxy<ndim>& part) { return gammam1*part.rho*part.u; }
  virtual FLOAT EntropicFunction(const EosParticleProxy<ndim>&) = 0;
  virtual FLOAT SoundSpeed(const EosParticleProxy<ndim>&) = 0;
  virtual FLOAT Temperature(const EosParticleProxy<ndim>&) = 0;
  virtual FLOAT SpecificInternalEnergy(const EosParticleProxy<ndim>&) = 0;

  virtual void set_nbody_data(Nbody<ndim> *) { } ;

  virtual FLOAT InternalEnergyFromPressure(const EosParticleProxy<ndim>&p){
    return p.p / (p.rho * gammam1);
  }
  // Construct State Vectors from primitive / conserved quantities.
  void ConstructStateVector(const EosParticleProxy<ndim>& p, StateVector<ndim>& state){
    state.Wprim.density  = p.rho ;
    state.Ucons.density  = p.rho ;
    state.Wprim.pressure = p.p ;
    state.Ucons.energy   = p.rho*p.u;

    for (int i=0; i < ndim; ++i) {
      state.Wprim.velocity[i] = p.v[i];
      state.Ucons.momentum[i] = p.rho*p.v[i];
      state.Ucons.energy += 0.5*p.rho*p.v[i]*p.v[i];
    }

    state.sound_speed = this->SoundSpeed(p);
  }


  const FLOAT gamma;                             ///< Ratio of specific heats
  const FLOAT gammam1;                           ///< gamma - 1

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
protected:
  using EOS<ndim>::gamma;
  using EOS<ndim>::gammam1;

 public:

  Isothermal(Parameters*, SimUnits *);
  virtual ~Isothermal();

  FLOAT EntropicFunction(const EosParticleProxy<ndim>&);
  FLOAT SoundSpeed(const EosParticleProxy<ndim>&);
  FLOAT Temperature(const EosParticleProxy<ndim>&);
  FLOAT SpecificInternalEnergy(const EosParticleProxy<ndim>&);

  const FLOAT temp0;
  const FLOAT mu_bar;

};



//=================================================================================================
//  Class Polytropic
/// \brief   Polytropic equation of state
/// \details Polytropic equation of state
/// \author  D. A. Hubber, G. Rosotti
/// \date    12/07/2016
//=================================================================================================
template <int ndim>
class Polytropic: public EOS<ndim>
{
  using EOS<ndim>::gamma;
  using EOS<ndim>::gammam1;

 public:

  Polytropic(Parameters*, SimUnits *);
  virtual ~Polytropic();

  FLOAT Pressure(const EosParticleProxy<ndim>&);
  FLOAT EntropicFunction(const EosParticleProxy<ndim>&);
  FLOAT SoundSpeed(const EosParticleProxy<ndim>&);
  FLOAT Temperature(const EosParticleProxy<ndim>&);
  FLOAT SpecificInternalEnergy(const EosParticleProxy<ndim>&);

  const FLOAT Kpoly;
  const FLOAT eta;

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

  Barotropic(Parameters*, SimUnits *);
  virtual ~Barotropic();

  FLOAT EntropicFunction(const EosParticleProxy<ndim>&);
  FLOAT SoundSpeed(const EosParticleProxy<ndim>&);
  FLOAT Temperature(const EosParticleProxy<ndim>&);
  FLOAT SpecificInternalEnergy(const EosParticleProxy<ndim>&);

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

  Barotropic2(Parameters*, SimUnits *);
  virtual ~Barotropic2();

  FLOAT EntropicFunction(const EosParticleProxy<ndim>&);
  FLOAT SoundSpeed(const EosParticleProxy<ndim>&);
  FLOAT Temperature(const EosParticleProxy<ndim>&);
  FLOAT SpecificInternalEnergy(const EosParticleProxy<ndim>&);

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


  Adiabatic(Parameters*, SimUnits *);
  virtual ~Adiabatic();

  FLOAT EntropicFunction(const EosParticleProxy<ndim>&);
  FLOAT SoundSpeed(const EosParticleProxy<ndim>&);
  FLOAT Temperature(const EosParticleProxy<ndim>&);
  FLOAT SpecificInternalEnergy(const EosParticleProxy<ndim>&);

  const FLOAT mu_bar;

};



//=================================================================================================
//  Class Isothermal
/// \brief   Isothermal equation of state
/// \details Isothermal equation of state
/// \author  D. A. Hubber, G. Rosotti
/// \date    03/04/2013
//=================================================================================================
template <int ndim>
class LocallyIsothermal: public Isothermal<ndim>
{
  using Isothermal<ndim>::gamma;
  using Isothermal<ndim>::gammam1;
  using Isothermal<ndim>::temp0;
  using Isothermal<ndim>::mu_bar;

 public:

  LocallyIsothermal(Parameters*, SimUnits *);
  virtual ~LocallyIsothermal();

  FLOAT SpecificInternalEnergy(const EosParticleProxy<ndim>&);
  FLOAT Temperature(const EosParticleProxy<ndim>&);

  virtual void set_nbody_data(Nbody<ndim>* nbody_aux) {
    nbody = nbody_aux;
  } ;

private:
  FLOAT templaw;
  FLOAT tempmin;

  Nbody<ndim>* nbody;
};


//=================================================================================================
//  Class DiscLocallyIsothermal
/// \brief   Locally isothermal equation of state for discs
/// \details Locally isothermal equation of state for discs
/// \author  G. Rosotti
/// \date    03/04/2013
//=================================================================================================
template <int ndim>
class DiscLocallyIsothermal: public Isothermal<ndim>
{
  using Isothermal<ndim>::gamma;
  using Isothermal<ndim>::gammam1;
  using Isothermal<ndim>::temp0;
  using Isothermal<ndim>::mu_bar;

 public:

  DiscLocallyIsothermal(Parameters*, SimUnits *);
  virtual ~DiscLocallyIsothermal();

  FLOAT SoundSpeed(const EosParticleProxy<ndim>&);
  FLOAT SpecificInternalEnergy(const EosParticleProxy<ndim>&);
  FLOAT Temperature(const EosParticleProxy<ndim>&);

  virtual void set_nbody_data(Nbody<ndim>* nbody_aux) {
    nbody = nbody_aux;
  } ;

private:
  FLOAT slope;
  FLOAT norm;
  FLOAT rin;

  Nbody<ndim>* nbody;
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

  IonisingRadiation(Parameters*, SimUnits *);
  virtual ~IonisingRadiation();

  FLOAT EntropicFunction(const EosParticleProxy<ndim>&);
  FLOAT SoundSpeed(const EosParticleProxy<ndim>&);
  FLOAT Temperature(const EosParticleProxy<ndim>&);
  FLOAT SpecificInternalEnergy(const EosParticleProxy<ndim>&);

  void set_nbody_data(Nbody<ndim>* nbody) {
    eos->set_nbody_data(nbody);
  }

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

 public:

  MCRadiationEOS(Parameters*, SimUnits *);
  virtual ~MCRadiationEOS();

  FLOAT EntropicFunction(const EosParticleProxy<ndim>&);
  FLOAT SoundSpeed(const EosParticleProxy<ndim>&);
  FLOAT Temperature(const EosParticleProxy<ndim>&);
  FLOAT SpecificInternalEnergy(const EosParticleProxy<ndim>&);

  void set_nbody_data(Nbody<ndim>* nbody) {
    eos->set_nbody_data(nbody);
  }

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

 public:

  Radws(Parameters*, SimUnits *);
  virtual ~Radws();

  FLOAT Pressure(const EosParticleProxy<ndim>&);
  FLOAT EntropicFunction(const EosParticleProxy<ndim>&);
  FLOAT SoundSpeed(const EosParticleProxy<ndim>&);
  FLOAT Temperature(const EosParticleProxy<ndim>&);
  FLOAT SpecificInternalEnergy(const EosParticleProxy<ndim>&);
  virtual FLOAT InternalEnergyFromPressure(const EosParticleProxy<ndim>&p);

  OpacityTable<ndim> *opacity_table;
};
#endif
