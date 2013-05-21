//=============================================================================
//  EOS.h
//  Contains (virtual) definitions for equation of state class.  
//  Also includes inherited class defintions for various equation of  
//  state options.
//=============================================================================


#ifndef _EOS_H_
#define _EOS_H_


#include "Precision.h"
#include "Constants.h"
#include "Parameters.h"
#include "SphParticle.h"



//=============================================================================
//  Class EOS
/// \brief   Main equation of state
/// \details Main equation of state class.  Only contains virtual function 
///          defintions. All functions must be defined by the inherited 
///          EOS classes.
/// \author  D. A. Hubber, G. Rosotti
/// \date    03/04/2013
//=============================================================================
template <int ndim>
class EOS
{
 public:

  EOS(FLOAT _gamma):
    gamma(_gamma),
    gammam1(gamma - (FLOAT) 1.0){};

  virtual FLOAT Pressure(SphParticle<ndim> &) = 0;
  virtual FLOAT EntropicFunction(SphParticle<ndim> &) = 0;
  virtual FLOAT SoundSpeed(SphParticle<ndim> &) = 0;
  virtual FLOAT Temperature(SphParticle<ndim> &) = 0;
  virtual FLOAT SpecificInternalEnergy(SphParticle<ndim> &) = 0;

  FLOAT gamma;
  FLOAT gammam1;

};



//=============================================================================
//  Class Isothermal
/// \brief   Isothermal equation of state
/// \details Isothermal equation of state
/// \author  D. A. Hubber, G. Rosotti
/// \date    03/04/2013
//=============================================================================
template <int ndim>
class Isothermal: public EOS<ndim>
{
  using EOS<ndim>::gamma;
  using EOS<ndim>::gammam1;

 public:

  Isothermal(FLOAT, FLOAT, FLOAT);
  ~Isothermal();

  FLOAT Pressure(SphParticle<ndim> &);
  FLOAT EntropicFunction(SphParticle<ndim> &);
  FLOAT SoundSpeed(SphParticle<ndim> &);
  FLOAT Temperature(SphParticle<ndim> &);
  FLOAT SpecificInternalEnergy(SphParticle<ndim> &);

  FLOAT temp0;
  FLOAT mu_bar;

};



//=============================================================================
//  Class Barotropic
/// \brief   Barotropic equation of state
/// \details Barotropic equation of state
/// \author  D. A. Hubber, G. Rosotti
/// \date    03/04/2013
//=============================================================================
template <int ndim>
class Barotropic: public EOS<ndim>
{
  using EOS<ndim>::gamma;
  using EOS<ndim>::gammam1;

 public:

  Barotropic(FLOAT, FLOAT, FLOAT, FLOAT);
  ~Barotropic();

  FLOAT Pressure(SphParticle<ndim> &);
  FLOAT EntropicFunction(SphParticle<ndim> &);
  FLOAT SoundSpeed(SphParticle<ndim> &);
  FLOAT Temperature(SphParticle<ndim> &);
  FLOAT SpecificInternalEnergy(SphParticle<ndim> &);

  FLOAT temp0;
  FLOAT mu_bar;
  FLOAT rho_bary;
  FLOAT invrho_bary;

};



//=============================================================================
//  Class Adiabatic
/// \brief   Adiabatic equation of state
/// \details Adiabatic equation of state.  Requires integrating the energy 
///          equation parallel to the main dynamical quantities.
/// \author  D. A. Hubber, G. Rosotti
/// \date    03/04/2013
//=============================================================================
template <int ndim>
class Adiabatic: public EOS<ndim>
{
  using EOS<ndim>::gamma;
  using EOS<ndim>::gammam1;

 public:

  Adiabatic(FLOAT, FLOAT, FLOAT);
  ~Adiabatic();

  FLOAT Pressure(SphParticle<ndim> &);
  FLOAT EntropicFunction(SphParticle<ndim> &);
  FLOAT SoundSpeed(SphParticle<ndim> &);
  FLOAT Temperature(SphParticle<ndim> &);
  FLOAT SpecificInternalEnergy(SphParticle<ndim> &);

  FLOAT temp0;
  FLOAT mu_bar;

};
#endif
