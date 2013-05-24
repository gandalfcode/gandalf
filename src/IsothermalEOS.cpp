//=============================================================================
//  IsothermalEOS.cpp
//  Contains all function definitions for the Isothermal Equation of state.
//=============================================================================


#include <math.h>
#include "EOS.h"
#include "Sph.h"


//=============================================================================
//  Isothermal::Isothermal()
/// Default constructor for isothermal EOS.  Passes and sets important 
/// thermal physics variables.
//=============================================================================
template <int ndim>
Isothermal<ndim>::Isothermal(FLOAT temp0aux, FLOAT mu_bar_aux, 
                             FLOAT gamma_aux, SimUnits *units):
  EOS<ndim>(gamma_aux)
{
  temp0 = temp0aux/units->temp.outscale;
  mu_bar = mu_bar_aux;
}



//=============================================================================
//  Isothermal::Isothermal()
//=============================================================================
template <int ndim>
Isothermal<ndim>::~Isothermal()
{
}



//=============================================================================
//  Isothermal::Pressure
/// Calculates and returns thermal pressure of referenced particle
//=============================================================================
template <int ndim>
FLOAT Isothermal<ndim>::Pressure(SphParticle<ndim> &part)
{
  return gammam1*part.rho*part.u;
}



//=============================================================================
//  Isothermal::EntropicFunction
/// Calculates and returns value of Entropic function (= P/rho^gamma) for 
/// referenced particle
//=============================================================================
template <int ndim>
FLOAT Isothermal<ndim>::EntropicFunction(SphParticle<ndim> &part)
{
  return gammam1*part.u*pow(part.rho,(FLOAT) 1.0 - gamma);
}



//=============================================================================
//  Isothermal::SoundSpeed
/// Returns isothermal sound speed of SPH particle
//=============================================================================
template <int ndim>
FLOAT Isothermal<ndim>::SoundSpeed(SphParticle<ndim> &part)
{
  return sqrt(gammam1*part.u);
}



//=============================================================================
//  Isothermal::SpecificInternalEnergy
//=============================================================================
template <int ndim>
FLOAT Isothermal<ndim>::SpecificInternalEnergy(SphParticle<ndim> &part)
{
  return temp0/gammam1/mu_bar;
}



//=============================================================================
//  Isothermal::Temperature
/// Return isothermal temperature of particle
//=============================================================================
template <int ndim>
FLOAT Isothermal<ndim>::Temperature(SphParticle<ndim> &part)
{
  return temp0;
}



template class Isothermal<1>;
template class Isothermal<2>;
template class Isothermal<3>;
