//=============================================================================
//  BarotropicEOS.cpp
//  Contains all function definitions for the Barotropic Equation of state.
//=============================================================================


#include <math.h>
#include "EOS.h"
#include "Sph.h"


//=============================================================================
//  Barotropic::Barotropic()
/// Default constructor for barotropic EOS.  Passes and sets important 
/// thermal physics variables.
//=============================================================================
template <int ndim>
Barotropic<ndim>::Barotropic(FLOAT temp0aux, FLOAT mu_bar_aux, 
			     FLOAT gamma_aux, FLOAT rho_bary_aux):
  EOS<ndim>(gamma_aux)
{
  temp0 = temp0aux;
  mu_bar = mu_bar_aux;
  rho_bary = rho_bary_aux;
  invrho_bary = 1.0/rho_bary_aux;
}



//=============================================================================
//  Barotropic::Barotropic()
//=============================================================================
template <int ndim>
Barotropic<ndim>::~Barotropic()
{
}



//=============================================================================
//  Barotropic::Pressure
/// Calculates and returns thermal pressure of referenced particle
//=============================================================================
template <int ndim>
FLOAT Barotropic<ndim>::Pressure(SphParticle<ndim> &part)
{
  return gammam1*part.rho*part.u;
}



//=============================================================================
//  Barotropic::EntropicFunction
/// Calculates and returns value of Entropic function (= P/rho^gamma) for 
/// referenced particle
//=============================================================================
template <int ndim>
FLOAT Barotropic<ndim>::EntropicFunction(SphParticle<ndim> &part)
{
  return gammam1*part.u*pow(part.rho,(FLOAT) 1.0 - gamma);
}



//=============================================================================
//  Barotropic::SoundSpeed
/// Returns isothermal sound speed of SPH particle
//=============================================================================
template <int ndim>
FLOAT Barotropic<ndim>::SoundSpeed(SphParticle<ndim> &part)
{
  return sqrt(gammam1*part.u);
}



//=============================================================================
//  Barotropic::SpecificInternalEnergy
//=============================================================================
template <int ndim>
FLOAT Barotropic<ndim>::SpecificInternalEnergy(SphParticle<ndim> &part)
{
  return temp0*(1.0 + pow(part.rho*invrho_bary,gammam1))/gammam1/mu_bar;
}



//=============================================================================
//  Barotropic::Temperature
/// Return isothermal temperature of particle
//=============================================================================
template <int ndim>
FLOAT Barotropic<ndim>::Temperature(SphParticle<ndim> &part)
{
  return temp0*(1.0 + pow(part.rho*invrho_bary,gammam1));
}



template class Barotropic<1>;
template class Barotropic<2>;
template class Barotropic<3>;

