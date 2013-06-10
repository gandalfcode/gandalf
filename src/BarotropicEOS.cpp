//=============================================================================
//  BarotropicEOS.cpp
//  Contains all function definitions for a barotropic Equation of state of 
//  the form T = temp0*(1 + (rho/rho_bary)^{gamma - 1}).
//  Used for star formation simulations to approximate the combined isothermal 
//  and optically-thich adiabatic regimes of the gas collapse phase.
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
Barotropic<ndim>::Barotropic(FLOAT temp0aux, FLOAT mu_bar_aux, FLOAT gamma_aux,
			     FLOAT rho_bary_aux, SimUnits *units):
  EOS<ndim>(gamma_aux)
{
  temp0 = temp0aux/units->temp.outscale;
  mu_bar = mu_bar_aux;
  rho_bary = rho_bary_aux/units->rho.outscale/units->rho.outcgs;
  invrho_bary = 1.0/rho_bary;
}



//=============================================================================
//  Barotropic::~Barotropic()
/// Barotropic EOS destructor
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
/// Returns specific internal energy
//=============================================================================
template <int ndim>
FLOAT Barotropic<ndim>::SpecificInternalEnergy(SphParticle<ndim> &part)
{
  return temp0*(1.0 + pow(part.rho*invrho_bary,gammam1))/gammam1/mu_bar;
}



//=============================================================================
//  Barotropic::Temperature
/// Returns temperature of particle.  Approximates gas in the isothermal 
/// regime (T = temp0 for rho << rho_bary) and in the optically thick 
/// adiabatic phase (T = const*rho^{gamma - 1} for rho >> rho_bary).
//=============================================================================
template <int ndim>
FLOAT Barotropic<ndim>::Temperature(SphParticle<ndim> &part)
{
  return temp0*(1.0 + pow(part.rho*invrho_bary,gammam1));
}



template class Barotropic<1>;
template class Barotropic<2>;
template class Barotropic<3>;

