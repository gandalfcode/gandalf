// ============================================================================
// AdiabaticEOS.cpp
// Contains all function definitions for the Adiabatic Equation of state.
// ============================================================================


#include <math.h>
#include "EOS.h"
#include "Sph.h"



// ============================================================================
// Adiabatic::Adiabatic()
// Default constructor for perfect gas EOS.  Passes and sets important 
// thermal physics variables.
// ============================================================================
template <int ndim>
Adiabatic<ndim>::Adiabatic(FLOAT temp0aux, FLOAT mu_bar_aux, FLOAT gamma_aux):
  EOS<ndim> (gamma_aux)
{
  temp0 = temp0aux;
  mu_bar = mu_bar_aux;
}



// ============================================================================
// Adiabatic::Adiabatic()
// ============================================================================
template <int ndim>
Adiabatic<ndim>::~Adiabatic()
{
}



// ============================================================================
// Adiabatic::Pressure
// Calculates and returns thermal pressure of referenced particle
// ============================================================================
template <int ndim>
FLOAT Adiabatic<ndim>::Pressure(SphParticle<ndim> &part)
{
  return gammam1*part.rho*part.u;
}



// ============================================================================
// Adiabatic::EntropicFunction
// Calculates and returns value of Entropic function (= P/rho^gamma) for 
// referenced particle
// ============================================================================
template <int ndim>
FLOAT Adiabatic<ndim>::EntropicFunction(SphParticle<ndim> &part)
{
  return gammam1*part.u*pow(part.rho,(FLOAT) 1.0 - gamma);
}



// ============================================================================
// Adiabatic::SoundSpeed
// Returns adiabatic sound speed of particle
// ============================================================================
template <int ndim>
FLOAT Adiabatic<ndim>::SoundSpeed(SphParticle<ndim> &part)
{
  return sqrt(gamma*gammam1*part.u);
}



// ============================================================================
// Adiabatic::SpecificInternalEnergy
// Returns specific internal energy of particle
// ============================================================================
template <int ndim>
FLOAT Adiabatic<ndim>::SpecificInternalEnergy(SphParticle<ndim> &part)
{
  return part.u;
}



// ============================================================================
// Adiabatic::Temperature
// Returns temperature of particle
// ============================================================================
template <int ndim>
FLOAT Adiabatic<ndim>::Temperature(SphParticle<ndim> &part)
{
  return gammam1*part.u;
}

template class Adiabatic<1>;
template class Adiabatic<2>;
template class Adiabatic<3>;
