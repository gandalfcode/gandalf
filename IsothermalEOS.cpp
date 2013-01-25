// ============================================================================
// IsothermalEOS.cpp
// ============================================================================


#include <math.h>
#include "EOS.h"
#include "Sph.h"
#include "Parameters.h"



// ============================================================================
// Isothermal::Isothermal()
// ============================================================================
Isothermal::Isothermal(float temp0aux, float mu_bar_aux, float gamma_aux)
{
  temp0 = temp0aux;
  mu_bar = mu_bar_aux;
  gamma = gamma_aux;
}



// ============================================================================
// Isothermal::Isothermal()
// ============================================================================
Isothermal::~Isothermal()
{
}


// ============================================================================
// Isothermal::Pressure
// ============================================================================
float Isothermal::Pressure(SphParticle &part)
{
  return (gamma - 1.0)*part.rho*part.u;
}



// ============================================================================
// Isothermal::EntropicFunction
// ============================================================================
float Isothermal::EntropicFunction(SphParticle &part)
{
  return (gamma - 1.0)*part.u*pow(part.rho,1.0 - gamma);
}



// ============================================================================
// Isothermal::SoundSpeed
// ============================================================================
float Isothermal::SoundSpeed(SphParticle &part)
{
  return sqrt((gamma - 1.0)*part.u);
}



// ============================================================================
// Isothermal::SoundSpeed
// ============================================================================
float Isothermal::SpecificInternalEnergy(SphParticle &part)
{
  return temp0/(gamma - 1.0)/mu_bar;
}



// ============================================================================
// Isothermal::Temperature
// ============================================================================
float Isothermal::Temperature(SphParticle &part)
{
  return temp0;
}
