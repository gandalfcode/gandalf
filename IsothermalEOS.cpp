// ============================================================================
// IsothermalEOS.cpp
// Contains all function definitions for the Isothermal Equation of state.
// ============================================================================


#include <math.h>
#include "EOS.h"
#include "Sph.h"


// ============================================================================
// Isothermal::Isothermal()
// Default constructor for isothermal EOS.  Passes and sets important 
// thermal physics variables.
// ============================================================================
Isothermal::Isothermal(float temp0aux, float mu_bar_aux, float gamma_aux)
{
  temp0 = temp0aux;
  mu_bar = mu_bar_aux;
  gamma = gamma_aux;
  gammam1 = gamma - 1.0;
}



// ============================================================================
// Isothermal::Isothermal()
// ============================================================================
Isothermal::~Isothermal()
{
}



// ============================================================================
// Isothermal::Pressure
// Calculates and returns thermal pressure of referenced particle
// ============================================================================
float Isothermal::Pressure(SphParticle &part)
{
  return (gamma - 1.0)*part.rho*part.u;
}



// ============================================================================
// Isothermal::EntropicFunction
// Calculates and returns value of Entropic function (= P/rho^gamma) for 
// referenced particle
// ============================================================================
float Isothermal::EntropicFunction(SphParticle &part)
{
  return (gamma - 1.0)*part.u*pow(part.rho,1.0 - gamma);
}



// ============================================================================
// Isothermal::SoundSpeed
// Returns isothermal sound speed of SPH particle
// ============================================================================
float Isothermal::SoundSpeed(SphParticle &part)
{
  return sqrt((gamma - 1.0)*part.u);
}



// ============================================================================
// Isothermal::SpecificInternalEnergy
// ============================================================================
float Isothermal::SpecificInternalEnergy(SphParticle &part)
{
  return temp0/(gamma - 1.0)/mu_bar;
}



// ============================================================================
// Isothermal::Temperature
// Return isothermal temperature of particle
// ============================================================================
float Isothermal::Temperature(SphParticle &part)
{
  return temp0;
}
