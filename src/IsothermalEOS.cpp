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
Isothermal::Isothermal(FLOAT temp0aux, FLOAT mu_bar_aux, FLOAT gamma_aux)
{
  temp0 = temp0aux;
  mu_bar = mu_bar_aux;
  gamma = gamma_aux;
  gammam1 = gamma - (FLOAT) 1.0;
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
FLOAT Isothermal::Pressure(SphParticle &part)
{
  return gammam1*part.rho*part.u;
}



// ============================================================================
// Isothermal::EntropicFunction
// Calculates and returns value of Entropic function (= P/rho^gamma) for 
// referenced particle
// ============================================================================
FLOAT Isothermal::EntropicFunction(SphParticle &part)
{
  return gammam1*part.u*pow(part.rho,(FLOAT) 1.0 - gamma);
}



// ============================================================================
// Isothermal::SoundSpeed
// Returns isothermal sound speed of SPH particle
// ============================================================================
FLOAT Isothermal::SoundSpeed(SphParticle &part)
{
  return sqrt(gammam1*part.u);
}



// ============================================================================
// Isothermal::SpecificInternalEnergy
// ============================================================================
FLOAT Isothermal::SpecificInternalEnergy(SphParticle &part)
{
  return temp0/gammam1/mu_bar;
}



// ============================================================================
// Isothermal::Temperature
// Return isothermal temperature of particle
// ============================================================================
FLOAT Isothermal::Temperature(SphParticle &part)
{
  return temp0;
}
