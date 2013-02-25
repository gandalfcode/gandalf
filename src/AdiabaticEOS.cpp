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
Adiabatic::Adiabatic(FLOAT temp0aux, FLOAT mu_bar_aux, FLOAT gamma_aux)
{
  temp0 = temp0aux;
  mu_bar = mu_bar_aux;
  gamma = gamma_aux;
  gammam1 = gamma - 1.0;
}



// ============================================================================
// Adiabatic::Adiabatic()
// ============================================================================
Adiabatic::~Adiabatic()
{
}



// ============================================================================
// Adiabatic::Pressure
// Calculates and returns thermal pressure of referenced particle
// ============================================================================
FLOAT Adiabatic::Pressure(SphParticle &part)
{
  return (gamma - 1.0)*part.rho*part.u;
}



// ============================================================================
// Adiabatic::EntropicFunction
// Calculates and returns value of Entropic function (= P/rho^gamma) for 
// referenced particle
// ============================================================================
FLOAT Adiabatic::EntropicFunction(SphParticle &part)
{
  return (gamma - 1.0)*part.u*pow(part.rho,1.0 - gamma);
}



// ============================================================================
// Adiabatic::SoundSpeed
// Returns adiabatic sound speed of particle
// ============================================================================
FLOAT Adiabatic::SoundSpeed(SphParticle &part)
{
  return sqrt(gamma*(gamma - 1.0)*part.u);
}



// ============================================================================
// Adiabatic::SpecificInternalEnergy
// Returns specific internal energy of particle
// ============================================================================
FLOAT Adiabatic::SpecificInternalEnergy(SphParticle &part)
{
  return part.u;
}



// ============================================================================
// Adiabatic::Temperature
// Returns temperature of particle
// ============================================================================
FLOAT Adiabatic::Temperature(SphParticle &part)
{
  return gammam1*part.u;
}
