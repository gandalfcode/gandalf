// ============================================================================
// EnergyPEC.cpp
// Contains functions for energy equation integration using a 
// Predict-Evalulate-Correct (PEC) scheme.
// N.B. this PEC scheme is the same as integrating the particle velocities 
// in the Leapfrog KDK scheme.
// ============================================================================


#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <math.h>
#include "Sph.h"
#include "SphKernel.h"
#include "SphIntegration.h"
#include "SphParticle.h"
#include "EOS.h"
#include "EnergyEquation.h"
#include "Debug.h"
using namespace std;



// ============================================================================
// EnergyEquation::EnergyEquation()
// ============================================================================
EnergyEquation::EnergyEquation(DOUBLE energy_mult_aux) :
  energy_mult(energy_mult_aux)
{
}

EnergyEquation::~EnergyEquation()
{
}



// ============================================================================
// EnergyPEC::EnergyPEC()
// ============================================================================
EnergyPEC::EnergyPEC(DOUBLE energy_mult_aux) :
  EnergyEquation(energy_mult_aux)
{
}



// ============================================================================
// EnergyPEC::~EnergyPEC()
// ============================================================================
EnergyPEC::~EnergyPEC()
{
}



// ============================================================================
// EnergyPEC::EnergyIntegration
// Integrate internal energy to first order from the beginning of the step to 
// the current simulation time, i.e. u(t+dt) = u(t) + dudt(t)*dt
// ============================================================================
void EnergyPEC::EnergyIntegration(int n, int level_step, int Nsph,
				  SphParticle *sph, FLOAT timestep)
{
  int i;                                // Particle counter
  int nstep;                            // Particle (integer) step size
  FLOAT dt;                             // Timestep since start of step

  debug2("[EnergyPEC::EnergyIntegration]");

  for (i=0; i<Nsph; i++) {
    nstep = pow(2,level_step - sph[i].level);
    if (n%nstep == 0) dt = timestep*(FLOAT) nstep;
    else dt = timestep*(FLOAT) (n%nstep);
    sph[i].u = sph[i].u0 + sph[i].dudt*dt;
  }

  return;
}
 


// ============================================================================
// EnergyPEC::CorrectionTerms
// Compute energy integration to second order at the end of the step by 
// adding a second order correction term.  The full integration becomes
// u(t+dt) = u(t) + 0.5*(dudt(t) + dudt(t+dt))*dt 
// ============================================================================
void EnergyPEC::EnergyCorrectionTerms(int n, int level_step, int Nsph, 
				      SphParticle *sph, FLOAT timestep)
{
  int i;                                // Particle counter
  int nstep;                            // Particle (integer) step size

  debug2("[EnergyPEC::EnergyCorrectionTerms]");

  for (i=0; i<Nsph; i++) {
    nstep = pow(2,level_step - sph[i].level);
    if (n%nstep == 0) 
      sph[i].u += 0.5*(sph[i].dudt - sph[i].dudt0)*timestep*(FLOAT) nstep;
  }

  return;
}



// ============================================================================
// EnergyPEC::EndTimestep
// Record all important thermal quantities at the end of the step for the 
// start of the new timestep.
// ============================================================================
void EnergyPEC::EndTimestep(int n, int level_step, int Nsph, SphParticle *sph)
{
  int i;                                // Particle counter
  int nstep;                            // Particle (integer) step size

  debug2("[EnergyPEC::EndTimestep]");

  for (i=0; i<Nsph; i++) {
    nstep = pow(2,level_step - sph[i].level);
    if (n%nstep == 0) {
      sph[i].u0 = sph[i].u;
      sph[i].dudt0 = sph[i].dudt;
    }
  }

  return;
}



// ============================================================================
// EnergyPEC::Timestep
// Compute explicit timestep such that u cannot change by a large fraction 
// in one step, i.e. dt = const*u/|dudt + epsilon| 
// where epsilon is to prevent the denominator becoming zero.
// ============================================================================
DOUBLE EnergyPEC::Timestep(SphParticle &part)
{
  return energy_mult*(DOUBLE) (part.u/(fabs(part.dudt) + small_number));
}
