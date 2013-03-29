// ============================================================================
// EnergyGodunovIntegration.cpp
// ..
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
// EnergyGodunovIntegration::EnergyGodunovIntegration()
// ============================================================================
EnergyGodunovIntegration::EnergyGodunovIntegration(DOUBLE energy_mult_aux) :
  EnergyEquation(energy_mult_aux)
{
}



// ============================================================================
// EnergyGodunovIntegration::~EnergyGodunovIntegration()
// ============================================================================
EnergyGodunovIntegration::~EnergyGodunovIntegration()
{
}



// ============================================================================
// EnergyGodunovIntegration::EnergyIntegration
// Integrate internal energy to first order from the beginning of the step to 
// the current simulation time, i.e. u(t+dt) = u(t) + dudt(t)*dt
// ============================================================================
void EnergyGodunovIntegration::EnergyIntegration(int n, int level_step, int Nsph,
				  SphParticle *sph, FLOAT timestep)
{
  int i;                                // Particle counter
  int nstep;                            // Particle (integer) step size
  FLOAT dt;                             // Timestep since start of step

  debug2("[EnergyGodunovIntegration::EnergyIntegration]");

  // --------------------------------------------------------------------------
#pragma omp parallel for default(shared) private(dt,nstep)
  for (i=0; i<Nsph; i++) {
    nstep = pow(2,level_step - sph[i].level);
    if (n%nstep == 0) dt = timestep*(FLOAT) nstep;
    else dt = timestep*(FLOAT) (n%nstep);
    sph[i].u = sph[i].u0 + sph[i].dudt*dt;
  }
  // --------------------------------------------------------------------------

  return;
}
 


// ============================================================================
// EnergyGodunovIntegration::CorrectionTerms
// Compute energy integration to second order at the end of the step by 
// adding a second order correction term.  The full integration becomes
// u(t+dt) = u(t) + 0.5*(dudt(t) + dudt(t+dt))*dt 
// ============================================================================
void EnergyGodunovIntegration::EnergyCorrectionTerms(int n, int level_step, 
						     int Nsph, 
						     SphParticle *sph, 
						     FLOAT timestep)
{
  return;
}



// ============================================================================
// EnergyGodunovIntegration::EndTimestep
// Record all important thermal quantities at the end of the step for the 
// start of the new timestep.
// ============================================================================
void EnergyGodunovIntegration::EndTimestep(int n, int level_step, int Nsph, SphParticle *sph)
{
  int i;                                // Particle counter
  int nstep;                            // Particle (integer) step size

  debug2("[EnergyGodunovIntegration::EndTimestep]");

  // --------------------------------------------------------------------------
#pragma omp parallel for default(shared) private(nstep)
  for (i=0; i<Nsph; i++) {
    nstep = pow(2,level_step - sph[i].level);
    if (n%nstep == 0) {
      sph[i].u0 = sph[i].u;
      sph[i].dudt0 = sph[i].dudt;
    }
  }
  // --------------------------------------------------------------------------

  return;
}



// ============================================================================
// EnergyGodunovIntegration::Timestep
// Compute explicit timestep such that u cannot change by a large fraction 
// in one step, i.e. dt = const*u/|dudt + epsilon| 
// where epsilon is to prevent the denominator becoming zero.
// ============================================================================
DOUBLE EnergyGodunovIntegration::Timestep(SphParticle &part)
{
  return energy_mult*(DOUBLE) (part.u/(fabs(part.dudt) + small_number));
}
