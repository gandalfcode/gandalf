// ============================================================================
// EnergyPEC.cpp
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
// EnergyPEC::EnergyPEC()
// ============================================================================
EnergyEquation::EnergyEquation(double energy_mult_aux) :
  energy_mult(energy_mult_aux)
{
}

EnergyEquation::~EnergyEquation()
{
}



// ============================================================================
// EnergyPEC::EnergyPEC()
// ============================================================================
EnergyPEC::EnergyPEC(double energy_mult_aux) :
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
// ============================================================================
void EnergyPEC::EnergyIntegration(int n, int level_step, 
				  int Nsph, SphParticle *sph, double dt)
{
  for (int i=0; i<Nsph; i++)
    sph[i].u = sph[i].u0 + sph[i].dudt*dt;

  return;
}
 


// ============================================================================
// EnergyPEC::CorrectionTerms
// ============================================================================
void EnergyPEC::EnergyCorrectionTerms(int n, int level_step, 
				      int Nsph, SphParticle *sph, double dt)
{
  int i,nstep;

  for (i=0; i<Nsph; i++) {
    nstep = pow(2,level_step - sph[i].level);
    if (n%nstep == 0) 
      sph[i].u = sph[i].u0 + 0.5*(sph[i].dudt + sph[i].dudt0)*dt;
  }

  return;
}



// ============================================================================
// EnergyPEC::EndTimestep
// ============================================================================
void EnergyPEC::EndTimestep(int n, int level_step, int Nsph, SphParticle *sph)
{
  int i,nstep;

  debug2("[EnergyPEC::EndTimestep]\n");

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
// ..
// ============================================================================
double EnergyPEC::Timestep(SphParticle &part)
{
  return energy_mult*part.u/(fabs(part.dudt) + small_number_dp);
}
