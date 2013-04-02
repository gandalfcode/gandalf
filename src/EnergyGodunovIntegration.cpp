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
#include "Exception.h"
#include "Debug.h"
using namespace std;




// ============================================================================
// EnergyGodunovIntegration::EnergyGodunovIntegration()
// ============================================================================
template <int ndim>
EnergyGodunovIntegration<ndim>::EnergyGodunovIntegration(DOUBLE energy_mult_aux) :
  EnergyEquation<ndim>(energy_mult_aux)
{
}



// ============================================================================
// EnergyGodunovIntegration::~EnergyGodunovIntegration()
// ============================================================================
template <int ndim>
EnergyGodunovIntegration<ndim>::~EnergyGodunovIntegration()
{
}



// ============================================================================
// EnergyGodunovIntegration::EnergyIntegration
// Integrate internal energy to first order from the beginning of the step to 
// the current simulation time, i.e. u(t+dt) = u(t) + dudt(t)*dt
// ============================================================================
template <int ndim>
void EnergyGodunovIntegration<ndim>::EnergyIntegration(int n, int level_step, int Nsph,
				  SphParticle<ndim> *sph, FLOAT timestep)
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
    //if (sph[i].u < small_number && sph[i].dudt < 0.0)
    //  sph[i].u = sph[i].u0*exp(dt*sph[i].dudt/sph[i].u0);
    if (sph[i].u != sph[i].u) {
      cout << "Something wrong with energy integration (NaN) : " << endl;
      cout << sph[i].u << "   " << sph[i].u0 << "   " << sph[i].dudt 
	   << "   " << dt << "   " << nstep << "    " << timestep << endl;
      exit(0);
    }
    if (sph[i].u < small_number) {
      cout << "Something wrong with energy integration (0) : " << endl;
      cout << sph[i].u << "   " << sph[i].u0 << "   " << sph[i].dudt 
	   << "   " << dt << "   " << nstep << "    " 
	   << sph[i].u0/sph[i].dudt << endl;
      cout << " dt_courant : " << sph[i].h/sph[i].sound << "   " 
	   << sph[i].u0/(sph[i].dudt + small_number) << endl;
      string message = "Problem with energy integration (0)";
      ExceptionHandler::getIstance().raise(message);
    }
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
template <int ndim>
void EnergyGodunovIntegration<ndim>::EnergyCorrectionTerms(int n, int level_step,
						     int Nsph, 
						     SphParticle<ndim> *sph,
						     FLOAT timestep)
{
  return;
}



// ============================================================================
// EnergyGodunovIntegration::EndTimestep
// Record all important thermal quantities at the end of the step for the 
// start of the new timestep.
// ============================================================================
template <int ndim>
void EnergyGodunovIntegration<ndim>::EndTimestep(int n, int level_step, int Nsph, SphParticle<ndim> *sph)
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
template <int ndim>
DOUBLE EnergyGodunovIntegration<ndim>::Timestep(SphParticle<ndim> &part)
{
  return this->energy_mult*(DOUBLE) (part.u/(fabs(part.dudt) + small_number));
}

template class EnergyGodunovIntegration<1>;
template class EnergyGodunovIntegration<2>;
template class EnergyGodunovIntegration<3>;
