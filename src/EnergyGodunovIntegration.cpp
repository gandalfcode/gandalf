//=============================================================================
//  EnergyGodunovIntegration.cpp
//  Energy integration scheme for Inutsuka (2002) Godunov SPH algorithm.
//  Conserves total energy (thermal + kinetic) to machine precision for 
//  direct sum forces and global timesteps.
//=============================================================================


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




//=============================================================================
//  EnergyGodunovIntegration::EnergyGodunovIntegration()
/// EnergyGodunovIntegration class constructor
//=============================================================================
template <int ndim>
EnergyGodunovIntegration<ndim>::EnergyGodunovIntegration(DOUBLE energy_mult_aux) :
  EnergyEquation<ndim>(energy_mult_aux)
{
}



//=============================================================================
//  EnergyGodunovIntegration::~EnergyGodunovIntegration()
/// EnergyGodunovIntegration class destructor
//=============================================================================
template <int ndim>
EnergyGodunovIntegration<ndim>::~EnergyGodunovIntegration()
{
}



//=============================================================================
//  EnergyGodunovIntegration::EnergyIntegration
/// Integrate internal energy to first order from the beginning of the step to 
/// the current simulation time, i.e. $u(t+dt) = u(t) + dudt(t)*dt$
//=============================================================================
template <int ndim>
void EnergyGodunovIntegration<ndim>::EnergyIntegration
(int n,                             ///< [in] Integer time in block time struct
 int Nsph,                          ///< [in] No. of SPH particles
 SphParticle<ndim> *sphdata,        ///< [inout] SPH particle data array
 FLOAT timestep)                    ///< [in] Base timestep value
{
  int dn;                           // Integer time since beginning of step
  int i;                            // Particle counter
  int nstep;                        // Particle (integer) step size
  FLOAT dt;                         // Timestep since start of step

  debug2("[EnergyGodunovIntegration::EnergyIntegration]");

  // --------------------------------------------------------------------------
#pragma omp parallel for default(shared) private(dt,nstep)
  for (i=0; i<Nsph; i++) {
    nstep = sphdata[i].nstep;
    dn = n - sphdata[i].nlast;
    dt = timestep*(FLOAT) dn;
    sphdata[i].u = sphdata[i].u0 + sphdata[i].dudt0*dt;

    //cout << "CHECKING UINT : " << i << "   " << nstep << "    " << dn << "    " << dt << "    " << sphdata[i].u << "    " << sphdata[i].u0 << "    " << sphdata[i].dudt0 << endl;

    //if (sphdata[i].u < small_number && sphdata[i].dudt < 0.0)
    //  sphdata[i].u = sphdata[i].u0*exp(dt*sphdata[i].dudt/sphdata[i].u0);
    if (sphdata[i].u != sphdata[i].u) {
      cout << "Something wrong with energy integration (NaN) : " << endl;
      cout << sphdata[i].u << "   " << sphdata[i].u0 << "  " << sphdata[i].dudt
	   << "   " << dt << "   " << nstep << "    " << timestep << endl;
      exit(0);
    }
    if (sphdata[i].u < small_number) {
      cout << "Something wrong with energy integration (0) : " << endl;
      cout << sphdata[i].u << "   " << sphdata[i].u0 << "  " << sphdata[i].dudt
	   << "   " << dt << "   " << nstep << "    " 
	   << sphdata[i].u0/sphdata[i].dudt << endl;
      cout << " dt_courant : " << sphdata[i].h/sphdata[i].sound << "   " 
	   << sphdata[i].u0/(sphdata[i].dudt + small_number) << endl;
      string message = "Problem with energy integration (0)";
      ExceptionHandler::getIstance().raise(message);
    }
  }
  // --------------------------------------------------------------------------

  return;
}
 


//=============================================================================
//  EnergyGodunovIntegration::CorrectionTerms
/// Empty function (no corrections needed)
//=============================================================================
template <int ndim>
void EnergyGodunovIntegration<ndim>::EnergyCorrectionTerms
(int n, int Nsph, SphParticle<ndim> *sph, FLOAT timestep)
{
  return;
}



//=============================================================================
//  EnergyGodunovIntegration::EndTimestep
/// Record all important thermal quantities at the end of the step for the 
/// start of the new timestep.
//=============================================================================
template <int ndim>
void EnergyGodunovIntegration<ndim>::EndTimestep
(int n,                             ///< [in] Integer time in block time struct
 int Nsph,                          ///< [in] No. of SPH particles
 SphParticle<ndim> *sphdata)        ///< [inout] SPH particle data array
{
  int dn;                           // Integer time since beginning of step
  int i;                            // Particle counter
  int nstep;                        // Particle (integer) step size

  debug2("[EnergyGodunovIntegration::EndTimestep]");

  // --------------------------------------------------------------------------
#pragma omp parallel for default(shared) private(nstep)
  for (i=0; i<Nsph; i++) {
    dn = n - sphdata[i].nlast;
    nstep = sphdata[i].nstep;
    if (n == sphdata[i].nlast) {
      sphdata[i].u0 = sphdata[i].u;
      sphdata[i].dudt0 = sphdata[i].dudt;
    }
    //cout << "CHECKING TEND : " << i << "   " << nstep << "    " << dn << "    " << n << "    " << sphdata[i].nlast << "    " << sphdata[i].u << "    " << sphdata[i].dudt << endl;
  }
  // --------------------------------------------------------------------------

  return;
}



//=============================================================================
//  EnergyGodunovIntegration::Timestep
/// Compute explicit timestep such that u cannot change by a large fraction 
/// in one step, i.e. dt = const*u/|dudt + epsilon| 
/// where epsilon is to prevent the denominator becoming zero.
//=============================================================================
template <int ndim>
DOUBLE EnergyGodunovIntegration<ndim>::Timestep(SphParticle<ndim> &part)
{
  return this->energy_mult*(DOUBLE) (part.u/(fabs(part.dudt) + small_number));
}



template class EnergyGodunovIntegration<1>;
template class EnergyGodunovIntegration<2>;
template class EnergyGodunovIntegration<3>;
