//=============================================================================
//  NbodyLeapfrogKDK.cpp
//  Contains functions for integrating star particle positions and velocities 
//  using the leapfrog kick-drift-kick (KDK) scheme.
//=============================================================================


#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <math.h>
#include "Precision.h"
#include "NbodyParticle.h"
#include "StarParticle.h"
#include "Parameters.h"
#include "Nbody.h"
#include "SphKernel.h"
#include "Debug.h"
#include "Exception.h"
#include "InlineFuncs.h"
using namespace std;



//=============================================================================
//  NbodyLeapfrogKDK::NbodyLeapfrogKDK()
/// N-body leapfrog KDK class constructor
//=============================================================================
template <int ndim, template<int> class kernelclass>
NbodyLeapfrogKDK<ndim, kernelclass>::NbodyLeapfrogKDK
(int nbody_softening_aux, int sub_systems_aux, 
 DOUBLE nbody_mult_aux, string KernelName) : 
  Nbody<ndim>(nbody_softening_aux, sub_systems_aux, 
              nbody_mult_aux, KernelName),
  kern(kernelclass<ndim>(KernelName))
{
}



//=============================================================================
//  NbodyLeapfrogKDK::~NbodyLeapfrog()
/// N-body leapfrog KDK class destructor
//=============================================================================
template <int ndim, template<int> class kernelclass>
NbodyLeapfrogKDK<ndim, kernelclass>::~NbodyLeapfrogKDK()
{
}



//=============================================================================
//  NbodyLeapfrogKDK::CalculateDirectGravForces
/// Calculate all star-star force contributions for active systems using 
/// direct summation with unsoftened gravity.
//=============================================================================
template <int ndim, template<int> class kernelclass>
void NbodyLeapfrogKDK<ndim, kernelclass>::CalculateDirectGravForces
(int N,                             ///< Number of stars
 NbodyParticle<ndim> **star)        ///< Array of stars/systems
{
  int i,j,k;                        // Star and dimension counters
  DOUBLE dr[ndim];                  // Relative position vector
  DOUBLE drsqd;                     // Distance squared
  DOUBLE invdrmag;                  // 1 / drmag

  debug2("[NbodyLeapfrogKDK::CalculateDirectGravForces]");

  // Loop over all (active) stars
  // --------------------------------------------------------------------------
  for (i=0; i<N; i++) {

    if (star[i]->active == 0) continue;

    star[i]->gpot = 0.0;
    for (k=0; k<ndim; k++) star[i]->a[k] = 0.0;
    for (k=0; k<ndim; k++) star[i]->adot[k] = 0.0;

    // Sum grav. contributions for all other stars (excluding star itself)
    // ------------------------------------------------------------------------
    for (j=0; j<N; j++) {
      if (i == j) continue;

      for (k=0; k<ndim; k++) dr[k] = star[j]->r[k] - star[i]->r[k];
      drsqd = DotProduct(dr,dr,ndim);
      invdrmag = 1.0/sqrt(drsqd);

      star[i]->gpot -= star[j]->m*invdrmag;
      for (k=0; k<ndim; k++) star[i]->a[k] += star[j]->m*dr[k]*pow(invdrmag,3);

    }
    // ------------------------------------------------------------------------

  }
  // --------------------------------------------------------------------------

  return;
}



//=============================================================================
//  NbodyLeapfrogKDK::CalculateDirectSPHForces
/// Calculate all ..
//=============================================================================
template <int ndim, template<int> class kernelclass>
void NbodyLeapfrogKDK<ndim, kernelclass>::CalculateDirectSPHForces
(int N,                             ///< Number of stars
 int Ngas,                          ///< Number of gas particles
 SphParticle<ndim> *sphdata,        ///< Array of SPH particles
 NbodyParticle<ndim> **star)        ///< Array of stars/systems
{
  int i,j,k;                        // Star and dimension counters
  DOUBLE dr[ndim];                  // Relative position vector
  DOUBLE drmag;                     // ..
  DOUBLE drsqd;                     // Distance squared
  DOUBLE invdrmag;                  // 1 / drmag
  DOUBLE paux;                      // ..
  DOUBLE gaux;                      // ..

  debug2("[NbodyLeapfrogKDK::CalculateDirectSPHForces]");

  // Loop over all (active) stars
  // --------------------------------------------------------------------------
  for (i=0; i<N; i++) {

    if (star[i]->active == 0) continue;

    // Sum grav. contributions for all other stars (excluding star itself)
    // ------------------------------------------------------------------------
    for (j=0; j<Ngas; j++) {

      for (k=0; k<ndim; k++) dr[k] = sphdata[j].r[k] - star[i]->r[k];
      drsqd = DotProduct(dr,dr,ndim);
      invdrmag = 1.0/sqrt(drsqd);

      paux = sphdata[j].invh*sphdata[j].invh*kern.wgrav(drmag*sphdata[j].invh)
	+ star[i]->invh*star[i]->invh*kern.wgrav(drmag*star[i]->invh);
      gaux = sphdata[j].invh*kern.wpot(drmag*sphdata[j].invh) + 
	star[i]->invh*kern.wpot(drmag*star[i]->invh);

      // Add ..
      for (k=0; k<ndim; k++) star[i]->a[k] += 0.5*sphdata[j].m*dr[k]*paux;
      star[i]->gpot += 0.5*sphdata[j].m*gaux;

    }
    // ------------------------------------------------------------------------

  }
  // --------------------------------------------------------------------------

  return;
}



//=============================================================================
//  NbodyLeapfrogKDK::AdvanceParticles
/// Integrate star positions to 2nd order, and star velocities to 1st
/// order from the beginning of the step to the current simulation time, i.e. 
/// r(t+dt) = r(t) + v(t)*dt + 0.5*a(t)*dt^2, 
/// v(t+dt) = v(t) + a(t)*dt.
/// Also set particles at the end of step as 'active' in order to compute 
/// the end-of-step force computation and velocity correction step.
//=============================================================================
template <int ndim, template<int> class kernelclass>
void NbodyLeapfrogKDK<ndim, kernelclass>::AdvanceParticles
(int n,                             ///< Integer time
 int N,                             ///< No. of stars/systems
 NbodyParticle<ndim> **star,        ///< Main star/system array
 DOUBLE timestep)                   ///< Smallest timestep value
{
  int i;                            // Particle counter
  int k;                            // Dimension counter
  int nstep;                        // Particle (integer) step size
  DOUBLE dt;                        // Timestep since start of step

  debug2("[NbodyLeapfrogKDK::AdvanceParticles]");

  // --------------------------------------------------------------------------
  for (i=0; i<N; i++) {

    // Compute time since beginning of step
    nstep = star[i]->nstep;
    if (n%nstep == 0) dt = timestep*(DOUBLE) nstep;
    else dt = timestep*(DOUBLE) (n%nstep);

    // Advance positions to second order and velocities to first order
    for (k=0; k<ndim; k++) star[i]->r[k] = star[i]->r0[k] +
			     star[i]->v0[k]*dt + 0.5*star[i]->a0[k]*dt*dt;
    for (k=0; k<vdim; k++)
      star[i]->v[k] = star[i]->v0[k] + star[i]->a0[k]*dt;

    // If at end of step, set system particle as active
    if (n%nstep == 0) star[i]->active = true;
  }
  // --------------------------------------------------------------------------

  return;
}
 


//=============================================================================
//  NbodyLeapfrogKDK::CorrectionTerms
/// Compute velocity integration to second order at the end of the step by 
/// adding a second order correction term, 
/// v(t+dt) -> v(t+dt) + 0.5*(a(t+dt) - a(t))*dt.
/// The full integration therefore becomes
/// v(t+dt) = v(t) + 0.5*(a(t) + a(t+dt))*dt.
//=============================================================================
template <int ndim, template<int> class kernelclass>
void NbodyLeapfrogKDK<ndim, kernelclass>::CorrectionTerms
(int n,                             ///< Integer time
 int N,                             ///< No. of stars/systems
 NbodyParticle<ndim> **star,        ///< Main star/system array
 DOUBLE timestep)                   ///< Smallest timestep value
{
  int i;                            // Particle counter
  int k;                            // Dimension counter
  int nstep;                        // Particle (integer) step size

  debug2("[NbodyLeapfrogKDK::CorrectionTerms]");

  // --------------------------------------------------------------------------
  for (i=0; i<N; i++) {
    nstep = star[i]->nstep;
    if (n%nstep == 0)
      for (k=0; k<ndim; k++) star[i]->v[k] +=
	0.5*(star[i]->a[k] - star[i]->a0[k])*timestep*(DOUBLE) nstep;
  }
  // --------------------------------------------------------------------------

  return;
}



//=============================================================================
//  NbodyLeapfrogKDK::EndTimestep
/// Record all important star particle quantities at the end of the step 
/// for the start of the new timestep.
//=============================================================================
template <int ndim, template<int> class kernelclass>
void NbodyLeapfrogKDK<ndim, kernelclass>::EndTimestep
(int n,                             ///< Integer time
 int N,                             ///< No. of stars/systems
 NbodyParticle<ndim> **star)        ///< Main star/system array
{
  int i;                            // Particle counter
  int k;                            // Dimension counter
  int nstep;                        // Particle (integer) step size

  debug2("[NbodyLeapfrogKDK::EndTimestep]");

  // --------------------------------------------------------------------------
  for (i=0; i<N; i++) {
    nstep = star[i]->nstep;
    //cout << "HERE : " << i << "   " << N << "   " << nstep << endl;
    //cout << "POS  : " << star[i]->r[0] << "    " << star[i]->r[1] << endl;
    if (n%nstep == 0) {
      for (k=0; k<ndim; k++) star[i]->r0[k] = star[i]->r[k];
      for (k=0; k<ndim; k++) star[i]->v0[k] = star[i]->v[k];
      for (k=0; k<ndim; k++) star[i]->a0[k] = star[i]->a[k];
      star[i]->active = false;
    }
  }
  // --------------------------------------------------------------------------

  return;
}



//=============================================================================
//  NbodyLeapfrogKDK::Timestep
/// Default timestep size for SPH particles.  Takes the minimum of : 
/// (i)  const*h/(sound_speed + h*|div_v|)    (Courant condition)
/// (ii) const*sqrt(h/|a|)                    (Acceleration condition)
//=============================================================================
template <int ndim, template<int> class kernelclass>
DOUBLE NbodyLeapfrogKDK<ndim, kernelclass>::Timestep
(NbodyParticle<ndim> *star)             ///< Reference to SPH particle
{
  DOUBLE timestep;                      // Minimum value of particle timesteps
  DOUBLE amag;                          // Magnitude of particle acceleration

  // Acceleration condition
  amag = sqrt(DotProduct(star->a,star->a,ndim));
  timestep = nbody_mult*sqrt(star->h/(amag + small_number_dp));

  //cout << "TIMESTEP : " << amag << "   " << star.h << "   " << timestep << endl;

  return timestep;
}



// Template class instances for each dimensionality value (1, 2 and 3)
template class NbodyLeapfrogKDK<1, M4Kernel>;
template class NbodyLeapfrogKDK<1, QuinticKernel>;
template class NbodyLeapfrogKDK<1, GaussianKernel>;
template class NbodyLeapfrogKDK<1, TabulatedKernel>;
template class NbodyLeapfrogKDK<2, M4Kernel>;
template class NbodyLeapfrogKDK<2, QuinticKernel>;
template class NbodyLeapfrogKDK<2, GaussianKernel>;
template class NbodyLeapfrogKDK<2, TabulatedKernel>;
template class NbodyLeapfrogKDK<3, M4Kernel>;
template class NbodyLeapfrogKDK<3, QuinticKernel>;
template class NbodyLeapfrogKDK<3, GaussianKernel>;
template class NbodyLeapfrogKDK<3, TabulatedKernel>;
