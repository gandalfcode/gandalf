//=============================================================================
//  NbodyHermite4.cpp
//  Contains functions for integrating star particle positions and velocities 
//  using the 4th-order Hermite scheme (Makino & Aarseth 1992).
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
//  NbodyHermite4::NbodyHermite4()
/// N-body 4th-order Hermite class constructor
//=============================================================================
template <int ndim, template<int> class kernelclass>
NbodyHermite4<ndim, kernelclass>::NbodyHermite4
(int nbody_softening_aux, int sub_systems_aux, 
 DOUBLE nbody_mult_aux, string KernelName) : 
  Nbody<ndim>(nbody_softening_aux, sub_systems_aux, 
              nbody_mult_aux, KernelName),
  kern(kernelclass<ndim>(KernelName))
{
}



//=============================================================================
//  NbodyHermite4::~NbodyHermite4()
/// N-body 4th-order Hermite class destructor
//=============================================================================
template <int ndim, template<int> class kernelclass>
NbodyHermite4<ndim, kernelclass>::~NbodyHermite4()
{
}



//=============================================================================
//  NbodyHermite4::CalculateDirectGravForces
/// Calculate all star-star force contributions for active systems using 
/// direct summation with unsoftened gravity.
//=============================================================================
template <int ndim, template<int> class kernelclass>
void NbodyHermite4<ndim, kernelclass>::CalculateDirectGravForces
(int N,                             ///< Number of stars
 NbodyParticle<ndim> **star)        ///< Array of stars/systems
{
  int i,j,k;                        // Star and dimension counters
  DOUBLE dr[ndim];                  // Relative position vector
  DOUBLE drdt;                      // Rate of change of distance
  DOUBLE drsqd;                     // Distance squared
  DOUBLE dv[ndim];                  // Relative velocity vector
  DOUBLE invdrmag;                  // 1 / drmag

  debug2("[NbodyHermite4::CalculateDirectGravForces]");

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
      for (k=0; k<ndim; k++) dv[k] = star[j]->v[k] - star[i]->v[k];
      drsqd = DotProduct(dr,dr,ndim);
      invdrmag = 1.0/sqrt(drsqd);
      drdt = DotProduct(dv,dr,ndim)*invdrmag;

      star[i]->gpot -= star[j]->m*invdrmag;
      for (k=0; k<ndim; k++) star[i]->a[k] += star[j]->m*dr[k]*pow(invdrmag,3);
      for (k=0; k<ndim; k++) star[i]->adot[k] +=
	star[j]->m*pow(invdrmag,3)*(dv[k] - 3.0*drdt*invdrmag*dr[k]);

    }
    // ------------------------------------------------------------------------

  }
  // --------------------------------------------------------------------------

  return;
}



//=============================================================================
//  NbodyHermite4::AdvanceParticles
/// Integrate particle positions to 3nd order, and particle velocities to 2nd
/// order from the beginning of the step to the current simulation time, i.e. 
/// r(t+dt) = r(t) + v(t)*dt + a(t)*dt^2/2 + adot(t)*dt^3/6, 
/// v(t+dt) = v(t) + a(t)*dt + adot(t)*dt^2/2.
/// Also set particles at the end of step as 'active' in order to compute 
/// the end-of-step force computation.
//=============================================================================
template <int ndim, template<int> class kernelclass>
void NbodyHermite4<ndim, kernelclass>::AdvanceParticles
(int n,                             ///< Integer time
 int N,                             ///< No. of stars/systems
 NbodyParticle<ndim> **star,        ///< Main star/system array
 DOUBLE timestep)                   ///< Smallest timestep value
{
  int i;                            // Particle counter
  int k;                            // Dimension counter
  int nstep;                        // Particle (integer) step size
  DOUBLE dt;                        // Timestep since start of step

  debug2("[NbodyHermite4::AdvanceParticles]");

  // Advance all systems
  // --------------------------------------------------------------------------
  for (i=0; i<Nsystem; i++) {

    // Compute time since beginning of step
    nstep = star[i]->nstep;
    if (n%nstep == 0) dt = timestep*(DOUBLE) nstep;
    else dt = timestep*(DOUBLE) (n%nstep);

    // Advance positions to third order and velocities to second order
    for (k=0; k<ndim; k++) star[i]->r[k] = star[i]->r0[k] +
      star[i]->v0[k]*dt + 0.5*star[i]->a0[k]*dt*dt +
      onesixth*star[i]->adot0[k]*dt*dt*dt;
    for (k=0; k<vdim; k++) star[i]->v[k] = star[i]->v0[k] +
      star[i]->a0[k]*dt + 0.5*star[i]->adot0[k]*dt*dt;

    // If at end of step, set system particle as active
    if (n%nstep == 0) star[i]->active = true;
  }
  // --------------------------------------------------------------------------

  return;
}
 


//=============================================================================
//  NbodyHermite4::CorrectionTerms
/// Compute 2nd and 3rd time derivatives of the acceleration using Hermite 
/// interpolation.  Finally correct positions to 5th order and velocities to 
/// 4th order using higher-order derivatives.
//=============================================================================
template <int ndim, template<int> class kernelclass>
void NbodyHermite4<ndim, kernelclass>::CorrectionTerms
(int n,                             ///< Integer time
 int N,                             ///< No. of stars/systems
 NbodyParticle<ndim> **star,        ///< Main star/system array
 DOUBLE timestep)                   ///< Smallest timestep value
{
  int i;                            // Particle counter
  int k;                            // Dimension counter
  int nstep;                        // Particle (integer) step size
  DOUBLE dt;                        // Physical time step size

  debug2("[NbodyHermite4::CorrectionTerms]");

  // Loop over all system particles
  // --------------------------------------------------------------------------
  for (i=0; i<N; i++) {
    nstep = star[i]->nstep;
    if (n%nstep == 0) continue;
    dt = timestep*(DOUBLE) nstep;
    
    for (k=0; k<ndim; k++) {
      star[i]->a2dot[k] = (-6.0*(star[i]->a0[k] - star[i]->a[k]) - dt*
			  (4.0*star[i]->adot0[k] + 2.0*star[i]->adot[k]))/dt/dt;
      
      star[i]->a3dot[k] = (12.0*(star[i]->a0[k] - star[i]->a[k]) + 6.0*dt*
			  (star[i]->adot0[k] + star[i]->adot[k]))/dt/dt/dt;
      star[i]->r[k] += star[i]->a2dot[k]*dt*dt*dt*dt/24.0 +
	star[i]->a3dot[k]*dt*dt*dt*dt*dt/120.0;
      star[i]->v[k] += star[i]->a2dot[k]*dt*dt*dt/6.0 +
	star[i]->a3dot[k]*dt*dt*dt*dt/24.0;
    }
    
  }
  // --------------------------------------------------------------------------

  return;
}



//=============================================================================
//  NbodyHermite4::EndTimestep
/// Record all important star particle quantities at the end of the step 
/// for the start of the new timestep.
//=============================================================================
template <int ndim, template<int> class kernelclass>
void NbodyHermite4<ndim, kernelclass>::EndTimestep
(int n,                             ///< Integer time
 int N,                             ///< No. of stars/systems
 NbodyParticle<ndim> **star)        ///< Main star/system array
{
  int i;                            // Particle counter
  int k;                            // Dimension counter
  int nstep;                        // Particle (integer) step size

  debug2("[NbodyHermite4::EndTimestep]");

  // Loop over all system particles
  // --------------------------------------------------------------------------
  for (i=0; i<Nsystem; i++) {
    nstep = star[i]->nstep;

    // If at end of the current step, set quantites for start of new step
    if (n%nstep == 0) {
      for (k=0; k<ndim; k++) star[i]->r0[k] = star[i]->r[k];
      for (k=0; k<ndim; k++) star[i]->v0[k] = star[i]->v[k];
      for (k=0; k<ndim; k++) star[i]->a0[k] = star[i]->a[k];
      for (k=0; k<ndim; k++) star[i]->adot0[k] = star[i]->adot[k];
      star[i]->active = false;
    }

  }
  // --------------------------------------------------------------------------

  return;
}



//=============================================================================
//  NbodyHermite4::Timestep
/// Default timestep size for SPH particles.  Takes the minimum of : 
/// (i)  const*h/(sound_speed + h*|div_v|)    (Courant condition)
/// (ii) const*sqrt(h/|a|)                    (Acceleration condition)
//=============================================================================
template <int ndim, template<int> class kernelclass>
DOUBLE NbodyHermite4<ndim, kernelclass>::Timestep
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
template class NbodyHermite4<1, M4Kernel>;
template class NbodyHermite4<1, QuinticKernel>;
template class NbodyHermite4<1, GaussianKernel>;
template class NbodyHermite4<1, TabulatedKernel>;
template class NbodyHermite4<2, M4Kernel>;
template class NbodyHermite4<2, QuinticKernel>;
template class NbodyHermite4<2, GaussianKernel>;
template class NbodyHermite4<2, TabulatedKernel>;
template class NbodyHermite4<3, M4Kernel>;
template class NbodyHermite4<3, QuinticKernel>;
template class NbodyHermite4<3, GaussianKernel>;
template class NbodyHermite4<3, TabulatedKernel>;

