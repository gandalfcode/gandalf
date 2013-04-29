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
#include "Nbody.h"
#include "Debug.h"
#include "InlineFuncs.h"
using namespace std;



//=============================================================================
//  NbodyHermite4::NbodyHermite4()
/// ..
//=============================================================================
template <int ndim, template<int> class kernelclass>
NbodyHermite4<ndim, kernelclass>::NbodyHermite4(int nbody_softening_aux, int sub_systems_aux, DOUBLE nbody_mult_aux, string KernelName) : 
  Nbody<ndim>(nbody_softening_aux, sub_systems_aux, nbody_mult_aux, KernelName),
  kern(kernelclass<ndim>(KernelName))
{
}



//=============================================================================
//  NbodyHermite4::~NbodyHermite4()
/// ..
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
void NbodyHermite4<ndim, kernelclass>::CalculateDirectGravForces(void)
{
  int i,j,k;                        // Star and dimension counters
  DOUBLE dr[ndim];                  // Relative position vector
  DOUBLE drsqd;                     // Distance squared
  DOUBLE invdrmag;                  // 1 / drmag
  StarParticle<ndim> stari;         // Local copy of star particle

  debug2("[NbodyHermite4::CalculateDirectGravForces]");

  // Loop over all (active) stars
  // --------------------------------------------------------------------------
  for (i=0; i<Nstar; i++) {
    if (stardata[i].active == 0) continue;

    stari = stardata[i];
    stari.gpot = 0.0;
    for (k=0; k<ndim; k++) stari.a[k] = 0.0;
    for (k=0; k<ndim; k++) stari.adot[k] = 0.0;

    // Sum grav. contributions for all other stars (excluding star itself)
    // ------------------------------------------------------------------------
    for (j=0; j<Nstar; j++) {
      if (i == j) continue;

      for (k=0; k<ndim; k++) dr[k] = stardata[j].r[k] - stari.r[k];
      drsqd = DotProduct(dr,dr,ndim);
      invdrmag = 1.0/sqrt(drsqd);

      stari.gpot -= stardata[j].m*invdrmag;
      for (k=0; k<ndim; k++) stari.a[k] += stari.m*dr[k]*pow(invdrmag,3);

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
void NbodyHermite4<ndim, kernelclass>::AdvanceParticles(int n, int level_step, int Nsystem,
				     StarParticle<ndim> *system, DOUBLE timestep)
{
  int i;                                // Particle counter
  int k;                                // Dimension counter
  int nstep;                            // Particle (integer) step size
  DOUBLE dt;                            // Timestep since start of step

  debug2("[NbodyHermite4::AdvanceParticles]");

  // Advance all systems
  // --------------------------------------------------------------------------
  for (i=0; i<Nsystem; i++) {

    // Compute time since beginning of step
    nstep = pow(2,level_step - system[i].level);
    if (n%nstep == 0) dt = timestep*(DOUBLE) nstep;
    else dt = timestep*(DOUBLE) (n%nstep);

    // Advance positions to third order and velocities to second order
    for (k=0; k<ndim; k++) system[i].r[k] = system[i].r0[k] + 
      system[i].v0[k]*dt + 0.5*system[i].a0[k]*dt*dt + 
      onesixth*system[i].adot0[k]*dt*dt*dt;
    for (k=0; k<vdim; k++) system[i].v[k] = system[i].v0[k] + 
      system[i].a0[k]*dt + 0.5*system[i].adot0[k]*dt*dt;

    // If at end of step, set system particle as active
    if (n%nstep == 0) system[i].active = true;
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
void NbodyHermite4<ndim, kernelclass>::CorrectionTerms(int n, int level_step, int Nsystem,
				    StarParticle<ndim> *system, DOUBLE timestep)
{
  int i;                                // Particle counter
  int k;                                // Dimension counter
  int nstep;                            // Particle (integer) step size

  debug2("[NbodyHermite4::CorrectionTerms]");

  // Loop over all system particles
  // --------------------------------------------------------------------------
  for (i=0; i<Nsystem; i++) {
    nstep = pow(2,level_step - system[i].level);

    // If at end of step, compute correction terms
    if (n%nstep == 0) {
      for (k=0; k<ndim; k++) {

	system[i].adot2[k] = 0.0;
	system[i].adot3[k] = 0.0;

	system[i].r[k] += 0.0;
	system[i].v[k] += 0.0;

      }
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
void NbodyHermite4<ndim, kernelclass>::EndTimestep(int n, int level_step, 
				int Nsystem, StarParticle<ndim> *system)
{
  int i;                                // Particle counter
  int k;                                // Dimension counter
  int nstep;                            // Particle (integer) step size

  debug2("[NbodyHermite4::EndTimestep]");

  // Loop over all system particles
  // --------------------------------------------------------------------------
  for (i=0; i<Nsystem; i++) {
    nstep = pow(2,level_step - system[i].level);

    // If at end of the current step, set quantites for start of new step
    if (n%nstep == 0) {
      for (k=0; k<ndim; k++) system[i].r0[k] = system[i].r[k];
      for (k=0; k<ndim; k++) system[i].v0[k] = system[i].v[k];
      for (k=0; k<ndim; k++) system[i].a0[k] = system[i].a[k];
      for (k=0; k<ndim; k++) system[i].adot0[k] = system[i].adot[k];
      system[i].active = false;
    }

  }
  // --------------------------------------------------------------------------

  return;
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

