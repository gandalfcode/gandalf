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
              nbody_mult_aux, KernelName, 1),
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

    //star[i]->gpot = 0.0;
    //for (k=0; k<ndim; k++) star[i]->a[k] = 0.0;
    //for (k=0; k<ndim; k++) star[i]->adot[k] = 0.0;

    // Sum grav. contributions for all other stars (excluding star itself)
    // ------------------------------------------------------------------------
    for (j=0; j<N; j++) {
      if (i == j) continue;

      for (k=0; k<ndim; k++) dr[k] = star[j]->r[k] - star[i]->r[k];
      drsqd = DotProduct(dr,dr,ndim);
      invdrmag = 1.0/sqrt(drsqd);

      // Add contribution to main star array
      for (k=0; k<ndim; k++) star[i]->a[k] += star[j]->m*dr[k]*pow(invdrmag,3);
      star[i]->gpot += star[j]->m*invdrmag;

    }
    // ------------------------------------------------------------------------

  }
  // --------------------------------------------------------------------------

  return;
}



//=============================================================================
//  NbodyLeapfrogKDK::CalculateDirectSPHForces
/// Calculate all SPH forces by direct summation.
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
  DOUBLE drmag;                     // Distance
  DOUBLE drsqd;                     // Distance squared
  DOUBLE invhmean;                  // ..
  DOUBLE invdrmag;                  // 1 / drmag
  DOUBLE paux;                      // Aux. force variable

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
      drmag = sqrt(drsqd);
      invdrmag = 1.0/drmag;
      invhmean = 2.0/(star[i]->h + sphdata[j].h);

      paux = sphdata[j].m*invhmean*invhmean*
        kern.wgrav(drmag*invhmean)*invdrmag;

      // Add contribution to main star array
      for (k=0; k<ndim; k++) star[i]->a[k] += dr[k]*paux;
      star[i]->gpot += sphdata[j].m*invhmean*kern.wpot(drmag*invhmean);

    }
    // ------------------------------------------------------------------------

  }
  // --------------------------------------------------------------------------

  return;
}



//=============================================================================
//  NbodyLeapfrogKDK::CalculateAllStartupQuantities
/// Empty function for Leapfrog-KDK (no additional quantities required).
//=============================================================================
template <int ndim, template<int> class kernelclass>
void NbodyLeapfrogKDK<ndim, kernelclass>::CalculateAllStartupQuantities
(int N,                             ///< Number of stars
 NbodyParticle<ndim> **star)        ///< Array of stars/systems
{
  return;
}



//=============================================================================
//  NbodyLeapfrogKDK::AdvanceParticles
/// Integrate star positions to 2nd order, and star velocities to 1st
/// order from the beginning of the step to the current simulation time, i.e. 
/// $r(t+dt) = r(t) + v(t)*dt + 0.5*a(t)*dt^2$, 
/// $v(t+dt) = v(t) + a(t)*dt$.
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

  // Advance positions and velocities of all star particles
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
/// $v(t+dt) -> v(t+dt) + 0.5*(a(t+dt) - a(t))*dt$.
/// The full integration therefore becomes
/// $v(t+dt) = v(t) + 0.5*(a(t) + a(t+dt))*dt$.
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

  // Loop over all star particles and calculate correction terms only for 
  // those at end of step.
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

  // Loop over all star particles and set values for those at end of step
  // --------------------------------------------------------------------------
  for (i=0; i<N; i++) {
    nstep = star[i]->nstep;
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
(NbodyParticle<ndim> *star)         ///< Reference to SPH particle
{
  DOUBLE timestep;                  // Minimum value of particle timesteps
  DOUBLE amag;                      // Magnitude of particle acceleration

  // Acceleration condition
  amag = sqrt(DotProduct(star->a,star->a,ndim));
  timestep = nbody_mult*sqrt(star->h/(amag + small_number_dp));

  timestep = min(timestep,star->dt_internal);

  return timestep;
}




//=============================================================================
//  NbodyLeapfrogKDK::IntegrateInternalMotion
/// This function integrates the internal motion of a system. First integrates
/// the internal motion of its sub-systems by recursively calling their method,
/// then integrates the COM of the sub-systems.
//=============================================================================
template <int ndim, template<int> class kernelclass>
void NbodyLeapfrogKDK<ndim, kernelclass>::IntegrateInternalMotion
(SystemParticle<ndim>* systemi,     ///< [inout] System that we wish to 
                                    ///<         integrate the internal motion
 DOUBLE tlocal_end)                 ///< [in]    Time to integrate the 
                                    ///<         internal motion for.
{
  int i;                                              // ..
  int it;                                             // Iteration counter
  int k;                                              // ..
  int Nchildren = systemi->Nchildren;                 // No. of child systems
  int nlocal_steps = 0;                               // ..
  DOUBLE dt;                                          // ..
  DOUBLE tlocal=0.0;                                  // ..
  DOUBLE rcom[ndim];                                  // ..
  DOUBLE vcom[ndim];                                  // ..
  DOUBLE acom[ndim];                                  // ..
  DOUBLE adotcom[ndim];                               // ..
  NbodyParticle<ndim>** children = systemi->children; // Child systems

  // Zero all COM summation variables
  for (k=0; k<ndim; k++) rcom[k] = 0.0;
  for (k=0; k<ndim; k++) vcom[k] = 0.0;
  for (k=0; k<ndim; k++) acom[k] = 0.0;

  // Make local copies of children and calculate COM properties
  // --------------------------------------------------------------------------
  for (i=0; i<Nchildren; i++) {
    for (k=0; k<ndim; k++) rcom[k] += children[i]->m*children[i]->r[k];
    for (k=0; k<ndim; k++) vcom[k] += children[i]->m*children[i]->v[k];
    for (k=0; k<ndim; k++) acom[k] += children[i]->m*children[i]->a[k];
  }

  // Normalise COM values
  for (k=0; k<ndim; k++) rcom[k] /= systemi->m;
  for (k=0; k<ndim; k++) vcom[k] /= systemi->m;
  for (k=0; k<ndim; k++) acom[k] /= systemi->m;


  // Now convert to COM frame
  // --------------------------------------------------------------------------
  for (i=0; i<Nchildren; i++) {
    for (k=0; k<ndim; k++) children[i]->r[k] -= rcom[k];
    for (k=0; k<ndim; k++) children[i]->r0[k] -= rcom[k];
    for (k=0; k<ndim; k++) children[i]->v[k] -= vcom[k];
    for (k=0; k<ndim; k++) children[i]->v0[k] -= vcom[k];
    for (k=0; k<ndim; k++) children[i]->a[k] -= acom[k];
    for (k=0; k<ndim; k++) children[i]->a0[k] -= acom[k];
    children[i]->active = true;
    children[i]->nstep = 1;
    children[i]->level = 0;
  }


  // Main time integration loop
  // ==========================================================================
  do {

    // Calculate time-step
    dt = std::min(big_number, tlocal_end - tlocal);
    for (i=0; i<Nchildren; i++) {
      dt = std::min(dt, Timestep(children[i]));
    }
    tlocal += dt;
    nlocal_steps +=1;

    // Advance position and velocities
    AdvanceParticles(nlocal_steps, Nchildren, children, dt);

    //Zero all acceleration terms
    for (i=0; i<Nchildren; i++) {
      children[i]->gpot = 0.0;
      children[i]->gpe_internal = 0.0;
      for (k=0; k<ndim; k++) children[i]->a[k] = 0.0;
    }
    
    // Calculate forces, derivatives and other terms
    CalculateDirectGravForces(Nchildren, children);
    
    // Apply correction terms
    CorrectionTerms(nlocal_steps, Nchildren, children, dt);

    // Now loop over children and, if they are systems, integrate
    // their internal motion
    // ------------------------------------------------------------------------
    for (i=0; i<Nchildren; i++) {

      if (children[i]->Ncomp > 1)
	// The cast is needed because the function is defined only in
	// SystemParticle, not in NbodyParticle.  
	// The safety of the cast relies on the correctness of the Ncomp value
	IntegrateInternalMotion(static_cast<SystemParticle<ndim>* > 
				(children[i]), dt);
    }

    // Set end-of-step variables
    EndTimestep(nlocal_steps, Nchildren, children);

  } while (tlocal < tlocal_end);
  // ==========================================================================


  // Copy children back to main coordinate system
  // --------------------------------------------------------------------------
  for (i=0; i<Nchildren; i++) {
    for (k=0; k<ndim; k++) children[i]->r[k] += systemi->r[k];
    for (k=0; k<ndim; k++) children[i]->r0[k] += systemi->r[k];
    for (k=0; k<ndim; k++) children[i]->v[k] += systemi->v[k];
    for (k=0; k<ndim; k++) children[i]->v0[k] += systemi->v[k];
    for (k=0; k<ndim; k++) children[i]->a[k] += systemi->a[k];
    for (k=0; k<ndim; k++) children[i]->a0[k] += systemi->a[k];
    children[i]->gpot = children[i]->gpot + systemi->gpot;
  }

  return;
}



// Template class instances for each dimensionality value (1, 2 and 3) and 
// employed kernel (M4, Quintic, Gaussian and tabulated).
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
