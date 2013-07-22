//=============================================================================
//  NbodyHermite4TS.cpp
//  Contains functions for integrating star particle positions and velocities 
//  using the 4th-order Hermite scheme (Makino & Aarseth 1992) using
//  time-symmetric iterations (Hut et al. 1995??).
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
//  NbodyHermite4TS::NbodyHermite4TS()
/// N-body 4th-order Hermite class constructor
//=============================================================================
template <int ndim, template<int> class kernelclass>
NbodyHermite4TS<ndim, kernelclass>::NbodyHermite4TS
(int nbody_softening_aux, int sub_systems_aux, 
 DOUBLE nbody_mult_aux, string KernelName, int Npec) :
  NbodyHermite4<ndim, kernelclass>(nbody_softening_aux, sub_systems_aux,
                      nbody_mult_aux, KernelName, Npec)
{
}



//=============================================================================
//  NbodyHermite4TS::~NbodyHermite4TS()
/// N-body 4th-order Hermite class destructor
//=============================================================================
template <int ndim, template<int> class kernelclass>
NbodyHermite4TS<ndim, kernelclass>::~NbodyHermite4TS()
{
}



//=============================================================================
//  NbodyHermite4TS::CorrectionTerms
/// Compute 2nd and 3rd time derivatives of the acceleration using Hermite 
/// interpolation.  Finally correct positions to 5th order and velocities to 
/// 4th order using higher-order derivatives.
//=============================================================================
template <int ndim, template<int> class kernelclass>
void NbodyHermite4TS<ndim, kernelclass>::CorrectionTerms
(int n,                             ///< Integer time
 int N,                             ///< No. of stars/systems
 NbodyParticle<ndim> **star,        ///< Main star/system array
 DOUBLE timestep)                   ///< Smallest timestep value
{
  int dn;                           // Integer time since beginning of step
  int i;                            // Particle counter
  int k;                            // Dimension counter
  int nstep;                        // Particle (integer) step size
  DOUBLE dt;                        // Physical time step size
  DOUBLE invdt;                     // 1 / dt

  debug2("[NbodyHermite4TS::CorrectionTerms]");

  // Loop over all system particles
  // --------------------------------------------------------------------------
  for (i=0; i<N; i++) {
    dn = n - star[i]->nlast;
    nstep = star[i]->nstep;

    if (dn == nstep) {
      dt = timestep*(DOUBLE) nstep;
      invdt = 1.0 / dt;
    
      for (k=0; k<ndim; k++) {
        star[i]->a2dot[k] = 
	  (-6.0*(star[i]->a0[k] - star[i]->a[k]) - dt*
	   (4.0*star[i]->adot0[k] + 2.0*star[i]->adot[k]))*invdt*invdt;
        star[i]->a3dot[k] = 
	  (12.0*(star[i]->a0[k] - star[i]->a[k]) + 6.0*dt*
	   (star[i]->adot0[k] + star[i]->adot[k]))*invdt*invdt*invdt;
      }

      for (k=0; k<ndim; k++) {
        star[i]->v[k] = star[i]->v0[k] 
	  + 0.5*(star[i]->a0[k] + star[i]->a[k])*dt
	  //  + 0.1*(star[i]->adot[k] - star[i]->adot0[k])*dt*dt;
	  - onetwelfth*(star[i]->adot[k] - star[i]->adot0[k])*dt*dt;
        star[i]->r[k] = star[i]->r0[k] 
	  + 0.5*(star[i]->v0[k] + star[i]->v[k])*dt
          //+ 0.1*(star[i]->a[k] - star[i]->a0[k])*dt*dt;
          - onetwelfth*(star[i]->a[k] - star[i]->a0[k])*dt*dt;
      }
    }
    
  }
  // --------------------------------------------------------------------------

  return;
}



//=============================================================================
//  NbodyHermite4TS::IntegrateInternalMotion
/// This function integrates the internal motion of a system. First integrates
/// the internal motion of its sub-systems by recursively calling their method,
/// then integrates the COM of the sub-systems.
//=============================================================================
template <int ndim, template<int> class kernelclass>
void NbodyHermite4TS<ndim, kernelclass>::IntegrateInternalMotion
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
  for (k=0; k<ndim; k++) adotcom[k] = 0.0;

  // Make local copies of children and calculate COM properties
  // --------------------------------------------------------------------------
  for (i=0; i<Nchildren; i++) {
    for (k=0; k<ndim; k++) rcom[k] += children[i]->m*children[i]->r[k];
    for (k=0; k<ndim; k++) vcom[k] += children[i]->m*children[i]->v[k];
    for (k=0; k<ndim; k++) acom[k] += children[i]->m*children[i]->a[k];
    for (k=0; k<ndim; k++) adotcom[k] += children[i]->m*children[i]->adot[k];
  }

  // Normalise COM values
  for (k=0; k<ndim; k++) rcom[k] /= systemi->m;
  for (k=0; k<ndim; k++) vcom[k] /= systemi->m;
  for (k=0; k<ndim; k++) acom[k] /= systemi->m;
  for (k=0; k<ndim; k++) adotcom[k] /= systemi->m;


  // Now convert to COM frame
  // --------------------------------------------------------------------------
  for (i=0; i<Nchildren; i++) {
    for (k=0; k<ndim; k++) children[i]->r[k] -= rcom[k];
    for (k=0; k<ndim; k++) children[i]->r0[k] -= rcom[k];
    for (k=0; k<ndim; k++) children[i]->v[k] -= vcom[k];
    for (k=0; k<ndim; k++) children[i]->v0[k] -= vcom[k];
    for (k=0; k<ndim; k++) children[i]->a[k] -= acom[k];
    for (k=0; k<ndim; k++) children[i]->a0[k] -= acom[k];
    for (k=0; k<ndim; k++) children[i]->adot[k] -= adotcom[k];
    for (k=0; k<ndim; k++) children[i]->adot0[k] -= adotcom[k];
    for (k=0; k<ndim; k++) children[i]->a2dot[k] -= 0.0;
    for (k=0; k<ndim; k++) children[i]->a3dot[k] -= 0.0;
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

    // Time-symmetric iteration loop
    // ------------------------------------------------------------------------
    for (it=0; it<Npec; it++) {

      //Zero all acceleration terms
      for (i=0; i<Nchildren; i++) {
        children[i]->gpot = 0.0;
	children[i]->gpe_internal = 0.0;
        for (k=0; k<ndim; k++) children[i]->a[k] = 0.0;
        for (k=0; k<ndim; k++) children[i]->adot[k] = 0.0;
        for (k=0; k<ndim; k++) children[i]->a2dot[k] = 0.0;
        for (k=0; k<ndim; k++) children[i]->a3dot[k] = 0.0;
      }

      // Calculate forces, derivatives and other terms
      CalculateDirectGravForces(Nchildren, children);
      
      // Apply correction terms
      CorrectionTerms(nlocal_steps, Nchildren, children, dt);
    }
    // ------------------------------------------------------------------------

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
    for (k=0; k<ndim; k++) children[i]->adot[k] += systemi->adot[k];
    for (k=0; k<ndim; k++) children[i]->adot0[k] += systemi->adot[k];
    children[i]->gpot = children[i]->gpot + systemi->gpot;
  }

  return;
}




// Template class instances for each dimensionality value (1, 2 and 3) and
// employed kernel (M4, Quintic, Gaussian and tabulated).
template class NbodyHermite4TS<1, M4Kernel>;
template class NbodyHermite4TS<1, QuinticKernel>;
template class NbodyHermite4TS<1, GaussianKernel>;
template class NbodyHermite4TS<1, TabulatedKernel>;
template class NbodyHermite4TS<2, M4Kernel>;
template class NbodyHermite4TS<2, QuinticKernel>;
template class NbodyHermite4TS<2, GaussianKernel>;
template class NbodyHermite4TS<2, TabulatedKernel>;
template class NbodyHermite4TS<3, M4Kernel>;
template class NbodyHermite4TS<3, QuinticKernel>;
template class NbodyHermite4TS<3, GaussianKernel>;
template class NbodyHermite4TS<3, TabulatedKernel>;
