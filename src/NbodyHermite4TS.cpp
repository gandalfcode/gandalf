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
  int i;                            // Particle counter
  int k;                            // Dimension counter
  int nstep;                        // Particle (integer) step size
  DOUBLE dt;                        // Physical time step size

  debug2("[NbodyHermite4TS::CorrectionTerms]");

  // Loop over all system particles
  // --------------------------------------------------------------------------
  for (i=0; i<N; i++) {
    nstep = star[i]->nstep;
    
    if (n%nstep == 0) {
      dt = timestep*(DOUBLE) nstep;

      for (k=0; k<ndim; k++) {
        star[i]->a2dot[k] = (-6.0*(star[i]->a0[k] - star[i]->a[k]) - dt*
		  	    (4.0*star[i]->adot0[k] + 2.0*star[i]->adot[k]))/dt/dt;
        star[i]->a3dot[k] = (12.0*(star[i]->a0[k] - star[i]->a[k]) + 6.0*dt*
			    (star[i]->adot0[k] + star[i]->adot[k]))/dt/dt/dt;

        star[i]->v[k] = star[i]->v0[k] + 0.5*(star[i]->a0[k] + star[i]->a[k])*dt
          + onetwelfth*(star[i]->adot[k] - star[i]->adot0[k])*dt*dt;
        star[i]->r[k] = star[i]->r0[k] + 0.5*(star[i]->v0[k] + star[i]->v[k])*dt
          + onetwelfth*(star[i]->a[k] - star[i]->a0[k])*dt*dt;
      }
    }
    
  }
  // --------------------------------------------------------------------------

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
