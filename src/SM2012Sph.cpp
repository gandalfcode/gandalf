// ============================================================================
// SM2012Sph.cpp
// Contains all functions for calculating SPH quantities using the method 
// proposed by Saitoh & Makino (2012) in the conservative formalation 
// of Hopkins (2012).
// ============================================================================


#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <math.h>
#include "Precision.h"
#include "Exception.h"
#include "Sph.h"
#include "SphKernel.h"
#include "SphParticle.h"
#include "Parameters.h"
#include "EOS.h"
#include "InlineFuncs.h"
using namespace std;



// ============================================================================
// SM2012Sph::SM2012Sph
// ============================================================================
template <int ndim, template<int> class kernelclass>
SM2012Sph<ndim, kernelclass >::SM2012Sph(int hydro_forces_aux,
	    int self_gravity_aux, FLOAT alpha_visc_aux, FLOAT beta_visc_aux,
	    FLOAT h_fac_aux, FLOAT h_converge_aux, aviscenum avisc_aux,
	    acondenum acond_aux, string gas_eos_aux, string KernelName):
  Sph<ndim>(hydro_forces_aux,
		    self_gravity_aux, alpha_visc_aux, beta_visc_aux,
		    h_fac_aux, h_converge_aux, avisc_aux,
		    acond_aux, gas_eos_aux, KernelName),
  kern(kernelclass<ndim>(KernelName))
{

  this->kernp = &kern;
  this->kernfac = (FLOAT) 1.0;
  this->kernfacsqd = (FLOAT) 1.0;
}



// ============================================================================
// SM2012Sph::~SM2012Sph
// ============================================================================
template <int ndim, template<int> class kernelclass>
SM2012Sph<ndim, kernelclass >::~SM2012Sph()
{
}



// ============================================================================
// SM2012Sph::ComputeH
// Compute the value of the smoothing length of particle 'i' by iterating  
// the relation : h = h_fac*(m/rho)^(1/ndim).
// Uses the previous value of h as a starting guess and then uses either 
// a Newton-Rhapson solver, or fixed-point iteration, to converge on the 
// correct value of h.  The maximum tolerance used for deciding whether the 
// iteration has converged is given by the 'h_converge' parameter.
// ============================================================================
template <int ndim, template<int> class kernelclass>
int SM2012Sph<ndim, kernelclass >::ComputeH
(int i,                                 // id of particle
 int Nneib,                             // No. of potential neighbours
 FLOAT *m,                              // Array of neib. masses
 FLOAT *mu,                             // Array of neib. internal energies
 FLOAT *drsqd,                          // Array of neib. distances (squared)
 SphParticle<ndim> &parti)                    // Particle i data
{
  int j;                                // Neighbour id
  int jj;                               // Aux. neighbour counter
  int k;                                // Dimension counter
  int iteration = 0;                    // h-rho iteration counter
  int iteration_max = 30;               // Max. no of iterations
  FLOAT h_max = big_number;             // Max. allowed value of h
  FLOAT h_lower_bound = 0.0;            // Lower bound on h
  FLOAT h_upper_bound = big_number;     // Upper bound on h
  FLOAT invhsqd;                        // (1 / h)^2
  FLOAT ssqd;                           // Kernel parameter squared, (r/h)^2


  // Main smoothing length iteration loop
  // ==========================================================================
  do {

    // Initialise all variables for this value of h
    iteration++;
    parti.invh = (FLOAT) 1.0/parti.h;
    parti.rho = (FLOAT) 0.0;
    parti.invomega = (FLOAT) 0.0;
    parti.q = (FLOAT) 0.0;
    parti.hfactor = pow(parti.invh,ndim);
    invhsqd = parti.invh*parti.invh;

    // Loop over all nearest neighbours in list to calculate 
    // density, omega and zeta.
    // ------------------------------------------------------------------------
    for (j=0; j<Nneib; j++) {
      ssqd = drsqd[j]*invhsqd;
      parti.rho += m[j]*parti.hfactor*kern.w0_s2(ssqd);
      parti.invomega += m[j]*parti.hfactor*parti.invh*kern.womega_s2(ssqd);
      parti.q += mu[j]*parti.hfactor*kern.w0_s2(ssqd);
    }
    // ------------------------------------------------------------------------

    if (parti.rho > (FLOAT) 0.0) parti.invrho = (FLOAT) 1.0/parti.rho;

    // If h changes below some fixed tolerance, exit iteration loop
    if (parti.rho > (FLOAT) 0.0 && parti.h > h_lower_bound &&
    		fabs(parti.h - h_fac*pow(parti.m*parti.invrho,
    				invndim)) < h_converge) break;

    // Use fixed-point iteration, i.e. h_new = h_fac*(m/rho_old)^(1/ndim), 
    // for now.  If this does not converge in a reasonable number of 
    // iterations (iteration_max), then assume something is wrong and switch 
    // to a bisection method, which should be guaranteed to converge, 
    // albeit much more slowly.  (N.B. will implement Newton-Raphson soon)
    // ------------------------------------------------------------------------
    if (iteration < iteration_max)
      parti.h = h_fac*pow(parti.m*parti.invrho,invndim);

    else if (iteration == iteration_max)
      parti.h = (FLOAT) 0.5*(h_lower_bound + h_upper_bound);

    else if (iteration < 5*iteration_max) {
      if (parti.rho < small_number ||
	  parti.rho*pow(parti.h,ndim) > pow(h_fac,ndim)*parti.m)
        h_upper_bound = parti.h;
      else 
        h_lower_bound = parti.h;
      parti.h = (FLOAT) 0.5*(h_lower_bound + h_upper_bound);
    }

    else {
      string message = "Problem with convergence of h-rho iteration";
      ExceptionHandler::getIstance().raise(message);
    }

    // If the smoothing length is too large for the neighbour list, exit 
    // routine and flag neighbour list error in order to generate a larger
    // neighbour list (not properly implemented yet).
    if (parti.h > h_max) return 0;
    
  } while (parti.h > h_lower_bound && parti.h < h_upper_bound);
  // ==========================================================================


  // Normalise all SPH sums correctly
  parti.h = h_fac*pow(parti.m*parti.invrho,invndim);
  parti.invh = (FLOAT) 1.0/parti.h;
  parti.invomega = (FLOAT) 1.0 + invndim*parti.h*parti.invomega*parti.invrho;
  parti.invomega = (FLOAT) 1.0/parti.invomega;
  parti.invq = (FLOAT) 1.0/parti.q;

  // Set important thermal variables here
  parti.u = eos->SpecificInternalEnergy(parti);
  parti.sound = eos->SoundSpeed(parti);
  parti.hfactor = pow(parti.invh,ndim+1);
  parti.pfactor = eos->Pressure(parti)*parti.invrho*parti.invq;
  parti.div_v = (FLOAT) 0.0;
  parti.dudt = (FLOAT) 0.0;
  

  return 1;
}



// ============================================================================
// SM2012Sph::ComputeSphNeibForces
// Compute SPH neighbour force pairs for 
// (i) All neighbour interactions of particle i with id j > i,
// (ii) Active neighbour interactions of particle j with id j > i
// (iii) All inactive neighbour interactions of particle i with id j < i.
// This ensures that all particle-particle pair interactions are only 
// computed once only for efficiency.
// ============================================================================
template <int ndim, template<int> class kernelclass>
void SM2012Sph<ndim, kernelclass >::ComputeSphNeibForces
(int i,                                 // id of particle
 int Nneib,                             // No. of neighbours in neibpart array
 int *neiblist,                         // id of gather neighbour in neibpart
 FLOAT *drmag,                          // Distances of gather neighbours
 FLOAT *invdrmag,                       // Inverse distances of gather neibs
 FLOAT *dr,                             // Position vector of gather neibs
 SphParticle<ndim> &parti,                    // Particle i data
 SphParticle<ndim> *neibpart)                 // Neighbour particle data
{
  int j;                                // Neighbour list id
  int jj;                               // Aux. neighbour counter
  int k;                                // Dimension counter
  FLOAT draux[ndim];                 // Relative position vector
  FLOAT dv[ndim];                    // Relative velocity vector
  FLOAT dvdr;                           // Dot product of dv and dr
  FLOAT wkerni;                         // Value of w1 kernel function
  FLOAT wkernj;                         // Value of w1 kernel function
  FLOAT vsignal;                        // Signal velocity
  FLOAT paux;                           // Aux. pressure force variable
  FLOAT uaux;                           // Aux. internal energy variable
  FLOAT winvrho;                        // 0.5*(wkerni + wkernj)*invrhomean


  // Compute hydro forces
  // ==========================================================================
  if (hydro_forces == 1) {

    // Loop over all potential neighbours in the list
    // ------------------------------------------------------------------------
    for (jj=0; jj<Nneib; jj++) {
      j = neiblist[jj];
      wkerni = parti.hfactor*kern.w1(drmag[jj]*parti.invh);
      wkernj = neibpart[j].hfactor*kern.w1(drmag[jj]*neibpart[j].invh);

      for (k=0; k<ndim; k++) draux[k] = dr[jj*ndim + k];
      for (k=0; k<ndim; k++) dv[k] = neibpart[j].v[k] - parti.v[k];
      dvdr = DotProduct(dv,draux,ndim);

      // Add contribution to velocity divergence
      parti.div_v -= neibpart[j].m*dvdr*wkerni;
      neibpart[j].div_v -= parti.m*dvdr*wkernj;

      // Main SPH pressure force term
      paux = (FLOAT) 0.5*(eos->gamma - 1.0)*parti.u*neibpart[j].u*
    		  (parti.invq + neibpart[j].invq)*(wkerni + wkernj);

      
      // Add dissipation terms (for approaching particle pairs)
      // ----------------------------------------------------------------------
      if (dvdr < (FLOAT) 0.0) {

    	winvrho = (FLOAT) 0.25*(wkerni + wkernj)*(parti.invrho + neibpart[j].invrho);
	
        // Artificial viscosity term
        if (avisc == mon97) {
          vsignal = parti.sound + neibpart[j].sound - beta_visc*dvdr;
          paux -= (FLOAT) alpha_visc*vsignal*dvdr*winvrho;
          uaux = (FLOAT) 0.5*alpha_visc*vsignal*dvdr*dvdr*winvrho;
          parti.dudt -= neibpart[j].m*uaux;
          neibpart[j].dudt -= parti.m*uaux;
        }

        // Artificial conductivity term
        if (acond == wadsley2008) {
	      uaux = (FLOAT) 0.5*dvdr*(neibpart[j].u - parti.u)*
	    		  (parti.invrho*wkerni + neibpart[j].invrho*wkernj);
	      parti.dudt += neibpart[j].m*uaux;
	      neibpart[j].dudt -= parti.m*uaux;
        }
        else if (acond == price2008) {
    	  vsignal = sqrt(fabs(eos->Pressure(parti) -
			      eos->Pressure(neibpart[j]))*0.5*
    			  (parti.invrho + neibpart[j].invrho));
          parti.dudt += 0.5*neibpart[j].m*vsignal*
            (parti.u - neibpart[j].u)*winvrho;
          neibpart[j].dudt -= 0.5*parti.m*vsignal*
            (parti.u - neibpart[j].u)*winvrho;
        }
	
      }
      // ----------------------------------------------------------------------

      uaux = 0.5*neibpart[j].m*neibpart[j].u*dvdr*(wkerni + wkernj);

      // Add total hydro contribution to acceleration for particle i
      for (k=0; k<ndim; k++) parti.a[k] += neibpart[j].m*draux[k]*paux;
      parti.dudt += 0.5*neibpart[j].m*neibpart[j].u*dvdr*(wkerni + wkernj)*
    		  parti.pfactor;
      
      // If neighbour is also active, add contribution to force here
      for (k=0; k<ndim; k++) neibpart[j].a[k] -= parti.m*draux[k]*paux;
      neibpart[j].dudt += 0.5*parti.m*parti.u*dvdr*(wkerni + wkernj)*
    		  neibpart[j].pfactor;

    }
    // ------------------------------------------------------------------------

  }
  // ==========================================================================


  return;
}



// ============================================================================
// SM2012Sph::ComputeSphNeibDudt
// Empty definition (
// ============================================================================
template <int ndim, template<int> class kernelclass>
void SM2012Sph<ndim, kernelclass >::ComputeSphNeibDudt
(int i, int Nneib, int *neiblist, FLOAT *drmag, 
 FLOAT *invdrmag, FLOAT *dr, SphParticle<ndim> &parti, SphParticle<ndim> *neibpart)
{
  return;
}



// ============================================================================
// SM2012Sph::ComputeSphDerivatives
// Empty definition
// ============================================================================
template <int ndim, template<int> class kernelclass>
void SM2012Sph<ndim, kernelclass >::ComputeSphDerivatives
(int i, int Nneib, int *neiblist, FLOAT *drmag, 
 FLOAT *invdrmag, FLOAT *dr, SphParticle<ndim> &parti, SphParticle<ndim> *neibpart)
{
  return;
}



// ============================================================================
// SM2012Sph::ComputePostHydroQuantities
// ..
// ============================================================================
template <int ndim, template<int> class kernelclass>
void SM2012Sph<ndim, kernelclass >::ComputePostHydroQuantities
(SphParticle<ndim> &parti)
{
  parti.div_v *= parti.invrho;
  return;
}



// ============================================================================
// SM2012Sph::ComputeGravForces
// Compute the contribution to the total gravitational force of particle 'i' 
// due to 'Nneib' neighbouring particles in the list 'neiblist'.
// ============================================================================
template <int ndim, template<int> class kernelclass>
void SM2012Sph<ndim, kernelclass >::ComputeDirectGravForces
(int i,                                 // id of particle
 int Ndirect,                           // No. of nearby 'gather' neighbours
 int *directlist,                       // id of gather neighbour in neibpart
 SphParticle<ndim> &parti,                    // Particle i data
 SphParticle<ndim> *sph)                      // Neighbour particle data
{
  return;
}



template class SM2012Sph<1, M4Kernel>;
template class SM2012Sph<1, QuinticKernel>;
template class SM2012Sph<1, GaussianKernel>;
template class SM2012Sph<1, TabulatedKernel>;
template class SM2012Sph<2, M4Kernel>;
template class SM2012Sph<2, QuinticKernel>;
template class SM2012Sph<2, GaussianKernel>;
template class SM2012Sph<2, TabulatedKernel>;
template class SM2012Sph<3, M4Kernel>;
template class SM2012Sph<3, QuinticKernel>;
template class SM2012Sph<3, GaussianKernel>;
template class SM2012Sph<3, TabulatedKernel>;
