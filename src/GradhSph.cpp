//=============================================================================
//  GradhSph.cpp
//  Contains all functions for calculating conservative 'grad-h' SPH quantities
//  (See Springel & Hernquist (2002) and Price & Monaghan (2007).
//=============================================================================


#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <math.h>
#include "Precision.h"
#include "Sph.h"
#include "SphParticle.h"
#include "Parameters.h"
#include "SphKernel.h"
#include "EOS.h"
#include "Debug.h"
#include "Exception.h"
#include "InlineFuncs.h"
using namespace std;



//=============================================================================
//  GradhSph::GradhSph
/// GradhSph class constructor
//=============================================================================
template <int ndim, template<int> class kernelclass>
GradhSph<ndim, kernelclass>::GradhSph(int hydro_forces_aux,
	    int self_gravity_aux, FLOAT alpha_visc_aux, FLOAT beta_visc_aux,
	    FLOAT h_fac_aux, FLOAT h_converge_aux, aviscenum avisc_aux,
	    acondenum acond_aux, string gas_eos_aux, string KernelName):
  Sph<ndim>(hydro_forces_aux,self_gravity_aux, alpha_visc_aux, beta_visc_aux,
	    h_fac_aux, h_converge_aux, avisc_aux,acond_aux, gas_eos_aux, 
            KernelName),
  kern(kernelclass<ndim>(KernelName))
{
  this->kernp = &kern;
  this->kernfac = (FLOAT) 1.0;
  this->kernfacsqd = (FLOAT) 1.0;
}



//=============================================================================
//  GradhSph::~GradhSph
/// GradhSph class destructor
//=============================================================================
template <int ndim, template<int> class kernelclass>
GradhSph<ndim, kernelclass>::~GradhSph()
{
}



//=============================================================================
//  GradhSph::ComputeH
/// Compute the value of the smoothing length of particle 'i' by iterating  
/// the relation : h = h_fac*(m/rho)^(1/ndim).
/// Uses the previous value of h as a starting guess and then uses either 
/// a Newton-Rhapson solver, or fixed-point iteration, to converge on the 
/// correct value of h.  The maximum tolerance used for deciding whether the 
/// iteration has converged is given by the 'h_converge' parameter.
//=============================================================================
template <int ndim, template<int> class kernelclass>
int GradhSph<ndim, kernelclass>::ComputeH
(int i,                             ///< [in] id of particle
 int Nneib,                         ///< [in] No. of potential neighbours
 FLOAT *m,                          ///< [in] Array of neib. masses
 FLOAT *mu,                         ///< [in] Array of m*u (not needed here)
 FLOAT *drsqd,                      ///< [in] Array of neib. distances squared
 SphParticle<ndim> &parti)          ///< [inout] Particle i data
{
  int j;                            // Neighbour id
  int jj;                           // Aux. neighbour counter
  int k;                            // Dimension counter
  int iteration = 0;                // h-rho iteration counter
  int iteration_max = 30;           // Max. no of iterations
  FLOAT h_max = big_number;         // Max. allowed value of h
  FLOAT h_lower_bound = 0.0;        // Lower bound on h
  FLOAT h_upper_bound = big_number; // Upper bound on h
  FLOAT invhsqd;                    // (1 / h)^2
  FLOAT ssqd;                       // Kernel parameter squared, (r/h)^2


  // Main smoothing length iteration loop
  // ==========================================================================
  do {

    // Initialise all variables for this value of h
    iteration++;
    parti.invh = (FLOAT) 1.0/parti.h;
    parti.rho = (FLOAT) 0.0;
    parti.invomega = (FLOAT) 0.0;
    parti.zeta = (FLOAT) 0.0;
    parti.hfactor = pow(parti.invh,ndim);
    invhsqd = parti.invh*parti.invh;

    // Loop over all nearest neighbours in list to calculate 
    // density, omega and zeta.
    // ------------------------------------------------------------------------
    for (j=0; j<Nneib; j++) {
      ssqd = drsqd[j]*invhsqd;
      parti.rho += m[j]*parti.hfactor*kern.w0_s2(ssqd);
      parti.invomega += m[j]*parti.hfactor*parti.invh*kern.womega_s2(ssqd);
      parti.zeta += m[j]*invhsqd*kern.wzeta_s2(ssqd);
    }
    // ------------------------------------------------------------------------

    if (parti.rho > (FLOAT) 0.0) parti.invrho = (FLOAT) 1.0/parti.rho;

    // If h changes below some fixed tolerance, exit iteration loop
    if (parti.rho > (FLOAT) 0.0 && parti.h > h_lower_bound &&
    		fabs(parti.h - h_fac*pow(parti.m*parti.invrho,
    				Sph<ndim>::invndim)) < h_converge) break;

    // Use fixed-point iteration, i.e. h_new = h_fac*(m/rho_old)^(1/ndim), 
    // for now.  If this does not converge in a reasonable number of 
    // iterations (iteration_max), then assume something is wrong and switch 
    // to a bisection method, which should be guaranteed to converge, 
    // albeit much more slowly.  (N.B. will implement Newton-Raphson soon)
    // ------------------------------------------------------------------------
    if (iteration < iteration_max)
      parti.h = h_fac*pow(parti.m*parti.invrho,Sph<ndim>::invndim);

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
      cout << "HERE : " << parti.h << "    " << parti.rho << "    " 
	   << h_upper_bound << "     " << h_lower_bound << "    " 
	   << parti.hfactor << "     " 
	   << parti.m*parti.hfactor*kern.w0(0.0) << endl;
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
  //parti.h = h_fac*pow(parti.m*parti.invrho,Sph<ndim>::invndim);
  parti.invh = (FLOAT) 1.0/parti.h;
  parti.invomega = (FLOAT) 1.0 + 
    Sph<ndim>::invndim*parti.h*parti.invomega*parti.invrho;
  parti.invomega = (FLOAT) 1.0/parti.invomega;
  parti.zeta = -Sph<ndim>::invndim*parti.h*parti.zeta*
    parti.invrho*parti.invomega;
  //parti.invomega = 1.0;
  //parti.zeta = 0.0;

  // Set important thermal variables here
  parti.u = eos->SpecificInternalEnergy(parti);
  parti.sound = eos->SoundSpeed(parti);
  parti.hfactor = pow(parti.invh,ndim+1);
  parti.pfactor = eos->Pressure(parti)*parti.invrho*
    parti.invrho*parti.invomega;
  parti.div_v = (FLOAT) 0.0;
  
  return 1;
}



//=============================================================================
//  GradhSph::ComputeSphNeibForces
/// Compute SPH neighbour force pairs for 
/// (i) All neighbour interactions of particle i with i.d. j > i,
/// (ii) Active neighbour interactions of particle j with i.d. j > i
/// (iii) All inactive neighbour interactions of particle i with i.d. j < i.
/// This ensures that all particle-particle pair interactions are only 
/// computed once only for efficiency.
//=============================================================================
template <int ndim, template<int> class kernelclass>
void GradhSph<ndim, kernelclass>::ComputeSphNeibForces
(int i,                             ///< [in] id of particle
 int Nneib,                         ///< [in] No. of neins in neibpart array
 int *neiblist,                     ///< [in] id of gather neibs in neibpart
 FLOAT *drmag,                      ///< [in] Distances of gather neighbours
 FLOAT *invdrmag,                   ///< [in] Inverse distances of gather neibs
 FLOAT *dr,                         ///< [in] Position vector of gather neibs
 SphParticle<ndim> &parti,          ///< [inout] Particle i data
 SphParticle<ndim> *neibpart)       ///< [inout] Neighbour particle data
{
  int j;                            // Neighbour list id
  int jj;                           // Aux. neighbour counter
  int k;                            // Dimension counter
  FLOAT draux[ndim];                // Relative position vector
  FLOAT dv[ndim];                   // Relative velocity vector
  FLOAT dvdr;                       // Dot product of dv and dr
  FLOAT wkerni;                     // Value of w1 kernel function
  FLOAT wkernj;                     // Value of w1 kernel function
  FLOAT vsignal;                    // Signal velocity
  FLOAT paux;                       // Aux. pressure force variable
  FLOAT uaux;                       // Aux. internal energy variable
  FLOAT winvrho;                    // 0.5*(wkerni + wkernj)*invrhomean


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
      paux = parti.pfactor*wkerni + neibpart[j].pfactor*wkernj;
      
      // Add dissipation terms (for approaching particle pairs)
      // ----------------------------------------------------------------------
      if (dvdr < (FLOAT) 0.0) {

    	winvrho = (FLOAT) 0.25*(wkerni + wkernj)*
	  (parti.invrho + neibpart[j].invrho);
	
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

      // Add total hydro contribution to acceleration for particle i
      for (k=0; k<ndim; k++) parti.a[k] += neibpart[j].m*draux[k]*paux;
      
      // If neighbour is also active, add contribution to force here
      for (k=0; k<ndim; k++) neibpart[j].a[k] -= parti.m*draux[k]*paux;

    }
    // ------------------------------------------------------------------------

  }
  // ==========================================================================

   // Compute gravitational contribution
    // ------------------------------------------------------------------------
    /*if (self_gravity == 1) {
      paux = (FLOAT) 0.5*
	(parti.invh*parti.invh*kern.wgrav(drmag[j]*parti.invh) +
	 parti.zeta*wkern - invdrmag[j]*invdrmag[j]);
      for (k=0; k<ndim; k++) parti.agrav[k] += neibpart[j].m*draux[k]*paux;
      parti.gpot += neibpart[j].m*parti.invh*wpot(drmag[j].parti.invh);
      }*/

  return;
}



//=============================================================================
//  GradhSph::ComputeSphNeibGravForces
/// Compute SPH neighbour force pairs for 
/// (i) All neighbour interactions of particle i with i.d. j > i,
/// (ii) Active neighbour interactions of particle j with i.d. j > i
/// (iii) All inactive neighbour interactions of particle i with i.d. j < i.
/// This ensures that all particle-particle pair interactions are only 
/// computed once only for efficiency.
//=============================================================================
template <int ndim, template<int> class kernelclass>
void GradhSph<ndim, kernelclass>::ComputeSphNeibGravForces
(int i,                             ///< [in] id of particle
 int Nneib,                         ///< [in] No. of neins in neibpart array
 int *neiblist,                     ///< [in] id of gather neibs in neibpart
 FLOAT *drmag,                      ///< [in] Distances of gather neighbours
 FLOAT *invdrmag,                   ///< [in] Inverse distances of gather neibs
 FLOAT *dr,                         ///< [in] Position vector of gather neibs
 SphParticle<ndim> &parti,          ///< [inout] Particle i data
 SphParticle<ndim> *neibpart)       ///< [inout] Neighbour particle data
{
  int j;                            // Neighbour list id
  int jj;                           // Aux. neighbour counter
  int k;                            // Dimension counter
  FLOAT draux[ndim];                // Relative position vector
  FLOAT dv[ndim];                   // Relative velocity vector
  FLOAT dvdr;                       // Dot product of dv and dr
  FLOAT paux;                       // Aux. pressure force variable
  FLOAT gaux;                       // Aux. internal energy variable


  // Compute hydro forces
  // ==========================================================================
  if (self_gravity == 1) {

    // Loop over all potential neighbours in the list
    // ------------------------------------------------------------------------
    for (jj=0; jj<Nneib; jj++) {
      j = neiblist[jj];
      for (k=0; k<ndim; k++) draux[k] = dr[jj*ndim + k];

      // Main SPH gravity terms
      paux = parti.invh*parti.invh*kern.wgrav(drmag[jj]*parti.invh) + 
	parti.zeta*parti.hfactor*kern.w1(drmag[jj]*parti.invh) + 
	neibpart[j].invh*neibpart[j].invh*
	kern.wgrav(drmag[jj]*neibpart[j].invh) + 
	neibpart[j].zeta*neibpart[j].hfactor*
	kern.w1(drmag[jj]*neibpart[j].invh);
      gaux = (parti.invh*kern.wpot(drmag[jj]*parti.invh) + 
	      neibpart[j].invh*kern.wpot(drmag[jj]*neibpart[j].invh));

      // Add total hydro contribution to acceleration for particle i
      for (k=0; k<ndim; k++) parti.agrav[k] += 0.5*neibpart[j].m*draux[k]*paux;
      parti.gpot += 0.5*neibpart[j].m*gaux;

      // If neighbour is also active, add contribution to force here
      //for (k=0; k<ndim; k++) neibpart[j].agrav[k] -= 0.5*parti.m*draux[k]*paux;
      //neibpart[j].gpot += 0.5*parti.m*gaux;

    }
    // ------------------------------------------------------------------------

  }
  // ==========================================================================


  return;
}



//=============================================================================
//  GradhSph::ComputeSphNeibDudt
/// Empty definition (not needed for grad-h SPH)
//=============================================================================
template <int ndim, template<int> class kernelclass>
void GradhSph<ndim, kernelclass >::ComputeSphNeibDudt
(int i, int Nneib, int *neiblist, FLOAT *drmag, FLOAT *invdrmag, 
 FLOAT *dr, SphParticle<ndim> &parti, SphParticle<ndim> *neibpart)
{
  return;
}



//=============================================================================
//  GradhSph::ComputeSphDerivatives
/// Empty definition
//=============================================================================
template <int ndim, template<int> class kernelclass>
void GradhSph<ndim, kernelclass>::ComputeSphDerivatives
(int i, int Nneib, int *neiblist, FLOAT *drmag, FLOAT *invdrmag, 
 FLOAT *dr, SphParticle<ndim> &parti, SphParticle<ndim> *neibpart)
{
  return;
}



//=============================================================================
//  GradhSph::ComputePostHydroQuantities
/// ..
//=============================================================================
template <int ndim, template<int> class kernelclass>
void GradhSph<ndim, kernelclass>::ComputePostHydroQuantities
(SphParticle<ndim> &parti)
{
  parti.div_v *= parti.invrho;
  parti.dudt -= eos->Pressure(parti)*parti.div_v*parti.invrho*parti.invomega;

  return;
}



//=============================================================================
//  GradhSph::ComputeGravForces
/// Compute the contribution to the total gravitational force of particle 'i' 
/// due to 'Nneib' neighbouring particles in the list 'neiblist'.
//=============================================================================
template <int ndim, template<int> class kernelclass>
void GradhSph<ndim, kernelclass>::ComputeDirectGravForces
(int i,                             // id of particle
 int Ndirect,                       // No. of nearby 'gather' neighbours
 int *directlist,                   // id of gather neighbour in neibpart
 SphParticle<ndim> &parti,          // Particle i data
 SphParticle<ndim> *sph)            // Neighbour particle data
{
  int j;                            // ..
  int jj;                           // ..
  int k;                            // ..
  FLOAT dr[ndim];                   // ..
  FLOAT drsqd;                      // ..
  FLOAT invdrmag;                   // ..

  // Loop over all neighbouring particles in list
  // --------------------------------------------------------------------------
  for (jj=0; jj<Ndirect; jj++) {
    j = directlist[jj];

    for (k=0; k<ndim; k++) dr[k] = sph[j].r[k] - parti.r[k];
    drsqd = DotProduct(dr,dr,ndim);
    invdrmag = 1.0/(sqrt(drsqd) + small_number);

    parti.gpot += sph[j].m*invdrmag;
    for (k=0; k<ndim; k++) 
      parti.agrav[k] += sph[j].m*dr[k]*pow(invdrmag,3);

  }
  // --------------------------------------------------------------------------

  return;
}



//=============================================================================
//  GradhSph::ComputeStarGravForces
/// ..
//=============================================================================
template <int ndim, template<int> class kernelclass>
void GradhSph<ndim, kernelclass>::ComputeStarGravForces
(int N,
 NbodyParticle<ndim> **nbodydata,
 SphParticle<ndim> &parti)
{
  int j;
  int k;
  FLOAT dr[ndim];
  FLOAT drmag;
  FLOAT drsqd;
  FLOAT invdrmag;
  FLOAT paux;
  FLOAT gaux;

  // --------------------------------------------------------------------------
  for (j=0; j<N; j++) {

    for (k=0; k<ndim; k++) dr[k] = nbodydata[j]->r[k] - parti.r[k];
    drsqd = DotProduct(dr,dr,ndim);
    drmag = sqrt(drsqd);

    paux = parti.invh*parti.invh*kern.wgrav(drmag*parti.invh) + 
      nbodydata[j]->invh*nbodydata[j]->invh*
      kern.wgrav(drmag*nbodydata[j]->invh);
    gaux = (parti.invh*kern.wpot(drmag*parti.invh) + 
	    nbodydata[j]->invh*kern.wpot(drmag*nbodydata[j]->invh));

    // Add total hydro contribution to acceleration for particle i
    for (k=0; k<ndim; k++) parti.agrav[k] += 0.5*nbodydata[j]->m*dr[k]*paux;
    parti.gpot += 0.5*nbodydata[j]->m*gaux;
    
  }
  // --------------------------------------------------------------------------

  return;
}



template class GradhSph<1, M4Kernel>;
template class GradhSph<1, QuinticKernel>;
template class GradhSph<1, GaussianKernel>;
template class GradhSph<1, TabulatedKernel>;
template class GradhSph<2, M4Kernel>;
template class GradhSph<2, QuinticKernel>;
template class GradhSph<2, GaussianKernel>;
template class GradhSph<2, TabulatedKernel>;
template class GradhSph<3, M4Kernel>;
template class GradhSph<3, QuinticKernel>;
template class GradhSph<3, GaussianKernel>;
template class GradhSph<3, TabulatedKernel>;
