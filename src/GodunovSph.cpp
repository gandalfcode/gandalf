// ============================================================================
// GodunovSph.cpp
// Contains all functions for calculating conservative 'grad-h' SPH quantities
// (See Springel & Hernquist (2002) and Price & Monaghan (2007).
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


static const int riemann_order = 1;
static const string riemann_solver = "hllc";


// ============================================================================
// GodunovSph::GodunovSph
// ============================================================================
template <typename kernelclass>
GodunovSph<kernelclass>::GodunovSph(int ndimaux, int vdimaux, int bdimaux, int hydro_forces_aux,
	    int self_gravity_aux, FLOAT alpha_visc_aux, FLOAT beta_visc_aux,
	    FLOAT h_fac_aux, FLOAT h_converge_aux, aviscenum avisc_aux,
	    acondenum acond_aux, string gas_eos_aux, string KernelName):
  Sph(ndimaux, vdimaux, bdimaux, hydro_forces_aux,
		    self_gravity_aux, alpha_visc_aux, beta_visc_aux,
		    h_fac_aux, h_converge_aux, avisc_aux,
		    acond_aux, gas_eos_aux, KernelName),
  kern(kernelclass(ndimaux, KernelName))
{
  allocated = false;
  Nsph = 0;
  Nsphmax = 0;
  kernp = &kern;
}



// ============================================================================
// GodunovSph::~GodunovSph
// ============================================================================
template <typename kernelclass>
GodunovSph<kernelclass>::~GodunovSph()
{
}



// ============================================================================
// GodunovSph::ComputeH
// Compute the value of the smoothing length of particle 'i' by iterating  
// the relation : h = h_fac*(m/rho)^(1/ndim).
// Uses the previous value of h as a starting guess and then uses either 
// a Newton-Rhapson solver, or fixed-point iteration, to converge on the 
// correct value of h.  The maximum tolerance used for deciding whether the 
// iteration has converged is given by the 'h_converge' parameter.
// ============================================================================
template <typename kernelclass>
int GodunovSph<kernelclass>::ComputeH
(int i,                                 // id of particle
 int Nneib,                             // No. of potential neighbours
 FLOAT *m,                              // Array of neib. masses
 FLOAT *mu,                             // Array of m*u (not needed here)
 FLOAT *drsqd,                          // Array of neib. distances (squared)
 SphParticle &parti)                    // Particle i data
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
  parti.zeta = -invndim*parti.h*parti.zeta*parti.invrho*parti.invomega;

  // Set important thermal variables here
  parti.u = eos->SpecificInternalEnergy(parti);
  parti.sound = eos->SoundSpeed(parti);
  parti.hfactor = pow(parti.invh,ndim+1);
  parti.press = eos->Pressure(parti);
  parti.pfactor = parti.press*parti.invrho*parti.invrho*parti.invomega;
  parti.div_v = (FLOAT) 0.0;


  return 1;
}



// ============================================================================
// GodunovSph::ComputeSphNeibForces
// Compute SPH neighbour force pairs for 
// (i) All neighbour interactions of particle i with i.d. j > i,
// (ii) Active neighbour interactions of particle j with i.d. j > i
// (iii) All inactive neighbour interactions of particle i with i.d. j < i.
// This ensures that all particle-particle pair interactions are only 
// computed once only for efficiency.
// ============================================================================
template <typename kernelclass>
void GodunovSph<kernelclass>::ComputeSphNeibForces
(int i,                                 // id of particle
 int Nneib,                             // No. of neighbours in neibpart array
 int *neiblist,                         // id of gather neighbour in neibpart
 FLOAT *drmag,                          // Distances of gather neighbours
 FLOAT *invdrmag,                       // Inverse distances of gather neibs
 FLOAT *dr,                             // Position vector of gather neibs
 SphParticle &parti,                    // Particle i data
 SphParticle *neibpart)                 // Neighbour particle data
{
  int j;                                // Neighbour list id
  int jj;                               // Aux. neighbour counter
  int k;                                // Dimension counter
  FLOAT draux[ndimmax];                 // Relative position vector
  FLOAT dv[ndimmax];                    // Relative velocity vector
  FLOAT dvdr;                           // Dot product of dv and dr
  FLOAT wkerni;                         // Value of w1 kernel function
  FLOAT wkernj;                         // Value of w1 kernel function
  FLOAT vsignal;                        // Signal velocity
  FLOAT paux;                           // Aux. pressure force variable
  FLOAT uaux;                           // Aux. internal energy variable
  FLOAT winvrho;                        // 0.5*(wkerni + wkernj)*invrhomean

  FLOAT Cij;                            // ..
  FLOAT Dij;                            // ..
  FLOAT Vsqdi;                          // ..
  FLOAT Vsqdj;                          // ..
  FLOAT Sij;                            // ..
  FLOAT pl,pr;
  FLOAT rhol,rhor;
  FLOAT vl,vr;
  FLOAT pstar;
  FLOAT vstar;
  FLOAT vtemp[ndimmax];
  FLOAT gradi;
  FLOAT gradj;

  static const FLOAT hconv = powf(invsqrttwo,ndim+1);


  // Compute hydro forces
  // ==========================================================================
  if (hydro_forces == 1) {

    // Loop over all potential neighbours in the list
    // ------------------------------------------------------------------------
    for (jj=0; jj<Nneib; jj++) {
      j = neiblist[jj];
      wkerni = hconv*parti.hfactor*kern.w1(invsqrttwo*drmag[jj]*parti.invh);
      wkernj = hconv*neibpart[j].hfactor*
	kern.w1(invsqrttwo*drmag[jj]*neibpart[j].invh);

      for (k=0; k<ndim; k++) draux[k] = dr[jj*ndim + k];
      for (k=0; k<ndim; k++) dv[k] = neibpart[j].v[k] - parti.v[k];
      dvdr = DotProduct(dv,draux,ndim);

      // Add contribution to velocity divergence
      parti.div_v -= neibpart[j].m*dvdr*wkerni;
      neibpart[j].div_v -= parti.m*dvdr*wkernj;

      // Linear interpolate quantites between left and right states
      Cij = (parti.invrho - neibpart[j].invrho)*invdrmag[jj];
      Dij = (FLOAT) 0.5*(parti.invrho + neibpart[j].invrho);
      Vsqdi = (FLOAT) 0.25*parti.h*parti.h*Cij*Cij + Dij*Dij;
      Vsqdj = (FLOAT) 0.25*neibpart[j].h*neibpart[j].h*Cij*Cij + Dij*Dij;
      Sij = (FLOAT) 0.25*Cij*Dij*(parti.h*parti.h/Vsqdi + 
				  neibpart[j].h*neibpart[j].h/Vsqdj);

      // Initialise the LHS and RHS of the Riemann problem
      InitialiseRiemannProblem(parti,neibpart[j],draux,drmag[jj],dvdr,
    		  parti.sound,neibpart[j].sound,pl,pr,rhol,rhor,vl,vr);

      // Now solve Riemann problem and return intermediate state variables
      if (riemann_solver == "hllc")
	HllcSolver("basic1",pl,pr,rhol,rhor,parti.sound,neibpart[j].sound,
		   vl,vr,eos->gamma,pstar,vstar);
      else if (riemann_solver == "mg")
	MgSolver("basic1",pl,pr,rhol,rhor,parti.sound,neibpart[j].sound,
		 vl,vr,eos->gamma,pstar,vstar);

      //FLOAT rhomean = 0.5*(rhol + rhor);
      //FLOAT soundmean = 0.5*(parti.sound + neibpart[j].sound);
      //pstar = 0.5*(pl + pr) - 0.5*dvdr*rhomean*soundmean;
      //vstar = 0.5*(vl + vr) - 0.5*(pr - pl)/rhomean/soundmean;

      // Main SPH pressure force term
      paux = pstar*(Vsqdi*wkerni + Vsqdj*wkernj);

      // Add total hydro contribution to acceleration for particle i
      for (k=0; k<ndim; k++) parti.a[k] += neibpart[j].m*draux[k]*paux;
      
      // If neighbour is also active, add contribution to force here
      for (k=0; k<ndim; k++) neibpart[j].a[k] -= parti.m*draux[k]*paux;

    }
    // ------------------------------------------------------------------------

  }
  // ==========================================================================


  return;
}


// ============================================================================
// GodunovSph::ComputeSphDerivatives
// ..
// ============================================================================
template <typename kernelclass>
void GodunovSph<kernelclass>::InitialiseRiemannProblem
(SphParticle partl,
 SphParticle partr,
 FLOAT draux[ndimmax],
 FLOAT drmag,
 FLOAT dvdr,
 FLOAT soundl,
 FLOAT soundr,
 FLOAT &pl,
 FLOAT &pr,
 FLOAT &rhol,
 FLOAT &rhor,
 FLOAT &vl,
 FLOAT &vr)
{
  FLOAT deltal;
  FLOAT deltar;
  FLOAT limiter;

  // ..
  pl = partl.press;
  pr = partr.press;
  rhol = partl.rho;
  rhor = partr.rho;
  vl = (FLOAT) 0.0;
  vr = dvdr;

  // --------------------------------------------------------------------------
  if (riemann_order == 2) {

	// Slope limiter for pressure gradients
	deltal = DotProduct(partl.gradP,draux,ndim)*drmag;
	deltar = DotProduct(partr.gradP,draux,ndim)*drmag;

	limiter = max3(0.0,deltal/(deltar + small_number),deltar/(deltal + small_number));
	limiter = 2.0*limiter/(limiter*limiter + 1);
	//limiter = 1.0;
	//deltal = min3(2.0*fabs(pl - pr),fabs(deltal),2.0*fabs(2.0*deltal + pr - pl))*sgn(deltal);
    //deltar = min3(2.0*fabs(pr - pl),fabs(deltar),2.0*fabs(2.0*deltar + pl - pr))*sgn(deltar);

	//if (sgn(pl - pr) != sgn(deltal) || sgn(deltal) == sgn(2.0*deltal + pr - pl)) deltal = 0.0;
	//if (sgn(pr - pl) != sgn(deltar) || sgn(deltar) == sgn(2.0*deltar + pl - pr)) deltar = 0.0;

	pl += limiter*(FLOAT) 0.5*deltal*(1.0 - partl.sound*partl.dt/drmag);
	pr -= limiter*(FLOAT) 0.5*deltar*(1.0 - partr.sound*partr.dt/drmag);


	// Slope limiter for pressure gradients
	deltal = DotProduct(partl.gradrho,draux,ndim)*drmag;
	deltar = DotProduct(partr.gradrho,draux,ndim)*drmag;
	limiter = max3(0.0,deltal/(deltar + small_number),deltar/(deltal + small_number));
	limiter = 2.0*limiter/(limiter*limiter + 1);
	//deltal = min3(2.0*fabs(rhol - rhor),fabs(deltal),2.0*fabs(2.0*deltal + pr - pl))*sgn(deltal);
    //deltar = min3(2.0*fabs(rhor - rhol),fabs(deltar),2.0*fabs(2.0*deltar + pl - pr))*sgn(deltar);
//limiter = 1.0;
	//if (sgn(rhol - rhor) != sgn(deltal) || sgn(deltal) == sgn(2.0*deltal + rhor - rhol)) deltal = 0.0;
	//if (sgn(rhor - rhol) != sgn(deltar) || sgn(deltar) == sgn(2.0*deltar + rhol - rhor)) deltar = 0.0;

	rhol += limiter*(FLOAT) 0.5*deltal*(1.0 - partl.sound*partl.dt/drmag);
	rhor -= limiter*(FLOAT) 0.5*deltar*(1.0 - partr.sound*partr.dt/drmag);


	// Slope limiter for pressure gradients
	deltal = DotProduct(partl.gradv[0],draux,ndim)*drmag;
	deltar = DotProduct(partr.gradv[0],draux,ndim)*drmag;
	limiter = max3(0.0,deltal/(deltar + small_number),deltar/(deltal + small_number));
	limiter = 2.0*limiter/(limiter*limiter + 1);
	//cout << "limiter : " << vl << "    " << vr << "    " << limiter << endl;
	//deltal = min3(2.0*fabs(vl - vr),fabs(deltal),2.0*fabs(2.0*deltal + pr - pl))*sgn(deltal);
    //deltar = min3(2.0*fabs(vr - vl),fabs(deltar),2.0*fabs(2.0*deltar + pl - pr))*sgn(deltar);
//limiter = 1.0;
	//if (sgn(vl - vr) != sgn(deltal) || sgn(deltal) == sgn(2.0*deltal + rhor - rhol)) deltal = 0.0;
	//if (sgn(vr - vl) != sgn(deltar) || sgn(deltar) == sgn(2.0*deltar + rhol - rhor)) deltar = 0.0;

	vl += limiter*(FLOAT) 0.5*deltal*(1.0 - partl.sound*partl.dt/drmag);
	vr -= limiter*(FLOAT) 0.5*deltar*(1.0 - partr.sound*partr.dt/drmag);

  }
  // --------------------------------------------------------------------------

  return;
}



// ============================================================================
// GodunovSph::ComputeSphDerivatives
// ..
// ============================================================================
template <typename kernelclass>
void GodunovSph<kernelclass>::ComputeSphDerivatives
(int i,                                 // id of particle
 int Nneib,                             // No. of neighbours in neibpart array
 int *neiblist,                         // id of gather neighbour in neibpart
 FLOAT *drmag,                          // Distances of gather neighbours
 FLOAT *invdrmag,                       // Inverse distances of gather neibs
 FLOAT *dr,                             // Position vector of gather neibs
 SphParticle &parti,                    // Particle i data
 SphParticle *neibpart)                 // Neighbour particle data
{
  int j;                                // Neighbour list id
  int jj;                               // Aux. neighbour counter
  int k;                                // Dimension counter
  int kk;
  FLOAT draux[ndimmax];                 // Relative position vector
  FLOAT dv[ndimmax];                    // Relative velocity vector
  FLOAT wkern;                          // Value of w1 kernel function

  for (k=0; k<ndim; k++) {
    parti.gradP[k] = (FLOAT) 0.0;
    parti.gradrho[k] = (FLOAT) 0.0;
    for (kk=0; kk<ndim; kk++) parti.gradv[k][kk] = (FLOAT) 0.0;
  }


  // Loop over all potential neighbours in the list
  // --------------------------------------------------------------------------
  for (jj=0; jj<Nneib; jj++) {
    j = neiblist[jj];
    wkern = parti.hfactor*kern.w1(drmag[jj]*parti.invh);
    
    for (k=0; k<ndim; k++) draux[k] = dr[jj*ndim + k];
    for (k=0; k<ndim; k++) dv[k] = neibpart[j].v[k] - parti.v[k];
    
    for (k=0; k<ndim; k++) parti.gradP[k] -= neibpart[j].m*
      (neibpart[j].press - parti.press)*wkern*parti.invrho*draux[k];
    for (k=0; k<ndim; k++) parti.gradrho[k] -= neibpart[j].m*
      (neibpart[j].rho - parti.rho)*wkern*parti.invrho*draux[k];
    for (k=0; k<ndim; k++) 
      for (kk=0; kk<ndim; kk++) parti.gradv[kk][k] -= neibpart[j].m*
	dv[kk]*wkern*parti.invrho*draux[k];
  }
  // --------------------------------------------------------------------------

  //cout << "Gradients : " << i << "   " << parti.v[0] << "   "
    //   << parti.gradv[0][0] << endl;

  return;
}



// ============================================================================
// GodunovSph::ComputeSphNeibDudt
// Compute SPH neighbour contributions to dudt for 
// (i) All neighbour interactions of particle i with i.d. j > i,
// (ii) Active neighbour interactions of particle j with i.d. j > i
// (iii) All inactive neighbour interactions of particle i with i.d. j < i.
// This ensures that all particle-particle pair interactions are only 
// computed once only for efficiency.
// ============================================================================
template <typename kernelclass>
void GodunovSph<kernelclass>::ComputeSphNeibDudt
(int i,                                 // id of particle
 int Nneib,                             // No. of neighbours in neibpart array
 int *neiblist,                         // id of gather neighbour in neibpart
 FLOAT *drmag,                          // Distances of gather neighbours
 FLOAT *invdrmag,                       // Inverse distances of gather neibs
 FLOAT *dr,                             // Position vector of gather neibs
 SphParticle &parti,                    // Particle i data
 SphParticle *neibpart)                 // Neighbour particle data
{
  int j;                                // Neighbour list id
  int jj;                               // Aux. neighbour counter
  int k;                                // Dimension counter
  FLOAT da[ndimmax];                    // ..
  FLOAT draux[ndimmax];                 // Relative position vector
  FLOAT dv[ndimmax];                    // Relative velocity vector
  FLOAT dvdr;                           // Dot product of dv and dr
  FLOAT wkerni;                         // Value of w1 kernel function
  FLOAT wkernj;                         // Value of w1 kernel function
  FLOAT vsignal;                        // Signal velocity
  FLOAT paux;                           // Aux. pressure force variable
  FLOAT uaux;                           // Aux. internal energy variable
  FLOAT winvrho;                        // 0.5*(wkerni + wkernj)*invrhomean

  FLOAT Cij;                            // ..
  FLOAT Dij;                            // ..
  FLOAT Vsqdi;                          // ..
  FLOAT Vsqdj;                          // ..
  FLOAT Sij;                            // ..
  FLOAT pl,pr;
  FLOAT rhol,rhor;
  FLOAT vl,vr;
  FLOAT pstar;
  FLOAT vstar;
  FLOAT vhalfi;
  FLOAT vhalfj;
  FLOAT dt;
  FLOAT vtemp[ndimmax];
  FLOAT gradi, gradj;

  static const FLOAT hconv = powf(invsqrttwo,ndim+1);


  // Compute hydro forces
  // ==========================================================================
  if (hydro_forces == 1) {

    // Loop over all potential neighbours in the list
    // ------------------------------------------------------------------------
    for (jj=0; jj<Nneib; jj++) {
      j = neiblist[jj];
      wkerni = hconv*parti.hfactor*kern.w1(invsqrttwo*drmag[jj]*parti.invh);
      wkernj = hconv*neibpart[j].hfactor*
	kern.w1(invsqrttwo*drmag[jj]*neibpart[j].invh);

      for (k=0; k<ndim; k++) draux[k] = dr[jj*ndim + k];
      for (k=0; k<ndim; k++) dv[k] = neibpart[j].v[k] - parti.v[k];
      for (k=0; k<ndim; k++) da[k] = neibpart[j].a[k] - parti.a[k];
      dvdr = DotProduct(dv,draux,ndim);
      vhalfi = (FLOAT) 0.5*DotProduct(parti.a,draux,ndim)*parti.dt;
      vhalfj = dvdr + (FLOAT) 0.5*DotProduct(neibpart[j].a,draux,ndim)*neibpart[j].dt;

      // Add contribution to velocity divergence
      parti.div_v -= neibpart[j].m*dvdr*wkerni;
      neibpart[j].div_v -= parti.m*dvdr*wkernj;

      // Linear interpolate quantites between left and right states
      Cij = (parti.invrho - neibpart[j].invrho)*invdrmag[jj];
      Dij = (FLOAT) 0.5*(parti.invrho + neibpart[j].invrho);
      Vsqdi = (FLOAT) 0.25*parti.h*parti.h*Cij*Cij + Dij*Dij;
      Vsqdj = (FLOAT) 0.25*neibpart[j].h*neibpart[j].h*Cij*Cij + Dij*Dij;
      Sij = (FLOAT) 0.25*Cij*Dij*(parti.h*parti.h/Vsqdi + 
				  neibpart[j].h*neibpart[j].h/Vsqdj);

      // Initialise the LHS and RHS of the Riemann problem
      InitialiseRiemannProblem(parti,neibpart[j],draux,drmag[jj],dvdr,
    		  parti.sound,neibpart[j].sound,pl,pr,rhol,rhor,vl,vr);

      // Now solve Riemann problem and return intermediate state variables
      if (riemann_solver == "hllc")
	HllcSolver("basic1",pl,pr,rhol,rhor,parti.sound,neibpart[j].sound,
		   vl,vr,eos->gamma,pstar,vstar);
      else if (riemann_solver == "mg")
	MgSolver("basic1",pl,pr,rhol,rhor,parti.sound,neibpart[j].sound,
		 vl,vr,eos->gamma,pstar,vstar);

      // Main SPH pressure force term
      uaux = pstar*(Vsqdi*wkerni + Vsqdj*wkernj);

      // Add total hydro contribution to acceleration for particle i
      for (k=0; k<ndim; k++) parti.dudt += neibpart[j].m*uaux*(vstar - vhalfi);
      
      // If neighbour is also active, add contribution to force here
      for (k=0; k<ndim; k++) neibpart[j].dudt -= parti.m*uaux*(vstar - vhalfj);

    }
    // ------------------------------------------------------------------------

  }
  // ==========================================================================

  return;
}



// ============================================================================
// GodunovSph::ComputePostHydroQuantities
// ..
// ============================================================================
template <typename kernelclass>
void GodunovSph<kernelclass>::ComputePostHydroQuantities
(SphParticle &parti)
{
  parti.div_v *= parti.invrho;
  //parti.dudt = (FLOAT) 0.0;
  //parti.dudt -= eos->Pressure(parti)*parti.div_v*parti.invrho*parti.invomega;

  return;
}



// ============================================================================
// GodunovSph::ComputeGravForces
// Compute the contribution to the total gravitational force of particle 'i' 
// due to 'Nneib' neighbouring particles in the list 'neiblist'.
// ============================================================================
template <typename kernelclass>
void GodunovSph<kernelclass>::ComputeDirectGravForces
(int i,                                 // id of particle
 int Ndirect,                           // No. of nearby 'gather' neighbours
 int *directlist,                       // id of gather neighbour in neibpart
 SphParticle &parti,                    // Particle i data
 SphParticle *sph)                      // Neighbour particle data
{
  return;
}



// ============================================================================
// HllcSolver
// Toro et al. (???) improved HLL solver.
// ============================================================================
template <typename kernelclass>
void GodunovSph<kernelclass>::HllcSolver
(string wave_speed,
 FLOAT pl,
 FLOAT pr,
 FLOAT rhol,
 FLOAT rhor,
 FLOAT soundl,
 FLOAT soundr,
 FLOAT vl,
 FLOAT vr,
 FLOAT gamma_eos,
 FLOAT &pstar,
 FLOAT &Sstar)
{
  FLOAT ql = 1.0;
  FLOAT qr = 1.0;
  FLOAT Sl;
  FLOAT Sr;
  FLOAT gamma_fac = 0.5f*(gamma_eos - 1.0f)/gamma_eos;

  // Compute wave speeds based on various options
  if (wave_speed == "basic1") { 
    Sl = vl - soundl;
    Sr = vr + soundr;
  }
  else if (wave_speed == "basic2") { 
    Sl = vl - soundl;
    Sr = vr + soundr;
  }
  else if (wave_speed == "pressure") {
    FLOAT rhomean = 0.5f*(rhol + rhor);
    FLOAT soundmean = 0.5f*(soundl + soundr);
    pstar = 0.5f*(pl + pr) - 0.5f*(vr - vl)*rhomean*soundmean;
    //Sstar = 0.5f*(vl + vr) - 0.5f*(pr - rl)/rhomean/soundmean;

    // Compute pressure-based speeds depending on pressure values
    if (pstar > pl) ql = sqrt(1.0f + gamma_fac*(pstar/pl - 1.0f));
    if (pstar > pr) qr = sqrt(1.0f + gamma_fac*(pstar/pr - 1.0f));
    Sl = vl - ql*soundl;
    Sr = vr + qr*soundr;
  }
  else {
    cout << "Unrecognised wave_speed option" << endl;
    Sstar = 0.0;
    pstar = 0.0;
    return;
  }

  // Compute intermediate ('star') velocity and pressure
  Sstar = (pr - pl + rhol*vl*(Sl - vl) - rhor*vr*(Sr - vr))/
    (rhol*(Sl - vl) - rhor*(Sr - vr));
  pstar = ((Sr - vr)*rhor*pl - (Sl - vl)*rhol*pr + rhol*rhor*(Sr - vr)*
	   (Sl - vl)*(vr - vl))/((Sr - vr)*rhor - (Sl - vl)*rhol);

  //pstar = pl + rhol*(Sl - vl)*(Sstar - vl);

  return;
}



// ============================================================================
// MgSolver
// Riemann solver implemented in MG (courtesy of S. Falle).
// ============================================================================
template <typename kernelclass>
void GodunovSph<kernelclass>::MgSolver
(string wave_speed,
 FLOAT pl,
 FLOAT pr,
 FLOAT rhol,
 FLOAT rhor,
 FLOAT soundl,
 FLOAT soundr,
 FLOAT vl,
 FLOAT vr,
 FLOAT gamma_eos,
 FLOAT &pstar,
 FLOAT &Sstar)
{
  FLOAT soundi,pressi,Wl,Wr;
  FLOAT tol = (FLOAT) 0.1;
  FLOAT gamma_aux1 = (FLOAT) 0.5*(gamma_eos - (FLOAT) 1.0)/gamma_eos;
  FLOAT gamma_aux2 = (FLOAT) 0.5*(gamma_eos + (FLOAT) 1.0)/gamma_eos;

  // Linear Riemann solver
  pstar = (soundr*rhor*pl + soundl*rhol*pr - rhor*soundr*rhol*soundl*
	   (vr - vl))/(rhor*soundr + rhol*soundl);
  Sstar = (soundr*rhor*vr + soundl*rhol*vl - (pr - pl))/
    (rhor*soundr + rhol*soundl);

  // If linear solver is not good enough an approximation, use iterative 
  // Riemann solver.
  // --------------------------------------------------------------------------
  if (fabs(pstar - pl) > tol*pl || fabs(pstar - pr) > tol*pr) {


    // Non-iterative solution for two rarefactions
    // ------------------------------------------------------------------------
    if (pstar < pl && pstar < pr) {
      soundi = soundl / (rhol*powf(pl,gamma_aux1));
      pstar = (0.5*(gamma_eos - 1.0)*(vr - vl) + soundl/rhol + soundr/rhor)/
	(soundi + soundr/(rhor*powf(pr,gamma_aux1)));
      pstar = fmax(small_number,powf(pstar,1.0/gamma_aux1));
      Wl = -soundl*rhol*wave(pstar,pl,gamma_aux1,gamma_aux2);
      Wr = soundr*rhor*wave(pstar,pr,gamma_aux1,gamma_aux2);
    }

    // Iterative solution for shock/rarefaction
    // ------------------------------------------------------------------------
    else {
      pressi = pstar;
      Wl = -soundl*rhol*wave(pstar,pl,gamma_aux1,gamma_aux2);
      Wr = soundr*rhor*wave(pstar,pr,gamma_aux1,gamma_aux2);
      int nit = 0;
      bool doit = true;

      // Iteration loop
      // ----------------------------------------------------------------------
      while (doit) {
	pstar = fmax((Wr*pl - Wl*pr + Wr*Wl*(vr - vl))/(Wr - Wl),small_number);

	if (fabs(pstar - pressi) > tol*pstar) {
	  pressi = pstar;
	  Wl = -soundl*rhol*wave(pstar,pl,gamma_aux1,gamma_aux2);
	  Wr = soundr*rhor*wave(pstar,pr,gamma_aux1,gamma_aux2);
	  nit++;

	  if (nit > 50) {
	    cout << "Convergence failure in Riemann solver" << endl;
	    pstar = (soundr*rhor*pl + soundl*rhol*pr - rhor*soundr*rhol*soundl*
		     (vr - vl))/(rhor*soundr + rhol*soundl);
	    Wl = -soundl;
	    Wr = soundr;
	    doit = false;
	  }
	}
	else 
	  doit = false;

      }
      // ----------------------------------------------------------------------

    }
    // ------------------------------------------------------------------------

    // Calculate resolved velocity
    Sstar = (Wr*vr - Wl*vl - pr + pl)/(Wr - Wl);

  }
  // --------------------------------------------------------------------------


  return;
}





template class GodunovSph<M4Kernel>;
template class GodunovSph<QuinticKernel>;
template class GodunovSph<GaussianKernel>;
template class GodunovSph<TabulatedKernel>;
