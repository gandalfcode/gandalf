//=================================================================================================
//  LV2008MFV.cpp
//  Contains all functions for calculating Meshless Finite-Volume Hydrodynamics quantities.
//
//  This file is part of GANDALF :
//  Graphical Astrophysics code for N-body Dynamics And Lagrangian Fluids
//  https://github.com/gandalfcode/gandalf
//  Contact : gandalfcode@gmail.com
//
//  Copyright (C) 2013  D. A. Hubber, G. Rosotti
//
//  GANDALF is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 2 of the License, or
//  (at your option) any later version.
//
//  GANDALF is distributed in the hope that it will be useful, but
//  WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  General Public License (http://www.gnu.org/licenses) for more details.
//=================================================================================================


#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cassert>
#include <iostream>
#include <math.h>
#include "Precision.h"
#include "MeshlessFV.h"
#include "Particle.h"
#include "Parameters.h"
#include "SmoothingKernel.h"
#include "EOS.h"
#include "Debug.h"
#include "Exception.h"
#include "InlineFuncs.h"
using namespace std;



//=================================================================================================
//  LV2008MFV::LV2008MFV
/// LV2008MFV class constructor.  Calls main SPH class constructor and also
/// sets additional kernel-related quantities
//=================================================================================================
template <int ndim, template<int> class kernelclass>
LV2008MFV<ndim, kernelclass>::LV2008MFV(int hydro_forces_aux, int self_gravity_aux, FLOAT h_fac_aux,
  FLOAT h_converge_aux, FLOAT gamma_aux, string gas_eos_aux, string KernelName, int size_sph):
  MeshlessFV<ndim>(hydro_forces_aux, self_gravity_aux, h_fac_aux, h_converge_aux,
                   gamma_aux, gas_eos_aux, KernelName, size_sph),
  kern(kernelclass<ndim>(KernelName))
{
  this->kernp      = &kern;
  this->kernfac    = (FLOAT) 1.0;
  this->kernfacsqd = (FLOAT) 1.0;
  this->kernrange  = this->kernp->kernrange;
}



//=================================================================================================
//  LV2008MFV::~LV2008MFV
/// LV2008MFV class destructor
//=================================================================================================
template <int ndim, template<int> class kernelclass>
LV2008MFV<ndim, kernelclass>::~LV2008MFV()
{
  //DeallocateMemory();
}



//=================================================================================================
//  LV2008MFV::ComputeH
/// Compute the value of the smoothing length of particle 'i' by iterating the relation :
/// h = h_fac*(m/rho)^(1/ndim).
/// Uses the previous value of h as a starting guess and then uses either a Newton-Rhapson solver,
/// or fixed-point iteration, to converge on the correct value of h.  The maximum tolerance used
/// for deciding whether the iteration has converged is given by the 'h_converge' parameter.
//=================================================================================================
template <int ndim, template<int> class kernelclass>
int LV2008MFV<ndim, kernelclass>::ComputeH
 (const int i,                         ///< [in] id of particle
  const int Nneib,                     ///< [in] No. of potential neighbours
  const FLOAT hmax,                    ///< [in] Max. h permitted by neib list
  FLOAT *m,                            ///< [in] Array of neib. masses
  FLOAT *mu,                           ///< [in] Array of m*u (not needed here)
  FLOAT *drsqd,                        ///< [in] Array of neib. distances squared
  FLOAT *gpot,                         ///< [in] Array of neib. grav. potentials
  MeshlessFVParticle<ndim> &part,      ///< [inout] Particle i data
  Nbody<ndim> *nbody)                  ///< [in] Main N-body object
{
  int j;                               // Neighbour id
  int k;                               // Dimension counter
  int iteration = 0;                   // h-rho iteration counter
  int iteration_max = 30;              // Max. no of iterations
  FLOAT dr[ndim];                      // Relative position vector
  FLOAT h_lower_bound = (FLOAT) 0.0;   // Lower bound on h
  FLOAT h_upper_bound = hmax;          // Upper bound on h
  FLOAT invhsqd;                       // (1 / h)^2
  FLOAT ssqd;                          // Kernel parameter squared, (r/h)^2


  // Main smoothing length iteration loop
  //===============================================================================================
  do {

    // Initialise all variables for this value of h
    iteration++;
    part.invh    = (FLOAT) 1.0/part.h;
    part.ndens   = (FLOAT) 0.0;
    part.hfactor = pow(part.invh,ndim);
    invhsqd      = part.invh*part.invh;

    // Loop over all nearest neighbours in list to calculate
    // density, omega and zeta.
    //---------------------------------------------------------------------------------------------
    for (j=0; j<Nneib; j++) {
      ssqd        = drsqd[j]*invhsqd;
      part.ndens += kern.w0_s2(ssqd);
    }
    //---------------------------------------------------------------------------------------------

    part.ndens *= part.hfactor;
    part.volume = (FLOAT) 1.0/part.ndens;
    part.rho    = part.m*part.ndens;


    if (part.rho > (FLOAT) 0.0) part.invrho = (FLOAT) 1.0/part.rho;

    // If h changes below some fixed tolerance, exit iteration loop
    if (part.rho > (FLOAT) 0.0 && part.h > h_lower_bound &&
        fabs(part.h - h_fac*pow(part.volume,MeshlessFV<ndim>::invndim)) < h_converge) break;

    // Use fixed-point iteration, i.e. h_new = h_fac*(m/rho_old)^(1/ndim), for now.  If this does
    // not converge in a reasonable number of iterations (iteration_max), then assume something is
    // wrong and switch to a bisection method, which should be guaranteed to converge, albeit much
    // more slowly.  (N.B. will implement Newton-Raphson soon)
    //---------------------------------------------------------------------------------------------
    if (iteration < iteration_max) {
      part.h = h_fac*pow(part.volume,MeshlessFV<ndim>::invndim);
    }
    else if (iteration == iteration_max) {
      part.h = (FLOAT) 0.5*(h_lower_bound + h_upper_bound);
    }
    else if (iteration < 5*iteration_max) {
      if (part.ndens < small_number || part.ndens*pow(part.h,ndim) > pow(h_fac,ndim))
        h_upper_bound = part.h;
      else
        h_lower_bound = part.h;
      part.h = (FLOAT) 0.5*(h_lower_bound + h_upper_bound);
    }
    else {
      cout << "H ITERATION : " << iteration << "    h : " << part.h
           << "   rho : " << part.rho << "   h_upper " << h_upper_bound << "    hmax :  " << hmax
           << "   h_lower : " << h_lower_bound << "    " << part.hfactor << "    m : " << part.m
           << "     " << part.m*part.hfactor*kern.w0(0.0) << "    " << Nneib << endl;
      string message = "Problem with convergence of h-rho iteration";
      ExceptionHandler::getIstance().raise(message);
    }

    // If the smoothing length is too large for the neighbour list, exit routine and flag neighbour
    // list error in order to generate a larger neighbour list (not properly implemented yet).
    if (part.h > hmax) return 0;

  } while (part.h > h_lower_bound && part.h < h_upper_bound);
  //===============================================================================================


  // Normalise all SPH sums correctly
  part.h         = max(h_fac*pow(part.volume,MeshlessFV<ndim>::invndim),h_lower_bound);
  part.invh      = (FLOAT) 1.0/part.h;
  part.hfactor   = pow(part.invh,ndim+1);
  part.hrangesqd = kernfacsqd*kern.kernrangesqd*part.h*part.h;
  part.div_v     = (FLOAT) 0.0;


  // Set important thermal variables here
  this->ComputeThermalProperties(part);


  // If h is invalid (i.e. larger than maximum h), then return error code (0)
  if (part.h <= hmax) return 1;
  else return -1;
}



//=================================================================================================
//  LV2008MFV::ComputeDerivatives
/// Compute SPH neighbour force pairs for
/// (i) All neighbour interactions of particle i with i.d. j > i,
/// (ii) Active neighbour interactions of particle j with i.d. j > i
/// (iii) All inactive neighbour interactions of particle i with i.d. j < i.
/// This ensures that all particle-particle pair interactions are only
/// computed once only for efficiency.
//=================================================================================================
template <int ndim, template<int> class kernelclass>
void LV2008MFV<ndim, kernelclass>::ComputePsiFactors
 (const int i,                                 ///< [in] id of particle
  const int Nneib,                             ///< [in] No. of neins in neibpart array
  int *neiblist,                               ///< [in] id of gather neibs in neibpart
  FLOAT *drmag,                                ///< [in] Distances of gather neighbours
  FLOAT *invdrmag,                             ///< [in] Inverse distances of gather neibs
  FLOAT *dr,                                   ///< [in] Position vector of gather neibs
  MeshlessFVParticle<ndim> &part,              ///< [inout] Particle i data
  MeshlessFVParticle<ndim> *neibpart)          ///< [inout] Neighbour particle data
{
  int j;                                       // Neighbour list id
  int jj;                                      // ..
  int k;                                       // Dimension counter
  FLOAT invdet;                                // Determinant
  FLOAT draux[ndim];                           // Relative position vector
  FLOAT drsqd;                                 // ..
  FLOAT E[ndim][ndim];                         // ..
  const FLOAT invhsqd = part.invh*part.invh;   // ..


  for (k=0; k<ndim; k++) {
    for (int kk=0; kk<ndim; kk++) {
      E[k][kk] = (FLOAT) 0.0;
      part.B[k][kk] = 0.0;
    }
  }


  // Loop over all potential neighbours in the list
  //-----------------------------------------------------------------------------------------------
  for (jj=0; jj<Nneib; jj++) {
    j = neiblist[jj];

    for (k=0; k<ndim; k++) draux[k] = neibpart[j].r[k] - part.r[k];
    drsqd = DotProduct(draux, draux, ndim);

    for (k=0; k<ndim; k++) {
      for (int kk=0; kk<ndim; kk++) {
        E[k][kk] += draux[k]*draux[kk]*part.hfactor*kern.w0_s2(drsqd*invhsqd)/part.ndens;
      }
    }
  }
  //-----------------------------------------------------------------------------------------------

  // Invert the matrix (depending on dimensionality)
  if (ndim == 1) {
    part.B[0][0] = 1.0/E[0][0];
  }
  else if (ndim == 2) {
    invdet = 1.0/(E[0][0]*E[1][1] - E[0][1]*E[1][0]);
    part.B[0][0] = invdet*E[1][1];
    part.B[0][1] = -invdet*E[0][1];
    part.B[1][0] = -invdet*E[1][0];
    part.B[1][1] = invdet*E[0][0];
  }


  return;
}



//=================================================================================================
//  LV2008MFV::ComputeDerivatives
/// Compute SPH neighbour force pairs for
/// (i) All neighbour interactions of particle i with i.d. j > i,
/// (ii) Active neighbour interactions of particle j with i.d. j > i
/// (iii) All inactive neighbour interactions of particle i with i.d. j < i.
/// This ensures that all particle-particle pair interactions are only
/// computed once only for efficiency.
//=================================================================================================
template <int ndim, template<int> class kernelclass>
void LV2008MFV<ndim, kernelclass>::ComputeGradients
 (const int i,                         ///< [in] id of particle
  const int Nneib,                     ///< [in] No. of neins in neibpart array
  int *neiblist,                       ///< [in] id of gather neibs in neibpart
  FLOAT *drmag,                        ///< [in] Distances of gather neighbours
  FLOAT *invdrmag,                     ///< [in] Inverse distances of gather neibs
  FLOAT *dr,                           ///< [in] Position vector of gather neibs
  MeshlessFVParticle<ndim> &part,      ///< [inout] Particle i data
  MeshlessFVParticle<ndim> *neibpart)  ///< [inout] Neighbour particle data
{
  int j;                               // Neighbour list id
  int jj;                              // Aux. neighbour counter
  int k;                               // Dimension counter
  int var;
  FLOAT alpha_mean;                    // Mean articial viscosity alpha value
  FLOAT draux[ndim];                   // Relative position vector
  FLOAT drsqd;
  FLOAT dvdr;                          // Dot product of dv and dr
  FLOAT wkerni;                        // Value of w1 kernel function
  FLOAT wkernj;                        // Value of w1 kernel function
  FLOAT vsignal;                       // Signal velocity
  FLOAT paux;                          // Aux. pressure force variable
  FLOAT winvrho;                       // 0.5*(wkerni + wkernj)*invrhomean
  FLOAT psitilda[ndim];
  const FLOAT invhsqd = part.invh*part.invh;


  // Zero gradients
  for (k=0; k<ndim; k++) {
    for (var=0; var<ndim+2; var++) {
      part.grad[k][var] = 0.0;
    }
  }
  for (var=0; var<nvar; var++) {
    part.Wmin[var] = part.Wprim[var];
    part.Wmax[var] = part.Wprim[var];
    part.Wmidmin[var] = big_number;
    part.Wmidmax[var] = -big_number;
  }


  // Loop over all potential neighbours in the list
  //-----------------------------------------------------------------------------------------------
  for (jj=0; jj<Nneib; jj++) {
    j = neiblist[jj];

    for (k=0; k<ndim; k++) draux[k] = neibpart[j].r[k] - part.r[k];
    drsqd = DotProduct(draux, draux, ndim);

    // Calculate psitilda values
    for (k=0; k<ndim; k++) psitilda[k] = 0.0;
    for (k=0; k<ndim; k++) {
      for (int kk=0; kk<ndim; kk++) {
        psitilda[k] += part.B[k][kk]*draux[kk]*part.hfactor*kern.w0_s2(drsqd*invhsqd)/part.ndens;
      }
    }

    // Calculate contribution to gradient from neighbour
    for (var=0; var<ndim+2; var++) {
      for (k=0; k<ndim; k++) {
        part.grad[k][var] += (neibpart[j].Wprim[var] - part.Wprim[var])*psitilda[k];
      }
    }


  }
  //-----------------------------------------------------------------------------------------------


  // Find all max and min values for meshless slope limiters
  //-----------------------------------------------------------------------------------------------
  for (jj=0; jj<Nneib; jj++) {
    j = neiblist[jj];

    for (k=0; k<ndim; k++) draux[k] = neibpart[j].r[k] - part.r[k];
    drsqd = DotProduct(draux, draux, ndim);

    // Calculate min and max values of primitives for slope limiters
    for (var=0; var<nvar; var++) {
      part.Wmin[var] = min(part.Wmin[var], neibpart[j].Wprim[var]);
      part.Wmax[var] = max(part.Wmax[var], neibpart[j].Wprim[var]);
      part.Wmidmin[var] = min(part.Wmidmin[var], part.Wprim[var] + 0.5*part.grad[var][0]*draux[0]);
      part.Wmidmax[var] = max(part.Wmidmax[var], part.Wprim[var] + 0.5*part.grad[var][0]*draux[0]);
      assert(part.Wmidmax[var] >= part.Wmidmin[var]);
      assert(part.Wmax[var] >= part.Wmin[var]);
    }

  }
  //-----------------------------------------------------------------------------------------------



  return;
}



//=================================================================================================
//  LV2008MFV::ComputeGodunovFlux
/// ...
//=================================================================================================
template <int ndim, template<int> class kernelclass>
void LV2008MFV<ndim, kernelclass>::ComputeGodunovFlux
 (const int i,                         ///< [in] id of particle
  const int Nneib,                     ///< [in] No. of neins in neibpart array
  const FLOAT timestep,                ///< ..
  int *neiblist,                       ///< [in] id of gather neibs in neibpart
  FLOAT *drmag,                        ///< [in] Distances of gather neighbours
  FLOAT *invdrmag,                     ///< [in] Inverse distances of gather neibs
  FLOAT *dr,                           ///< [in] Position vector of gather neibs
  MeshlessFVParticle<ndim> &part,      ///< [inout] Particle i data
  MeshlessFVParticle<ndim> *neibpart)  ///< [inout] Neighbour particle data
{
  int j;                               // Neighbour list id
  int jj;                              // Aux. neighbour counter
  int k;                               // Dimension counter
  int var;                             // ..
  FLOAT Aij[ndim];                     // ..
  FLOAT alpha;                         // ..
  FLOAT draux[ndim];                   // Position vector of part relative to neighbour
  FLOAT dr_unit[ndim];                 // Unit vector from neighbour to part
  FLOAT drsqd;                         // ..
  FLOAT dvdr;                          // Dot product of dv and dr
  FLOAT invdrmagaux;
  FLOAT psitildai[ndim];
  FLOAT psitildaj[ndim];
  FLOAT rface[ndim];
  FLOAT vface[ndim];
  FLOAT flux[nvar][ndim];
  FLOAT Wleft[nvar];
  FLOAT Wright[nvar];
  FLOAT Wdot[nvar];
  FLOAT gradW[nvar][ndim];
  FLOAT dW[nvar];


  // Initialise all particle flux variables
  for (var=0; var<nvar; var++) part.dQdt[var] = 0.0;


  // Loop over all potential neighbours in the list
  //-----------------------------------------------------------------------------------------------
  for (jj=0; jj<Nneib; jj++) {
    j = neiblist[jj];

    // Calculate position and velocity of the face
    for (k=0; k<ndim; k++) rface[k] = (FLOAT) 0.5*(part.r[k] + neibpart[j].r[k]);
    for (k=0; k<ndim; k++) vface[k] = (FLOAT) 0.5*(part.v[k] + neibpart[j].v[k]);

    // Compute slope-limited values for LHS
    for (k=0; k<ndim; k++) draux[k] = rface[k] - part.r[k];
    limiter.ComputeLimitedSlopes(part, neibpart[j], draux, gradW, dW);
    for (var=0; var<nvar; var++) Wleft[var] = part.Wprim[var] + dW[var];
    for (k=0; k<ndim; k++) Wleft[k] -= vface[k];

    // Time-integrate LHS state to half-timestep value
    this->CalculatePrimitiveTimeDerivative(Wleft, gradW, Wdot);
    for (var=0; var<nvar; var++) Wleft[var] -= (FLOAT) 0.5*timestep*Wdot[var];


    // Compute slope-limited values for RHS
    for (k=0; k<ndim; k++) draux[k] = rface[k] - neibpart[j].r[k];
    limiter.ComputeLimitedSlopes(neibpart[j], part, draux, gradW, dW);
    for (var=0; var<nvar; var++) Wright[var] = neibpart[j].Wprim[var] + dW[var];
    for (k=0; k<ndim; k++) Wright[k] -= vface[k];

    // Time-integrate RHS state to half-timestep value
    this->CalculatePrimitiveTimeDerivative(Wright, gradW, Wdot);
    for (var=0; var<nvar; var++) Wright[var] -= (FLOAT) 0.5*timestep*Wdot[var];

    for (k=0; k<ndim; k++) draux[k] = part.r[k] - neibpart[j].r[k];
    drsqd = DotProduct(draux, draux, ndim);
    invdrmagaux = 1.0/sqrt(drsqd + small_number);
    for (k=0; k<ndim; k++) dr_unit[k] = draux[k]*invdrmagaux;


    // Calculate Godunov flux using the selected Riemann solver
    riemann->ComputeFluxes(Wright, Wleft, dr_unit, vface, flux);
    /*if (part.r[0] > neibpart[j].r[0]) {
      riemann->ComputeFluxes(Wright, Wleft, dr_unit, vface, flux);
    }
    else {
      riemann->ComputeFluxes(Wleft, Wright, dr_unit, vface, flux);
    }*/


    // Calculate psitilda values
    for (k=0; k<ndim; k++) {
      psitildai[k] = (FLOAT) 0.0;
      psitildaj[k] = (FLOAT) 0.0;
      for (int kk=0; kk<ndim; kk++) {
        psitildai[k] += neibpart[j].B[k][kk]*draux[kk]*neibpart[j].hfactor*
          kern.w0_s2(drsqd*neibpart[j].invh*neibpart[j].invh)/neibpart[j].ndens;
        psitildaj[k] -= part.B[k][kk]*draux[kk]*part.hfactor*
          kern.w0_s2(drsqd*part.invh*part.invh)/part.ndens;
      }
      Aij[k] = part.volume*psitildaj[k] - neibpart[j].volume*psitildai[k];
    }


    // Now boost back into main coordinate frame
    /*for (k=0; k<ndim; k++) {
      flux[ietot][k] += (FLOAT) 0.5*DotProduct(vface, vface, ndim)*flux[irho][k] + vface[k]*flux[ivx][k];
      flux[ivx][k] += vface[k]*flux[irho][k];
    }*/
    

    // Finally calculate flux terms for all quantities based on Lanson & Vila gradient operators
    for (var=0; var<nvar; var++) {
      part.dQdt[var] -= DotProduct(flux[var], Aij, ndim);
    }

  }
  //-----------------------------------------------------------------------------------------------


  return;
}



//=================================================================================================
//  Ghosts::CopySphDataToGhosts
/// Copy any newly calculated data from original SPH particles to ghosts.
//=================================================================================================
template <int ndim, template<int> class kernelclass>
void LV2008MFV<ndim, kernelclass>::CopyDataToGhosts
 (DomainBox<ndim> &simbox,
  MeshlessFVParticle<ndim> *partdata)  ///< [inout] Neighbour particle data
{
  int i;                               // Particle id
  int iorig;                           // Original (real) particle id
  int itype;                           // Ghost particle type
  int j;                               // Ghost particle counter

  debug2("[LV2008MFV::CopySphDataToGhosts]");


  //---------------------------------------------------------------------------
//#pragma omp parallel for default(none) private(i,iorig,itype,j) shared(simbox,sph,partdata)
  for (j=0; j<this->NPeriodicGhost; j++) {
    i = this->Nhydro + j;
    iorig = partdata[i].iorig;
    itype = partdata[i].itype;

    partdata[i]        = partdata[iorig];
    partdata[i].iorig  = iorig;
    partdata[i].itype  = itype;
    partdata[i].active = false;

    // Modify ghost position based on ghost type
    if (itype == x_lhs_periodic) {
      partdata[i].r[0] += simbox.boxsize[0];
    }
    else if (itype == x_rhs_periodic) {
      partdata[i].r[0] -= simbox.boxsize[0];
    }
    else if (itype == y_lhs_periodic && ndim > 1) {
      partdata[i].r[1] += simbox.boxsize[1];
    }
    else if (itype == y_rhs_periodic && ndim > 1) {
      partdata[i].r[1] -= simbox.boxsize[1];
    }
    else if (itype == z_lhs_periodic && ndim == 3) {
      partdata[i].r[2] += simbox.boxsize[2];
    }
    else if (itype == z_rhs_periodic && ndim == 3) {
      partdata[i].r[2] -= simbox.boxsize[2];
    }

  }
  //---------------------------------------------------------------------------

  return;
}





template class LV2008MFV<1, M4Kernel>;
template class LV2008MFV<2, M4Kernel>;
template class LV2008MFV<3, M4Kernel>;
template class LV2008MFV<1, QuinticKernel>;
template class LV2008MFV<2, QuinticKernel>;
template class LV2008MFV<3, QuinticKernel>;
template class LV2008MFV<1, GaussianKernel>;
template class LV2008MFV<2, GaussianKernel>;
template class LV2008MFV<3, GaussianKernel>;
template class LV2008MFV<1, TabulatedKernel>;
template class LV2008MFV<2, TabulatedKernel>;
template class LV2008MFV<3, TabulatedKernel>;
