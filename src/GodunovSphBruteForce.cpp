//=============================================================================
//  GodunovSphBruteForce.cpp
//  Contains all routines for generating SPH neighbour lists using 
//  brute-force (i.e. direct summation over all particles).
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
//=============================================================================



#include <iostream>
#include <math.h>
#include "SphNeighbourSearch.h"
#include "Sph.h"
#include "Parameters.h"
#include "SphParticle.h"
#include "Debug.h"
#include "InlineFuncs.h"
#include "SphKernel.h"
#if defined MPI_PARALLEL
#include "MpiNode.h"
#endif
using namespace std;



//=============================================================================
//  GodunovSphBruteForce::GodunovSphBruteForce
/// GodunovSphBruteForce class constructor
//=============================================================================
template <int ndim, template<int> class ParticleType>
GodunovSphBruteForce<ndim,ParticleType>::GodunovSphBruteForce
(FLOAT kernrangeaux,
 DomainBox<ndim> *boxaux,
 SphKernel<ndim> *kernaux,
 CodeTiming *timingaux):
  BruteForceSearch<ndim,ParticleType>(kernrangeaux,boxaux,kernaux,timingaux)
{
}



//=============================================================================
//  GodunovSphBruteForce::~GodunovSphBruteForce
/// GodunovSphBruteForce class destructor
//=============================================================================
template <int ndim, template<int> class ParticleType>
GodunovSphBruteForce<ndim,ParticleType>::~GodunovSphBruteForce()
{
}



//=============================================================================
//  GodunovSphBruteForce::UpdateAllSphProperties
/// Routine for computing SPH properties (smoothing lengths, densities and 
/// forces) for all active SPH particle using neighbour lists generated 
/// using brute force (i.e. direct summation).
//=============================================================================
template <int ndim, template<int> class ParticleType>
void GodunovSphBruteForce<ndim,ParticleType>::UpdateAllSphProperties
(int Nsph,                          ///< [in] No. of SPH particles
 int Ntot,                          ///< [in] No. of SPH + ghost particles
 SphParticle<ndim> *sph_gen,        ///< [inout] Pointer to SPH ptcl array
 Sph<ndim> *sph,                    ///< [in] Pointer to SPH object
 Nbody<ndim> *nbody)                ///< [in] Pointer to N-body object
{
  int i,j,jj,k;                     // Particle and dimension counters
  int Nneib = 0;                    // No. of (non-dead) neighbours
  int okflag;                       // Flag valid smoothing length
  int *neiblist;                    // List of neighbours
  FLOAT dr[ndim];                   // Relative distance vector
  FLOAT rp[ndim];                   // Position of current particle
  FLOAT *drsqd;                     // Distance squared
  FLOAT *gpot;                      // Array of neib. grav. potentials
  FLOAT *m;                         // Array of neib. position vectors
  FLOAT *mu;                        // Array of neib. mass*u values
  ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sph_gen);

  debug2("[GodunovSphBruteForce::UpdateAllSphProperties]");

  // Store masses in separate array
  gpot = new FLOAT[Ntot];
  m = new FLOAT[Ntot];
  mu = new FLOAT[Ntot];
  neiblist = new int[Ntot];
  for (i=0; i<Ntot; i++) {
    if (sphdata[i].itype == dead) continue;
    neiblist[Nneib] = i;
    gpot[Nneib] = sphdata[i].gpot;
    m[Nneib] = sphdata[i].m;
    mu[Nneib] = sphdata[i].m*sphdata[i].u;
    Nneib++;
  }

  // Create parallel threads
  //===========================================================================
#pragma omp parallel default(none) private(dr,drsqd,i,j,jj,k,okflag,rp)	\
  shared(gpot,m,mu,nbody,neiblist,Nneib,Nsph,Ntot,sph,sphdata)
  {
    drsqd = new FLOAT[Ntot];

    // Compute smoothing lengths of all SPH particles
    //-------------------------------------------------------------------------
#pragma omp for
    for (i=0; i<Nsph; i++) {

      // Skip over inactive particles
      if (!sphdata[i].active || sphdata[i].itype == dead) continue;

      for (k=0; k<ndim; k++) rp[k] = sphdata[i].r[k];

      // Compute distances and the reciprical between the current particle 
      // and all neighbours here
      //-----------------------------------------------------------------------
      for (jj=0; jj<Nneib; jj++) { 
	j = neiblist[jj];
    	for (k=0; k<ndim; k++) dr[k] = sphdata[j].r[k] - rp[k];
    	drsqd[jj] = DotProduct(dr,dr,ndim);
      }
      //-----------------------------------------------------------------------

      // Compute all SPH gather properties
      okflag = sph->ComputeH(i,Nneib,big_number,m,mu,drsqd,
                             gpot,sphdata[i],nbody);
  
    }
    //-------------------------------------------------------------------------

    delete[] drsqd;

  }
  //===========================================================================

  delete[] neiblist;
  delete[] mu;
  delete[] m;
  delete[] gpot;

  return;
}




//=============================================================================
//  GodunovSphBruteForce::UpdateAllSphDerivatives
/// Compute all SPH derivatives required for 2nd-order Riemann solver in 
/// Godunov SPH method.
//=============================================================================
template <int ndim, template<int> class ParticleType>
void GodunovSphBruteForce<ndim,ParticleType>::UpdateAllSphDerivatives
(int Nsph,                            ///< [in] ..
 int Ntot,                            ///< [in] ..
 SphParticle<ndim> *sph_gen, 
 Sph<ndim> *sph)                      ///< [inout] Pointer to SPH object
{
  int i,j,k;                          // Particle and dimension counters
  int Nneib;                          // No. of neighbours
  int *neiblist;                      // List of neighbour ids
  FLOAT draux[ndim];                  // Relative distance vector
  FLOAT drsqd;                        // Distance squared
  FLOAT hrangesqd;                    // Kernel extent (squared)
  FLOAT rp[ndim];                     // Position of current particle
  FLOAT *dr;                          // Array of neib. position vectors
  FLOAT *drmag;                       // Array of neib. distances
  FLOAT *invdrmag;                    // Array of neib. inverse distances
  SphParticle<ndim> *neibpart;        // Local copies of neib. particles
  ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sph_gen);

  debug2("[GodunovSphBruteForce::UpdateAllSphForces]");

  // The potential number of neighbours is given by ALL the particles
  Nneib = Ntot;

  // Allocate memory for storing neighbour ids and position data
  neiblist = new int[Ntot];
  dr = new FLOAT[ndim*Ntot];
  drmag = new FLOAT[Ntot];
  invdrmag = new FLOAT[Ntot];
  neibpart = new SphParticle<ndim>[Ntot];

  // Record local copies of (all) neighbour properties
  for (j=0; j<Ntot; j++) neibpart[j] = sphdata[j];


  // Compute smoothing lengths of all SPH particles
  //---------------------------------------------------------------------------
  for (i=0; i<Nsph; i++) {
    for (k=0; k<ndim; k++) rp[k] = sphdata[i].r[k];
    hrangesqd = pow(kernp->kernrange*sphdata[i].h,2);
    Nneib = 0;

    // Compute distances and the reciprical between the current particle 
    // and all neighbours here
    //-------------------------------------------------------------------------
    for (j=0; j<Ntot; j++) {
      for (k=0; k<ndim; k++) draux[k] = sphdata[j].r[k] - rp[k];
      drsqd = DotProduct(draux,draux,ndim);
      if (drsqd < hrangesqd) {
    	neiblist[Nneib] = j;
    	drmag[Nneib] = sqrt(drsqd);
    	invdrmag[Nneib] = (FLOAT) 1.0/(drmag[Nneib] + small_number);
    	for (k=0; k<ndim; k++) dr[Nneib*ndim + k] = draux[k]*invdrmag[Nneib];
    	Nneib++;
      }
    }
    //-------------------------------------------------------------------------

    // Compute all SPH hydro forces
    sph->ComputeSphDerivatives(i,Nneib,neiblist,drmag,invdrmag,dr,
			       sphdata[i],neibpart);

  }
  //---------------------------------------------------------------------------

  delete[] neibpart;
  delete[] invdrmag;
  delete[] drmag;
  delete[] dr;
  delete[] neiblist;

  return;
}



//=============================================================================
//  GodunovSphBruteForce::UpdateAllSphDudt
/// Compute the compressional heating rate (dudt) for all active particles.
//=============================================================================
template <int ndim, template<int> class ParticleType>
void GodunovSphBruteForce<ndim,ParticleType>::UpdateAllSphDudt
(int Nsph,                            ///< [in] ..
 int Ntot,                            ///< [in] ..
 SphParticle<ndim> *sph_gen,          ///< ..
 Sph<ndim> *sph)                      ///< [inout] Pointer to SPH object
{
  int i,j,k;                          // Particle and dimension counters
  int Nneib;                          // No. of neighbours
  int *neiblist;                      // List of neighbour ids
  FLOAT draux[ndim];                  // Relative distance vector
  FLOAT drsqd;                        // Distance squared
  FLOAT hrangesqdi;                   // Gather kernel range (squared)
  FLOAT hrangesqdj;                   // Scatter kernel range (squared)
  FLOAT rp[ndim];                     // Position of current particle
  FLOAT *dr;                          // Array of neib. position vectors
  FLOAT *drmag;                       // Array of neib. distances
  FLOAT *invdrmag;                    // Array of neib. inverse distances
  struct SphParticle<ndim> *neibpart; // Local copies of neighbour particles
  ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sph_gen);

  debug2("[GodunovSphBruteForce::UpdateAllSphDudt]");

  // The potential number of neighbours is given by ALL the particles
  Nneib = Ntot;

  // Allocate memory for storing neighbour ids and position data
  neiblist = new int[Ntot];
  dr = new FLOAT[ndim*Ntot];
  drmag = new FLOAT[Ntot];
  invdrmag = new FLOAT[Ntot];
  neibpart = new SphParticle<ndim>[Ntot];

  for (j=0; j<Ntot; j++) neibpart[j] = sphdata[j];
  for (j=0; j<Ntot; j++) neibpart[j].dudt = (FLOAT) 0.0;


  // Compute smoothing lengths of all SPH particles
  //---------------------------------------------------------------------------
  for (i=0; i<Nsph; i++) {
    for (k=0; k<ndim; k++) rp[k] = sphdata[i].r[k];
    hrangesqdi = pow(kernfac*kernp->kernrange*sphdata[i].h,2);
    Nneib = 0;

    // Compute distances and the reciprical between the current particle 
    // and all neighbours here
    //-------------------------------------------------------------------------
    for (j=0; j<Ntot; j++) {
      hrangesqdj = pow(kernfac*kernp->kernrange*sphdata[j].h,2);
      for (k=0; k<ndim; k++) draux[k] = sphdata[j].r[k] - rp[k];
      drsqd = DotProduct(draux,draux,ndim);
      if ((drsqd < hrangesqdi || drsqd < hrangesqdj) &&
	  ((j < i && !sphdata[j].active) || j > i)) {
    	neiblist[Nneib] = j;
    	drmag[Nneib] = sqrt(drsqd);
    	invdrmag[Nneib] = (FLOAT) 1.0/(drmag[Nneib] + small_number);
    	for (k=0; k<ndim; k++) dr[Nneib*ndim + k] = draux[k]*invdrmag[Nneib];
    	Nneib++;
      }
    }
    //-------------------------------------------------------------------------

    // Compute all SPH hydro forces
    sph->ComputeSphNeibDudt(i,Nneib,neiblist,drmag,invdrmag,dr,
			    sphdata[i],neibpart);

  }
  //---------------------------------------------------------------------------


  // Now add all active neighbour contributions to the main arrays
  for (j=0; j<Ntot; j++) {
    if (neibpart[j].active) {
      sphdata[j].dudt += neibpart[j].dudt;
    }
  }

  delete[] neibpart;
  delete[] invdrmag;
  delete[] drmag;
  delete[] dr;
  delete[] neiblist;


  return;
}



template class GodunovSphBruteForce<1,GodunovSphParticle>;
template class GodunovSphBruteForce<2,GodunovSphParticle>;
template class GodunovSphBruteForce<3,GodunovSphParticle>;
