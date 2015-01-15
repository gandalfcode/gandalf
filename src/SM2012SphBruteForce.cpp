//=================================================================================================
//  SM2012SphBruteForce.cpp
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
//=================================================================================================



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



//=================================================================================================
//  SM2012SphBruteForce::SM2012SphBruteForce
/// SM2012SphBruteForce class constructor
//=================================================================================================
template <int ndim, template<int> class ParticleType>
SM2012SphBruteForce<ndim,ParticleType>::SM2012SphBruteForce
(FLOAT kernrangeaux,
 DomainBox<ndim> *boxaux,
 SphKernel<ndim> *kernaux,
 CodeTiming *timingaux):
  BruteForceSearch<ndim,ParticleType>(kernrangeaux,boxaux,kernaux,timingaux)
{
}



//=================================================================================================
//  SM2012SphBruteForce::~SM2012SphBruteForce
/// SM2012SphBruteForce class destructor
//=================================================================================================
template <int ndim, template<int> class ParticleType>
SM2012SphBruteForce<ndim,ParticleType>::~SM2012SphBruteForce()
{
}



//=================================================================================================
//  SM2012SphBruteForce::UpdateAllSphProperties
/// Routine for computing SPH properties (smoothing lengths, densities and forces) for all active
/// SPH particle using neighbour lists generated using brute force (i.e. direct summation).
//=================================================================================================
template <int ndim, template<int> class ParticleType>
void SM2012SphBruteForce<ndim,ParticleType>::UpdateAllSphProperties
(int Nsph,                          ///< [in] No. of SPH particles
 int Ntot,                          ///< [in] No. of SPH + ghost particles
 SphParticle<ndim> *sph_gen,        ///< [inout] Pointer to SPH ptcl array
 Sph<ndim> *sph,                    ///< [in] Pointer to SPH object
 Nbody<ndim> *nbody)                ///< [in] Pointer to N-body object
{
  int i,j,jj,k;                     // Particle and dimension counters
  int Nneib = 0;                    // No. of (non-dead) neighbours
  //int okflag;                       // Flag valid smoothing length
  int *neiblist;                    // List of neighbours
  FLOAT dr[ndim];                   // Relative distance vector
  FLOAT rp[ndim];                   // Position of current particle
  FLOAT *drsqd;                     // Distance squared
  FLOAT *gpot;                      // Array of neib. grav. potentials
  FLOAT *m;                         // Array of neib. position vectors
  FLOAT *mu;                        // Array of neib. mass*u values
  ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sph_gen);

  debug2("[SM2012SphBruteForce::UpdateAllSphProperties]");

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
  //===============================================================================================
#pragma omp parallel default(none) private(dr,drsqd,i,j,jj,k,rp)	\
  shared(gpot,m,mu,nbody,neiblist,Nneib,Nsph,Ntot,sph,sphdata)
  {
    drsqd = new FLOAT[Ntot];

    // Compute smoothing lengths of all SPH particles
    //---------------------------------------------------------------------------------------------
#pragma omp for
    for (i=0; i<Nsph; i++) {

      // Skip over inactive particles
      if (!sphdata[i].active || sphdata[i].itype == dead) continue;

      for (k=0; k<ndim; k++) rp[k] = sphdata[i].r[k];

      // Compute distances and the reciprical between the current particle
      // and all neighbours here
      //-------------------------------------------------------------------------------------------
      for (jj=0; jj<Nneib; jj++) {
        j = neiblist[jj];
        for (k=0; k<ndim; k++) dr[k] = sphdata[j].r[k] - rp[k];
        drsqd[jj] = DotProduct(dr,dr,ndim);
      }
      //-------------------------------------------------------------------------------------------

      // Compute all SPH gather properties
      //okflag =
      sph->ComputeH(i,Nneib,big_number,m,mu,drsqd,gpot,sphdata[i],nbody);

    }
    //---------------------------------------------------------------------------------------------

    delete[] drsqd;

  }
  //===============================================================================================

  delete[] neiblist;
  delete[] mu;
  delete[] m;
  delete[] gpot;

  return;
}



template class SM2012SphBruteForce<1,SM2012SphParticle>;
template class SM2012SphBruteForce<2,SM2012SphParticle>;
template class SM2012SphBruteForce<3,SM2012SphParticle>;
