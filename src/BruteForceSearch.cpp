//=============================================================================
//  BruteForceSearch.cpp
//  Contains all routines for generating SPH neighbour lists using 
//  brute-force (i.e. direct summation over all particles).
//=============================================================================



#include <iostream>
#include <math.h>
#include "SphNeighbourSearch.h"
#include "Sph.h"
#include "Parameters.h"
#include "SphParticle.h"
#include "Debug.h"
#include "InlineFuncs.h"
using namespace std;



//=============================================================================
//  BruteForceSearch::BruteForceSearch
/// BruteForceSearch class constructor
//=============================================================================
template <int ndim>
BruteForceSearch<ndim>::BruteForceSearch()
{
}



//=============================================================================
//  BruteForceSearch::~BruteForceSearch
/// BruteForceSearch class destructor
//=============================================================================
template <int ndim>
BruteForceSearch<ndim>::~BruteForceSearch()
{
}



//=============================================================================
//  BruteForceSearch::UpdateTree
/// For Brute Force neighbour searching, there is no tree to construct.
//=============================================================================
template <int ndim>
void BruteForceSearch<ndim>::UpdateTree(Sph<ndim> *sph, Parameters &simparams)
{
  return;
}



//=============================================================================
//  BruteForceSearch::UpdateAllSphProperties
/// Routine for computing SPH properties (smoothing lengths, densities and 
/// forces) for all active SPH particle using neighbour lists generated 
/// using brute force (i.e. direct summation).
//=============================================================================
template <int ndim>
void BruteForceSearch<ndim>::UpdateAllSphProperties
(Sph<ndim> *sph                     ///< [inout] Pointer to SPH object
 )
{
  int i,j,k;                        // Particle and dimension counters
  int okflag;                       // Flag valid smoothing length
  int Nneib;                        // No. of neighbours
  int Nfar;                         // No. of 'far' neighbours
  int *neiblist;                    // List of neighbour ids
  FLOAT dr[ndim];                   // Relative distance vector
  FLOAT drsqdaux;                   // Distance squared
  FLOAT rp[ndim];                   // Position of current particle
  FLOAT *drsqd;                     // Array of neib. distances (sqd)
  FLOAT *m;                         // Array of neib. position vectors
  FLOAT *mu;                        // Array of neib. mass*u values

  debug2("[BruteForceSearch::UpdateAllSphProperties]");

  // Store masses in separate array
  m = new FLOAT[sph->Ntot];
  mu = new FLOAT[sph->Ntot];
  for (i=0; i<sph->Ntot; i++) m[i] = sph->sphdata[i].m;
  for (i=0; i<sph->Ntot; i++) mu[i] = sph->sphdata[i].m*sph->sphdata[i].u;

  // Create parallel threads
  // ==========================================================================
#pragma omp parallel default(shared) private(dr,drsqd,i,j,k,neiblist,okflag,rp)
  {
    neiblist = new int[sph->Ntot];
    drsqd = new FLOAT[sph->Ntot];

    // Compute smoothing lengths of all SPH particles
    // ------------------------------------------------------------------------
#pragma omp for
    for (i=0; i<sph->Nsph; i++) {
      for (k=0; k<ndim; k++) rp[k] = sph->sphdata[i].r[k];

      // Compute distances and the reciprical between the current particle 
      // and all neighbours here
      // ----------------------------------------------------------------------
      for (j=0; j<sph->Ntot; j++) { 
    	neiblist[j] = j;
    	for (k=0; k<ndim; k++) dr[k] = sph->sphdata[j].r[k] - rp[k];
    	drsqd[j] = DotProduct(dr,dr,ndim);
      }
      // ----------------------------------------------------------------------

      // Compute all SPH gather properties
      okflag = sph->ComputeH(i,sph->Ntot,big_number,
                             m,mu,drsqd,sph->sphdata[i]);
  
    }
    // ------------------------------------------------------------------------

    delete[] drsqd;
    delete[] neiblist;

  }
  // ==========================================================================

  delete[] mu;
  delete[] m;

  return;
}



//=============================================================================
//  BruteForceSearch::UpdateAllSphHydroForces
/// Routine for computing SPH properties (smoothing lengths, densities and 
/// forces) for all active SPH particle using neighbour lists generated 
/// using brute force (i.e. direct summation).
//=============================================================================
template <int ndim>
void BruteForceSearch<ndim>::UpdateAllSphHydroForces
(Sph<ndim> *sph                       ///< [inout] Pointer to SPH object
)
{
  int i,j,k;                          // Particle and dimension counters
  int okflag;                         // Flag valid smoothing length
  int Nneib;                          // No. of neighbours
  int *neiblist;                      // List of neighbour ids
  FLOAT draux[ndim];                  // Relative distance vector
  FLOAT drsqd;                        // Distance squared
  FLOAT hrangesqdi;                   // Gather kernel extent (squared)
  FLOAT hrangesqdj;                   // Scatter kernel extent (squared)
  FLOAT rp[ndim];                     // Position of current particle
  FLOAT *dr;                          // Array of neib. position vectors
  FLOAT *drmag;                       // Array of neib. distances
  FLOAT *invdrmag;                    // Array of neib. inverse distances
  struct SphParticle<ndim> *neibpart; // Local copies of neib. particles

  debug2("[BruteForceSearch::UpdateAllSphForces]");

  // The potential number of neighbours is given by ALL the particles
  Nneib = sph->Ntot;

  // Allocate memory for storing neighbour ids and position data
  neiblist = new int[sph->Ntot];
  dr = new FLOAT[ndim*sph->Ntot];
  drmag = new FLOAT[sph->Ntot];
  invdrmag = new FLOAT[sph->Ntot];
  neibpart = new SphParticle<ndim>[sph->Ntot];

  for (j=0; j<sph->Ntot; j++) neibpart[j] = sph->sphdata[j];

  // Compute smoothing lengths of all SPH particles
  // --------------------------------------------------------------------------
  for (i=0; i<sph->Nsph; i++) {
    for (k=0; k<ndim; k++) rp[k] = sph->sphdata[i].r[k];
    hrangesqdi = pow(sph->kernfac*sph->kernp->kernrange*sph->sphdata[i].h,2);
    Nneib = 0;

    // Make local copies of all potential neighbours
    for (j=0; j<sph->Ntot; j++) {
      neibpart[j].div_v = (FLOAT) 0.0;
      neibpart[j].dudt = (FLOAT) 0.0;
      for (k=0; k<ndim; k++) neibpart[j].a[k] = (FLOAT) 0.0;
    }

    // Compute distances and the reciprical between the current particle 
    // and all neighbours here
    // ------------------------------------------------------------------------
    for (j=0; j<sph->Ntot; j++) {
      hrangesqdj = pow(sph->kernfac*sph->kernp->kernrange*sph->sphdata[j].h,2);
      for (k=0; k<ndim; k++) draux[k] = sph->sphdata[j].r[k] - rp[k];
      drsqd = DotProduct(draux,draux,ndim);
      if ((drsqd < hrangesqdi || drsqd < hrangesqdj) &&
	  ((j < i && !sph->sphdata[j].active) || j > i)) {
    	neiblist[Nneib] = j;
    	drmag[Nneib] = sqrt(drsqd);
    	invdrmag[Nneib] = (FLOAT) 1.0/(drmag[Nneib] + small_number);
    	for (k=0; k<ndim; k++) dr[Nneib*ndim + k] = draux[k]*invdrmag[Nneib];
    	Nneib++;
      }
    }
    // ------------------------------------------------------------------------

    // Compute all SPH hydro forces
    sph->ComputeSphHydroForces(i,Nneib,neiblist,drmag,invdrmag,dr,
			      sph->sphdata[i],neibpart);

    // Now add all active neighbour contributions to the main arrays
    for (j=0; j<sph->Ntot; j++) {
      if (neibpart[j].active) {
        for (k=0; k<ndim; k++) sph->sphdata[j].a[k] += neibpart[j].a[k];
        sph->sphdata[j].dudt += neibpart[j].dudt;
        sph->sphdata[j].div_v += neibpart[j].div_v;
      }
    }

  }
  // --------------------------------------------------------------------------

  delete[] neibpart;
  delete[] invdrmag;
  delete[] drmag;
  delete[] dr;
  delete[] neiblist;


  // Compute other important SPH quantities after hydro forces are computed
  if (sph->hydro_forces == 1) {
    for (i=0; i<sph->Nsph; i++) {
      if (sph->sphdata[i].active) 
        sph->ComputePostHydroQuantities(sph->sphdata[i]);
    }
  }

  return;
}



//=============================================================================
//  BruteForceSearch::UpdateAllSphGravityProperties
/// Empty function for now
//=============================================================================
template <int ndim>
void BruteForceSearch<ndim>::UpdateAllSphForces
(Sph<ndim> *sph                     ///< Pointer to SPH object
)
{
  return;
}



//=============================================================================
//  BruteForceSearch::UpdateAllSphGravForces
/// Routine for computing SPH properties (smoothing lengths, densities and 
/// forces) for all active SPH particle using neighbour lists generated 
/// using brute force (i.e. direct summation).
//=============================================================================
template <int ndim>
void BruteForceSearch<ndim>::UpdateAllSphGravForces
(Sph<ndim> *sph                       ///< [inout] Pointer to SPH object
)
{
  int i,j,k;                          // Particle and dimension counters
  int okflag;                         // Flag valid smoothing length
  int Nneib;                          // No. of neighbours
  int *neiblist;                      // List of neighbour ids
  FLOAT draux[ndim];                  // Relative distance vector
  FLOAT drsqd;                        // Distance squared
  FLOAT hrangesqdi;                   // Gather kernel extent (squared)
  FLOAT hrangesqdj;                   // Scatter kernel extent (squared)
  FLOAT rp[ndim];                     // Position of current particle
  FLOAT *dr;                          // Array of neib. position vectors
  FLOAT *drmag;                       // Array of neib. distances
  FLOAT *invdrmag;                    // Array of neib. inverse distances
  FLOAT *m;                           // Array of neib. masses
  struct SphParticle<ndim> *neibpart; // Local copies of neib. particles

  debug2("[BruteForceSearch::UpdateAllSphGravForces]");

  // The potential number of neighbours is given by ALL the particles
  Nneib = sph->Ntot;

  // Allocate memory for storing neighbour ids and position data
  neiblist = new int[sph->Ntot];
  dr = new FLOAT[ndim*sph->Ntot];
  drmag = new FLOAT[sph->Ntot];
  invdrmag = new FLOAT[sph->Ntot];
  m = new FLOAT[sph->Ntot];
  neibpart = new SphParticle<ndim>[sph->Ntot];

  for (j=0; j<sph->Ntot; j++) neibpart[j] = sph->sphdata[j];

  // Compute smoothing lengths of all SPH particles
  // --------------------------------------------------------------------------
  for (i=0; i<sph->Nsph; i++) {
    for (k=0; k<ndim; k++) rp[k] = sph->sphdata[i].r[k];
    hrangesqdi = pow(sph->kernfac*sph->kernp->kernrange*sph->sphdata[i].h,2);
    Nneib = 0;

    // Make local copies of all potential neighbours
    for (j=0; j<sph->Ntot; j++) {
      neibpart[j].gpot = 0.0;
      for (k=0; k<ndim; k++) neibpart[j].agrav[k] = (FLOAT) 0.0;
    }

    // Compute distances and the reciprical between the current particle 
    // and all neighbours here
    // ------------------------------------------------------------------------
    for (j=0; j<sph->Ntot; j++) {

      for (k=0; k<ndim; k++) draux[k] = sph->sphdata[j].r[k] - rp[k];
      drsqd = DotProduct(draux,draux,ndim);
      neiblist[Nneib] = j;
      drmag[Nneib] = sqrt(drsqd);
      invdrmag[Nneib] = (FLOAT) 1.0/(drmag[Nneib] + small_number);
      for (k=0; k<ndim; k++) dr[Nneib*ndim + k] = draux[k]*invdrmag[Nneib];
      Nneib++;

    }
    // ------------------------------------------------------------------------

    // Compute all SPH hydro forces
    //sph->ComputeSphNeibGravForces(i,Nneib,neiblist,drmag,invdrmag,dr,
    //				  sph->sphdata[i],neibpart);

    // Now add all active neighbour contributions to the main arrays
    //for (j=0; j<sph->Ntot; j++) {
    //  if (neibpart[j].active) {
    //sph->sphdata[i].gpot += neibpart[j].gpot;
    //  for (k=0; k<ndim; k++) 
    //  sph->sphdata[j].agrav[k] += neibpart[j].agrav[k];
    //}
    //}

  }
  // --------------------------------------------------------------------------

  delete[] neibpart;
  delete[] invdrmag;
  delete[] drmag;
  delete[] dr;
  delete[] neiblist;

  return;
}



//=============================================================================
//  BruteForceSearch::UpdateAllSphDerivatives
/// ..
//=============================================================================
template <int ndim>
void BruteForceSearch<ndim>::UpdateAllSphDerivatives
(Sph<ndim> *sph                       ///< Pointer to SPH object
)
{
  int i,j,k;                          // Particle and dimension counters
  int okflag;                         // Flag valid smoothing length
  int Nneib;                          // No. of neighbours
  int Nfar;                           // No. of 'far' neighbours
  int *neiblist;                      // List of neighbour ids
  FLOAT draux[ndim];                  // Relative distance vector
  FLOAT drsqd;                        // Distance squared
  FLOAT hrangesqd;                    // Kernel extent (squared)
  FLOAT rp[ndim];                     // Position of current particle
  FLOAT *dr;                          // Array of neib. position vectors
  FLOAT *drmag;                       // Array of neib. distances
  FLOAT *invdrmag;                    // Array of neib. inverse distances
  struct SphParticle<ndim> *neibpart; // Local copies of neib. particles

  debug2("[BruteForceSearch::UpdateAllSphForces]");

  // The potential number of neighbours is given by ALL the particles
  Nneib = sph->Ntot;

  // Allocate memory for storing neighbour ids and position data
  neiblist = new int[sph->Ntot];
  dr = new FLOAT[ndim*sph->Ntot];
  drmag = new FLOAT[sph->Ntot];
  invdrmag = new FLOAT[sph->Ntot];
  neibpart = new SphParticle<ndim>[sph->Ntot];

  // Record local copies of (all) neighbour properties
  for (j=0; j<sph->Ntot; j++) neibpart[j] = sph->sphdata[j];

  // Compute smoothing lengths of all SPH particles
  // --------------------------------------------------------------------------
  for (i=0; i<sph->Nsph; i++) {
    for (k=0; k<ndim; k++) rp[k] = sph->sphdata[i].r[k];
    hrangesqd = pow(sph->kernp->kernrange*sph->sphdata[i].h,2);
    Nneib = 0;

    // Compute distances and the reciprical between the current particle 
    // and all neighbours here
    // ------------------------------------------------------------------------
    for (j=0; j<sph->Ntot; j++) {
      for (k=0; k<ndim; k++) draux[k] = sph->sphdata[j].r[k] - rp[k];
      drsqd = DotProduct(draux,draux,ndim);
      if (drsqd < hrangesqd) {
    	neiblist[Nneib] = j;
    	drmag[Nneib] = sqrt(drsqd);
    	invdrmag[Nneib] = (FLOAT) 1.0/(drmag[Nneib] + small_number);
    	for (k=0; k<ndim; k++) dr[Nneib*ndim + k] = draux[k]*invdrmag[Nneib];
    	Nneib++;
      }
    }
    // ------------------------------------------------------------------------

    // Compute all SPH hydro forces
    sph->ComputeSphDerivatives(i,Nneib,neiblist,drmag,invdrmag,dr,
			       sph->sphdata[i],neibpart);

  }
  // --------------------------------------------------------------------------

  delete[] neibpart;
  delete[] invdrmag;
  delete[] drmag;
  delete[] dr;
  delete[] neiblist;

  return;
}



//=============================================================================
//  BruteForceSearch::UpdateAllSphDudt
/// Compute the compressional heating rate (dudt) for all active particles.
//=============================================================================
template <int ndim>
void BruteForceSearch<ndim>::UpdateAllSphDudt
(Sph<ndim> *sph                       ///< Pointer to SPH object
)
{
  int i,j,k;                          // Particle and dimension counters
  int okflag;                         // Flag valid smoothing length
  int Nneib;                          // No. of neighbours
  int Nfar;                           // No. of 'far' neighbours
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

  debug2("[BruteForceSearch::UpdateAllSphForces]");

  // The potential number of neighbours is given by ALL the particles
  Nneib = sph->Ntot;

  // Allocate memory for storing neighbour ids and position data
  neiblist = new int[sph->Ntot];
  dr = new FLOAT[ndim*sph->Ntot];
  drmag = new FLOAT[sph->Ntot];
  invdrmag = new FLOAT[sph->Ntot];
  neibpart = new SphParticle<ndim>[sph->Ntot];

  for (j=0; j<sph->Ntot; j++) neibpart[j] = sph->sphdata[j];

  // Compute smoothing lengths of all SPH particles
  // --------------------------------------------------------------------------
  for (i=0; i<sph->Nsph; i++) {
    for (k=0; k<ndim; k++) rp[k] = sph->sphdata[i].r[k];
    hrangesqdi = pow(sph->kernfac*sph->kernp->kernrange*sph->sphdata[i].h,2);
    Nneib = 0;

    // Make local copies of all potential neighbours
     for (j=0; j<sph->Ntot; j++) {
       neibpart[j].dudt = (FLOAT) 0.0;
     }

    // Compute distances and the reciprical between the current particle 
    // and all neighbours here
    // ------------------------------------------------------------------------
    for (j=0; j<sph->Ntot; j++) {
      hrangesqdj = pow(sph->kernfac*sph->kernp->kernrange*sph->sphdata[j].h,2);
      for (k=0; k<ndim; k++) draux[k] = sph->sphdata[j].r[k] - rp[k];
      drsqd = DotProduct(draux,draux,ndim);
      if ((drsqd < hrangesqdi || drsqd < hrangesqdj) &&
	  ((j < i && !sph->sphdata[j].active) || j > i)) {
    	neiblist[Nneib] = j;
    	drmag[Nneib] = sqrt(drsqd);
    	invdrmag[Nneib] = (FLOAT) 1.0/(drmag[Nneib] + small_number);
    	for (k=0; k<ndim; k++) dr[Nneib*ndim + k] = draux[k]*invdrmag[Nneib];
    	Nneib++;
      }
    }
    // ------------------------------------------------------------------------

    // Compute all SPH hydro forces
    sph->ComputeSphNeibDudt(i,Nneib,neiblist,drmag,invdrmag,dr,
			    sph->sphdata[i],neibpart);

    // Now add all active neighbour contributions to the main arrays
    for (j=0; j<sph->Ntot; j++) {
      if (neibpart[j].active) {
        sph->sphdata[j].dudt += neibpart[j].dudt;
      }
    }

  }
  // --------------------------------------------------------------------------

  delete[] neibpart;
  delete[] invdrmag;
  delete[] drmag;
  delete[] dr;
  delete[] neiblist;


  return;
}



template class BruteForceSearch<1>;
template class BruteForceSearch<2>;
template class BruteForceSearch<3>;
