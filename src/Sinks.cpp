//=============================================================================
//  Sinks.cpp
//  All routines for creating new sinks and accreting gas and updating all 
//  sink particle propterties.
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
#include "Sph.h"
#include "Nbody.h"
#include "Sinks.h"
#include "Debug.h"
#include "Exception.h"
#include "InlineFuncs.h"
using namespace std;



//=============================================================================
//  Sinks::Sinks()
/// Sinks class constructor
//=============================================================================
template <int ndim>
Sinks<ndim>::Sinks()
{
  allocated_memory = false;
  Nsink = 0;
}



//=============================================================================
//  Sinks::~Sinks()
/// Sinks class destructor
//=============================================================================
template <int ndim>
Sinks<ndim>::~Sinks()
{
}



//=============================================================================
//  Sinks::AllocateMemory
/// Allocate all memory required for storing sink particle data.
//=============================================================================
template <int ndim>
void Sinks<ndim>::AllocateMemory(int N)
{
  if (N > Nsinkmax) {
    if (allocated_memory) DeallocateMemory();
    Nsinkmax = N;
    sink = new SinkParticle<ndim>[Nsinkmax];
    allocated_memory = true;
  }

  return;
}



//=============================================================================
//  Sinks::DeallocateMemory
/// Deallocate all sink particle arrays
//=============================================================================
template <int ndim>
void Sinks<ndim>::DeallocateMemory(void)
{
  if (allocated_memory) delete[] sink;
  allocated_memory = false;

  return;
}



//=============================================================================
//  Sinks::SearchForNewSinkParticles
/// Searches through all SPH particles for new sink particle candidates, and 
/// if a particle satisfies all tests, then a sink is created.
//=============================================================================
template <int ndim>
void Sinks<ndim>::SearchForNewSinkParticles
(int n,                             ///< [in] Current integer time
 Sph<ndim> *sph,                    ///< [inout] Object containing SPH ptcls
 Nbody<ndim> *nbody)                ///< [inout] Object containing star ptcls
{
  bool sink_flag;                   // Flag if particle is to become a sink
  int i;                            // Particle counter
  int isink;                        // i.d. of SPH particle to form sink from
  int k;                            // Dimension counter
  int s;                            // Sink counter
  FLOAT dr[ndim];                   // Relative position vector
  FLOAT drsqd;                      // Distance squared
  FLOAT rho_max = 0.0;              // Maximum density of sink candidates

  debug2("[Sinks::SearchForNewSinkParticles]");


  // Continuous loop to search for new sinks.  If a new sink is found, then 
  // repeat entire process to search for other sinks on current timestep. 
  // If no sinks are found, then exit and return to main program.
  // ==========================================================================
  do {
    isink = -1;

    // Loop over all SPH particles finding the particle with the highest 
    // density that obeys all of the formation criteria, if any do.
    // ------------------------------------------------------------------------
    for (i=0; i<sph->Nsph; i++) {
      sink_flag = true;

      // Only consider SPH particles located at a local potential minimum
      if (!sph->sphdata[i].potmin) continue;

      //cout << "Found candidate : " << i << "    " << sph->sphdata[i].r[0] 
      //   << "    " << sph->sphdata[i].r[1] << endl;

      // If density of SPH particle is too low, skip to next particle
      if (sph->sphdata[i].rho < rho_sink) continue;

      // Make sure candidate particle is at the end of its current timestep
      //if (n%sph->sphdata[i].nstep != 0) continue;

      // If SPH particle neighbours a nearby sink, skip to next particle
      for (s=0; s<Nsink; s++) {
        for (k=0; k<ndim; k++) 
	  dr[k] = sph->sphdata[i].r[k] - sink[s].star->r[k];
        drsqd = DotProduct(dr,dr,ndim);
        if (drsqd < pow(sink_radius*sph->sphdata[i].h + sink[s].radius,2))
          sink_flag = false;
      }

      // If candidate particle has passed all the tests, then check if it is 
      // the most dense candidate.  If yes, record the particle id and density
      if (sink_flag && sph->sphdata[i].rho > rho_max) {
        isink = i;
        rho_max = sph->sphdata[i].rho;
      }

    }
    // ------------------------------------------------------------------------


    // If all conditions have been met, then create a new sink particle
    if (isink != -1) {
      CreateNewSinkParticle(isink,sph,nbody);
      cout << "Found sink particle : " << isink << "    " 
           << sph->sphdata[isink].rho << "     " << rho_sink << endl;
      //exit(0);
    }


  } while (isink != -1);
  // ==========================================================================

  return;
}



//=============================================================================
//  Sinks::CreateNewSinkParticle
/// Create a new sink particle from specified SPH particle 'isink'.
//=============================================================================
template <int ndim>
void Sinks<ndim>::CreateNewSinkParticle
(int isink,                         ///< [in] i.d. of SPH ptcl
 Sph<ndim> *sph,                    ///< [inout] Object containing SPH ptcls
 Nbody<ndim> *nbody)                ///< [inout] Object containing star ptcls
{
  int k;                            // Dimension counter

  debug2("[Sinks::CreateNewSinkParticle]");

  // If we've reached the maximum number of sinks, then throw exception
  if (Nsink == Nsinkmax || nbody->Nstar == nbody->Nstarmax) {
    cout << "Run out of memory : " << Nsink << "    " << Nsinkmax << endl;
    exit(0);
  }

  // First create new star and set pointer
  sink[Nsink].star = &nbody->stardata[nbody->Nstar];

  // If we have space in main arrays, then create sink
  sink[Nsink].star->m = sph->sphdata[isink].m;
  sink[Nsink].star->radius = sph->kernp->kernrange*sph->sphdata[isink].h;
  sink[Nsink].radius = sph->kernp->kernrange*sph->sphdata[isink].h;
  sink[Nsink].star->h = sph->sphdata[isink].h;
  sink[Nsink].star->invh = 1.0/sph->sphdata[isink].h;
  sink[Nsink].star->nstep = sph->sphdata[isink].nstep;
  sink[Nsink].star->level = sph->sphdata[isink].level;
  for (k=0; k<ndim; k++) sink[Nsink].star->r[k] = sph->sphdata[isink].r[k];
  for (k=0; k<ndim; k++) sink[Nsink].star->v[k] = sph->sphdata[isink].v[k];
  for (k=0; k<ndim; k++) sink[Nsink].star->a[k] = sph->sphdata[isink].a[k];
  for (k=0; k<ndim; k++) sink[Nsink].star->adot[k] = 0.0;
  for (k=0; k<ndim; k++) sink[Nsink].star->r0[k] = sph->sphdata[isink].r0[k];
  for (k=0; k<ndim; k++) sink[Nsink].star->v0[k] = sph->sphdata[isink].v0[k];
  for (k=0; k<ndim; k++) sink[Nsink].star->a0[k] = sph->sphdata[isink].a0[k];
  for (k=0; k<ndim; k++) sink[Nsink].star->adot0[k] = 0.0;
  for (k=0; k<3; k++) sink[Nsink].angmom[k] = 0.0;

  // Increment star and sink counters
  nbody->Nstar++;
  Nsink++;

  return;
}



//=============================================================================
//  Sinks::AcceteMassToSinks
/// Identify all SPH particles inside sinks and accrete some fraction (or all) 
/// of the gas mass to the sinks if selected accretion criteria are satisfied.
//=============================================================================
template <int ndim>
void Sinks<ndim>::AccreteMassToSinks
(Sph<ndim> *sph,                    ///< [inout] Object containing SPH ptcls
 Nbody<ndim> *nbody)                ///< [inout] Object containing star ptcls
{
  int i,j,k;                        // Particle and dimension counters
  int Ndead = 0;                    // No. of 'dead' (i.e. accreted) particles
  int Nlist = 0;                    // Max. no of gas particles inside sink
  int Nlisttot = 0;                 // Total number of gas ptcls inside sinks
  int Nneib;                        // ..
  int s;                            // Sink counter
  int saux;                         // Aux. sink i.d.
  int *accrete_list;                // ..
  int *deadlist;                    // List of 'dead' particles
  int *ilist;                       // List of particle ids
  FLOAT asqd;                       // ..
  FLOAT dr[ndim];                   // Relative position vector
  FLOAT drmag;                      // Distance
  FLOAT drsqd;                      // Distance squared
  FLOAT dt;                         // Sink/star timestep
  FLOAT dv[ndim];                   // Relative velocity vector
  FLOAT efrac;                      // ..
  FLOAT macc;                       // ..
  FLOAT macc_temp;                  // ..
  FLOAT mold;                       // ..
  FLOAT mtemp;                      // ..
  FLOAT rold[ndim];                 // ..
  FLOAT timestep;                   // ..
  FLOAT vold[ndim];                 // ..
  FLOAT wnorm;                      // Kernel normalisation factor
  FLOAT *rsqdlist;                  // Array of particle-sink distances
  SinkParticle<ndim> s1;            // Local reference to sink

  debug2("[Sinks::AccreteMassToSinks]");

  // Allocate local memory and initialise values
  accrete_list = new int[sph->Ntot];
  for (i=0; i<sph->Ntot; i++) accrete_list[i] = -1;
  for (s=0; s<Nsinkmax; s++) sink[s].Ngas = 0;
  

  // Determine which sink each SPH particle accretes to.  If none, flag -1
  // --------------------------------------------------------------------------
  for (i=0; i<sph->Nsph; i++) {
    saux = -1;
    for (s=0; s<Nsink; s++) {
      for (k=0; k<ndim; k++) dr[k] = sph->sphdata[i].r[k] - sink[s].star->r[k];
      drsqd = DotProduct(dr,dr,drsqd);
      if (drsqd <= sink[s].radius*sink[s].radius) saux = s;
    }
    accrete_list[i] = saux;
    if (saux != -1) {
      sink[saux].Ngas++;
      Nlisttot++;
      Nlist = max(Nlist,sink[saux].Ngas);
    }
  }

  // If there are no particles inside any sink, deallocate memory and return
  if (Nlist == 0) {
    delete[] accrete_list;
    return;
  }

  // Otherwise, allocate additional memory and proceed to accrete mass
  deadlist = new int[Nlisttot];
  ilist = new int[Nlist];
  rsqdlist = new FLOAT[Nlist];


  // Calculate the accretion timescale and the total mass accreted from all 
  // particles for each sink.
  // ==========================================================================
  for (s=0; s<Nsink; s++) {
    s1 = sink[s];

    // Skip sink particle unless it's at the beginning of its current step
    if (s1.Ngas == 0) continue;

    // Initialise all variables for current sink
    s1.menc  = 0.0;
    s1.trad  = 0.0;
    s1.tvisc = 1.0;
    s1.ketot = 0.0;
    s1.gpetot = 0.0;
    Nneib = 0;
    wnorm = 0.0;

    // Calculate distances (squared) from sink to all neighbouring particles
    for (i=0; i<sph->Nsph; i++) {
      if (accrete_list[i] == s) {
	for (k=0; k<ndim; k++) dr[k] = sph->sphdata[i].r[k] - s1.star->r[k];
        drsqd = DotProduct(dr,dr,ndim);
	if (drsqd < s1.radius*s1.radius) {
	  ilist[Nneib] = i;
	  rsqdlist[Nneib] = drsqd;
	  Nneib++;
	}
      }
    }

    // Double-check that numbers add up here
    if (Nneib != s1.Ngas) cout << "WTF?? : " << Nneib << "   " << s1.Ngas << endl;

    // Sort particle ids by increasing distance from the sink
    Heapsort(Nneib,ilist,rsqdlist);


    // Calculate all important quantities (e.g. energy contributions) due to 
    // all particles inside the sink
    // ------------------------------------------------------------------------
    for (j=0; j<Nneib; j++) {
      i = ilist[j];

      for (k=0; k<ndim; k++) dr[k] = sph->sphdata[i].r[k] - s1.star->r[k];
      drsqd = DotProduct(dr,dr,ndim);
      drmag = sqrt(drsqd) + small_number;
      for (k=0; k<ndim; k++) dr[k] /= drmag;

      s1.menc += sph->sphdata[i].m;
      wnorm += sph->sphdata[i].m*sph->kernp->w0(drmag*s1.star->invh)*
	pow(s1.star->invh,ndim)*sph->sphdata[i].invrho;

      // Sum total grav. potential energy of all particles inside sink
      s1.gpetot += 0.5*sph->sphdata[i].m*(s1.star->m + s1.menc)*
	s1.star->invh*sph->kernp->wpot(drmag*s1.star->invh);

      // Compute rotational component of kinetic energy
      for (k=0; k<ndim; k++) dv[k] = sph->sphdata[i].v[k] - s1.star->v[k];
      for (k=0; k<ndim; k++) dv[k] -= DotProduct(dv,dr,ndim)*dr[k];
      s1.rotketot += sph->sphdata[i].m*DotProduct(dv,dv,ndim)*
	sph->kernp->w0(drmag*s1.star->invh)*pow(s1.star->invh,ndim)*sph->sphdata[i].invrho;

      // Compute total kinetic energy
      for (k=0; k<ndim; k++) dv[k] = sph->sphdata[i].v[k] - s1.star->v[k];
      s1.ketot += sph->sphdata[i].m*DotProduct(dv,dv,ndim)*
	sph->kernp->w0(drmag*s1.star->invh)*pow(s1.star->invh,ndim)*sph->sphdata[i].invrho;

      // Add contributions to average timescales from particles
      s1.tvisc *= pow(sqrt(drmag)/sph->sphdata[i].sound/
		      sph->sphdata[i].sound,sph->sphdata[i].m);
      s1.trad += fabs(4.0*pi*drsqd*sph->sphdata[i].m*DotProduct(dv,dr,ndim)*
		      sph->kernp->w0(drmag*s1.star->invh)*pow(s1.star->invh,ndim));  
    }


    // Normalise sums
    s1.ketot *= 0.5*s1.menc/wnorm;
    s1.rotketot *= 0.5*s1.menc/wnorm;


    // Calculate the sink accretion timescale and the total amount of mass 
    // accreted by sink s this timestep.  If the contained mass is greater 
    // than the maximum allowed, accrete the excess mass.  Otherwise, 
    // accrete a small amount based on freefall/viscous timescale.
    // ------------------------------------------------------------------------
    if (smooth_accretion == 1) {
      efrac    = min(2.0*s1.rotketot/s1.gpetot,1.0);
      s1.tvisc = (sqrt(s1.star->m + s1.menc)*pow(s1.tvisc,1.0/s1.menc))
	/alpha_ss;
      s1.trad  = s1.menc / s1.trad;
      s1.trot  = twopi*sqrt(pow(s1.radius,3)/(s1.menc + s1.star->m));
          
      // Finally calculate accretion timescale and mass accreted
      s1.taccrete = pow(s1.trad,1.0 - efrac)*pow(s1.tvisc,efrac);
      if (s1.menc > s1.mmax) s1.taccrete *= pow(s1.mmax/s1.menc,2);
      dt = (FLOAT) s1.star->nstep*timestep;
      macc = s1.menc*max(1.0 - exp(-dt/s1.taccrete),0.0);
    }
    else {
      macc = s1.menc;
    }


    // Now accrete particles to sink
    // ------------------------------------------------------------------------
    macc_temp = macc;
    for (k=0; k<ndim; k++) rold[k] = s1.star->r[k];
    for (k=0; k<ndim; k++) vold[k] = s1.star->v[k];
    mold = s1.star->m;

    for (k=0; k<ndim; k++) s1.star->r[k] *= s1.star->m;
    for (k=0; k<ndim; k++) s1.star->v[k] *= s1.star->m;
    for (k=0; k<ndim; k++) s1.star->a[k] *= s1.star->m;

    // Loop over all neighbouring particles
    // ------------------------------------------------------------------------
    for (j=0; j<Nneib; j++) {
      i = ilist[j];
      mtemp = min(sph->sphdata[i].m,macc_temp);
      dt = sph->sphdata[i].dt;

      // Special conditions for total particle accretion
      if (smooth_accretion == 0 || 
	  sph->sphdata[i].m - mtemp < smooth_accrete_frac*sph->mmean ||
	  dt < smooth_accrete_dt*s1.trot)
	mtemp = sph->sphdata[i].m;
      macc_temp -= mtemp;

      // Now accrete COM quantities to sink particle
      s1.star->m += mtemp;
      for (k=0; k<ndim; k++) s1.star->r[k] += mtemp*sph->sphdata[i].r[k];
      for (k=0; k<ndim; k++) s1.star->v[k] += mtemp*sph->sphdata[i].v[k];
      for (k=0; k<ndim; k++) s1.star->a[k] += mtemp*sph->sphdata[i].a[k];
      s1.utot += mtemp*sph->sphdata[i].u;

      // If we've reached/exceeded the mass limit, do not include more ptcls
      if (macc_temp < small_number) break;

    }

    // Normalise COM quantities
    for (k=0; k<ndim; k++) s1.star->r[k] /= s1.star->m;
    for (k=0; k<ndim; k++) s1.star->v[k] /= s1.star->m;
    for (k=0; k<ndim; k++) s1.star->a[k] /= s1.star->m;

    // Calculate angular momentum of old COM around new COM
    for (k=0; k<ndim; k++) dr[k] = rold[k] - s1.star->r[k];
    for (k=0; k<ndim; k++) dv[k] = vold[k] - s1.star->v[k];
    s1.angmom[2] += mold*(dr[0]*dv[1] - dr[1]*dv[0]);
    if (ndim == 3) {
      s1.angmom[0] += mold*(dr[1]*dv[2] - dr[2]*dv[1]);
      s1.angmom[1] += mold*(dr[2]*dv[0] - dr[0]*dv[2]);
    }

    // Now add angular momentum contribution of individual SPH particles
    // ------------------------------------------------------------------------
    for (j=0; j<Nneib; j++) {
      i = ilist[j];
      mtemp = min(sph->sphdata[i].m,macc_temp);
      dt = sph->sphdata[i].dt;

      // Special conditions for total particle accretion
      if (smooth_accretion == 0 || 
	  sph->sphdata[i].m - mtemp < smooth_accrete_frac*sph->mmean ||
	  dt < smooth_accrete_dt*s1.trot)
	mtemp = sph->sphdata[i].m;
      macc -= mtemp;

      // Calculate angular momentum of old COM around new COM
      for (k=0; k<ndim; k++) dr[k] = sph->sphdata[i].r[k] - s1.star->r[k];
      for (k=0; k<ndim; k++) dv[k] = sph->sphdata[i].v[k] - s1.star->v[k];
      s1.angmom[2] += mtemp*(dr[0]*dv[1] - dr[1]*dv[0]);
      if (ndim == 3) {
	s1.angmom[0] += mtemp*(dr[1]*dv[2] - dr[2]*dv[1]);
	s1.angmom[1] += mtemp*(dr[2]*dv[0] - dr[0]*dv[2]);
      }

      // If we've reached/exceeded the mass limit, do not include more ptcls
      if (macc < small_number) break;

    }
    
    // Calculate internal sink timestep here
    asqd = DotProduct(s1.star->a,s1.star->a,ndim);
    s1.star->dt_internal = 0.4*sqrt(s1.radius/(sqrt(asqd) + small_number));

    // Copy back sink
    sink[s] = s1;

  }
  // ==========================================================================


  return;
}



// Create template class instances of the main SphSimulation object for
// each dimension used (1, 2 and 3)
template class Sinks<1>;
template class Sinks<2>;
template class Sinks<3>;
