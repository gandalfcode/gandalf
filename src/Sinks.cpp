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
           << sph->sphdata[isink].rho << endl;
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
  int i;                            // SPH particle counter
  int k;                            // Dimension counter
  int deadlist[1];                  // List of 'dead' particles
  FLOAT dr[ndim];                   // Relative position vector
  FLOAT drsqd;                      // Distance squared

  debug2("[Sinks::CreateNewSinkParticle]");

  // If we've reached the maximum number of sinks, then throw exception
  if (Nsink == Nsinkmax || nbody->Nstar == nbody->Nstarmax) {
    cout << "Run out of memory : " << Nsink << "    " << Nsinkmax << endl;
    exit(0);
  }

  // First create new star and set pointer
  sink[Nsink].star = &(nbody->stardata[nbody->Nstar]);
  nbody->nbodydata[nbody->Nnbody] = &(nbody->stardata[nbody->Nstar]);

  // If we have space in main arrays, then create sink
  sink[Nsink].star->m = sph->sphdata[isink].m;
  sink[Nsink].radius = sph->kernp->kernrange*sph->sphdata[isink].h;
  sink[Nsink].star->h = sph->sphdata[isink].h;
  sink[Nsink].star->invh = 1.0/sph->sphdata[isink].h;
  sink[Nsink].star->hfactor = sph->sphdata[i].hfactor;
  sink[Nsink].star->radius = sph->kernp->kernrange*sph->sphdata[isink].h;
  sink[Nsink].star->gpot = sph->sphdata[i].gpot;
  sink[Nsink].star->gpe_internal = 0.0;
  sink[Nsink].star->dt = sph->sphdata[i].dt;
  sink[Nsink].star->nstep = sph->sphdata[isink].nstep;
  sink[Nsink].star->level = sph->sphdata[isink].level;
  sink[Nsink].star->active = sph->sphdata[isink].active;
  sink[Nsink].star->Ncomp = 1;
  for (k=0; k<ndim; k++) sink[Nsink].star->r[k] = sph->sphdata[isink].r[k];
  for (k=0; k<ndim; k++) sink[Nsink].star->v[k] = sph->sphdata[isink].v[k];
  for (k=0; k<ndim; k++) sink[Nsink].star->a[k] = sph->sphdata[isink].a[k];
  for (k=0; k<ndim; k++) sink[Nsink].star->adot[k] = 0.0;
  for (k=0; k<ndim; k++) sink[Nsink].star->r0[k] = sph->sphdata[isink].r0[k];
  for (k=0; k<ndim; k++) sink[Nsink].star->v0[k] = sph->sphdata[isink].v0[k];
  for (k=0; k<ndim; k++) sink[Nsink].star->a0[k] = sph->sphdata[isink].a0[k];
  for (k=0; k<ndim; k++) sink[Nsink].star->adot0[k] = 0.0;
  for (k=0; k<3; k++) sink[Nsink].angmom[k] = 0.0;

  // Calculate total mass inside sink
  sink[Nsink].mmax = 0.0;
  for (i=0; i<sph->Nsph; i++) {
    for (k=0; k<ndim; k++) 
      dr[k] = sink[Nsink].star->r[k] - sph->sphdata[i].r[k];
    drsqd = DotProduct(dr,dr,ndim);
    if (drsqd < pow(sink[Nsink].radius,2)) 
      sink[Nsink].mmax += sph->sphdata[i].m;
  }

  cout << "Created new sink particle : " << isink << "    " << Nsink << endl;
  cout << "radius : " << sink[Nsink].radius << endl;
  cout << "m : " << sink[Nsink].star->m << "    mmax : " << sink[Nsink].mmax << endl;
  cout << "r : " << sink[Nsink].star->r[0] << "   " << sink[Nsink].star->r[0] << endl;

  // Remove SPH particle from main arrays
  deadlist[0] = isink;
  sph->DeleteParticles(1,deadlist);
  
  // Increment star and sink counters
  nbody->Nstar++;
  nbody->Nnbody++;
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
 Nbody<ndim> *nbody,                ///< [inout] Object containing star ptcls
 int n,                             ///< [in] Integer timestep
 DOUBLE timestep,                   ///< [in] Minimum timestep level
 int ssink)
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
  int *ilist2;                      // List of particle ids
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
  FLOAT vold[ndim];                 // ..
  FLOAT wnorm;                      // Kernel normalisation factor
  FLOAT *rsqdlist;                  // Array of particle-sink distances
  SinkParticle<ndim> *s1;           // Local reference to sink
  int *insink;

  debug2("[Sinks::AccreteMassToSinks]");

  insink = new int[sph->Ntot];
  for (i=0; i<sph->Nsph; i++) insink[i] = 0;

  // Allocate local memory and initialise values
  accrete_list = new int[sph->Ntot];
  for (i=0; i<sph->Ntot; i++) accrete_list[i] = -1;
  for (s=0; s<Nsinkmax; s++) sink[s].Ngas = 0;


  FLOAT mcom = 0.0;
  FLOAT rcom[ndim];
  FLOAT vcom[ndim];
  FLOAT acom[ndim];
  for (k=0; k<ndim; k++) rcom[k] = 0.0;
  for (k=0; k<ndim; k++) vcom[k] = 0.0;
  for (k=0; k<ndim; k++) acom[k] = 0.0;
  for (i=0; i<sph->Nsph; i++) {
    mcom += sph->sphdata[i].m;
    for (k=0; k<ndim; k++) rcom[k] += sph->sphdata[i].m*sph->sphdata[i].r[k];
    for (k=0; k<ndim; k++) vcom[k] += sph->sphdata[i].m*sph->sphdata[i].v[k];
    for (k=0; k<ndim; k++) acom[k] += sph->sphdata[i].m*sph->sphdata[i].a[k];
  }
  for (s=0; s<Nsink; s++) {
    mcom += sink[s].star->m;
    for (k=0; k<ndim; k++) rcom[k] += sink[s].star->m*sink[s].star->r[k];
    for (k=0; k<ndim; k++) vcom[k] += sink[s].star->m*sink[s].star->v[k];
    for (k=0; k<ndim; k++) acom[k] += sink[s].star->m*sink[s].star->a[k];
  }
  for (k=0; k<ndim; k++) rcom[k] /= mcom;
  for (k=0; k<ndim; k++) vcom[k] /= mcom;
  for (k=0; k<ndim; k++) acom[k] /= mcom;
  cout << "mcom (before) : " << mcom << endl;
  cout << "rcom (before) : " << rcom[0] << "   " << rcom[1] << "    " << rcom[2] << endl;
  cout << "vcom (before) : " << vcom[0] << "   " << vcom[1] << "    " << vcom[2] << endl;
  cout << "acom (before) : " << acom[0] << "   " << acom[1] << "    " << acom[2] << endl;
  //return;


  // Determine which sink each SPH particle accretes to.  If none, flag -1
  // --------------------------------------------------------------------------
  for (i=0; i<sph->Nsph; i++) {
    saux = -1;
    for (s=0; s<Nsink; s++) {
      for (k=0; k<ndim; k++) dr[k] = sph->sphdata[i].r[k] - sink[s].star->r[k];
      drsqd = DotProduct(dr,dr,ndim);
      if (drsqd <= sink[s].radius*sink[s].radius) saux = s;
      //cout << "WTF?? : " << i << "   " << sqrt(drsqd) << "    " << sink[s].radius << endl;
    }
    accrete_list[i] = saux;
    if (saux != -1) {
      sink[saux].Ngas++;
      Nlisttot++;
      Nlist = max(Nlist,sink[saux].Ngas);
    }
  }

  cout << "Accretion Nlist : " << Nlist << endl;

  // If there are no particles inside any sink, deallocate memory and return
  if (Nlist == 0) {
    delete[] accrete_list;
    return;
  }

  // Otherwise, allocate additional memory and proceed to accrete mass
  deadlist = new int[Nlisttot];
  ilist = new int[Nlisttot];
  ilist2 = new int[Nlisttot];
  rsqdlist = new FLOAT[Nlisttot];


  // Calculate the accretion timescale and the total mass accreted from all 
  // particles for each sink.
  // ==========================================================================
  s = ssink; //for (s=0; s<Nsink; s++) {
    //if (s == 0) continue;
    cout << "Accreting sink " << s << endl;
    cout << "Ngas : " << sink[s].Ngas << endl;
    // Skip sink particle unless it's at the beginning of its current step
    //if (sink[s].Ngas == 0) continue;

    // Initialise all variables for current sink
    sink[s].menc  = 0.0;
    sink[s].trad  = 0.0;
    sink[s].tvisc = 1.0;
    sink[s].ketot = 0.0;
    sink[s].rotketot = 0.0;
    sink[s].gpetot = 0.0;
    Nneib = 0;
    wnorm = 0.0;

    // Calculate distances (squared) from sink to all neighbouring particles
    for (i=0; i<sph->Nsph; i++) {
      if (accrete_list[i] == s) {
	for (k=0; k<ndim; k++) dr[k] = sph->sphdata[i].r[k] - sink[s].star->r[k];
        drsqd = DotProduct(dr,dr,ndim);
	ilist[Nneib] = i;
	rsqdlist[Nneib] = drsqd;
	Nneib++;
	insink[i]++;
	if (insink[i] > 1) {
	  cout << "WTF??  Neighbour more than 1 sink : " << i << "   " << s << endl;
	  exit(0);
	}
      }
    }

    // Double-check that numbers add up here
    if (Nneib != sink[s].Ngas) cout << "WTF?? : " << Nneib << "   " << sink[s].Ngas << endl;
    if (Nneib > Nlist) cout << "ERROR!!!" << Nneib << "   " << Nlist << endl;

    // Sort particle ids by increasing distance from the sink
    //Heapsort(Nneib,ilist2,rsqdlist);


    // Calculate all important quantities (e.g. energy contributions) due to 
    // all particles inside the sink
    // ------------------------------------------------------------------------
    for (j=0; j<Nneib; j++) {
      i = ilist[j]; //i = ilist[ilist2[j]];

      for (k=0; k<ndim; k++) dr[k] = sph->sphdata[i].r[k] - sink[s].star->r[k];
      drsqd = DotProduct(dr,dr,ndim);
      drmag = sqrt(drsqd) + small_number;
      for (k=0; k<ndim; k++) dr[k] /= drmag;

      sink[s].menc += sph->sphdata[i].m;
      wnorm += sph->sphdata[i].m*sph->kernp->w0(drmag*sink[s].star->invh)*
	pow(sink[s].star->invh,ndim)*sph->sphdata[i].invrho;

      continue;

      // Sum total grav. potential energy of all particles inside sink
      sink[s].gpetot += 0.5*sph->sphdata[i].m*(sink[s].star->m + sink[s].menc)*
	sink[s].star->invh*sph->kernp->wpot(drmag*sink[s].star->invh);

      // Compute rotational component of kinetic energy
      for (k=0; k<ndim; k++) dv[k] = sph->sphdata[i].v[k] - sink[s].star->v[k];
      for (k=0; k<ndim; k++) dv[k] -= DotProduct(dv,dr,ndim)*dr[k];
      sink[s].rotketot += sph->sphdata[i].m*DotProduct(dv,dv,ndim)*
	sph->kernp->w0(drmag*sink[s].star->invh)*pow(sink[s].star->invh,ndim)*sph->sphdata[i].invrho;

      // Compute total kinetic energy
      for (k=0; k<ndim; k++) dv[k] = sph->sphdata[i].v[k] - sink[s].star->v[k];
      sink[s].ketot += sph->sphdata[i].m*DotProduct(dv,dv,ndim)*
	sph->kernp->w0(drmag*sink[s].star->invh)*pow(sink[s].star->invh,ndim)*sph->sphdata[i].invrho;

      // Add contributions to average timescales from particles
      sink[s].tvisc *= pow(sqrt(drmag)/sph->sphdata[i].sound/
		      sph->sphdata[i].sound,sph->sphdata[i].m);
      sink[s].trad += fabs(4.0*pi*drsqd*sph->sphdata[i].m*DotProduct(dv,dr,ndim)*
		      sph->kernp->w0(drmag*sink[s].star->invh)*pow(sink[s].star->invh,ndim));  

      //cout << "Adding values : " << j << "    " << i << "   " << drmag 
      //   << "    " << sqrt(rsqdlist[j]) << "     " 
      //   << sink[s].ketot << "    " << "     " << sink[s].rotketot << "    " 
      //   << sink[s].gpetot << "    " << sink[s].tvisc << "    " << sink[s].trad << endl;

    }


    // Normalise sums
    sink[s].ketot *= 0.5*sink[s].menc/wnorm;
    sink[s].rotketot *= 0.5*sink[s].menc/wnorm;


    // Calculate the sink accretion timescale and the total amount of mass 
    // accreted by sink s this timestep.  If the contained mass is greater 
    // than the maximum allowed, accrete the excess mass.  Otherwise, 
    // accrete a small amount based on freefall/viscous timescale.
    // ------------------------------------------------------------------------
    if (smooth_accretion == 1) {
      efrac    = min(2.0*sink[s].rotketot/sink[s].gpetot,1.0);
      sink[s].tvisc = (sqrt(sink[s].star->m + sink[s].menc)*pow(sink[s].tvisc,1.0/sink[s].menc))
	/alpha_ss;
      sink[s].trad  = sink[s].menc / sink[s].trad;
      sink[s].trot  = twopi*sqrt(pow(sink[s].radius,3)/(sink[s].menc + sink[s].star->m));
          
      // Finally calculate accretion timescale and mass accreted
      sink[s].taccrete = pow(sink[s].trad,1.0 - efrac)*pow(sink[s].tvisc,efrac);
      if (sink[s].menc > sink[s].mmax) sink[s].taccrete *= pow(sink[s].mmax/sink[s].menc,2);
      dt = (FLOAT) sink[s].star->nstep*timestep;
      macc = sink[s].menc*max(1.0 - exp(-dt/sink[s].taccrete),0.0);
    }
    else {
      macc = sink[s].menc;
    }

    cout << "Accreting sink " << s << endl;
    cout << "smooth_accretion : " << smooth_accretion << endl;
    cout << "Mass enclosed " << sink[s].menc << "    " << sink[s].mmax << endl;
    cout << "Timestep : " << timestep << "    " << dt << "    efrac : " << efrac <<endl;
    cout << "tvisc : " << sink[s].tvisc << "    trad : " << sink[s].trad 
	 << "   trot : " << sink[s].trot << endl;
    cout << "tacc : " << sink[s].taccrete << "    macc : " << macc << endl;
    //exit(0);


    // Now accrete particles to sink
    // ------------------------------------------------------------------------
    macc_temp = macc;
    for (k=0; k<ndim; k++) rold[k] = sink[s].star->r[k];
    for (k=0; k<ndim; k++) vold[k] = sink[s].star->v[k];
    mold = sink[s].star->m;

    for (k=0; k<ndim; k++) sink[s].star->r[k] *= sink[s].star->m;
    for (k=0; k<ndim; k++) sink[s].star->v[k] *= sink[s].star->m;
    for (k=0; k<ndim; k++) sink[s].star->a[k] *= sink[s].star->m;

    // Loop over all neighbouring particles
    // ------------------------------------------------------------------------
    for (j=0; j<Nneib; j++) {
      i = ilist[j]; //ilist[ilist2[j]];
      mtemp = min(sph->sphdata[i].m,macc_temp);
      dt = sph->sphdata[i].dt;

      // Special conditions for total particle accretion
      if (smooth_accretion == 0 || 
	  sph->sphdata[i].m - mtemp < smooth_accrete_frac*sph->mmean ||
	  dt < smooth_accrete_dt*sink[s].trot) {
	mtemp = sph->sphdata[i].m;
      }
      macc_temp -= mtemp;

      // Now accrete COM quantities to sink particle
      sink[s].star->m += mtemp;
      for (k=0; k<ndim; k++) sink[s].star->r[k] += mtemp*sph->sphdata[i].r[k];
      for (k=0; k<ndim; k++) sink[s].star->v[k] += mtemp*sph->sphdata[i].v[k];
      for (k=0; k<ndim; k++) sink[s].star->a[k] += mtemp*sph->sphdata[i].a[k];
      sink[s].utot += mtemp*sph->sphdata[i].u;

      // If we've reached/exceeded the mass limit, do not include more ptcls
      if (macc_temp < small_number) break;
    }

    // Normalise COM quantities
    for (k=0; k<ndim; k++) sink[s].star->r[k] /= sink[s].star->m;
    for (k=0; k<ndim; k++) sink[s].star->v[k] /= sink[s].star->m;
    for (k=0; k<ndim; k++) sink[s].star->a[k] /= sink[s].star->m;

    for (k=0; k<ndim; k++) sink[s].star->r0[k] = sink[s].star->r[k];
    for (k=0; k<ndim; k++) sink[s].star->v0[k] = sink[s].star->v[k];
    for (k=0; k<ndim; k++) sink[s].star->a0[k] = sink[s].star->a[k];

    //continue;

    // Calculate angular momentum of old COM around new COM
    for (k=0; k<ndim; k++) dr[k] = rold[k] - sink[s].star->r[k];
    for (k=0; k<ndim; k++) dv[k] = vold[k] - sink[s].star->v[k];
    sink[s].angmom[2] += mold*(dr[0]*dv[1] - dr[1]*dv[0]);
    if (ndim == 3) {
      sink[s].angmom[0] += mold*(dr[1]*dv[2] - dr[2]*dv[1]);
      sink[s].angmom[1] += mold*(dr[2]*dv[0] - dr[0]*dv[2]);
    }

    // Now add angular momentum contribution of individual SPH particles
    // ------------------------------------------------------------------------
    for (j=0; j<Nneib; j++) {
      i = ilist[j]; //ilist[ilist2[j]];
      mtemp = min(sph->sphdata[i].m,macc);
      dt = sph->sphdata[i].dt;

      // Special conditions for total particle accretion
      if (smooth_accretion == 0 || 
	  sph->sphdata[i].m - mtemp < smooth_accrete_frac*sph->mmean ||
	  dt < smooth_accrete_dt*sink[s].trot) {
	mtemp = sph->sphdata[i].m;
	deadlist[Ndead++] = i;
      }
      else sph->sphdata[i].m -= mtemp;
      macc -= mtemp;

      // Calculate angular momentum of old COM around new COM
      for (k=0; k<ndim; k++) dr[k] = sph->sphdata[i].r[k] - sink[s].star->r[k];
      for (k=0; k<ndim; k++) dv[k] = sph->sphdata[i].v[k] - sink[s].star->v[k];
      sink[s].angmom[2] += mtemp*(dr[0]*dv[1] - dr[1]*dv[0]);
      if (ndim == 3) {
	sink[s].angmom[0] += mtemp*(dr[1]*dv[2] - dr[2]*dv[1]);
	sink[s].angmom[1] += mtemp*(dr[2]*dv[0] - dr[0]*dv[2]);
      }

      // If we've reached/exceeded the mass limit, do not include more ptcls
      if (macc < small_number) break;

    }
    
    // Calculate internal sink timestep here
    asqd = DotProduct(sink[s].star->a,sink[s].star->a,ndim);
    sink[s].star->dt_internal = 0.4*sqrt(sink[s].radius/(sqrt(asqd) + small_number));

    //}
  // ==========================================================================


  // If particles have been accreted, need to reorder arrays
  if (Ndead > 0) cout << "Accreting dead particles : " << Ndead << endl;
  if (Ndead > 0) sph->DeleteParticles(Ndead,deadlist);
  if (Ndead > 0) cout << "Ngas : " << sph->Nsph << "    Nstar : " 
		      << nbody->Nstar << "   Nsink : " << Nsink << endl;

  // Free memory
  delete[] rsqdlist;
  delete[] ilist2;
  delete[] ilist;
  delete[] deadlist;
  delete[] accrete_list;

  mcom = 0.0;
  for (k=0; k<ndim; k++) rcom[k] = 0.0;
  for (k=0; k<ndim; k++) vcom[k] = 0.0;
  for (k=0; k<ndim; k++) acom[k] = 0.0;
  for (i=0; i<sph->Nsph; i++) {
    mcom += sph->sphdata[i].m;
    for (k=0; k<ndim; k++) rcom[k] += sph->sphdata[i].m*sph->sphdata[i].r[k];
    for (k=0; k<ndim; k++) vcom[k] += sph->sphdata[i].m*sph->sphdata[i].v[k];
    for (k=0; k<ndim; k++) acom[k] += sph->sphdata[i].m*sph->sphdata[i].a[k];
  }
  for (s=0; s<Nsink; s++) {
    mcom += sink[s].star->m;
    for (k=0; k<ndim; k++) rcom[k] += sink[s].star->m*sink[s].star->r[k];
    for (k=0; k<ndim; k++) vcom[k] += sink[s].star->m*sink[s].star->v[k];
    for (k=0; k<ndim; k++) acom[k] += sink[s].star->m*sink[s].star->a[k];
  }
  for (k=0; k<ndim; k++) rcom[k] /= mcom;
  for (k=0; k<ndim; k++) vcom[k] /= mcom;
  for (k=0; k<ndim; k++) acom[k] /= mcom;
  cout << "mcom (after) : " << mcom << endl;
  cout << "rcom (after) : " << rcom[0] << "   " << rcom[1] << "    " << rcom[2] << endl;
  cout << "vcom (after) : " << vcom[0] << "   " << vcom[1] << "    " << vcom[2] << endl;
  cout << "acom (after) : " << acom[0] << "   " << acom[1] << "    " << acom[2] << endl;

  return;
}



// Create template class instances of the main SphSimulation object for
// each dimension used (1, 2 and 3)
template class Sinks<1>;
template class Sinks<2>;
template class Sinks<3>;
