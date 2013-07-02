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
/// ..
//=============================================================================
template <int ndim>
Sinks<ndim>::Sinks()
{
  allocated_memory = false;
  Nsink = 0;
}



//=============================================================================
//  Sinks::~Sinks()
/// ..
//=============================================================================
template <int ndim>
Sinks<ndim>::~Sinks()
{
}



//=============================================================================
//  Sinks::AllocateMemory
/// ..
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
/// ..
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
/// ..
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
      
      // If density of SPH particle is too low, skip to next particle
      if (sph->sphdata[i].rho < rho_sink) continue;

      // Only consider SPH particles located at a local potential minimum
      if (sph->sphdata[i].gpotmin <= 0.9999999*sph->sphdata[i].gpot) continue;

      // Make sure candidate particle is at the end of its current timestep
      if (n%sph->sphdata[i].nstep != 0) continue;

      // If SPH particle neighbours a nearby sink, skip to next particle
      for (s=0; s<Nsink; s++) {
        for (k=0; k<ndim; k++) 
	  dr[k] = sph->sphdata[i].r[k] - sink[s].star->r[k];
        drsqd = DotProduct(dr,dr,ndim);
        if (drsqd < pow(sink_radius*sph->sphdata[i].h + sink[s].radius,2))
          sink_flag = false;
      }

      cout << "SINK?? : " << i << "    " << sph->sphdata[i].rho << "    " << rho_sink << endl;

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
      exit(0);
    }


  } while (isink != -1);
  // ==========================================================================

  return;
}



//=============================================================================
//  Sinks::CreateNewSinkParticle
/// ..
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

  cout << "CHECKING : " << sink[Nsink].star << "    " << &nbody->stardata[nbody->Nstar] << endl;
  cout << "CHECKING2 : " << sph->sphdata[isink].m << endl;

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
  for (k=0; k<ndim; k++) sink[Nsink].star->adot[k] = 0.0; //sph->sphdata[isink].a[k];
  for (k=0; k<ndim; k++) sink[Nsink].star->r0[k] = sph->sphdata[isink].r0[k];
  for (k=0; k<ndim; k++) sink[Nsink].star->v0[k] = sph->sphdata[isink].v0[k];
  for (k=0; k<ndim; k++) sink[Nsink].star->a0[k] = sph->sphdata[isink].a0[k];
  for (k=0; k<ndim; k++) sink[Nsink].star->adot0[k] = 0.0;//sph->sphdata[isink].a[k];
  for (k=0; k<3; k++) sink[Nsink].angmom[k] = 0.0;

  // Increment star and sink counters
  nbody->Nstar++;
  Nsink++;

  return;
}



//=============================================================================
//  Sinks::AcceteMassToSinks
/// ..
//=============================================================================
template <int ndim>
void Sinks<ndim>::AccreteMassToSinks
(Sph<ndim> *sph,                    ///< [inout] Object containing SPH ptcls
 Nbody<ndim> *nbody)                ///< [inout] Object containing star ptcls
{
  debug2("[Sinks::AccreteMassToSinks]");

  return;
}



// Create template class instances of the main SphSimulation object for
// each dimension used (1, 2 and 3)
template class Sinks<1>;
template class Sinks<2>;
template class Sinks<3>;
