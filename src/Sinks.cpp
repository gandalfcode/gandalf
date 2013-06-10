//=============================================================================
//  Sinks.cpp
//  ..
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
void Sinks<ndim>::AllocateMemory(void)
{
  return;
}



//=============================================================================
//  Sinks::DeallocateMemory
/// ..
//=============================================================================
template <int ndim>
void Sinks<ndim>::DeallocateMemory(void)
{
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
  FLOAT rho_max;                    // Maximum density of sink candidates

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
      if (sph->sphdata[i].rho > rho_sink) continue;

      // Only consider SPH particles located at a local potential minimum
      if (!sph->sphdata[i].potmin) continue;

      // Make sure candidate particle is at the end of its current timestep
      if (n%sph[i]->sphdata[i].nstep != 0) continue;

      // If SPH particle neighbours a nearby sink, skip to next particle
      for (s=0; s<Nsink; s++) {
        for (k=0; k<ndim; k++) dr[k] = sph->sphdata[i].r[k] - sink[s].r[k];
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
      CreateNewSinkParticle(isink);
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
void Sinks<ndim>::CreateNewSinkParticle(int isink)
{
  debug2("[Sinks::CreateNewSinkParticle]");

  // If we've reached the maximum number of sinks, then throw exception
  if (Nsink == Nsinkmax) {
    

  }

  return;
}



//=============================================================================
//  Sinks::AcceteMassToSinks
/// ..
//=============================================================================
template <int ndim>
void Sinks<ndim>::AccreteMassToSinks(void)
{
  debug2("[Sinks::AccreteMassToSinks]");

  return;
}



//=============================================================================
//  Sinks::UpdateSinkProperties
/// ..
//=============================================================================
template <int ndim>
void Sinks<ndim>::UpdateSinkProperties(void)
{
  debug2("[Sinks::UpdateSinkProperties]");

  return;
}



//=============================================================================
//  Sinks::UpdateStarProperties
/// ..
//=============================================================================
template <int ndim>
void Sinks<ndim>::UpdateStarProperties(void)
{
  debug2("[Sinks::UpdateStarProperties]");
  return;
}
