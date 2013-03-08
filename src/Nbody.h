// ============================================================================
// Nbody.h
// ..
// ============================================================================


#ifndef _NBODY_H_
#define _NBODY_H_


#include <string>
#include "Constants.h"
#include "Dimensions.h"
#include "Parameters.h"
using namespace std;


// ============================================================================
// Class Nbody
// Main N-body class.
// ============================================================================
class Nbody
{
 public:

#if !defined(SWIG) && !defined(FIXED_DIMENSIONS)
 Nbody(int ndimaux, int vdimaux):
  ndim(ndimaux),
    vdim(vdimaux)
    {};
#endif

  // N-body array memory allocation functions
  // --------------------------------------------------------------------------
  void AllocateMemory(int);
  void DeallocateMemory(void);

  // Other functions
  // --------------------------------------------------------------------------
  void CalculateDirectGravForces(int, StarParticle &, StarParticle *);

  // N-body counters and main data arrays
  // --------------------------------------------------------------------------
  bool allocated;                        // Is N-body memory allocated
  int Nstar;                             // No. of star particles
  int Nsystem;                           // No. of system particles

  struct StarParticle *stardata;         // Main star particle data array
  struct StarParticle *system;           // Main system particle array

};
