// ============================================================================
// Sph.h
// ============================================================================


#ifndef _SPH_H_
#define _SPH_H_


#include <string>
#include "Constants.h"
#include "Dimensions.h"
#include "SphParticle.h"
#include "SphKernel.h"
#include "Parameters.h"
#include "EOS.h"
using namespace std;


// ============================================================================
// Class Sph
// Main parent Sph class.  Different SPH implementations 
// (e.g. grad-h SPH, Saitoh & Makino 2012) are derived from this class.
// Each implementation requires defining its own version of each function 
// (e.g. ComputeH for its own method of computing smoothing lengths).
// ============================================================================
class Sph
{
 public:

#if !defined(FIXED_DIMENSIONS)
  Sph(int ndimaux, int vdimaux, int bdimaux):
    ndim(ndimaux),
    vdim(vdimaux),
    bdim(bdimaux),
    invndim (1.0/(FLOAT)ndim)
    {};
#endif

  // SPH functions for computing SPH sums with neighbouring particles 
  // (fully coded in each separate SPH implementation, and not in Sph.cpp)
  // --------------------------------------------------------------------------
  virtual int ComputeH(int,SphParticle &,int,
		       SphParticle *,FLOAT *,FLOAT *,FLOAT *) = 0;
  virtual void ComputeHydroForces(int,SphParticle &,int,
				  SphParticle *,FLOAT *,FLOAT *,FLOAT *) = 0;
  virtual void ComputeGravForces(int,int,SphParticle *) = 0;
  virtual void ComputeMeanhZeta(int,int,int *) = 0;

  // SPH array memory allocation functions
  // --------------------------------------------------------------------------
  void AllocateMemory(int);
  void DeallocateMemory(void);
  void SphBoundingBox(FLOAT *, FLOAT *, int);
  void InitialSmoothingLengthGuess(void);

  // SPH particle counters and main particle data array
  // --------------------------------------------------------------------------
  bool allocated;                       // Is SPH memory allocated?
  int Nsph;                             // No. of SPH particles in simulation
  int Nghost;                           // No. of ghost SPH particles
  int Ntot;                             // No. of real + ghost particles
  int Nsphmax;                          // Max. no. of SPH particles in array
  int Nghostmax;                        // ..

  FLOAT alpha_visc;                     // alpha artificial viscosity parameter
  FLOAT beta_visc;                      // beta artificial viscosity parameter
  FLOAT h_fac;                          // ..
  FLOAT h_converge;                     // ..
  string avisc;                         // ..
  string acond;                         // ..
  string gas_eos;                       // ..
  int hydro_forces;                     // ..
  int self_gravity;                     // ..

  struct SphParticle *sphdata;          // Main SPH particle data array
  EOS *eos;                             // Equation-of-state
  SphKernel* kernp;

#if !defined(FIXED_DIMENSIONS)
  const int ndim;
  const int vdim;
  const int bdim;
  const FLOAT invndim;
#endif

};



// ============================================================================
// Class GradhSph
// Class definition for conservative 'grad-h' SPH simulations (as derived 
// from the parent Sph class).  Full code for each of these class functions 
// written in 'GradhSph.cpp'.
// ============================================================================
template <class kernelclass>
class GradhSph: public Sph
{
 public:

  kernelclass kern;                      // SPH kernel

  GradhSph(int,int,int);
  ~GradhSph();

  int ComputeH(int,SphParticle &,int,SphParticle *,FLOAT *,FLOAT *,FLOAT *);
  void ComputeHydroForces(int,SphParticle &,int,SphParticle *,FLOAT *,FLOAT *,FLOAT *);
  void ComputeGravForces(int,int,SphParticle *);
  void ComputeMeanhZeta(int,int,int *);

};


#endif
