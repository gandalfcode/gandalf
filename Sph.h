// ============================================================================
// Sph.h
// ============================================================================


#ifndef _SPH_H_
#define _SPH_H_


#include "Constants.h"
#include "Dimensions.h"
#include "SphParticle.h"
#include "SphKernel.h"
#include "Parameters.h"
#include "EOS.h"


// ============================================================================
// CLASS Sph
// ============================================================================
class Sph
{
 public:

  // SPH functions for computing SPH sums with neighbouring particles 
  // (to be coded in each separate SPH implementation)
  // --------------------------------------------------------------------------
  virtual void ComputeH(int,int,int *,Parameters &) = 0;
  virtual void ComputeSphProperties(int,int,int *,Parameters &) = 0;
  virtual void ComputeHydroForces(int,int,int *,Parameters &) = 0;
  virtual void ComputeGravForce(int,int,float *) = 0;
  virtual void RandomBox(void) = 0;


  // SPH array memory allocation functions
  // --------------------------------------------------------------------------
  void AllocateMemory(int);
  void DeallocateMemory(void);
  void SphBoundingBox(float *, float *);
  void InitialSmoothingLengthGuess(void);


#if !defined(FIXED_DIMENSIONS)
  int ndim;
  int vdim;
  int bdim;
#endif


  // SPH particle counters and main particle data array
  // --------------------------------------------------------------------------
  int Nsph;                             // No. of SPH particles in simulation
  int Nsphmax;                          // Max. no. of SPH particles in array
  struct SphParticle *sphdata;          // Main SPH particle data array

  SphKernel *kern;
  EOS *eos;

  double alpha_visc;
  double beta_visc;

};



// ============================================================================
// Class GradhSph
// ============================================================================
class GradhSph: public Sph
{
 public:

  GradhSph(int,int,int);
  ~GradhSph();

  void ComputeH(int,int,int *,Parameters &);
  void ComputeSphProperties(int,int,int *,Parameters &);
  void ComputeHydroForces(int,int,int *, Parameters &);
  void ComputeGravForce(int,int,float *);

  void RandomBox(void);

};


#endif
