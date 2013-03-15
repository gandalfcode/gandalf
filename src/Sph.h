// ============================================================================
// Sph.h
// Contains main parent virtual class plus child classes for various SPH 
// algorithms that are implemented.
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


#if !defined(SWIG) && defined(FIXED_DIMENSIONS)
  Sph(int ndimaux, int vdimaux, int bdimaux, int hydro_forces_aux,
    int self_gravity_aux, FLOAT alpha_visc_aux, FLOAT beta_visc_aux,
    FLOAT h_fac_aux, FLOAT h_converge_aux, string avisc_aux,
    string acond_aux, string gas_eos_aux, string KernelName):
	hydro_forces(hydro_forces_aux),
	self_gravity(self_gravity_aux),
	alpha_visc(alpha_visc_aux),
	beta_visc(beta_visc_aux),
	h_fac(h_fac_aux),
	h_converge(h_converge_aux),
	avisc(avisc_aux),
	acond(acond_aux),
	gas_eos(gas_eos_aux),
    kerntab(TabulatedKernel(ndimaux, KernelName))
      {};
#elif !defined(SWIG) && !defined(FIXED_DIMENSIONS)
  Sph(int ndimaux, int vdimaux, int bdimaux, int hydro_forces_aux,
    int self_gravity_aux, FLOAT alpha_visc_aux, FLOAT beta_visc_aux,
    FLOAT h_fac_aux, FLOAT h_converge_aux, string avisc_aux,
    string acond_aux, string gas_eos_aux, string KernelName):
	hydro_forces(hydro_forces_aux),
	self_gravity(self_gravity_aux),
	alpha_visc(alpha_visc_aux),
	beta_visc(beta_visc_aux),
	h_fac(h_fac_aux),
	h_converge(h_converge_aux),
	avisc(avisc_aux),
	acond(acond_aux),
	gas_eos(gas_eos_aux),
    ndim(ndimaux), 
    vdim(vdimaux), 
    bdim(bdimaux), 
    invndim(1.0/(FLOAT)ndimaux),
    kerntab(TabulatedKernel(ndimaux, KernelName))
      {};
#endif

  // SPH functions for computing SPH sums with neighbouring particles 
  // (fully coded in each separate SPH implementation, and not in Sph.cpp)
  // --------------------------------------------------------------------------
  virtual int ComputeH(int, int, FLOAT *, FLOAT *, SphParticle &) = 0;
  virtual void ComputeSphNeibForces(int, int, int *, 
				    FLOAT *, FLOAT *, FLOAT *, 
				    SphParticle &, SphParticle *) = 0;
  virtual void ComputeDirectGravForces(int, int, int *,
				       SphParticle &, SphParticle *) = 0;

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
  int Nghostmax;                        // Max. allowed no. of ghost particles

  const FLOAT alpha_visc;               // alpha artificial viscosity parameter
  const FLOAT beta_visc;                // beta artificial viscosity parameter
  const FLOAT h_fac;                    // Smoothing length-density factor
  const FLOAT h_converge;               // h-rho iteration tolerance
  const string avisc;                   // Artificial viscosity option
  const string acond;                   // Artificial conductivity option
  const string gas_eos;                 // Gas EOS option
  const int hydro_forces;               // Compute hydro forces?
  const int self_gravity;               // Compute gravitational forces?

  struct SphParticle *sphdata;          // Main SPH particle data array
  EOS *eos;                             // Equation-of-state
  SphKernel *kernp;                     // Pointer to chosen kernel object
  TabulatedKernel kerntab;        // Tabulated version of the chosen kernel

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

  GradhSph(int, int, int, int, int, FLOAT, FLOAT, FLOAT, FLOAT,
		  string, string, string, string, string);
  ~GradhSph();

  int ComputeH(int, int, FLOAT *, FLOAT *, SphParticle &);
  void ComputeSphNeibForces(int, int, int *, FLOAT *, FLOAT *, 
			    FLOAT *, SphParticle &, SphParticle *);
  void ComputeDirectGravForces(int, int, int *, SphParticle &, SphParticle *);

};


#endif
