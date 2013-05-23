//=============================================================================
//  Nbody.h
//  Main N-body class
//=============================================================================


#ifndef _NBODY_H_
#define _NBODY_H_


#include <string>
#include "Precision.h"
#include "Constants.h"
#include "Parameters.h"
#include "SphKernel.h"
#include "NbodyParticle.h"
#include "StarParticle.h"
#include "SystemParticle.h"
#include "SphParticle.h"
using namespace std;



//=============================================================================
//  Class Nbody
/// \brief   Main N-body class.
/// \details Main N-body class for computing forces between stars/systems 
///          for N-body dynamics.
/// \author  D. A. Hubber
/// \date    15/04/2013
//=============================================================================
template <int ndim>
class Nbody
{
 public:

  // Constructor and destructor functions
  // --------------------------------------------------------------------------
  Nbody(int nbody_softening_aux, int sub_systems_aux, DOUBLE nbody_mult_aux,
	string KernelName);


  // N-body array memory allocation functions
  // --------------------------------------------------------------------------
  void AllocateMemory(int);
  void DeallocateMemory(void);


  // Other functions
  // --------------------------------------------------------------------------
  //void CalculateDirectSoftenedGravForces(void);
  virtual void CalculateDirectGravForces(int,NbodyParticle<ndim> **) = 0;
  virtual void CalculateDirectSPHForces(int,int,SphParticle<ndim> *,
					NbodyParticle<ndim> **) = 0;
  virtual void AdvanceParticles(int,int,NbodyParticle<ndim> **,DOUBLE) = 0;
  virtual void CorrectionTerms(int,int,NbodyParticle<ndim> **,DOUBLE) = 0;
  virtual void EndTimestep(int,int,NbodyParticle<ndim> **) = 0;
  virtual DOUBLE Timestep(NbodyParticle<ndim> *) = 0;


  // N-body counters and main data arrays
  // --------------------------------------------------------------------------
  bool allocated;                       ///< Is N-body memory allocated
  int Nnbody;                           ///< ..
  int Nnbodymax;                        ///< ..
  int Nstar;                            ///< No. of star particles
  int Nstarmax;                         ///< Max. no. of star particles
  int Nsystem;                          ///< No. of system particles
  int Nsystemmax;                       ///< No. of system particles

  const int nbody_softening;            ///< Use softened-gravity for stars?
  const int sub_systems;                ///< Create sub-systems?
  const DOUBLE nbody_mult;              ///< N-body timestep multiplier
  static const int vdim=ndim;           ///< Local copy of vdim
  static const FLOAT invndim=1./ndim;   ///< Copy of 1/ndim

  SphKernel<ndim> *kernp;               ///< Pointer to chosen kernel object
  TabulatedKernel<ndim> kerntab;        ///< Tabulated version of chosen kernel
  struct NbodyParticle<ndim> **nbodydata;  ///< Generic N-body array of ptrs
  struct StarParticle<ndim> *stardata;  ///< Main star particle data array
  struct SystemParticle<ndim> *system;  ///< Main system particle array

};



//=============================================================================
//  Class NbodyLeapfrogKDK
/// Class definition for N-body Leapfrog kick-drift-kick integration scheme.
//=============================================================================
#if !defined(SWIG)
template <int ndim, template<int> class kernelclass>
class NbodyLeapfrogKDK: public Nbody<ndim>
{
  using Nbody<ndim>::allocated;
  using Nbody<ndim>::nbody_mult;
  using Nbody<ndim>::Nstar;
  using Nbody<ndim>::Nsystem;
  using Nbody<ndim>::stardata;
  using Nbody<ndim>::system;

 public:

  NbodyLeapfrogKDK(int, int, DOUBLE, string);
  ~NbodyLeapfrogKDK();

  void CalculateDirectGravForces(int,NbodyParticle<ndim> **);
  void CalculateDirectSPHForces(int,int,SphParticle<ndim> *,
				NbodyParticle<ndim> **);
  void AdvanceParticles(int,int,NbodyParticle<ndim> **,DOUBLE);
  void CorrectionTerms(int,int,NbodyParticle<ndim> **,DOUBLE);
  void EndTimestep(int,int,NbodyParticle<ndim> **);
  DOUBLE Timestep(NbodyParticle<ndim> *);

  static const int vdim=ndim;           ///< Local copy of vdim
  static const FLOAT invndim=1./ndim;   ///< Copy of 1/ndim
  kernelclass<ndim> kern;               ///< SPH kernel

};



//=============================================================================
//  Class NbodyHermite4
/// Class definition for N-body 4th-order Hermite integration scheme.
//=============================================================================
template <int ndim, template<int> class kernelclass>
class NbodyHermite4: public Nbody<ndim>
{
  using Nbody<ndim>::allocated;
  using Nbody<ndim>::nbody_mult;
  using Nbody<ndim>::Nstar;
  using Nbody<ndim>::Nsystem;
  using Nbody<ndim>::stardata;
  using Nbody<ndim>::system;

 public:

  NbodyHermite4(int, int, DOUBLE, string);
  ~NbodyHermite4();

  void CalculateDirectGravForces(int,NbodyParticle<ndim> **);
  void CalculateDirectSPHForces(int,int,SphParticle<ndim> *,
				NbodyParticle<ndim> **);
  void AdvanceParticles(int,int,NbodyParticle<ndim> **,DOUBLE);
  void CorrectionTerms(int,int,NbodyParticle<ndim> **,DOUBLE);
  void EndTimestep(int,int,NbodyParticle<ndim> **);
  DOUBLE Timestep(NbodyParticle<ndim> *);

  static const int vdim=ndim;           ///< Local copy of vdim
  static const FLOAT invndim=1./ndim;   ///< Copy of 1/ndim
  kernelclass<ndim> kern;               ///< SPH kernel

};
#endif

#endif
