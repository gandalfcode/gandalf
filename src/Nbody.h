//=============================================================================
//  Nbody.h
//  Main N-body class
//
//  This file is part of GANDALF :
//  Graphical Astrophysics code for N-body Dynamics and Lagrangian Fluids
//  https://github.com/gandalfcode/gandalf
//  Contact : gandalfcode@gmail.com
//
//  Copyright (C) 2013  D. A. Hubber, G Rosotti
//
//  GANDALF is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 2 of the License, or
//  (at your option) any later version.
//
//  GANDALF is distributed in the hope that it will be useful, but
//  WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  General Public License (http://www.gnu.org/licenses) for more details.
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
	string KernelName, int);


  // N-body array memory allocation functions
  // --------------------------------------------------------------------------
  void AllocateMemory(int);
  void DeallocateMemory(void);


  // Other functions
  // --------------------------------------------------------------------------
  virtual void AdvanceParticles(int, int, NbodyParticle<ndim> **,DOUBLE) = 0;
  virtual void CalculateAllStartupQuantities(int, NbodyParticle<ndim> **) = 0;
  virtual void CalculateDirectGravForces(int, NbodyParticle<ndim> **) = 0;
  virtual void CalculateDirectSPHForces(int, int, SphParticle<ndim> *,
                                        NbodyParticle<ndim> **) = 0;
  virtual void CalculatePerturberForces(int, int, NbodyParticle<ndim> **,
                                        NbodyParticle<ndim> *,
                                        DOUBLE *, DOUBLE *) = 0;
  virtual void PerturberCorrectionTerms(int, int, NbodyParticle<ndim> **, DOUBLE) = 0;
  virtual void CorrectionTerms(int, int, NbodyParticle<ndim> **,DOUBLE) = 0;
  virtual void CorrectPerturbedChildStars(SystemParticle<ndim>* system,
                                          int, DOUBLE, DOUBLE tend) = 0;
  virtual void EndTimestep(int, int, NbodyParticle<ndim> **) = 0;
  virtual DOUBLE Timestep(NbodyParticle<ndim> *) = 0;
  virtual void IntegrateInternalMotion(SystemParticle<ndim>* system,
                                       int, DOUBLE, DOUBLE tend);


  // N-body counters and main data arrays
  // --------------------------------------------------------------------------
  bool allocated;                       ///< Is N-body memory allocated
  int Nnbody;                           ///< No. of N-body particles
  int Nnbodymax;                        ///< Max. no. of N-body particles
  int Nstar;                            ///< No. of star particles
  int Nstarmax;                         ///< Max. no. of star particles
  int Nsystem;                          ///< No. of system particles
  int Nsystemmax;                       ///< No. of system particles
  int reset_tree;                       ///< Reset all star properties for tree
  int perturbers;                       ///< Use perturbers or not

  const int nbody_softening;            ///< Use softened-gravity for stars?
  const int sub_systems;                ///< Create sub-systems?
  const int Npec;                       ///< Number of iterations (if using  
                                        ///< a time-symmetric integrator)
  const DOUBLE nbody_mult;              ///< N-body timestep multiplier
 
  static const int vdim=ndim;           ///< Local copy of vdim
  static const FLOAT invndim=1./ndim;   ///< Copy of 1/ndim

  SphKernel<ndim> *kernp;               ///< Pointer to chosen kernel object
  TabulatedKernel<ndim> kerntab;        ///< Tabulated version of chosen kernel
  struct NbodyParticle<ndim> **nbodydata; ///< Generic N-body array of ptrs
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
public:
  using Nbody<ndim>::allocated;
  using Nbody<ndim>::nbody_mult;
  using Nbody<ndim>::Nstar;
  using Nbody<ndim>::Nsystem;
  using Nbody<ndim>::stardata;
  using Nbody<ndim>::system;

  NbodyLeapfrogKDK(int, int, DOUBLE, string);
  ~NbodyLeapfrogKDK();

  void AdvanceParticles(int, int, NbodyParticle<ndim> **,DOUBLE);
  void CalculateAllStartupQuantities(int, NbodyParticle<ndim> **);
  void CalculateDirectGravForces(int, NbodyParticle<ndim> **);
  void CalculateDirectSPHForces(int, int, SphParticle<ndim> *,
				NbodyParticle<ndim> **);
  void CalculatePerturberForces(int, int, NbodyParticle<ndim> **,
				NbodyParticle<ndim> *, DOUBLE *, DOUBLE *);
  void CorrectionTerms(int, int, NbodyParticle<ndim> **, DOUBLE);
  void CorrectPerturbedChildStars(SystemParticle<ndim>* system,
                                  int, DOUBLE, DOUBLE tend);
  void EndTimestep(int, int, NbodyParticle<ndim> **);
  void IntegrateInternalMotion(SystemParticle<ndim>* system,
                               int, DOUBLE, DOUBLE tend);
  void PerturberCorrectionTerms(int, int, NbodyParticle<ndim> **, DOUBLE);
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
public:

  using Nbody<ndim>::allocated;
  using Nbody<ndim>::nbody_mult;
  using Nbody<ndim>::Nstar;
  using Nbody<ndim>::Nsystem;
  using Nbody<ndim>::stardata;
  using Nbody<ndim>::system;

  NbodyHermite4(int, int, DOUBLE, string, int Npec=1);
  ~NbodyHermite4();

  void AdvanceParticles(int, int, NbodyParticle<ndim> **,DOUBLE);
  void CalculateAllStartupQuantities(int, NbodyParticle<ndim> **);
  void CalculateDirectGravForces(int, NbodyParticle<ndim> **);
  void CalculateDirectSPHForces(int, int, SphParticle<ndim> *,
				NbodyParticle<ndim> **);
  void CalculatePerturberForces(int, int, NbodyParticle<ndim> **,
				NbodyParticle<ndim> *, DOUBLE *, DOUBLE *);
  void CorrectionTerms(int, int, NbodyParticle<ndim> **,DOUBLE);
  void EndTimestep(int, int, NbodyParticle<ndim> **);
  void IntegrateInternalMotion(SystemParticle<ndim>* system,
                               int, DOUBLE, DOUBLE tend);
  void CorrectPerturbedChildStars(SystemParticle<ndim>* system,
                                  int, DOUBLE, DOUBLE tend);
  void PerturberCorrectionTerms(int, int, NbodyParticle<ndim> **, DOUBLE);
  DOUBLE Timestep(NbodyParticle<ndim> *);

  static const int vdim=ndim;           ///< Local copy of vdim
  static const FLOAT invndim=1./ndim;   ///< Copy of 1/ndim
  kernelclass<ndim> kern;               ///< SPH kernel

};



//=============================================================================
//  Class NbodyHermite4TS
/// Class definition for N-body 4th-order Hermite integration scheme using
/// time-symmetric iteration.
//=============================================================================
template <int ndim, template<int> class kernelclass>
class NbodyHermite4TS: public NbodyHermite4<ndim, kernelclass>
{
public:

  using NbodyHermite4<ndim, kernelclass>::allocated;
  using Nbody<ndim>::nbody_mult;
  using Nbody<ndim>::Npec;
  using Nbody<ndim>::Nstar;
  using Nbody<ndim>::Nsystem;
  using Nbody<ndim>::stardata;
  using Nbody<ndim>::system;
  using Nbody<ndim>::Timestep;
  using Nbody<ndim>::AdvanceParticles;
  using Nbody<ndim>::EndTimestep;
  using Nbody<ndim>::CalculateDirectGravForces;
  using Nbody<ndim>::perturbers;

  NbodyHermite4TS(int, int, DOUBLE, string, int);
  ~NbodyHermite4TS();

  void CorrectionTerms(int, int, NbodyParticle<ndim> **,DOUBLE);
  void IntegrateInternalMotion(SystemParticle<ndim>* system,
                               int, DOUBLE, DOUBLE tend);

};
#endif

#endif
