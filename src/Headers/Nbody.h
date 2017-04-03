//=================================================================================================
//  Nbody.h
//  Main N-body class
//
//  This file is part of GANDALF :
//  Graphical Astrophysics code for N-body Dynamics And Lagrangian Fluids
//  https://github.com/gandalfcode/gandalf
//  Contact : gandalfcode@gmail.com
//
//  Copyright (C) 2013  D. A. Hubber, G. Rosotti
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
//=================================================================================================


#ifndef _NBODY_H_
#define _NBODY_H_


#include <string>
#include "Precision.h"
#include "CodeTiming.h"
#include "Constants.h"
#include "DomainBox.h"
#include "Ewald.h"
#include "ExternalPotential.h"
#include "Hydrodynamics.h"
#include "Parameters.h"
#include "SmoothingKernel.h"
#include "NbodyParticle.h"
#include "StarParticle.h"
#include "SystemParticle.h"
#include "Particle.h"
#include "SimUnits.h"
using namespace std;

template <int ndim>
class Hydrodynamics;


//=================================================================================================
//  Class Nbody
/// \brief   Main N-body class.
/// \details Main N-body class for computing forces between stars/systems
///          for N-body dynamics.
/// \author  D. A. Hubber
/// \date    15/04/2013
//=================================================================================================
template <int ndim>
class Nbody
{
protected:

#if defined _OPENMP
  static const int maxNbodyPerThread = 16;     ///< Max. no. of N-body particles per OpenMP thread
                                               ///< if conditional OpenMP is employed.
  const int maxNbodyOpenMp;                    ///< Max. total of N-body particles for OpenMP
                                               ///< if conditional OpenMP is employed.
#endif


public:

  // Constructor and destructor functions
  //-----------------------------------------------------------------------------------------------
  Nbody(int, int, int, DOUBLE, string, int);


  // N-body array memory allocation functions
  //-----------------------------------------------------------------------------------------------
  void AllocateMemory(int);
  void DeallocateMemory(void);
  void LoadStellarPropertiesTable(SimUnits *);
  void UpdateStellarProperties(void);


  // N-body gravitational acceleration routines
  //-----------------------------------------------------------------------------------------------
  void CheckBoundaries(int, int, FLOAT, FLOAT, DomainBox<ndim> &, NbodyParticle<ndim> **);
  virtual void CalculateDirectGravForces(int, NbodyParticle<ndim> **,
                                         DomainBox<ndim> &, Ewald<ndim> *);


  // Other functions
  //-----------------------------------------------------------------------------------------------
  virtual void AdvanceParticles(int, int, FLOAT, FLOAT, NbodyParticle<ndim> **) = 0;
  virtual void CalculateAllStartupQuantities(int, NbodyParticle<ndim> **,
                                             DomainBox<ndim> &, Ewald<ndim> *) = 0;
  virtual void CalculateDirectSmoothedGravForces(int, NbodyParticle<ndim> **,
                                                 DomainBox<ndim> &, Ewald<ndim> *) = 0;
  virtual void CalculateDirectHydroForces(NbodyParticle<ndim> *, int, int, int *, int *,
                                          Hydrodynamics<ndim> *, DomainBox<ndim> &, Ewald<ndim> *) = 0;
  virtual void CalculatePerturberForces(int, int, NbodyParticle<ndim> **, NbodyParticle<ndim> *,
                                        DomainBox<ndim> &, Ewald<ndim> *, FLOAT *, FLOAT *);
  virtual void PerturberCorrectionTerms(int, int, FLOAT, FLOAT, NbodyParticle<ndim> **) = 0;
  virtual void CorrectionTerms(int, int, FLOAT, FLOAT, NbodyParticle<ndim> **) = 0;
  virtual void UpdateChildStars(SystemParticle<ndim> *);
  virtual void EndTimestep(int, int, FLOAT, FLOAT, NbodyParticle<ndim> **) = 0;
  virtual DOUBLE Timestep(NbodyParticle<ndim> *) = 0;
  virtual void IntegrateInternalMotion(SystemParticle<ndim>* system, const int, const FLOAT,
                                       const FLOAT, DomainBox<ndim> &, Ewald<ndim> *);


  // N-body counters and main data arrays
  //----------------------------------------------------------------------------------------------
  bool allocated;                       ///< Is N-body memory allocated
  int Nnbody;                           ///< No. of N-body particles
  int Nnbodymax;                        ///< Max. no. of N-body particles
  int Nstar;                            ///< No. of star particles
  int Nstarmax;                         ///< Max. no. of star particles
  int Nsystem;                          ///< No. of system particles
  int Nsystemmax;                       ///< No. of system particles
  int reset_tree;                       ///< Reset all star properties for tree

  const int nbody_softening;            ///< Use softened-gravity for stars?
  const int perturbers;                 ///< Use perturbers or not
  const int sub_systems;                ///< Create sub-systems?
  const int Npec;                       ///< No. of iterations if using time-symmetric integrator
  const DOUBLE nbody_mult;              ///< N-body timestep multiplier
  static const int vdim=ndim;           ///< Local copy of vdim
  //static const DOUBLE invndim=1.0/ndim; ///< Copy of 1/ndim
  static const FLOAT invndim;

  NbodyParticle<ndim> **nbodydata;      ///< Generic N-body array of ptrs
  StarParticle<ndim> *stardata;         ///< Main star particle data array
  SystemParticle<ndim> *system;         ///< Main system particle array

  CodeTiming *timing;                   ///< Pointer to code timing object
  ExternalPotential<ndim> *extpot;      ///< Pointer to external potential object
  SmoothingKernel<ndim> *kernp;         ///< Pointer to chosen kernel object
  TabulatedKernel<ndim> kerntab;        ///< Tabulated version of chosen kernel


  // Data structures and array for storing stellar property data
  //-----------------------------------------------------------------------------------------------
  struct StellarTableElement {
    FLOAT mass;
    FLOAT luminosity;
    FLOAT NLyC;
    FLOAT Teff;
    FLOAT mdot;
    FLOAT vwind;
  };
  int Nstellartable;
  struct StellarTableElement *stellartable;

};


// Declare invndim constant here (prevents warnings with some compilers)
template <int ndim>
const FLOAT Nbody<ndim>::invndim = 1.0/ndim;



//=================================================================================================
//  Class NbodyLeapfrogKDK
/// Class definition for N-body Leapfrog kick-drift-kick integration scheme.
//=================================================================================================
#if !defined(SWIG)
template <int ndim, template<int> class kernelclass>
class NbodyLeapfrogKDK : public Nbody<ndim>
{
public:
  using Nbody<ndim>::allocated;
  using Nbody<ndim>::invndim;
  using Nbody<ndim>::nbody_mult;
  using Nbody<ndim>::Nstar;
  using Nbody<ndim>::Nsystem;
  using Nbody<ndim>::stardata;
  using Nbody<ndim>::system;
  using Nbody<ndim>::timing;

  NbodyLeapfrogKDK(int, int, int, DOUBLE, string);
  ~NbodyLeapfrogKDK();

  void AdvanceParticles(int, int, FLOAT, FLOAT, NbodyParticle<ndim> **);
  void CalculateAllStartupQuantities(int, NbodyParticle<ndim> **, DomainBox<ndim> &, Ewald<ndim> *) {};
  void CalculateDirectSmoothedGravForces(int, NbodyParticle<ndim> **, DomainBox<ndim> &, Ewald<ndim> *);
  void CalculateDirectHydroForces(NbodyParticle<ndim> *, int, int, int *, int *,
                                 Hydrodynamics<ndim> *, DomainBox<ndim> &, Ewald<ndim> *);
  void CorrectionTerms(int, int, FLOAT, FLOAT, NbodyParticle<ndim> **);
  void EndTimestep(int, int, FLOAT, FLOAT, NbodyParticle<ndim> **);
  void PerturberCorrectionTerms(int, int, FLOAT, FLOAT, NbodyParticle<ndim> **);
  DOUBLE Timestep(NbodyParticle<ndim> *);

  static const int vdim=ndim;             ///< Local copy of vdim
  //static const FLOAT invndim=1.0/ndim;   ///< Copy of 1/ndim
  kernelclass<ndim> kern;                 ///< SPH kernel

};



//=================================================================================================
//  Class NbodyLeapfrogDKD
/// Class definition for N-body Leapfrog drift-kick-drift integration scheme.
//=================================================================================================
template <int ndim, template<int> class kernelclass>
class NbodyLeapfrogDKD : public Nbody<ndim>
{
public:
  using Nbody<ndim>::allocated;
  using Nbody<ndim>::invndim;
  using Nbody<ndim>::nbody_mult;
  using Nbody<ndim>::Nstar;
  using Nbody<ndim>::Nsystem;
  using Nbody<ndim>::stardata;
  using Nbody<ndim>::system;
  using Nbody<ndim>::timing;

  NbodyLeapfrogDKD(int, int, int, DOUBLE, string);
  ~NbodyLeapfrogDKD();

  void AdvanceParticles(int, int, FLOAT, FLOAT, NbodyParticle<ndim> **);
  void CalculateAllStartupQuantities(int, NbodyParticle<ndim> **, DomainBox<ndim> &, Ewald<ndim> *) {};
  void CalculateDirectSmoothedGravForces(int, NbodyParticle<ndim> **, DomainBox<ndim> &, Ewald<ndim> *);
  void CalculateDirectHydroForces(NbodyParticle<ndim> *, int, int, int *, int *,
                                  Hydrodynamics<ndim> *, DomainBox<ndim> &, Ewald<ndim> *);
  void CorrectionTerms(int, int, FLOAT, FLOAT, NbodyParticle<ndim> **) {};
  void EndTimestep(int, int, FLOAT, FLOAT, NbodyParticle<ndim> **);
  void PerturberCorrectionTerms(int, int, FLOAT, FLOAT, NbodyParticle<ndim> **) {};
  DOUBLE Timestep(NbodyParticle<ndim> *);

  static const int vdim=ndim;             ///< Local copy of vdim
  //static const FLOAT invndim=1.0/ndim;   ///< Copy of 1/ndim
  kernelclass<ndim> kern;                 ///< SPH kernel

};



//=================================================================================================
//  Class NbodyHermite4
/// Class definition for N-body 4th-order Hermite integration scheme.
//=================================================================================================
template <int ndim, template<int> class kernelclass>
class NbodyHermite4 : public Nbody<ndim>
{
public:

  using Nbody<ndim>::allocated;
  using Nbody<ndim>::invndim;
  using Nbody<ndim>::nbody_mult;
  using Nbody<ndim>::Nstar;
  using Nbody<ndim>::Nsystem;
  using Nbody<ndim>::stardata;
  using Nbody<ndim>::system;
  using Nbody<ndim>::perturbers;
  using Nbody<ndim>::timing;

  NbodyHermite4(int, int, int, DOUBLE, string, int Npec=1);
  ~NbodyHermite4();

  void AdvanceParticles(int, int, FLOAT, FLOAT, NbodyParticle<ndim> **);
  void CalculateAllStartupQuantities(int, NbodyParticle<ndim> **, DomainBox<ndim> &, Ewald<ndim> *);
  void CalculateDirectSmoothedGravForces(int, NbodyParticle<ndim> **, DomainBox<ndim> &, Ewald<ndim> *);
  void CalculateDirectHydroForces(NbodyParticle<ndim> *, int, int, int *, int *,
                                  Hydrodynamics<ndim> *, DomainBox<ndim> &, Ewald<ndim> *);
  void CorrectionTerms(int, int, FLOAT, FLOAT, NbodyParticle<ndim> **);
  void EndTimestep(int, int, FLOAT, FLOAT, NbodyParticle<ndim> **);
  void PerturberCorrectionTerms(int, int, FLOAT, FLOAT, NbodyParticle<ndim> **);
  DOUBLE Timestep(NbodyParticle<ndim> *);

  static const int vdim=ndim;             ///< Local copy of vdim
  //static const FLOAT invndim=1.0/ndim;   ///< Copy of 1/ndim
  kernelclass<ndim> kern;                 ///< SPH kernel

};



//=================================================================================================
//  Class NbodyHermite4TS
/// Class definition for N-body 4th-order Hermite scheme using time-symmetric iteration.
//=================================================================================================
template <int ndim, template<int> class kernelclass>
class NbodyHermite4TS : public NbodyHermite4<ndim, kernelclass>
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
  using Nbody<ndim>::timing;

  NbodyHermite4TS(int, int, int, DOUBLE, string, int);
  ~NbodyHermite4TS();

  void CorrectionTerms(int, int, FLOAT, FLOAT, NbodyParticle<ndim> **);
  //void IntegrateInternalMotion(SystemParticle<ndim>* system, int, FLOAT, FLOAT);

};



//=================================================================================================
//  Class NbodyHermite6TS
/// Class definition for N-body 6th-order Hermite scheme using time-symmetric iteration.
//=================================================================================================
template <int ndim, template<int> class kernelclass>
class NbodyHermite6TS : public Nbody<ndim>
{
public:

  using Nbody<ndim>::allocated;
  using Nbody<ndim>::invndim;
  using Nbody<ndim>::nbody_mult;
  using Nbody<ndim>::Nstar;
  using Nbody<ndim>::Nsystem;
  using Nbody<ndim>::stardata;
  using Nbody<ndim>::system;
  using Nbody<ndim>::perturbers;
  using Nbody<ndim>::timing;

  NbodyHermite6TS(int, int, int, DOUBLE, string, int);
  ~NbodyHermite6TS();

  void CalculateDirectGravForces(int, NbodyParticle<ndim> **, DomainBox<ndim> &, Ewald<ndim> *);
  void AdvanceParticles(int, int, FLOAT, FLOAT, NbodyParticle<ndim> **);
  void CalculateAllStartupQuantities(int, NbodyParticle<ndim> **, DomainBox<ndim> &, Ewald<ndim> *);
  void CalculateDirectSmoothedGravForces(int, NbodyParticle<ndim> **, DomainBox<ndim> &, Ewald<ndim> *);
  void CalculateDirectHydroForces(NbodyParticle<ndim> *, int, int, int *, int *,
                                  Hydrodynamics<ndim> *, DomainBox<ndim> &, Ewald<ndim> *);
  void CorrectionTerms(int, int, FLOAT, FLOAT, NbodyParticle<ndim> **);
  void EndTimestep(int, int, FLOAT, FLOAT, NbodyParticle<ndim> **);
  void PerturberCorrectionTerms(int, int, FLOAT, FLOAT, NbodyParticle<ndim> **);
  DOUBLE Timestep(NbodyParticle<ndim> *);

  static const int vdim=ndim;             ///< Local copy of vdim
  //static const FLOAT invndim=1.0/ndim;   ///< Copy of 1/ndim
  kernelclass<ndim> kern;                 ///< SPH kernel

};
#endif
#endif
