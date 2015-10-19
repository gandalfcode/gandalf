//=================================================================================================
//  MeshlessFV.h
//  Contains main parent virtual class plus child classes for various MeshlessFV
//  algorithms that are implemented.
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


#ifndef _MESHLESS_FV_H_
#define _MESHLESS_FV_H_


#include <assert.h>
#include <string>
#include "Constants.h"
#include "DomainBox.h"
#include "EOS.h"
#include "ExternalPotential.h"
#include "FV.h"
#include "Hydrodynamics.h"
#include "NbodyParticle.h"
#include "Nbody.h"
#include "Parameters.h"
#include "Particle.h"
#include "Precision.h"
#include "RiemannSolver.h"
#include "SlopeLimiter.h"
#include "SmoothingKernel.h"
#if defined _OPENMP
#include "omp.h"
#endif
using namespace std;


//=================================================================================================
//  Class MeshlessFV
/// \brief   Main parent MeshlessFV class.
/// \details
/// \author  D. A. Hubber, J. Ngoumou
/// \date    19/02/2015
//=================================================================================================
//template <int ndim>
template <int ndim>
class MeshlessFV : public FV<ndim>
{
public:

  using Hydrodynamics<ndim>::allocated;
  using Hydrodynamics<ndim>::eos;
  using Hydrodynamics<ndim>::h_fac;
  using Hydrodynamics<ndim>::hydrodata_unsafe;
  using Hydrodynamics<ndim>::hydro_forces;
  using Hydrodynamics<ndim>::invndim;
  using Hydrodynamics<ndim>::iorder;
  using Hydrodynamics<ndim>::kernfac;
  using Hydrodynamics<ndim>::kernfacsqd;
  using Hydrodynamics<ndim>::kernp;
  using Hydrodynamics<ndim>::kernrange;
  using Hydrodynamics<ndim>::mmean;
  using Hydrodynamics<ndim>::Ngather;
  using Hydrodynamics<ndim>::Nghost;
  using Hydrodynamics<ndim>::Nghostmax;
  using Hydrodynamics<ndim>::NImportedParticles;
  using Hydrodynamics<ndim>::Nhydro;
  using Hydrodynamics<ndim>::Nhydromax;
  using Hydrodynamics<ndim>::Nmpighost;
  using Hydrodynamics<ndim>::NPeriodicGhost;
  using Hydrodynamics<ndim>::Ntot;
  using Hydrodynamics<ndim>::riemann;
  using Hydrodynamics<ndim>::self_gravity;
  using Hydrodynamics<ndim>::size_hydro_part;
  using FV<ndim>::gamma_eos;
  using FV<ndim>::gammam1;

  static const int nvar = ndim + 2;
  static const int ivx = 0;
  static const int ivy = 1;
  static const int ivz = 2;
  static const int irho = ndim;
  static const int ietot = ndim + 1;
  static const int ipress = ndim + 1;


  // Constructor
  //-----------------------------------------------------------------------------------------------
  MeshlessFV(int hydro_forces_aux, int self_gravity_aux, FLOAT _accel_mult,
             FLOAT _courant_mult, FLOAT h_fac_aux, FLOAT h_converge_aux,
             FLOAT gamma_aux, string gas_eos_aux, string KernelName, int size_mfv_part);
  ~MeshlessFV();


  virtual void AllocateMemory(int);
  virtual void DeallocateMemory(void);
  virtual void DeleteDeadParticles(void);
  virtual void ReorderParticles(void);


  // MeshlessFV functions for computing MeshlessFV sums with neighbouring particles
  // (fully coded in each separate MeshlessFV implementation, and not in MeshlessFV.cpp)
  //-----------------------------------------------------------------------------------------------
  virtual int ComputeH(const int, const int, const FLOAT, FLOAT *, FLOAT *, FLOAT *, FLOAT *,
                       MeshlessFVParticle<ndim> &, Nbody<ndim> *) = 0;
  virtual void ComputePsiFactors(const int, const int, int *, FLOAT *, FLOAT *, FLOAT *,
                                 MeshlessFVParticle<ndim> &, MeshlessFVParticle<ndim> *) = 0;
  virtual void ComputeGradients(const int, const int, int *, FLOAT *, FLOAT *, FLOAT *,
                                MeshlessFVParticle<ndim> &, MeshlessFVParticle<ndim> *) = 0;
  virtual void ComputeGodunovFlux(const int, const int, const FLOAT, int *, FLOAT *, FLOAT *, FLOAT *,
                                  MeshlessFVParticle<ndim> &, MeshlessFVParticle<ndim> *) = 0;
  virtual void ComputeSmoothedGravForces(const int, const int, int *, MeshlessFVParticle<ndim> &,
                                         MeshlessFVParticle<ndim> *) = 0;
  virtual void ComputeDirectGravForces(const int, const int, int *, MeshlessFVParticle<ndim> &,
                                       MeshlessFVParticle<ndim> *) = 0;
  virtual void ComputeStarGravForces(const int, NbodyParticle<ndim> **, MeshlessFVParticle<ndim> &) = 0;
  virtual void CopyDataToGhosts(DomainBox<ndim> &, MeshlessFVParticle<ndim> *) = 0;


  // MeshlessFV array memory allocation functions
  //-----------------------------------------------------------------------------------------------
  void InitialSmoothingLengthGuess(void);
  void ComputeThermalProperties(MeshlessFVParticle<ndim> &);

  FLOAT Timestep(MeshlessFVParticle<ndim> &);
  void UpdatePrimitiveVector(MeshlessFVParticle<ndim> &);
  void UpdateArrayVariables(MeshlessFVParticle<ndim> &);
  void IntegrateConservedVariables(MeshlessFVParticle<ndim> &, FLOAT);
  void EndTimestep(const int, const int, const FLOAT, const FLOAT, MeshlessFVParticle<ndim> *);


  // Functions needed to hide some implementation details
  //-----------------------------------------------------------------------------------------------
  MeshlessFVParticle<ndim>& GetMeshlessFVParticlePointer(const int i) {
    return *((MeshlessFVParticle<ndim>*)((unsigned char*) hydrodata_unsafe + i*size_hydro_part));
  }
  MeshlessFVParticle<ndim>* GetMeshlessFVParticleArray() {return hydrodata;}
  virtual Particle<ndim>* GetParticleArray() {return hydrodata;};


  // Const variables (read in from parameters file)
  //-----------------------------------------------------------------------------------------------
  const bool staticParticles;          ///< ..
  const FLOAT accel_mult;              ///< ..
  const FLOAT courant_mult;            ///< ..
  //const FLOAT gamma_eos;               ///< ..
  //const FLOAT gammam1;                 ///< ..
  const FLOAT h_converge;              ///< h-rho iteration tolerance


  // MeshlessFV particle counters and main particle data array
  //-----------------------------------------------------------------------------------------------
  int create_sinks;                    ///< Create new sink particles?
  int fixed_sink_mass;                 ///< Fix masses of sink particles
  FLOAT msink_fixed;                   ///< Fixed sink mass value
  FLOAT hmin_sink;                     ///< Minimum smoothing length of sinks
  string riemann_solver;               ///< Selected Riemann solver
  string slope_limiter;                ///< Selected slope limiter
  //MeshlessFVType MeshlessFVtype[Nhydrotypes];        ///< Array of MeshlessFV types

  SlopeLimiter<ndim,MeshlessFVParticle> *limiter;
  MeshlessFVParticle<ndim> *hydrodata;

};



//=================================================================================================
//  Class MfvMuscl
/// \brief   Meshless Finite-Volume scheme described by Lanson & Vila (2008)
/// \details Meshless Finite-Volume scheme described by Lanson & Vila (2008).
/// \author  D. A. Hubber, J. Ngoumou
/// \date    25/02/2015
//=================================================================================================
//template <int ndim>
template <int ndim, template<int> class kernelclass>
class MfvMuscl : public MeshlessFV<ndim>
{
 public:

  using MeshlessFV<ndim>::allocated;
  using MeshlessFV<ndim>::gamma_eos;
  using MeshlessFV<ndim>::gammam1;
  using MeshlessFV<ndim>::h_converge;
  using MeshlessFV<ndim>::h_fac;
  using MeshlessFV<ndim>::hydrodata;
  using MeshlessFV<ndim>::hydrodata_unsafe;
  using MeshlessFV<ndim>::invndim;
  using MeshlessFV<ndim>::kernfac;
  using MeshlessFV<ndim>::kernfacsqd;
  using MeshlessFV<ndim>::kernp;
  using MeshlessFV<ndim>::kernrange;
  using MeshlessFV<ndim>::limiter;
  using MeshlessFV<ndim>::mmean;
  using MeshlessFV<ndim>::Ngather;
  using MeshlessFV<ndim>::Nghost;
  using MeshlessFV<ndim>::Nghostmax;
  using MeshlessFV<ndim>::NImportedParticles;
  using MeshlessFV<ndim>::Nhydro;
  using MeshlessFV<ndim>::Nhydromax;
  using MeshlessFV<ndim>::Nmpighost;
  using MeshlessFV<ndim>::NPeriodicGhost;
  using MeshlessFV<ndim>::Ntot;
  using MeshlessFV<ndim>::riemann;
  using MeshlessFV<ndim>::size_hydro_part;
  using MeshlessFV<ndim>::staticParticles;

  static const int nvar = ndim + 2;
  static const int ivx = 0;
  static const int ivy = 1;
  static const int ivz = 2;
  static const int irho = ndim;
  static const int ietot = ndim + 1;
  static const int ipress = ndim + 1;



  // Constructor
  //-----------------------------------------------------------------------------------------------
  MfvMuscl(int hydro_forces_aux, int self_gravity_aux, FLOAT _accel_mult, FLOAT _courant_mult,
            FLOAT h_fac_aux, FLOAT h_converge_aux,
            FLOAT gamma_aux, string gas_eos_aux, string KernelName, int size_MeshlessFV_part);
  ~MfvMuscl();



  // MeshlessFV functions for computing MeshlessFV sums with neighbouring particles
  // (fully coded in each separate MeshlessFV implementation, and not in MeshlessFV.cpp)
  //-----------------------------------------------------------------------------------------------
  int ComputeH(const int, const int, const FLOAT, FLOAT *, FLOAT *, FLOAT *, FLOAT *,
               MeshlessFVParticle<ndim> &, Nbody<ndim> *);
  void ComputePsiFactors(const int, const int, int *, FLOAT *, FLOAT *, FLOAT *,
                         MeshlessFVParticle<ndim> &, MeshlessFVParticle<ndim> *);
  void ComputeGradients(const int, const int, int *, FLOAT *, FLOAT *, FLOAT *,
                                    MeshlessFVParticle<ndim> &, MeshlessFVParticle<ndim> *);
  void ComputeGodunovFlux(const int, const int, const FLOAT, int *, FLOAT *, FLOAT *, FLOAT *,
                          MeshlessFVParticle<ndim> &, MeshlessFVParticle<ndim> *);
  void CopyDataToGhosts(DomainBox<ndim> &, MeshlessFVParticle<ndim> *);
  void ComputeSmoothedGravForces(const int, const int, int *,
                                 MeshlessFVParticle<ndim> &, MeshlessFVParticle<ndim> *);
  void ComputeDirectGravForces(const int, const int, int *,
                               MeshlessFVParticle<ndim> &, MeshlessFVParticle<ndim> *);
  void ComputeStarGravForces(const int, NbodyParticle<ndim> **, MeshlessFVParticle<ndim> &);

  kernelclass<ndim> kern;                  ///< SPH kernel

};



//=================================================================================================
//  Class MfvRungeKutta
/// \brief   Meshless Finite-Volume scheme described by Lanson & Vila (2008)
/// \details Meshless Finite-Volume scheme described by Lanson & Vila (2008).
/// \author  D. A. Hubber, J. Ngoumou
/// \date    25/02/2015
//=================================================================================================
//template <int ndim>
template <int ndim, template<int> class kernelclass>
class MfvRungeKutta : public MeshlessFV<ndim>
{
 public:

  using MeshlessFV<ndim>::allocated;
  using MeshlessFV<ndim>::gamma_eos;
  using MeshlessFV<ndim>::gammam1;
  using MeshlessFV<ndim>::h_converge;
  using MeshlessFV<ndim>::h_fac;
  using MeshlessFV<ndim>::hydrodata;
  using MeshlessFV<ndim>::hydrodata_unsafe;
  using MeshlessFV<ndim>::invndim;
  using MeshlessFV<ndim>::kernfac;
  using MeshlessFV<ndim>::kernfacsqd;
  using MeshlessFV<ndim>::kernp;
  using MeshlessFV<ndim>::kernrange;
  using MeshlessFV<ndim>::limiter;
  using MeshlessFV<ndim>::mmean;
  using MeshlessFV<ndim>::Ngather;
  using MeshlessFV<ndim>::Nghost;
  using MeshlessFV<ndim>::Nghostmax;
  using MeshlessFV<ndim>::NImportedParticles;
  using MeshlessFV<ndim>::Nhydro;
  using MeshlessFV<ndim>::Nhydromax;
  using MeshlessFV<ndim>::Nmpighost;
  using MeshlessFV<ndim>::NPeriodicGhost;
  using MeshlessFV<ndim>::Ntot;
  using MeshlessFV<ndim>::riemann;
  using MeshlessFV<ndim>::size_hydro_part;
  using MeshlessFV<ndim>::staticParticles;

  static const int nvar = ndim + 2;
  static const int ivx = 0;
  static const int ivy = 1;
  static const int ivz = 2;
  static const int irho = ndim;
  static const int ietot = ndim + 1;
  static const int ipress = ndim + 1;


  // Constructor
  //-----------------------------------------------------------------------------------------------
  MfvRungeKutta(int hydro_forces_aux, int self_gravity_aux, FLOAT _accel_mult, FLOAT _courant_mult,
                FLOAT h_fac_aux, FLOAT h_converge_aux,
                FLOAT gamma_aux, string gas_eos_aux, string KernelName, int size_MeshlessFV_part);
  ~MfvRungeKutta();



  // MeshlessFV functions for computing MeshlessFV sums with neighbouring particles
  // (fully coded in each separate MeshlessFV implementation, and not in MeshlessFV.cpp)
  //-----------------------------------------------------------------------------------------------
  int ComputeH(const int, const int, const FLOAT, FLOAT *, FLOAT *, FLOAT *, FLOAT *,
               MeshlessFVParticle<ndim> &, Nbody<ndim> *);
  void ComputePsiFactors(const int, const int, int *, FLOAT *, FLOAT *, FLOAT *,
                         MeshlessFVParticle<ndim> &, MeshlessFVParticle<ndim> *);
  void ComputeGradients(const int, const int, int *, FLOAT *, FLOAT *, FLOAT *,
                                    MeshlessFVParticle<ndim> &, MeshlessFVParticle<ndim> *);
  void ComputeGodunovFlux(const int, const int, const FLOAT, int *, FLOAT *, FLOAT *, FLOAT *,
                          MeshlessFVParticle<ndim> &, MeshlessFVParticle<ndim> *);
  void CopyDataToGhosts(DomainBox<ndim> &, MeshlessFVParticle<ndim> *);
  void ComputeSmoothedGravForces(const int, const int, int *,
                                 MeshlessFVParticle<ndim> &, MeshlessFVParticle<ndim> *);
  void ComputeDirectGravForces(const int, const int, int *,
                               MeshlessFVParticle<ndim> &, MeshlessFVParticle<ndim> *);
  void ComputeStarGravForces(const int, NbodyParticle<ndim> **, MeshlessFVParticle<ndim> &);

  kernelclass<ndim> kern;                  ///< SPH kernel

};
#endif
