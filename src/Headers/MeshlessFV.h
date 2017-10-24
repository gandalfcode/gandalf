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
#include "NeighbourManager.h"
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
/// \details Main parent MeshlessFV class.
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
  using Hydrodynamics<ndim>::kernp;
  using Hydrodynamics<ndim>::kernrange;
  using Hydrodynamics<ndim>::mmean;
  using Hydrodynamics<ndim>::Ngather;
  using Hydrodynamics<ndim>::Nghost;
  using Hydrodynamics<ndim>::NImportedParticles;
  using Hydrodynamics<ndim>::Nhydro;
  using Hydrodynamics<ndim>::Nhydromax;
  using Hydrodynamics<ndim>::Nmpighost;
  using Hydrodynamics<ndim>::NPeriodicGhost;
  using Hydrodynamics<ndim>::Ntot;
  using Hydrodynamics<ndim>::self_gravity;
  using Hydrodynamics<ndim>::size_hydro_part;
  using Hydrodynamics<ndim>::types;

  static const int nvar = ndim + 2;
  static const int ivx = 0;
  static const int ivy = 1;
  static const int ivz = 2;
  static const int irho = ndim;
  static const int ietot = ndim + 1;
  static const int ipress = ndim + 1;

  typedef MeshlessFVParticle<ndim> ParticleType ;
  typedef typename ParticleType::FluxParticle     FluxNeib ;
  typedef typename ParticleType::GradientParticle GradNeib ;
  typedef typename ParticleType::GravParticle     GravNeib ;
  typedef typename GravityNeighbourLists<GravNeib>::DirectType DirectNeib ;

  // Constructor
  //-----------------------------------------------------------------------------------------------
  MeshlessFV(int hydro_forces_aux, int self_gravity_aux, FLOAT _accel_mult, FLOAT _courant_mult,
             FLOAT h_fac_aux, FLOAT h_converge_aux, FLOAT gamma_aux, string gas_eos_aux,
             string KernelName, int size_mfv_part, SimUnits &units, Parameters *params);
 virtual ~MeshlessFV();

  virtual void AllocateMemory(int);
  virtual void DeallocateMemory(void);
  virtual int DeleteDeadParticles(void) {
    return this->template DoDeleteDeadParticles<MeshlessFVParticle>() ;
  }
  virtual void AccreteMassFromParticle(const FLOAT dm, Particle<ndim> &part) {
    MeshlessFVParticle<ndim>& mfvpart = static_cast<MeshlessFVParticle<ndim>& > (part);
    mfvpart.m -= dm;
    mfvpart.dQ[irho] -= dm;
    //mfvpart.dQ[ietot] -= dm*part.u;
    for (int k=0; k<ndim; k++) mfvpart.dQ[k] -= dm*mfvpart.v[k];
  }

  virtual void ZeroAccelerations();



  // MeshlessFV functions for computing MeshlessFV sums with neighbouring particles
  // (fully coded in each separate MeshlessFV implementation, and not in MeshlessFV.cpp)
  //-----------------------------------------------------------------------------------------------
  virtual int ComputeH(const int, const int, const FLOAT, FLOAT *, FLOAT *, FLOAT *, FLOAT *,
                       MeshlessFVParticle<ndim> &, Nbody<ndim> *) = 0;
  virtual void ComputeGradients(MeshlessFVParticle<ndim>&, NeighbourList<GradNeib>&) = 0;
  virtual void ComputeGodunovFlux(MeshlessFVParticle<ndim>&, NeighbourList<FluxNeib>&, FLOAT) = 0;
  virtual void ComputeSmoothedGravForces(MeshlessFVParticle<ndim>&, NeighbourList<GravNeib>&) = 0;
  virtual void ComputeDirectGravForces(MeshlessFVParticle<ndim>&, NeighbourList<DirectNeib>&) = 0;
  virtual void ComputeStarGravForces(const int, NbodyParticle<ndim> **, MeshlessFVParticle<ndim> &) = 0;


  // Other functions.
  //-----------------------------------------------------------------------------------------------
  void ComputeThermalProperties(MeshlessFVParticle<ndim> &);
  void InitialSmoothingLengthGuess(void);
  void UpdatePrimitiveVector(MeshlessFVParticle<ndim> &);
  void UpdateArrayVariables(MeshlessFVParticle<ndim> &, FLOAT [nvar]);


  // Functions needed to hide some implementation details
  //-----------------------------------------------------------------------------------------------
  MeshlessFVParticle<ndim>& GetMeshlessFVParticlePointer(const int i) {
	  long int numBytes = (long int) i * (long int) size_hydro_part;
    return *((MeshlessFVParticle<ndim>*)((unsigned char*) hydrodata_unsafe + numBytes));
  }
  MeshlessFVParticle<ndim>* GetMeshlessFVParticleArray() {return hydrodata;}


  // Const variables (read in from parameters file)
  //-----------------------------------------------------------------------------------------------
  const bool staticParticles;          ///< ..
  const FLOAT accel_mult;              ///< ..
  const FLOAT courant_mult;            ///< ..
  const FLOAT h_converge;              ///< h-rho iteration tolerance


  // MeshlessFV particle counters and main particle data array
  //-----------------------------------------------------------------------------------------------
  int create_sinks;                    ///< Create new sink particles?
  int fixed_sink_mass;                 ///< Fix masses of sink particles
  FLOAT msink_fixed;                   ///< Fixed sink mass value
  FLOAT hmin_sink;                     ///< Minimum smoothing length of sinks
  string riemann_solver;               ///< Selected Riemann solver
  string slope_limiter;                ///< Selected slope limiter
  string timestep_limiter;            ///< Time of limiter used for block timesteps

  MeshlessFVParticle<ndim> *hydrodata;
  CodeTiming* timing;

};



//=================================================================================================
//  Class MfvCommon
/// \brief   Intermediate class for Meshless FV scheme containing all common functions.
/// \details Intermediate class for Meshless FV scheme containing all common functions.
/// \author  D. A. Hubber, J. Ngoumou
/// \date    25/02/2015
//=================================================================================================
//template <int ndim>
template <int ndim, template<int> class kernelclass, class SlopeLimiterType>
class MfvCommon : public MeshlessFV<ndim>
{
 public:

  using MeshlessFV<ndim>::allocated;
  using MeshlessFV<ndim>::h_converge;
  using MeshlessFV<ndim>::h_fac;
  using MeshlessFV<ndim>::hydrodata;
  using MeshlessFV<ndim>::hydrodata_unsafe;
  using MeshlessFV<ndim>::invndim;
  using MeshlessFV<ndim>::kernp;
  using MeshlessFV<ndim>::kernrange;
  using MeshlessFV<ndim>::mmean;
  using MeshlessFV<ndim>::Ngather;
  using MeshlessFV<ndim>::Nghost;
  using MeshlessFV<ndim>::NImportedParticles;
  using MeshlessFV<ndim>::Nhydro;
  using MeshlessFV<ndim>::Nhydromax;
  using MeshlessFV<ndim>::Nmpighost;
  using MeshlessFV<ndim>::NPeriodicGhost;
  using MeshlessFV<ndim>::Ntot;
  using MeshlessFV<ndim>::size_hydro_part;
  using MeshlessFV<ndim>::staticParticles;
  using Hydrodynamics<ndim>::create_sinks;
  using Hydrodynamics<ndim>::hmin_sink;
  using Hydrodynamics<ndim>::types;

  static const int nvar = ndim + 2;
  static const int ivx = 0;
  static const int ivy = 1;
  static const int ivz = 2;
  static const int irho = ndim;
  static const int ietot = ndim + 1;
  static const int ipress = ndim + 1;

  typedef MeshlessFVParticle<ndim> ParticleType ;
  typedef typename MeshlessFV<ndim>::FluxNeib     FluxNeib ;
  typedef typename MeshlessFV<ndim>::GradNeib     GradNeib ;
  typedef typename MeshlessFV<ndim>::GravNeib     GravNeib ;
  typedef typename MeshlessFV<ndim>::DirectNeib   DirectNeib ;

  // Constructor
  //-----------------------------------------------------------------------------------------------
  MfvCommon(int hydro_forces_aux, int self_gravity_aux, FLOAT _accel_mult, FLOAT _courant_mult,
           FLOAT h_fac_aux, FLOAT h_converge_aux, FLOAT gamma_aux, string gas_eos_aux,
           string KernelName, int size_MeshlessFV_part, SimUnits &units, Parameters *params);
  virtual ~MfvCommon();

  // MeshlessFV functions for computing MeshlessFV sums with neighbouring particles
  // (fully coded in each separate MeshlessFV implementation, and not in MeshlessFV.cpp)
  //-----------------------------------------------------------------------------------------------
  int ComputeH(const int, const int, const FLOAT, FLOAT *, FLOAT *, FLOAT *, FLOAT *,
               MeshlessFVParticle<ndim> &, Nbody<ndim> *);
  virtual void ComputeGradients(MeshlessFVParticle<ndim>&, NeighbourList<GradNeib>&);
  virtual void ComputeSmoothedGravForces(MeshlessFVParticle<ndim>&, NeighbourList<GravNeib>&);
  virtual void ComputeDirectGravForces(MeshlessFVParticle<ndim>&, NeighbourList<DirectNeib>&);
  void ComputeStarGravForces(const int, NbodyParticle<ndim> **, MeshlessFVParticle<ndim> &);

  kernelclass<ndim> kern;                  ///< SPH kernel
  SlopeLimiterType limiter;
  ExactRiemannSolver<ndim> riemannExact ;
  HllcRiemannSolver<ndim> riemannHLLC ;
  int RiemannSolverType ;
  ViscousFlux<ndim> viscosity;
  bool need_viscosity;


};



//=================================================================================================
//  Class MfvMuscl
/// \brief   Meshless Finite-Volume scheme described by Lanson & Vila (2008)
/// \details Meshless Finite-Volume scheme described by Lanson & Vila (2008).
/// \author  D. A. Hubber, J. Ngoumou
/// \date    25/02/2015
//=================================================================================================
//template <int ndim>
template <int ndim, template<int> class kernelclass, class SlopeLimiterType>
class MfvMuscl : public MfvCommon<ndim,kernelclass,SlopeLimiterType>
{
 public:

  using Hydrodynamics<ndim>::create_sinks;
  using Hydrodynamics<ndim>::hmin_sink;
  using MeshlessFV<ndim>::allocated;
  using MeshlessFV<ndim>::h_converge;
  using MeshlessFV<ndim>::h_fac;
  using MeshlessFV<ndim>::hydrodata;
  using MeshlessFV<ndim>::hydrodata_unsafe;
  using MeshlessFV<ndim>::invndim;
  using MeshlessFV<ndim>::kernp;
  using MeshlessFV<ndim>::kernrange;
  using MeshlessFV<ndim>::mmean;
  using MeshlessFV<ndim>::Ngather;
  using MeshlessFV<ndim>::Nghost;
  using MeshlessFV<ndim>::NImportedParticles;
  using MeshlessFV<ndim>::Nhydro;
  using MeshlessFV<ndim>::Nhydromax;
  using MeshlessFV<ndim>::Nmpighost;
  using MeshlessFV<ndim>::NPeriodicGhost;
  using MeshlessFV<ndim>::Ntot;
  using MeshlessFV<ndim>::size_hydro_part;
  using MeshlessFV<ndim>::staticParticles;
  using MfvCommon<ndim,kernelclass,SlopeLimiterType>::limiter;
  using MfvCommon<ndim,kernelclass,SlopeLimiterType>::kern;
  using MfvCommon<ndim,kernelclass,SlopeLimiterType>::riemannExact;
  using MfvCommon<ndim,kernelclass,SlopeLimiterType>::riemannHLLC;
  using MfvCommon<ndim,kernelclass,SlopeLimiterType>::RiemannSolverType;
  using MfvCommon<ndim,kernelclass,SlopeLimiterType>::viscosity;
  using MfvCommon<ndim,kernelclass,SlopeLimiterType>::need_viscosity;

  static const int nvar = ndim + 2;
  static const int ivx = 0;
  static const int ivy = 1;
  static const int ivz = 2;
  static const int irho = ndim;
  static const int ietot = ndim + 1;
  static const int ipress = ndim + 1;

  typedef MeshlessFVParticle<ndim> ParticleType ;
  typedef typename MeshlessFV<ndim>::FluxNeib     FluxNeib ;
  typedef typename MeshlessFV<ndim>::GradNeib     GradNeib ;
  typedef typename MeshlessFV<ndim>::GravNeib     GravNeib ;
  typedef typename MeshlessFV<ndim>::DirectNeib   DirectNeib ;



  // Constructor
  //-----------------------------------------------------------------------------------------------
  MfvMuscl(int hydro_forces_aux, int self_gravity_aux, FLOAT _accel_mult, FLOAT _courant_mult,
           FLOAT h_fac_aux, FLOAT h_converge_aux, FLOAT gamma_aux, string gas_eos_aux,
           string KernelName, int size_MeshlessFV_part, SimUnits &units, Parameters *params);
  ~MfvMuscl();


  //-----------------------------------------------------------------------------------------------
  virtual void ComputeGodunovFlux(MeshlessFVParticle<ndim>&, NeighbourList<FluxNeib>&, FLOAT);


};



//=================================================================================================
//  Class MfvRungeKutta
/// \brief   Meshless Finite-Volume scheme described by Lanson & Vila (2008)
/// \details Meshless Finite-Volume scheme described by Lanson & Vila (2008).
/// \author  D. A. Hubber, J. Ngoumou
/// \date    25/02/2015
//=================================================================================================
//template <int ndim>
template <int ndim, template<int> class kernelclass, class SlopeLimiterType>
class MfvRungeKutta : public MfvCommon<ndim,kernelclass,SlopeLimiterType>
{
 public:

  using Hydrodynamics<ndim>::create_sinks;
  using Hydrodynamics<ndim>::hmin_sink;
  using MeshlessFV<ndim>::allocated;
  using MeshlessFV<ndim>::h_converge;
  using MeshlessFV<ndim>::h_fac;
  using MeshlessFV<ndim>::hydrodata;
  using MeshlessFV<ndim>::hydrodata_unsafe;
  using MeshlessFV<ndim>::invndim;
  using MeshlessFV<ndim>::kernp;
  using MeshlessFV<ndim>::kernrange;
  using MeshlessFV<ndim>::mmean;
  using MeshlessFV<ndim>::Ngather;
  using MeshlessFV<ndim>::Nghost;
  using MeshlessFV<ndim>::NImportedParticles;
  using MeshlessFV<ndim>::Nhydro;
  using MeshlessFV<ndim>::Nhydromax;
  using MeshlessFV<ndim>::Nmpighost;
  using MeshlessFV<ndim>::NPeriodicGhost;
  using MeshlessFV<ndim>::Ntot;
  using MeshlessFV<ndim>::size_hydro_part;
  using MeshlessFV<ndim>::staticParticles;
  using MfvCommon<ndim,kernelclass,SlopeLimiterType>::limiter;
  using MfvCommon<ndim,kernelclass,SlopeLimiterType>::kern;
  using MfvCommon<ndim,kernelclass,SlopeLimiterType>::riemannExact;
  using MfvCommon<ndim,kernelclass,SlopeLimiterType>::riemannHLLC;
  using MfvCommon<ndim,kernelclass,SlopeLimiterType>::RiemannSolverType;
  using MfvCommon<ndim,kernelclass,SlopeLimiterType>::viscosity;
  using MfvCommon<ndim,kernelclass,SlopeLimiterType>::need_viscosity;

  static const int nvar = ndim + 2;
  static const int ivx = 0;
  static const int ivy = 1;
  static const int ivz = 2;
  static const int irho = ndim;
  static const int ietot = ndim + 1;
  static const int ipress = ndim + 1;

  typedef MeshlessFVParticle<ndim> ParticleType ;
  typedef typename MeshlessFV<ndim>::FluxNeib     FluxNeib ;
  typedef typename MeshlessFV<ndim>::GradNeib     GradNeib ;
  typedef typename MeshlessFV<ndim>::GravNeib     GravNeib ;
  typedef typename MeshlessFV<ndim>::DirectNeib   DirectNeib ;

  // Constructor
  //-----------------------------------------------------------------------------------------------
  MfvRungeKutta(int hydro_forces_aux, int self_gravity_aux, FLOAT _accel_mult, FLOAT _courant_mult,
                FLOAT h_fac_aux, FLOAT h_converge_aux, FLOAT gamma_aux, string gas_eos_aux,
                string KernelName, int size_MeshlessFV_part, SimUnits &units, Parameters *params);
  ~MfvRungeKutta() {};


  // ..
  //-----------------------------------------------------------------------------------------------
  virtual void ComputeGodunovFlux(MeshlessFVParticle<ndim>&, NeighbourList<FluxNeib>&,FLOAT);


};


//=================================================================================================
//  MeshlessFV::UpdatePrimitiveVector
/// Updates the primitive vector from particle quantities.
//=================================================================================================
template <int ndim>
inline void MeshlessFV<ndim>::UpdatePrimitiveVector(MeshlessFVParticle<ndim> &part)
{
  for (int k=0; k<ndim; k++) part.Wprim[k] = part.v[k];
  part.Wprim[irho] = part.rho;
  part.Wprim[ipress] = part.pressure;
}



#endif
