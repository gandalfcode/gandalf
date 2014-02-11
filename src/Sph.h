//=============================================================================
//  Sph.h
//  Contains main parent virtual class plus child classes for various SPH 
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
//=============================================================================


#ifndef _SPH_H_
#define _SPH_H_


#include <assert.h>
#include <string>
#include "Precision.h"
#include "Constants.h"
#include "SphParticle.h"
#include "SphKernel.h"
#include "NbodyParticle.h"
#include "Nbody.h"
#include "Parameters.h"
#include "EOS.h"
#include "RiemannSolver.h"
using namespace std;
#if defined _OPENMP
#include "omp.h"
#endif

enum aviscenum{noneav, mon97, mon97td};
enum acondenum{noneac, wadsley2008, price2008};



//=============================================================================
//  Class Sph
/// \brief   Main parent Sph class.
/// \details Different SPH implementations (e.g. grad-h SPH, Saitoh &
///          Makino 2012) are derived from this class.  Each implementation
///          requires defining its own version of each function (e.g. ComputeH
///          for its own method of computing smoothing lengths).
/// \author  D. A. Hubber, G. Rosotti
/// \date    03/04/2013
//=============================================================================
template <int ndim>
class Sph
{
private:

#if defined _OPENMP
  void InitParticleLocks();
  void DestroyParticleLocks();
#endif

 public:

  const aviscenum avisc;
  const acondenum acond;

  // Constructor
  //---------------------------------------------------------------------------
  Sph(int hydro_forces_aux, int self_gravity_aux, FLOAT alpha_visc_aux, 
      FLOAT beta_visc_aux, FLOAT h_fac_aux, FLOAT h_converge_aux, 
      aviscenum avisc_aux, acondenum acond_aux, string gas_eos_aux, 
      string KernelName);


  // SPH functions for computing SPH sums with neighbouring particles 
  // (fully coded in each separate SPH implementation, and not in Sph.cpp)
  //---------------------------------------------------------------------------
  virtual int ComputeH(int, int, FLOAT, FLOAT *, FLOAT *, FLOAT *, FLOAT *,
                       SphParticle<ndim> &, Nbody<ndim> *) = 0;
  virtual void ComputeSphHydroForces(int, int, int *, FLOAT *, FLOAT *, 
                                     FLOAT *, SphParticle<ndim> &,
                                     SphParticle<ndim> *) = 0;
  virtual void ComputeSphHydroGravForces(int, int, int *, SphParticle<ndim> &, 
					 SphParticle<ndim> *) = 0;
  virtual void ComputeSphGravForces(int, int, int *, SphParticle<ndim> &,
				    SphParticle<ndim> *) = 0;
  virtual void ComputeDirectGravForces(int, int, int *, SphParticle<ndim> &, 
                                       SphParticle<ndim> *) = 0;
  virtual void ComputeSphNeibDudt(int, int, int *, FLOAT *, FLOAT *,
				  FLOAT *, SphParticle<ndim> &, 
                                  SphParticle<ndim> *) = 0;
  virtual void ComputeSphDerivatives(int, int, int *, FLOAT *, FLOAT *,
				     FLOAT *, SphParticle<ndim> &, 
                                     SphParticle<ndim> *) = 0;
  virtual void ComputeStarGravForces(int, NbodyParticle<ndim> **, 
				     SphParticle<ndim> &) = 0;


  // SPH array memory allocation functions
  //---------------------------------------------------------------------------
  void AllocateMemory(int);
  void DeallocateMemory(void);
  void DeleteParticles(int, int *);
  void ReorderParticles(void);
  void SphBoundingBox(FLOAT *, FLOAT *, int);
  void InitialSmoothingLengthGuess(void);


  // Functions needed to hide some implementation details
  //---------------------------------------------------------------------------
  SphParticle<ndim>* GetParticleIPointer(int i) {return &sphdata[i];};
#if defined _OPENMP
  omp_lock_t& GetParticleILock(int i) {return locks[i];};
  omp_lock_t* locks;
#endif


  // SPH particle counters and main particle data array
  //---------------------------------------------------------------------------
  bool allocated;                     ///< Is SPH memory allocated?
  int Ngather;                        ///< Average no. of gather neighbours
  int Nsph;                           ///< No. of SPH particles in simulation
  int Nghost;                         ///< No. of ghost SPH particles
  int NPeriodicGhost;                 ///< No. of periodic ghost particles
  int Ntot;                           ///< No. of real + ghost particles
  int Nsphmax;                        ///< Max. no. of SPH particles in array
  int Nghostmax;                      ///< Max. allowed no. of ghost particles

  const FLOAT alpha_visc;             ///< alpha artificial viscosity parameter
  const FLOAT beta_visc;              ///< beta artificial viscosity parameter
  const FLOAT h_fac;                  ///< Smoothing length-density factor
  const FLOAT h_converge;             ///< h-rho iteration tolerance
  const string gas_eos;               ///< Gas EOS option
  const int hydro_forces;             ///< Compute hydro forces?
  const int self_gravity;             ///< Compute gravitational forces?
  static const FLOAT invndim=1./ndim; ///< Copy of 1/ndim
  int create_sinks;                   ///< Create new sink particles?
  int time_dependent_avisc;           ///< Use time-dependent viscosity?
  FLOAT mmean;                        ///< Mean SPH particle mass
  FLOAT hmin_sink;                    ///< Minimum smoothing length of sinks

  string riemann_solver;              ///< Selected Riemann solver
  string slope_limiter;               ///< Selected slope limiter
  int riemann_order;                  ///< Order of Riemann solver
  FLOAT alpha_visc_min;               ///< Min. time-dependent viscosity alpha
  FLOAT kernfac;                      ///< Kernel range neighbour fraction
  FLOAT kernfacsqd;                   ///< Kernel range neib. fraction squared

  int *iorder;                        ///< Array containing particle ordering
  FLOAT *rsph;                        ///< Position array (for efficiency)

  SphIntParticle<ndim>* sphintdata;   ///< Pointer to particle integration data
  SphParticle<ndim> *sphdata;         ///< Pointer to particle data
  SphKernel<ndim> *kernp;             ///< Pointer to chosen kernel object
  TabulatedKernel<ndim> kerntab;      ///< Tabulated version of chosen kernel
  EOS<ndim> *eos;                     ///< Equation-of-state
  RiemannSolver *riemann;             ///< Riemann solver

};



//=============================================================================
//  Class GradhSph
/// \brief   Class definition for conservative 'grad-h' SPH simulations.
/// \details Class definition for conservative 'grad-h' SPH simulations 
///          (as derived from the parent Sph class).  Full code for each of 
///          these class functions written in 'GradhSph.cpp'.
/// \author  D. A. Hubber, G. Rosotti
/// \date    03/04/2013
//=============================================================================
#if !defined(SWIG)
template <int ndim, template<int> class kernelclass>
class GradhSph: public Sph<ndim>
{
  using Sph<ndim>::allocated;
  using Sph<ndim>::Nsph;
  using Sph<ndim>::Ntot;
  using Sph<ndim>::sphdata;
  using Sph<ndim>::eos;
  using Sph<ndim>::h_fac;
  using Sph<ndim>::kernfacsqd;
  using Sph<ndim>::invndim;
  using Sph<ndim>::h_converge;
  using Sph<ndim>::hydro_forces;
  using Sph<ndim>::self_gravity;
  using Sph<ndim>::avisc;
  using Sph<ndim>::beta_visc;
  using Sph<ndim>::alpha_visc;
  using Sph<ndim>::alpha_visc_min;
  using Sph<ndim>::acond;
  using Sph<ndim>::create_sinks;
  using Sph<ndim>::hmin_sink;
  using Sph<ndim>::time_dependent_avisc;

 public:

  GradhSph(int, int, FLOAT, FLOAT, FLOAT, FLOAT,
           aviscenum, acondenum, string, string);
  ~GradhSph();

  int ComputeH(int, int, FLOAT, FLOAT *, FLOAT *, FLOAT *, FLOAT *,
               SphParticle<ndim> &, Nbody<ndim> *);
  void ComputeSphGravForces(int, int, int *, SphParticle<ndim> &,
			    SphParticle<ndim> *);
  void ComputeSphHydroGravForces(int, int, int *, SphParticle<ndim> &, 
				 SphParticle<ndim> *);
  void ComputeSphHydroForces(int, int, int *, FLOAT *, FLOAT *, FLOAT *,
			     SphParticle<ndim> &, SphParticle<ndim> *);
  void ComputeSphNeibDudt(int, int, int *, FLOAT *, FLOAT *,
			  FLOAT *, SphParticle<ndim> &, SphParticle<ndim> *);
  void ComputeSphDerivatives(int, int, int *, FLOAT *, FLOAT *, FLOAT *, 
			     SphParticle<ndim> &, SphParticle<ndim> *);
  void ComputeDirectGravForces(int, int, int *, 
			       SphParticle<ndim> &, SphParticle<ndim> *);
  void ComputeStarGravForces(int, NbodyParticle<ndim> **, SphParticle<ndim> &);

  kernelclass<ndim> kern;                  ///< SPH kernel

};



//=============================================================================
//  Class SM2012Sph
/// \brief   Class definition for Saitoh & Makino (2012) SPH simulations
/// \details Class definition for Saitoh & Makino (2012) SPH simulations 
///          (as derived from the parent Sph class).  Full code for each of 
///          these class functions written in 'SM2012Sph.cpp'.
/// \author  D. A. Hubber, G. Rosotti
/// \date    03/04/2013
//=============================================================================
template <int ndim, template<int> class kernelclass>
class SM2012Sph: public Sph<ndim>
{
  using Sph<ndim>::allocated;
  using Sph<ndim>::Nsph;
  using Sph<ndim>::Ntot;
  using Sph<ndim>::eos;
  using Sph<ndim>::h_fac;
  using Sph<ndim>::kernfacsqd;
  using Sph<ndim>::invndim;
  using Sph<ndim>::h_converge;
  using Sph<ndim>::hydro_forces;
  using Sph<ndim>::avisc;
  using Sph<ndim>::beta_visc;
  using Sph<ndim>::alpha_visc;
  using Sph<ndim>::alpha_visc_min;
  using Sph<ndim>::acond;
  using Sph<ndim>::create_sinks;
  using Sph<ndim>::hmin_sink;
  using Sph<ndim>::time_dependent_avisc;

 public:

  SM2012Sph(int, int, FLOAT, FLOAT, FLOAT, FLOAT,
            aviscenum, acondenum, string, string);
  ~SM2012Sph();

  int ComputeH(int, int, FLOAT, FLOAT *, FLOAT *, FLOAT *, FLOAT *,
               SphParticle<ndim> &, Nbody<ndim> *);
  void ComputeSphHydroForces(int, int, int *, FLOAT *, FLOAT *, FLOAT *,
			     SphParticle<ndim> &, SphParticle<ndim> *);
  void ComputeSphHydroGravForces(int, int, int *, SphParticle<ndim> &, 
				 SphParticle<ndim> *);
  void ComputeSphGravForces(int, int, int *,
			    SphParticle<ndim> &, SphParticle<ndim> *);
  void ComputeSphNeibDudt(int, int, int *, FLOAT *, FLOAT *,
			  FLOAT *, SphParticle<ndim> &, SphParticle<ndim> *);
  void ComputeSphDerivatives(int, int, int *, FLOAT *, FLOAT *, FLOAT *, 
			     SphParticle<ndim> &, SphParticle<ndim> *);
  void ComputeDirectGravForces(int, int, int *,
                               SphParticle<ndim> &, SphParticle<ndim> *);
  void ComputeStarGravForces(int, NbodyParticle<ndim> **, SphParticle<ndim> &);


  kernelclass<ndim> kern;                  ///< SPH kernel

};



//=============================================================================
//  Class GodunovSph
/// Class definition for Godunov SPH (Inutsuka 2002) algorithm.
/// Full code for each of these class functions
/// written in 'GodunovSph.cpp'.
//=============================================================================
template <int ndim, template<int> class kernelclass>
class GodunovSph: public Sph<ndim>
{
  using Sph<ndim>::allocated;
  using Sph<ndim>::Nsph;
  using Sph<ndim>::Ntot;
  using Sph<ndim>::eos;
  using Sph<ndim>::h_fac;
  using Sph<ndim>::kernfacsqd;
  using Sph<ndim>::invndim;
  using Sph<ndim>::h_converge;
  using Sph<ndim>::hydro_forces;
  using Sph<ndim>::avisc;
  using Sph<ndim>::beta_visc;
  using Sph<ndim>::alpha_visc;
  using Sph<ndim>::acond;
  using Sph<ndim>::riemann;
  using Sph<ndim>::riemann_solver;
  using Sph<ndim>::riemann_order;
  using Sph<ndim>::slope_limiter;
  using Sph<ndim>::create_sinks;
  using Sph<ndim>::hmin_sink;

 public:

  GodunovSph(int, int, FLOAT, FLOAT, FLOAT, FLOAT,
             aviscenum, acondenum, string, string);
  ~GodunovSph();

  int ComputeH(int, int, FLOAT, FLOAT *, FLOAT *, FLOAT *, FLOAT *,
               SphParticle<ndim> &, Nbody<ndim> *);
  void ComputeSphHydroForces(int, int, int *, FLOAT *, FLOAT *, FLOAT *, 
			     SphParticle<ndim> &, SphParticle<ndim> *);
  void ComputeSphHydroGravForces(int, int, int *, SphParticle<ndim> &, 
				 SphParticle<ndim> *);
  void ComputeSphGravForces(int, int, int *,
			    SphParticle<ndim> &, SphParticle<ndim> *);
  void ComputeSphNeibDudt(int, int, int *, FLOAT *, FLOAT *,
  			  FLOAT *, SphParticle<ndim> &, SphParticle<ndim> *);
  void ComputeSphDerivatives(int, int, int *, FLOAT *, FLOAT *, FLOAT *, 
			     SphParticle<ndim> &, SphParticle<ndim> *);
  void ComputeDirectGravForces(int, int, int *, 
                               SphParticle<ndim> &, SphParticle<ndim> *);
  void InitialiseRiemannProblem(SphParticle<ndim>, SphParticle<ndim>, FLOAT *,
                                FLOAT, FLOAT, FLOAT, FLOAT, FLOAT, FLOAT &, 
                                FLOAT &, FLOAT &, FLOAT &, FLOAT &, FLOAT &);
  void ComputeStarGravForces(int, NbodyParticle<ndim> **, SphParticle<ndim> &);

  kernelclass<ndim> kern;                 ///< SPH kernel

};
#endif


#endif
