//=================================================================================================
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
//=================================================================================================


#ifndef _SPH_H_
#define _SPH_H_


#include <assert.h>
#include <string>
#include "Constants.h"
#include "DomainBox.h"
#include "EOS.h"
#include "ExternalPotential.h"
#include "Hydrodynamics.h"
#include "Particle.h"
#include "NbodyParticle.h"
#include "Nbody.h"
#include "Parameters.h"
#include "Precision.h"
#include "RiemannSolver.h"
#include "SimUnits.h"
#include "SmoothingKernel.h"
#if defined _OPENMP
#include "omp.h"
#endif
using namespace std;


enum aviscenum{noav, mon97, mon97mm97, mon97cd2010};
enum acondenum{noac, wadsley2008, price2008};
enum tdaviscenum{notdav, mm97, cd2010};


//=================================================================================================
//  Class Sph
/// \brief   Main parent Sph class.
/// \details Different SPH implementations (e.g. grad-h SPH, Saitoh & Makino 2012) are derived from
///          this class.  Each implementation requires defining its own version of each function
///          (e.g. ComputeH for its own method of computing smoothing lengths).
/// \author  D. A. Hubber, G. Rosotti
/// \date    03/04/2013
//=================================================================================================
template <int ndim>
class Sph : public Hydrodynamics<ndim>
{
 private:
  const int size_sph_part;

 protected:
  void* sphdata_unsafe;

 public:
  using Hydrodynamics<ndim>::allocated;
  using Hydrodynamics<ndim>::eos;
  using Hydrodynamics<ndim>::h_fac;
  using Hydrodynamics<ndim>::hydrodata_unsafe;
  using Hydrodynamics<ndim>::invndim;
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
  using Hydrodynamics<ndim>::types;


  // Constructor
  //-----------------------------------------------------------------------------------------------
  Sph(int hydro_forces_aux, int self_gravity_aux, FLOAT alpha_visc_aux, FLOAT beta_visc_aux,
      FLOAT h_fac_aux, FLOAT h_converge_aux, aviscenum avisc_aux, acondenum acond_aux,
      tdaviscenum tdavisc_aux, string gas_eos_aux, string KernelName, int size_sph_part,
      SimUnits &units, Parameters *params);

  virtual ~Sph() {} ;

  virtual void AllocateMemory(int) = 0;
  virtual void DeallocateMemory(void) = 0;
  virtual void DeleteDeadParticles(void) = 0;
  virtual void ReorderParticles(void) = 0;


  // SPH functions for computing SPH sums with neighbouring particles
  // (fully coded in each separate SPH implementation, and not in Sph.cpp)
  //-----------------------------------------------------------------------------------------------
  virtual int ComputeH(const int, const int, const FLOAT, FLOAT *, FLOAT *, FLOAT *, FLOAT *,
                       SphParticle<ndim> &, Nbody<ndim> *) = 0;
  virtual void ComputeThermalProperties(SphParticle<ndim> &) = 0;
  virtual void ComputeSphHydroForces(const int, const int, const int *, const FLOAT *,
                                     const FLOAT *, const FLOAT *, SphParticle<ndim> &,
                                     SphParticle<ndim> *) = 0;
  virtual void ComputeSphHydroGravForces(const int, const int, int *, SphParticle<ndim> &,
                                         SphParticle<ndim> *) = 0;
  virtual void ComputeSphGravForces(const int, const int, int *,
                                    SphParticle<ndim> &, SphParticle<ndim> *) = 0;
  virtual void ComputeDirectGravForces(const int, const int, int *, SphParticle<ndim> &,
                                       SphParticle<ndim> *) = 0;
  virtual void ComputeSphNeibDudt(const int, const int, int *, FLOAT *, FLOAT *, FLOAT *,
                                  SphParticle<ndim> &, SphParticle<ndim> *) = 0;
  virtual void ComputeSphDerivatives(const int, const int, int *, FLOAT *, FLOAT *, FLOAT *,
                                     SphParticle<ndim> &, SphParticle<ndim> *) = 0;
  virtual void ComputeStarGravForces(const int, NbodyParticle<ndim> **, SphParticle<ndim> &) = 0;


  // SPH array memory allocation functions
  //-----------------------------------------------------------------------------------------------
  void InitialSmoothingLengthGuess(void);


  // Functions needed to hide some implementation details
  //-----------------------------------------------------------------------------------------------
  SphParticle<ndim>& GetSphParticlePointer(const int i) {
    long int numBytes = (long int) i * (long int) size_sph_part;
    return *((SphParticle<ndim>*)((unsigned char*) sphdata_unsafe + numBytes));
  };
  virtual SphParticle<ndim>* GetSphParticleArray() = 0;


  // Const variables (read in from parameters file)
  //-----------------------------------------------------------------------------------------------
  const acondenum acond;               ///< Artificial conductivity enum
  const aviscenum avisc;               ///< Artificial viscosity enum
  const tdaviscenum tdavisc;           ///< Time-dependent art. viscosity enum
  const FLOAT alpha_visc;              ///< alpha artificial viscosity parameter
  const FLOAT beta_visc;               ///< beta artificial viscosity parameter
  const FLOAT h_converge;              ///< h-rho iteration tolerance


  // SPH particle counters and main particle data array
  //-----------------------------------------------------------------------------------------------
  int create_sinks;                    ///< Create new sink particles?
  int fixed_sink_mass;                 ///< Fix masses of sink particles
  //int Ngather;                         ///< Average no. of gather neighbours
  FLOAT alpha_visc_min;                ///< Min. time-dependent viscosity alpha
  FLOAT msink_fixed;                   ///< Fixed sink mass value
  FLOAT hmin_sink;                     ///< Minimum smoothing length of sinks
  string riemann_solver;               ///< Selected Riemann solver
  string slope_limiter;                ///< Selected slope limiter

};



//=================================================================================================
//  Class GradhSph
/// \brief   Class definition for conservative 'grad-h' SPH simulations.
/// \details Class definition for conservative 'grad-h' SPH simulations (as derived from the parent
///          Sph class).  Full code for each of these class functions written in 'GradhSph.cpp'.
/// \author  D. A. Hubber, G. Rosotti
/// \date    03/04/2013
//=================================================================================================
#if !defined(SWIG)
template <int ndim, template<int> class kernelclass>
class GradhSph: public Sph<ndim>
{
public:
  using Sph<ndim>::allocated;
  using Sph<ndim>::Nhydro;
  using Sph<ndim>::Ntot;
  using Sph<ndim>::eos;
  using Sph<ndim>::fixed_sink_mass;
  using Sph<ndim>::hydrodata_unsafe;
  using Sph<ndim>::h_fac;
  using Sph<ndim>::kernp;
  using Sph<ndim>::kernfac;
  using Sph<ndim>::kernfacsqd;
  using Sph<ndim>::invndim;
  using Sph<ndim>::h_converge;
  using Sph<ndim>::hydro_forces;
  using Sph<ndim>::msink_fixed;
  using Sph<ndim>::self_gravity;
  using Sph<ndim>::avisc;
  using Sph<ndim>::beta_visc;
  using Sph<ndim>::alpha_visc;
  using Sph<ndim>::alpha_visc_min;
  using Sph<ndim>::acond;
  using Sph<ndim>::create_sinks;
  using Sph<ndim>::hmin_sink;
  using Sph<ndim>::Nhydromax;
  using Sph<ndim>::iorder;
  using Sph<ndim>::sphdata_unsafe;

 //public:

  GradhSph(int, int, FLOAT, FLOAT, FLOAT, FLOAT,
           aviscenum, acondenum, tdaviscenum, string, string, SimUnits &, Parameters *);
  virtual ~GradhSph();

  virtual Particle<ndim>* GetParticleArray() {return sphdata;};
  virtual SphParticle<ndim>* GetSphParticleArray() {return sphdata;};

  virtual void AllocateMemory(int);
  virtual void DeallocateMemory(void);
  virtual void DeleteDeadParticles(void);
  virtual void ReorderParticles(void);

  int ComputeH(const int, const int, const FLOAT, FLOAT *, FLOAT *, FLOAT *, FLOAT *,
               SphParticle<ndim> &, Nbody<ndim> *);
  void ComputeThermalProperties(SphParticle<ndim> &);
  void ComputeSphGravForces(const int, const int, int *, SphParticle<ndim> &, SphParticle<ndim> *);
  void ComputeSphHydroGravForces(const int, const int, int *,
                                 SphParticle<ndim> &, SphParticle<ndim> *);
  void ComputeSphHydroForces(const int, const int, const int *, const FLOAT *, const FLOAT *,
                             const FLOAT *, SphParticle<ndim> &, SphParticle<ndim> *);
  void ComputeSphNeibDudt(const int, const int, int *, FLOAT *, FLOAT *, FLOAT *,
                          SphParticle<ndim> &, SphParticle<ndim> *) {};
  void ComputeSphDerivatives(const int, const int, int *, FLOAT *, FLOAT *, FLOAT *,
                             SphParticle<ndim> &, SphParticle<ndim> *) {};
  void ComputeDirectGravForces(const int, const int, int *,
                               SphParticle<ndim> &, SphParticle<ndim> *);
  void ComputeStarGravForces(const int, NbodyParticle<ndim> **, SphParticle<ndim> &);

  kernelclass<ndim> kern;                  ///< SPH kernel
  GradhSphParticle<ndim> *sphdata;         ///< Pointer to particle data

};



//=================================================================================================
//  Class SM2012Sph
/// \brief   Class definition for Saitoh & Makino (2012) SPH simulations
/// \details Class definition for Saitoh & Makino (2012) SPH simulations (as derived from parent
///          Sph class).  Full code for each class functions is written in 'SM2012Sph.cpp'.
/// \author  D. A. Hubber, G. Rosotti
/// \date    03/04/2013
//=================================================================================================
template <int ndim, template<int> class kernelclass>
class SM2012Sph: public Sph<ndim>
{
  using Sph<ndim>::allocated;
  using Sph<ndim>::Nhydro;
  using Sph<ndim>::Ntot;
  using Sph<ndim>::eos;
  using Sph<ndim>::h_fac;
  using Sph<ndim>::kernfacsqd;
  using Sph<ndim>::invndim;
  using Sph<ndim>::h_converge;
  using Sph<ndim>::hydrodata_unsafe;
  using Sph<ndim>::hydro_forces;
  using Sph<ndim>::avisc;
  using Sph<ndim>::beta_visc;
  using Sph<ndim>::alpha_visc;
  using Sph<ndim>::alpha_visc_min;
  using Sph<ndim>::acond;
  using Sph<ndim>::create_sinks;
  using Sph<ndim>::hmin_sink;
  using Sph<ndim>::Nhydromax;
  using Sph<ndim>::kernp;
  using Sph<ndim>::iorder;
  using Sph<ndim>::sphdata_unsafe;

 public:

  SM2012Sph(int, int, FLOAT, FLOAT, FLOAT, FLOAT,
            aviscenum, acondenum, tdaviscenum, string, string, SimUnits &, Parameters *);
  virtual ~SM2012Sph();

  virtual Particle<ndim>* GetParticleArray() {return sphdata;};
  virtual SphParticle<ndim>* GetSphParticleArray() {return sphdata;};

  virtual void AllocateMemory(int);
  virtual void DeallocateMemory(void);
  virtual void DeleteDeadParticles(void);
  virtual void ReorderParticles(void);

  int ComputeH(const int, const int, const FLOAT, FLOAT *, FLOAT *, FLOAT *, FLOAT *,
               SphParticle<ndim> &, Nbody<ndim> *);
  void ComputeThermalProperties(SphParticle<ndim> &);
  void ComputeSphHydroForces(const int, const int, const int *, const FLOAT *, const FLOAT *,
                             const FLOAT *, SphParticle<ndim> &, SphParticle<ndim> *);
  void ComputeSphHydroGravForces(const int, const int, int *,
                                 SphParticle<ndim> &, SphParticle<ndim> *) {};
  void ComputeSphGravForces(const int, const int, int *,
                            SphParticle<ndim> &, SphParticle<ndim> *) {};
  void ComputeSphNeibDudt(const int, const int, int *, FLOAT *, FLOAT *, FLOAT *,
                          SphParticle<ndim> &, SphParticle<ndim> *) {};
  void ComputeSphDerivatives(const int, const int, int *, FLOAT *, FLOAT *, FLOAT *,
                             SphParticle<ndim> &, SphParticle<ndim> *) {};
  void ComputeDirectGravForces(const int, const int, int *,
                               SphParticle<ndim> &, SphParticle<ndim> *) {};
  void ComputeStarGravForces(const int, NbodyParticle<ndim> **, SphParticle<ndim> &) {};

  kernelclass<ndim> kern;                   ///< SPH kernel
  SM2012SphParticle<ndim> *sphdata;         ///< Pointer to particle data

};



//=================================================================================================
//  Class NullSph
/// \brief   Class definition for empty SPH class (needed in NbodySimulation).
/// \details Class definition for empty SPH class (needed in NbodySimulation).
/// \author  D. A. Hubber, G. Rosotti
/// \date    03/04/2013
//=================================================================================================
template <int ndim>
class NullSph: public Sph<ndim>
{
  using Sph<ndim>::allocated;
  using Sph<ndim>::Nhydro;
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
  using Sph<ndim>::create_sinks;
  using Sph<ndim>::hmin_sink;
  using Sph<ndim>::Nhydromax;
  using Sph<ndim>::kernp;
  using Sph<ndim>::iorder;
  using Sph<ndim>::sphdata_unsafe;

 public:

  NullSph(int hydro_forces_aux, int self_gravity_aux, FLOAT alpha_visc_aux,
          FLOAT beta_visc_aux, FLOAT h_fac_aux, FLOAT h_converge_aux,
          aviscenum avisc_aux, acondenum acond_aux, tdaviscenum tdavisc_aux,
          string gas_eos_aux, string KernelName, int size_sph_part, SimUnits &units, Parameters *):
    Sph<ndim>(hydro_forces_aux, self_gravity_aux, alpha_visc_aux,
              beta_visc_aux, h_fac_aux, h_converge_aux, avisc_aux, acond_aux,
              tdavisc_aux, gas_eos_aux, KernelName, size_sph_part, units),
   sphdata(NULL) {};

  virtual ~NullSph() ;

  virtual Particle<ndim>* GetParticleArray() {return sphdata;};
  virtual SphParticle<ndim>* GetSphParticleArray() {return sphdata;};

  virtual void AllocateMemory(int) {};
  virtual void DeallocateMemory(void) {};
  virtual void DeleteDeadParticles(void) {};
  virtual void ReorderParticles(void) {};

  int ComputeH(const int, const int, const FLOAT, FLOAT *, FLOAT *, FLOAT *, FLOAT *,
               SphParticle<ndim> &, Nbody<ndim> *) {return -1;}
  void ComputeThermalProperties(SphParticle<ndim> &) {};
  void ComputeSphHydroForces(const int, const int, const int *, const FLOAT *, const FLOAT *,
                             const FLOAT *, SphParticle<ndim> &, SphParticle<ndim> *) {};
  void ComputeSphHydroGravForces(const int, const int, int *,
                                 SphParticle<ndim> &, SphParticle<ndim> *) {};
  void ComputeSphGravForces(const int, const int, int *,
                            SphParticle<ndim> &, SphParticle<ndim> *) {};
  void ComputeSphNeibDudt(const int, const int, int *, FLOAT *, FLOAT *, FLOAT *,
                          SphParticle<ndim> &, SphParticle<ndim> *) {};
  void ComputeSphDerivatives(const int, const int, int *, FLOAT *, FLOAT *, FLOAT *,
                             SphParticle<ndim> &, SphParticle<ndim> *) {};
  void ComputeDirectGravForces(const int, const int, int *,
                               SphParticle<ndim> &, SphParticle<ndim> *) {};
  void ComputeStarGravForces(int, NbodyParticle<ndim> **, SphParticle<ndim> &) {};

  //kernelclass<ndim> kern;               ///< SPH kernel
  SphParticle<ndim> *sphdata;           ///< Pointer to particle data

};
#endif


#endif
