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
#include "NeighbourManager.h"
#include "Parameters.h"
#include "Precision.h"
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
  using Hydrodynamics<ndim>::create_sinks;
  using Hydrodynamics<ndim>::eos;
  using Hydrodynamics<ndim>::h_fac;
  using Hydrodynamics<ndim>::hmin_sink;
  using Hydrodynamics<ndim>::hydrodata_unsafe;
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
  using Hydrodynamics<ndim>::types;

  typedef typename SphParticle<ndim>::DensityParticle DensityParticle;

  // Constructor
  //-----------------------------------------------------------------------------------------------
  Sph(int hydro_forces_aux, int self_gravity_aux, FLOAT alpha_visc_aux, FLOAT beta_visc_aux,
      FLOAT h_fac_aux, FLOAT h_converge_aux, aviscenum avisc_aux, acondenum acond_aux,
      tdaviscenum tdavisc_aux, string gas_eos_aux, string KernelName, int size_sph_part,
      SimUnits &units, Parameters *params);

  virtual ~Sph() {} ;

  virtual void AllocateMemory(int) = 0;
  virtual void DeallocateMemory(void) = 0;
  virtual int DeleteDeadParticles(void) = 0;
  virtual void AccreteMassFromParticle(const FLOAT dm, Particle<ndim> &part) {part.m -= dm;}

  virtual void ZeroAccelerations() ;



  // SPH functions for computing SPH sums with neighbouring particles
  // (fully coded in each separate SPH implementation, and not in Sph.cpp)
  //-----------------------------------------------------------------------------------------------
  virtual int ComputeH(SphParticle<ndim> &, FLOAT, const vector<DensityParticle> &,
                       Nbody<ndim> *) = 0;
  virtual void ComputeThermalProperties(SphParticle<ndim> &) = 0;
  virtual void ComputeStarGravForces(const int, NbodyParticle<ndim> **, SphParticle<ndim> &) = 0;

#if !defined(SWIG)
  template<template<int> class Kernel>
  void ComputeCullenAndDehnenViscosity(SphParticle<ndim> &, const vector<DensityParticle> &,
                                       Kernel<ndim>&);
#endif

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
  int conservative_sph_star_gravity;   ///< Use Hubber et al. (2013) conservative SPH-star gravity
  int fixed_sink_mass;                 ///< Fix masses of sink particles
  FLOAT alpha_visc_min;                ///< Min. time-dependent viscosity alpha
  FLOAT msink_fixed;                   ///< Fixed sink mass value
  string riemann_solver;               ///< Selected Riemann solver
  string slope_limiter;                ///< Selected slope limiter

};

#if !defined(SWIG)

template <int ndim>
class GradhSphBase: public Sph<ndim>
{
public:
  typedef GradhSphParticle<ndim> ParticleType ;
  typedef typename ParticleType::HydroForcesParticle HydroNeib ;
  typedef typename GravityNeighbourLists<HydroNeib>::DirectType DirectNeib ;


  GradhSphBase(int hydro_forces_aux, int self_gravity_aux,
  FLOAT alpha_visc_aux, FLOAT beta_visc_aux, FLOAT h_fac_aux, FLOAT h_converge_aux,
  aviscenum avisc_aux, acondenum acond_aux, tdaviscenum tdavisc_aux,
  string gas_eos_aux, string KernelName, SimUnits &units, Parameters *params):
    Sph<ndim>(hydro_forces_aux, self_gravity_aux, alpha_visc_aux, beta_visc_aux,
              h_fac_aux, h_converge_aux, avisc_aux, acond_aux, tdavisc_aux,
              gas_eos_aux, KernelName, sizeof(GradhSphParticle<ndim>), units, params) {};


  virtual void ComputeSphGravForces(GradhSphParticle<ndim>&, NeighbourList<HydroNeib>&) = 0;
  virtual void ComputeSphHydroGravForces(GradhSphParticle<ndim>&, NeighbourList<HydroNeib>&) = 0;
  virtual void ComputeSphHydroForces(GradhSphParticle<ndim>&, NeighbourList<HydroNeib>&) = 0;
  void ComputeDirectGravForces(GradhSphParticle<ndim>&, NeighbourList<DirectNeib>&) ;


};

//=================================================================================================
//  Class GradhSph
/// \brief   Class definition for conservative 'grad-h' SPH simulations.
/// \details Class definition for conservative 'grad-h' SPH simulations (as derived from the parent
///          Sph class).  Full code for each of these class functions written in 'GradhSph.cpp'.
/// \author  D. A. Hubber, G. Rosotti
/// \date    03/04/2013
//=================================================================================================
template <int ndim, template<int> class kernelclass>
class GradhSph: public GradhSphBase<ndim>
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
  using Sph<ndim>::sphdata_unsafe;
  using Sph<ndim>::tdavisc;

  using Hydrodynamics<ndim>::rho_sink;
  using Hydrodynamics<ndim>::sink_particles;

  typedef typename GradhSphBase<ndim>::HydroNeib HydroNeib;
  typedef typename GradhSphBase<ndim>::DirectNeib  DirectNeib;
  typedef typename SphParticle<ndim>::DensityParticle DensityParticle;


 //public:

  GradhSph(int, int, FLOAT, FLOAT, FLOAT, FLOAT,
           aviscenum, acondenum, tdaviscenum, string, string, SimUnits &, Parameters *);
  virtual ~GradhSph();

  virtual SphParticle<ndim>* GetSphParticleArray() {return sphdata;};

  virtual void AllocateMemory(int);
  virtual void DeallocateMemory(void);
  virtual int DeleteDeadParticles(void) {
    return this->template DoDeleteDeadParticles<GradhSphParticle>() ;
  }

  virtual int ComputeH(SphParticle<ndim> &, FLOAT, const vector<DensityParticle> &, Nbody<ndim> *);
  void ComputeThermalProperties(SphParticle<ndim> &);
  virtual void ComputeSphGravForces(GradhSphParticle<ndim>&, NeighbourList<HydroNeib>&);
  virtual void ComputeSphHydroGravForces(GradhSphParticle<ndim>&, NeighbourList<HydroNeib>&);
  virtual void ComputeSphHydroForces(GradhSphParticle<ndim>&, NeighbourList<HydroNeib>&);
  void ComputeStarGravForces(const int, NbodyParticle<ndim> **, SphParticle<ndim> &);
#if defined MPI_PARALLEL
  virtual void FinishReturnExport ();
#endif

  inline FLOAT h_rho_func(const FLOAT m, const FLOAT rho) const
  {
    return h_fac*pow(m/rho, invndim);
    //return h_fac*pow((FLOAT) 0.5*m*((FLOAT) 1.0/rho + (FLOAT) 1.0/rho_sink), invndim);
  }
  inline FLOAT h_rho_deriv(const FLOAT h, const FLOAT m, const FLOAT rho) const
  {
    return -invndim*h/rho;
    //return -h*invndim*pow(rho, -(FLOAT) 2.0)/((FLOAT) 1.0/rho + (FLOAT) 1.0/rho_sink);
  }


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
  using Sph<ndim>::sphdata_unsafe;

 public:
  typedef typename SphParticle<ndim>::DensityParticle DensityParticle;

  SM2012Sph(int, int, FLOAT, FLOAT, FLOAT, FLOAT,
            aviscenum, acondenum, tdaviscenum, string, string, SimUnits &, Parameters *);
  virtual ~SM2012Sph();

  virtual SphParticle<ndim>* GetSphParticleArray() {return sphdata;};

  virtual void AllocateMemory(int);
  virtual void DeallocateMemory(void);
  virtual int DeleteDeadParticles(void) {
    return this->template DoDeleteDeadParticles<SM2012SphParticle>() ;
  }
  virtual int ComputeH(SphParticle<ndim> &, FLOAT, const vector<DensityParticle> &, Nbody<ndim> *);
  void ComputeThermalProperties(SphParticle<ndim> &);
  void ComputeSphHydroForces(const int, const int, const int *, const FLOAT *, const FLOAT *,
                             const FLOAT *, SphParticle<ndim> &, SphParticle<ndim>* );
  void ComputeSphHydroGravForces(const int, const int, int *,
                                 SphParticle<ndim> &, SphParticle<ndim> *) {};
  void ComputeSphGravForces(const int, const int, int *,
                            SphParticle<ndim> &, SphParticle<ndim> *) {};
  void ComputeDirectGravForces(const int, const int, int *,
                               SphParticle<ndim> &, SphParticle<ndim> *) {};
  void ComputeStarGravForces(const int, NbodyParticle<ndim> **, SphParticle<ndim> &) {};

  kernelclass<ndim> kern;                   ///< SPH kernel
  SM2012SphParticle<ndim> *sphdata;         ///< Pointer to particle data

};


//=================================================================================================
//  CurlVelSqd
/// Square of the curl in 1,2 and 3 dimensions
//=================================================================================================
inline FLOAT CurlVelSqd(FLOAT gradv[1][1]) {
  return 0;
}
inline FLOAT CurlVelSqd(FLOAT gradv[2][2]) {
  FLOAT curl = (gradv[1][0] - gradv[0][1]) ;
  return curl*curl ;
}
inline FLOAT CurlVelSqd(FLOAT gradv[3][3]) {
  FLOAT curl[3] = {
      gradv[1][2] - gradv[2][1],
      gradv[2][0] - gradv[0][2],
      gradv[0][1] - gradv[1][0] } ;

  return DotProduct(curl, curl, 3) ;
}

//=================================================================================================
//  Sph::ComputeCullenAndDehnenViscosity
/// Compute the viscosity limiter for the Cullen & Dehnen switch
//=================================================================================================
template <int ndim>
template<template<int> class Kernel>
void Sph<ndim>::ComputeCullenAndDehnenViscosity
(SphParticle<ndim> & parti,                         ///< [inout] Particle to compute the switch for
const vector<DensityParticle> &ngbs,                ///< [in] List of neighbours
Kernel<ndim>& kern)                                 ///< [in] Kernel
{

  // Buffers for gradients of vectors
  FLOAT dv[ndim][ndim];
  FLOAT da[ndim][ndim];
  FLOAT rr[ndim][ndim];
  FLOAT dvdx[ndim][ndim];
  FLOAT dadx[ndim][ndim];

  for (int i=0; i<ndim; ++i)
    for (int j=0; j<ndim; ++j) {
      rr[i][j] = da[i][j] = dv[i][j] = dadx[i][j] = dvdx[i][j] = 0 ;
    }

  FLOAT invh = 1 / parti.h;
  FLOAT hfac = invh * parti.hfactor / parti.rho ;
  int Nneib = ngbs.size() ;
  for (int i=0; i < Nneib; ++i) {

    FLOAT dr[ndim];
    for (int j=0; j < ndim; j++) dr[j] = ngbs[i].r[j] - parti.r[j] ;
    FLOAT w = ngbs[i].m * hfac * kern.w1(invh * sqrt(DotProduct(dr, dr, ndim)));

    for (int j=0; j < ndim; j++)
      for (int k=0; k < ndim; k++) {
        rr[j][k] += w * dr[j] * dr[k] ;
        dv[j][k] += w * dr[j] * (ngbs[i].v[k] - parti.v[k]) ;
        da[j][k] += w * dr[j] * (ngbs[i].a[k] - parti.a[k]) ;
      }
  }

  // Invert the rr matrix and compute the gradients
  FLOAT T[ndim][ndim] ;
  InvertMatrix(rr, T) ;

  // Check the accuracy of the integral gradients (using the square of the condition number),
  // if it's bad, we'll set alpha_loc to alpha_max.
  double modR(0), modT(0) ;
  for (int i=0; i<ndim; i++)
    for (int j=0; j<ndim; j++){
      modR += rr[i][j]*rr[i][j];
      modT +=  T[i][j]* T[i][j];
    }
  double sqd_condition_number = modR*modT / (ndim*ndim) ;

  FLOAT alpha_loc = 0;
  if (sqd_condition_number > 1e4) {
    // Bad gradients
    alpha_loc = alpha_visc ;
  } else {
    // Ok gradients
    for (int i=0; i<ndim; i++)
      for (int j=0; j<ndim; j++)
        for (int k=0; k<ndim; k++) {
          dvdx[i][j] += T[j][k] * dv[k][i];
          dadx[i][j] += T[j][k] * da[k][i];
        }

    // Now compute the components needed for the limiter:
    FLOAT ddivdt = 0;
    FLOAT divv2 = 0;
    for (int i=0; i<ndim; ++i) {
      ddivdt += dadx[i][i] ;
      for (int j=0; j<ndim; ++j)
        ddivdt -= dvdx[i][j]*dvdx[j][i];

      divv2 += dvdx[i][i] ;
    }

    divv2 *= divv2;
    FLOAT curlv2 = CurlVelSqd(dvdx) ;

    FLOAT f_balsara = 1 ;
    if (curlv2 > 0)
      f_balsara = (divv2 / (divv2 + curlv2)) ;

    if (ddivdt < 0) {
      alpha_loc = (10 * parti.h*parti.h / (parti.sound*parti.sound)) * f_balsara * (-ddivdt) ;
      alpha_loc = min(alpha_loc, alpha_visc) ;
    }
  }

  if (alpha_loc > parti.alpha)
    parti.alpha = alpha_loc ;

  parti.dalphadt = (FLOAT)0.1*parti.sound*(max(alpha_visc_min, alpha_loc) - parti.alpha)*invh;
}



#endif


#endif
