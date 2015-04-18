//=================================================================================================
//  FV.h
//  Contains main parent virtual class plus child classes for various FV
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


#ifndef _FV_H_
#define _FV_H_


#include <assert.h>
#include <string>
#include "Precision.h"
#include "Constants.h"
#include "Hydrodynamics.h"
#include "Particle.h"
#include "SmoothingKernel.h"
#include "NbodyParticle.h"
#include "Nbody.h"
#include "Parameters.h"
#include "DomainBox.h"
#include "EOS.h"
#include "RiemannSolver.h"
#include "SlopeLimiter.h"
#include "ExternalPotential.h"
#if defined _OPENMP
#include "omp.h"
#endif
using namespace std;


//=================================================================================================
//  Class FV
/// \brief   Main parent FV class.
/// \details
/// \author  D. A. Hubber, J. Ngoumou
/// \date    19/02/2015
//=================================================================================================
//template <int ndim>
template <int ndim>
class FV : public Hydrodynamics<ndim>
{
public:

  using Hydrodynamics<ndim>::allocated;
  using Hydrodynamics<ndim>::eos;
  using Hydrodynamics<ndim>::h_fac;
  using Hydrodynamics<ndim>::hydrodata_unsafe;
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
  using Hydrodynamics<ndim>::size_hydro_part;

  static const FLOAT invndim=1./ndim;  ///< Copy of 1/ndim
  static const int nvar = ndim + 2;
  static const int ivx = 0;
  static const int ivy = 1;
  static const int ivz = 2;
  static const int irho = ndim;
  static const int ietot = ndim + 1;
  static const int ipress = ndim + 1;


  // Constructor
  //-----------------------------------------------------------------------------------------------
  FV(int hydro_forces_aux, int self_gravity_aux, FLOAT _accel_mult, FLOAT _courant_mult,
     FLOAT h_fac_aux, FLOAT h_converge_aux,
     FLOAT gamma_aux, string gas_eos_aux, string KernelName, int size_mfv_part);
  ~FV();


  virtual void AllocateMemory(int) = 0;
  virtual void DeallocateMemory(void) = 0;
  virtual void DeleteDeadParticles(void) = 0;
  virtual void ReorderParticles(void) = 0;

  void CalculateFluxVectorFromPrimitive(FLOAT Wprim[nvar], FLOAT flux[nvar][ndim]);
  void CalculatePrimitiveTimeDerivative(FLOAT Wprim[nvar], FLOAT gradW[nvar][ndim], FLOAT Wdot[nvar]);
  void ConvertQToConserved(const FLOAT, const FLOAT Qcons[nvar], FLOAT Ucons[nvar]);
  void ConvertConservedToQ(const FLOAT, const FLOAT Ucons[nvar], FLOAT Qcons[nvar]);
  void ConvertConservedToPrimitive(const FLOAT Ucons[nvar], FLOAT Wprim[nvar]);
  void ConvertPrimitiveToConserved(const FLOAT Wprim[nvar], FLOAT Ucons[nvar]);
  //void CalculateConservedFluxFromConserved(int k, FLOAT Ucons[nvar], FLOAT flux[nvar]);


  // Functions needed to hide some implementation details
  //-----------------------------------------------------------------------------------------------
  /*FVParticle<ndim>& GetFVParticlePointer(const int i) {
    return *((FVParticle<ndim>*)((unsigned char*) hydrodata_unsafe + i*size_hydro_part));
  }
  FVParticle<ndim>* GetFVParticleArray() {return hydrodata;}
  virtual Particle<ndim>* GetParticleArray() {return hydrodata;};*/


};
#endif
