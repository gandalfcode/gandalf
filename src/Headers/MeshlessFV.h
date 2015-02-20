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
#include "ExternalPotential.h"
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
template <int ndim, template<int> class kernelclass>
class MeshlessFV : public Hydrodynamics<ndim>
{
 public:

  using Hydrodynamics<ndim>::allocated;
  using Hydrodynamics<ndim>::h_fac;
  using Hydrodynamics<ndim>::hydrodata_unsafe;
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
  using Hydrodynamics<ndim>::size_hydro_part;

  // Constructor
  //-----------------------------------------------------------------------------------------------
  MeshlessFV(int hydro_forces_aux, int self_gravity_aux, FLOAT h_fac_aux, FLOAT h_converge_aux,
             string gas_eos_aux, string KernelName, int size_MeshlessFV_part);


  virtual void AllocateMemory(int) = 0;
  virtual void DeallocateMemory(void) = 0;
  virtual void DeleteDeadParticles(void) = 0;
  virtual void ReorderParticles(void) = 0;


  // MeshlessFV functions for computing MeshlessFV sums with neighbouring particles
  // (fully coded in each separate MeshlessFV implementation, and not in MeshlessFV.cpp)
  //-----------------------------------------------------------------------------------------------
  int ComputeH(const int, const int, const FLOAT, FLOAT *, FLOAT *, FLOAT *, FLOAT *,
               MeshlessFVParticle<ndim> &, Nbody<ndim> *) = 0;
  void ComputeThermalProperties(MeshlessFVParticle<ndim> &) = 0;
  void ComputeNormalisedDerivatives(const int, const int, int *, FLOAT *, FLOAT *, FLOAT *,
                                    MeshlessFVParticle<ndim> &, MeshlessFVParticle<ndim> *);
  //virtual void ComputeDirectGravForces(const int, const int, int *, MeshlessFVParticle<ndim> &,
  //                                     MeshlessFVParticle<ndim> *) = 0;
  //virtual void ComputeStarGravForces(const int, NbodyParticle<ndim> **, MeshlessFVParticle<ndim> &) = 0;


  // MeshlessFV array memory allocation functions
  //-----------------------------------------------------------------------------------------------
  void InitialSmoothingLengthGuess(void);


  // Functions needed to hide some implementation details
  //-----------------------------------------------------------------------------------------------
  MeshlessFVParticle<ndim>& GetMeshlessFVParticlePointer(const int i) {
    return *((MeshlessFVParticle<ndim>*)((unsigned char*)hydrodata_unsafe + i*size_hydro_part));
  }
  MeshlessFVParticle<ndim>* GetMeshlessFVParticleArray() {return hydrodata;}


  // Const variables (read in from parameters file)
  //-----------------------------------------------------------------------------------------------
  const FLOAT h_converge;              ///< h-rho iteration tolerance
  static const FLOAT invndim=1./ndim;  ///< Copy of 1/ndim


  // MeshlessFV particle counters and main particle data array
  //-----------------------------------------------------------------------------------------------
  int create_sinks;                    ///< Create new sink particles?
  int fixed_sink_mass;                 ///< Fix masses of sink particles
  int riemann_order;                   ///< Order of Riemann solver
  FLOAT msink_fixed;                   ///< Fixed sink mass value
  FLOAT hmin_sink;                     ///< Minimum smoothing length of sinks
  string riemann_solver;               ///< Selected Riemann solver
  string slope_limiter;                ///< Selected slope limiter
  //MeshlessFVType MeshlessFVtype[Nhydrotypes];        ///< Array of MeshlessFV types

  MeshlessFVParticle<ndim> *hydrodata;

};
#endif
