//=================================================================================================
//  Hydrodynamics.h
//  Contains main parent virtual class plus child classes for various hydrodynamics
//  algorithms that are implemented (e.g. SPH, Meshless Finite-Volume)
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


#ifndef _HYDRODYNAMICS_H_
#define _HYDRODYNAMICS_H_


#include <assert.h>
#include <string>
#include "Precision.h"
#include "Constants.h"
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

template <int ndim>
class Hydrodynamics;

template <int ndim>
class EOS;

static const FLOAT ghost_range = 1.6;


//=================================================================================================
//  Class Hydrodynamics
/// \brief   Main parent hydrodynamics class
/// \details Different hydrodynamics implementations (e.g. SPH, Meshless Finite-Volume)
///          are derived from this class.  Each implementation requires defining its own version
///          of each function (e.g. for computing hydrodynamical fluxes/forces).
/// \author  D. A. Hubber, G. Rosotti
/// \date    02/02/2015
//=================================================================================================
template <int ndim>
class Hydrodynamics
{
protected:
  const int size_hydro_part;
  void* hydrodata_unsafe;

public:

  // Constructor
  //-----------------------------------------------------------------------------------------------
  Hydrodynamics(int hydro_forces_aux, int self_gravity_aux, FLOAT h_fac_aux,
                string gas_eos_aux, string KernelName, int size_hydro_part);


  // SPH array memory allocation functions
  //-----------------------------------------------------------------------------------------------
  virtual void AllocateMemory(int) = 0;
  virtual void DeallocateMemory(void) = 0;
  virtual void DeleteDeadParticles(void) = 0;
  virtual void ReorderParticles(void) = 0;
  void ComputeBoundingBox(FLOAT *, FLOAT *, const int);
  void CheckXBoundaryGhostParticle(const int, const FLOAT, const DomainBox<ndim> &);
  void CheckYBoundaryGhostParticle(const int, const FLOAT, const DomainBox<ndim> &);
  void CheckZBoundaryGhostParticle(const int, const FLOAT, const DomainBox<ndim> &);
  void CreateBoundaryGhostParticle(const int, const int, const int, const FLOAT, const FLOAT);


  // Functions needed to hide some implementation details
  //-----------------------------------------------------------------------------------------------
  Particle<ndim>& GetParticlePointer(const int i) {
    return *((Particle<ndim>*)((unsigned char*) hydrodata_unsafe + i*size_hydro_part));
  };
  virtual Particle<ndim>* GetParticleArray() = 0;


  // Const variables (read in from parameters file)
  //-----------------------------------------------------------------------------------------------
  const int hydro_forces;              ///< Compute hydro forces?
  const int self_gravity;              ///< Compute gravitational forces?
  const string gas_eos;                ///< Gas EOS option
  const FLOAT h_fac;                   ///< Smoothing length-density factor
  static const FLOAT invndim;          ///< Copy of 1/ndim


  // SPH particle counters and main particle data array
  //-----------------------------------------------------------------------------------------------
  bool allocated;                      ///< Is memory allocated?
  int Ngather;                         ///< No. of gather neighbours
  int Nghost;                          ///< No. of ghost particles
  int Nghostmax;                       ///< Max. allowed no. of ghost particles
  int NImportedParticles;              ///< No. of imported particles
                                       ///< (to compute forces on behalf of other processors)
  int Nmpighost;                       ///< No. of MPI ghost particles
  int NPeriodicGhost;                  ///< No. of periodic ghost particles
  int Nhydro;                          ///< No. of hydro particles in simulation
  int Nhydromax;                       ///< Max. no. of hydro particles in array
  int Ntot;                            ///< No. of real + ghost particles
  FLOAT hmin_sink;                     ///< Minimum smoothing length for sinks
  FLOAT kernfac;                       ///< Kernel range neighbour fraction
  FLOAT kernfacsqd;                    ///< Kernel range neib. fraction squared
  FLOAT kernrange;                     ///< Kernel range
  FLOAT mmean;                         ///< Mean SPH particle mass


  int *iorder;                         ///< Array containing particle ordering
  EOS<ndim> *eos;                      ///< Equation-of-state
  SmoothingKernel<ndim> *kernp;        ///< Pointer to chosen kernel object
  TabulatedKernel<ndim> kerntab;       ///< Tabulated version of chosen kernel
  RiemannSolver<ndim> *riemann;        ///< Riemann solver
  ExternalPotential<ndim> *extpot;     ///< Pointer to external potential object

};


// Declare invndim constant here (prevents warnings with some compilers)
template <int ndim>
const FLOAT Hydrodynamics<ndim>::invndim = 1.0/ndim;

#endif
