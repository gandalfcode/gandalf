//=================================================================================================
//  Sinks.h
//  Main sink particle class
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


#ifndef _SINKS_H_
#define _SINKS_H_


#include <string>
#include "Precision.h"
#include "CodeTiming.h"
#include "Constants.h"
#include "Hydrodynamics.h"
#include "Nbody.h"
#include "NbodyParticle.h"
#include "NeighbourSearch.h"
#include "Parameters.h"
#include "Particle.h"
#include "SmoothingKernel.h"
#include "SinkParticle.h"
#include "SystemParticle.h"
#if defined MPI_PARALLEL
#include "MpiControl.h"
#endif
using namespace std;




//=================================================================================================
//  Class Sinks
/// \brief   Main sink particle class.
/// \details Main sink particle class for searching for and creating new sinks, and for
///          controlling the accretion of SPH particles onto sink particles.
/// \author  D. A. Hubber
/// \date    08/06/2013
//=================================================================================================
template<int ndim>
class Sinks
{
#if defined MPI_PARALLEL
 MpiControl<ndim>* mpicontrol;
#endif
 public:

  // Constructor and destructor
  //-----------------------------------------------------------------------------------------------
  Sinks(NeighbourSearch<ndim> *);
  ~Sinks();

#if defined MPI_PARALLEL
  void SetMpiControl(MpiControl<ndim>* mpicontrol_aux) {mpicontrol=mpicontrol_aux;};
#endif

  // Function prototypes
  //-----------------------------------------------------------------------------------------------
  void AllocateMemory(int);
  void DeallocateMemory(void);
  void SearchForNewSinkParticles(const int, const FLOAT,
                                 Hydrodynamics<ndim> *, Nbody<ndim> *);
  void CreateNewSinkParticle(const int, const FLOAT, Particle<ndim>&,
                             Hydrodynamics<ndim> *, Nbody<ndim> *);
  void AccreteMassToSinks(const int, const FLOAT,
                          Hydrodynamics<ndim> *, Nbody<ndim> *);


  // Local class variables
  //-----------------------------------------------------------------------------------------------
  bool allocated_memory;               ///< Has sink memory been allocated?
  int Nsink;                           ///< No. of sink particles
  int Nsinkmax;                        ///< Max. no. of sink particles
  int sink_particles;                  ///< Using sink particles?
  int create_sinks;                    ///< Create new sink particles?
  int smooth_accretion;                ///< Use smooth accretion?
  FLOAT alpha_ss;                      ///< Shakura-Sunyaev alpha viscosity
  FLOAT rho_sink;                      ///< Sink formation density
  FLOAT sink_radius;                   ///< New sink radius (in units of h)
  FLOAT smooth_accrete_frac;           ///< Minimum particle mass fraction
  FLOAT smooth_accrete_dt;             ///< Minimum particle timestep fraction
  string sink_radius_mode;             ///< Sink radius mode

  SinkParticle<ndim> *sink;            ///< Main sink particle array
  CodeTiming *timing;                  ///< Pointer to code timing object
  NeighbourSearch<ndim> *neibsearch;   ///< Pointer to neighbour search object

};
#endif
