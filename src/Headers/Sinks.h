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
#include "StarParticle.h"
#include "SystemParticle.h"
#if defined MPI_PARALLEL
#include "MpiControl.h"
#endif
using namespace std;



//=================================================================================================
//  Class SinkParticle
/// \brief   Individual sink particle data structure
/// \details Main parent N-body particle data structure.  All main other N-body particle types
///          (e.g. stars, systems, sinks) are derived from this class.
/// \author  D. A. Hubber
/// \date    15/04/2013
//=================================================================================================
template <int ndim>
class SinkParticle
{
 public:

  StarParticle<ndim> *star;            ///< Pointer to connected star particle
  int istar;                           ///< i.d. of connected star particle
  int Ngas;                            ///< No. of gas particles inside sink
  FLOAT macctot;                      ///< Total accreted mass
  FLOAT radius;                       ///< Softening/sink radius of particle
  FLOAT dt;                           ///< Particle timestep
  FLOAT dmdt;                         ///< Accretion rate
  FLOAT menc;                         ///< Gas mass enclosed within sink radius
  FLOAT mmax;                         ///< Max. mass before increasing dmdt
  FLOAT racc;                         ///< Accretion radius
  FLOAT ketot;                        ///< Internal kinetic energy
  FLOAT gpetot;                       ///< Internal grav. pot. energy
  FLOAT rotketot;                     ///< Internal rotational kinetic energy
  FLOAT utot;                         ///< Internal energy accrted by sink
  FLOAT taccrete;                     ///< Accretion timescale
  FLOAT trad;                         ///< Radial accretion timescale
  FLOAT trot;                         ///< Rotational period at sink radius
  FLOAT tvisc;                        ///< Viscous accretion timescale
  FLOAT angmom[3];                    ///< Internal sink angular momentum


  // Star particle constructor to initialise all values
  //-----------------------------------------------------------------------------------------------
  SinkParticle()
  {
    star     = 0;
    istar    = -1;
    Ngas     = 0;
    macctot  = 0.0;
    radius   = 0.0;
    dt       = 0.0;
    dmdt     = 0.0;
    menc     = 0.0;
    mmax     = 0.0;
    racc     = 0.0;
    ketot    = 0.0;
    gpetot   = 0.0;
    rotketot = 0.0;
    utot     = 0.0;
    taccrete = 0.0;
    trad     = 0.0;
    trot     = 0.0;
    tvisc    = 0.0;
    for (int k=0; k<3; k++) angmom[k] = 0.0;
  }

};



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
  int Nsinkfixed;                      ///< Fixed maximum no. of sinks for testing purposes
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
