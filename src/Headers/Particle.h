//=================================================================================================
//  Particle.h
//  Main particle data structures
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


#ifndef _PARTICLE_H_
#define _PARTICLE_H_


#include "Precision.h"
#include "Constants.h"
#ifdef MPI_PARALLEL
#include <stddef.h>
#include "mpi.h"
#include "Exception.h"
#endif


enum ptype{gas, icm, boundary, cdm, dust, dead,
           x_lhs_periodic, x_lhs_mirror, x_rhs_periodic, x_rhs_mirror,
           y_lhs_periodic, y_lhs_mirror, y_rhs_periodic, y_rhs_mirror,
           z_lhs_periodic, z_lhs_mirror, z_rhs_periodic, z_rhs_mirror,
           Nhydrotypes};

enum parttype{gas_type, icm_type, cdm_type, dust_type, Ntypes} ;

static const int parttype_converter[] = { gas, icm, cdm, dust } ;
static const int parttype_reverse_converter[] =
{ gas_type, icm_type, gas_type, cdm_type, dust_type} ;

typedef bool Typemask[Ntypes];

//static Typemask truemask = {true};


//=================================================================================================
//  Structure ParticleType
/// \brief  Structure containing particle type information
/// \author D. A. Hubber
/// \date   14/10/2015
//=================================================================================================
struct ParticleTypeInfo
{
  int N;                               ///< Current no. of particles
  bool hydro_forces;                   ///< Does particle experience hydro forces?
  bool self_gravity;                   ///< Does particle experience gravitational forces?
  bool drag_forces ;                   ///< Does particle experience drag forces?
  Typemask hmask;                      ///< Neighbour mask for computing smoothing lengths
  Typemask hydromask;                  ///< Neighbour mask for computing hydro forces
  Typemask gravmask;                   ///< Neighbour mask for computing gravity forces
  Typemask dragmask;                   ///< Neighbour mask for computing drag forces

  ParticleTypeInfo() {
    N = 0;
    hydro_forces = false;
    self_gravity = false;
    drag_forces  = false ;
    for (int k=0; k<Ntypes; k++) hmask[k] = false;
    for (int k=0; k<Ntypes; k++) hydromask[k] = false;
    for (int k=0; k<Ntypes; k++) gravmask[k] = false;
    for (int k=0; k<Ntypes; k++) dragmask[k] = false;
  }
};



//=================================================================================================
//  Structure Particle
/// \brief  Main base particle data structure.
/// \author D. A. Hubber, G. Rosotti
/// \date   01/10/2013
//=================================================================================================
template <int ndim>
struct Particle
{
  bool active;                      ///< Flag if active (i.e. recompute step)
  bool potmin;                      ///< Is particle at a potential minima?
  int iorig;                        ///< Original particle i.d.
  int itype;                        ///< SPH particle type (eg boundary/dead)
  int ptype;                        ///< SPH particle type (gas/cdm/dust)
  int sinkid;                       ///< i.d. of sink particle
  int levelneib;                    ///< Min. timestep level of neighbours
  int nstep;                        ///< Integer step-size of particle
  int nlast;                        ///< Integer time at beginning of step
  int level;                        ///< Current timestep level of particle
  FLOAT r[ndim];                    ///< Position
  FLOAT v[ndim];                    ///< Velocity
  FLOAT a[ndim];                    ///< Total acceleration
  FLOAT r0[ndim];                   ///< Position at beginning of step
  FLOAT v0[ndim];                   ///< Velocity at beginning of step
  FLOAT a0[ndim];                   ///< Acceleration at beginning of step
  FLOAT agrav[ndim];                ///< Gravitational acceleration
  FLOAT a_dust[ndim];                ///< Gravitational acceleration
  FLOAT m;                          ///< Particle mass
  FLOAT h;                          ///< SPH smoothing length
  FLOAT h_dust ;                    ///< Gas Smoothing length for dust
  FLOAT hrangesqd;                  ///< Kernel extent (squared)
  FLOAT invh;                       ///< 1 / h
  FLOAT hfactor;                    ///< invh^(ndim + 1)
  FLOAT invrho;                     ///< 1 / rho
  FLOAT sound;                      ///< Sound speed
  FLOAT rho;                        ///< Density
  FLOAT u;                          ///< Specific internal energy
  FLOAT u0;                         ///< u at beginning of step
  FLOAT dudt0;                      ///< dudt at beginning of step
  FLOAT dudt;                       ///< Compressional heating rate
  FLOAT gpot;                       ///< Gravitational potential
  DOUBLE dt;                        ///< Particle timestep
  DOUBLE tlast;                     ///< Time at beginning of current step
  FLOAT ionfrac;                    ///< Ionisation fraction
  FLOAT Xion;                       ///< Ionisation fraciton (from tree)
  FLOAT mu_bar;                     ///< mean molecular weight
  FLOAT ueq;                        ///< equilibrium internal energy
  FLOAT dt_therm;                   ///< thermalization time scale
  FLOAT rad_pres[ndim];             ///< Acceleration from radiation pressure cmscott
  int ionstate;                     ///< States current ionisation state of the particle
                                    ///< (0 is neutral, 1 is smoothed and 2 is ionised)

  Particle() {
    active = false;
    potmin = false;
    iorig = -1;
    itype = gas;
    ptype = gas_type;
    level = 0;
    nstep = 0;
    nlast = 0;
    sinkid = -1;
    for (int k=0; k<ndim; k++) r[k] = (FLOAT) 0.0;
    for (int k=0; k<ndim; k++) v[k] = (FLOAT) 0.0;
    for (int k=0; k<ndim; k++) a[k] = (FLOAT) 0.0;
    for (int k=0; k<ndim; k++) r0[k] = (FLOAT) 0.0;
    for (int k=0; k<ndim; k++) v0[k] = (FLOAT) 0.0;
    for (int k=0; k<ndim; k++) a0[k] = (FLOAT) 0.0;
    for (int k=0; k<ndim; k++) agrav[k] = (FLOAT) 0.0;
    for (int k=0; k<ndim; k++) a_dust[k] = (FLOAT) 0.0;
    m         = (FLOAT) 0.0;
    h         = (FLOAT) 0.0;
    h_dust    = (FLOAT) 0.0;
    hrangesqd = (FLOAT) 0.0;
    invh      = (FLOAT) 0.0;
    hfactor   = (FLOAT) 0.0;
    invrho    = (FLOAT) 0.0;
    rho       = (FLOAT) 0.0;
    sound     = (FLOAT) 0.0;
    u         = (FLOAT) 0.0;
    u0        = (FLOAT) 0.0;
    dudt      = (FLOAT) 0.0;
    dudt0     = (FLOAT) 0.0;
    gpot      = (FLOAT) 0.0;
    dt        = (DOUBLE) 0.0;
    tlast     = (DOUBLE) 0.0;
    ionfrac   = (FLOAT) 0.999;
    Xion      = (FLOAT) 0.999;
    mu_bar    = (FLOAT) 1.0;
  }

};



//=================================================================================================
//  Structure SphParticle
/// \brief  Main base SPH particle data structure.
/// \author D. A. Hubber, G. Rosotti
/// \date   01/10/2013
//=================================================================================================
template <int ndim>
struct SphParticle : public Particle<ndim>
{
  FLOAT pfactor;                    ///< Pressure factor in SPH EOM
  FLOAT div_v;                      ///< Velocity divergence
  FLOAT alpha;                      ///< Artificial viscosity alpha value
  FLOAT dalphadt;                   ///< Rate of change of alpha


  SphParticle()
  {
    pfactor  = (FLOAT) 0.0;
    div_v    = (FLOAT) 0.0;
    alpha    = (FLOAT) 0.0;
    dalphadt = (FLOAT) 0.0;
  }

};



//=================================================================================================
//  Structure GradhSphParticle
/// \brief  `grad-h' SPH particle data structure.
/// \author D. A. Hubber, G. Rosotti
/// \date   01/10/2013
//=================================================================================================
template <int ndim>
struct GradhSphParticle : public SphParticle<ndim>
{
  FLOAT invomega;                   ///< grad-h omega/f correction term
  FLOAT zeta;                       ///< grad-h gravity correction term
  //FLOAT chi;                        ///< grad-h star-gravity correction term

  GradhSphParticle () {
    invomega = (FLOAT) 1.0;
    zeta = (FLOAT) 0.0;
    //chi = (FLOAT) 0.0;
  }

#ifdef MPI_PARALLEL
  static MPI_Datatype CreateMpiDataType() {
    MPI_Datatype particle_type;
    MPI_Datatype types[1] = {MPI_BYTE};
    MPI_Aint offsets[1] = {0};
    int blocklen[1] = {sizeof(GradhSphParticle<ndim>)};

    MPI_Type_create_struct(1,blocklen,offsets,types,&particle_type);

    return particle_type;
  }
#endif

};



//=================================================================================================
//  Structure SM2012SphParticle
/// \brief  Saitoh & Makino (2012) SPH particle data structure.
/// \author D. A. Hubber, G. Rosotti
/// \date   01/10/2013
//=================================================================================================
template <int ndim>
struct SM2012SphParticle : public SphParticle<ndim>
{
  FLOAT q;                          ///< Internal energy density
  FLOAT invq;                       ///< 1 / q

  SM2012SphParticle () {
    q = (FLOAT) 0.0;
    invq = (FLOAT) 0.0;
  }

#ifdef MPI_PARALLEL
  static MPI_Datatype CreateMpiDataType() {
    MPI_Datatype particle_type;
    MPI_Datatype types[1] = {MPI_BYTE};
    MPI_Aint offsets[1] = {0};
    int blocklen[1] = {sizeof(SM2012SphParticle<ndim>)};

    MPI_Type_create_struct(1,blocklen,offsets,types,&particle_type);

    return particle_type;
  }
#endif

};



//=================================================================================================
//  Structure MeshlessFVParticle
/// \brief  Main base Meshless Finite-Volume particle data structure.
/// \author D. A. Hubber, J. Ngoumou
/// \date   19/02/2015
//=================================================================================================
template <int ndim>
struct MeshlessFVParticle : public Particle<ndim>
{
  FLOAT invh;                          ///< 1 / h
  FLOAT hfactor;                       ///< invh^(ndim + 1)
  FLOAT invrho;                        ///< 1 / rho
  FLOAT invomega;                      ///< ..
  FLOAT press;                         ///< Thermal pressure
  FLOAT pfactor;                       ///< Pressure factor in SPH EOM
  FLOAT sound;                         ///< Sound speed
  FLOAT div_v;                         ///< Velocity divergence
  FLOAT vsig_max;                      ///< Maximum signal velocity to all neighbours
  FLOAT ndens;                         ///< Number density of neighbours
  FLOAT volume;                        ///< 'Volume' of particle
  FLOAT zeta;                          ///< ..
  FLOAT B[ndim][ndim];                 ///< Inverse matrix for gradient calculations
  FLOAT Wprim[ndim+2];                 ///< ..
  FLOAT Wmin[ndim+2];                  ///< ..
  FLOAT Wmax[ndim+2];                  ///< ..
  FLOAT Wmidmax[ndim+2];               ///< ..
  FLOAT Wmidmin[ndim+2];               ///< ..
  FLOAT Qcons[ndim+2];                 ///< ..
  FLOAT Qcons0[ndim+2];                ///< ..
  FLOAT grad[ndim+2][ndim];            ///< ..
  FLOAT dQ[ndim+2];                    ///< ..
  FLOAT dQdt[ndim+2];                  ///< Time derivative of conserved variables
  FLOAT alpha_slope[ndim+2];           ///< ..
  FLOAT Utot;                          ///< ..
  FLOAT rdmdt[ndim];                   ///< ..
  FLOAT rdmdt0[ndim];                  ///< ..
  FLOAT vreg[ndim];                    ///< ..


  // SPH particle constructor to initialise all values
  //-----------------------------------------------------------------------------------------------
  MeshlessFVParticle()
  {
    invh      = (FLOAT) 0.0;
    hfactor   = (FLOAT) 0.0;
    invrho    = (FLOAT) 0.0;
    invomega  = (FLOAT) 1.0;
    press     = (FLOAT) 0.0;
    pfactor   = (FLOAT) 0.0;
    sound     = (FLOAT) 0.0;
    div_v     = (FLOAT) 0.0;
    ndens     = (FLOAT) 0.0;
    volume    = (FLOAT) 0.0;
    vsig_max  = (FLOAT) 0.0;
    zeta      = (FLOAT) 0.0;
  }

};
#endif
