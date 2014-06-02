//=============================================================================
//  SphParticle.h
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
//=============================================================================


#ifndef _SPH_PARTICLE_H_
#define _SPH_PARTICLE_H_


#include "Precision.h"
#include "Constants.h"
#ifdef MPI_PARALLEL
#include <stddef.h>
#include "mpi.h"
#include "Exception.h"
#endif


enum ptype{gas, icm, boundary, cdm, dead, 
           x_lhs_periodic, x_lhs_mirror, x_rhs_periodic, x_rhs_mirror,
           y_lhs_periodic, y_lhs_mirror, y_rhs_periodic, y_rhs_mirror,
           z_lhs_periodic, z_lhs_mirror, z_rhs_periodic, z_rhs_mirror,
	   Nsphtypes};
enum eosenum{isothermal, barotropic, barotropic2,
	     energy_eqn, constant_temp, nspheos};



//=============================================================================
//  Structure Particle
/// \brief  Main base particle data structure.
/// \author D. A. Hubber, G. Rosotti
/// \date   01/10/2013
//=============================================================================
template <int ndim>
struct Particle
{
  bool active;                      ///< Flag if active (i.e. recompute step)
  int iorig;                        ///< Original particle i.d.
  int itype;                        ///< SPH particle type
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
  FLOAT m;                          ///< Particle mass
  FLOAT h;                          ///< SPH smoothing length
  FLOAT hrangesqd;                  ///< Kernel extent (squared)
  FLOAT rho;                        ///< Density
  FLOAT u;                          ///< Specific internal energy
  FLOAT u0;                         ///< u at beginning of step
  FLOAT dudt0;                      ///< dudt at beginning of step
  FLOAT dudt;                       ///< Compressional heating rate
  FLOAT gpot;                       ///< Gravitational potential
  FLOAT gpe;                        ///< Gravitational potential energy
  DOUBLE dt;                        ///< Particle timestep

  Particle() {
    active = false;
    iorig = -1;
    itype = gas;
    level = 0;
    nstep = 0;
    nlast = 0;
    for (int k=0; k<ndim; k++) r[k] = (FLOAT) 0.0;
    for (int k=0; k<ndim; k++) v[k] = (FLOAT) 0.0;
    for (int k=0; k<ndim; k++) a[k] = (FLOAT) 0.0;
    for (int k=0; k<ndim; k++) r0[k] = (FLOAT) 0.0;
    for (int k=0; k<ndim; k++) v0[k] = (FLOAT) 0.0;
    for (int k=0; k<ndim; k++) a0[k] = (FLOAT) 0.0;
    for (int k=0; k<ndim; k++) agrav[k] = (FLOAT) 0.0;
    m = (FLOAT) 0.0;
    h = (FLOAT) 0.0;
    hrangesqd = (FLOAT) 0.0;
    rho = (FLOAT) 0.0;
    u = (FLOAT) 0.0;
    u0 = (FLOAT) 0.0;
    dudt = (FLOAT) 0.0;
    dudt0 = (FLOAT) 0.0;
    gpot = (FLOAT) 0.0;
    gpe = (FLOAT) 0.0;
    dt = (DOUBLE) 0.0;
  }

};



//=============================================================================
//  Structure SphParticle
/// \brief  Main base SPH particle data structure.
/// \author D. A. Hubber, G. Rosotti
/// \date   01/10/2013
//=============================================================================
template <int ndim>
struct SphParticle : public Particle<ndim>
{
  bool potmin;                      ///< Is particle at a potential minima?
  int levelneib;                    ///< Min. timestep level of neighbours
  int sinkid;                       ///< i.d. of sink particle
  FLOAT invh;                       ///< 1 / h
  FLOAT hfactor;                    ///< invh^(ndim + 1)
  FLOAT invrho;                     ///< 1 / rho
  FLOAT press;                      ///< Thermal pressure
  FLOAT pfactor;                    ///< Pressure factor in SPH EOM
  FLOAT sound;                      ///< Sound speed
  FLOAT div_v;                      ///< Velocity divergence
  FLOAT alpha;                      ///< Artificial viscosity alpha value
  FLOAT dalphadt;                   ///< Rate of change of alpha
  FLOAT adot[ndim];                 ///< 1st time derivative of acceleration


  // SPH particle constructor to initialise all values
  //---------------------------------------------------------------------------
  SphParticle()
  {
    potmin = false;
    levelneib = 0;
    sinkid = -1;
    invh = (FLOAT) 0.0;
    hfactor = (FLOAT) 0.0;
    invrho = (FLOAT) 0.0;
    press = (FLOAT) 0.0;
    pfactor = (FLOAT) 0.0;
    sound = (FLOAT) 0.0;
    div_v = (FLOAT) 0.0;
    alpha = (FLOAT) 0.0;
    dalphadt = (FLOAT) 0.0;
    for (int k=0; k<ndim; k++) adot[k] = (FLOAT) 0.0;
  }

#ifdef MPI_PARALLEL
  static MPI_Datatype CreateMpiDataType() {
    MPI_Datatype particle_type;
    MPI_Datatype types[1] = {MPI_BYTE};
    MPI_Aint offsets[1] = {0};
    int blocklen[1] = {sizeof(SphParticle<ndim>)};
    
    MPI_Type_create_struct(1,blocklen,offsets,types,&particle_type);
    
    return particle_type;
  }
#endif

};



//=============================================================================
//  Structure GradhSphParticle
/// \brief  `grad-h' SPH particle data structure.
/// \author D. A. Hubber, G. Rosotti
/// \date   01/10/2013
//=============================================================================
template <int ndim>
struct GradhSphParticle : public SphParticle<ndim>
{
  FLOAT invomega;                   ///< grad-h omega/f correction term
  FLOAT zeta;                       ///< grad-h gravity correction term
  FLOAT chi;                        ///< grad-h star-gravity correction term

  GradhSphParticle () {
    invomega = (FLOAT) 1.0;
    zeta = (FLOAT) 0.0;
    chi = (FLOAT) 0.0;
  }

};



//=============================================================================
//  Structure SM2012SphParticle
/// \brief  Saitoh & Makino (2012) SPH particle data structure.
/// \author D. A. Hubber, G. Rosotti
/// \date   01/10/2013
//=============================================================================
template <int ndim>
struct SM2012SphParticle : public SphParticle<ndim>
{
  FLOAT q;                          ///< Internal energy density
  FLOAT invq;                       ///< 1 / q

  SM2012SphParticle () {
    q = (FLOAT) 0.0;
    invq = (FLOAT) 0.0;
  }

};



//=============================================================================
//  Structure GodunovSphParticle
/// \brief  Godunov SPH particle data structure.
/// \author D. A. Hubber, G. Rosotti
/// \date   01/10/2013
//=============================================================================
template <int ndim>
struct GodunovSphParticle : public SphParticle<ndim>
{
  FLOAT gradrho[ndim];              ///< Density gradient
  FLOAT gradP[ndim];                ///< Pressure gradient
  FLOAT gradv[ndim][ndim];          ///< Velocity gradient matrix

  GodunovSphParticle() {
    for (int k=0; k<ndim; k++) gradrho[k] = (FLOAT) 0.0;
    for (int k=0; k<ndim; k++) gradP[k] = (FLOAT) 0.0;
    for (int k=0; k<ndim; k++)
      for (int kk=0; kk<ndim; kk++)
        gradv[k][kk] = (FLOAT) 0.0;
  }

};



//=============================================================================
//  Structure SphType
/// \brief  ..
/// \author D. A. Hubber, G. Rosotti
/// \date   10/02/2014
//=============================================================================
struct SphType
{
  eosenum eos;
};
#endif
