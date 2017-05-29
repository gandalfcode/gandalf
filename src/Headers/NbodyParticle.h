//=================================================================================================
//  NbodyParticle.h
//  Main N-body particle data structure.
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


#ifndef _NBODY_PARTICLE_H_
#define _NBODY_PARTICLE_H_


#include "Precision.h"
#include "Constants.h"
#include "Flags.h"


//=================================================================================================
//  Class NbodyParticle
/// \brief   N-body particle data structure
/// \details Main parent N-body particle data structure.  All main other N-body particle types
///          (e.g. stars, systems, sinks) are derived from this class.
/// \author  D. A. Hubber
/// \date    15/04/2013
//=================================================================================================
template <int ndim>
class NbodyParticle
{
 public:

  type_flag flags;                     ///< Flags for active / end of time step
  int istar;                           ///< Internal i.d.
  int Ncomp;                           ///< No. of internal components
  int level;                           ///< Current timestep level
  int nstep;                           ///< Integer step-size of particle
  int nlast;                           ///< Integer time at beginning of step
  FLOAT r[ndim];                       ///< Position
  FLOAT v[ndim];                       ///< Velocity
  FLOAT a[ndim];                       ///< Acceleration
  FLOAT adot[ndim];                    ///< Time derivative of acceleration (jerk)
  FLOAT a2dot[ndim];                   ///< 2nd time derivative of acceleration
  FLOAT a3dot[ndim];                   ///< 3rd time derivative of acceleration
  FLOAT r0[ndim];                      ///< Position at beginning of step
  FLOAT v0[ndim];                      ///< Velocity at beginning of step
  FLOAT a0[ndim];                      ///< Acceleration at beginning of step
  FLOAT adot0[ndim];                   ///< Jerk at beginning of step
  FLOAT a2dot0[ndim];                  ///< 2nd time derivative at beginning of step
  FLOAT apert[ndim];                   ///< Acceleration due to perturbers
  FLOAT adotpert[ndim];                ///< Jerk due to perturbers
  FLOAT m;                             ///< Star mass
  FLOAT h;                             ///< Smoothing length
  FLOAT invh;                          ///< 1 / h
  FLOAT radius;                        ///< Softening/sink radius of particle
  //FLOAT hfactor;                       ///< invh^(ndim + 1)
  FLOAT gpot;                          ///< Gravitational potential
  FLOAT gpe;                           ///< Gravitational potential energy
  FLOAT gpe_internal;                  ///< Internal grav. potential energy
  FLOAT gpe_pert;                      ///< Perturber grav. potential energy
  DOUBLE dt;                           ///< Particle timestep
  DOUBLE dt_next;                      ///< Particle timestep for next step
  DOUBLE dt_internal;                  ///< Internal timestep (e.g. due to sub-systems)
  DOUBLE tlast;                        ///< Time at beginning of last step
  DOUBLE NLyC;                         ///< No. of ionising photons per second


  // Star particle constructor to initialise all values
  //---------------------------------------------------------------------------
  NbodyParticle()
  {
    level  = 0;
    nstep  = 0;
    nlast  = 0;
    Ncomp  = 1;
    for (int k=0; k<ndim; k++) r[k]        = (FLOAT) 0.0;
    for (int k=0; k<ndim; k++) v[k]        = (FLOAT) 0.0;
    for (int k=0; k<ndim; k++) a[k]        = (FLOAT) 0.0;
    for (int k=0; k<ndim; k++) adot[k]     = (FLOAT) 0.0;
    for (int k=0; k<ndim; k++) a2dot[k]    = (FLOAT) 0.0;
    for (int k=0; k<ndim; k++) a3dot[k]    = (FLOAT) 0.0;
    for (int k=0; k<ndim; k++) r0[k]       = (FLOAT) 0.0;
    for (int k=0; k<ndim; k++) v0[k]       = (FLOAT) 0.0;
    for (int k=0; k<ndim; k++) a0[k]       = (FLOAT) 0.0;
    for (int k=0; k<ndim; k++) adot0[k]    = (FLOAT) 0.0;
    for (int k=0; k<ndim; k++) apert[k]    = (FLOAT) 0.0;
    for (int k=0; k<ndim; k++) adotpert[k] = (FLOAT) 0.0;
    m            = 0;
    h            = 0;
    invh         = (FLOAT) 0.0;
    //hfactor      = (FLOAT) 0.0;
    gpot         = (FLOAT) 0.0;
    gpe          = (FLOAT) 0.0;
    gpe_internal = (FLOAT) 0.0;
    gpe_pert     = (FLOAT) 0.0;
    dt           = 0.0;
    dt_internal  = big_number;
    tlast        = 0.0;
    NLyC         = 0.0;
  }

};
#endif
