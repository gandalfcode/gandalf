//=================================================================================================
//  SinkParticle.h
//  Sink particle structure
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


#ifndef _SINKPARTICLE_H_
#define _SINKPARTICLE_H_


#include "Precision.h"
#include "StarParticle.h"

//=================================================================================================
//  Class SinkParticle
/// \brief   Individual sink particle data structure
/// \details Contains information relative to the sink and a pointer to the parent star
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
#endif
