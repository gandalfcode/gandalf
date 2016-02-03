//=================================================================================================
//  Supernova.h
//  Contains all definitions for classes related to creating Supernova via adding new particles
//  with a given thermal or kinetic energy input.
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


#ifndef _SUPERNOVA_H_
#define _SUPERNOVA_H_


#include <assert.h>
#include <string>
#include "Precision.h"
#include "Constants.h"
#include "Hydrodynamics.h"
#include "NeighbourSearch.h"
#include "Particle.h"


//=================================================================================================
//  Class Supernova
/// \brief   Main parent Supernova class
/// \details Injects a single Supernova; adds accelerated particles within injection radius
/// \author  S. Walch, D. A. Hubber
/// \date    19/01/2016
//=================================================================================================
template <int ndim>
class Supernova
{
public:

  Supernova();
  ~Supernova();

  // pass SN position, Einj, Etherm/Ekin, Minj, Rinj, SNid, hydro
  void SupernovaInjection(int, int, int, FLOAT, FLOAT *, FLOAT, FLOAT, FLOAT, FLOAT,
                          Hydrodynamics<ndim> *, NeighbourSearch<ndim> *, RandomNumber *);

};
#endif
