//=============================================================================
//  BinaryOrbit.h
//  Contains definitions for binary orbit data structure.
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


#ifndef _BINARY_ORBIT_H_
#define _BINARY_ORBIT_H_


#include <iostream>
#include <string>
#include "Precision.h"
#include "Constants.h"
using namespace std;



//=============================================================================
//  Structure BinaryOrbit
/// \brief   Class definition for individual binary orbit.
/// \details Class definition for individual binary orbit.
/// \author  D. A. Hubber, G. Rosotti
/// \date    01/08/2013
//=============================================================================
struct BinaryOrbit {
  int ichild1;                     ///< id of first component (system or star)
  int ichild2;                     ///< id of second component (system or star)
  int Nstar;                       ///< Total no. of stars in system
  FLOAT r[3];                      ///< Position of COM of binary
  FLOAT v[3];                      ///< Velocity of COM of binary
  FLOAT m;                         ///< Total mass of binary
  FLOAT angmom[3];                 ///< Angular momentum of binary
  FLOAT binen;                     ///< Specific binding energy
  FLOAT sma;                       ///< Semi-major axis
  FLOAT ecc;                       ///< Orbital eccentricity
  FLOAT period;                    ///< Orbital period
  FLOAT q;                         ///< Mass ratio
  string systemtype;               ///< Type of system (binary, triple, etc..)
};
#endif
