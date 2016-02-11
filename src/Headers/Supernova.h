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



//=================================================================================================
//  Class SupernovaDriver
/// \brief   Main parent Supernova class
/// \details Injects a single Supernova; adds accelerated particles within injection radius
/// \author  S. Walch, D. A. Hubber
/// \date    19/01/2016
//=================================================================================================
template <int ndim>
class SupernovaDriver
{
public:

  int Nsupernova;                      ///< ..

  //SupernovaDriver(Parameters *params, SimUnits &units, bool restart, FLOAT time);
  SupernovaDriver() {Nsupernova = 0;}
  ~SupernovaDriver() {};

  void Update(FLOAT, Hydrodynamics *) = 0;

};



//=================================================================================================
//  Class SedovTestDriver
/// \brief   Simple Sedov-explosion driver for testing fidelity of explosions
/// \details Simple Sedov-explosion driver for testing fidelity of explosions
/// \author  D. A. Hubber
/// \date    10/02/2016
//=================================================================================================
template <int ndim>
class SedovTestDriver : public SupernovaDriver<ndim>
{
public:

  Supernova supernova;                 ///< ..
  FLOAT tsupernova;

  SedovTestDriver(Parameters *params, SimUnits &units, bool restart, FLOAT time);
  ~SedovTestDriver();

  void Update(FLOAT, Hydrodynamics<ndim> *);

};



//=================================================================================================
//  Class SilccSupernovaDriver
/// \brief   Main parent SupernovaDriving class
/// \details ..
/// \author  D. A. Hubber, S. Walch, T. Balduin, P. Rohde
/// \date    19/01/2016
//=================================================================================================
/*template <int ndim>
class SilccSupernovaDriving
{
public:
  // SN ID (defined in constructor)
  int SNid;

  // Global SN parameters that could go into the parameters file
  FLOAT R_therm_kin;  //ratio of thermal and kinetic E to be injected
  FLOAT Minj;         //ejecta mass (default 7 Msun)
  FLOAT Rinj;         //injection radius (default 1pc)

  int nSN;
  FLOAT *SNtime;
  FLOAT *SNposx;
  FLOAT *SNposy;
  FLOAT *SNposz;
  FLOAT *SNEinj;


  // Constructor
  SupernovaDriving(Parameters *params, SimUnits &units, bool restart, FLOAT time);
  // update pass current simulation time, time step, hydro
  void Update(FLOAT, Hydrodynamics*);

};*/
#endif
