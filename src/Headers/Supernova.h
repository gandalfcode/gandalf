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
#include "DomainBox.h"
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
  void SupernovaInjection(const int, const int, const int, const int, const FLOAT,
                          FLOAT *, FLOAT, FLOAT, FLOAT, FLOAT,
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

  virtual void Update(const int, const int, const int, const FLOAT, Hydrodynamics<ndim> *, NeighbourSearch<ndim> *, RandomNumber *) = 0;

};



//=================================================================================================
//  Class NullSupernovaDriver
/// \brief   ...
/// \details ...
/// \author  D. A. Hubber
/// \date    18/02/2016
//=================================================================================================
template <int ndim>
class NullSupernovaDriver : public SupernovaDriver<ndim>
{
public:

  NullSupernovaDriver() {};
  ~NullSupernovaDriver() {};

  virtual void Update(const int, const int, const int, const FLOAT, Hydrodynamics<ndim> *,
                      NeighbourSearch<ndim> *, RandomNumber *) {};

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
  using SupernovaDriver<ndim>::Nsupernova;

  Supernova<ndim>  supernova;                 ///< ..
  FLOAT tsupernova;

  SedovTestDriver(Parameters *params, SimUnits &units);
  ~SedovTestDriver();

  virtual void Update(const int, const int, const int, const FLOAT,
                      Hydrodynamics<ndim> *, NeighbourSearch<ndim> *, RandomNumber *);

};



//=================================================================================================
//  Class RandomSedovTestDriver
/// \brief   Random supernova driver to test multiple explosion with block timesteps.
/// \details Random supernova driver to test multiple explosion with block timesteps.
/// \author  D. A. Hubber
/// \date    06/09/2016
//=================================================================================================
template <int ndim>
class RandomSedovTestDriver : public SupernovaDriver<ndim>
{
public:
  using SupernovaDriver<ndim>::Nsupernova;

  FLOAT tnext;                                   ///< Time for next supernova explosion
  FLOAT tsupernova;                              ///< Average time inbetween supernovae
  Supernova<ndim> supernova;                     ///< Instance of supernova class
  DomainBox<ndim> *simbox;                       ///< Pointer to simulation domain box

  RandomSedovTestDriver(Parameters *params, SimUnits &units, DomainBox<ndim> &);
  ~RandomSedovTestDriver();

  virtual void Update(const int, const int, const int, const FLOAT,
                      Hydrodynamics<ndim> *, NeighbourSearch<ndim> *, RandomNumber *);

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
