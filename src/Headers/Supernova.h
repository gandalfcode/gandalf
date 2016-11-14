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
	SimulationBase* sim; ///< Pointer to simulation object
public:

  Supernova(SimulationBase* _sim): sim(_sim) {};
  ~Supernova() {};

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
protected:
  int Nsupernova;                      ///< Number of supernovae
  Supernova<ndim>  supernova;                 ///< Supernova object


public:

  SupernovaDriver(SimulationBase* sim): supernova(sim) {Nsupernova = 0;}
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

  NullSupernovaDriver(SimulationBase* sim): SupernovaDriver<ndim>(sim) {};
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
  using SupernovaDriver<ndim>::supernova;

  FLOAT tsupernova;

  SedovTestDriver(SimulationBase* sim, Parameters *params, SimUnits &units);
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
  using SupernovaDriver<ndim>::supernova;


  FLOAT tnext;                                   ///< Time for next supernova explosion
  FLOAT tsupernova;                              ///< Time between supernovae
  DomainBox<ndim> *simbox;                       ///< Pointer to simulation domain box

  RandomSedovTestDriver(SimulationBase* sim, Parameters *params, SimUnits &units, DomainBox<ndim> &);
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
template <int ndim>
class SilccSupernovaDriver : public SupernovaDriver<ndim>
{
public:
  using SupernovaDriver<ndim>::Nsupernova;
  using SupernovaDriver<ndim>::supernova;

  int SNid;
  FLOAT tsupernova;
  FLOAT Minj;
  FLOAT Rinj;
  FLOAT R_therm_kin;
  FLOAT *SNtime;
  FLOAT *SNposx;
  FLOAT *SNposy;
  FLOAT *SNposz;
  FLOAT *SNEinj;
  //Supernova<ndim>  supernova;                 ///< ..

  SilccSupernovaDriver(SimulationBase* sim, Parameters *params, SimUnits &units);
  ~SilccSupernovaDriver();

  virtual void Update(const int, const int, const int, const FLOAT,
                      Hydrodynamics<ndim> *, NeighbourSearch<ndim> *, RandomNumber *);

};
#endif
