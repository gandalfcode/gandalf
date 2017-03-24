//=================================================================================================
//  SupernovaDriver.cpp
//  Routines for creating supernova explosions in GANDALF simulations.
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


#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include "Constants.h"
#include "Debug.h"
#include "Exception.h"
#include "InlineFuncs.h"
#include "Parameters.h"
#include "Particle.h"
#include "Precision.h"
#include "RandomNumber.h"
#include "Supernova.h"
using namespace std;



//=================================================================================================
//  SedovTestDriver::SedovTestDriver
/// ...
//=================================================================================================
template <int ndim>
SedovTestDriver<ndim>::SedovTestDriver
 (SimulationBase* sim, Parameters *params, SimUnits &units) : SupernovaDriver<ndim>(sim)
{
  tsupernova = 1.0;
}



//=================================================================================================
//  SedovTestDriver::Update
/// ...
//=================================================================================================
template <int ndim>
void SedovTestDriver<ndim>::Update
 (const int n,                               ///< Current integer time
  const int level_step,                      ///< ..
  const int level_max,                       ///< Max. block timestep level
  const FLOAT t,                             ///< Physical simulation time
  Hydrodynamics<ndim> *hydro,                ///< Pointer to hydrodynamics object
  NeighbourSearch<ndim> *neibsearch,         ///< Pointer to neighbour search object
  RandomNumber *randnumb)                    ///< Pointer to random number generator
{

  // Create supernova at requested time (and ensure that only one is made)
  //-----------------------------------------------------------------------------------------------
  if (Nsupernova == 0 && t >= tsupernova) {
#ifdef OUTPUT_ALL
    cout << "ADDING SUPERNOVA!!" << endl;
#endif

    FLOAT SNpos[ndim];
    FLOAT Einj        = (FLOAT) 0.01;
    FLOAT R_therm_kin = (FLOAT) 100000.0;
    FLOAT Minj        = (FLOAT) 0.005;
    FLOAT Rinj        = hydro->GetParticlePointer(0).h; //(FLOAT) 0.0;
    for (int k=0; k<ndim; k++) SNpos[k] = (FLOAT) 0.0;

    supernova.SupernovaInjection(n, level_step, level_max, Nsupernova, t, SNpos, Einj,
                                 R_therm_kin, Minj, Rinj, hydro, neibsearch, randnumb);
    Nsupernova++;
  }
  //-----------------------------------------------------------------------------------------------

  return;
}



//=================================================================================================
//  RandomSedovTestDriver::SedovTestDriver
/// ...
//=================================================================================================
template <int ndim>
RandomSedovTestDriver<ndim>::RandomSedovTestDriver
 (SimulationBase* sim, Parameters *params, SimUnits &units, DomainBox<ndim> &_simbox) :
  SupernovaDriver<ndim>(sim)
{
  tsupernova = 0.5;
  tnext      = 0.5*tsupernova;
  simbox     = &(_simbox);
}



//=================================================================================================
//  RandomTestDriver::Update
/// ...
//=================================================================================================
template <int ndim>
void RandomSedovTestDriver<ndim>::Update
 (const int n,                               ///< Current integer time
  const int level_step,                      ///< ..
  const int level_max,                       ///< Max. block timestep level
  const FLOAT t,                             ///< Physical simulation time
  Hydrodynamics<ndim> *hydro,                ///< Pointer to hydrodynamics object
  NeighbourSearch<ndim> *neibsearch,         ///< Pointer to neighbour search object
  RandomNumber *randnumb)                    ///< Pointer to random number generator
{
  // Use random number generator to decide the position of the supernova
  //-----------------------------------------------------------------------------------------------
  if (t >= tnext) {

    FLOAT SNpos[ndim];
    FLOAT Einj        = (FLOAT) 0.01;
    FLOAT R_therm_kin = (FLOAT) 100000.0;
    FLOAT Minj        = (FLOAT) 0.005;
    FLOAT Rinj        = hydro->GetParticlePointer(0).h; //(FLOAT) 0.0;
    for (int k=0; k<ndim; k++) SNpos[k] = simbox->min[k] + randnumb->floatrand()*simbox->size[k];

    supernova.SupernovaInjection(n, level_step, level_max, Nsupernova, t, SNpos, Einj,
                                 R_therm_kin, Minj, Rinj, hydro, neibsearch, randnumb);
    Nsupernova++;

#ifdef OUTPUT_ALL
    cout << "ADDING SUPERNOVA!   Nsupernova : " << Nsupernova << endl;
#endif
    tnext = ((FLOAT) Nsupernova + 0.5)*tsupernova;
  }
  //-----------------------------------------------------------------------------------------------

  return;
}



//=================================================================================================
//  SilccSupernovaDriver::SilccSupernovaDriver
/// ...
//=================================================================================================
template <int ndim>
SilccSupernovaDriver<ndim>::SilccSupernovaDriver
 (SimulationBase* sim, Parameters *params, SimUnits &units) : SupernovaDriver<ndim>(sim)
{
  // Local references to parameter variables for brevity
  map<string, double> &floatparams = params->floatparams;
  map<string, string> &stringparams = params->stringparams;


  // Parameters to be read from the file
  ifstream file;
  ifstream SNinfile;
  string SNfile;
  string line;
  int i;
  int nSN;

  Minj = floatparams["Minj"]/(units.m.outscale);
  Rinj = floatparams["Rinj"]/(units.r.outscale);
  SNfile = stringparams["SNfile_name"];
  R_therm_kin = floatparams["R_therm_kin"];
  SNid = 0;

  cout << "Minj: " << Minj << endl;
  cout << "Rinj: " << Rinj << endl;

  file.open(SNfile.c_str(), ios::in);


  //-----------------------------------------------------------------------------------------------
  if (file.good()) {
    i = 0;
    getline(file, line);
    getline(file, line);
    istringstream istr(line);
    istr >> nSN;

    SNtime = new FLOAT[nSN];
    SNposx = new FLOAT[nSN];
    SNposy = new FLOAT[nSN];
    SNposz = new FLOAT[nSN];
    SNEinj = new FLOAT[nSN];

    // read table
    //---------------------------------------------------------------------------------------------
    while (getline(file, line)) {
      istringstream istr(line);

      if (istr >> SNtime[i] >> SNposx[i] >> SNposy[i] >> SNposz[i] >> SNEinj[i]) {

        SNtime[i] = SNtime[i]/(units.t.outscale*units.t.outcgs);
        SNposx[i] = SNposx[i]/(units.r.outscale*units.r.outcgs);
        SNposy[i] = SNposy[i]/(units.r.outscale*units.r.outcgs);
        SNposz[i] = SNposz[i]/(units.r.outscale*units.r.outcgs);
        SNEinj[i] = SNEinj[i]/(units.E.outscale*units.E.outcgs);
        i++;
      }
    }
    //---------------------------------------------------------------------------------------------

  }
  //-----------------------------------------------------------------------------------------------
  else {
    cout << "Dateifehler" << endl;
  }
  //-----------------------------------------------------------------------------------------------

  SNinfile.close();
}



//=================================================================================================
//  SilccSupernovaDriver::Update
/// ...
//=================================================================================================
template <int ndim>
void SilccSupernovaDriver<ndim>::Update
 (const int n,                               ///< Current integer time
  const int level_step,                      ///< ..
  const int level_max,                       ///< Max. block timestep level
  const FLOAT t,                             ///< Physical simulation time
  Hydrodynamics<ndim> *hydro,                ///< Pointer to hydrodynamics object
  NeighbourSearch<ndim> *neibsearch,         ///< Pointer to neighbour search object
  RandomNumber *randnumb)                    ///< Pointer to random number generator
{

  // Create supernova at requested time (and ensure that only one is made)
  //-----------------------------------------------------------------------------------------------
  if (t >= SNtime[SNid]) {

    cout << "ADDING SUPERNOVA " << SNid <<  " !!" << endl;

    FLOAT SNpos[ndim];
    //FLOAT Rinj        = hydro->GetParticlePointer(0).h; //(FLOAT) 0.0;
    SNpos[0] = SNposx[SNid];
    SNpos[1] = SNposy[SNid];
    SNpos[2] = SNposz[SNid];

    supernova.SupernovaInjection(n, level_step, level_max, SNid, t, SNpos, SNEinj[SNid],
                                 R_therm_kin, Minj, Rinj, hydro, neibsearch, randnumb);
    SNid++;
  }
  //-----------------------------------------------------------------------------------------------

  return;
}



template class SedovTestDriver<1>;
template class RandomSedovTestDriver<1>;
template class SilccSupernovaDriver<1>;

template class SedovTestDriver<2>;
template class RandomSedovTestDriver<2>;
template class SilccSupernovaDriver<2>;

template class SedovTestDriver<3>;
template class RandomSedovTestDriver<3>;
template class SilccSupernovaDriver<3>;
