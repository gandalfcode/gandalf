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
#include <math.h>
#include "Constants.h"
#include "Precision.h"
#include "NbodyParticle.h"
#include "StarParticle.h"
#include "Parameters.h"
#include "Sph.h"
#include "Nbody.h"
#include "Sinks.h"
#include "Debug.h"
#include "Exception.h"
#include "InlineFuncs.h"
using namespace std;



//=================================================================================================
//  ...
/// SupernovaDriving class constructor
// For now this class treats the Supernova driving from an external file
//=================================================================================================
template <int ndim>
SedovTestDriver<ndim>::SedovTestDriver
 (Parameters *params, SimUnits &units, bool restart, FLOAT time) : SupernovaDriver()
{
  // Local references to parameter variables for brevity
  map<string, int> &intparams = params->intparams;
  map<string, double> &floatparams = params->floatparams;
  map<string, string> &stringparams = params->stringparams;

  tsupernova = 1.0;
}



//=================================================================================================
//  SedovTestDriver::Update
/// ...
//=================================================================================================
template <int ndim>
void SedovTestDriver<ndim>::Update
 (int n,                                     ///< Current integer time
  int level_max,                             ///< Max. block timestep level
  FLOAT t,                                   ///< Physical simulation time
  Hydrodynamics<ndim> *hydro,                ///< Pointer to hydrodynamics object
  NeighbourSearch<ndim> *neibsearch,         ///< Pointer to neighbour search object
  RandomNumber *randnumb)                    ///< Pointer to random number generator
{

  FLOAT SNpos[ndim];  //current Supernova position
  FLOAT Einj;         //total energy of SN (default 1e51 erg)

    // Helper variables (local)
  FLOAT t_ratio1;


  if (SNid == -1) {
    cout << "Welcome, this is your first Supernova!";
    SNid =0;
  }

  // get last SN from table
  t_ratio1 = SNtime[SNid]/time;


  // check if a Supernova should be injected:
  while (t_ratio1<= 1.){
    // inject SN now!...
    cout << "Injecting SN number "<< SNid;

    // Fill entries to pass to SupernovaInjection
    SNpos[0]=SNposx[SNid];
    SNpos[1]=SNposy[SNid];
    SNpos[2]=SNposz[SNid];
    Einj = SNEinj[SNid];

  supernova.SupernovaInjection(SNpos, Einj, R_therm_kin, Minj, Rinj, SNid, hydro);


  //return
}
