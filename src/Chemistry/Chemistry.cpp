//=================================================================================================
//  Chemistry.cpp
//  All routines for creating new sinks and accreting gas and updating all
//  sink particle propterties.
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
#include <math.h>
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
/// Chemistry class constructor
//=================================================================================================
template <int ndim>
Chemistry<ndim>::Chemistry(Parameters *params)
{
  // Local references to parameter variables for brevity
  map<string, int> &intparams = params->intparams;
  map<string, double> &floatparams = params->floatparams;
  map<string, string> &stringparams = params->stringparams;


  initialAbundances[AB_H2] = floatparams["ab_h2"];
  initialAbundances[AB_HP] = floatparams["ab_hp"];

}



//=================================================================================================
//  ...
/// Chemistry class destructor
//=================================================================================================
template <int ndim>
Chemistry<ndim>::~Chemistry()
{
}



//=================================================================================================
//  Chemistry::EvolveAbundances
/// ...
//=================================================================================================
template <int ndim>
void SCOChemistry<ndim>::EvolveAbundances(ChemistryParticle &part, FLOAT dt, SimUnits &units)
{
  FLOAT dt_cgs = dt*units.t.outscale*units.t.outcgs;
  FLOAT rho_cgs = part.hydropart->rho*units.rho.outscale*units.rho.outcgs;
  FLOAT u_cgs = part.hydropart->u*units.u.outscale*units.u.outcgs;

  FLOAT u_new = ExternalEvolveAbundances(...);

  u_new = u_new/units.u.outscale/units.u.outcgs;
  part.hydropart->u = u_new;

  return;
}
