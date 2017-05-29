//=================================================================================================
//  LocalIsotherm.cpp
//  Contains functions for a locally isothermal EOS.
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
 
 
 
#include <math.h>

#include "../Headers/Integration.h"
#include "Debug.h"
#include "Constants.h"
#include "EOS.h"
#include "Hydrodynamics.h"
#include "Sinks.h"
#include "Parameters.h"
using namespace std;
 
 
 
//=================================================================================================
//  LocalIsotherm::LocalIsotherm()
/// LocalIsotherm class constructor
//=================================================================================================
template <int ndim>
LocallyIsothermal<ndim>::LocallyIsothermal(Parameters* simparams, SimUnits *units):
  Isothermal<ndim>(simparams, units),
  templaw(simparams->floatparams["templaw"]),
  tempmin(simparams->floatparams["tempmin"]/units->temp.outscale),
  nbody(0)
{
}
 

//=================================================================================================
//  LocalIsotherm::LocalIsotherm()
/// LocalIsotherm class desstructor
//=================================================================================================
template <int ndim>
LocallyIsothermal<ndim>::~LocallyIsothermal()
{

}

//=================================================================================================
//  LocalIsotherm::Temperature
/// Set the temperature given distance to nearest star
//=================================================================================================
template <int ndim>
FLOAT LocallyIsothermal<ndim>::Temperature(Particle<ndim> & part)
{
  FLOAT stardistmin=1e30;
  StarParticle<ndim>* star = nbody->stardata;

  // Compute distance to closest star
  for(int i=0; i<nbody->Nstar; i++){
    FLOAT dr[ndim] ;
    for (int j=0; j<ndim;j++) dr[j] = part.r[j]-star[i].r[j];
    FLOAT stardist = DotProduct(dr,dr,ndim);

    stardistmin = std::min(stardist, stardistmin) ;
  }
  stardistmin = sqrt(stardistmin) ;

  // Compute temperature
  FLOAT temp = max(temp0*pow(stardistmin,-templaw), tempmin);
  return temp;
}

//=================================================================================================
//  LocalIsotherm::SpecificInternalEnergy
/// Set the internal energt given distance to nearest star
//=================================================================================================
template <int ndim>
FLOAT LocallyIsothermal<ndim>::SpecificInternalEnergy(Particle<ndim> & part)
{
  FLOAT temp = LocallyIsothermal<ndim>::Temperature(part);
  return temp/gammam1/mu_bar;
}
 
template class LocallyIsothermal<1>;
template class LocallyIsothermal<2>;
template class LocallyIsothermal<3>;

