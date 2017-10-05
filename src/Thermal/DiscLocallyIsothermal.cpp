//=================================================================================================
//  DiscLocallyIsotherm.cpp
//  Contains functions for a locally isothermal EOS used in discs.
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
#include "Debug.h"
#include "Constants.h"
#include "EOS.h"
#include "Hydrodynamics.h"
#include "Sinks.h"
#include "Integration.h"
#include "Parameters.h"
using namespace std;



//=================================================================================================
//  DiscLocallyIsothermal::DiscLocallyIsothermal()
/// DiscLocallyIsothermal class constructor
//=================================================================================================
template <int ndim>
DiscLocallyIsothermal<ndim>::DiscLocallyIsothermal(Parameters* simparams, SimUnits *units):
  Isothermal<ndim>(simparams, units),
  slope(simparams->floatparams["DiscIcQ"]),
  norm(simparams->floatparams["DiscIcHr"]*sqrt(1./simparams->floatparams["DiscIcRin"])),
  rin(simparams->floatparams["DiscIcRin"]),
  nbody(0)
{
}


//=================================================================================================
//  DiscLocallyIsothermal::DiscLocallyIsothermal()
/// DiscLocallyIsothermal class desstructor
//=================================================================================================
template <int ndim>
DiscLocallyIsothermal<ndim>::~DiscLocallyIsothermal()
{

}

//=================================================================================================
//  DiscLocallyIsothermal::SoundSpeed
/// Set the sound speed given distance to the star
//=================================================================================================
template <int ndim>
FLOAT DiscLocallyIsothermal<ndim>::SoundSpeed(Particle<ndim> & part)
{
  StarParticle<ndim>* star = nbody->stardata;

  // Compute distance to star
  FLOAT dr[ndim] ;
  for (int j=0; j<ndim;j++) dr[j] = part.r[j]-star[0].r[j];
  FLOAT r2 = DotProduct(dr,dr,ndim);
  FLOAT r = sqrt(r2);

  // Compute sound speed
  FLOAT cs = norm*pow(r/rin,-slope);;
  return cs;
}

//=================================================================================================
//  DiscLocallyIsothermal::SpecificInternalEnergy
/// Set the internal energy given the position of the particle
//=================================================================================================
template <int ndim>
FLOAT DiscLocallyIsothermal<ndim>::SpecificInternalEnergy(Particle<ndim> & part)
{
  FLOAT cs = SoundSpeed(part);
  return cs*cs/gammam1;
}


//=================================================================================================
//  DiscLocallyIsothermal::Temperature
/// Set the temperature given the position of the particle
//=================================================================================================
template <int ndim>
FLOAT DiscLocallyIsothermal<ndim>::Temperature(Particle<ndim> & part)
{
  FLOAT cs = SoundSpeed(part);
  return mu_bar*cs*cs;
}



template class DiscLocallyIsothermal<1>;
template class DiscLocallyIsothermal<2>;
template class DiscLocallyIsothermal<3>;

