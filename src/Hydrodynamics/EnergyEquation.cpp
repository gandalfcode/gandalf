//=================================================================================================
//  EnergyEquation.cpp
//  Contains basic (empty) functions for Energy Equation with templates.
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
#include "Sph.h"
#include "SmoothingKernel.h"
#include "Particle.h"
#include "EOS.h"
#include "EnergyEquation.h"

#include "../Headers/Integration.h"
#include "Debug.h"
using namespace std;



//=================================================================================================
//  EnergyEquation::EnergyEquation()
/// EnergyEquation constructor
//=================================================================================================
template <int ndim>
EnergyEquation<ndim>::EnergyEquation(DOUBLE energy_mult_aux) :
  energy_mult(energy_mult_aux)
{
}



//=================================================================================================
//  EnergyEquation::~EnergyEquation()
/// EnergyEquation destructor
//=================================================================================================
template <int ndim>
EnergyEquation<ndim>::~EnergyEquation()
{
}



// Class instances for each dimensionality (1, 2 and 3)
template class EnergyEquation<1>;
template class EnergyEquation<2>;
template class EnergyEquation<3>;
template class NullEnergy<1>;
template class NullEnergy<2>;
template class NullEnergy<3>;
