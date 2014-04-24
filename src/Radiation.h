//=============================================================================
//  Radiation.h
//  Contains definitions for all classes that control the transport of 
//  radiation through the computational domain.
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


#ifndef _RADIATION_H_
#define _RADIATION_H_


#include <map>
#include <string>
#include <list>
#include "CodeTiming.h"
#include "EnergyEquation.h"
#include "DomainBox.h"
#include "Nbody.h"
#include "Precision.h"
#include "Parameters.h"
#include "SimUnits.h"
#include "Sinks.h"
#include "SphKernel.h"
#include "SphNeighbourSearch.h"
using namespace std;



//=============================================================================
//  Class Radiation
/// \brief   Main base radiation class
/// \details Main base radiation class from which child classes containing 
///          implementations are inherited.
/// \author  D. A. Hubber
/// \date    21/04/2014
//=============================================================================
template <int ndim>
class Radiation
{
 public:

  Radiation();
  ~Radiation();

  virtual void UpdateRadiationField(void) = 0;

};



//=============================================================================
//  Class VoronoiMonteCarloRadiation
/// \brief   Monte-Carlo radiation transport using Voronoi tessellation
/// \details ..
/// \author  D. A. Hubber
/// \date    21/04/2014
//=============================================================================
template <int ndim>
class Radiation
{
 public:

  VoronoiMonteCarloRadiation();
  ~VoronoiMonteCarloRadiation();
  
  virtual void UpdateRadiationField(void);

};
#endif
