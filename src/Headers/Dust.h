//=================================================================================================
//  Dust.h
//  Contains main parent virtual class plus child classes for various dust
//  algorithms that are implemented.
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
// Guard against multiple inclusions
#ifndef __GANDALF_DUST_H__
#define __GANDALF_DUST_H__

#include "Parameters.h"
#include "Particle.h"
#include "Precision.h"
#include "Constants.h"
#include "CodeTiming.h"
#include "InlineFuncs.h"
#include "Tree.h"



//=================================================================================================
//  Class DustBase
/// \brief   DustBase class definition.
/// \details  Base class that defines the interface for computing drag forces.
/// \author  R. A. Booth
/// \date    17/10/2015
//=================================================================================================
template<int ndim>
class DustBase
{
protected:
	DustBase() { } ;
public:
	virtual ~DustBase() {} ;
	virtual void UpdateAllDragForces(int , int, Particle<ndim> *) = 0 ;
};

//=================================================================================================
//  Class DustFactory
/// \brief   DustFactory function definition
/// \details  Selects the appropriate dust implementation based on the parameters
/// \author  R. A. Booth
/// \date    23/10/2015
//=================================================================================================
template<int ndim, template<int> class ParticleType>
class DustFactory
{
public:
  static DustBase<ndim>* ProcessParameters(Parameters * params,
							               TreeBase<ndim>* t, TreeBase<ndim>* ghost,
							               TreeBase<ndim>* mpi_tree)  ;
} ;

#endif//__GANDALF_DUST_H__
