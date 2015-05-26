//=================================================================================================
//  SystemParticle.h
//  Main system particle data structure
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


#ifndef _SYSTEM_PARTICLE_H_
#define _SYSTEM_PARTICLE_H_


#include "Precision.h"
#include "Constants.h"
#include "NbodyParticle.h"


static const int Ncompmax = 4;
static const int Npertmax = 4;



//=================================================================================================
//  Structure SystemParticle
/// \brief  System particle data structure
/// \author D. A. Hubber, G. Rosotti
/// \date   10/05/2013
//=================================================================================================
template <int ndim>
class SystemParticle: public NbodyParticle<ndim>
{
public:
  int inode;                                ///< NN-tree node id
  int Nchildren;                            ///< Number of nbody children
  int Npert;                                ///< Number of perturbers
  NbodyParticle<ndim>* children[Ncompmax];  ///< Array of ptrs to children
  NbodyParticle<ndim>* perturber[Npertmax]; ///< Array of ptrs to perturbers

};
#endif
