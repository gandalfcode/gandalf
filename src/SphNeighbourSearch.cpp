//=============================================================================
//  SphNeighbourSearch.cpp
//  ..
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


#include <iostream>
#include "SphNeighbourSearch.h"
#include "Sph.h"
#include "Parameters.h"
#include "SphParticle.h"
#include "Debug.h"
#include "InlineFuncs.h"
using namespace std;


//=============================================================================
// SphNeighbourSearch::SphNeighbourSearch
// Empty default constructor
//=============================================================================
template <int ndim>
SphNeighbourSearch<ndim>::SphNeighbourSearch
(FLOAT kernrangeaux,
 DomainBox<ndim> *boxaux,
 SphKernel<ndim> *kernaux,
 CodeTiming *timingaux):
  kernrange(kernrangeaux),
  kernrangesqd(kernrangeaux*kernrangeaux),
  box(boxaux),
  kernp(kernaux),
  timing(timingaux)
{
}



//=============================================================================
// SphNeighbourSearch::~SphNeighbourSearch
// Empty default destructor
//=============================================================================
template <int ndim>
SphNeighbourSearch<ndim>::~SphNeighbourSearch()
{
}



template class SphNeighbourSearch<1>;
template class SphNeighbourSearch<2>;
template class SphNeighbourSearch<3>;
