//=================================================================================================
//  QuinticKernel.cpp
//  Quintic Kernel set-up functions
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
#include <iostream>
#include "Constants.h"
#include "SmoothingKernel.h"
using namespace std;



//=================================================================================================
//  QuinticKernel::QuinticKernel
//=================================================================================================
template <int ndim>
QuinticKernel<ndim>::QuinticKernel(string KernelName):
  SmoothingKernel<ndim>()
{
  kernrange = (FLOAT) 3.0;
  invkernrange = onethird;
  kernrangesqd = (FLOAT) 9.0;
  if (ndim == 1) {
	  kernnorm = (FLOAT) (1.0/120.0);
	  kernnormdrag = 2.0 ;
	  }
  else if (ndim == 2) {
	  kernnorm = invpi*(FLOAT) (7.0/478.0);
	  kernnormdrag = 2868. /2771;
  }
  else if (ndim == 3) {
	  kernnorm = invpi*(FLOAT) (1/120.);
	  kernnormdrag = 5 / 7.;
  }
}



//=================================================================================================
//  QuinticKernel::~QuinticKernel
//=================================================================================================
template <int ndim>
QuinticKernel<ndim>::~QuinticKernel()
{
}



template class QuinticKernel<1>;
template class QuinticKernel<2>;
template class QuinticKernel<3>;
