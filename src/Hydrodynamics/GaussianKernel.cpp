//=============================================================================
//  GaussianKernel.cpp
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


#include <math.h>
#include <iostream>
#include "Constants.h"
#include "SmoothingKernel.h"
using namespace std;



//=============================================================================
//  GaussianKernel::GaussianKernel
/// GaussianKernel contructor to initialise all kernel constants
//=============================================================================
template <int ndim>
GaussianKernel<ndim>::GaussianKernel(string kernelname):
  SmoothingKernel<ndim>()
{
  kernrange = (FLOAT) 3.0;
  invkernrange = onethird;
  kernrangesqd = (FLOAT) 9.0;
  if (ndim == 1) {
	  kernnorm = sqrt(invpi);
	  kernnormdrag = 2.0 ;
  }
  else if (ndim == 2) {
	  kernnorm = invpi;
	  kernnormdrag = 1.0 ;
  }
  else if (ndim == 3){
	  kernnorm = invpi*sqrt(invpi);
	  kernnormdrag = 2. / 3 ;
  }
}



//=============================================================================
//  GaussianKernel::~GaussianKernel
/// GaussianKernel destructor
//=============================================================================
template <int ndim>
GaussianKernel<ndim>::~GaussianKernel()
{
}



template class GaussianKernel<1>;
template class GaussianKernel<2>;
template class GaussianKernel<3>;

