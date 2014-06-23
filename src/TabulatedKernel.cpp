//=============================================================================
//  TabulatedKernel.cpp
//  Contains functions for creating tabulated kernels for all required 
//  kernel functions and derivatives.
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


#include "SphKernel.h"


//=============================================================================
//  SphKernel::KernelFactory
/// Create and return new kernel object selected in parameters file.
//=============================================================================
template <int ndim>
SphKernel<ndim>* SphKernel<ndim>::KernelFactory(string KernelName) 
{
  if (KernelName == "m4")
    return new M4Kernel<ndim>(KernelName);
  else if (KernelName == "quintic")
    return new QuinticKernel<ndim>(KernelName);
  else if (KernelName == "gaussian")
    return new GaussianKernel<ndim>(KernelName);
  else {
    string message = "Unrecognised kernel: " + KernelName;
    ExceptionHandler::getIstance().raise(message);
  }
  return NULL;
}



//=============================================================================
//  TabulatedKernel::TabulatedKernel
/// Create and return new kernel object selected in parameters file.
//=============================================================================
template <int ndim>
TabulatedKernel<ndim>::TabulatedKernel(string KernelName, int resaux):
  SphKernel<ndim>()
{
  res = resaux;

  kernel = SphKernel<ndim>::KernelFactory (KernelName);

  this->kernrange    = kernel->kernrange;
  this->kernrangesqd = kernel->kernrangesqd;
  this->invkernrange = kernel->invkernrange;
  this->kernnorm     = kernel->kernnorm;
  resinvkernrange    = res/this->kernrange;
  resinvkernrangesqd = res/this->kernrangesqd;

  // Allocate memory
  tableW0        = new FLOAT[res];
  tableW1        = new FLOAT[res];
  tableWomega    = new FLOAT[res];
  tableWzeta     = new FLOAT[res];
  tableWgrav     = new FLOAT[res];
  tableWpot      = new FLOAT[res];
  tableW0_s2     = new FLOAT[res];
  tableWomega_s2 = new FLOAT[res];
  tableWzeta_s2  = new FLOAT[res];
  tableLOS       = new FLOAT[res];

  // Initialize the tables
  initializeTable(tableW0,&SphKernel<ndim>::w0);
  initializeTable(tableW1,&SphKernel<ndim>::w1);
  initializeTable(tableWomega,&SphKernel<ndim>::womega);
  initializeTable(tableWzeta,&SphKernel<ndim>::wzeta);
  initializeTable(tableWgrav,&SphKernel<ndim>::wgrav);
  initializeTable(tableWpot,&SphKernel<ndim>::wpot);
  initializeTableSqd(tableW0_s2,&SphKernel<ndim>::w0);
  initializeTableSqd(tableWomega_s2,&SphKernel<ndim>::womega);
  initializeTableSqd(tableWzeta_s2,&SphKernel<ndim>::wzeta);
  initializeTableLOS();

  // Delete kernel object now that we don't need it anymore
  delete kernel;
}



//=============================================================================
//  TabulatedKernel::initializeTableLOS
/// Create and return new kernel object selected in parameters file.
//=============================================================================
template <int ndim>
void TabulatedKernel<ndim>::initializeTableLOS() 
{
  int i;                                    // Kernel table element counter
  int j;                                    // Integration step counter
  FLOAT dist;                               // ..
  FLOAT impactparametersqd;                 // Kernel impact parameter squared
  FLOAT intstep;                            // ..
  FLOAT position;                           // ..
  FLOAT sum;                                // Integration sum
  FLOAT s;                                  // Distance from the center
  const FLOAT step = kernel->kernrange/res; // Step in the tabulated variable
  const int intsteps = 4000;                // No. of steps per integration

  //---------------------------------------------------------------------------
  for (i=0; i<res; i++) {
    impactparametersqd = pow((FLOAT) i*step,2);
    sum = 0.0;

    // Half-length of the integration path
    dist = sqrt(this->kernrangesqd - impactparametersqd); 
    intstep = dist/intsteps;

    // Now numerically integrate through kernel
    for (j=0; j<intsteps; j++) {
      position = intstep*j;

      // Compute distance from the center
      s = sqrt(position*position + impactparametersqd);
      sum += kernel->w0(s)*intstep;
    }
    
    // Multiply by 2 because we integrated only along half of the path
    tableLOS[i] = 2.0*sum; 
  }
  //---------------------------------------------------------------------------

  return;
}



template class TabulatedKernel<1>;
template class TabulatedKernel<2>;
template class TabulatedKernel<3>;
