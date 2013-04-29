//=============================================================================
//  GaussianKernel.cpp
//=============================================================================


#include <math.h>
#include <iostream>
#include "Constants.h"
#include "SphKernel.h"
using namespace std;



//=============================================================================
//  GaussianKernel::GaussianKernel
/// GaussianKernel contructor to initialise all kernel constants
//=============================================================================
template <int ndim>
GaussianKernel<ndim>::GaussianKernel(string kernelname):
  SphKernel<ndim>()
{
  this->kernrange = (FLOAT) 3.0;
  this->invkernrange = onethird;
  this->kernrangesqd = (FLOAT) 9.0;
  if (ndim == 1) this->kernnorm = sqrt(invpi);
  else if (ndim == 2) this->kernnorm = invpi;
  else if (ndim == 3) this->kernnorm = invpi*sqrt(invpi);
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

