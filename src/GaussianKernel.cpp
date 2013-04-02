// ============================================================================
// GaussianKernel.cpp
// ============================================================================


#include <math.h>
#include <iostream>
#include "Constants.h"
#include "Dimensions.h"
#include "SphKernel.h"
using namespace std;



// ============================================================================
// GaussianKernel::GaussianKernel
// ============================================================================
template <int ndim>
GaussianKernel<ndim>::GaussianKernel(string kernelname):
  SphKernel<ndim>()
{
//#if !defined(FIXED_DIMENSIONS)
//  ndim = ndimaux;
//  ndimpr = (FLOAT) ndim;
//#endif
  this->kernrange = (FLOAT) 3.0;
  this->invkernrange = onethird;
  this->kernrangesqd = (FLOAT) 9.0;
  if (ndim == 1) this->kernnorm = sqrt(invpi);
  else if (ndim == 2) this->kernnorm = invpi;
  else if (ndim == 3) this->kernnorm = invpi*sqrt(invpi);
}



// ============================================================================
// GaussianKernel::~GaussianKernel
// ============================================================================
template <int ndim>
GaussianKernel<ndim>::~GaussianKernel()
{
}


template class GaussianKernel<1>;
template class GaussianKernel<2>;
template class GaussianKernel<3>;

