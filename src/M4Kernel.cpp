// ============================================================================
// M4Kernel.cpp
// ============================================================================


#include <math.h>
#include <iostream>
#include "Constants.h"
#include "SphKernel.h"
using namespace std;



// ============================================================================
// M4Kernel::M4Kernel
// ============================================================================
template <int ndim>
M4Kernel<ndim>::M4Kernel(string kernelname):
  SphKernel<ndim>()
{
//#if !defined(FIXED_DIMENSIONS)
//  ndim = ndimaux;
//  ndimpr = (FLOAT) ndim;
//#endif
  this->kernrange = (FLOAT) 2.0;
  this->invkernrange = (FLOAT) 0.5;
  this->kernrangesqd = (FLOAT) 4.0;
  if (ndim == 1) this->kernnorm = twothirds;
  else if (ndim == 2) this->kernnorm = invpi*(FLOAT) (10.0/7.0);
  else if (ndim == 3) this->kernnorm = invpi;
}



// ============================================================================
// M4Kernel::~M4Kernel
// ============================================================================
template <int ndim>
M4Kernel<ndim>::~M4Kernel()
{
}



template class M4Kernel<1>;
template class M4Kernel<2>;
template class M4Kernel<3>;
