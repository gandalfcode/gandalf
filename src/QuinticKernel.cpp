// ============================================================================
// QuinticKernel.cpp
// ============================================================================


#include <math.h>
#include <iostream>
#include "Constants.h"
#include "Dimensions.h"
#include "SphKernel.h"
using namespace std;



// ============================================================================
// QuinticKernel::QuinticKernel
// ============================================================================
template <int ndim>
QuinticKernel<ndim>::QuinticKernel(string KernelName):
  SphKernel<ndim>()
{
//#if !defined(FIXED_DIMENSIONS)
//  ndim = ndimaux;
//  ndimpr = (FLOAT) ndim;
//#endif
  this->kernrange = (FLOAT) 3.0;
  this->invkernrange = onethird;
  this->kernrangesqd = (FLOAT) 9.0;
  if (ndim == 1) this->kernnorm = (FLOAT) (1.0/120.0);
  else if (ndim == 2) this->kernnorm = invpi*(FLOAT) (7.0/478.0);
  else if (ndim == 3) this->kernnorm = invpi*(FLOAT) (3.0/359.0);
}



// ============================================================================
// QuinticKernel::~QuinticKernel
// ============================================================================
template <int ndim>
QuinticKernel<ndim>::~QuinticKernel()
{
}



template class QuinticKernel<1>;
template class QuinticKernel<2>;
template class QuinticKernel<3>;



