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
GaussianKernel::GaussianKernel(int ndimaux, string kernelname)
{
#if !defined(FIXED_DIMENSIONS)
  ndim = ndimaux;
  ndimpr = (FLOAT) ndimaux;
#endif
  kernrange = (FLOAT) 3.0;
  invkernrange = onethird;
  kernrangesqd = (FLOAT) 9.0;
  if (ndim == 1) kernnorm = sqrt(invpi);
  else if (ndim == 2) kernnorm = invpi;
  else if (ndim == 3) kernnorm = invpi*sqrt(invpi);
}



// ============================================================================
// GaussianKernel::~GaussianKernel
// ============================================================================
GaussianKernel::~GaussianKernel()
{
}



