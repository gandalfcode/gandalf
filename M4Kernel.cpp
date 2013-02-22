// ============================================================================
// M4Kernel.cpp
// ============================================================================


#include <math.h>
#include <iostream>
#include "Constants.h"
#include "Dimensions.h"
#include "SphKernel.h"
using namespace std;



// ============================================================================
// M4Kernel::M4Kernel
// ============================================================================
M4Kernel::M4Kernel(int ndimaux)
{
#if !defined(FIXED_DIMENSIONS)
  ndim = ndimaux;
  ndimpr = (float) ndimaux;
#endif
  kernrange = 2.0;
  invkernrange = 0.5;
  kernrangesqd = 4.0;
  if (ndim == 1) kernnorm = twothirds;
  else if (ndim == 2) kernnorm = invpi*10.0/7.0;
  else if (ndim == 3) kernnorm = invpi;
}



// ============================================================================
// M4Kernel::~M4Kernel
// ============================================================================
M4Kernel::~M4Kernel()
{
}



