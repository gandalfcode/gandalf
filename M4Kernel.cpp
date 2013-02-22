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



// ============================================================================
// M4Kernel::w0
// ============================================================================
inline float M4Kernel::w0(float s)
{
  if (s < 1.0)
    return kernnorm*(1.0 - 1.5*s*s + 0.75*s*s*s);
  else if (s < 2.0)
    return 0.25*kernnorm*powf(2.0 - s,3.0);
  else
    return 0.0;
}



// ============================================================================
// M4Kernel::w1
// ============================================================================
inline float M4Kernel::w1(float s)
{
  if (s < 1.0)
    return kernnorm*(-3.0*s + 2.25*s*s);
  else if (s < 2.0)
    return -0.75*kernnorm*(2.0 - s)*(2.0 - s);
  else
    return 0.0;
}



// ============================================================================
// M4Kernel::womega
// ============================================================================
inline float M4Kernel::womega(float s)
{
  if (s < 1.0)
    return kernnorm*(-ndimpr + 1.5*(ndimpr + 2.0)*s*s - 
		     0.75*(ndimpr + 3.0)*pow(s,3));
  else if (s < 2.0)
    return kernnorm*(-2.0*ndimpr + 3.0*(ndimpr + 1.0)*s - 1.50*
		     (ndimpr + 2.0)*s*s + 0.25*(ndimpr + 3.0)*pow(s,3));
  else
    return 0.0;
}



// ============================================================================
// M4Kernel::wzeta
// ============================================================================
inline float M4Kernel::wzeta(float s)
{
  if (s < 1.0)
    return 1.4 - 2.0*s*s + 1.5*pow(s,4) - 0.6*pow(s,5);
  else if (s < 2.0)
    return 1.6 - 4.0*s*s + 4.0*pow(s,3) - 1.5*pow(s,4) + 0.2*pow(s,5);
  else
    return 0.0;
}



// ============================================================================
// M4Kernel::wgrav
// ============================================================================
inline float M4Kernel::wgrav(float s)
{
  if (s < 1.0)
    return 1.33333333333333*s - 1.2*pow(s,3) + 0.5*pow(s,4);
  else if (s < 2.0)
    return 2.66666666666667*s - 3.0*s*s + 1.2*pow(s,3) -
      0.1666666666666667*pow(s,4) - 0.06666666666667/(s*s);
  else
    return 0.0;
}



// ============================================================================
// M4Kernel::wpot
// ============================================================================
inline float M4Kernel::wpot(float s)
{
  if (s < 1.0)
    return 1.4 - 0.666666666666666*s*s + 0.3*pow(s,4) - 0.1*pow(s,5);
  else if (s < 2.0)
    return -1.0/(15.0*s) + 1.6 - 1.333333333333333*s*s + 
      pow(s,3) - 0.3*pow(s,4) + (1.0/30.0)*pow(s,5);
  else
    return 0.0;
}
