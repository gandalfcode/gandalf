// ============================================================================
// QuinticKernel.Cpp
// ============================================================================


#include <cmath>
#include <iostream>
#include "Constants.h"
#include "Dimensions.h"
#include "SphKernel.h"
using namespace std;



// ============================================================================
// QuinticKernel::QuinticKernel
// ============================================================================
QuinticKernel::QuinticKernel(int ndimaux)
{
#if !defined(FIXED_DIMENSIONS)
  ndim = ndimaux;
  ndimpr = (float) ndimaux;
#endif
  kernrange = 3.0;
  invkernrange = onethird;
  kernrangesqd = 9.0;
  if (ndim == 1) kernnorm = 1.0/120.0;
  else if (ndim == 2) kernnorm = invpi*7.0/478.0;
  else if (ndim == 3) kernnorm = invpi*3.0/359.0;
}



// ============================================================================
// QuinticKernel::~QuinticKernel
// ============================================================================
QuinticKernel::~QuinticKernel()
{
}



// ============================================================================
// QuinticKernel::w0
// ============================================================================
float QuinticKernel::w0(float s)
{
  if (s < 1.0)
    return kernnorm*(66.0 - 60.0*s*s + 30.0*pow(s,4) - 10.0*pow(s,5));
  else if (s < 2.0)
    return kernnorm*(51.0 + 75.0*s - 210.0*s*s + 150.0*pow(s,3) - 
		     45.0*pow(s,4) + 5.0*pow(s,5));
  else if (s < 3.0)
    return kernnorm*(243.0 - 405*s + 270.0*s*s - 90.0*pow(s,3) +
		     15.0*pow(s,4) - pow(s,5));
  else
    return 0.0;
}



// ============================================================================
// QuinticKernel::w1
// ============================================================================
float QuinticKernel::w1(float s)
{
  if (s < 1.0)
    return kernnorm*(-120.0*s + 120.0*pow(s,3) - 50.0*pow(s,4));
  else if (s < 2.0)
    return kernnorm*(75.0 - 420.0*s + 450.0*s*s - 
		     180.0*pow(s,3) + 25.0*pow(s,4));
  else if (s < 2.0)
    return kernnorm*(-405.0 + 540.0*s - 270.0*s*s + 
		     60.0*pow(s,3) - 5.0*pow(s,4));
  else
    return 0.0;
}



// ============================================================================
// QuinticKernel::womega
// ============================================================================
float QuinticKernel::womega(float s)
{
  if (s < 1.0)
    return kernnorm*(-66.0*ndimpr + 60.0*(ndimpr + 2.0)*s*s - 
		     30.0*(ndimpr + 4.0)*pow(s,4) + 
		     10.0*(ndimpr + 5.0)*pow(s,5));
  else if (s < 2.0)
    return kernnorm*(-51.0*ndimpr - 75.0*(ndimpr + 1.0)*s +
		     210.0*(ndimpr + 2.0)*s*s - 
		     150.0*(ndimpr + 3.0)*pow(s,3) +
		     45.0*(ndimpr + 4.0)*pow(s,4) - 
		     5.0*(ndimpr + 5.0)*pow(s,5));
  else if (s < 3.0)
    return kernnorm*(-243.0*ndimpr + 405.0*(ndimpr + 1.0)*s -
		     270.0*(ndimpr + 2.0)*s*s + 
		     90.0*(ndimpr + 3.0)*pow(s,3) - 
		     15.0*(ndimpr + 4.0)*pow(s,4) + 
		     (ndimpr + 5.0)*pow(s,5));
  else
    return 0.0;
}



// ============================================================================
// QuinticKernel::wzeta
// ============================================================================
float QuinticKernel::wzeta(float s)
{
  if (s < 1.0)
    return 1.4 - 2.0*s*s + 1.5*pow(s,4) - 0.6*pow(s,5);
  else if (s < 2.0)
    return 1.6 - 4.0*s*s + 4.0*pow(s,3) - 1.5*pow(s,4) + 0.2*pow(s,5);
  else
    return 0.0;
}



// ============================================================================
// QuinticKernel::wgrav
// ============================================================================
float QuinticKernel::wgrav(float s)
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
// QuinticKernel::wpot
// ============================================================================
float QuinticKernel::wpot(float s)
{
  if (s < 1.0)
    return 1.4 - 0.666666666666666*s*s + 0.3*pow(s,4) - 0.1*pow(s,5);
  else if (s < 2.0)
    return -1.0/(15.0*s) + 1.6 - 1.333333333333333*s*s + 
      pow(s,3) - 0.3*pow(s,4) + (1.0/30.0)*pow(s,5);
  else
    return 0.0;
}

