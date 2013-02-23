// ============================================================================
// SphKernel.h
// ============================================================================


#ifndef _SPH_KERNEL_H_
#define _SPH_KERNEL_H_

#include <math.h>
#include "Dimensions.h"


// ============================================================================
// Class SphKernel
// ============================================================================
class SphKernel
{
 public:

  virtual float w0(float) = 0;
  virtual float w1(float) = 0;
  virtual float womega(float) = 0;
  virtual float wzeta(float) = 0;
  virtual float wgrav(float) = 0;
  virtual float wpot(float) = 0;

  float kernrange;
  float invkernrange;
  float kernrangesqd;
  float kernnorm;
#if !defined(FIXED_DIMENSIONS)
  int ndim;
  float ndimpr;
#endif

};



// ============================================================================
// Class M4Kernel
// ============================================================================
class M4Kernel: public SphKernel
{
 public:

  M4Kernel(int);
  ~M4Kernel();

  // M4 kernel function prototypes
  // --------------------------------------------------------------------------
  float w0(float);
  float w1(float);
  float womega(float);
  float wzeta(float);
  float wgrav(float);
  float wpot(float);

};


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



// ============================================================================
// Class QuinticKernel
// ============================================================================
class QuinticKernel: public SphKernel
{
 public:

  QuinticKernel(int);
  ~QuinticKernel();

  // M4 kernel function prototypes
  // --------------------------------------------------------------------------
  float w0(float);
  float w1(float);
  float womega(float);
  float wzeta(float);
  float wgrav(float);
  float wpot(float);

};


// ============================================================================
// QuinticKernel::w0
// ============================================================================
inline float QuinticKernel::w0(float s)
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
inline float QuinticKernel::w1(float s)
{
  if (s < 1.0)
    return kernnorm*(-120.0*s + 120.0*pow(s,3) - 50.0*pow(s,4));
  else if (s < 2.0)
    return kernnorm*(75.0 - 420.0*s + 450.0*s*s -
             180.0*pow(s,3) + 25.0*pow(s,4));
  else if (s < 3.0)
    return kernnorm*(-405.0 + 540.0*s - 270.0*s*s +
             60.0*pow(s,3) - 5.0*pow(s,4));
  else
    return 0.0;
}



// ============================================================================
// QuinticKernel::womega
// ============================================================================
inline float QuinticKernel::womega(float s)
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
inline float QuinticKernel::wzeta(float s)
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
inline float QuinticKernel::wgrav(float s)
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
inline float QuinticKernel::wpot(float s)
{
  if (s < 1.0)
    return 1.4 - 0.666666666666666*s*s + 0.3*pow(s,4) - 0.1*pow(s,5);
  else if (s < 2.0)
    return -1.0/(15.0*s) + 1.6 - 1.333333333333333*s*s +
      pow(s,3) - 0.3*pow(s,4) + (1.0/30.0)*pow(s,5);
  else
    return 0.0;
}


#endif
