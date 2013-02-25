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

  virtual FLOAT w0(FLOAT) = 0;
  virtual FLOAT w1(FLOAT) = 0;
  virtual FLOAT womega(FLOAT) = 0;
  virtual FLOAT wzeta(FLOAT) = 0;
  virtual FLOAT wgrav(FLOAT) = 0;
  virtual FLOAT wpot(FLOAT) = 0;

  FLOAT kernrange;
  FLOAT invkernrange;
  FLOAT kernrangesqd;
  FLOAT kernnorm;
#if !defined(FIXED_DIMENSIONS)
  int ndim;
  FLOAT ndimpr;
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
  FLOAT w0(FLOAT);
  FLOAT w1(FLOAT);
  FLOAT womega(FLOAT);
  FLOAT wzeta(FLOAT);
  FLOAT wgrav(FLOAT);
  FLOAT wpot(FLOAT);

};


// ============================================================================
// M4Kernel::w0
// ============================================================================
inline FLOAT M4Kernel::w0(FLOAT s)
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
inline FLOAT M4Kernel::w1(FLOAT s)
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
inline FLOAT M4Kernel::womega(FLOAT s)
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
inline FLOAT M4Kernel::wzeta(FLOAT s)
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
inline FLOAT M4Kernel::wgrav(FLOAT s)
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
inline FLOAT M4Kernel::wpot(FLOAT s)
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
  FLOAT w0(FLOAT);
  FLOAT w1(FLOAT);
  FLOAT womega(FLOAT);
  FLOAT wzeta(FLOAT);
  FLOAT wgrav(FLOAT);
  FLOAT wpot(FLOAT);

};


// ============================================================================
// QuinticKernel::w0
// ============================================================================
inline FLOAT QuinticKernel::w0(FLOAT s)
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
inline FLOAT QuinticKernel::w1(FLOAT s)
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
inline FLOAT QuinticKernel::womega(FLOAT s)
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
inline FLOAT QuinticKernel::wzeta(FLOAT s)
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
inline FLOAT QuinticKernel::wgrav(FLOAT s)
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
inline FLOAT QuinticKernel::wpot(FLOAT s)
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
