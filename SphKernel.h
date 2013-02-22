// ============================================================================
// SphKernel.h
// ============================================================================


#ifndef _SPH_KERNEL_H_
#define _SPH_KERNEL_H_

#include <math.h>

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

#endif
