// ============================================================================
// SphKernel.h
// ============================================================================


#ifndef _SPH_KERNEL_H_
#define _SPH_KERNEL_H_


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

#endif
