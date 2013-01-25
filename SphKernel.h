// ============================================================================
// KERNEL.H
// ============================================================================


#ifndef _KERNEL_H_
#define _KERNEL_H_




// ============================================================================
// CLASS KERNEL
// ============================================================================
class SphKernel
{
 public:

  virtual float w0(float) = 0;
  virtual float w1(float) = 0;
  virtual void Setup(int) = 0;

  float kernrange;
  float invkernrange;
  float kernrangesqd;
  float kernnorm;

};



// ============================================================================
// CLASS M4
// ============================================================================
class m4: public SphKernel
{
 public:

  m4();
  ~m4();

  // M4-kernel function prototypes
  // --------------------------------------------------------------------------
  float w0(float);
  float w1(float);
  //float w1_tc(float);
  void Setup(int);

};

#endif
