// ============================================================================
// SphKernel.h
// ============================================================================


#ifndef _SPH_KERNEL_H_
#define _SPH_KERNEL_H_

#include <math.h>
#include "Exception.h"
#include <iostream>
using namespace std;

// ============================================================================
// Class SphKernel
// ============================================================================
template <int ndim>
class SphKernel
{
 public:


  SphKernel(): ndimpr(ndim) {};

  static SphKernel<ndim>* KernelFactory (string KernelName);

  virtual FLOAT w0(FLOAT) = 0;
 //For the versions using the squared distance, the default behaviour
 //is to call the standard one with the square root
  virtual inline FLOAT w0_s2(FLOAT s) {return this->w0(sqrt(s));};
  virtual FLOAT w1(FLOAT) = 0;
  virtual FLOAT womega(FLOAT) = 0;
  virtual inline FLOAT womega_s2(FLOAT s) {return this->womega(sqrt(s));};
  virtual FLOAT wzeta(FLOAT) = 0;
  virtual inline FLOAT wzeta_s2(FLOAT s) {return this->wzeta(sqrt(s));};
  virtual FLOAT wgrav(FLOAT) = 0;
  virtual FLOAT wpot(FLOAT) = 0;
  virtual FLOAT wLOS(FLOAT) {
    //We do not provide a default behaviour
    string message = "Using a non-tabulated kernel, cannot use the column integrated kernel!";
    ExceptionHandler::getIstance().raise(message);
    return 0;
  };
  virtual ~SphKernel(){};
  FLOAT kernrange;
  FLOAT invkernrange;
  FLOAT kernrangesqd;
  FLOAT kernnorm;
//#if !defined(FIXED_DIMENSIONS)
//  int ndim;
  FLOAT ndimpr;
//#endif

};



// ============================================================================
// Class M4Kernel
// ============================================================================
template <int ndim>
class M4Kernel: public SphKernel<ndim>
{

 public:

  M4Kernel(string);
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
template <int ndim>
inline FLOAT M4Kernel<ndim>::w0(FLOAT s)
{
  if (s < (FLOAT) 1.0)
    return (this->kernnorm)*((FLOAT) 1.0 - (FLOAT) 1.5*s*s + (FLOAT) 0.75*s*s*s);
  else if (s < (FLOAT) 2.0)
    return (FLOAT) 0.25*(this->kernnorm)*pow((FLOAT) 2.0 - s,3);
  else
    return (FLOAT) 0.0;
}



// ============================================================================
// M4Kernel::w1
// ============================================================================
template <int ndim>
inline FLOAT M4Kernel<ndim>::w1(FLOAT s)
{
  if (s < (FLOAT) 1.0)
    return (this->kernnorm)*(-(FLOAT) 3.0*s + (FLOAT) 2.25*s*s);
  else if (s < (FLOAT) 2.0)
    return -(FLOAT) 0.75*(this->kernnorm)*((FLOAT) 2.0 - s)*((FLOAT) 2.0 - s);
  else
    return (FLOAT) 0.0;
}



// ============================================================================
// M4Kernel::womega
// ============================================================================
template <int ndim>
inline FLOAT M4Kernel<ndim>::womega(FLOAT s)
{
  if (s < (FLOAT) 1.0)
    return (this->kernnorm)*(-(this->ndimpr) + (FLOAT) 1.5*((this->ndimpr) + (FLOAT) 2.0)*s*s -
		     (FLOAT) 0.75*((this->ndimpr) + (FLOAT) 3.0)*pow(s,3));
  else if (s < (FLOAT) 2.0)
    return (this->kernnorm)*(-(FLOAT) 2.0*(this->ndimpr) +
		     (FLOAT) 3.0*((this->ndimpr) + (FLOAT) 1.0)*s - (FLOAT) 1.50*
		     ((this->ndimpr) + (FLOAT) 2.0)*s*s +
		     (FLOAT) 0.25*((this->ndimpr) + (FLOAT) 3.0)*pow(s,3));
  else
    return (FLOAT) 0.0;
}



// ============================================================================
// M4Kernel::wzeta
// ============================================================================
template <int ndim>
inline FLOAT M4Kernel<ndim>::wzeta(FLOAT s)
{
  if (s < (FLOAT) 1.0)
    return (FLOAT) 1.4 - (FLOAT) 2.0*s*s + (FLOAT) 1.5*pow(s,4) 
      - (FLOAT) 0.6*pow(s,5);
  else if (s < (FLOAT) 2.0)
    return (FLOAT) 1.6 - (FLOAT) 4.0*s*s + (FLOAT) 4.0*pow(s,3) 
      - (FLOAT) 1.5*pow(s,4) + (FLOAT) 0.2*pow(s,5);
  else
    return 0.0;
}



// ============================================================================
// M4Kernel::wgrav
// ============================================================================
template <int ndim>
inline FLOAT M4Kernel<ndim>::wgrav(FLOAT s)
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
template <int ndim>
inline FLOAT M4Kernel<ndim>::wpot(FLOAT s)
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
template <int ndim>
class QuinticKernel: public SphKernel<ndim>
{
 public:

  QuinticKernel(string);
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
template <int ndim>
inline FLOAT QuinticKernel<ndim>::w0(FLOAT s)
{
  if (s < 1.0)
    return (this->kernnorm)*(66.0 - 60.0*s*s + 30.0*pow(s,4) - 10.0*pow(s,5));
  else if (s < 2.0)
    return (this->kernnorm)*(51.0 + 75.0*s - 210.0*s*s + 150.0*pow(s,3) -
             45.0*pow(s,4) + 5.0*pow(s,5));
  else if (s < 3.0)
    return (this->kernnorm)*(243.0 - 405*s + 270.0*s*s - 90.0*pow(s,3) +
             15.0*pow(s,4) - pow(s,5));
  else
    return 0.0;
}



// ============================================================================
// QuinticKernel::w1
// ============================================================================
template <int ndim>
inline FLOAT QuinticKernel<ndim>::w1(FLOAT s)
{
  if (s < 1.0)
    return (this->kernnorm)*(-120.0*s + 120.0*pow(s,3) - 50.0*pow(s,4));
  else if (s < 2.0)
    return (this->kernnorm)*(75.0 - 420.0*s + 450.0*s*s -
             180.0*pow(s,3) + 25.0*pow(s,4));
  else if (s < 3.0)
    return (this->kernnorm)*(-405.0 + 540.0*s - 270.0*s*s +
             60.0*pow(s,3) - 5.0*pow(s,4));
  else
    return 0.0;
}



// ============================================================================
// QuinticKernel::womega
// ============================================================================
template <int ndim>
inline FLOAT QuinticKernel<ndim>::womega(FLOAT s)
{
  if (s < 1.0)
    return (this->kernnorm)*(-66.0*(this->ndimpr) + 60.0*((this->ndimpr) + 2.0)*s*s -
             30.0*((this->ndimpr) + 4.0)*pow(s,4) +
             10.0*((this->ndimpr) + 5.0)*pow(s,5));
  else if (s < 2.0)
    return (this->kernnorm)*(-51.0*(this->ndimpr) - 75.0*((this->ndimpr) + 1.0)*s +
             210.0*((this->ndimpr) + 2.0)*s*s -
             150.0*((this->ndimpr) + 3.0)*pow(s,3) +
             45.0*((this->ndimpr) + 4.0)*pow(s,4) -
             5.0*((this->ndimpr) + 5.0)*pow(s,5));
  else if (s < 3.0)
    return (this->kernnorm)*(-243.0*(this->ndimpr) + 405.0*((this->ndimpr) + 1.0)*s -
             270.0*((this->ndimpr) + 2.0)*s*s +
             90.0*((this->ndimpr) + 3.0)*pow(s,3) -
             15.0*((this->ndimpr) + 4.0)*pow(s,4) +
             ((this->ndimpr) + 5.0)*pow(s,5));
  else
    return 0.0;
}



// ============================================================================
// QuinticKernel::wzeta
// ============================================================================
template <int ndim>
inline FLOAT QuinticKernel<ndim>::wzeta(FLOAT s)
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
template <int ndim>
inline FLOAT QuinticKernel<ndim>::wgrav(FLOAT s)
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
template <int ndim>
inline FLOAT QuinticKernel<ndim>::wpot(FLOAT s)
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
// Class GaussianKernel
// ============================================================================
template <int ndim>
class GaussianKernel: public SphKernel<ndim>
{
 public:

  GaussianKernel(string);
  ~GaussianKernel();

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
// GaussianKernel::w0
// ============================================================================
template <int ndim>
inline FLOAT GaussianKernel<ndim>::w0(FLOAT s)
{
  if (s < (FLOAT) 3.0)
    return (this->kernnorm)*exp(-s*s);
  else
    return (FLOAT) 0.0;
}



// ============================================================================
// GaussianKernel::w1
// ============================================================================
template <int ndim>
inline FLOAT GaussianKernel<ndim>::w1(FLOAT s)
{
  if (s < (FLOAT) 3.0)
    return -(FLOAT) 2.0*(this->kernnorm)*s*exp(-s*s);
  else
    return (FLOAT) 0.0;
}



// ============================================================================
// GaussianKernel::womega
// ============================================================================
template <int ndim>
inline FLOAT GaussianKernel<ndim>::womega(FLOAT s)
{
  if (s < (FLOAT) 3.0)
    return (this->kernnorm)*((FLOAT) 2.0*s*exp(-s*s) - (this->ndimpr)*exp(-s*s));
  else
    return (FLOAT) 0.0;
}



// ============================================================================
// GaussianKernel::wzeta
// ============================================================================
template <int ndim>
inline FLOAT GaussianKernel<ndim>::wzeta(FLOAT s)
{
  if (s < (FLOAT) 3.0)
    return 0.0;
  else
    return 0.0;
}



// ============================================================================
// GaussianKernel::wgrav
// ============================================================================
template <int ndim>
inline FLOAT GaussianKernel<ndim>::wgrav(FLOAT s)
{
  if (s < (FLOAT) 3.0)
    return 0.0;
  else
    return 0.0;
}



// ============================================================================
// GaussianKernel::wpot
// ============================================================================
template <int ndim>
inline FLOAT GaussianKernel<ndim>::wpot(FLOAT s)
{
  if (s < (FLOAT) 3.0)
    return 0.0;
  else
    return 0.0;
}

template <int ndim>
class TabulatedKernel: public SphKernel<ndim> {
private:
  SphKernel<ndim>* kernel;
  int res;
  FLOAT resinvkernrange;
  FLOAT resinvkernrangesqd;
  FLOAT* tableW0;
  FLOAT* tableW1;
  FLOAT* tableWomega;
  FLOAT* tableWzeta;
  FLOAT* tableWgrav;
  FLOAT* tableWpot;
  FLOAT* tableW0_s2;
  FLOAT* tableWomega_s2;
  FLOAT* tableWzeta_s2;
  FLOAT* tableLOS;

  void initializeTable(FLOAT* table, FLOAT (SphKernel<ndim>::*function) (FLOAT s)) {
    const FLOAT step = kernel->kernrange/res;
    for (int i=0; i< res; i++) {
      table[i] = (kernel->*function)(step*i);
    }
  }

  void initializeTableSqd (FLOAT* table, FLOAT (SphKernel<ndim>::*function) (FLOAT s)) {
    const FLOAT step = kernel->kernrangesqd/res;
    for (int i=0; i< res; i++) {
      table[i] = (kernel->*function)(sqrt(step*i));
    }
  }

  void initializeTableLOS();


  FLOAT tableLookup (FLOAT* table, FLOAT s) {
    if (s>=(this->kernrange))
      return 0;
    FLOAT indexf = s*resinvkernrange;
    int index = (int)indexf;
    return table[index];
  }

  FLOAT tableLookupSqd(FLOAT* table, FLOAT s) {
    if (s>=(this->kernrangesqd))
      return 0;
    FLOAT indexf = s*resinvkernrangesqd;
    int index = (int)indexf;
    return table[index];
  }

public:
  TabulatedKernel(string KernelName, int resaux=1000);

  FLOAT w0(FLOAT s);
  FLOAT w0_s2(FLOAT s);
  FLOAT w1(FLOAT s);
  FLOAT womega(FLOAT s);
  FLOAT womega_s2(FLOAT s);
  FLOAT wzeta(FLOAT s);
  FLOAT wzeta_s2(FLOAT s);
  FLOAT wgrav(FLOAT s);
  FLOAT wpot(FLOAT s);
  FLOAT wLOS(FLOAT s);


  ~TabulatedKernel() {
    delete[] tableW0;
    delete[] tableW1;
    delete[] tableWomega;
    delete[] tableWzeta;
    delete[] tableWgrav;
    delete[] tableWpot;
    delete[] tableW0_s2;
    delete[] tableWomega_s2;
    delete[] tableWzeta_s2;
    delete[] tableLOS;
  }


};

template <int ndim>
inline FLOAT TabulatedKernel<ndim>::w0 (FLOAT s) {
  return tableLookup(tableW0, s);
}

template <int ndim>
inline FLOAT TabulatedKernel<ndim>::w0_s2 (FLOAT s2) {
  return tableLookupSqd(tableW0_s2, s2);
}

template <int ndim>
inline FLOAT TabulatedKernel<ndim>::w1 (FLOAT s) {
  return tableLookup(tableW1, s);
}

template <int ndim>
inline FLOAT TabulatedKernel<ndim>::womega (FLOAT s) {
  return tableLookup(tableWomega, s);
}

template <int ndim>
inline FLOAT TabulatedKernel<ndim>::womega_s2 (FLOAT s2) {
  return tableLookupSqd(tableWomega_s2, s2);
}

template <int ndim>
inline FLOAT TabulatedKernel<ndim>::wzeta (FLOAT s) {
  return tableLookup(tableWzeta, s);
}

template <int ndim>
inline FLOAT TabulatedKernel<ndim>::wzeta_s2 (FLOAT s2) {
  return tableLookupSqd(tableWzeta_s2, s2);
}

template <int ndim>
inline FLOAT TabulatedKernel<ndim>::wgrav (FLOAT s) {
  return tableLookup(tableWgrav, s);
}

template <int ndim>
inline FLOAT TabulatedKernel<ndim>::wpot (FLOAT s) {
  return tableLookup(tableWpot, s);
}

template <int ndim>
inline FLOAT TabulatedKernel<ndim>::wLOS (FLOAT s) {
  return tableLookup(tableLOS, s);
}

#endif
