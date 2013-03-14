// ============================================================================
// SphKernel.h
// ============================================================================


#ifndef _SPH_KERNEL_H_
#define _SPH_KERNEL_H_

#include <math.h>
#include "Dimensions.h"
#include "Exception.h"
#include <iostream>
using namespace std;

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
  virtual FLOAT wLOS(FLOAT) {};
  virtual ~SphKernel(){};
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

  M4Kernel(int, string);
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
  if (s < (FLOAT) 1.0)
    return kernnorm*((FLOAT) 1.0 - (FLOAT) 1.5*s*s + (FLOAT) 0.75*s*s*s);
  else if (s < (FLOAT) 2.0)
    return (FLOAT) 0.25*kernnorm*pow((FLOAT) 2.0 - s,3);
  else
    return (FLOAT) 0.0;
}



// ============================================================================
// M4Kernel::w1
// ============================================================================
inline FLOAT M4Kernel::w1(FLOAT s)
{
  if (s < (FLOAT) 1.0)
    return kernnorm*(-(FLOAT) 3.0*s + (FLOAT) 2.25*s*s);
  else if (s < (FLOAT) 2.0)
    return -(FLOAT) 0.75*kernnorm*((FLOAT) 2.0 - s)*((FLOAT) 2.0 - s);
  else
    return (FLOAT) 0.0;
}



// ============================================================================
// M4Kernel::womega
// ============================================================================
inline FLOAT M4Kernel::womega(FLOAT s)
{
  if (s < (FLOAT) 1.0)
    return kernnorm*(-ndimpr + (FLOAT) 1.5*(ndimpr + (FLOAT) 2.0)*s*s -
		     (FLOAT) 0.75*(ndimpr + (FLOAT) 3.0)*pow(s,3));
  else if (s < (FLOAT) 2.0)
    return kernnorm*(-(FLOAT) 2.0*ndimpr + 
		     (FLOAT) 3.0*(ndimpr + (FLOAT) 1.0)*s - (FLOAT) 1.50*
		     (ndimpr + (FLOAT) 2.0)*s*s + 
		     (FLOAT) 0.25*(ndimpr + (FLOAT) 3.0)*pow(s,3));
  else
    return (FLOAT) 0.0;
}



// ============================================================================
// M4Kernel::wzeta
// ============================================================================
inline FLOAT M4Kernel::wzeta(FLOAT s)
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

  QuinticKernel(int, string);
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

static SphKernel* KernelFactory (int ndimaux, string KernelName) {
  if (KernelName == "m4")
    return new M4Kernel(ndimaux, KernelName);
  else if (KernelName == "quintic")
    return new QuinticKernel(ndimaux,KernelName);
  else {
    string message = "Unrecognised kernel: " + KernelName;
    ExceptionHandler::getIstance().raise(message);
  }
  return NULL;
}

class TabulatedKernel: public SphKernel {
private:
  SphKernel* kernel;
  int res;
  FLOAT resf;
  FLOAT* tableW0;
  FLOAT* tableW1;
  FLOAT* tableWomega;
  FLOAT* tableWzeta;
  FLOAT* tableWgrav;
  FLOAT* tableWpot;
  FLOAT* tableLOS;

  void initializeTable(FLOAT* table, FLOAT (SphKernel::*function) (FLOAT s)) {
    const FLOAT step = kernel->kernrange/res;
    for (int i=0; i< res; i++) {
      table[i] = (kernel->*function)(step*i);
    }
  }

  // -------------------------------
  // Lookup table for the line of sight integrated kernel
  // ------------------------------
  void initializeTableLOS() {
    const FLOAT step = kernel->kernrange/res; //step in the tabulated variable
    const int intsteps = 4000; //how many steps for each integration
    FLOAT sum;
    FLOAT s; //distance from the center
    for (int i=0; i<res; i++) {
      FLOAT impactparameter = i*step;
      FLOAT impactparametersqd = impactparameter*impactparameter;
      sum=0;
      FLOAT dist = sqrt(kernrangesqd-impactparametersqd); //half-length of the integration path
      FLOAT intstep = dist/intsteps;
      for (int j=0; j<intsteps; j++) {
        FLOAT position = intstep*j;
        s = sqrt(position*position+impactparametersqd); //compute distance from the center
        sum += kernel->w0(s)*intstep;
      }
      tableLOS[i] = 2*sum; //multiply by 2 because we integrated only along half of the path
    }
  };


  FLOAT tableLookup (FLOAT* table, FLOAT s) {
    if (s>kernrange)
      return 0;
    FLOAT indexf = s*resf/kernrange;
    int index = (int)indexf;
    return table[index];
  }

public:
  TabulatedKernel(int ndimaux, string KernelName, int resaux=1000)
    {
    cout << "Using tabulated kernel" << endl;
    res = resaux;
    resf = (int)res;

    kernel = KernelFactory (ndimaux, KernelName);

    kernrange = kernel->kernrange;
    kernrangesqd = kernel->kernrangesqd;
    invkernrange = kernel->invkernrange;
    kernnorm = kernel->kernnorm;
#if !defined (FIXED_DIMENSIONS)
    ndim = kernel->ndim;
    ndimpr = kernel->ndimpr;
#endif

    //allocate memory
    tableW0 = new FLOAT[res];
    tableW1 = new FLOAT[res];
    tableWomega = new FLOAT[res];
    tableWzeta = new FLOAT[res];
    tableWgrav = new FLOAT[res];
    tableWpot = new FLOAT[res];
    tableLOS = new FLOAT[res];


    //initialize the tables
    initializeTable(tableW0,&SphKernel::w0);
    initializeTable(tableW1,&SphKernel::w1);
    initializeTable(tableWomega,&SphKernel::womega);
    initializeTable(tableWzeta,&SphKernel::wzeta);
    initializeTable(tableWgrav,&SphKernel::wgrav);
    initializeTable(tableWpot,&SphKernel::wpot);
    initializeTableLOS();

    //deallocates the kernel now that we don't need it anymore
    delete kernel;
  }

  FLOAT w0(FLOAT s);
  FLOAT w1(FLOAT s);
  FLOAT womega(FLOAT s);
  FLOAT wzeta(FLOAT s);
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
    delete[] tableLOS;
  }


};

inline FLOAT TabulatedKernel::w0 (FLOAT s) {
  return tableLookup(tableW0, s);
}

inline FLOAT TabulatedKernel::w1 (FLOAT s) {
  return tableLookup(tableW1, s);
}

inline FLOAT TabulatedKernel::womega (FLOAT s) {
  return tableLookup(tableWomega, s);
}

inline FLOAT TabulatedKernel::wzeta (FLOAT s) {
  return tableLookup(tableWzeta, s);
}

inline FLOAT TabulatedKernel::wgrav (FLOAT s) {
  return tableLookup(tableWgrav, s);
}

inline FLOAT TabulatedKernel::wpot (FLOAT s) {
  return tableLookup(tableWpot, s);
}

inline FLOAT TabulatedKernel::wLOS (FLOAT s) {
  return tableLookup(tableLOS, s);
}


#endif
