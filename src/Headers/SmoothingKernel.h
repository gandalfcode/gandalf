//=================================================================================================
//  SmoothingKernel.h
//  Contains all routines for computing kernel functions, calculatiung and
//  storing kernel tables for quick look-up of values.
//
//  This file is part of GANDALF :
//  Graphical Astrophysics code for N-body Dynamics And Lagrangian Fluids
//  https://github.com/gandalfcode/gandalf
//  Contact : gandalfcode@gmail.com
//
//  Copyright (C) 2013  D. A. Hubber, G. Rosotti
//
//  GANDALF is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 2 of the License, or
//  (at your option) any later version.
//
//  GANDALF is distributed in the hope that it will be useful, but
//  WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  General Public License (http://www.gnu.org/licenses) for more details.
//=================================================================================================


#ifndef _SMOOTHING_KERNEL_H_
#define _SMOOTHING_KERNEL_H_


#include <iostream>
#include <math.h>
#include "Exception.h"
#include "Precision.h"
using namespace std;



//=================================================================================================
//  Class SmoothingKernel
/// \brief   Parent class for all SPH kernel options.
/// \details Contains all (virtual) SPH kernel functions.  All functions
///          are fully defined in inheritated child class implementations
///          (i.e. M4, Quintic or Gaussian kernels).
/// \author  D. A. Hubber, G. Rosotti
/// \date    03/04/2013
//=================================================================================================
template <int ndim>
class SmoothingKernel
{
 public:

  SmoothingKernel() {};
  virtual ~SmoothingKernel(){};

  static SmoothingKernel<ndim>* KernelFactory (string KernelName);

  // Main kernel function and associated derivative and integrated forms.
  //---------------------------------------------------------------------------
  virtual FLOAT w0(const FLOAT) = 0;
  virtual FLOAT w1(const FLOAT) = 0;
  virtual FLOAT womega(const FLOAT) = 0;
  virtual FLOAT wzeta(const FLOAT) = 0;
  virtual FLOAT wgrav(const FLOAT) = 0;
  virtual FLOAT wpot(const FLOAT) = 0;
  virtual FLOAT wdrag(const FLOAT q) {
	  return kernnormdrag * q*q * w0(q) ;
  }
  virtual FLOAT wLOS(const FLOAT) {
    //We do not provide a default behaviour
    string message = "Using a non-tabulated kernel, cannot use the column integrated kernel!";
    ExceptionHandler::getIstance().raise(message);
    return 0;
  };


  // For the versions using the squared distance, the default behaviour
  // is to call the standard one with the square root
  //---------------------------------------------------------------------------
  virtual inline FLOAT w0_s2(const FLOAT s) {return w0(sqrt(s));};
  virtual inline FLOAT womega_s2(const FLOAT s) {return womega(sqrt(s));};
  virtual inline FLOAT wzeta_s2(const FLOAT s) {return wzeta(sqrt(s));};


  // Kernel variables
  //---------------------------------------------------------------------------
  FLOAT kernrange;                  ///< Maximum extent of kernel
  FLOAT invkernrange;               ///< 1/kernrange
  FLOAT kernrangesqd;               ///< kernrange^2
  FLOAT kernnorm;                   ///< Kernel normalisation constant
  FLOAT kernnormdrag;               ///< Normalization factor for drag kernel
};



//=================================================================================================
//  Class M4Kernel
/// \brief   M4 kernel class including functions
/// \details Contains all SPH kernel function definitions for M4 kernel.
/// \author  D. A. Hubber, G. Rosotti
/// \date    03/04/2013
//=================================================================================================
template <int ndim>
class M4Kernel: public SmoothingKernel<ndim>
{
	using SmoothingKernel<ndim>::kernnorm;                 ///< Kernel normalisation constant
	using SmoothingKernel<ndim>::kernnormdrag;             ///< Normalization factor for drag kernel
public:
	using SmoothingKernel<ndim>::kernrange;                ///< Maximum extent of kernel
	using SmoothingKernel<ndim>::invkernrange;             ///< 1/kernrange
	using SmoothingKernel<ndim>::kernrangesqd;             ///< kernrange^2


  M4Kernel(string);
  ~M4Kernel();

  // M4 kernel function prototypes
  //---------------------------------------------------------------------------
  FLOAT w0(const FLOAT);
  FLOAT w1(const FLOAT);
  FLOAT womega(const FLOAT);
  FLOAT wzeta(const FLOAT);
  FLOAT wgrav(const FLOAT);
  FLOAT wpot(const FLOAT);
};



//=================================================================================================
//  M4Kernel<ndim>::w0
/// Main SPH smoothing kernel function, $W(s=r/h)$, for M4 kernel.
//=================================================================================================
template <int ndim>
inline FLOAT M4Kernel<ndim>::w0(const FLOAT s)  ///< [in] Kernel parameter, r/h
{
  if (s < (FLOAT) 1.0)
    return (kernnorm)*((FLOAT) 1.0 - (FLOAT) 1.5*s*s + (FLOAT) 0.75*s*s*s);
  else if (s < (FLOAT) 2.0)
    return (FLOAT) 0.25*(kernnorm)*pow((FLOAT) 2.0 - s,3);
  else
    return (FLOAT) 0.0;
}



//=================================================================================================
//  M4Kernel<ndim>::w1
/// First spatial derivative of main smoothing kernel, dWdr, for M4 kernel.
//=================================================================================================
template <int ndim>
inline FLOAT M4Kernel<ndim>::w1(const FLOAT s)  ///< [in] Kernel parameter, r/h
{
  if (s < (FLOAT) 1.0)
    return (kernnorm)*(-(FLOAT) 3.0*s + (FLOAT) 2.25*s*s);
  else if (s < (FLOAT) 2.0)
    return -(FLOAT) 0.75*(kernnorm)*((FLOAT) 2.0 - s)*((FLOAT) 2.0 - s);
  else
    return (FLOAT) 0.0;
}



//=================================================================================================
//  M4Kernel<ndim>::womega
/// Partial derivative of kernel with respect to smoothing length, $dW/dh$,
/// to compute omega/f correction term for M4 kernel.
//=================================================================================================
template <int ndim>
inline FLOAT M4Kernel<ndim>::womega(const FLOAT s)  ///< [in] Kernel parameter, r/h
{
  if (s < (FLOAT) 1.0)
    return (kernnorm)*
      (-(ndim) + (FLOAT) 1.5*((ndim) + (FLOAT) 2.0)*s*s -
       (FLOAT) 0.75*((ndim) + (FLOAT) 3.0)*pow(s,3));
  else if (s < (FLOAT) 2.0)
    return (kernnorm)*
      (-(FLOAT) 2.0*(ndim) +
       (FLOAT) 3.0*((ndim) + (FLOAT) 1.0)*s - (FLOAT) 1.50*
       ((ndim) + (FLOAT) 2.0)*s*s +
       (FLOAT) 0.25*((ndim) + (FLOAT) 3.0)*pow(s,3));
  else
    return (FLOAT) 0.0;
}



//=================================================================================================
//  M4Kernel<ndim>::wzeta
/// Partial derivative of potential kernel with respect to h.  Used to compute
/// gravitational zeta correction term when using grad-h SPH and M4 kernel.
//=================================================================================================
template <int ndim>
inline FLOAT M4Kernel<ndim>::wzeta(const FLOAT s)  ///< [in] Kernel parameter, r/h
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



//=================================================================================================
//  M4Kernel<ndim>::wgrav
/// Volume intergated kernel function to compute kernel-softened gravitational
/// force for M4 kernel.
//=================================================================================================
template <int ndim>
inline FLOAT M4Kernel<ndim>::wgrav(const FLOAT s)  ///< [in] Kernel parameter, r/h
{
  if (s < (FLOAT) 1.0)
    return (FLOAT) 1.333333333333333333333*s - (FLOAT) 1.2*pow(s,3) + (FLOAT) 0.5*pow(s,4);
  else if (s < (FLOAT) 2.0)
    return (FLOAT) 2.6666666666666666667*s - (FLOAT) 3.0*s*s + (FLOAT) 1.2*pow(s,3) -
      (FLOAT) 0.166666666666666666667*pow(s,4) - (FLOAT) 0.06666666666666666667/(s*s);
  else
    return (FLOAT) 1.0/(s*s);
}



//=================================================================================================
//  M4Kernel<ndim>::wpot
/// Volume integrated kernel function for computing softened gravitational
/// potential using M4 kernel.
//=================================================================================================
template <int ndim>
inline FLOAT M4Kernel<ndim>::wpot(const FLOAT s)  ///< [in] Kernel parameter, r/h
{
  if (s < (FLOAT) 1.0)
    return (FLOAT) 1.4 - (FLOAT) 0.666666666666666666666666*s*s +
      (FLOAT) 0.3*pow(s,4) - (FLOAT) 0.1*pow(s,5);
  else if (s < (FLOAT) 2.0)
    return -(FLOAT) 1.0/((FLOAT) 15.0*s) + (FLOAT) 1.6 - (FLOAT) 1.33333333333333333333333333*s*s
      + pow(s,3) - (FLOAT) 0.3*pow(s,4) + (FLOAT) (1.0/30.0)*pow(s,5);
  else
    return (FLOAT) 1.0/s;
}



//=================================================================================================
//  Class QuinticKernel
/// \brief   Quintic kernel class including functions
/// \details Contains all SPH kernel function definitions for Quintic kernel.
/// \author  D. A. Hubber, G. Rosotti
/// \date    03/04/2013
//=================================================================================================
template <int ndim>
class QuinticKernel: public SmoothingKernel<ndim>
{
	using SmoothingKernel<ndim>::kernnorm;                 ///< Kernel normalisation constant
	using SmoothingKernel<ndim>::kernnormdrag;             ///< Normalization factor for drag kernel
public:
	using SmoothingKernel<ndim>::kernrange;                ///< Maximum extent of kernel
	using SmoothingKernel<ndim>::invkernrange;             ///< 1/kernrange
	using SmoothingKernel<ndim>::kernrangesqd;             ///< kernrange^2

  QuinticKernel(string);
  ~QuinticKernel();

  // Quintic kernel function prototypes
  //---------------------------------------------------------------------------
  FLOAT w0(const FLOAT);
  FLOAT w1(const FLOAT);
  FLOAT womega(const FLOAT);
  FLOAT wzeta(const FLOAT);
  FLOAT wgrav(const FLOAT);
  FLOAT wpot(const FLOAT);

};



//=================================================================================================
//  QuinticKernel<ndim>::w0
/// Main SPH smoothing kernel function, $W(s=r/h)$, for quintic kernel.
//=================================================================================================
template <int ndim>
inline FLOAT QuinticKernel<ndim>::w0(const FLOAT s)
{
  if (s < 1.0)
    return (kernnorm)*(66.0 - 60.0*s*s + 30.0*pow(s,4) - 10.0*pow(s,5));
  else if (s < 2.0)
    return (kernnorm)*(51.0 + 75.0*s - 210.0*s*s + 150.0*pow(s,3) -
             45.0*pow(s,4) + 5.0*pow(s,5));
  else if (s < 3.0)
    return (kernnorm)*(243.0 - 405*s + 270.0*s*s - 90.0*pow(s,3) +
             15.0*pow(s,4) - pow(s,5));
  else
    return 0.0;
}



//=================================================================================================
//  QuinticKernel<ndim>::w1
/// First spatial derivative of smoothing kernel, dWdr, for Quintic kernel.
//=================================================================================================
template <int ndim>
inline FLOAT QuinticKernel<ndim>::w1(const FLOAT s)
{
  if (s < 1.0)
    return (kernnorm)*(-120.0*s + 120.0*pow(s,3) - 50.0*pow(s,4));
  else if (s < 2.0)
    return (kernnorm)*(75.0 - 420.0*s + 450.0*s*s - 180.0*pow(s,3) + 25.0*pow(s,4));
  else if (s < 3.0)
    return (kernnorm)*(-405.0 + 540.0*s - 270.0*s*s + 60.0*pow(s,3) - 5.0*pow(s,4));
  else
    return 0.0;
}



//=================================================================================================
//  QuinticKernel<ndim>::womega
/// Derivative of main kernel function w.r.t the smoothing length.
//=================================================================================================
template <int ndim>
inline FLOAT QuinticKernel<ndim>::womega(const FLOAT s)
{
  if (s < 1.0)
    return (kernnorm)*
      (-66.0*(ndim) + 60.0*((ndim) + 2.0)*s*s -
       30.0*((ndim) + 4.0)*pow(s,4) + 10.0*((ndim) + 5.0)*pow(s,5));
  else if (s < 2.0)
    return (kernnorm)*
      (-51.0*(ndim) - 75.0*((ndim) + 1.0)*s + 210.0*((ndim) + 2.0)*s*s -
       150.0*((ndim) + 3.0)*pow(s,3) + 45.0*((ndim) + 4.0)*pow(s,4) -
       5.0*((ndim) + 5.0)*pow(s,5));
  else if (s < 3.0)
    return (kernnorm)*
      (-243.0*(ndim) + 405.0*((ndim) + 1.0)*s - 270.0*((ndim) + 2.0)*s*s +
       90.0*((ndim) + 3.0)*pow(s,3) - 15.0*((ndim) + 4.0)*pow(s,4) +
       ((ndim) + 5.0)*pow(s,5));
  else
    return 0.0;
}



//=================================================================================================
//  QuinticKernel<ndim>::wzeta
/// Derivative of potential kernel w.r.t. smoothing length.
//=================================================================================================
template <int ndim>
inline FLOAT QuinticKernel<ndim>::wzeta(const FLOAT s)
{
  if (s < (FLOAT) 1.0)
    return (FLOAT) 33.0*s*s - (FLOAT) 15.0*pow(s,4) + (FLOAT) 5.0*pow(s,6) -
      (FLOAT) 1.42857142857*pow(s,7) - (FLOAT) 34.14285714;
  else if (s < 2.0)
    return (FLOAT) 25.5*s*s + (FLOAT) 25.0*pow(s,3) - (FLOAT) 52.5*pow(s,4)
      + (FLOAT) 30.0*pow(s,5) - (FLOAT) 7.5*pow(s,6)
      + (FLOAT) 0.7142857143*pow(s,7) - 33.785714286;
  else if (s < (FLOAT) 3.0)
    return (FLOAT) 121.5*s*s - (FLOAT) 135.0*pow(s,3) +
      (FLOAT) 67.5*pow(s,4) - (FLOAT) 18.0*pow(s,5) +
      (FLOAT) 2.5*pow(s,6) - (FLOAT) 0.142857143*pow(s,7) -
      (FLOAT) 52.07142857;
  else
    return 0.0;
}



//=================================================================================================
//  QuinticKernel<ndim>::wgrav
/// Gravitational force kernel.
//=================================================================================================
template <int ndim>
inline FLOAT QuinticKernel<ndim>::wgrav(const FLOAT s)
{
  if (s < 1.0)
    return (12.0/359.0)*(22.0*s - 12.0*pow(s,3) + (30.0/7.0)*pow(s,5) - (5.0/4.0)*pow(s,6));
  else if (s < 2.0)
    return (12.0/359.0)*(17.0*s + (75.0/4.0)*s*s - 42.0*pow(s,3) +
      25.0*pow(s,4) - (45.0/7.0)*pow(s,5) + (5.0/8.0)*pow(s,6) + (5.0/56.0)/(s*s));
  else if (s < 3.0)
    return (12.0/359.0)*(81.0*s - (405.0/4.0)*pow(s,2) + 54.0*pow(s,3) -
      15.0*pow(s,4) + (15.0/7.0)*pow(s,5) - (1.0/8.0)*pow(s,6) - (507.0/56.0)/(s*s));
  else
    return 1.0/(s*s);
}



//=================================================================================================
//  QuinticKernel<ndim>::wpot
/// Gravitational potential kernel
//=================================================================================================
template <int ndim>
inline FLOAT QuinticKernel<ndim>::wpot(const FLOAT s)
{
  if (s < 1.0)
    return (12.0/359.0)*(-11.0*s*s + 3.0*pow(s,4) - (5.0/7.0)*pow(s,6) +
      (5.0/28.0)*pow(s,7) + (478.0/14.0));
  else if (s < 2.0)
    return (12.0/359.0)*(-(17.0/2.0)*s*s - (25.0/4.0)*pow(s,3) + (21.0/2.0)*pow(s,4) -
      5.0*pow(s,5) + (15.0/14.0)*pow(s,6) - (5.0/56.0)*pow(s,7) + (473.0/14.0) + (5.0/56.0)/s);
  else if (s < 3.0)
    return (12.0/359.0)*(-(81.0/2.0)*s*s + (135.0/4.0)*pow(s,3) - (27.0/2.0)*pow(s,4) +
      3.0*pow(s,5) - (5.0/14.0)*pow(s,6) + (1.0/56.0)*pow(s,7) + (729.0/14.0) - (507.0/56.0)/s);
  else
    return 1.0/s;
}



//=================================================================================================
//  Class GaussianKernel
/// \brief   Gaussian kernel class including functions
/// \details Contains all SPH kernel function definitions for Gaussian kernel.
///          Gaussian kernel formally extends to infinity, so for
///          computational efficiency, this kernel only extends to 3h.
/// \author  D. A. Hubber, G. Rosotti
/// \date    03/04/2013
//=================================================================================================
template <int ndim>
class GaussianKernel: public SmoothingKernel<ndim>
{
	using SmoothingKernel<ndim>::kernnorm;                 ///< Kernel normalisation constant
	using SmoothingKernel<ndim>::kernnormdrag;             ///< Normalization factor for drag kernel
public:
	using SmoothingKernel<ndim>::kernrange;                ///< Maximum extent of kernel
	using SmoothingKernel<ndim>::invkernrange;             ///< 1/kernrange
	using SmoothingKernel<ndim>::kernrangesqd;             ///< kernrange^2

  GaussianKernel(string);
  ~GaussianKernel();

  // Gaussian kernel function prototypes
  //---------------------------------------------------------------------------
  FLOAT w0(const FLOAT);
  FLOAT w1(const FLOAT);
  FLOAT womega(const FLOAT);
  FLOAT wzeta(const FLOAT);
  FLOAT wgrav(const FLOAT);
  FLOAT wpot(const FLOAT);

};



//=================================================================================================
//  GaussianKernel<ndim>::w0
/// Main SPH smoothing kernel function, $W(s=r/h)$, for Gaussian kernel.
//=================================================================================================
template <int ndim>
inline FLOAT GaussianKernel<ndim>::w0(const FLOAT s)
{
  if (s < kernrange)
    return (kernnorm)*exp(-s*s);
  else
    return (FLOAT) 0.0;
}



//=================================================================================================
//  GaussianKernel<ndim>::w1
/// First spatial derivative of smoothing kernel, dWdr, for Gaussian kernel.
//=================================================================================================
template <int ndim>
inline FLOAT GaussianKernel<ndim>::w1(const FLOAT s)
{
  if (s < kernrange)
    return -(FLOAT) 2.0*(kernnorm)*s*exp(-s*s);
  else
    return (FLOAT) 0.0;
}



//=================================================================================================
//  GaussianKernel<ndim>::womega
/// Derivative of main SPH kernel w.r.t. smoothing length.
//=================================================================================================
template <int ndim>
inline FLOAT GaussianKernel<ndim>::womega(const FLOAT s)
{
  if (s < kernrange)
    return (kernnorm)*
      ((FLOAT) 2.0*s*exp(-s*s) - (ndim)*exp(-s*s));
  else
    return (FLOAT) 0.0;
}



//=================================================================================================
//  GaussianKernel<ndim>::wzeta
/// Derivative of gravitational potential kernel w.r.t. smoothing length.
//=================================================================================================
template <int ndim>
inline FLOAT GaussianKernel<ndim>::wzeta(const FLOAT s)
{
  if (s < kernrange)
    return 0.0;
  else
    return 0.0;
}



//=================================================================================================
//  GaussianKernel<ndim>::wgrav
/// Gravitational force kernel
//=================================================================================================
template <int ndim>
inline FLOAT GaussianKernel<ndim>::wgrav(const FLOAT s)
{
  if (s < kernrange)
    return 0.0;
  else
    return 0.0;
}



//=================================================================================================
//  GaussianKernel<ndim>::wpot
/// Gravitational potential kernel.
//=================================================================================================
template <int ndim>
inline FLOAT GaussianKernel<ndim>::wpot(const FLOAT s)
{
  if (s < kernrange)
    return 0.0;
  else
    return 0.0;
}



//=================================================================================================
//  Class TabulatedKernel
/// \brief   Class to tabulate arbitrary kernel function for quick look-up.
/// \details Tabulates all functions of any arbitrary kernel choice.
///          Allows quick look-up of kernel values rather than computing
///          relatively expensive polynomial (with/without divisions).
/// \author  G. Rosotti
/// \date    03/04/2013
//=================================================================================================
template <int ndim>
class TabulatedKernel: public SmoothingKernel<ndim>
{
	using SmoothingKernel<ndim>::kernnorm;                 ///< Kernel normalisation constant
	using SmoothingKernel<ndim>::kernnormdrag;             ///< Normalization factor for drag kernel
public:
	using SmoothingKernel<ndim>::kernrange;                ///< Maximum extent of kernel
	using SmoothingKernel<ndim>::invkernrange;             ///< 1/kernrange
	using SmoothingKernel<ndim>::kernrangesqd;               ///< kernrange^2
 private:

  SmoothingKernel<ndim>* kernel;       ///< Pointer to kernel object

  int res;                             ///< 'Resolution' of kernel table
  FLOAT resinvkernrange;               ///< ??
  FLOAT resinvkernrangesqd;            ///< ??
  FLOAT* tableW0;                      ///< Tabulated W (main kernel)
  FLOAT* tableW1;                      ///< Tabulated dW/dr kernel
  FLOAT* tableWomega;                  ///< Tabulated dW/dh kernel
  FLOAT* tableWzeta;                   ///< Tabulated zeta kernel
  FLOAT* tableWgrav;                   ///< Tabulated smoothed gravity kernel
  FLOAT* tableWpot;                    ///< Tabulated smoothed potential kernel
  FLOAT* tableWdrag;                   ///< Tabulated Drag Kernel (double hump)
  FLOAT* tableW0_s2;                   ///< Tabulated W with ssqd argument
  FLOAT* tableWomega_s2;               ///< Tabulated Womega with ssqd argument
  FLOAT* tableWzeta_s2;                ///< Tabulated Wzeta with ssqd argument
  FLOAT* tableLOS;                     ///< Tabulated Line-of-sight kernel

  void initializeTableLOS();



  //===============================================================================================
  //  TabulatedKernel<ndim>::initializeTable
  /// Initialise kernel table
  //===============================================================================================
  void initializeTable(FLOAT* table, FLOAT (SmoothingKernel<ndim>::*function) (const FLOAT s)) {
    const FLOAT step = kernel->kernrange/res;
    for (int i=0; i< res; i++) {
      table[i] = (kernel->*function)(step*i);
    }
  }

  //===============================================================================================
  //  TabulatedKernel<ndim>::duplicateTable
  /// Create a copy of a given table
  //===============================================================================================
  FLOAT* duplicateTable(FLOAT* table) const
  {
	  FLOAT * newtable = new FLOAT[res] ;
	  for (int i=0; i<res; i++){
		 newtable[i] = table[i] ;
	  }
	  return newtable ;
  }

  //===============================================================================================
  //  TabulatedKernel<ndim>::initializeTable
  /// Initialise kernel table
  //===============================================================================================
  void initializeTableSqd(FLOAT* table, FLOAT (SmoothingKernel<ndim>::*function) (const FLOAT s)) {
    const FLOAT step = kernel->kernrangesqd/res;
    for (int i=0; i< res; i++) {
      table[i] = (kernel->*function)(sqrt(step*i));
    }
  }


  //===============================================================================================
  //  TabulatedKernel<ndim>::tableLookup
  /// ..
  //===============================================================================================
  FLOAT tableLookup(FLOAT* table, const FLOAT s) {
    if (s >= (kernrange)) return (FLOAT) 0.0;
    FLOAT indexf = s*resinvkernrange;
    int index = (int) indexf;
    return table[index];
  }



  //===============================================================================================
  //  TabulatedKernel<ndim>::tableLookupSqd
  /// ..
  //===============================================================================================
  FLOAT tableLookupSqd(FLOAT* table, const FLOAT s) {
    if (s >= (kernrangesqd)) return (FLOAT) 0.0;
    FLOAT indexf = s*resinvkernrangesqd;
    int index = (int) indexf;
    return table[index];
  }



  //===============================================================================================
  //  TabulatedKernel<ndim>::GravTableLookup
  /// ..
  //===============================================================================================
  FLOAT GravTableLookup(FLOAT* table, const FLOAT s) {
    if (s >= (kernrange)) return (FLOAT) 1.0/(s*s);
    FLOAT indexf = s*resinvkernrange;
    int index = (int) indexf;
    return table[index];
  }



  //===============================================================================================
  //  TabulatedKernel<ndim>::GravPotTableLookup
  /// ..
  //===============================================================================================
  FLOAT GravPotTableLookup(FLOAT* table, const FLOAT s) {
    if (s >= (kernrange)) return (FLOAT) 1.0/s;
    FLOAT indexf = s*resinvkernrange;
    int index = (int) indexf;
    return table[index];
  }



 public:
  TabulatedKernel(string KernelName, int resaux=1000);
  TabulatedKernel(const TabulatedKernel<ndim>&);

  ~TabulatedKernel() {
    delete[] tableW0;
    delete[] tableW1;
    delete[] tableWomega;
    delete[] tableWzeta;
    delete[] tableWgrav;
    delete[] tableWpot;
    delete[] tableWdrag;
    delete[] tableW0_s2;
    delete[] tableWomega_s2;
    delete[] tableWzeta_s2;
    delete[] tableLOS;
  }

  FLOAT w0(const FLOAT s);
  FLOAT w0_s2(const FLOAT s);
  FLOAT w1(const FLOAT s);
  FLOAT womega(const FLOAT s);
  FLOAT womega_s2(const FLOAT s);
  FLOAT wzeta(const FLOAT s);
  FLOAT wzeta_s2(const FLOAT s);
  FLOAT wgrav(const FLOAT s);
  FLOAT wpot(const FLOAT s);
  FLOAT wdrag(const FLOAT s);
  FLOAT wLOS(const FLOAT s);

};



// Templated functions for tabulated kernel look-up
//-------------------------------------------------------------------------------------------------
template <int ndim>
inline FLOAT TabulatedKernel<ndim>::w0 (const FLOAT s) {
  return tableLookup(tableW0, s);
}

template <int ndim>
inline FLOAT TabulatedKernel<ndim>::w0_s2 (FLOAT s2) {
  return tableLookupSqd(tableW0_s2, s2);
}

template <int ndim>
inline FLOAT TabulatedKernel<ndim>::w1 (const FLOAT s) {
  return tableLookup(tableW1, s);
}

template <int ndim>
inline FLOAT TabulatedKernel<ndim>::womega (const FLOAT s) {
  return tableLookup(tableWomega, s);
}

template <int ndim>
inline FLOAT TabulatedKernel<ndim>::womega_s2 (FLOAT s2) {
  return tableLookupSqd(tableWomega_s2, s2);
}

template <int ndim>
inline FLOAT TabulatedKernel<ndim>::wzeta (const FLOAT s) {
  return tableLookup(tableWzeta, s);
}

template <int ndim>
inline FLOAT TabulatedKernel<ndim>::wzeta_s2 (FLOAT s2) {
  return tableLookupSqd(tableWzeta_s2, s2);
}

template <int ndim>
inline FLOAT TabulatedKernel<ndim>::wgrav (const FLOAT s) {
  return GravTableLookup(tableWgrav, s);
}

template <int ndim>
inline FLOAT TabulatedKernel<ndim>::wpot (const FLOAT s) {
  return GravPotTableLookup(tableWpot, s);
}

template <int ndim>
inline FLOAT TabulatedKernel<ndim>::wdrag (const FLOAT s) {
  return GravPotTableLookup(tableWdrag, s);
}

template <int ndim>
inline FLOAT TabulatedKernel<ndim>::wLOS (const FLOAT s) {
  return tableLookup(tableLOS, s);
}

#endif
