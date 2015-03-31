//=================================================================================================
//  Ewald.h
//  Contains class definition for computing periodic gravity in 1D, 2D and 3D.
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


#ifndef _EWALD_H_
#define _EWALD_H_


#include <iostream>
#include <string>
#include <cmath>
#include <math.h>
#include "Precision.h"
#include "DomainBox.h"
#include "Constants.h"
#include "CodeTiming.h"
using namespace std;
#ifdef GANDALF_GSL
#include <gsl/gsl_sf.h>
#endif



//=================================================================================================
//  Class Ewald
/// \brief   ...
/// \details ...
/// \author  F. Dinnbier & D. A. Hubber
/// \date    09/10/2014
//=================================================================================================
template <int ndim>
class Ewald
{
 public:

  // Constructor and destructor
  //-----------------------------------------------------------------------------------------------
  Ewald(DomainBox<ndim> &, int, int, int, DOUBLE, DOUBLE, DOUBLE, DOUBLE, CodeTiming *);
  ~Ewald();


  // Other functions
  //-----------------------------------------------------------------------------------------------
  void CalculatePeriodicCorrection(FLOAT, FLOAT *, FLOAT *, FLOAT &);
  DOUBLE erfcx(DOUBLE);
  DOUBLE SimpsonInt (DOUBLE (Ewald<ndim>::*f)(DOUBLE, int, DOUBLE, DOUBLE),
                     int, DOUBLE, DOUBLE, DOUBLE, DOUBLE, int);
  DOUBLE GravInt2p1i (int, int, DOUBLE, DOUBLE, DOUBLE, DOUBLE, DOUBLE, DOUBLE);
  DOUBLE DerGravInt2p1i(int, int, DOUBLE, DOUBLE, DOUBLE, DOUBLE, DOUBLE, DOUBLE);
  DOUBLE IntFPot(DOUBLE, int, DOUBLE, DOUBLE);
  DOUBLE IntFAcc(DOUBLE, int, DOUBLE, DOUBLE);
  DOUBLE AccShort(DOUBLE, DOUBLE, DOUBLE);
  DOUBLE PotLong1p2i(int, DOUBLE, DOUBLE, DOUBLE, DOUBLE, DOUBLE);
  DOUBLE AccLong1p2iPer(int, DOUBLE, DOUBLE, DOUBLE, DOUBLE, DOUBLE);
  DOUBLE AccLong1p2iIso(int, DOUBLE, DOUBLE, DOUBLE, DOUBLE, DOUBLE);
  DOUBLE PotLong2p1i(DOUBLE, DOUBLE, DOUBLE, DOUBLE,
                     DOUBLE, DOUBLE, int, int, DOUBLE, DOUBLE);
  DOUBLE AccLong2p1iPer(DOUBLE, DOUBLE, DOUBLE, DOUBLE, DOUBLE,
                         DOUBLE, int, int, DOUBLE, DOUBLE);
  DOUBLE AccLong2p1iIso(DOUBLE, DOUBLE, DOUBLE, DOUBLE, DOUBLE,
                         DOUBLE, int, int, DOUBLE, DOUBLE);
  DOUBLE AccLong3pPer(int, DOUBLE, DOUBLE, DOUBLE, DOUBLE);


  // Ewald class variables
  //-----------------------------------------------------------------------------------------------
  /*static const int ewald_periodicity = 7;
  static const int gr_bhewaldseriesn = 10;
  static const DOUBLE lx_per = 2*5.6415e+18;
  static const DOUBLE ly_per = 3*5.6415e+18;
  static const DOUBLE lz_per = 4*5.6415e+18;
  static const double ewald_mult = 1.0;
  static const DOUBLE ixmin=1.0e-8;
  static const DOUBLE ixmax=10.0;
  static const int in=3000;*/
  int ewald_periodicity;
  int one_component;
  int Ngrid[3]; 				// Number of gridpoints in x, y, z
  DOUBLE* ewald_field;
  DOUBLE* ewald_fields;
  DOUBLE* ewald_fieldl;
//  DOUBLE* ewald_coord;
  DOUBLE dI[3];
  DOUBLE accPlane;
  DOUBLE potC1p2i;
  const int gr_bhewaldseriesn;
  const int in;
  const DOUBLE ewald_mult;
  const DOUBLE ixmin;
  const DOUBLE ixmax;
  const DOUBLE EFratio;
  const DOUBLE lx_per;
  const DOUBLE ly_per;
  const DOUBLE lz_per;
  const int nEwaldGrid;

  CodeTiming* timing;

};
#endif
