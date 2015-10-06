//=================================================================================================
//  Constants.h
//  Definitions for all astrophysical and dimensionless constants used
//  throughout the code.
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


#ifndef _CONSTANTS_H_
#define _CONSTANTS_H_


#include "Precision.h"


// Physical constants in SI units (unless stated otherwise)
//-------------------------------------------------------------------------------------------------
static const DOUBLE r_pc            = 3.08568025E16;
static const DOUBLE r_au            = 1.49597870E11;
static const DOUBLE r_sun           = 6.955E8;
static const DOUBLE r_earth         = 6.371E6;
static const DOUBLE pc_au           = 206264.9855910382;
static const DOUBLE km_cm           = 1.E5;
static const DOUBLE m_sun           = 1.98892E30;
static const DOUBLE m_jup           = 1.8986E27;
static const DOUBLE m_earth         = 5.9736E24;
static const DOUBLE myr             = 3.1556952E13;
static const DOUBLE yr              = 3.1556952E7;
static const DOUBLE day             = 8.64E4;
static const DOUBLE amu             = 1.660538782E-27;
static const DOUBLE m_hydrogen      = 1.66054E-27;
static const DOUBLE invm_hydrogen   = 6.02217E+26;
static const DOUBLE c_light         = 2.99792458E8;
static const DOUBLE G_const         = 6.67384E-11;
static const DOUBLE k_boltzmann     = 1.3806503E-23;
static const DOUBLE stefboltz       = 5.67037321E-8;
static const DOUBLE e_charge        = 1.6021765E-19;
static const DOUBLE mu_0            = 1.25663706144E-6;
static const DOUBLE kappa_const     = 2.09E-4;
static const DOUBLE L_sun           = 3.839E26;


// Dimensionless numerical constants
//-------------------------------------------------------------------------------------------------
static const FLOAT pi               = (FLOAT) 3.14159265358979;
static const FLOAT twopi            = (FLOAT) 6.28318530717959;
static const FLOAT invpi            = (FLOAT) 0.31830988618379;
static const FLOAT invlogetwo       = (FLOAT) 1.44269504088896;
static const FLOAT invlog10two      = (FLOAT) 3.321928095;
static const FLOAT invsqrttwo       = (FLOAT) 0.707106781;
static const FLOAT sqrttwo          = (FLOAT) 1.414213562;
static const FLOAT onethird         = (FLOAT) 0.33333333333333333333333;
static const FLOAT onesixth         = (FLOAT) 0.16666666666666666666666;
static const FLOAT twothirds        = (FLOAT) 0.66666666666666666666666;
static const FLOAT onetwelfth       = (FLOAT) 0.08333333333333333333333;
static const FLOAT big_number       = (FLOAT) 9.9e20;
static const FLOAT small_number     = (FLOAT) 1.0e-20;


// Dimensionless numerical constants in double precision
//-------------------------------------------------------------------------------------------------
static const DOUBLE pi_dp           = 3.14159265358979;
static const DOUBLE twopi_dp        = 6.28318530717959;
static const DOUBLE invpi_dp        = 0.31830988618379;
static const DOUBLE invlogetwo_dp   = 1.44269504088896;
static const DOUBLE invlog10two_dp  = 3.321928095;
static const DOUBLE invsqrttwo_dp   = 0.707106781;
static const DOUBLE sqrttwo_dp      = 1.414213562;
static const DOUBLE onethird_dp     = 0.3333333333333333333333;
static const DOUBLE onesixth_dp     = 0.1666666666666666666666;
static const DOUBLE twothirds_dp    = 0.6666666666666666666666;
static const DOUBLE onetwelfth_dp   = 0.0833333333333333333333;
static const DOUBLE big_number_dp   = 9.9e50;
static const DOUBLE small_number_dp = 1.0e-50;


#endif
