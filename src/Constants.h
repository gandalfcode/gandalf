// ============================================================================
// Constants.h
// Definitions for all astrophysical and dimensionless constants used 
// throughout the code.
// ============================================================================


#include "Precision.h"


#ifndef _CONSTANTS_H_
#define _CONSTANTS_H_


// Physical constants in SI units (unless stated otherwise)
static const DOUBLE r_pc    = 3.08568025E16;
static const DOUBLE r_au    = 1.49597870E11;
static const DOUBLE r_sun   = 6.955E8;
static const DOUBLE r_earth = 6.371E6;
static const DOUBLE pc_au   = 206264.9855910382;
static const DOUBLE km_cm   = 1.E5;
static const DOUBLE m_sun   = 1.98892E30;
static const DOUBLE m_jup   = 1.8986E27;
static const DOUBLE m_earth = 5.9736E24;
static const DOUBLE myr     = 3.1556952E13;
static const DOUBLE yr      = 3.1556952E7;
static const DOUBLE day     = 8.64E4;
static const DOUBLE amu     = 1.660538782E-27;
static const DOUBLE m_hydrogen  = 1.66054E-27;
static const DOUBLE G_const     = 6.67384E-11;
static const DOUBLE k_boltzmann = 1.3806503E-23;
static const DOUBLE stefboltz   = 5.67037321E-8;
static const DOUBLE e_charge    = 1.6021765E-19;
static const DOUBLE mu_0        = 1.25663706144E-6;
static const DOUBLE kappa_const = 2.09E-4;
static const DOUBLE L_sun       = 3.839E26;

// Dimensionless numerical constants
static const FLOAT pi = 3.1415926536;
static const FLOAT twopi = 6.283185307;
static const FLOAT invpi = 0.318309886;
static const FLOAT invlogetwo = 1.442695041;
static const FLOAT invlog10two = 3.321928095;
static const FLOAT invsqrttwo = 0.707106781;
static const FLOAT sqrttwo = 1.414213562;
static const FLOAT onethird = 0.333333333333;
static const FLOAT onesixth = 0.16666666666666;
static const FLOAT twothirds = 0.66666666666666;
static const FLOAT big_number = 9.9e20;
static const FLOAT small_number = 1.0e-20;

static const DOUBLE pi_dp = 3.1415926536;
static const DOUBLE twopi_dp = 6.283185307;
static const DOUBLE invpi_dp = 0.318309886;
static const DOUBLE invlogetwo_dp = 1.442695041;
static const DOUBLE invlog10two_dp = 3.321928095;
static const DOUBLE invsqrttwo_dp = 0.707106781;
static const DOUBLE sqrttwo_dp = 1.414213562;
static const DOUBLE onethird_dp = 0.333333333333333;
static const DOUBLE onesixth_dp = 0.166666666666666;
static const DOUBLE twothirds_dp = 0.666666666666666;
static const DOUBLE big_number_dp = 9.9e50;
static const DOUBLE small_number_dp = 1.0e-50;


#endif
