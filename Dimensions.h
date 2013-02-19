// ============================================================================
// Dimensions.h
// Definitions of all static dimensionality variables.  
// If FIXED_DIMENSIONS macro is defined, then all dimensionality variables 
// are set statically here.  Otherwise, some variables are set locally in 
// each class definition.
// ============================================================================


#include "Precision.h"
#include "Constants.h"


#ifndef _DIMENSIONS_H_
#define _DIMENSIONS_H_


// Maximum allowed no. of dimensions chosen in Makefile
// ----------------------------------------------------------------------------
#if NDIM==1
static const int ndimmax = 1;
static const int bdimmax = 1;
#elif NDIM==2
static const int ndimmax = 2;
static const int bdimmax = 2;
#elif NDIM==3
static const int ndimmax = 3;
static const int bdimmax = 3;
#endif
#if VDIM==1
static const int vdimmax = 1;
#elif VDIM==2
static const int vdimmax = 2;
#elif VDIM==3
static const int vdimmax = 3;
#endif


// Dimensions are fixed at compilation-time with parameters
// ----------------------------------------------------------------------------
#if defined(FIXED_DIMENSIONS)
#if NDIM==1
static const int ndim = 1;
static const int bdim = 1;
static const FLOAT ndimpr = 1.0;
static const FLOAT invndim = 1.0;
static const FLOAT ndimplus1 = 2.0;
#elif NDIM==2
static const int ndim = 2;
static const int bdim = 2;
static const FLOAT ndimpr = 2.0;
static const FLOAT invndim = 0.5;
static const FLOAT ndimplus1 = 3.0;
#elif NDIM==3
static const int ndim = 3;
static const int bdim = 3;
static const FLOAT ndimpr = 3.0;
static const FLOAT invndim = onethird;
static const FLOAT ndimplus1 = 4.0;
#endif
#if VDIM==1
static const int vdim = 1;
#elif VDIM==2
static const int vdim = 2;
#elif VDIM==3
static const int vdim = 3;
#endif
static const int ndimp1 = ndim + 1;
static const int ndimm1 = ndim - 1;
static const int vdimp1 = vdim + 1;
static const int vdimm1 = vdim - 1;
#endif
// ----------------------------------------------------------------------------


#endif
