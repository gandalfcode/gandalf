//=============================================================================
//  Precision.h
//  Contains macro definitions for floating point precision in all routines.
//=============================================================================


#ifndef _PRECISION_H_
#define _PRECISION_H_

#if defined(SINGLE_PRECISION)
typedef float FLOAT;
typedef double DOUBLE;
#elif defined(DOUBLE_PRECISION)
typedef double FLOAT;
typedef double DOUBLE;
#endif

#endif
