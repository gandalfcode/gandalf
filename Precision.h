// ============================================================================
// Precision.h
// Contains macro definitions for floating point precision in all routines.
// ============================================================================

#ifndef _PRECISION_H_
#define _PRECISION_H_

#if defined(SINGLE_PRECISION)
//#define FLOAT float
typedef float FLOAT;
typedef double DOUBLE;
//#define DOUBLE double
#elif defined(DOUBLE_PRECISION)
typedef double FLOAT;
typedef double DOUBLE;
//#define FLOAT double
//#define DOUBLE double
#endif

#endif
