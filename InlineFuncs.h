// ============================================================================
// InlineFuncs.h
// ============================================================================


#ifndef _INLINE_FUNCS_H_
#define _INLINE_FUNCS_H_


#include "Dimensions.h"


// ============================================================================
// DotProduct
// Calculates the dot product between two vectors, v1 and v2.
// For optimisation reasons, assuming both vectors are of length NDIM.
// ============================================================================
template <typename T>
static inline T DotProduct(T *v1, T *v2)
{
#if NDIM == 1
  return v1[0]*v2[0];
#elif NDIM == 2
  return v1[0]*v2[0] + v1[1]*v2[1];
#elif NDIM == 3
  return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
#endif
}


#endif
