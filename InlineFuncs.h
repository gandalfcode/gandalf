// ============================================================================
// InlineFuncs.h
// Contains definitions of any useful small utility functions that can be 
// inlined to improve readability/performance of the code.
// ============================================================================


#ifndef _INLINE_FUNCS_H_
#define _INLINE_FUNCS_H_


#include <string>
#include "Dimensions.h"
using namespace std;

static const float kernnorm = invpi*7.0/478.0;


// ============================================================================
// DotProduct
// Calculates the dot product between two vectors, v1 and v2.
// For optimisation reasons, assumes both vectors are of length NDIM.
// ============================================================================
template <typename T>
static inline T DotProduct(T *v1, T *v2, int ndim)
{
  if (ndim == 1)
    return v1[0]*v2[0];
  else if (ndim == 2)
    return v1[0]*v2[0] + v1[1]*v2[1];
  else if (ndim == 3)
    return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
}



// ============================================================================
// ..
// ============================================================================
template <typename T>
static inline void PrintArray(string message, int Tsize, T *array)
{
  cout << message;
  for (int i=0; i<Tsize; i++) cout << array[i] << "  ";
  cout << endl;
  return;
}

inline float w0(float s)
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
inline float w1(float s)
{
  if (s < 1.0)
    return kernnorm*(-120.0*s + 120.0*pow(s,3) - 50.0*pow(s,4));
  else if (s < 2.0)
    return kernnorm*(75.0 - 420.0*s + 450.0*s*s -
             180.0*pow(s,3) + 25.0*pow(s,4));
  else if (s < 2.0)
    return kernnorm*(-405.0 + 540.0*s - 270.0*s*s +
             60.0*pow(s,3) - 5.0*pow(s,4));
  else
    return 0.0;
}



// ============================================================================
// QuinticKernel::womega
// ============================================================================
inline float womega(float s)
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

#endif
