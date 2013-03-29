// ============================================================================
// InlineFuncs.h
// Contains definitions of any useful small utility functions that can be 
// inlined to improve readability/performance of the code.
// ============================================================================


#ifndef _INLINE_FUNCS_H_
#define _INLINE_FUNCS_H_


#include <string>
#include <math.h>
#include "Precision.h"
#include "Constants.h"
using namespace std;


// ============================================================================
// DotProduct
// Calculates the dot product between two vectors, v1 and v2, 
// of given length 'ndim'
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
// PrintArray
// Print values of a given array to standard output
// ============================================================================
template <typename T>
static inline void PrintArray(string message, int Tsize, T *array)
{
  cout << message;
  for (int i=0; i<Tsize; i++) cout << array[i] << "  ";
  cout << endl;
  return;
}


template <typename T>
static inline T min3(T v1, T v2, T v3)
{
   T vmin = v1;
   if (v2 < vmin) vmin = v2;
   if (v3 < vmin) vmin = v2;
   return vmin;
}


template <typename T>
static inline T max3(T v1, T v2, T v3)
{
   T vmax = v1;
   if (v2 > vmax) vmax = v2;
   if (v3 > vmax) vmax = v2;
   return vmax;
}


template <typename T> int sgn(T val) {
	return (T(0) < val) - (val < T(0));
}


//file wave.c
/*
+ wave - Nonlinear wave speeds.
     Description :
         This function calculates the wave speed for a wave connecting
         states with pressures pi, p ahead an behind respectively.
*/
static inline FLOAT wave(FLOAT p, FLOAT pi, FLOAT g3, FLOAT g4)
{
  FLOAT x, w;
  x = p/pi;
  if (fabs(x - 1.0) < 1.0e-03)
    /* Use linear expression */
    w = 1.0 + 0.5*g3*(x - 1.0);
  else{
    /* Use non-linear expression */
    if (x >= 1.0)
      /* Shock */
      w = sqrt(1.0 + g3*(x - 1.0));
    else
      /* Rarefaction */
      w = g4*(1.0 - x)/(1.0 - (FLOAT)pow((DOUBLE)x, (DOUBLE)g4));
  }
  return(w);
}
//end wave.c




#endif
