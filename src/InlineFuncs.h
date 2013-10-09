//=============================================================================
//  InlineFuncs.h
//  Contains definitions of any useful small utility functions that can be 
//  inlined to improve readability/performance of the code.
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
//=============================================================================


#ifndef _INLINE_FUNCS_H_
#define _INLINE_FUNCS_H_


#include <string>
#include <math.h>
#include "Precision.h"
#include "Constants.h"
using namespace std;


//=============================================================================
//  DotProduct
//  Calculates the dot product between two vectors, v1 and v2, 
//  of given length 'ndim'
//=============================================================================
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



//=============================================================================
//  PrintArray
//  Print values of a given array to standard output
//=============================================================================
template <typename T>
static inline void PrintArray(string message, int Tsize, T *array)
{
  cout << message;
  for (int i=0; i<Tsize; i++) cout << array[i] << "  ";
  cout << endl;
  return;
}



//=============================================================================
//  min3
//  Return minimum of 3 given values.
//=============================================================================
template <typename T>
static inline T min3(T v1, T v2, T v3)
{
   T vmin = v1;
   if (v2 < vmin) vmin = v2;
   if (v3 < vmin) vmin = v2;
   return vmin;
}



//=============================================================================
//  max3
//  Return maximum of 3 given values.
//=============================================================================
template <typename T>
static inline T max3(T v1, T v2, T v3)
{
   T vmax = v1;
   if (v2 > vmax) vmax = v2;
   if (v3 > vmax) vmax = v2;
   return vmax;
}



//=============================================================================
//  sgn
//  Sign function.  Returns (a) -1 if T < 0, (b) 0 if T = 0, (c) +1 if T > 0.
//=============================================================================
template <typename T> 
static inline int sgn(T val)
{
  return (T(0) < val) - (val < T(0));
}



//=============================================================================
//  Heapsort
//  Sorts a 1D array of values using the Heapsort algorithm.
//  (Courtesy of Anthony Whitworth - 18/04/2013)
//=============================================================================
template <typename T>
static inline void Heapsort
(int q_TOT,                         ///< No. of values to be sorted
 int *qV,                           ///< Sorted ids of q-values
 T *V)                              ///< (Templated) array of values to sort
{
  int q,qq;                         // Dummy ids
  int q1;                           // Dummy id (for comparison)
  int q2;                           // Buffer for ranked id

  // Give V-value arbitrary rank
  for (q=0; q<q_TOT; q++) qV[q] = q_TOT - q - 1;

  // Build heap
  for (q=1; q<q_TOT; q++) {
    qq = q;  q1 = qq/2;

    while(V[qV[qq]] > V[qV[q1]]) {
      q2 = qV[qq];  qV[qq] = qV[q1];   qV[q1] = q2;
      if (q1 == 0) break;
      qq = q1;   q1 = qq/2;
    }
  }

  // Check local ordering
#if defined(VERIFY_ALL)
  for (q=1; q<q_TOT; q++)
    if (V[qV[q]] > V[qV[q/2]]) cout << "Tree not locally hierarchical" << endl;
#endif

  // Invert heap
  for (q=q_TOT-1; q>0; q--) {
    q2 = qV[q];  qV[q] = qV[0];   qV[0] = q2;
    if (q == 1) break;
    qq = 0; q1 = 1;
    if (V[qV[1]] < V[qV[2]] && q > 2) q1 = 2;
    while (V[qV[qq]] < V[qV[q1]]) {
      q2 = qV[qq];  qV[qq] = qV[q1];   qV[q1] = q2;
      qq = q1;   q1 = 2*qq;
      if (q1 >= q) break;
      if (V[qV[q1]] < V[qV[q1+1]] && q1+1 < q) q1++;
    }
  }

#if defined(VERIFY_ALL)
  for (q=1; q<q_TOT; q++)
    if (V[qV[q]] < V[qV[q-1]]) 
      cout << "Not properly ranked : "
	   << q << "   " 
	   << q_TOT << "   "
	   << V[qV[q-2]] << "   "
	   << V[qV[q-1]] << "   "
	   << V[qV[q]] << "   " 
	   << V[qV[q+1]] << endl;
#endif

  return;
}



//=============================================================================
//  InsertionSortIds
//  ..
//=============================================================================
static inline void InsertionSortIds
(int Nsort,                         ///< No. of values to be sorted
 int *ids)                          ///< List of particle ids
{
  int i,iaux,j;

  for (j=1; j<Nsort; j++) {
    iaux = ids[j];
    for (i=j-1; i>=0; i--) {
      if (ids[i] <= iaux) break;
      ids[i+1] = ids[i];
    }
    ids[i+1] = iaux;
  }

  return;
}
#endif
