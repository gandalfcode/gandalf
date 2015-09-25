//=================================================================================================
//  InlineFuncs.h
//  Contains definitions of any useful small utility functions that can be
//  templated and/or inlined to improve readability/performance of the code.
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


#ifndef _INLINE_FUNCS_H_
#define _INLINE_FUNCS_H_


#include <assert.h>
#include <string>
#include <math.h>
#include "Precision.h"
#include "Constants.h"
#include "Particle.h"
#include "DomainBox.h"
using namespace std;


//=================================================================================================
//  DotProduct
//  Calculates the dot product between two vectors, v1 and v2,
//  of given length 'ndim'
//=================================================================================================
template <typename T>
static inline T DotProduct(T *v1, T *v2, int ndim)
{
  if (ndim == 1)
    return v1[0]*v2[0];
  else if (ndim == 2)
    return v1[0]*v2[0] + v1[1]*v2[1];
  else
    return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
}



//=================================================================================================
//  PrintArray
//  Print values of a given array to standard output
//=================================================================================================
template <typename T>
static inline void PrintArray(string message, int Tsize, T *array)
{
  cout << message;
  for (int i=0; i<Tsize; i++) cout << array[i] << "  ";
  cout << endl;
  return;
}



//=================================================================================================
//  min3
//  Return minimum of 3 given values.
//=================================================================================================
template <typename T>
static inline T min3(T v1, T v2, T v3)
{
   T vmin = v1;
   if (v2 < vmin) vmin = v2;
   if (v3 < vmin) vmin = v3;
   return vmin;
}



//=================================================================================================
//  max3
//  Return maximum of 3 given values.
//=================================================================================================
template <typename T>
static inline T max3(T v1, T v2, T v3)
{
   T vmax = v1;
   if (v2 > vmax) vmax = v2;
   if (v3 > vmax) vmax = v3;
   return vmax;
}



//=================================================================================================
//  sgn
//  Sign function.  Returns (a) -1 if T < 0, (b) 0 if T = 0, (c) +1 if T > 0.
//=================================================================================================
template <typename T>
static inline int sgn(T val)
{
  return (T(0) < val) - (val < T(0));
}



//=================================================================================================
//  Heapsort
//  Sorts a 1D array of values using the Heapsort algorithm.
//  (Courtesy of Anthony Whitworth - 18/04/2013)
//=================================================================================================
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



//=================================================================================================
//  InsertionSort
//  Sort list of quantities into ascending order.
//=================================================================================================
template <typename T>
static inline void InsertionSort
(int Nsort,                         ///< No. of values to be sorted
 T *r)                              ///< List of values to be sorted
{
  int i,j;
  T raux;

  for (j=1; j<Nsort; j++) {
    raux = r[j];
    for (i=j-1; i>=0; i--) {
      if (r[i] <= raux) break;
      r[i+1] = r[i];
    }
    r[i+1] = raux;
  }

  return;
}



//=================================================================================================
//  InsertionSortIds
//  Sort ...
//=================================================================================================
template <typename T>
static inline void InsertionSortIds
 (int Nsort,                           ///< No. of values to be sorted
  int *ids,                            ///< List of particle ids
  T *r)                                ///< ...
{
  int i,iaux,j;
  T raux;

  for (j=1; j<Nsort; j++) {
    raux = r[j];
    iaux = ids[j];
    for (i=j-1; i>=0; i--) {
      if (r[i] <= raux) break;
      r[i+1] = r[i];
      ids[i+1] = ids[i];
    }
    r[i+1] = raux;
    ids[i+1] = iaux;
  }

  return;
}



//=================================================================================================
//  EulerAngleRotation
//  Rotate given vector around specified Euler angles
//=================================================================================================
template <typename T>
static inline void EulerAngleMatrix
 (T phi,
  T theta,
  T psi,
  T Arot[3][3])
{
  Arot[0][0] = cos(theta)*cos(psi);
  Arot[1][0] = cos(phi)*sin(psi) + sin(phi)*sin(theta)*cos(psi);
  Arot[2][0] = sin(phi)*sin(psi) - cos(phi)*sin(theta)*cos(psi);
  Arot[0][1] = -cos(theta)*sin(psi);
  Arot[1][1] = cos(phi)*cos(psi) - sin(phi)*sin(theta)*sin(psi);
  Arot[2][1] = sin(phi)*cos(psi) + cos(phi)*sin(theta)*sin(psi);
  Arot[0][2] = sin(theta);
  Arot[1][2] = -sin(phi)*cos(theta);
  Arot[2][2] = cos(phi)*cos(theta);

  return;
}



//=================================================================================================
//  EulerAngleRotation
//  Rotate given vector around specified Euler angles
//=================================================================================================
template <typename T>
static inline void EulerAngleRotation
 (T phi,
  T theta,
  T psi,
  T vec[3])
{
  int k;
  T Arot[3][3];
  T vecaux[3];

  EulerAngleMatrix(phi,theta,psi,Arot);

  for (k=0; k<3; k++) vecaux[k] = vec[k];

  for (k=0; k<3; k++) {
    vec[k] = Arot[0][k]*vecaux[0] + Arot[1][k]*vecaux[1] + Arot[2][k]*vecaux[2];
  }

  cout << "rot angles : " << phi << "    " << theta << "    " << psi << endl;
  cout << "vec orig : " << vecaux[0] << "   " << vecaux[1] << "   " << vecaux[2] << endl;
  cout << "vec rot  : " << vec[0] << "   " << vec[1] << "   " << vec[2] << endl;

  return;
}



//=================================================================================================
//  EulerAngleArrayRotation
//  Rotate given array of vectors around specified Euler angles
//=================================================================================================
template <typename T>
static inline void EulerAngleArrayRotation
 (int N,
  T phi,
  T theta,
  T psi,
  T *vec)
{
  int i;
  int k;
  T Arot[3][3];
  T vecaux[3];

  EulerAngleMatrix(phi,theta,psi,Arot);

  for (i=0; i<N; i++) {
    for (k=0; k<3; k++) vecaux[k] = vec[3*i + k];

    for (k=0; k<3; k++) {
      vec[3*i + k] = Arot[0][k]*vecaux[0] + Arot[1][k]*vecaux[1] + Arot[2][k]*vecaux[2];
    }

  }

  return;
}



//=================================================================================================
//  clamp
//  ...
//=================================================================================================
static inline FLOAT clamp (FLOAT value, FLOAT min, FLOAT max)
{
  bool smaller = value < min;
  if (smaller) return min;
  bool bigger = value > max;
  if (bigger) return max;
  return value;
}




//=================================================================================================
//  BoxOverlap
/// Check if two bounding boxes overlap.  If yes, then returns true.
//=================================================================================================
static inline bool BoxOverlap
 (const int ndim,
  FLOAT *box1min,                      ///< Minimum extent of box 1
  FLOAT *box1max,                      ///< Maximum extent of box 1
  FLOAT *box2min,                      ///< Minimum extent of box 2
  FLOAT *box2max)                      ///< Maximum extent of box 2
{
  if (ndim == 1) {
    if (box1min[0] > box2max[0]) return false;
    if (box2min[0] > box1max[0]) return false;
    return true;
  }
  else if (ndim == 2) {
    if (box1min[0] > box2max[0]) return false;
    if (box2min[0] > box1max[0]) return false;
    if (box1min[1] > box2max[1]) return false;
    if (box2min[1] > box1max[1]) return false;
    return true;
  }
  else {
    if (box1min[0] > box2max[0]) return false;
    if (box2min[0] > box1max[0]) return false;
    if (box1min[1] > box2max[1]) return false;
    if (box2min[1] > box1max[1]) return false;
    if (box1min[2] > box2max[2]) return false;
    if (box2min[2] > box1max[2]) return false;
    return true;
  }
}



//=================================================================================================
//  BoxOverlap
/// Check if two bounding boxes overlap.  If yes, then returns true.
//=================================================================================================
static inline bool BoxOverlap
 (const int ndim,
  const FLOAT *box1min,                ///< Minimum extent of box 1
  const FLOAT *box1max,                ///< Maximum extent of box 1
  const FLOAT *box2min,                ///< Minimum extent of box 2
  const FLOAT *box2max)                ///< Maximum extent of box 2
{
  if (ndim == 1) {
    if (box1min[0] > box2max[0]) return false;
    if (box2min[0] > box1max[0]) return false;
    return true;
  }
  else if (ndim == 2) {
    if (box1min[0] > box2max[0]) return false;
    if (box2min[0] > box1max[0]) return false;
    if (box1min[1] > box2max[1]) return false;
    if (box2min[1] > box1max[1]) return false;
    return true;
  }
  else {
    if (box1min[0] > box2max[0]) return false;
    if (box2min[0] > box1max[0]) return false;
    if (box1min[1] > box2max[1]) return false;
    if (box2min[1] > box1max[1]) return false;
    if (box1min[2] > box2max[2]) return false;
    if (box2min[2] > box1max[2]) return false;
    return true;
  }
}



//=================================================================================================
//  FractionalBoxOverlap
/// Returns what fraction of box 1 overlaps box 2
//=================================================================================================
static inline FLOAT FractionalBoxOverlap
 (const int ndim,
  FLOAT *box1min,                      ///< Minimum extent of box 1
  FLOAT *box1max,                      ///< Maximum extent of box 1
  FLOAT *box2min,                      ///< Minimum extent of box 2
  FLOAT *box2max)                      ///< Maximum extent of box 2
{
  int k;                             // ..
  FLOAT area    = (FLOAT) 1.0;       // ..
  FLOAT overlap = (FLOAT) 1.0;       // ..
  FLOAT frac    = (FLOAT) 0.0;       // ..

  for (k=0; k<ndim; k++) assert(box1min[k] < box1max[k]);
  for (k=0; k<ndim; k++) assert(box2min[k] < box2max[k]);

  for (k=0; k<ndim; k++) {
    area *= (box1max[k] - box1min[k]);
    if (box1max[k] <= box2min[k] || box1min[k] >= box2max[k]) {
      overlap *= (FLOAT) 0.0;
    }
    else if (box1min[k] > box2min[k] && box1max[k] < box2max[k]) {
      overlap *= (box1max[k] - box1min[k]);
    }
    else if (box1min[k] < box2min[k] && box1max[k] > box2max[k]) {
      overlap *= (box2max[k] - box2min[k]);
    }
    else if (box1min[k] < box2min[k] && box1max[k] < box2max[k]) {
      overlap *= (box1max[k] - box2min[k]);
    }
    else if (box1min[k] > box2min[k] && box1max[k] > box2max[k]) {
      overlap *= (box2max[k] - box1min[k]);
    }
    else {
      cout << "Should never get here!! " << overlap << "    area : " << area << endl;
    }
  }

  if (area > 0.0) frac = overlap/area;
  assert(frac >= 0.0 && frac <= 1.0);

  return frac;
}



//=================================================================================================
//  FractionalBoxOverlap
/// Returns what fraction of box 1 overlaps box 2
//=================================================================================================
static inline FLOAT FractionalBoxOverlap
 (const int ndim,                      ///< Dimensionality
  const FLOAT *box2min,                ///< Minimum extent of box 1
  const FLOAT *box2max,                ///< Maximum extent of box 1
  const FLOAT *box1min,                ///< Minimum extent of box 2
  const FLOAT *box1max)                ///< Maximum extent of box 2
{
  int k;                               // Dimension index counter
  FLOAT area    = (FLOAT) 1.0;         // Area/volume of box 1
  FLOAT overlap = (FLOAT) 1.0;         // Overlap area/volume of boxes
  FLOAT frac    = (FLOAT) 0.0;         // Fraction of box1 that is overlapped by box 2

  for (k=0; k<ndim; k++) assert(box1min[k] < box1max[k]);
  for (k=0; k<ndim; k++) assert(box2min[k] < box2max[k]);

  for (k=0; k<ndim; k++) {
    area *= (box1max[k] - box1min[k]);
    if (box1max[k] <= box2min[k] || box1min[k] >= box2max[k]) {
      overlap *= (FLOAT) 0.0;
    }
    else if (box1min[k] > box2min[k] && box1max[k] < box2max[k]) {
      overlap *= (box1max[k] - box1min[k]);
    }
    else if (box1min[k] < box2min[k] && box1max[k] > box2max[k]) {
      overlap *= (box2max[k] - box2min[k]);
    }
    else if (box1min[k] < box2min[k] && box1max[k] < box2max[k]) {
      overlap *= (box1max[k] - box2min[k]);
    }
    else if (box1min[k] > box2min[k] && box1max[k] > box2max[k]) {
      overlap *= (box2max[k] - box1min[k]);
    }
    else {
      cout << "Should never get here!! " << overlap << "    area : " << area << endl;
    }
  }

  if (area > 0.0) frac = overlap/area;
  assert(frac >= 0.0 && frac <= 1.0);

  return frac;
}



//=================================================================================================
//  ParticleBoxOverlap
//  ...
//=================================================================================================
template <int ndim>
static inline bool ParticleBoxOverlap (SphParticle<ndim>& part, Box<ndim>& box)
{
  // Find the closest point to the circle within the rectangle
  FLOAT closest_coord[ndim];
  for (int i=0; i<ndim; i++) {
    closest_coord[i] = clamp(part.r[i],box.boxmin[i],box.boxmax[i]);
  }

  // Calculate the distance between the circle's center and this closest point
  FLOAT distanceSquared = 0.;
  for (int i=0; i<ndim; i++) {
    FLOAT distance_coord = closest_coord[i] - part.r[i];
    distanceSquared += distance_coord*distance_coord;
  }

  // If the distance is less than the circle's radius, an intersection occurs
  return distanceSquared < part.hrangesqd;
}



//=================================================================================================
//  ParticleInBox
//  ..
//=================================================================================================
template <int ndim, template<int> class Particle >
static inline bool ParticleInBox (Particle<ndim>& part, Box<ndim>& box)
{
  for (int k=0; k<ndim; k++) {
    if (part.r[k] < box.boxmin[k] || part.r[k] > box.boxmax[k]) return false;
  }
  return true;
}



//=================================================================================================
//  compute_displs
/// Given a vector of counts, compute the displacement of the elements and save it
/// in displs
//=================================================================================================
static inline void compute_displs (std::vector<int>& displs, std::vector<int>& counts) {

  const int size = displs.size();

  if (size==0)
    return;

  displs[0]=0;

  if (size==1)
    return;

  for (int i=1; i<size; i++) {
    displs[i] = counts[i-1]+displs[i-1];
  }

}



//=================================================================================================
//  VerifyUniqueIds
/// ...
//=================================================================================================
static inline bool VerifyUniqueIds(const int N, const int Nrange, const int *values)
{
  int i,j;
  int *counter = new int[Nrange];
  for (i=0; i<Nrange; i++) counter[i] = 0;

  for (int j=0; j<N; j++) {
    i = values[j];
    if (!(i >= 0 && i < Nrange)) {
      cout << "i : " << i << "     Nrange : " << Nrange << "    N : " << N << endl;
      delete[] counter;
      return false;
    }
    counter[i]++;
  }

  for (j=0; j<N; j++) {
    if (!(counter[i] == 0 || counter[i] == 1)) {
      cout << "Invalid counter : " << i << "    count : " << counter[i] << endl;
      delete[] counter;
      return false;
    }
  }

  delete[] counter;

  return true;
}


//=================================================================================================
//  ComputeTimestepLevel
/// Compute the timestep level of dt by a block hierarchy defined by a maximum timestep dt_max.
/// Also ensure the timestep level is never negative (since all timesteps greater than the
/// maximum timestep are automatically anchored to the bottom timestep level).
//=================================================================================================
static inline int ComputeTimestepLevel
 (const DOUBLE dt,                     ///< [in] Timestep to compute level for
  const DOUBLE dt_max)                 ///< [in] Maximum timestep, corresponding to minimum level
{
  assert(dt_max > 0.0);
  assert(dt > 0.0);
  int level = max((int) (invlogetwo*log(dt_max/dt)) + 1, 0);
  return (int) level;
}

#endif
