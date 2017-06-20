//=================================================================================================
//  Multipole.h
//  Header file containing class definitions multipole gravitational forces
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

#ifndef _MULTIPOLE_H_
#define _MULTIPOLE_H_

#include <math.h>
#include "TreeCell.h"
using namespace std;

enum multipole_method { monopole, quadrupole, fast_monopole, fast_quadrupole } ;


//=================================================================================================
// Gravity force functions
///
/// These functions compute the force on particles due to distant tree cells using one of a
/// variety of approximations:
///   Monopole or Quadropole forces from each cell summed per particle
///   Fast Monopole, per cell, Taylor expanded to the location of each particle.
//=================================================================================================

//=================================================================================================
//  ComputeCellMonopoleForces
/// Compute the force on particle 'parti' due to all cells obtained in the
/// gravity tree walk.  Uses only monopole moments (i.e. COM) of the cell.
//=================================================================================================
template<int ndim>
void ComputeCellMonopoleForces(FLOAT &gpot, ///< [inout] Grav. potential
FLOAT agrav[ndim],                   ///< [inout] Acceleration array
FLOAT rp[ndim],                      ///< [in] Position of point
int Ngravcell,                       ///< [in] No. of tree cells in list
MultipoleMoment<ndim> *gravcell)     ///< [in] List of tree cell ids
{
  FLOAT dr[ndim];                      // Relative position vector

  // Loop over all neighbouring particles in list
  //-----------------------------------------------------------------------------------------------
  for (int cc=0; cc<Ngravcell; cc++) {
    MultipoleMoment<ndim>& cell = gravcell[cc];

    FLOAT mc = cell.m;
    for (int k=0; k<ndim; k++) dr[k] = cell.r[k] - rp[k];
    FLOAT drsqd    = DotProduct(dr,dr,ndim) + small_number;
    FLOAT invdrsqd = (FLOAT) 1.0/drsqd;
    FLOAT invdrmag = sqrt(invdrsqd);
    FLOAT invdr3   = invdrsqd*invdrmag;

    gpot += mc*invdrmag;
    for (int k=0; k<ndim; k++) agrav[k] += mc*dr[k]*invdr3;

  }
  //-----------------------------------------------------------------------------------------------

  return;
}


//=================================================================================================
//  ComputeQuadrupole
/// Compute the quadropole part of the force, depending on dimension.
//=================================================================================================
inline void ComputeQuadropole(const MultipoleMoment<1>& cell, const FLOAT rp[1],
                              FLOAT agrav[1], FLOAT& gpot) {

  FLOAT dr[1] ;
  for (int k=0; k<1; k++) dr[k] = rp[k] - cell.r[k];
  FLOAT drsqd    = DotProduct(dr,dr,1) + small_number;
  FLOAT invdrsqd = (FLOAT) 1.0/drsqd;
  FLOAT invdrmag = sqrt(invdrsqd);
  FLOAT invdr5 = invdrsqd*invdrsqd*invdrmag ;

  // First add monopole term for acceleration
  for (int k=0; k<1; k++) agrav[k] -= cell.m*dr[k]*invdrsqd*invdrmag;

  FLOAT qscalar = cell.q[0]*dr[0]*dr[0];
  FLOAT qfactor = 2.5*qscalar*invdr5*invdrsqd;

  agrav[0] += (cell.q[0]*dr[0])*invdr5 - qfactor*dr[0];
  gpot += cell.m*invdrmag + 0.5*qscalar*invdr5;
}
inline void ComputeQuadropole(const MultipoleMoment<2>& cell, const FLOAT rp[2],
                              FLOAT agrav[2], FLOAT& gpot) {
  FLOAT dr[2] ;
  for (int k=0; k<2; k++) dr[k] = rp[k] - cell.r[k];
  FLOAT drsqd    = DotProduct(dr,dr,2) + small_number;
  FLOAT invdrsqd = (FLOAT) 1.0/drsqd;
  FLOAT invdrmag = sqrt(invdrsqd);
  FLOAT invdr5 = invdrsqd*invdrsqd*invdrmag ;

  // First add monopole term for acceleration
  for (int k=0; k<2; k++) agrav[k] -= cell.m*dr[k]*invdrsqd*invdrmag;

  FLOAT qscalar = cell.q[0]*dr[0]*dr[0] + cell.q[2]*dr[1]*dr[1] +
    2.0*cell.q[1]*dr[0]*dr[1];
  FLOAT qfactor = 2.5*qscalar*invdr5*invdrsqd;

  agrav[0] += (cell.q[0]*dr[0] + cell.q[1]*dr[1])*invdr5 - qfactor*dr[0];
  agrav[1] += (cell.q[1]*dr[0] + cell.q[2]*dr[1])*invdr5 - qfactor*dr[1];
  gpot += cell.m*invdrmag + 0.5*qscalar*invdr5;
}
inline void ComputeQuadropole(const MultipoleMoment<3>& cell, const FLOAT rp[3],
                              FLOAT agrav[3], FLOAT& gpot) {
  FLOAT dr[3] ;
  for (int k=0; k<3; k++) dr[k] = rp[k] - cell.r[k]; //cell.r[k] - rp[k];
  FLOAT drsqd    = DotProduct(dr,dr,3) + small_number;
  FLOAT invdrsqd = (FLOAT) 1.0/drsqd;
  FLOAT invdrmag = sqrt(invdrsqd);
  FLOAT invdr5 = invdrsqd*invdrsqd*invdrmag ;

  // First add monopole term for acceleration
  for (int k=0; k<3; k++) agrav[k] -= cell.m*dr[k]*invdrsqd*invdrmag;

  FLOAT qscalar = cell.q[0]*dr[0]*dr[0] + cell.q[2]*dr[1]*dr[1] -
          (cell.q[0] + cell.q[2])*dr[2]*dr[2] +
           2.0*(cell.q[1]*dr[0]*dr[1] + cell.q[3]*dr[0]*dr[2] + cell.q[4]*dr[1]*dr[2]);

  FLOAT qfactor = 2.5*qscalar*invdr5*invdrsqd;

  agrav[0] +=
      (cell.q[0]*dr[0] + cell.q[1]*dr[1] + cell.q[3]*dr[2])*invdr5 - qfactor*dr[0];
  agrav[1] +=
      (cell.q[1]*dr[0] + cell.q[2]*dr[1] + cell.q[4]*dr[2])*invdr5 - qfactor*dr[1];
  agrav[2] +=
      (cell.q[3]*dr[0] + cell.q[4]*dr[1] - (cell.q[0] + cell.q[2])*dr[2])*invdr5 - qfactor*dr[2];

  gpot += cell.m*invdrmag + 0.5*qscalar*invdr5;
}


//=================================================================================================
//  ComputeCellQuadrupoleForces
/// Compute the force on particle 'parti' due to all cells obtained in the
/// gravity tree walk including the quadrupole moment correction term.
//=================================================================================================
template <int ndim>
void ComputeCellQuadrupoleForces
 (FLOAT &gpot,                         ///< [inout] Grav. potential
  FLOAT agrav[ndim],                   ///< [inout] Acceleration array
  FLOAT rp[ndim],                      ///< [in] Position of point
  int Ngravcell,                       ///< [in] No. of tree cells in list
  MultipoleMoment<ndim> *gravcell)     ///< [in] List of tree cell ids
{

  // Loop over all neighbouring particles in list
  //-----------------------------------------------------------------------------------------------
  for (int cc=0; cc<Ngravcell; cc++) {
    ComputeQuadropole(gravcell[cc], rp, agrav, gpot) ;
  }
  //-----------------------------------------------------------------------------------------------


  return;
}

//=================================================================================================
//  class FastMultipoleForces
/// \brief Class for computing the gravitational forces using the fast monopole expansion
//=================================================================================================
template<int ndim>
class FastMultipoleForces {
public:

  FastMultipoleForces() {} ;
  FastMultipoleForces(const FLOAT r[ndim])
  {
    set_target_cell(r) ;
  }

  void set_target_cell(const FLOAT r[ndim]) {
    for (int k=0; k<ndim; k++) rc[k]   = r[k];
    reset() ;
  }

  inline void AddMonopoleContribution(const MultipoleMoment<ndim>& cell) ;
  inline void AddQuadrupoleContribution(const MultipoleMoment<ndim>& cell) ;

  inline void ApplyForcesTaylor(const FLOAT r[ndim], FLOAT agrav[ndim], FLOAT& gpot) const ;

private:
  void reset() {
    for (int k=0; k<ndim; k++) ac[k]   = 0;
    for (int k=0; k<ndim; k++) dphi[k] = 0;
    for (int k=0; k<(ndim*(ndim+1))/2; k++) q[k] = 0;
    for (int k=0; k<(ndim*(ndim+1)*(ndim+2))/6; k++) q2[k] = 0;

    pot = 0 ;
  }


  static const bool HoT = false ; // Include 2nd order terms in taylor expansion

  FLOAT rc[ndim] ;
  FLOAT ac[ndim]  ;
  FLOAT dphi[ndim] ;
  FLOAT q[(ndim*(ndim+1))/2];
  FLOAT q2[(ndim*(ndim+1)*(ndim+2))/6];
  FLOAT pot ;
};

//=================================================================================================
//  AddCellContribution
/// Compute the fast monopole cell-cell terms depending on dimension.
//=================================================================================================
template<>
inline void FastMultipoleForces<1>::AddMonopoleContribution(const MultipoleMoment<1>& cell) {
  FLOAT dr[1] ;
  for (int k=0; k<1; k++) dr[k] = cell.r[k] - rc[k];
  FLOAT invdrmag = sqrt((FLOAT) 1.0/DotProduct(dr,dr,1));
  FLOAT invdrsqd = invdrmag*invdrmag;
  FLOAT invdr3   = invdrsqd*invdrmag;

  FLOAT mc = cell.m;
  pot += mc*invdrmag;

  mc *= invdr3;
  for (int k=0; k<1; k++) ac[k] += mc*dr[k];
  for (int k=0; k<1; k++) dphi[k] += mc*dr[k];
  q[0] += mc*(3.0*dr[0]*dr[0]*invdrsqd - 1);

  if (HoT) {
    //                     (5  x_i   x_j   x_k  / r^2    - (x_i d_j,k + x_j d_i,k + x_k d_i,j))
    q2[0] += 3*mc*invdrsqd*(5*dr[0]*dr[0]*dr[0]*invdrsqd - (  dr[0]   +   dr[0]   +   dr[0]  ));
  }
}
template<>
inline void FastMultipoleForces<2>::AddMonopoleContribution(const MultipoleMoment<2>& cell) {
  FLOAT dr[2] ;
  for (int k=0; k<2; k++) dr[k] = cell.r[k] - rc[k];
  FLOAT invdrmag = sqrt((FLOAT) 1.0/DotProduct(dr,dr,2));
  FLOAT invdrsqd = invdrmag*invdrmag;
  FLOAT invdr3   = invdrsqd*invdrmag;

  FLOAT mc = cell.m;
  pot += mc*invdrmag;

  mc *= invdr3;
  for (int k=0; k<2; k++) ac[k] += mc*dr[k];
  for (int k=0; k<2; k++) dphi[k] += mc*dr[k];
  q[0] += mc*(3.0*dr[0]*dr[0]*invdrsqd - 1);
  q[1] += mc*(3.0*dr[0]*dr[1]*invdrsqd);
  q[2] += mc*(3.0*dr[1]*dr[1]*invdrsqd - 1);

  if (HoT) {
    //                     (5  x_i   x_j   x_k  / r^2    - (x_i d_j,k + x_j d_i,k + x_k d_i,j))
    q2[0] += 3*mc*invdrsqd*(5*dr[0]*dr[0]*dr[0]*invdrsqd - (  dr[0]   +   dr[0]   +   dr[0]  ));
    q2[1] += 3*mc*invdrsqd*(5*dr[0]*dr[0]*dr[1]*invdrsqd - (    0     +     0     +   dr[1]  ));
    q2[2] += 3*mc*invdrsqd*(5*dr[0]*dr[1]*dr[1]*invdrsqd - (  dr[0]   +     0     +     0    ));
    q2[3] += 3*mc*invdrsqd*(5*dr[1]*dr[1]*dr[1]*invdrsqd - (  dr[1]   +   dr[1]   +   dr[1]  ));
  }
}
template<>
inline void FastMultipoleForces<3>::AddMonopoleContribution(const MultipoleMoment<3>& cell)  {
  FLOAT dr[3] ;
  for (int k=0; k<3; k++) dr[k] = cell.r[k] - rc[k];
  FLOAT invdrmag = sqrt((FLOAT) 1.0/DotProduct(dr,dr,3));
  FLOAT invdrsqd = invdrmag*invdrmag;
  FLOAT invdr3   = invdrsqd*invdrmag;

  FLOAT mc = cell.m;
  pot += mc*invdrmag;

  mc *= invdr3;
  for (int k=0; k<3; k++) ac[k] += mc*dr[k];
  for (int k=0; k<3; k++) dphi[k] += mc*dr[k];
  q[0] += mc*(3.0*dr[0]*dr[0]*invdrsqd - 1);
  q[1] += mc*(3.0*dr[0]*dr[1]*invdrsqd);
  q[2] += mc*(3.0*dr[1]*dr[1]*invdrsqd - 1);
  q[3] += mc*(3.0*dr[2]*dr[0]*invdrsqd);
  q[4] += mc*(3.0*dr[2]*dr[1]*invdrsqd);
  q[5] += mc*(3.0*dr[2]*dr[2]*invdrsqd - 1);

  if (HoT) {
    //                     (5  x_i   x_j   x_k  / r^2    - (x_i d_j,k + x_j d_i,k + x_k d_i,j))
    q2[0] += 3*mc*invdrsqd*(5*dr[0]*dr[0]*dr[0]*invdrsqd - (  dr[0]   +   dr[0]   +   dr[0]  ));
    q2[1] += 3*mc*invdrsqd*(5*dr[0]*dr[0]*dr[1]*invdrsqd - (    0     +     0     +   dr[1]  ));
    q2[2] += 3*mc*invdrsqd*(5*dr[0]*dr[1]*dr[1]*invdrsqd - (  dr[0]   +     0     +     0    ));
    q2[3] += 3*mc*invdrsqd*(5*dr[1]*dr[1]*dr[1]*invdrsqd - (  dr[1]   +   dr[1]   +   dr[1]  ));
    q2[4] += 3*mc*invdrsqd*(5*dr[0]*dr[0]*dr[2]*invdrsqd - (    0     +     0     +   dr[2]  ));
    q2[5] += 3*mc*invdrsqd*(5*dr[0]*dr[1]*dr[2]*invdrsqd - (    0     +     0     +     0    ));
    q2[6] += 3*mc*invdrsqd*(5*dr[0]*dr[2]*dr[2]*invdrsqd - (  dr[0]   +     0     +     0    ));
    q2[7] += 3*mc*invdrsqd*(5*dr[1]*dr[1]*dr[2]*invdrsqd - (    0     +     0     +   dr[2]  ));
    q2[8] += 3*mc*invdrsqd*(5*dr[1]*dr[2]*dr[2]*invdrsqd - (  dr[1]   +     0     +     0    ));
    q2[9] += 3*mc*invdrsqd*(5*dr[2]*dr[2]*dr[2]*invdrsqd - (  dr[2]   +   dr[2]   +   dr[2]  ));
  }
}

template<>
inline void FastMultipoleForces<1>::AddQuadrupoleContribution(const MultipoleMoment<1>& cell) {
  AddMonopoleContribution(cell);

  FLOAT dr[1] ;
  for (int k=0; k<1; k++) dr[k] = rc[k] - cell.r[k];
  FLOAT drsqd    = DotProduct(dr,dr,1) + small_number;
  FLOAT invdrsqd = (FLOAT) 1.0/drsqd;
  FLOAT invdrmag = sqrt(invdrsqd);
  FLOAT invdr5 = invdrsqd*invdrsqd*invdrmag ;

  FLOAT qscalar = cell.q[0]*dr[0]*dr[0] ;
  FLOAT qfactor = 2.5*qscalar*invdr5*invdrsqd;

  FLOAT qx[1];
  qx[0] = (cell.q[0]*dr[0])*invdr5;

  pot += 0.5*qscalar*invdr5;

  for (int k=0; k<1; k++) ac[k]   += qx[k] - qfactor*dr[k];
  for (int k=0; k<1; k++) dphi[k] += qx[k] - qfactor*dr[k];
  for (int k=0; k<1; k++) qx[k] *= 5*invdrsqd;

  q[0] += qfactor*(7*dr[0]*dr[0]*invdrsqd - 1);

  q[0] -= qx[0]*dr[0] + qx[0]*dr[0] - cell.q[0]*invdr5;
}
template<>
inline void FastMultipoleForces<2>::AddQuadrupoleContribution(const MultipoleMoment<2>& cell) {
  AddMonopoleContribution(cell);

  FLOAT dr[2] ;
  for (int k=0; k<2; k++) dr[k] = rc[k] - cell.r[k];
  FLOAT drsqd    = DotProduct(dr,dr,2) + small_number;
  FLOAT invdrsqd = (FLOAT) 1.0/drsqd;
  FLOAT invdrmag = sqrt(invdrsqd);
  FLOAT invdr5 = invdrsqd*invdrsqd*invdrmag ;

  FLOAT qscalar = cell.q[0]*dr[0]*dr[0] + cell.q[2]*dr[1]*dr[1] +
      2.0*cell.q[1]*dr[0]*dr[1];

  FLOAT qfactor = 2.5*qscalar*invdr5*invdrsqd;

  FLOAT qx[2];
  qx[0] = (cell.q[0]*dr[0] + cell.q[1]*dr[1])*invdr5;
  qx[1] = (cell.q[1]*dr[0] + cell.q[2]*dr[1])*invdr5;

  pot += 0.5*qscalar*invdr5;

  for (int k=0; k<2; k++) ac[k]   += qx[k] - qfactor*dr[k];
  for (int k=0; k<2; k++) dphi[k] += qx[k] - qfactor*dr[k];
  for (int k=0; k<2; k++) qx[k] *= 5*invdrsqd;

  q[0] += qfactor*(7*dr[0]*dr[0]*invdrsqd - 1);
  q[1] += qfactor*(7*dr[0]*dr[1]*invdrsqd);
  q[2] += qfactor*(7*dr[1]*dr[1]*invdrsqd - 1);

  q[0] -= qx[0]*dr[0] + qx[0]*dr[0] - cell.q[0]*invdr5;
  q[1] -= qx[0]*dr[1] + qx[1]*dr[0] - cell.q[1]*invdr5;
  q[2] -= qx[1]*dr[1] + qx[1]*dr[1] - cell.q[2]*invdr5;
}
template<>
inline void FastMultipoleForces<3>::AddQuadrupoleContribution(const MultipoleMoment<3>& cell) {
  AddMonopoleContribution(cell);

  FLOAT dr[3] ;
  for (int k=0; k<3; k++) dr[k] = rc[k] - cell.r[k];
  FLOAT drsqd    = DotProduct(dr,dr,3) + small_number;
  FLOAT invdrsqd = (FLOAT) 1.0/drsqd;
  FLOAT invdrmag = sqrt(invdrsqd);
  FLOAT invdr5 = invdrsqd*invdrsqd*invdrmag ;

  FLOAT qscalar = cell.q[0]*dr[0]*dr[0] + cell.q[2]*dr[1]*dr[1] -
      (cell.q[0] + cell.q[2])*dr[2]*dr[2] +
      2.0*(cell.q[1]*dr[0]*dr[1] + cell.q[3]*dr[0]*dr[2] + cell.q[4]*dr[1]*dr[2]);

  FLOAT qfactor = 2.5*qscalar*invdr5*invdrsqd;

  FLOAT qx[3];
  qx[0] = (cell.q[0]*dr[0] + cell.q[1]*dr[1] + cell.q[3]*dr[2])*invdr5;
  qx[1] = (cell.q[1]*dr[0] + cell.q[2]*dr[1] + cell.q[4]*dr[2])*invdr5;
  qx[2] = (cell.q[3]*dr[0] + cell.q[4]*dr[1] -(cell.q[0]+cell.q[2])*dr[2])*invdr5;

  pot += 0.5*qscalar*invdr5;

  for (int k=0; k<3; k++) ac[k]   += qx[k] - qfactor*dr[k];
  for (int k=0; k<3; k++) dphi[k] += qx[k] - qfactor*dr[k];
  for (int k=0; k<3; k++) qx[k] *= 5*invdrsqd;

  q[0] += qfactor*(7*dr[0]*dr[0]*invdrsqd - 1);
  q[1] += qfactor*(7*dr[0]*dr[1]*invdrsqd);
  q[2] += qfactor*(7*dr[1]*dr[1]*invdrsqd - 1);
  q[3] += qfactor*(7*dr[0]*dr[2]*invdrsqd);
  q[4] += qfactor*(7*dr[1]*dr[2]*invdrsqd);
  q[5] += qfactor*(7*dr[2]*dr[2]*invdrsqd - 1);

  q[0] -= qx[0]*dr[0] + qx[0]*dr[0] - cell.q[0]*invdr5;
  q[1] -= qx[0]*dr[1] + qx[1]*dr[0] - cell.q[1]*invdr5;
  q[2] -= qx[1]*dr[1] + qx[1]*dr[1] - cell.q[2]*invdr5;
  q[3] -= qx[0]*dr[2] + qx[2]*dr[0] - cell.q[3]*invdr5;
  q[4] -= qx[1]*dr[2] + qx[2]*dr[1] - cell.q[4]*invdr5;
  q[5] -= qx[2]*dr[2] + qx[2]*dr[2] + (cell.q[0] + cell.q[2])*invdr5;
}
//=================================================================================================
//  ApplyMonopoleForces
/// Apply the fast monopole forces at a given location.
//=================================================================================================
template<>
inline void FastMultipoleForces<1>::ApplyForcesTaylor(const FLOAT r[1], FLOAT agrav[1], FLOAT& gpot) const {
  FLOAT dr[1];
  for (int k=0; k<1; k++) dr[k] = r[k] - rc[k];

  agrav[0] += ac[0] + q[0]*dr[0];
  gpot += pot + dphi[0]*dr[0];

  // 2nd order terms:
  if (HoT) {
    agrav[0] += 0.5*q2[0]*dr[0]*dr[0];
    gpot += 0.5*q[0]*dr[0]*dr[0];
  }
}
template<>
inline void FastMultipoleForces<2>::ApplyForcesTaylor(const FLOAT r[2], FLOAT agrav[2], FLOAT& gpot) const {
  FLOAT dr[2];
  for (int k=0; k<2; k++) dr[k] = r[k] - rc[k];

  agrav[0] += ac[0] + q[0]*dr[0] + q[1]*dr[1];
  agrav[1] += ac[1] + q[1]*dr[0] + q[2]*dr[1];
  gpot += pot + dphi[0]*dr[0] + dphi[1]*dr[1];

  // 2nd order terms:
  if (HoT) {
    agrav[0] += 0.5*(q2[0]*dr[0]*dr[0] +  q2[2]*dr[1]*dr[1] + 2*q2[1]*dr[0]*dr[1]);
    agrav[1] += 0.5*(q2[1]*dr[0]*dr[0] +  q2[3]*dr[1]*dr[1] + 2*q2[2]*dr[0]*dr[1]);
    gpot += 0.5*(q[0]*dr[0]*dr[0] + q[2]*dr[1]*dr[1] + 2*q[1]*dr[0]*dr[1]);
  }
}
template<>
inline void FastMultipoleForces<3>::ApplyForcesTaylor(const FLOAT r[3], FLOAT agrav[3], FLOAT& gpot) const {
  FLOAT dr[3];
  for (int k=0; k<3; k++) dr[k] = r[k] - rc[k];

  agrav[0] += ac[0] + q[0]*dr[0] + q[1]*dr[1] + q[3]*dr[2];
  agrav[1] += ac[1] + q[1]*dr[0] + q[2]*dr[1] + q[4]*dr[2];
  agrav[2] += ac[2] + q[3]*dr[0] + q[4]*dr[1] + q[5]*dr[2];
  gpot += pot + dphi[0]*dr[0] + dphi[1]*dr[1] + dphi[2]*dr[2];

  // 2nd order terms:
  if (HoT) {
    agrav[0] += 0.5*(q2[0]*dr[0]*dr[0] +  q2[2]*dr[1]*dr[1] + q2[6]*dr[2]*dr[2] +
        2*(q2[1]*dr[0]*dr[1] + q2[4]*dr[0]*dr[2] + q2[5]*dr[1]*dr[2]));

    agrav[1] += 0.5*(q2[1]*dr[0]*dr[0] +  q2[3]*dr[1]*dr[1] + q2[8]*dr[2]*dr[2] +
        2*(q2[2]*dr[0]*dr[1] + q2[5]*dr[0]*dr[2] + q2[7]*dr[1]*dr[2]));

    agrav[2] += 0.5*(q2[4]*dr[0]*dr[0] +  q2[7]*dr[1]*dr[1] + q2[9]*dr[2]*dr[2] +
        2*(q2[5]*dr[0]*dr[1] + q2[6]*dr[0]*dr[2] + q2[8]*dr[1]*dr[2]));

    gpot += 0.5*(q[0]*dr[0]*dr[0] + q[2]*dr[1]*dr[1] + q[5]*dr[2]*dr[2] +
        2*(q[1]*dr[0]*dr[1] + q[3]*dr[0]*dr[2] + q[4]*dr[1]*dr[2]));
  }
}


//=================================================================================================
//  ComputeFastMonopoleForces
/// Compute the force on particle 'parti' due to all cells obtained in the
/// gravity tree walk.  Uses only monopole moments (i.e. COM) of the cell.
//=================================================================================================
template <int ndim, template<int> class ParticleType>
void ComputeFastMonopoleForces
 (int Nactive,                         ///< [in] No. of active particles
  int Ngravcell,                       ///< [in] No. of tree cells in list
  MultipoleMoment<ndim> *gravcell,     ///< [in] List of tree cell ids
  TreeCellBase<ndim> &cell,            ///< [in] Current cell pointer
  ParticleType<ndim> *activepart,      ///< [inout] Active Hydrodynamics particle array
  const ParticleTypeRegister& types)   ///< [in] Flags specifying which particles need grav forces
{

  FastMultipoleForces<ndim> monopole(cell.r) ;

  //-----------------------------------------------------------------------------------------------
  for (int cc=0; cc<Ngravcell; cc++) {
#ifndef MPI_PARALLEL
    assert(cell.id != gravcell[cc].id);
#endif
    monopole.AddMonopoleContribution(gravcell[cc]);
  }

  for (int j=0; j<Nactive; j++)
    if (types[activepart[j].ptype].self_gravity)
      monopole.ApplyForcesTaylor(activepart[j].r, activepart[j].atree, activepart[j].gpot) ;
  //-----------------------------------------------------------------------------------------------

  return;
}

//=================================================================================================
//  ComputeFastQuadrupoleForces
/// Compute the force on particle 'parti' due to all cells obtained in the
/// gravity tree walk.  Uses quadrupole moments (i.e. COM) of the cell.
//=================================================================================================
template <int ndim, template<int> class ParticleType>
void ComputeFastQuadrupoleForces
 (int Nactive,                         ///< [in] No. of active particles
  int Ngravcell,                       ///< [in] No. of tree cells in list
  MultipoleMoment<ndim> *gravcell,     ///< [in] List of tree cell ids
  TreeCellBase<ndim> &cell,            ///< [in] Current cell pointer
  ParticleType<ndim> *activepart,      ///< [inout] Active Hydrodynamics particle array
  const ParticleTypeRegister& types)   ///< [in] Flags specifying which particles need grav forces
{

  FastMultipoleForces<ndim> monopole(cell.r) ;

  //-----------------------------------------------------------------------------------------------
  for (int cc=0; cc<Ngravcell; cc++) {
#ifndef MPI_PARALLEL
    assert(cell.id != gravcell[cc].id);
#endif
    monopole.AddQuadrupoleContribution(gravcell[cc]);
  }

  for (int j=0; j<Nactive; j++)
    if (types[activepart[j].ptype].self_gravity)
      monopole.ApplyForcesTaylor(activepart[j].r, activepart[j].atree, activepart[j].gpot) ;
  //-----------------------------------------------------------------------------------------------

  return;
}


#endif // _MULTIPOLE_H_
