//=================================================================================================
//  SlopeLimiter.h
//  Contains all routines for ...
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


#ifndef _SLOPE_LIMITER_H_
#define _SLOPE_LIMITER_H_


#include <assert.h>
#include <iostream>
#include <math.h>
#include "Exception.h"
#include "InlineFuncs.h"
#include "Precision.h"
using namespace std;



//=================================================================================================
//  Class SlopeLimiter
/// \brief   Parent class for all slope limiters
/// \details ...
/// \author  D. A. Hubber
/// \date    23/03/2015
//=================================================================================================
template <int ndim, template<int> class ParticleType>
class SlopeLimiter
{
 public:

  SlopeLimiter() {};
  ~SlopeLimiter() {};

  virtual void ComputeLimitedSlopes(ParticleType<ndim> &parti, ParticleType<ndim> &partj,
                                    FLOAT draux[ndim], FLOAT gradW[ndim+2][ndim], FLOAT dW[ndim+2]) = 0;
  virtual void ComputeAlphaSlope(ParticleType<ndim> &parti, ParticleType<ndim> &partj, FLOAT draux[ndim]) {};
};



//=================================================================================================
//  Class ZeroSlopeLimiter
/// \brief   Imposes 1st-order (Godunov method) by zeroing slopes.
/// \details Imposes 1st-order (Godunov method) by zeroing slopes.
/// \author  D. A. Hubber
/// \date    28/08/2015
//=================================================================================================
template <int ndim, template<int> class ParticleType>
class ZeroSlopeLimiter : public SlopeLimiter<ndim,ParticleType>
{
 public:

  ZeroSlopeLimiter() {};
  ~ZeroSlopeLimiter() {};

  //===============================================================================================
  void ComputeLimitedSlopes(ParticleType<ndim> &parti, ParticleType<ndim> &partj,
                            FLOAT draux[ndim], FLOAT gradW[ndim+2][ndim], FLOAT dW[ndim+2])
  {
    for (int var=0; var<ndim+2; var++) {
      dW[var] = (FLOAT) 0.0;
      for (int k=0; k<ndim; k++) gradW[var][k] = (FLOAT) 0.0;
    }
  }

};



//=================================================================================================
//  Class NullLimiter
/// \brief   Null slope limiter.  Extrapolates variables fully without limiting their values.
/// \details Null slope limiter.  Extrapolates variables fully without limiting their values.
/// \author  D. A. Hubber
/// \date    23/03/2015
//=================================================================================================
template <int ndim, template<int> class ParticleType>
class NullLimiter : public SlopeLimiter<ndim,ParticleType>
{
 public:

  NullLimiter() {};
  ~NullLimiter() {};

  //===============================================================================================
  void ComputeLimitedSlopes(ParticleType<ndim> &parti, ParticleType<ndim> &partj,
                            FLOAT draux[ndim], FLOAT gradW[ndim+2][ndim], FLOAT dW[ndim+2])
  {
    for (int var=0; var<ndim+2; var++) {
      dW[var] = DotProduct(parti.grad[var], draux, ndim);
      for (int k=0; k<ndim; k++) gradW[var][k] = parti.grad[var][k];
    }
  }

};



//=================================================================================================
//  Class Springel2009Limiter
/// \brief   ...
/// \details ...
/// \author  D. A. Hubber
/// \date    23/03/2015
//=================================================================================================
template <int ndim, template<int> class ParticleType>
class Springel2009Limiter : public SlopeLimiter<ndim,ParticleType>
{
 public:

  Springel2009Limiter() {};
  ~Springel2009Limiter() {};

  //===============================================================================================
  void ComputeLimitedSlopes(ParticleType<ndim> &parti, ParticleType<ndim> &partj,
                            FLOAT draux[ndim], FLOAT gradW[ndim+2][ndim], FLOAT dW[ndim+2])
  {
    for (int var=0; var<ndim+2; var++) {
      dW[var] = DotProduct(parti.grad[var], draux, ndim);
      dW[var] = parti.alpha_slope[var]*dW[var];
      for (int k=0; k<ndim; k++) gradW[var][k] = parti.alpha_slope[var]*parti.grad[var][k];
    }
  }

  virtual void ComputeAlphaSlope(ParticleType<ndim> &parti, ParticleType<ndim> &partj, FLOAT draux[ndim])
  {
    FLOAT psi;
    FLOAT dW[ndim+2];

    for (int var=0; var<ndim+2; var++) {
      dW[var] = DotProduct(parti.grad[var], draux, ndim);
      if (dW[var] > (FLOAT) 0.0) psi = (parti.Wmax[var] - parti.Wprim[var])/dW[var];
      else if (dW[var] < (FLOAT) 0.0) psi = (parti.Wmin[var] - parti.Wprim[var])/dW[var];
      else psi = (FLOAT) 1.0;
      parti.alpha_slope[var] = min(parti.alpha_slope[var], psi);
    }
  }

};



//=================================================================================================
//  Class TESS2011Limiter
/// \brief   ...
/// \details ...
/// \author  D. A. Hubber
/// \date    23/03/2015
//=================================================================================================
template <int ndim, template<int> class ParticleType>
class TESS2011Limiter : public SlopeLimiter<ndim,ParticleType>
{
 public:

  TESS2011Limiter() {};
  ~TESS2011Limiter() {};

  //===============================================================================================
  void ComputeLimitedSlopes(ParticleType<ndim> &parti, ParticleType<ndim> &partj,
                            FLOAT draux[ndim], FLOAT gradW[ndim+2][ndim], FLOAT dW[ndim+2])
  {
    for (int var=0; var<ndim+2; var++) {
      dW[var] = DotProduct(parti.grad[var], draux, ndim);
      dW[var] = parti.alpha_slope[var]*dW[var];
      for (int k=0; k<ndim; k++) gradW[var][k] = parti.alpha_slope[var]*parti.grad[var][k];
    }
  }

  virtual void ComputeAlphaSlope(ParticleType<ndim> &parti, ParticleType<ndim> &partj, FLOAT draux[ndim])
  {
    const FLOAT theta = (FLOAT) 0.5;
    FLOAT psi;
    FLOAT dW[ndim+2];

    for (int var=0; var<ndim+2; var++) {
      dW[var] = DotProduct(parti.grad[var], draux, ndim);
      if (dW[var] > (FLOAT) 0.0) psi = max(theta*(partj.Wprim[var] - parti.Wprim[var])/dW[var], (FLOAT) 0.0);
      else if (dW[var] < 0.0) psi = max(theta*(partj.Wprim[var] - parti.Wprim[var])/dW[var], (FLOAT) 0.0);
      else psi = (FLOAT) 1.0;
      parti.alpha_slope[var] = min(parti.alpha_slope[var], psi);
    }
  }

};



//=================================================================================================
//  Class Balsara2004Limiter
/// \brief   ...
/// \details ...
/// \author  D. A. Hubber
/// \date    23/03/2015
//=================================================================================================
template <int ndim, template<int> class ParticleType>
class Balsara2004Limiter : public SlopeLimiter<ndim,ParticleType>
{
 public:

  Balsara2004Limiter() {};
  ~Balsara2004Limiter() {};

  //===============================================================================================
  void ComputeLimitedSlopes(ParticleType<ndim> &parti, ParticleType<ndim> &partj,
                            FLOAT draux[ndim], FLOAT gradW[ndim+2][ndim], FLOAT dW[ndim+2])
  {
    int var;
    FLOAT alpha = (FLOAT) 0.0;
    const FLOAT beta = (FLOAT) 1.0;

    for (var=0; var<ndim+2; var++) {

      dW[var] = DotProduct(parti.grad[var], draux, ndim);
      alpha = min((FLOAT) 1.0, beta*min((parti.Wmax[var] - parti.Wprim[var])/(parti.Wmidmax[var] - parti.Wprim[var]),
                                        (parti.Wprim[var] - parti.Wmin[var])/(parti.Wprim[var] - parti.Wmidmin[var])));
      alpha = max((FLOAT) 0.0, alpha);

      dW[var] = alpha*dW[var];
      for (int k=0; k<ndim; k++) gradW[var][k] = alpha*parti.grad[var][k];

      assert(alpha >= (FLOAT) 0.0 && alpha <= (FLOAT) 1.0);
    }
  }

};



//=================================================================================================
//  Class GizmoLimiter
/// \brief   ...
/// \details ...
/// \author  D. A. Hubber
/// \date    23/03/2015
//=================================================================================================
template <int ndim, template<int> class ParticleType>
class GizmoLimiter : public SlopeLimiter<ndim,ParticleType>
{
 public:

  GizmoLimiter() {};
  ~GizmoLimiter() {};

  //===============================================================================================
  void ComputeLimitedSlopes(ParticleType<ndim> &parti, ParticleType<ndim> &partj,
                            FLOAT draux[ndim], FLOAT gradW[ndim+2][ndim], FLOAT dW[ndim+2])
  {
    int var;
    FLOAT alpha = (FLOAT) 0.0;
    FLOAT dr[ndim];
    FLOAT phiminus;
    FLOAT phiplus;
    FLOAT phimid;
    const FLOAT beta = (FLOAT) 2.0;
    const FLOAT psi1 = (FLOAT) 0.5;
    const FLOAT psi2 = (FLOAT) 0.375;  //0.25;
    FLOAT gradPerp[ndim];
    FLOAT dr_unit[ndim];

    const FLOAT alim = 0.25;
    const FLOAT stol = 0.0;


    //---------------------------------------------------------------------------------------------
    for (var=0; var<ndim+2; var++) {

      for (int k=0; k<ndim; k++) dr[k] = partj.r[k] - parti.r[k];
      FLOAT drmag = sqrt(DotProduct(dr, dr, ndim));
      for (int k=0; k<ndim; k++) dr_unit[k] = dr[k]/drmag;

      dW[var] = DotProduct(parti.grad[var], draux, ndim);
      for (int k=0; k<ndim; k++) gradPerp[k] = parti.grad[var][k] - DotProduct(parti.grad[var], dr_unit, ndim)*dr_unit[k];




      /*FLOAT d_abs = 0.0;
      d_abs += DotProduct(parti.grad[var], parti.grad[var], ndim);

      if (d_abs > 0.0) {
        FLOAT cfac = 1.0 / (alim*sqrt(parti.hrangesqd*d_abs));
        FLOAT fabs_max = fabs(parti.Wmax[var]);
        FLOAT fabs_min = fabs(parti.Wmin[var]);
        FLOAT abs_min = min(fabs_max, fabs_min);
        if (stol > 0.0) {

        }
        else {
          cfac *= abs_min;
        }
        if (cfac < 1.0) {
          for (int k=0; k<ndim; k++) gradW[var][k] = parti.grad[var][k]*cfac;
          dW[var] = DotProduct(gradW[var], draux, ndim);
        }
        else {
          for (int k=0; k<ndim; k++) gradW[var][k] = parti.grad[var][k];
          dW[var] = DotProduct(gradW[var], draux, ndim);
        }
      }*/

      alpha = min((FLOAT) 1.0, beta*min((parti.Wmax[var] - parti.Wprim[var])/(parti.Wmidmax[var] - parti.Wprim[var]),
                                        (parti.Wprim[var] - parti.Wmin[var])/(parti.Wprim[var] - parti.Wmidmin[var])));
      alpha = max((FLOAT) 0.0, alpha);
      alpha = 1.0;

      dW[var] = alpha*dW[var];
      for (int k=0; k<ndim; k++) gradW[var][k] = alpha*parti.grad[var][k];



      const FLOAT delta1 = psi1*fabs(parti.Wprim[var] - partj.Wprim[var]);
      const FLOAT delta2 = psi2*fabs(parti.Wprim[var] - partj.Wprim[var]);
      const FLOAT phimin = min(parti.Wprim[var], partj.Wprim[var]);
      const FLOAT phimax = max(parti.Wprim[var], partj.Wprim[var]);
      const FLOAT phibar = parti.Wprim[var] + (partj.Wprim[var] - parti.Wprim[var])*
        sqrt(DotProduct(draux, draux, ndim))/drmag;
      const FLOAT phimid0 = parti.Wprim[var] + dW[var]; //DotProduct(gradW[var], draux, ndim);

      if (sgn(phimin - delta1) == sgn(phimin)) {
        phiminus = phimin - delta1;
      }
      else {
        phiminus = phimin / ((FLOAT) 1.0 + delta1/fabs(phimin));
      }

      if (sgn(phimax + delta1) == sgn(phimax)) {
        phiplus = phimax + delta1;
      }
      else {
        phiplus = phimax / ((FLOAT) 1.0 + delta1/fabs(phimax));
      }

      if (parti.Wprim[var] < partj.Wprim[var]) {
        phimid = max(phiminus, min(phibar + delta2, phimid0));
      }
      else if (parti.Wprim[var] > partj.Wprim[var]) {
        phimid = min(phiplus, max(phibar - delta2, phimid0));
      }
      else {
        phimid = parti.Wprim[var];
      }

      FLOAT drsqd = DotProduct(draux, draux, ndim);
      dW[var] = phimid - parti.Wprim[var];
      for (int k=0; k<ndim; k++) gradW[var][k] = gradPerp[k] + dW[var]*draux[k]/drsqd;

    }
    //---------------------------------------------------------------------------------------------
  }

};




//=================================================================================================
//  Class GizmoLimiter
/// \brief   ...
/// \details ...
/// \author  D. A. Hubber
/// \date    23/03/2015
//=================================================================================================
template <int ndim, template<int> class ParticleType>
class Gizmo2Limiter : public SlopeLimiter<ndim,ParticleType>
{
 public:

  Gizmo2Limiter() {};
  ~Gizmo2Limiter() {};

  //===============================================================================================
  void ComputeLimitedSlopes(ParticleType<ndim> &parti, ParticleType<ndim> &partj,
                            FLOAT draux[ndim], FLOAT gradW[ndim+2][ndim], FLOAT dW[ndim+2])
  {
    int var;
    FLOAT alpha = (FLOAT) 0.0;
    FLOAT dr[ndim];
    FLOAT phiminus;
    FLOAT phiplus;
    FLOAT phimid;
    const FLOAT beta = (FLOAT) 2.0;
    const FLOAT psi1 = (FLOAT) 0.5;
    const FLOAT psi2 = (FLOAT) 0.375;  //0.25;
    FLOAT gradPerp[ndim];
    FLOAT dr_unit[ndim];

    const FLOAT alim = 0.25;
    const FLOAT stol = 0.0;


    //---------------------------------------------------------------------------------------------
    for (var=0; var<ndim+2; var++) {

      for (int k=0; k<ndim; k++) dr[k] = partj.r[k] - parti.r[k];
      FLOAT drmag = sqrt(DotProduct(dr, dr, ndim));
      for (int k=0; k<ndim; k++) dr_unit[k] = dr[k]/drmag;

      dW[var] = DotProduct(parti.grad[var], draux, ndim);
      for (int k=0; k<ndim; k++) gradPerp[k] = 0.0; //parti.grad[var][k] - DotProduct(parti.grad[var], dr_unit, ndim)*dr_unit[k];



      alpha = min((FLOAT) 1.0, beta*min((parti.Wmax[var] - parti.Wprim[var])/(parti.Wmidmax[var] - parti.Wprim[var]),
                                        (parti.Wprim[var] - parti.Wmin[var])/(parti.Wprim[var] - parti.Wmidmin[var])));
      alpha = max((FLOAT) 0.0, alpha);
      //alpha = 1.0;

      dW[var] = alpha*dW[var];
      for (int k=0; k<ndim; k++) gradW[var][k] = alpha*parti.grad[var][k];



      const FLOAT delta1 = psi1*fabs(parti.Wprim[var] - partj.Wprim[var]);
      const FLOAT delta2 = psi2*fabs(parti.Wprim[var] - partj.Wprim[var]);
      const FLOAT phimin = min(parti.Wprim[var], partj.Wprim[var]);
      const FLOAT phimax = max(parti.Wprim[var], partj.Wprim[var]);
      const FLOAT phibar = parti.Wprim[var] + (partj.Wprim[var] - parti.Wprim[var])*
        sqrt(DotProduct(draux, draux, ndim))/drmag;
      const FLOAT phimid0 = parti.Wprim[var] + dW[var]; //DotProduct(gradW[var], draux, ndim);

      if (sgn(phimin - delta1) == sgn(phimin)) {
        phiminus = phimin - delta1;
      }
      else {
        phiminus = phimin / ((FLOAT) 1.0 + delta1/fabs(phimin));
      }

      if (sgn(phimax + delta1) == sgn(phimax)) {
        phiplus = phimax + delta1;
      }
      else {
        phiplus = phimax / ((FLOAT) 1.0 + delta1/fabs(phimax));
      }

      if (parti.Wprim[var] < partj.Wprim[var]) {
        phimid = max(phiminus, min(phibar + delta2, phimid0));
      }
      else if (parti.Wprim[var] > partj.Wprim[var]) {
        phimid = min(phiplus, max(phibar - delta2, phimid0));
      }
      else {
        phimid = parti.Wprim[var];
      }

      FLOAT drsqd = DotProduct(draux, draux, ndim);
      dW[var] = phimid - parti.Wprim[var];
      for (int k=0; k<ndim; k++) gradW[var][k] = gradPerp[k] + dW[var]*draux[k]/drsqd;

    }
    //---------------------------------------------------------------------------------------------
  }

};



//=================================================================================================
//  Class MinModLimiter
/// \brief   Null slope limiter.  Extrapolates variables fully without limiting their values.
/// \details Null slope limiter.  Extrapolates variables fully without limiting their values.
/// \author  D. A. Hubber
/// \date    23/03/2015
//=================================================================================================
template <int ndim, template<int> class ParticleType>
class MinModLimiter : public SlopeLimiter<ndim,ParticleType>
{
 public:

  MinModLimiter() {};
  ~MinModLimiter() {};

  //===============================================================================================
  void ComputeLimitedSlopes(ParticleType<ndim> &parti, ParticleType<ndim> &partj,
                            FLOAT draux[ndim], FLOAT gradW[ndim+2][ndim], FLOAT dW[ndim+2])
  {
    for (int var=0; var<ndim+2; var++) {

      FLOAT dr[ndim], dr_unit[ndim];
      for (int k=0; k<ndim; k++) dr[k] = partj.r[k] - parti.r[k];
      FLOAT drmag = sqrt(DotProduct(dr, dr, ndim));
      for (int k=0; k<ndim; k++) dr_unit[k] = dr[k] / drmag;

      FLOAT D = DotProduct(parti.grad[var], dr_unit, ndim);

      if (D == (FLOAT) 0.0 && partj.Wprim[var] == parti.Wprim[var]) {
        dW[var] = (FLOAT) 0.0;
        for (int k=0; k<ndim; k++) gradW[var][k] = (FLOAT) 0.0;
        continue;
      }

      FLOAT r = (partj.Wprim[var] - parti.Wprim[var]) /
        (parti.Wprim[var] - partj.Wprim[var] + 2.0*D*drmag);
      r = max(r, (FLOAT) 0.0);

      //FLOAT phi = min(2.0/(1.0 + r), 2.0*r/(1.0 + r));
      //FLOAT phi = min(1.0, min(4.0/(1.0 + r), 4.0*r/(1.0 + r)));
      //FLOAT phi = 2.0*r/(r*r + 1.0);
      FLOAT phi = (FLOAT) 4.0*r/pow(r + (FLOAT) 1.0, 2);

      dW[var] = (FLOAT) 0.5*phi*D*drmag;
      for (int k=0; k<ndim; k++) gradW[var][k] = dW[var]*dr_unit[k];

    }
  }

};
#endif
