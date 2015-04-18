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
    int var;
    FLOAT alpha;

    for (var=0; var<ndim+2; var++) {
      dW[var] = DotProduct(parti.grad[var], draux, ndim);
      if (dW[var] > 0.0) alpha = (parti.Wmax[var] - parti.Wprim[var])/dW[var];
      else if (dW[var] < 0.0) alpha = (parti.Wmin[var] - parti.Wprim[var])/dW[var];
      else alpha = 1.0;
      alpha = min(1.0, alpha);
      //alpha = 0.0;
      dW[var] = alpha*dW[var];
      for (int k=0; k<ndim; k++) gradW[var][k] = alpha*parti.grad[var][k];
      assert(alpha >= 0.0);
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
    int var;
    FLOAT alpha;
    const FLOAT theta = 0.5;

    for (var=0; var<ndim+2; var++) {
      dW[var] = DotProduct(parti.grad[var], draux, ndim);
      if (dW[var] > 0.0) alpha = max(theta*(partj.Wprim[var] - parti.Wprim[var])/dW[var], 0.0);
      else if (dW[var] < 0.0) alpha = max(theta*(partj.Wprim[var] - parti.Wprim[var])/dW[var], 0.0);
      else alpha = 1.0;
      alpha = min(1.0, alpha);
      //alpha = 0.0;
      dW[var] = alpha*dW[var];
      for (int k=0; k<ndim; k++) gradW[var][k] = alpha*parti.grad[var][k];
      assert(alpha >= 0.0);
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
    FLOAT alpha = 0.0;
    const FLOAT beta = 1.0;

    for (var=0; var<ndim+2; var++) {

      dW[var] = DotProduct(parti.grad[var], draux, ndim);
      alpha = min(1.0, beta*min((parti.Wmax[var] - parti.Wprim[var])/(parti.Wmidmax[var] - parti.Wprim[var]),
                                (parti.Wprim[var] - parti.Wmin[var])/(parti.Wprim[var] - parti.Wmidmin[var])));
      alpha = max(0.0, alpha);

      dW[var] = alpha*dW[var];
      for (int k=0; k<ndim; k++) gradW[var][k] = alpha*parti.grad[var][k];

      assert(alpha >= 0.0 && alpha <= 1.0);
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
    FLOAT alpha = 0.0;
    FLOAT dr[ndim];
    FLOAT phiminus;
    FLOAT phiplus;
    FLOAT phimid;
    const FLOAT beta = 1.0;
    const FLOAT psi1 = 0.5;
    const FLOAT psi2 = 0.25;

    //---------------------------------------------------------------------------------------------
    for (var=0; var<ndim+2; var++) {

      dW[var] = DotProduct(parti.grad[var], draux, ndim);
      alpha = min(1.0, beta*min((parti.Wmax[var] - parti.Wprim[var])/(parti.Wmidmax[var] - parti.Wprim[var]),
                                (parti.Wprim[var] - parti.Wmin[var])/(parti.Wprim[var] - parti.Wmidmin[var])));
      alpha = max(0.0, alpha);

      dW[var] = alpha*dW[var];
      for (int k=0; k<ndim; k++) gradW[var][k] = alpha*parti.grad[var][k];


      for (int k=0; k<ndim; k++) dr[k] = partj.r[k] - parti.r[k];
      FLOAT drmag = sqrt(DotProduct(dr, dr, ndim));

      const FLOAT delta1 = psi1*fabs(parti.Wprim[var] - partj.Wprim[var]);
      const FLOAT delta2 = psi2*fabs(parti.Wprim[var] - partj.Wprim[var]);
      const FLOAT phimin = min(parti.Wprim[var], partj.Wprim[var]);
      const FLOAT phimax = max(parti.Wprim[var], partj.Wprim[var]);
      const FLOAT phibar = parti.Wprim[var] + (partj.Wprim[var] - parti.Wprim[var])*
        sqrt(DotProduct(draux, draux, ndim))/drmag;
      const FLOAT phimid0 = parti.Wprim[var] + DotProduct(gradW[var], draux, ndim);

      if (sgn(phimin - delta1) == sgn(phimin)) {
        phiminus = phimin - delta1;
      }
      else {
        phiminus = phimin / (1.0 + delta1/fabs(phimin));
      }

      if (sgn(phimax + delta1) == sgn(phimax)) {
        phiplus = phimax + delta1;
      }
      else {
        phiplus = phimax / (1.0 + delta1/fabs(phimax));
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
      for (int k=0; k<ndim; k++) gradW[var][k] = dW[var]*draux[k]/drsqd;

    }
    //---------------------------------------------------------------------------------------------
  }

};
#endif
