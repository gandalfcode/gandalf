//=================================================================================================
//  ExternalPotential.h
//  Class definitions for all external potential fields.
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


#ifndef _EXTERNAL_POTENTIAL_H_
#define _EXTERNAL_POTENTIAL_H_


#include <string>
#include "Precision.h"
#include "Constants.h"
#include "InlineFuncs.h"
using namespace std;



//=================================================================================================
//  Class ExternalPotential
/// \brief   Class to compute and return all terms of external potential fields
/// \details Class to compute and return all terms of external potential fields
/// \author  D. A. Hubber
/// \date    10/03/2014
//=================================================================================================
template <int ndim>
class ExternalPotential
{
 public:

  ExternalPotential() {};
  ~ExternalPotential() {};

  virtual void AddExternalPotential(DOUBLE *, DOUBLE *, DOUBLE *, DOUBLE *, DOUBLE &) = 0;

};



//=================================================================================================
//  Class NullPotential
/// \brief   Null class when there is no external potential field to add
/// \details Null class when there is no external potential field to add
/// \author  D. A. Hubber
/// \date    10/03/2014
//=================================================================================================
template <int ndim>
class NullPotential : public ExternalPotential<ndim>
{
 public:

  NullPotential() {};
  ~NullPotential() {};

  void AddExternalPotential(DOUBLE *, DOUBLE *, DOUBLE *, DOUBLE *, DOUBLE &) {}

};



//=================================================================================================
//  Class VerticalPotential
/// \brief   Add simple constant gravitational field potential
/// \details ...
/// \author  D. A. Hubber
/// \date    14/03/2015
//=================================================================================================
template <int ndim>
class VerticalPotential : public ExternalPotential<ndim>
{
 public:

  VerticalPotential(int _kgrav, DOUBLE _avert, DOUBLE _rzero) :
    kgrav(_kgrav), avert(_avert), rzero(_rzero) {}
  ~VerticalPotential();

  const int kgrav;                     ///< ..
  const DOUBLE avert;                  ///< ..
  const DOUBLE rzero;                  ///< ..


  void AddExternalPotential
   (DOUBLE rp[ndim],                   ///< Position of particle
    DOUBLE vp[ndim],                   ///< Velocity of particle
    DOUBLE ap[ndim],                   ///< Acceleration of particle
    DOUBLE adotp[ndim],                ///< 'Jerk' of particle
    DOUBLE &potp)                      ///< Potential of particle
  {
    ap[kgrav]    += avert;
    adotp[kgrav] += 0.0;
    potp         += (rp[kgrav] - rzero)* avert;

    return;
  }

};



//=================================================================================================
//  Class PlummerPotential
/// \brief   Add potential, acceleration and jerk for background Plummer potential.
/// \details Add potential, acceleration and jerk for background Plummer potential.
/// \author  D. A. Hubber
/// \date    10/03/2014
//=================================================================================================
template <int ndim>
class PlummerPotential : public ExternalPotential<ndim>
{
 public:

  PlummerPotential(DOUBLE mplummeraux, DOUBLE rplummeraux) :
    mplummer(mplummeraux), rplummer(rplummeraux) {}
  ~PlummerPotential();


  const DOUBLE mplummer;               ///< Mass of Plummer sphere
  const DOUBLE rplummer;               ///< Core radius of Plummer sphere


  void AddExternalPotential
   (DOUBLE rp[ndim],                   ///< Position of particle
    DOUBLE vp[ndim],                   ///< Velocity of particle
    DOUBLE ap[ndim],                   ///< Acceleration of particle
    DOUBLE adotp[ndim],                ///< 'Jerk' of particle
    DOUBLE &potp)                      ///< Potential of particle
  {
    int k;                             // Dimension counter
    DOUBLE drsqd;                      // Distance squared
    DOUBLE dvdr;                       // Dot product of velocity and position

    drsqd = DotProduct(rp,rp,ndim);
    dvdr = DotProduct(rp,vp,ndim);
    for (k=0; k<ndim; k++) ap[k] -= mplummer*rp[k]*pow(drsqd + rplummer*rplummer,-1.5);
    for (k=0; k<ndim; k++) adotp[k] += 3.0*mplummer*pow(drsqd + rplummer*rplummer,-2.5)*dvdr*rp[k]
      - mplummer*pow(drsqd + rplummer*rplummer,-1.5)*vp[k];
    potp += 2.0*mplummer*pow(drsqd + rplummer*rplummer,-0.5);

    return;
  }

};
#endif
