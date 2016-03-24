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
#include "SimUnits.h"
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

  virtual void AddExternalPotential(const FLOAT *, const FLOAT *, FLOAT *, FLOAT *, FLOAT &) = 0;

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

  void AddExternalPotential(const FLOAT *, const FLOAT *, FLOAT *, FLOAT *, FLOAT &) {}

};



//=================================================================================================
//  Class VerticalPotential
/// \brief   Add simple constant gravitational field potential
/// \details Add simple constant gravitational field potential
/// \author  D. A. Hubber
/// \date    14/03/2015
//=================================================================================================
template <int ndim>
class VerticalPotential : public ExternalPotential<ndim>
{
public:

  const int kgrav;                     ///< 'Direction' of grav. field
  const FLOAT avert;                   ///< Size (plus sign) of grav. field
  const FLOAT rzero;                   ///< 'Zero' -height of potential


  VerticalPotential(int _kgrav, FLOAT _avert, FLOAT _rzero) :
    kgrav(_kgrav), avert(_avert), rzero(_rzero) {}
  ~VerticalPotential();

  void AddExternalPotential
   (const FLOAT rp[ndim],              ///< Position of particle
    const FLOAT vp[ndim],              ///< Velocity of particle
    FLOAT ap[ndim],                    ///< Acceleration of particle
    FLOAT adotp[ndim],                 ///< 'Jerk' of particle
    FLOAT &potp)                       ///< Potential of particle
  {
    ap[kgrav]    += avert;
    adotp[kgrav] += (FLOAT) 0.0;
    potp         += (rp[kgrav] - rzero)*avert;

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

  const FLOAT mplummer;                ///< Mass of Plummer sphere
  const FLOAT rplummer;                ///< Core radius of Plummer sphere


  PlummerPotential(FLOAT mplummeraux, FLOAT rplummeraux) :
    mplummer(mplummeraux), rplummer(rplummeraux) {}
  ~PlummerPotential();


  void AddExternalPotential
   (const FLOAT rp[ndim],              ///< Position of particle
    const FLOAT vp[ndim],              ///< Velocity of particle
    FLOAT ap[ndim],                    ///< Acceleration of particle
    FLOAT adotp[ndim],                 ///< 'Jerk' of particle
    FLOAT &potp)                       ///< Potential of particle
  {
    int k;                             // Dimension counter
    FLOAT drsqd;                       // Distance squared
    FLOAT dvdr;                        // Dot product of velocity and position

    drsqd = DotProduct(rp,rp,ndim);
    dvdr = DotProduct(rp,vp,ndim);
    for (k=0; k<ndim; k++) ap[k] -= mplummer*rp[k]*pow(drsqd + rplummer*rplummer,-(FLOAT) 1.5);
    for (k=0; k<ndim; k++) adotp[k] += (FLOAT) 3.0*mplummer*
      pow(drsqd + rplummer*rplummer, -(FLOAT) 2.5)*dvdr*rp[k]
      - mplummer*pow(drsqd + rplummer*rplummer, -(FLOAT) 1.5)*vp[k];
    potp += (FLOAT) 2.0*mplummer*pow(drsqd + rplummer*rplummer, -(FLOAT) 0.5);

    return;
  }

};



//=================================================================================================
//  Class SilccPotential
/// \brief   Add potential, acceleration and jerk for Silcc simulations.
/// \details Add potential, acceleration and jerk for Silcc simulations.
/// \author  D. A. Hubber
/// \date    10/03/2014
//=================================================================================================
template <int ndim>
class SilccPotential : public ExternalPotential<ndim>
{
public:

  FLOAT rho_star;
  FLOAT sigma_star;
  FLOAT z_d;

  SilccPotential(FLOAT _sigma_star, FLOAT _z_d, SimUnits &simunits)
  {
    sigma_star = _sigma_star/simunits.sigma.outscale;
    z_d = _z_d/simunits.r.outscale,
    rho_star = (FLOAT) 0.25*sigma_star/z_d;
  }
  ~SilccPotential();


  void AddExternalPotential
   (const FLOAT rp[ndim],              ///< Position of particle
    const FLOAT vp[ndim],              ///< Velocity of particle
    FLOAT ap[ndim],                    ///< Acceleration of particle
    FLOAT adotp[ndim],                 ///< 'Jerk' of particle
    FLOAT &potp)                       ///< Potential of particle
  {
    return;
  }

};
#endif
