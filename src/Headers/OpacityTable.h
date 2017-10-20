//=================================================================================================
//  OpacityTable.h
//  Class definitions of radiative feedback methods.
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

#ifndef _OPACITY_TABLE_H_
#define _OPACITY_TABLE_H_

#include "Particle.h"
#include "SimUnits.h"

//=================================================================================================
//  Class OpacityTable
/// \brief   Container for opacity table data with helper functions.
/// \details Allows the reading of an opacity table (e.g. eos.bell.cc.dat) via the constructor.
///          The values of gamma and mu_bar can then be accessed using an EOS class to calculate
///          thermodynamical properties.
/// \author  A. P. Mercer
/// \date    18/10/2017
//=================================================================================================
template <int ndim>
class OpacityTable
{
 public:

  OpacityTable(string, SimUnits *);
  ~OpacityTable();

  int GetIDens(const FLOAT);
  int GetITemp(const FLOAT);
  int GetIEner(const FLOAT, const FLOAT);
  int GetIEner(const FLOAT, const int);
  void GetKappa(const int, const int, FLOAT &, FLOAT &, FLOAT &);
  FLOAT GetEnergy(const int, const int);
  FLOAT GetMuBar(Particle<ndim> &part);
  FLOAT GetGamma(Particle<ndim> &part);
  FLOAT GetGamma1(Particle<ndim> &part);
  //-----------------------------------------------------------------------------------------------

  int ndens;
  int ntemp;
  FLOAT fcol;
  FLOAT *eos_dens;
  FLOAT *eos_temp ;
  FLOAT **eos_energy;
  FLOAT **eos_mu;
  FLOAT **kappa_table;
  FLOAT **kappar_table;
  FLOAT **kappap_table;
  FLOAT **eos_gamma;
  FLOAT **eos_gamma1;   // First adiabatic index
};

#endif
