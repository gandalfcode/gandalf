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

#include "InlineFuncs.h"
#include "SimUnits.h"

template<int> class EosParticleProxy;

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

  OpacityTable(Parameters*, SimUnits *);
  ~OpacityTable();

  int GetIDens(const FLOAT);
  int GetITemp(const FLOAT);
  int GetIEner(const FLOAT, const FLOAT);
  int GetIEner(const FLOAT, const int);
  void GetKappa(const int, const int, FLOAT &, FLOAT &, FLOAT &);
  FLOAT GetEnergy(const int, const int);
  FLOAT GetMuBar(const EosParticleProxy<ndim> &part);
  FLOAT GetGamma(const EosParticleProxy<ndim> &part);
  FLOAT GetGamma1(const EosParticleProxy<ndim> &part);
  FLOAT GetEnergyFromPressure(const FLOAT, const FLOAT);
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


//=================================================================================================
//  OpacityTable::GetIDens()
/// GetIDens returns table index for density
//=================================================================================================
template <int ndim>
inline int OpacityTable<ndim>::GetIDens
 (const FLOAT rho)
{
  return getClosestIndex(eos_dens, eos_dens + ndens, rho);
}



//=================================================================================================
//  OpacityTable::GetITemp()
/// GetITemp returns table index for temperature
//=================================================================================================
template <int ndim>
inline int OpacityTable<ndim>::GetITemp
 (const FLOAT log_temp)
{
  return getClosestIndex(eos_temp, eos_temp + ntemp , log_temp);
}



//=================================================================================================
//  OpacityTable::GetIEner()
/// GetIEner returns table index for specific internal energy
//=================================================================================================
template <int ndim>
inline int OpacityTable<ndim>::GetIEner
 (const FLOAT u, const FLOAT rho)
{
  int idens = GetIDens(log10(rho));
  return getClosestIndex(eos_energy[idens], eos_energy[idens] + ntemp , u);
}
template <int ndim>
inline int OpacityTable<ndim>::GetIEner
 (const FLOAT u, const int idens)
{
  return getClosestIndex(eos_energy[idens], eos_energy[idens] + ntemp , u);
}


//=================================================================================================
//  OpacityTable::GetKappa()
/// GetKappa returns Kappa  for index of density and temp
//=================================================================================================
template <int ndim>
inline void OpacityTable<ndim>::GetKappa
 (const int idens,
  const int itemp,
  FLOAT &kappa,
  FLOAT &kappar,
  FLOAT &kappap)
{
  kappa = kappa_table[idens][itemp];
  kappap = kappap_table[idens][itemp];
  kappar = kappap;
}


//=================================================================================================
//  OpacityTable::GetEnergy()
/// GetEnergy returns Energy for index of density and temp
//=================================================================================================

template <int ndim>
inline FLOAT OpacityTable<ndim>::GetEnergy(const int idens, const int itemp)
{
  return eos_energy[idens][itemp];
}



//=================================================================================================
//  OpacityTable::GetMuBar()
/// GetMuBar returns MuBar for index of density and temp
//=================================================================================================
template <int ndim>
inline FLOAT OpacityTable<ndim>::GetMuBar
 (const EosParticleProxy<ndim> &part)
{
  int idens = GetIDens(part.rho);
  int iener = GetIEner(part.u, idens);
  return eos_mu[idens][iener];
}



//=================================================================================================
//  OpacityTable::GetGamma()
/// GetGamma returns gamma from P = (gamma - 1)*rho*u for index of density and u
//=================================================================================================
template <int ndim>
inline FLOAT OpacityTable<ndim>::GetGamma
 (const EosParticleProxy<ndim> &part)
{
  int idens = GetIDens(part.rho);
  int iener = GetIEner(part.u, idens);
  return eos_gamma[idens][iener];
}




//=================================================================================================
//  OpacityTable::GetGamma1()
/// GetGamma1 returns the first adiabatic index i.e. dln(P)/dln(rho) at constant entropy for index
/// of density and u
//=================================================================================================
template <int ndim>
inline FLOAT OpacityTable<ndim>::GetGamma1
 (const EosParticleProxy<ndim> &part)
{
  int idens = GetIDens(part.rho);
  int iener = GetIEner(part.u, idens);
  return eos_gamma1[idens][iener];
}

#endif
