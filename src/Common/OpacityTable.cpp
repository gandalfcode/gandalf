//=================================================================================================
//  OpacityTable.cpp
//  Contains functions for Rad-WS radiation cooling scheme (Stamatellos et al. 2007) class.
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


#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <math.h>

#include "../Headers/Integration.h"
#include "Constants.h"
#include "Debug.h"
#include "OpacityTable.h"
#include "Exception.h"
#include "Parameters.h"
#include "Particle.h"
#include "SimUnits.h"
using namespace std;



//=================================================================================================
//  OpacityTable::OpacityTable()
/// OpacityTable class constructor
//=================================================================================================
template <int ndim>
OpacityTable<ndim>::OpacityTable
 (string radws_table,
  SimUnits *simunits)
{
  int i, j, l;
  DOUBLE eos_dens_, eos_temp_, eos_energy_, eos_mu_, kappa_, kappar_, kappap_, eos_gamma_, eos_gamma1_;
  string line;

  // Check for correct opacity units
  if (simunits->kappa.outunit != "cm2_g") {
    ExceptionHandler::getIstance().raise("Error: Wrong units for opacity, use cm2_g");
  }

  ifstream file;
  file.open(radws_table.c_str(), ios::in);
  if (!file.is_open()) {
    ExceptionHandler::getIstance().raise("Error: Cannot open opacity table!");
  }

  //-----------------------------------------------------------------------------------------------
  if (file.good()) {
    do {
      getline(file, line);
    } while (line[0] == '#');
    istringstream istr(line);
    istr >> ndens >> ntemp >> fcol;

    eos_dens     = new FLOAT[ndens];
    eos_temp     = new FLOAT[ntemp];
    eos_energy   = new FLOAT*[ndens];
    eos_mu       = new FLOAT*[ndens];
    kappa_table  = new FLOAT*[ndens];
    kappar_table = new FLOAT*[ndens];
    kappap_table = new FLOAT*[ndens];
    eos_gamma    = new FLOAT*[ndens];
    eos_gamma1   = new FLOAT*[ndens];

    for (i = 0; i < ndens; ++i){
      eos_energy[i]   = new FLOAT[ntemp];
      eos_mu[i]       = new FLOAT[ntemp];
      kappa_table[i]  = new FLOAT[ntemp];
      kappar_table[i] = new FLOAT[ntemp];
      kappap_table[i] = new FLOAT[ntemp];
      eos_gamma[i]    = new FLOAT[ntemp];
      eos_gamma1[i]   = new FLOAT[ntemp];
    }

    // read table
    i = 0;
    l = 0;
    j = 0;

    //---------------------------------------------------------------------------------------------
    while (getline(file, line)) {
      istringstream istr(line);

      if (istr >> eos_dens_ >> eos_temp_ >> eos_energy_ >> eos_mu_>> kappa_ >> kappar_ >> kappap_ >> eos_gamma_ >> eos_gamma1_) {
        eos_energy[i][j]   = eos_energy_/(simunits->u.outscale * simunits-> u.outcgs);
        eos_mu[i][j]       = eos_mu_;
        kappa_table[i][j]  = kappa_/(simunits->kappa.outscale * simunits-> kappa.outcgs);
        kappar_table[i][j] = kappar_/(simunits->kappa.outscale * simunits-> kappa.outcgs);
        kappap_table[i][j] = kappap_/(simunits->kappa.outscale * simunits-> kappa.outcgs);
        eos_gamma[i][j]    = eos_gamma_;
        eos_gamma1[i][j]    = eos_gamma1_;

        if (l < ntemp) {
          eos_temp[l] = log10(eos_temp_/(simunits->temp.outscale * simunits->temp.outcgs));
        }

        ++l;
        ++j;

        if (!(l % ntemp)) {
          eos_dens[i] = log10(eos_dens_/(simunits->rho.outscale * simunits-> rho.outcgs));
          ++i;
          j = 0;
        }

      }
      else {
        ExceptionHandler::getIstance().raise("Error: The specified opacity file cannot be read!");
      }

    }
    file.close();
  }
}



//=================================================================================================
//  OpacityTable::~OpacityTable()
/// OpacityTable class destructor
//=================================================================================================
template <int ndim>
OpacityTable<ndim>::~OpacityTable()
{
  for (int i = 0; i < ndens; i++) {
    delete[] eos_energy[i];
    delete[] eos_mu[i];
    delete[] kappa_table[i];
    delete[] kappar_table[i];
    delete[] kappap_table[i];
    delete[] eos_gamma[i];
    delete[] eos_gamma1[i];
  }

  delete[] eos_energy;
  delete[] eos_mu;
  delete[] kappa_table;
  delete[] kappar_table;
  delete[] kappap_table;
  delete[] eos_gamma;
  delete[] eos_gamma1;
}


// ------------------------------------------------------------------------------------------//
//			find  closest index in list for value: level
// ------------------------------------------------------------------------------------------//
template <typename BidirectionalIterator, typename T>
BidirectionalIterator getClosest
 (BidirectionalIterator first,
  BidirectionalIterator last,
  const T &value)
{
  BidirectionalIterator before = std::lower_bound(first, last, value);

  if (before == first) return first;
  if (before == last)  return --last;

  BidirectionalIterator after = before;
  --before;

  return (*after - value) < (value - *before) ? after : before;
}


template <typename BidirectionalIterator, typename T>
std::size_t getClosestIndex(BidirectionalIterator first,
                            BidirectionalIterator last,
                            const T &value)
{
  return std::distance(first, getClosest(first, last, value));
}



//=================================================================================================
//  OpacityTable::GetIDens()
/// GetIDens returns table index for density
//=================================================================================================
template <int ndim>
int OpacityTable<ndim>::GetIDens
 (const FLOAT rho)
{
  return getClosestIndex(eos_dens, eos_dens + ndens, rho);
}



//=================================================================================================
//  OpacityTable::GetITemp()
/// GetITemp returns table index for temperature
//=================================================================================================
template <int ndim>
int OpacityTable<ndim>::GetITemp
 (const FLOAT temp)
{
  return getClosestIndex(eos_temp, eos_temp + ntemp , temp);
}



//=================================================================================================
//  OpacityTable::GetIEner()
/// GetIEner returns table index for specific internal energy
//=================================================================================================
template <int ndim>
int OpacityTable<ndim>::GetIEner
 (const FLOAT u, const FLOAT rho)
{
  int idens = GetIDens(log10(rho));
  return getClosestIndex(eos_energy[idens], eos_energy[idens] + ntemp , u);
}
template <int ndim>
int OpacityTable<ndim>::GetIEner
 (const FLOAT u, const int idens)
{
  return getClosestIndex(eos_energy[idens], eos_energy[idens] + ntemp , u);
}


//=================================================================================================
//  OpacityTable::GetKappa()
/// GetKappa returns Kappa  for index of density and temp
//=================================================================================================
template <int ndim>
void OpacityTable<ndim>::GetKappa
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
FLOAT OpacityTable<ndim>::GetEnergy(const int idens, const int itemp)
{
  return eos_energy[idens][itemp];
}



//=================================================================================================
//  OpacityTable::GetMuBar()
/// GetMuBar returns MuBar for index of density and temp
//=================================================================================================
template <int ndim>
FLOAT OpacityTable<ndim>::GetMuBar
 (Particle<ndim> &part)
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
FLOAT OpacityTable<ndim>::GetGamma
 (Particle<ndim> &part)
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
FLOAT OpacityTable<ndim>::GetGamma1
 (Particle<ndim> &part)
{
  int idens = GetIDens(part.rho);
  int iener = GetIEner(part.u, idens);
  return eos_gamma1[idens][iener];
}

template class OpacityTable<1>;
template class OpacityTable<2>;
template class OpacityTable<3>;
