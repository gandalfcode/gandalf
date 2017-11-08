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
#include "EOS.h"
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
 (Parameters* simparams,
  SimUnits *simunits)
{
  int i, j, l;
  DOUBLE eos_dens_, eos_temp_, eos_energy_, eos_mu_, kappa_, kappar_, kappap_, eos_gamma_, eos_gamma1_;
  string line;

  // Check for correct opacity units
  if (simunits->kappa.outunit != "cm2_g" && simparams->stringparams["energy_integration"] == "radws") {
    ExceptionHandler::getIstance().raise("Error: Wrong units for opacity, use cm2_g");
  }

  string radws_table = simparams->stringparams["radws_table"];

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



//=================================================================================================
//  OpacityTable::GetEnergyFromPressure()
/// GetEnergyFromPressure computes the internal energy starting from the pressure
//=================================================================================================
template <int ndim>
FLOAT OpacityTable<ndim>::GetEnergyFromPressure(const FLOAT rho, const FLOAT P) {

  int idens = GetIDens(rho);

  FLOAT* p_gamma  = eos_gamma[idens];
  FLOAT* p_energy = eos_energy[idens];

  // Solve for P/rho == u*(gamma-1)

  int l=0, u=ntemp-1;


  if ((p_gamma[u]-1) * rho * p_energy[u] < P) return P / (rho*(p_gamma[u]-1)) ;
  if ((p_gamma[l]-1) * rho * p_energy[l] > P) return P / (rho*(p_gamma[l]-1)) ;


  // Get the value l that gives a equal or smaller internal energy than the desired value:
  while (l + 1 < u) {
    int c = (l + u)/2;

    FLOAT Pi = (p_gamma[c]-1) * rho * p_energy[c];

    if (Pi > P)
      u = c;
    else
      l = c;
  }

  assert(l+1 < ntemp);

  // Check whether the lower or upper bound is closer:
  if (((p_gamma[l+1]-1)*rho*p_energy[l+1] - P) < (P - (p_gamma[l]-1)*rho*p_energy[l]))
      l++ ;


 return P / (rho*(p_gamma[l]-1)) ;

}

template class OpacityTable<1>;
template class OpacityTable<2>;
template class OpacityTable<3>;
