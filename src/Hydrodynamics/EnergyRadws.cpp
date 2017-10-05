//=================================================================================================
//  EnergyRadws.cpp
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
#include "EnergyEquation.h"
#include "EOS.h"
#include "Exception.h"
#include "Hydrodynamics.h"
#include "Parameters.h"
#include "Particle.h"
#include "RadiativeFB.h"
#include "SimUnits.h"
#include "SmoothingKernel.h"
using namespace std;



//=================================================================================================
//  EnergyRadws::EnergyRadws()
/// EnergyRadws class constructor
//=================================================================================================
template <int ndim, template <int> class ParticleType>
EnergyRadws<ndim,ParticleType>::EnergyRadws
 (DOUBLE energy_mult_, string radws_table, FLOAT temp_ambient_, int lombardi_,
  SimUnits *simunits, EOS<ndim> *eos_, RadiativeFB<ndim> *radfb_) :
  EnergyEquation<ndim>(energy_mult_)
{
  // Set lombardi method flag
  lombardi = lombardi_;

  // Set radiative feedback
  radfb = radfb_;

  int i, j, l;
  // int ndens, ntemp; defined in EnergyEquation
  DOUBLE eos_dens_, eos_temp_, eos_energy_, eos_mu_, kappa_, kappar_, kappap_, eos_gamma_;
  DOUBLE num, denom, tempunit, fcol;
  string line;

  eos = eos_;

  // compute rad_const = Stefan-Boltzmann in code units
  num       = pow(simunits->r.outscale*simunits->r.outSI,2)*simunits->t.outscale*simunits->t.outSI;
  denom     = simunits->E.outscale * simunits->E.outSI;
  tempunit  = simunits->temp.outscale * simunits->temp.outSI;
  rad_const = stefboltz*(num*pow(tempunit,4.0))/denom;
  temp_ambient0 = temp_ambient_ / tempunit;
  temp_min = 5.0 / tempunit;

  // Check for correct opacity units
  if (simunits ->kappa.outunit != "cm2_g") {
    ExceptionHandler::getIstance().raise("Error: Wrong units for opacity, use cm2_g");
  }

  ifstream file;
  file.open(radws_table.c_str(), ios::in);

  //-----------------------------------------------------------------------------------------------
  if (file.good()) {
    getline(file, line);
    istringstream istr(line);
    istr >> ndens >> ntemp >> fcol;

    eos_dens     = new FLOAT[ndens];
    eos_temp     = new FLOAT[ntemp];
    eos_energy   = new FLOAT*[ndens];
    eos_mu       = new FLOAT*[ndens];
    eos_gamma    = new FLOAT*[ndens];
    kappa_table  = new FLOAT*[ndens];
    kappar_table = new FLOAT*[ndens];
    kappap_table = new FLOAT*[ndens];

    for (i = 0; i < ndens; ++i){
      eos_energy[i]   = new FLOAT[ntemp];
      eos_mu[i]       = new FLOAT[ntemp];
      eos_gamma[i]    = new FLOAT[ntemp];
      kappa_table[i]  = new FLOAT[ntemp];
      kappar_table[i] = new FLOAT[ntemp];
      kappap_table[i] = new FLOAT[ntemp];
    }

    // read table
    i = 0;
    l = 0;
    j = 0;

    //---------------------------------------------------------------------------------------------
    while (getline(file, line)) {
      istringstream istr(line);

      if (istr >> eos_dens_ >> eos_temp_ >> eos_energy_ >> eos_mu_>> kappa_ >> kappar_ >> kappap_ >> eos_gamma_) {

        eos_energy[i][j]   = eos_energy_/(simunits->u.outscale * simunits-> u.outcgs);
        eos_mu[i][j]       = eos_mu_;
        eos_gamma[i][j]    = eos_gamma_;
        kappa_table[i][j]  = kappa_/(simunits->kappa.outscale * simunits-> kappa.outcgs);
        kappar_table[i][j] = kappar_/(simunits->kappa.outscale * simunits-> kappa.outcgs);
        kappap_table[i][j] = kappap_/(simunits->kappa.outscale * simunits-> kappa.outcgs);

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

      } else cout << "Dateifehler" << endl;

    }
    //---------------------------------------------------------------------------------------------

    file.close();
  }
  //-----------------------------------------------------------------------------------------------

  // Set fcol2
  if (lombardi) {
    fcol2 = fcol * fcol * 4.0 * pi;
  }
  else {
    fcol2 = fcol * fcol;
  }
}



//=================================================================================================
//  EnergyRadws::~EnergyRadws()
/// EnergyRadws class destructor
//=================================================================================================
template <int ndim, template <int> class ParticleType>
EnergyRadws<ndim,ParticleType>::~EnergyRadws()
{
  for (int i = 0; i < ndens; i++) {
    delete[] eos_energy[i];
    delete[] eos_mu[i];
    delete[] kappa_table[i];
    delete[] kappar_table[i];
    delete[] kappap_table[i];
    delete[] eos_gamma[i];
  }

  delete[] eos_energy;
  delete[] eos_mu;
  delete[] kappa_table;
  delete[] kappar_table;
  delete[] kappap_table;
  delete[] eos_gamma;
}



//=================================================================================================
//  EnergyRadws::EnergyIntegration
/// Integrate internal energy to first order from the beginning of the step to
/// the current simulation time, i.e. u(t+dt) = u(t) + dudt(t)*dt .
//=================================================================================================
template <int ndim, template <int> class ParticleType>
void EnergyRadws<ndim,ParticleType>::EnergyIntegration
 (const int n,                         ///< [in] Integer time in block time struct
  const FLOAT t,                       ///< [in] Current simulation time
  const FLOAT timestep,                ///< [in] Base timestep value
  Hydrodynamics<ndim>* hydro)
{
  int i;                               // Particle counter
  FLOAT dt;                            // Timestep since start of step
  ParticleType<ndim>* partdata = hydro->template GetParticleArray<ParticleType>();

  debug2("[EnergyRadws::EnergyIntegration]");
  CodeTiming::BlockTimer timer = timing->StartNewTimer("ENERGY_RADWS_INTEGRATION");


  //-----------------------------------------------------------------------------------------------
#pragma omp parallel for default(none) private(dt,i) shared(partdata, hydro)
  for (i=0; i<hydro->Nhydro; i++) {
    ParticleType<ndim>& part = partdata[i];
    if (part.flags.is_dead()) continue;

    // Compute time since beginning of current step
    dt = t - part.tlast;

    if (part.dt_therm <= small_number) {
      part.u = part.u0;
    }
    else if (dt < (FLOAT) 40.0 * part.dt_therm) {
      part.u = part.u0 * exp(-dt / part.dt_therm)
             + part.ueq * ((FLOAT) 1.0 - exp(-dt / part.dt_therm));
    }
    else if (dt >= (FLOAT) 40.0 * part.dt_therm) {
      part.u = part.ueq;
    }

  }
  //-----------------------------------------------------------------------------------------------

  return;
}



//=================================================================================================
//  EnergyRadws::EndTimestep
/// Record all important thermal quantities at the end of the step for start of the new timestep.
//=================================================================================================
template <int ndim, template <int> class ParticleType>
void EnergyRadws<ndim,ParticleType>::EndTimestep
 (const int n,                         ///< [in] Integer time in block time struct
  const FLOAT t,                       ///< [in] Current simulation time
  const FLOAT timestep,                ///< [in] Base timestep value
  Hydrodynamics<ndim>* hydro)
{
  int dn;                              // Integer time since beginning of step
  int i;                               // Particle counter
  FLOAT temp;                          // Particle temperature
  FLOAT temp_amb = temp_ambient0;      // Ambient particle temperature
  FLOAT col2;                          // RadWS or Lombardi metric

  ParticleType<ndim>* partdata = hydro->template GetParticleArray<ParticleType>();

  debug2("[EnergyRadws::EndTimestep]");
  CodeTiming::BlockTimer timer = timing->StartNewTimer("ENERGY_RADWS_END_TIMESTEP");


  //-----------------------------------------------------------------------------------------------
#pragma omp parallel for default(none) private(dn, i, temp, col2) shared(partdata, hydro, temp_amb)
  for (i=0; i<hydro->Nhydro; i++) {
    ParticleType<ndim> &part = partdata[i];
    if (part.flags.is_dead()) continue;
    dn = n - part.nlast;

    if (part.flags.check(end_timestep)) {

      temp = GetTemp(part);
      if (radfb) temp_amb = radfb->AmbientTemp(part);
      col2 = GetCol2(part);

      EnergyFindEqui(part.rho, temp, temp_amb, col2, part.u, part.dudt,
                     part.ueq, part.dt_therm);


      part.u0 = part.u;
      part.dudt0 = part.dudt;
    }
  }
  //-----------------------------------------------------------------------------------------------

  return;
}



//=================================================================================================
//  EnergyRadws::~EnergyFindEqui()
/// Computes the thermal equilibrium state of the particle (i.e. its equilibrium temperature and
/// internal energy) including the thermal timescale to reach this equilibrium.
//=================================================================================================
template <int ndim, template <int> class ParticleType>
void EnergyRadws<ndim,ParticleType>:: EnergyFindEqui
 (const FLOAT rho,                             ///< [in] ..
  const FLOAT temp,                            ///< [in] ..
  const FLOAT temp_ambient,                    ///< [in] ..
  const FLOAT col2,                            ///< [in] ..
  const FLOAT u,                               ///< [in] ..
  const FLOAT dudt,                            ///< [in] ..
  FLOAT &ueq_p,                                ///< [out] ..
  FLOAT &dt_thermal)                           ///< [out] ..
{
  const FLOAT logrho   = log10(rho);           // ..
  const int idens      = GetIDens(logrho);     // ..
  int itemp;                                   // ..
  FLOAT logtemp;                               // ..
  FLOAT Tequi;                                 // ..
  FLOAT dudt_eq;                               // ..
  FLOAT dudt_rad;                              // ..
  FLOAT dudt_tot;                              // ..
  FLOAT kappa;                                 // ..
  FLOAT kappar;                                // ..
  FLOAT kappap;                                // ..

  logtemp = log10(temp);
  itemp   = GetITemp(logtemp);

  assert(idens >= 0 && idens <= ndens - 1);
  assert(itemp >= 0 && itemp <= ntemp - 1);

  GetKappa(idens, itemp, logrho, logtemp, kappa, kappar, kappap);
  dudt_rad = ebalance((FLOAT) 0.0, temp_ambient, temp, kappa, kappap, col2);

  // Calculate equilibrium temperature using implicit scheme described in Stamatellos et al. (2007)
  EnergyFindEquiTemp(idens, rho, temp, temp_ambient, col2, dudt, Tequi);
  assert(Tequi >= temp_min);

  // Get ueq_p and dudt_eq from Tequi
  logtemp = log10(Tequi);
  itemp = GetITemp(logtemp);
  GetKappa(idens, itemp, logrho, logtemp, kappa, kappar, kappap);
  ueq_p = GetEnergy(idens, itemp, logrho, logtemp);
  dudt_eq = ebalance((FLOAT) 0.0, temp_ambient, Tequi, kappa, kappap, col2);

  // Thermalization time scale
  dt_thermal = (ueq_p - u) / (dudt + dudt_rad);
  if (dudt_eq == dudt_rad) dt_thermal = 1E30;

  return;
}



//=================================================================================================
//  EnergyRadws::EnergyFindEquiTemp
/// EnergyFindEquiTemp returns equilibrium temperature
//=================================================================================================
template <int ndim, template <int> class ParticleType>
void EnergyRadws<ndim,ParticleType>::EnergyFindEquiTemp
 (const int idens,                         ///< ..
  const FLOAT rho,                         ///< ..
  const FLOAT temp,                        ///< ..
  const FLOAT temp_ambient,                ///< ..
  const FLOAT col2,                        ///< ..
  const FLOAT dudt,                        ///< ..
  FLOAT &Tequi)
{
  int itemp;
  int itemplow, itemphigh;

  FLOAT accuracy = 0.001; // not shure what is a good accuracy (Paul)
  FLOAT balance, minbalance;
  FLOAT balanceLow, balanceHigh ;
  FLOAT kappa, kappar, kappap;
  FLOAT minkappa, minkappar, minkappap;
  FLOAT kappaLow, kappaHigh, kappapLow, kappapHigh;
  FLOAT logtemp, logrho;
  FLOAT Tlow, Thigh, Tlow_log, Thigh_log, Tequi_log;
  FLOAT dtemp;

  logrho = log10(rho);

  logtemp = log10(temp);
  itemp = GetITemp(logtemp);

  // Rosseland and Planck opacities at rho_p and temperature temp(p)
  GetKappa(idens, itemp, logrho, logtemp, kappa, kappar, kappap);
  balance = ebalance((FLOAT) 0.0, temp_ambient, temp, kappa, kappap, col2);

  // Get min. kappa and min balance
  logtemp = log10(temp_min);
  itemp = GetITemp(logtemp);
  GetKappa(idens, itemp, logrho, logtemp, minkappa, minkappar, minkappap);
  minbalance = ebalance((FLOAT) 0.0, temp_ambient, temp_min, minkappa, minkappap, col2);
  // Find equilibrium temperature
  // ------------------------------------------------------------------------------------------------
  if (dudt <= -balance) {
    if (dudt <= -minbalance) {
      Tequi = temp_min;
      return;
    }

    if (itemp <= 1) {
      cout << "Reached itemp = 1, returning" << endl;
      Tequi = pow(10.0, eos_temp[1]);
      return;
    }

    Tlow_log  = eos_temp[itemp - 1];
    Tlow      = pow(10.0, Tlow_log);
    Thigh_log = eos_temp[itemp + 1];
    Thigh     = pow(10.0, Thigh_log);
    itemplow  = itemp - 1;

    GetKappa(idens, itemplow, logrho, Tlow_log, kappaLow, kappar, kappapLow);
    balanceLow = ebalance(dudt, temp_ambient, Tlow, kappaLow, kappapLow, col2);

    itemphigh = itemp;
    GetKappa(idens, itemphigh, logrho, Thigh_log, kappaHigh, kappar, kappapHigh);
    balanceHigh = ebalance(dudt, temp_ambient ,Thigh, kappaHigh, kappapHigh, col2);

    while (balanceLow * balanceHigh > 0.0){
      if (itemp <= 1) {
        Tequi = pow(10.0, eos_temp[1]);
        return;
      }

      Thigh       = Tlow;
      kappaHigh   = kappaLow;
      kappapHigh  = kappapLow;
      balanceHigh = balanceLow;
      itemplow    = itemplow - 1;
      Tlow_log    = eos_temp[itemplow];
      Tlow        = pow(10.0, Tlow_log);

      GetKappa(idens, itemplow, logrho, Tlow_log, kappaLow, kappar, kappapLow);
      balanceLow = ebalance(dudt, temp_ambient , Tlow, kappaLow, kappapLow, col2);
    }

    itemphigh = itemplow + 1;
  } else {
    if (itemp >= ntemp - 2) {
      Tequi = pow(10.0, eos_temp[ntemp - 2]);
      return;
    }

    Tlow_log  = eos_temp[itemp - 1];
    Tlow      = pow(10.0, Tlow_log);
    Thigh_log = eos_temp[itemp + 1];
    Thigh     = pow(10.0, Thigh_log);

    itemplow = itemp;
    GetKappa(idens, itemplow, logrho, Tlow_log, kappaLow, kappar, kappapLow);
    balanceLow = ebalance(dudt, temp_ambient, Tlow, kappaLow, kappapLow, col2);

    itemphigh = itemp + 1;
    GetKappa(idens, itemphigh, logrho, Thigh_log, kappaHigh, kappar, kappapHigh);
    balanceHigh = ebalance(dudt, temp_ambient, Thigh, kappaHigh, kappapHigh, col2);

    while (balanceLow * balanceHigh > 0.0) {
      if (itemphigh >= ntemp - 2) {
        Tequi = pow(10.0, eos_temp[ntemp - 2]);
        return;
      }

      Tlow       = Thigh;
      kappaLow   = kappaHigh;
      kappapLow  = kappapHigh;
      balanceLow = balanceHigh;
      itemphigh  = itemphigh + 1;
      Thigh_log  = eos_temp[itemphigh];
      Thigh      = pow(10.0, Thigh_log);

      GetKappa(idens, itemphigh, logrho, Thigh_log, kappaHigh, kappar, kappapHigh);
      balanceHigh = ebalance(dudt, temp_ambient, Thigh, kappaHigh, kappapHigh, col2);
    }

    itemplow = itemphigh - 1;
  }

  Tequi = 0.5 * (Tlow + Thigh);
  dtemp = Thigh - Tlow;

  // Refine the search in between Thigh and Tlow
  //-----------------------------------------------------------------------------------------------
  while (dtemp != 0.0 && fabs(2.0 * dtemp / (Thigh + Tlow)) > accuracy) {
    Tequi_log = log10(Tequi);
    GetKappa(idens, itemplow, logrho, Tlow_log, kappaLow, kappar, kappapLow);
    balance = ebalance(dudt, temp_ambient , Tequi, kappa, kappap, col2);

    if (balance == 0.0) {
      Tequi = 0.5 * (Thigh + Tlow);
      return;
    }

    if (balanceLow * balance < 0.0){
      Thigh       = Tequi;
      balanceHigh = balance;
      kappaHigh   = kappa;
      kappapHigh  = kappap;
    }
    else {
      Tlow       = Tequi;
      balanceLow = balance;
      kappaLow   = kappa;
      kappapLow  = kappap;
    }

    Tequi = 0.5 * (Thigh + Tlow);
    kappa = (kappaHigh + kappaLow) / 2;
    kappap = (kappapHigh + kappapLow) / 2;
    dtemp = Thigh - Tlow;
  }

  if (Tequi < temp_ambient) Tequi = temp_ambient;

  return;
}



// ------------------------------------------------------------------------------------------//
//			find  closest index in list for value: level
// ------------------------------------------------------------------------------------------//

template <typename BidirectionalIterator, typename T>
BidirectionalIterator getClosest(BidirectionalIterator first,
                                 BidirectionalIterator last,
                                 const T &value)
{
  BidirectionalIterator before = std::lower_bound(first, last, value);

  if (before == first) return first;
  if (before == last)  return --last; // iterator must be bidirectional

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
//  EnergyRadws::GetIDens()
/// GetIDens returns table index for density
//=================================================================================================
template <int ndim, template <int> class ParticleType>
int EnergyRadws<ndim,ParticleType>::GetIDens
 (const FLOAT rho)
{
  int idens = getClosestIndex(eos_dens, eos_dens + ndens, rho);

  if (rho < eos_dens[idens]) {
    idens = max(idens - 1, 0);
  }
  return idens;
}



//=================================================================================================
//  EnergyRadws::GetITemp()
/// GetITemp returns table index for temperature
//=================================================================================================
template <int ndim, template <int> class ParticleType>
int EnergyRadws<ndim,ParticleType>::GetITemp
 (const FLOAT temp)
{
  int itemp = getClosestIndex(eos_temp, eos_temp + ntemp , temp);

  if (temp <= eos_temp[itemp]) {
    itemp = max(itemp - 1, 0);
  }

  return itemp;
}



//=================================================================================================
//  EnergyRadws::GetKappa()
/// GetKappa returns Kappa  for index of density and temp
//=================================================================================================
template <int ndim, template <int> class ParticleType>
void EnergyRadws<ndim,ParticleType>::GetKappa
 (const int idens,
  const int itemp,
  const FLOAT logrho,
  const FLOAT logtemp,
  FLOAT &kappa,
  FLOAT &kappar,
  FLOAT &kappap)
{
  const FLOAT epsilon = (logtemp - eos_temp[itemp])/(eos_temp[itemp+1]-eos_temp[itemp]);
  const FLOAT delta = (logrho - eos_dens[idens])/(eos_dens[idens+1]-eos_dens[idens]);

  kappa = kappa_table[idens+1][itemp+1]*epsilon*delta +
    kappa_table[idens+1][itemp]*(1-epsilon)*delta +
    kappa_table[idens][itemp+1]*epsilon*(1-delta) +
    kappa_table[idens][itemp]*(1-epsilon)*(1-delta);

  kappap = kappap_table[idens+1][itemp+1]*epsilon*delta +
    kappap_table[idens+1][itemp]*(1-epsilon)*delta +
    kappap_table[idens][itemp+1]*epsilon*(1-delta) +
    kappap_table[idens][itemp]*(1-epsilon)*(1-delta);

  kappar = kappap;

  return;
}



//=================================================================================================
//  EnergyRadws::eBalance()
//  Calculates net heating rate due to hydro (i.e. expansion/contraction)
//  and radiative effects (i.e heating/cooling)
//=================================================================================================
template <int ndim, template <int> class ParticleType>
FLOAT EnergyRadws<ndim,ParticleType>::ebalance
 (const FLOAT dudt,
  const FLOAT temp_ex,
  const FLOAT temp,
  const FLOAT kappa,
  const FLOAT kappap,
  const FLOAT col2)
{
  return dudt - 4.0*rad_const*(pow(temp,4) - pow(temp_ex,4))/((col2*kappa) + (1.0/kappap));
}


//=================================================================================================
//  EnergyRadws::GetEnergy()
/// GetEnergy returns Energy  for index of density and temp
//=================================================================================================

template <int ndim, template <int> class ParticleType>
FLOAT EnergyRadws<ndim,ParticleType>:: GetEnergy(int idens, int itemp, FLOAT logrho, FLOAT logtemp)
{
  const FLOAT epsilon = (logtemp - eos_temp[itemp])/(eos_temp[itemp+1] - eos_temp[itemp]);
  const FLOAT delta = (logrho - eos_dens[idens])/(eos_dens[idens+1] - eos_dens[idens]);

  return eos_energy[idens+1][itemp+1]*epsilon*delta +
    eos_energy[idens+1][itemp]*(1 - epsilon)*delta +
    eos_energy[idens][itemp+1]*epsilon*(1 - delta) +
    eos_energy[idens][itemp]*(1 - epsilon)*(1 - delta);
 }



//=================================================================================================
//  EnergyRadws::GetMuBar()
/// GetMuBar returns MuBar  for index of density and temp
//=================================================================================================
template <int ndim, template <int> class ParticleType>
FLOAT EnergyRadws<ndim,ParticleType>::GetMuBar
 (const int idens,
  const int itemp,
  const FLOAT logrho,
  const FLOAT logtemp)
{
  const FLOAT epsilon = (logtemp - eos_temp[itemp])/(eos_temp[itemp+1] - eos_temp[itemp]);
  const FLOAT delta = (logrho - eos_dens[idens])/(eos_dens[idens+1] - eos_dens[idens]);

  return eos_mu[idens+1][itemp+1]*epsilon*delta +
    eos_mu[idens+1][itemp]*(1 - epsilon)*delta +
    eos_mu[idens][itemp+1]*epsilon*(1 - delta) +
    eos_mu[idens][itemp]*(1 - epsilon)*(1 - delta);
}



//=================================================================================================
//  EnergyRadws::GetTemp()
/// GetTemp returns the temperature of particle form it's density and specific internal energy.
/// The value is taken from the EoS table, not via a direct calculation.
//=================================================================================================
template <int ndim, template <int> class ParticleType>
FLOAT EnergyRadws<ndim,ParticleType>::GetTemp
 (Particle<ndim> &part)
{
  FLOAT result = temp_min;

  FLOAT logrho = log10(part.rho);
  int idens = GetIDens(logrho);

  if (part.u == 0.0) return result;

  int temp_index = -1;
  FLOAT du = 1E20;
  for (int i = 0; i < ntemp; ++i) {
    if (abs(part.u - eos_energy[idens][i]) < du) {
      du = abs(part.u - eos_energy[idens][i]);
      temp_index = i;
    }
  }

  result = pow(10.0, eos_temp[temp_index]);

  // Limit temperature to 5K (CMB)
  if (result < temp_min) {
    result = temp_min;
    temp_index = GetITemp(temp_min);
  }

  part.mu_bar = eos_mu[idens][temp_index];
  part.gamma = eos_gamma[idens][temp_index];

  return result;
}

//=================================================================================================
//  EnergyRadws::GetCol2()
/// GetTemp returns the square of the mass-weighted column density average of a particle
/// pseudocloud. The standard RadWS (Stamatellos et al. 2007) method uses the gravitational
/// potential and the density of the particle. The Lombardi et al. (2015) method instead
/// uses the hydrodynamical acceleration of a particle as well as its pressure.
//=================================================================================================
template <int ndim, template <int> class ParticleType>
FLOAT EnergyRadws<ndim,ParticleType>::GetCol2
 (Particle<ndim> &part)
{
  if (!lombardi) {
    return fcol2 * part.gpot_hydro * part.rho;
  }
  else {
    FLOAT mag_da = 0.0, P = 0.0;

    FLOAT tot_a = 0.0, tot_atree = 0.0;
    for (int k = 0; k < ndim; ++k) {
      tot_a += part.a[k] * part.a[k];
      tot_atree += part.atree[k] * part.atree[k];
    }
    FLOAT mag_a = sqrt(tot_a);
    FLOAT mag_atree = sqrt(tot_atree);
    mag_da = mag_a - mag_atree;
    P = (part.gamma - 1.0) * part.rho * part.u;

    return (fcol2 * P * P) / (mag_da * mag_da);
  }
}

template class EnergyRadws<1, GradhSphParticle>;
template class EnergyRadws<2, GradhSphParticle>;
template class EnergyRadws<3, GradhSphParticle>;
template class EnergyRadws<1, SM2012SphParticle>;
template class EnergyRadws<2, SM2012SphParticle>;
template class EnergyRadws<3, SM2012SphParticle>;
