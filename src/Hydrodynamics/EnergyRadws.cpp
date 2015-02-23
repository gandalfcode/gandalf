//=============================================================================
//  EnergyRadws.cpp
//  ...
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
//=============================================================================


#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <math.h>
#include "Sph.h"
#include "SmoothingKernel.h"
#include "SphIntegration.h"
#include "Particle.h"
#include "EOS.h"
#include "EnergyEquation.h"
#include "Exception.h"
#include "Debug.h"
#include "SimUnits.h"
#include "Constants.h"
using namespace std;




//=============================================================================
//  EnergyRadws::EnergyRadws()
/// EnergyRadws class constructor
//=============================================================================
template <int ndim, template <int> class ParticleType>
EnergyRadws<ndim,ParticleType>::EnergyRadws(DOUBLE energy_mult_aux, string radws_table,
                                            FLOAT temp_ambient_, SimUnits *simunits) :
  EnergyEquation<ndim>(energy_mult_aux)
{

  int i, j, l;
    // int ndens, ntemp; defined in EnergyEquation
  double eos_dens_, eos_temp_, eos_energy_, eos_mu_, kappa_, kappar_, kappap_;
  string line;
  float num, denom, tunit;

  // compute rad_const = Stefan-Boltzmann in code units
  num = simunits ->r.outscale * simunits ->r.outSI;
  denom = simunits -> E.outscale * simunits ->E.outSI;
  tunit = simunits -> t.outscale * simunits -> t.outSI;
  rad_const = stefboltz*num*num*tunit/denom ;

  temp_ambient = temp_ambient_ / (simunits->temp.outscale * simunits-> temp.outSI);

  // check that user wants cm2_g for opacity units
  if (simunits ->kappa.outunit != "cm2_g") {
    cout << "ERROR! Selected wrong unit for opacity. Use cm2_g";
    exit(0);
  }


  ifstream file;
  file.open(radws_table.c_str(),ios::in) ;
  // allocate table entries


  cout << "EnergyRadws:" << radws_table << "  "<< file.good() << endl;

  if (file.good()){
    getline(file, line);
    cout << line << endl;
    istringstream istr(line);
    istr >> ndens >> ntemp ;



    cout << "ndens, ntemp= " << ndens << "  "<< ntemp << endl;

    eos_dens = new FLOAT[ndens];
    eos_temp = new FLOAT[ntemp];

    eos_energy = new FLOAT*[ndens];
    eos_mu = new FLOAT*[ndens];
    kappa_table = new FLOAT*[ndens];
    kappar_table = new FLOAT*[ndens];
    kappap_table = new FLOAT*[ndens];

    for (i=0; i<ndens; i++){

      eos_energy[i] = new FLOAT[ntemp];
      eos_mu[i] = new FLOAT[ntemp];
      kappa_table[i] = new FLOAT[ntemp];
      kappar_table[i] = new FLOAT[ntemp];
      kappap_table[i] = new FLOAT[ntemp];

    }

    // read table
    i=0;
    l=0;
    j=0;

    while (getline(file, line)) {
      istringstream istr(line);


      if (istr >> eos_dens_>> eos_temp_>> eos_energy_>> eos_mu_>> kappa_>>
	  kappar_>> kappap_) {


	eos_energy[i][j] = eos_energy_/(simunits->E.outscale * simunits-> E.outcgs);
	eos_mu[i][j] = eos_mu_;
	kappa_table[i][j] = kappa_/(simunits->kappa.outscale * simunits-> kappa.outcgs);
	kappar_table[i][j] = kappar_/(simunits->kappa.outscale * simunits-> kappa.outcgs);
	kappap_table[i][j] = kappap_/(simunits->kappa.outscale * simunits-> kappa.outcgs);

      if (l < ntemp) {
	eos_temp[l] = log10(eos_temp_/(simunits->temp.outscale * simunits-> temp.outcgs));
      };

      l += 1;
      j += 1;

      if (l % ntemp == 0) {
	eos_dens[i] = log10(eos_dens_/(simunits->rho.outscale * simunits-> rho.outcgs));
	i += 1;
	j=0;
      };



    } else cout << "Dateifehler" << endl;
  }

  file.close();
  }
  //cout << getPosition_dens(test_dens) << endl;
  //cout << getPosition_temp(test_temp) << endl;
  //cout << kappap[getPosition_dens(test_dens)][getPosition_temp(test_temp)] << endl;

  //////////////////// SAVE eos_dens and eos_temp in log10 ////////
  //  eos_dens = log10(eos_dens);
  // eos_temp = log10(eos_temp);

  cout << "eos_dens=" <<eos_dens[0] << "  "<<eos_dens[ndens-1] <<endl ;
  cout << "eos_temp="<< eos_temp[0]<< "  "<<eos_temp[ntemp-1] <<endl;

}



//=============================================================================
//  EnergyRadws::~EnergyRadws()
/// EnergyRadws class destructor
//=============================================================================
template <int ndim, template <int> class ParticleType>
EnergyRadws<ndim,ParticleType>::~EnergyRadws()
{

  int i;
  //deallocate arrays

  for (i=0; i<ndens; i++){

    delete[] eos_energy[i];
    delete[] eos_mu[i];
    delete[] kappa_table[i];
    delete[] kappar_table[i];
    delete[] kappap_table[i];

  }

  delete[] eos_energy;
  delete[] eos_mu;
  delete[] kappa_table;
  delete[] kappar_table;
  delete[] kappap_table;


}



//=================================================================================================
//  EnergyRadws::EnergyIntegration
/// Integrate internal energy to first order from the beginning of the step to
/// the current simulation time, i.e. u(t+dt) = u(t) + dudt(t)*dt
//=================================================================================================
template <int ndim, template <int> class ParticleType>
void EnergyRadws<ndim,ParticleType>::EnergyIntegration
 (const int n,                         ///< [in] Integer time in block time struct
  const int Npart,                     ///< [in] Number of particles
  const FLOAT t,                       ///< [in] Current simulation time
  const FLOAT timestep,                ///< [in] Base timestep value
  SphParticle<ndim>* sph_gen)          ///< [inout] Pointer to SPH particle array
{
  int i;                               // Particle counter
  FLOAT dt;                            // Timestep since start of step
  ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sph_gen);

  debug2("[EnergyRadws::EnergyIntegration]");
  timing->StartTimingSection("ENERGY_RADWS_INTEGRATION");

  //-----------------------------------------------------------------------------------------------
#pragma omp parallel for default(none) private(dt,i) shared(sphdata)
  for (i=0; i<Npart; i++) {
    SphParticle<ndim>& part = sphdata[i];
    if (part.itype == dead) continue;

    // Compute time since beginning of current step
    dt = t - part.tlast;
    part.u = part.u0 + part.dudt0*dt;
  }
  //-----------------------------------------------------------------------------------------------

  timing->EndTimingSection("ENERGY_RADWS_INTEGRATION");

  return;
}



//=================================================================================================
//  EnergyRadws::EndTimestep
/// Record all important thermal quantities at the end of the step for start of the new timestep.
//=================================================================================================
template <int ndim, template <int> class ParticleType>
void EnergyRadws<ndim,ParticleType>::EndTimestep
 (const int n,                         ///< [in] Integer time in block time struct
  const int Npart,                     ///< [in] Number of particles
  const FLOAT t,                       ///< [in] Current simulation time
  const FLOAT timestep,                ///< [in] Base timestep value
  SphParticle<ndim>* sph_gen)          ///< [inout] Pointer to SPH particle array
{
  int dn;                           // Integer time since beginning of step
  int i;                            // Particle counter
  ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sph_gen);

  debug2("[EnergyRadws::EndTimestep]");
  timing->StartTimingSection("ENERGY_RADWS_END_TIMESTEP");

  //-----------------------------------------------------------------------------------------------
#pragma omp parallel for default(none) private(dn,i) shared(sphdata)
  for (i=0; i<Npart; i++) {
    SphParticle<ndim>& part = sphdata[i];
    if (part.itype == dead) continue;
    dn = n - part.nlast;

    if (dn == part.nstep) {
      part.u     += 0.5*(part.dudt - part.dudt0)*(t - part.tlast); //timestep*(FLOAT) nstep;
      part.u0    = part.u;
      part.dudt0 = part.dudt;
    }

  }
  //-----------------------------------------------------------------------------------------------

  timing->EndTimingSection("ENERGY_RADWS_END_TIMESTEP");

  return;
}



//=============================================================================
//  EnergyRadws::~EnergyFindEqui()
/// EnergyFindEqui returns equilibrium temperature
//=============================================================================
template <int ndim, template <int> class ParticleType>
void EnergyRadws<ndim,ParticleType>:: EnergyFindEqui(FLOAT rho, FLOAT temp, FLOAT gpot, FLOAT u, FLOAT dudt, FLOAT &ueq_p, FLOAT &dt_thermal)
{
  FLOAT fcolumn2 = 0.010816;

  int idens, itemp;

  FLOAT col2;
  FLOAT Tequi;
  FLOAT logtemp, logrho;
  FLOAT dudt_eq, dudt_rad, dudt_tot ;
  FLOAT kappa, kappar, kappap;

  col2 = fcolumn2 * gpot * rho;

  logrho = log10(rho);

  idens = GetIDens(rho);
  itemp = GetITemp(temp);
  logtemp = log10(temp);

  GetKappa(idens,  itemp,  logrho, logtemp , kappa,  kappar, kappap);
  dudt_rad = ebalance((FLOAT)0.0, temp_ambient , temp, kappa, kappap, col2);

  //////////// GET TEQUI

  EnergyFindEquiTemp(idens,  rho,  temp, col2, dudt, Tequi);

  itemp = GetITemp(Tequi);
  logtemp = log10(Tequi);

  GetKappa(idens,  itemp,  logrho, logtemp , kappa,  kappar, kappap);
  ueq_p = GetEnergy(idens, itemp, logrho, logtemp);

  dudt_eq = ebalance((FLOAT)0.0, temp_ambient , Tequi, kappa, kappap, col2);

  //// GET CHANGE IN DUDT

  dudt_tot = -(dudt_eq - dudt_rad);

  /// Thermalization time scale

  if (dudt_tot == 0.0){
    dt_thermal = 1.e20;
  }
  else {
    dt_thermal = (ueq_p - u)/(dudt + dudt_rad);
  }

  ////DONE!!!
  return;

}
//=============================================================================
//  EnergyRadws::~EnergyFindEqui()
/// EnergyFindEquiTemp returns equilibrium temperature
//=============================================================================
template <int ndim, template <int> class ParticleType>
void EnergyRadws<ndim,ParticleType>:: EnergyFindEquiTemp(int idens, FLOAT rho, FLOAT temp, FLOAT col2, FLOAT dudt, FLOAT &Tequi)
{
  int itemp;
  int itemplow, itemphigh;

  FLOAT accuracy = 0.001;
  FLOAT balance, minbalance ;
  FLOAT balanceLow, balanceHigh ;
  FLOAT kappa, kappar, kappap;
  FLOAT kappaLow, kappaHigh, kappapLow, kappapHigh;
  FLOAT logtemp, logrho , tempmin;
  FLOAT Tlow, Thigh, Tlow_log, Thigh_log, Tequi_log;
  FLOAT dtemp;

  itemp = GetITemp(temp);

  // Rosseland and Planck opacities at rho_p and temperature temp(p)
  logrho = log10(rho);
  logtemp = log10(temp);

  GetKappa(idens,  itemp,  logrho,  logtemp, kappa,  kappar, kappap);

  // calculate radiative heating/cooling rate
  balance = ebalance((FLOAT)0.0, temp_ambient , temp, kappa, kappap, col2);

// Find equilibrium temperature
// ----------------------------------------------------------------------------

  // Rosseland and Planck opacities at rho_p and mintemperature 5K
  tempmin = max(5.0,temp_ambient);
  logtemp = log10(tempmin);
  GetKappa(idens,  itemp,  logrho,  logtemp, kappa,  kappar, kappap);
  // calculate minradiative heating/cooling rate
  minbalance = ebalance((FLOAT)0.0, temp_ambient , tempmin, kappa, kappap, col2);

  if (dudt <= -minbalance){
    Tequi = tempmin;
    return;
  }
  else if (dudt <= -balance){
    if (itemp <= 1) {
      Tequi = pow(10.0, eos_temp[1]);
      return;
    }
    Tlow_log = eos_temp[itemp-1];
    Tlow = pow(10.0, Tlow_log);
    Thigh_log = eos_temp[itemp+1];
    Thigh = pow(10.0, Thigh_log);

    itemplow= itemp-1;
    GetKappa(idens,  itemplow,  logrho,  Tlow_log, kappaLow,  kappar, kappapLow);
    balanceLow = ebalance(dudt, temp_ambient , Tlow, kappaLow, kappapLow, col2);

    itemphigh = itemp +1;
    GetKappa(idens,  itemphigh,  logrho,  Thigh_log, kappaHigh,  kappar, kappapHigh);
    balanceHigh = ebalance(dudt, temp_ambient , Thigh, kappaHigh, kappapHigh, col2);


    while (balanceLow*balanceHigh > 0.0){
      if (itemp <= 1) {
	Tequi = pow(10.0, eos_temp[1]);
	return;
      }

      Thigh= Tlow ;
      kappaHigh   = kappaLow;
      kappapHigh  = kappapLow;
      balanceHigh = balanceLow;
      itemplow = itemplow -1;
      Tlow_log = eos_temp[itemplow];
      Tlow = pow(10.0, Tlow_log);

      GetKappa(idens,  itemplow,  logrho,  Tlow_log, kappaLow,  kappar, kappapLow);
      balanceLow = ebalance(dudt, temp_ambient , Tlow, kappaLow, kappapLow, col2);

    }

    dtemp = Thigh-Tlow;
  }
  else {
    // dudt > -balance; search higher temperatures

    if (itemp >= ntemp-2) {
      Tequi = pow(10.0, eos_temp[ntemp-2]);
      return;
    }
    Tlow_log = eos_temp[itemp-1];
    Tlow = pow(10.0, Tlow_log);
    Thigh_log = eos_temp[itemp+1];
    Thigh = pow(10.0, Thigh_log);

    itemplow= itemp-1;
    GetKappa(idens,  itemplow,  logrho,  Tlow_log, kappaLow,  kappar, kappapLow);
    balanceLow = ebalance(dudt, temp_ambient , Tlow, kappaLow, kappapLow, col2);

    itemphigh = itemp +1;
    GetKappa(idens,  itemphigh , logrho,  Thigh_log, kappaHigh,  kappar, kappapHigh);
    balanceHigh = ebalance(dudt, temp_ambient , Thigh, kappaHigh, kappapHigh, col2);

    while (balanceLow*balanceHigh > 0.0){
      if (itemp >= ntemp-2) {
	Tequi = pow(10.0, eos_temp[ntemp-2]);
	return;
      }

      Tlow= Thigh ;
      kappaLow   = kappaHigh;
      kappapLow  = kappapHigh;
      balanceLow = balanceHigh;
      itemphigh = itemphigh +1;
      Thigh_log = eos_temp[itemphigh];
      Thigh = pow(10.0, Thigh_log);

      GetKappa(idens,  itemphigh,  logrho,  Thigh_log, kappaHigh,  kappar, kappapHigh);
      balanceHigh = ebalance(dudt, temp_ambient , Thigh, kappaHigh, kappapHigh, col2);

    }
    itemplow = itemphigh-1;
    dtemp = Thigh-Tlow;

  }
  ///===========================================================================


  Tequi = 0.5*(Tlow + Thigh);
  dtemp = Thigh-Tlow;

  ///===========================================================================
  // Refine the search in between Thigh and Tlow

  while(fabs(2.*dtemp/(Thigh+Tlow)) > accuracy){
    Tequi_log = log10(Tequi);
    GetKappa(idens,  itemplow,  logrho,  Tequi_log, kappa,  kappar, kappap);
    balance = ebalance(dudt, temp_ambient , Tequi, kappa, kappap, col2);

    if (balance == 0.0) {
      // LEAVE ROUTINE
      Tequi = 0.5*(Thigh + Tlow);
      return;
    }

    if (balanceLow*balance < 0.0){
      Thigh = Tequi;
      balanceHigh = balance;
      kappaHigh = kappa;
      kappapHigh = kappap;
    }
    else{
      Tlow = Tequi;
      balanceLow = balance;
      kappaLow = kappa;
      kappapLow = kappap;
    }

    Tequi = 0.5*(Thigh + Tlow);
    dtemp = Thigh-Tlow;

  }



  return;

}
// ------------------------------------------------------------------------------------------//
//			find  closest index in list for value: level
// ------------------------------------------------------------------------------------------//

template <typename BidirectionalIterator, typename T>
BidirectionalIterator getClosest(BidirectionalIterator first,
                                 BidirectionalIterator last,
                                 const T & value)
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
                            const T & value)
{
    return std::distance(first, getClosest(first, last, value));
}



//=============================================================================
//  EnergyRadws::GetIDens()
/// GetIDens returns table index for density
//=============================================================================
template <int ndim, template <int> class ParticleType>
int EnergyRadws<ndim,ParticleType>:: GetIDens(FLOAT rho)
{
  int idens ;
  idens = getClosestIndex(eos_dens, eos_dens + ndens , rho);

  if (rho < eos_dens[idens]) {
    idens = max(idens-1, 0);
  }

  return idens;

}
//=============================================================================
//=============================================================================
//  EnergyRadws::GetITemp()
/// GetITemp returns table index for temperature
//=============================================================================
template <int ndim, template <int> class ParticleType>
int EnergyRadws<ndim,ParticleType>:: GetITemp(FLOAT temp)
{
  int itemp;

  itemp= getClosestIndex(eos_temp, eos_temp + ntemp , temp);

  if (temp < eos_temp[itemp]) {
    itemp = max(itemp-1, 0);
  }

  return itemp;

}

//=============================================================================
//  EnergyRadws::GetKappa()
/// GetKappa returns Kappa  for index of density and temp
//=============================================================================
template <int ndim, template <int> class ParticleType>
  void EnergyRadws<ndim,ParticleType>:: GetKappa(int idens, int itemp, FLOAT logrho, FLOAT logtemp,
				    FLOAT &kappa, FLOAT &kappar, FLOAT &kappap)
{
  FLOAT epsilon;
  FLOAT delta;


  epsilon = (logtemp - eos_temp[itemp])/(eos_temp[itemp+1]-eos_temp[itemp]);
  delta = (logrho - eos_dens[idens])/(eos_dens[idens+1]-eos_dens[idens]);


  kappa = kappa_table[idens+1][itemp+1]*epsilon*delta + kappa_table[idens+1][itemp]*(1-epsilon)*delta +
    kappa_table[idens][itemp+1]*epsilon*(1-delta) +kappa_table[idens][itemp]*(1-epsilon)*(1-delta);

  //  kappar_loc = kappar[idens+1,itemp+1]*epsilon*delta + kappar[idens+1,itemp]*(1-epsilon)*delta +
  //  kappar[idens,itemp+1]*epsilon*(1-delta) +kappar[idens,itemp]*(1-epsilon)*(1-delta);

  kappap = kappap_table[idens+1][itemp+1]*epsilon*delta + kappap_table[idens+1][itemp]*(1-epsilon)*delta +
    kappap_table[idens][itemp+1]*epsilon*(1-delta) +kappap_table[idens][itemp]*(1-epsilon)*(1-delta);

  return;

}

//=============================================================================
//  EnergyRadws::eBalance()
//  Calculates net heating rate due to hydro (i.e. expansion/contraction)
//  and radiative effects (i.e heating/cooling)
//=============================================================================
template <int ndim, template <int> class ParticleType>
FLOAT EnergyRadws<ndim,ParticleType>:: ebalance(FLOAT dudt, FLOAT temp_ex, FLOAT temp,FLOAT kappa, FLOAT kappap,FLOAT col2)
{
  return dudt - 4.0*rad_const*(pow(temp,4) - pow(temp_ex,4))/((col2*kappa) + (1.0/kappap));
}


//=============================================================================
//  EnergyRadws::GetEnergy()
/// GetEnergy returns Energy  for index of density and temp
//=============================================================================
template <int ndim, template <int> class ParticleType>
FLOAT EnergyRadws<ndim,ParticleType>:: GetEnergy(int idens, int itemp, FLOAT logrho, FLOAT logtemp)
{
  FLOAT epsilon;
  FLOAT delta;


  epsilon = (logtemp - eos_temp[itemp])/(eos_temp[itemp+1]-eos_temp[itemp]);
  delta = (logrho - eos_dens[idens])/(eos_dens[idens+1]-eos_dens[idens]);


  return eos_energy[idens+1][itemp+1]*epsilon*delta + eos_energy[idens+1][itemp]*(1-epsilon)*delta +
    eos_energy[idens][itemp+1]*epsilon*(1-delta) +eos_energy[idens][itemp]*(1-epsilon)*(1-delta);

 }
//=============================================================================

/// 1D/2D/3D


template class EnergyRadws<1, GradhSphParticle>;
template class EnergyRadws<2, GradhSphParticle>;
template class EnergyRadws<3, GradhSphParticle>;
template class EnergyRadws<1, SM2012SphParticle>;
template class EnergyRadws<2, SM2012SphParticle>;
template class EnergyRadws<3, SM2012SphParticle>;
