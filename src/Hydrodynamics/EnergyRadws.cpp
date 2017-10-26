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
template <int ndim>
EnergyRadwsBase<ndim>::EnergyRadwsBase
 (Parameters* simparams,
  SimUnits *simunits,
  Radws<ndim> *eos_,
  RadiativeFB<ndim> *radfb_) :
  EnergyEquation<ndim>(simparams->floatparams["energy_mult"])
{
  radfb = radfb_;
  eos = eos_;
  lombardi = simparams->intparams["lombardi_method"];
  table = eos->opacity_table;
  ndens = table->ndens;
  ntemp = table->ntemp;

  // compute rad_const = Stefan-Boltzmann in code units
  FLOAT num = 0.0, denom = 0.0, tempunit = 0.0;
  num       = pow(simunits->r.outscale*simunits->r.outSI,2)*simunits->t.outscale*simunits->t.outSI;
  denom     = simunits->E.outscale * simunits->E.outSI;
  tempunit  = simunits->temp.outscale * simunits->temp.outSI;
  rad_const = stefboltz*(num*pow(tempunit,4.0))/denom;
  temp_ambient0 = simparams->floatparams["temp_ambient"] / tempunit;
  temp_min = 5.0 / tempunit;

  // Set fcol2
  FLOAT fcol = table->fcol;
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
template <int ndim>
EnergyRadwsBase<ndim>::~EnergyRadwsBase()
{

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

      temp = eos->Temperature(part);
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
//  EnergyRadws::EndTimestep
/// Compute the cooling for the old and new time-steps for the meshless
//=================================================================================================
template <int ndim>
void EnergyRadws<ndim,MeshlessFVParticle>::EndTimestep
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

  MeshlessFVParticle<ndim>* partdata = hydro->template GetParticleArray<MeshlessFVParticle>();
  MeshlessFV<ndim>* mfv = reinterpret_cast<MeshlessFV<ndim>*>(hydro);

  debug2("[EnergyRadws::EndTimestep]");
  CodeTiming::BlockTimer timer = timing->StartNewTimer("ENERGY_RADWS_END_TIMESTEP");


  //-----------------------------------------------------------------------------------------------
#pragma omp parallel for default(none) private(dn, i, temp, col2) shared(partdata, mfv, temp_amb)
  for (i=0; i<mfv->Nhydro; i++) {
    MeshlessFVParticle<ndim> &part = partdata[i];
    if (part.flags.is_dead()) continue;

    if (part.flags.check(end_timestep)) {

      // Compute the primitive variables
      FLOAT Qcons[MeshlessFV<ndim>::nvar] ;
      for (int k=0; k<MeshlessFV<ndim>::nvar; k++) Qcons[k] = part.Qcons0[k] + part.dQ[k];

      mfv->UpdateArrayVariables(part, Qcons);
      mfv->ComputeThermalProperties(part);
      mfv->UpdatePrimitiveVector(part);

      // Compute an estimate of dudt:
      FLOAT dudt = part.dQ[MeshlessFV<ndim>::ietot] -
          0.5*(DotProduct(Qcons,Qcons,ndim)/Qcons[MeshlessFV<ndim>::irho] -
              DotProduct(part.Qcons0,part.Qcons0,ndim)/part.Qcons0[MeshlessFV<ndim>::irho]);
      dudt /= part.dt ;

      temp = eos->Temperature(part);
      if (radfb) temp_amb = radfb->AmbientTemp(part);
      col2 = GetCol2(part);

      EnergyFindEqui(part.rho, temp, temp_amb, col2, part.u, dudt,
          part.ueq, part.dt_therm);


      // Work out the cooling rate of the previous and new timestep.
      if (part.dt < 1e-4*part.dt_therm) {
        FLOAT Xi_old = part.dt / part.dt_therm ;
        part.dQ[MeshlessFV<ndim>::ietot] +=
            (part.ueq - part.u) * Xi_old*(1 - Xi_old/2);
      }
      else {
        part.dQ[MeshlessFV<ndim>::ietot] +=
            (part.ueq - part.u) * (1 - exp( - part.dt / part.dt_therm));
      }

      FLOAT Xi_new;
      if (part.dt_next < 1e-4*part.dt_therm) {
        Xi_new = part.dt_next/part.dt_therm;
        Xi_new = (1 - Xi_new/2) / part.dt_therm;
      }
      else {
        Xi_new = (1 - exp( - part.dt_next / part.dt_therm)) / part.dt_next ;
      }

      // Average cooling rate over the next time-step
      part.cooling = - (part.u - part.ueq)*exp(-part.dt/part.dt_therm)*Xi_new ;

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
template <int ndim>
void EnergyRadwsBase<ndim>::EnergyFindEqui
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
  const int idens = table->GetIDens(logrho);
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
  itemp   = table->GetITemp(logtemp);

  assert(idens >= 0 && idens <= ndens - 1);
  assert(itemp >= 0 && itemp <= ntemp - 1);

  table->GetKappa(idens, itemp, kappa, kappar, kappap);
  dudt_rad = ebalance((FLOAT) 0.0, temp_ambient, temp, kappa, kappap, col2);

  // Calculate equilibrium temperature using implicit scheme described in Stamatellos et al. (2007)
  EnergyFindEquiTemp(idens, rho, temp, temp_ambient, col2, dudt, Tequi);
  assert(Tequi >= temp_min);

  // Get ueq_p and dudt_eq from Tequi
  logtemp = log10(Tequi);
  itemp = table->GetITemp(logtemp);
  table->GetKappa(idens, itemp, kappa, kappar, kappap);
  ueq_p = table->GetEnergy(idens, itemp);
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
template <int ndim>
void EnergyRadwsBase<ndim>::EnergyFindEquiTemp
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

  FLOAT accuracy = 0.001;
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
  itemp = table->GetITemp(logtemp);

  // Rosseland and Planck opacities at rho_p and temperature temp(p)
  table->GetKappa(idens, itemp, kappa, kappar, kappap);
  balance = ebalance((FLOAT) 0.0, temp_ambient, temp, kappa, kappap, col2);

  // Get min. kappa and min balance
  logtemp = log10(temp_min);
  itemp = table->GetITemp(logtemp);
  table->GetKappa(idens, itemp, minkappa, minkappar, minkappap);
  minbalance = ebalance((FLOAT) 0.0, temp_ambient, temp_min, minkappa, minkappap, col2);

  // Find equilibrium temperature
  // ------------------------------------------------------------------------------------------------
  if (dudt <= -balance) {
    if (dudt <= -minbalance) {
      Tequi = temp_min;
      return;
    }

    Tlow_log  = table->eos_temp[itemp - 1];
    Tlow      = pow(10.0, Tlow_log);
    Thigh_log = table->eos_temp[itemp + 1];
    Thigh     = pow(10.0, Thigh_log);
    itemplow  = itemp - 1;

    table->GetKappa(idens, itemplow, kappaLow, kappar, kappapLow);
    balanceLow = ebalance(dudt, temp_ambient, Tlow, kappaLow, kappapLow, col2);

    itemphigh = itemp;
    table->GetKappa(idens, itemphigh, kappaHigh, kappar, kappapHigh);
    balanceHigh = ebalance(dudt, temp_ambient ,Thigh, kappaHigh, kappapHigh, col2);

    while (balanceLow * balanceHigh > 0.0){
      assert(itemp >= 1);

      Thigh       = Tlow;
      kappaHigh   = kappaLow;
      kappapHigh  = kappapLow;
      balanceHigh = balanceLow;
      itemplow    = itemplow - 1;
      Tlow_log    = table->eos_temp[itemplow];
      Tlow        = pow(10.0, Tlow_log);

      table->GetKappa(idens, itemplow, kappaLow, kappar, kappapLow);
      balanceLow = ebalance(dudt, temp_ambient , Tlow, kappaLow, kappapLow, col2);
    }

    itemphigh = itemplow + 1;
  } else {
    Tlow_log  = table->eos_temp[itemp - 1];
    Tlow      = pow(10.0, Tlow_log);
    Thigh_log = table->eos_temp[itemp + 1];
    Thigh     = pow(10.0, Thigh_log);

    itemplow = itemp;
    table->GetKappa(idens, itemplow, kappaLow, kappar, kappapLow);
    balanceLow = ebalance(dudt, temp_ambient, Tlow, kappaLow, kappapLow, col2);

    itemphigh = itemp + 1;
    table->GetKappa(idens, itemphigh, kappaHigh, kappar, kappapHigh);
    balanceHigh = ebalance(dudt, temp_ambient, Thigh, kappaHigh, kappapHigh, col2);

    while (balanceLow * balanceHigh > 0.0) {
      assert(itemp <= ntemp - 2);

      Tlow       = Thigh;
      kappaLow   = kappaHigh;
      kappapLow  = kappapHigh;
      balanceLow = balanceHigh;
      itemphigh  = itemphigh + 1;
      Thigh_log  = table->eos_temp[itemphigh];
      Thigh      = pow(10.0, Thigh_log);

      table->GetKappa(idens, itemphigh, kappaHigh, kappar, kappapHigh);
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
    table->GetKappa(idens, itemplow, kappaLow, kappar, kappapLow);
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



//=================================================================================================
//  EnergyRadws::eBalance()
//  Calculates net heating rate due to hydro (i.e. expansion/contraction)
//  and radiative effects (i.e heating/cooling)
//=================================================================================================
template <int ndim>
FLOAT EnergyRadwsBase<ndim>::ebalance
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
//  EnergyRadws::GetCol2()
/// GetTemp returns the square of the mass-weighted column density average of a particle
/// pseudocloud. The standard RadWS (Stamatellos et al. 2007) method uses the gravitational
/// potential and the density of the particle. The Lombardi et al. (2015) method instead
/// uses the hydrodynamical acceleration of a particle as well as its pressure.
//=================================================================================================
template <int ndim>
FLOAT EnergyRadwsBase<ndim>::GetCol2
 (Particle<ndim> &part)
{
  if (!lombardi) {
    assert(part.gpot_hydro <= part.gpot);
    return fcol2 * part.gpot_hydro * part.rho;
  }
  else {
    FLOAT P = 0.0, ahydro = 0.0;
    for (int k = 0; k < ndim; ++k) {
      ahydro += pow(part.a[k] - part.atree[k], 2.0);
    }
    P = eos->Pressure(part);

    return (fcol2 * P * P) / (ahydro + small_number);
  }
}

template class EnergyRadws<1, GradhSphParticle>;
template class EnergyRadws<2, GradhSphParticle>;
template class EnergyRadws<3, GradhSphParticle>;
template class EnergyRadws<1, SM2012SphParticle>;
template class EnergyRadws<2, SM2012SphParticle>;
template class EnergyRadws<3, SM2012SphParticle>;
template class EnergyRadws<1, MeshlessFVParticle>;
template class EnergyRadws<2, MeshlessFVParticle>;
template class EnergyRadws<3, MeshlessFVParticle>;
