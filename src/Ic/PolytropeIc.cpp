//=================================================================================================
//  PolytropeIc.cpp
//  Subroutines used to generate initial conditions with a Polytrope density field.
//  Can either solve for the general of isothermal Lane-Emden equations.
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


#include <fstream>
#include <sstream>
#include "Precision.h"
#include "Debug.h"
#include "Ic.h"
using namespace std;



//=================================================================================================
//  PolytropeIc::PolytropeIc
/// Constructor for PolytropicIc class
//=================================================================================================
template <int ndim>
PolytropeIc<ndim>::PolytropeIc(Simulation<ndim>* _sim, FLOAT _invndim) :
  Ic<ndim>(_sim, _invndim)
{
  int Ntablemax = 1000;
  xiArray    = new FLOAT[Ntablemax];
  psiArray   = new FLOAT[Ntablemax];
  phiArray   = new FLOAT[Ntablemax];
  muArray    = new FLOAT[Ntablemax];
  pressArray = new FLOAT[Ntablemax];
  rhoArray   = new FLOAT[Ntablemax];
  thetaArray = new FLOAT[Ntablemax];
  massArray  = new FLOAT[Ntablemax];
}



//=================================================================================================
//  PolytropeIc::~PolytropeIc
/// Constructor for PolytropicIc class
//=================================================================================================
template <int ndim>
PolytropeIc<ndim>::~PolytropeIc()
{
}



//=================================================================================================
//  PolytropeIc::Generate
/// Create the basic particle distribution before stretching to the correct density profile
//=================================================================================================
template <int ndim>
void PolytropeIc<ndim>::Generate(void)
{
  // Only compile for 3-dimensional case
  //-----------------------------------------------------------------------------------------------
  if (ndim == 3) {

    int i;                               // Particle counter
    int k;                               // Dimension counter
    FLOAT mp;                            // Mass of individual particle
    //FLOAT z;                             // z-position of newly inserted particle

    // Create local copies of initial conditions parameters
    int Npart      = simparams->intparams["Nhydro"];
    FLOAT eta_eos  = simparams->floatparams["eta_eos"];
    string eos     = simparams->stringparams["gas_eos"];

    debug2("[PolytropeIc::Generate]");

    // Check parameters are valid
    if (eta_eos < (FLOAT) 1.0) {
      ExceptionHandler::getIstance().raise("Invalid value for eta_eos (< 1)");
    }

    // Allocate local and main particle memory
    hydro->Nhydro = Npart;
    sim->AllocateParticleMemory();
    mp = (FLOAT) 0.0;
    //mp = m_box / (FLOAT) Npart;

    // Calculate cloud parameters depending on the chosen equation of state
    if (eos == "isothermal" || fabs(eta_eos - 1.0) < (FLOAT) 0.0001) {
      std::cout << "Setting up isothermal polytrope" << std::endl;
    }
    else {
      std::cout << "Setting up general polytrope" << std::endl;
    }


    std::cout << "Nhydro : " << Npart << std::endl;

    // Create box of particles


    // Record particle properties in main memory
    // Stretch box to occupy correct volume specified in parameters file
    for (i=0; i<hydro->Nhydro; i++) {
      Particle<ndim>& part = hydro->GetParticlePointer(i);

      //for (int k=0; k<ndim; k++)
      for (k=0; k<ndim; k++) part.v[k] = (FLOAT) 0.0;
      part.m = mp;
      part.u = (FLOAT) 1.5;
      //part.h = hydro->h_fac*powf(mp/rho,invndim);
      part.ptype = gas_type;

    }

  }
  //-----------------------------------------------------------------------------------------------

  return;
}



//=================================================================================================
//  PolytropeIc::GetValue
/// ...
//=================================================================================================
template <int ndim>
FLOAT PolytropeIc<ndim>::GetValue
 (const std::string var,
  const FLOAT r[ndim])
{
  if (var == "x") {
    return r[0];
  }
  else if (ndim > 1 && var == "y") {
    return r[1];
  }
  else if (ndim > 2 && var == "z") {
    return r[2];
  }
  else if (var == "rho") {
    return (FLOAT) 0.0;
  }
  else {
    std::cout << "Invalid string variable for PolytropeIc::GetValue" << std::endl;
    return 0.0;
  }
}



//=================================================================================================
//  PolytropeIc::ComputeIsothermalLaneEmdenEquation
/// Create initial conditions for binary accretion simulation.
//=================================================================================================
template <int ndim>
void PolytropeIc<ndim>::ComputeIsothermalLaneEmdenSolution
 (const int Nmax,                            ///< ..
  const FLOAT delta_xi,                      ///< ..
  FLOAT *muArray,                            ///< ..
  FLOAT *phiArray,                           ///< ..
  FLOAT *psiArray,                           ///< ..
  FLOAT *xiArray)                            ///< ..
{
  int i;
  FLOAT k1_phi, k1_psi;
  FLOAT k2_phi, k2_psi;
  FLOAT k3_phi, k3_psi;
  FLOAT k4_phi, k4_psi;
  FLOAT phi;
  FLOAT psi;
  FLOAT xi;

  debug2("[PolytropeIc::ComputeIsothermalLaneEmdenSolution]");

  // Tabulate central values using boundary conditions
  xiArray[0]  = (FLOAT) 0.0;
  psiArray[0] = (FLOAT) 0.0;
  phiArray[0] = (FLOAT) 0.0;
  muArray[0]  = (FLOAT) 0.0;

  // Use first few terms of series solution for first step
  // (due to singularity in differential equation at xi = 0)
  xi  = delta_xi;
  psi = onesixth*pow(xi,2) - pow(xi,4)/(FLOAT) 120.0 + pow(xi,6)/(FLOAT) 1890.0;
  phi = onethird*xi - pow(xi,3)/(FLOAT) 30.0 + pow(xi,5)/(FLOAT) 315.0;
  xiArray[1]  = xi;
  psiArray[1] = psi;
  phiArray[1] = phi;
  muArray[1]  = phi*xi*xi;


  // Now loop over all over integration points
  //-----------------------------------------------------------------------------------------------
  for (i=2; i<Nmax; i++) {

    // Solve using 4th order Runge-Kutta method
    k1_phi = delta_xi*(exp(-psi) - (FLOAT) 2.0*phi/xi);
    k1_psi = delta_xi*phi;

    k2_phi = delta_xi*(exp(-psi - (FLOAT) 0.5*k1_psi) -
      (FLOAT) 2.0*(phi + (FLOAT) 0.5*k1_phi)/(xi + (FLOAT) 0.5*delta_xi));
    k2_psi = delta_xi*(phi + (FLOAT) 0.5*k1_phi);

    k3_phi = delta_xi*(exp(-psi - (FLOAT) 0.5*k2_psi) -
      (FLOAT) 2.0*(phi + (FLOAT) 0.5*k2_phi)/(xi + (FLOAT) 0.5*delta_xi));
    k3_psi = delta_xi*(phi + (FLOAT) 0.5*k2_phi);

    k4_phi = delta_xi*(exp(-psi - k3_psi) - (FLOAT) 2.0*(phi + k3_phi)/(xi + delta_xi));
    k4_psi = delta_xi*(phi + k3_phi);

    phi = phi + onesixth*(k1_phi + k4_phi) + onethird*(k2_phi + k3_phi);
    psi = psi + onesixth*(k1_psi + k4_psi) + onethird*(k2_psi + k3_psi);
    xi = (FLOAT) i*delta_xi;

    // Tabulate values
    xiArray[i]  = xi;
    psiArray[i] = psi;
    phiArray[i] = phi;
    muArray[i]  = phi*xi*xi;
  }
  //-----------------------------------------------------------------------------------------------


  return;
}



//=================================================================================================
//  PolytropeIc::ComputeLaneEmdenEquation
/// Create initial conditions for binary accretion simulation.
//=================================================================================================
template <int ndim>
void PolytropeIc<ndim>::ComputeLaneEmdenSolution
 (const int Nmax,                            ///< ..
  const FLOAT delta_xi,                      ///< ..
  const FLOAT npoly,                         ///< ..
  FLOAT *muArray,                            ///< ..
  FLOAT *phiArray,                           ///< ..
  FLOAT *psiArray,                           ///< ..
  FLOAT *xiArray)                            ///< ..
{
  int i;
  FLOAT k1_phi, k1_psi;
  FLOAT k2_phi, k2_psi;
  FLOAT k3_phi, k3_psi;
  FLOAT k4_phi, k4_psi;
  FLOAT phi;
  FLOAT psi;
  FLOAT xi;

  debug2("[PolytropeIc::ComputeLaneEmdenSolution]");

  // Tabulate central values using boundary conditions
  xiArray[0]  = (FLOAT) 0.0;
  psiArray[0] = (FLOAT) 0.0;
  phiArray[0] = (FLOAT) 0.0;
  muArray[0]  = (FLOAT) 0.0;

  // Use first few terms of series solution for first step
  // (due to singularity in differential equation at xi = 0)
  xi  = delta_xi;
  psi = onesixth*pow(xi,2) - pow(xi,4)/(FLOAT) 120.0 + pow(xi,6)/(FLOAT) 1890.0;
  phi = onethird*xi - pow(xi,3)/(FLOAT) 30.0 + pow(xi,5)/(FLOAT) 315.0;
  xiArray[1]  = xi;
  psiArray[1] = psi;
  phiArray[1] = phi;
  muArray[1]  = phi*xi*xi;


  // Now loop over all over integration points
  //-----------------------------------------------------------------------------------------------
  for (i=2; i<Nmax; i++) {

    // Solve using 4th order Runge-Kutta method
    k1_phi = delta_xi*(exp(-psi) - (FLOAT) 2.0*phi/xi);
    k1_psi = delta_xi*phi;

    k2_phi = delta_xi*(exp(-psi - (FLOAT) 0.5*k1_psi) -
      (FLOAT) 2.0*(phi + (FLOAT) 0.5*k1_phi)/(xi + (FLOAT) 0.5*delta_xi));
    k2_psi = delta_xi*(phi + (FLOAT) 0.5*k1_phi);

    k3_phi = delta_xi*(exp(-psi - (FLOAT) 0.5*k2_psi) -
      (FLOAT) 2.0*(phi + (FLOAT) 0.5*k2_phi)/(xi + (FLOAT) 0.5*delta_xi));
    k3_psi = delta_xi*(phi + (FLOAT) 0.5*k2_phi);

    k4_phi = delta_xi*(exp(-psi - k3_psi) - (FLOAT) 2.0*(phi + k3_phi)/(xi + delta_xi));
    k4_psi = delta_xi*(phi + k3_phi);

    phi = phi + onesixth*(k1_phi + k4_phi) + onethird*(k2_phi + k3_phi);
    psi = psi + onesixth*(k1_psi + k4_psi) + onethird*(k2_psi + k3_psi);
    xi = (FLOAT) i*delta_xi;

    // Tabulate values
    xiArray[i]  = xi;
    psiArray[i] = psi;
    phiArray[i] = phi;
    muArray[i]  = phi*xi*xi;
  }
  //-----------------------------------------------------------------------------------------------


  return;
}




template class PolytropeIc<1>;
template class PolytropeIc<2>;
template class PolytropeIc<3>;
