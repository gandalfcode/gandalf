//=================================================================================================
//  Sph.cpp
//  Contains important default routines for Sph class.
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
#include <math.h>
#include "Precision.h"
#include "Sph.h"
#include "SmoothingKernel.h"
#include "Particle.h"
#include "Parameters.h"
#include "EOS.h"
#include "Debug.h"
#include "InlineFuncs.h"
using namespace std;



//=================================================================================================
//  Sph::Sph
/// Constructor for parent SPH class.  Initialises important variables and
/// sets important parameters using initialialisation lists.
//=================================================================================================
template <int ndim>
Sph<ndim>::Sph(int _hydro_forces, int _self_gravity, FLOAT _alpha_visc, FLOAT _beta_visc,
               FLOAT _h_fac, FLOAT _h_converge, aviscenum _avisc, acondenum _acond,
               tdaviscenum _tdavisc, string _gas_eos, string _KernelName, int _size_sph,
               SimUnits &units, Parameters *params):
  Hydrodynamics<ndim>(_hydro_forces, _self_gravity, _h_fac, _gas_eos, _KernelName, _size_sph),
  size_sph_part(_size_sph),
  acond(_acond),
  avisc(_avisc),
  tdavisc(_tdavisc),
  alpha_visc(_alpha_visc),
  beta_visc(_beta_visc),
  h_converge(_h_converge),
  create_sinks(0),
  fixed_sink_mass(0)
{
  // Local references to parameter variables for brevity
  map<string, int> &intparams = params->intparams;
  map<string, double> &floatparams = params->floatparams;
  map<string, string> &stringparams = params->stringparams;
  string gas_radiation = stringparams["radiation"];

  // Set other class variables here
  Ngather = 0;
  hmin_sink = big_number;


  // Thermal physics object.  If energy equation is chosen, also initiate
  // the energy integration object.
  //-----------------------------------------------------------------------------------------------
  if ((_gas_eos == "energy_eqn" || _gas_eos == "constant_temp" ||
       _gas_eos == "isothermal" || _gas_eos == "barotropic" ||
       _gas_eos == "barotropic2") && gas_radiation == "ionisation") {
    eos = new IonisingRadiation<ndim>
      (_gas_eos, floatparams["temp0"], floatparams["mu_bar"],
       floatparams["gamma_eos"], floatparams["rho_bary"], &units);
  }
  else if ((_gas_eos == "energy_eqn" || _gas_eos == "constant_temp" ||
            _gas_eos == "isothermal" || _gas_eos == "barotropic" ||
            _gas_eos == "barotropic2") && gas_radiation == "monoionisation") {
    eos = new MCRadiationEOS<ndim>
      (_gas_eos, floatparams["temp0"], floatparams["temp_ion"], floatparams["mu_bar"],
       floatparams["mu_ion"], floatparams["gamma_eos"], floatparams["rho_bary"], &units);
  }
  else if (_gas_eos == "energy_eqn" || _gas_eos == "constant_temp") {
    eos = new Adiabatic<ndim>
      (floatparams["temp0"], floatparams["mu_bar"], floatparams["gamma_eos"]);
  }
  else if (_gas_eos == "isothermal") {
    eos = new Isothermal<ndim>
      (floatparams["temp0"], floatparams["mu_bar"], floatparams["gamma_eos"], &units);
  }
  else if (_gas_eos == "barotropic") {
    eos = new Barotropic<ndim>(floatparams["temp0"], floatparams["mu_bar"],
                               floatparams["gamma_eos"], floatparams["rho_bary"], &units);
  }
  else if (_gas_eos == "barotropic2") {
    eos = new Barotropic2<ndim>(floatparams["temp0"], floatparams["mu_bar"],
                                     floatparams["gamma_eos"], floatparams["rho_bary"], &units);
  }
  else if (_gas_eos == "rad_ws" || _gas_eos == "radws") {
    eos = new Radws<ndim>(floatparams["temp0"], floatparams["mu_bar"], floatparams["gamma_eos"]);
  }
  else {
    string message = "Unrecognised parameter : gas_eos = " + _gas_eos;
    ExceptionHandler::getIstance().raise(message);
  }


}



//=================================================================================================
//  Sph::InitialSmoothingLengthGuess
/// Perform initial guess of smoothing.  In the abscence of more sophisticated techniques, we guess
/// the smoothing length assuming a uniform density medium with the same volume and total mass.
//=================================================================================================
template <int ndim>
void Sph<ndim>::InitialSmoothingLengthGuess(void)
{
  int i;                           // Particle counter
  FLOAT h_guess;                   // Global guess of smoothing length
  FLOAT volume;                    // Volume of global bounding box
  FLOAT rmin[ndim];                // Min. extent of bounding box
  FLOAT rmax[ndim];                // Max. extent of bounding box

  debug2("[Sph::InitialSmoothingLengthGuess]");

  // Calculate bounding box containing all SPH particles
  this->ComputeBoundingBox(rmax,rmin,Nhydro);

  // Depending on the dimensionality, calculate the average smoothing
  // length assuming a uniform density distribution filling the bounding box.
  //-----------------------------------------------------------------------------------------------
  if (ndim == 1) {
    Ngather = (int) (2.0*kernp->kernrange*h_fac);
    volume = rmax[0] - rmin[0];
    h_guess = (volume*(FLOAT) Ngather)/(4.0*(FLOAT) Nhydro);
  }
  //-----------------------------------------------------------------------------------------------
  else if (ndim == 2) {
    Ngather = (int) (pi*pow(kernp->kernrange*h_fac,2));
    volume = (rmax[0] - rmin[0])*(rmax[1] - rmin[1]);
    h_guess = sqrtf((volume*(FLOAT) Ngather)/(4.0*(FLOAT) Nhydro));
  }
  //-----------------------------------------------------------------------------------------------
  else if (ndim == 3) {
    Ngather = (int) (4.0*pi*pow(kernp->kernrange*h_fac,3)/3.0);
    volume = (rmax[0] - rmin[0])*(rmax[1] - rmin[1])*(rmax[2] - rmin[2]);
    h_guess = powf((3.0*volume*(FLOAT) Ngather)/(32.0*pi*(FLOAT) Nhydro),onethird);
  }
  //-----------------------------------------------------------------------------------------------

  // Set all smoothing lengths equal to average value
  for (i=0; i<Nhydro; i++) {
    SphParticle<ndim>& part = GetSphParticlePointer(i);
    part.h = h_guess;
    part.invh = 1.0/h_guess;
    part.hrangesqd = kernfacsqd*kernp->kernrangesqd*part.h*part.h;
  }

  return;
}



template class Sph<1>;
template class Sph<2>;
template class Sph<3>;
