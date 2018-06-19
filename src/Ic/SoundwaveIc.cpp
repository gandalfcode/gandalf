//=================================================================================================
//  SoundwaveIc.cpp
//  Class for generating initial conditions for simple turbulent core simulations.
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
//  SoundwaveIc::SoundwaveIc
/// Set-up SILCC-type simulation initial conditions.
//=================================================================================================
template <int ndim>
SoundwaveIc<ndim>::SoundwaveIc(Simulation<ndim>* _sim, FLOAT _invndim) :
  Ic<ndim>(_sim, _invndim)
{
  if (simparams->intparams["dimensionless"] == 0) {
    ExceptionHandler::getIstance().raise("dimensionless units required");
  }
}



//=================================================================================================
//  Silcc::Generate
/// Set-up SILCC-type simulation initial conditions.
//=================================================================================================
template <int ndim>
void SoundwaveIc<ndim>::Generate(void)
{
  int i,k;                                       // Particle and dimension counters
  int Nlattice[ndim];                            // Lattice size
  int Npart;                                     // No. of particles
  FLOAT csound;                                  // Sound speed
  FLOAT lambda;                                  // Wavelength of perturbation
  FLOAT kwave;                                   // Wave number of perturbing sound wave
  FLOAT ugas;                                    // Internal energy of gas
  FLOAT volume;                                  // Volume of simulation box
  FLOAT *r;                                      // Particle positions

  // Make local copies of parameters for setting up problem
  Nlattice[0]     = simparams->intparams["Nlattice1[0]"];
  Nlattice[1]     = simparams->intparams["Nlattice1[1]"];
  Nlattice[2]     = simparams->intparams["Nlattice1[2]"];
  FLOAT rhofluid1 = simparams->floatparams["rhofluid1"];
  FLOAT press1    = simparams->floatparams["press1"];
  FLOAT gamma     = simparams->floatparams["gamma_eos"];
  FLOAT gammaone  = gamma - (FLOAT) 1.0;
  FLOAT amp       = simparams->floatparams["amp"];
  FLOAT temp0     = simparams->floatparams["temp0"];
  FLOAT mu_bar    = simparams->floatparams["mu_bar"];
  std::string particle_dist = simparams->stringparams["particle_distribution"];

  debug2("[SoundwaveIc::Generate]");

  if (hydro->gas_eos == "isothermal") {
    ugas   = temp0/gammaone/mu_bar;
    press1 = gammaone*rhofluid1*ugas;
    csound = sqrt(press1/rhofluid1);
  }
  else {
    ugas   = press1/rhofluid1/gammaone;
    csound = sqrt(gamma*press1/rhofluid1);
  }
  lambda = icBox.max[0] - icBox.min[0];
  kwave = twopi/lambda;

  if (ndim == 1) {
    volume = icBox.max[0] - icBox.min[0];
    Npart  = simparams->intparams["Nhydro"];
    Nlattice[0] = Npart;
  }
  else if (ndim == 2) {
    volume = (icBox.max[0] - icBox.min[0])*(icBox.max[1] - icBox.min[1]);
    Npart = Nlattice[0]*Nlattice[1];
  }
  else if (ndim == 3) {
    volume = (icBox.max[0] - icBox.min[0])*
      (icBox.max[1] - icBox.min[1])*(icBox.max[2] - icBox.min[2]);
    Npart = Nlattice[0]*Nlattice[1]*Nlattice[2];
  }
  FLOAT mp = rhofluid1*volume/(FLOAT) Npart;

  hydro->Nhydro = Npart;
  bool dusty_wave = simparams->stringparams["dust_forces"] != "none" ;
  if (dusty_wave) hydro->Nhydro *= 2;

  sim->AllocateParticleMemory();
  r = new FLOAT[ndim*Npart];

  // Add regular distribution of SPH particles
  Ic<ndim>::AddCubicLattice(Npart, Nlattice, icBox, false, r);

  // Add sinusoidal density perturbation to particle distribution
  Ic<ndim>::AddSinusoidalDensityPerturbation(Npart, amp, lambda, r);

  // Set all other particle quantities
  for (i=0; i<Npart; i++) {
    Particle<ndim>& part = hydro->GetParticlePointer(i);
    for (k=0; k<ndim; k++) part.r[k] = r[ndim*i + k];
    for (k=0; k<ndim; k++) part.v[k] = (FLOAT) 0.0;
    part.v[0] = csound*amp*sin(kwave*r[ndim*i]);
    part.m    = mp;
    part.h    = hydro->h_fac*pow(part.m/rhofluid1,invndim);
    part.u    = ugas;
    part.rho  = rhofluid1*(1.0 + amp*sin(twopi*part.r[0]/lambda));
  }

  if (dusty_wave) {
    FLOAT d2g = simparams->floatparams["dust_mass_factor"] ;
    for (i = 0; i < Npart; ++i) {
      Particle<ndim>& Pg = hydro->GetParticlePointer(i) ;
      Particle<ndim>& Pd = hydro->GetParticlePointer(i+Npart) ;
      Pd = Pg ;
      Pd.m *= d2g ;
      Pd.h_dust = Pd.h ;
      Pd.u = 0 ;
      Pg.ptype = gas_type ;
      Pd.ptype = dust_type ;
    }

    sim->initial_h_provided = true;
  }

  delete[] r;

  return;
}



template class SoundwaveIc<1>;
template class SoundwaveIc<2>;
template class SoundwaveIc<3>;
