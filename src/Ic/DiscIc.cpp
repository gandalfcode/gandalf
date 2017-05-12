//=================================================================================================
//  DiscIc.cpp
//  Class for generating initial conditions for accretion discs
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

#include <cstdlib>
#include <cmath>
#include "Precision.h"
#include "Debug.h"
#include "Ic.h"


//=================================================================================================
//  DiscIc::DiscIc
/// Class constructor
//=================================================================================================
template <int ndim>
DiscIc<ndim>::DiscIc(Simulation<ndim>* _sim, FLOAT _invndim) :
  Ic<ndim>(_sim, _invndim)
{
  if (simparams->intparams["dimensionless"] == 0) {
    ExceptionHandler::getIstance().raise("dimensionless units required");
  }
  if (simparams->intparams["ndim"] == 1){
    ExceptionHandler::getIstance().raise("The disc IC require at least 2 dimensions!");
  }
  if (simparams->intparams["sink_particles"] == 0) {
    ExceptionHandler::getIstance().raise("You can't run a disc simulation without sink accretion!");
  }
}

//=================================================================================================
//  DiscIc::Generate
/// Generate an accretion disc
//=================================================================================================
template <int ndim>
void DiscIc<ndim>::Generate(void)
{
  debug2("[DiscIc::Generate");

  // Local copy of important parameters
  const int Npart = simparams->intparams["Nhydro"];
  const FLOAT gammaone = simparams->floatparams["gamma_eos"] - 1.0;
  const FLOAT mass = simparams->floatparams["DiscIcMass"];
  const FLOAT p = simparams->floatparams["DiscIcP"];
  const FLOAT q = simparams->floatparams["DiscIcQ"];
  const FLOAT rin = simparams->floatparams["DiscIcRin"];
  const FLOAT rout = simparams->floatparams["DiscIcRout"];
  const FLOAT H_r = simparams->floatparams["DiscIcHr"];

  // Assume units where GM_\ast=1
  const FLOAT GStar_M = 1;


  // Allocate global and local memory for all particles
  hydro->Nhydro = Npart;
  if (simparams->intparams["DiscIcPlanet"]==0) {
    sim->nbody->Nstar = 1;
  }
  else {
    sim->nbody->Nstar = 2;
  }

  sim->AllocateParticleMemory();

  // Compute the particle mass
  const FLOAT massi = mass/Npart;

  // Normalisation of sound speed at rin
  const FLOAT cs0 = H_r*std::sqrt(GStar_M/rin);

  // Normalization of the surface density at rin
  const FLOAT sig0=(2-p)*mass*std::pow(rin,-p)/(2*pi) / (pow(rout,2-p)-pow(rin,2-p));

  FLOAT f_max;
  if (p <=1 ) {
    f_max = pow(rout/rin,-(p-1));
  }
  else {
    f_max = 1;
  }

  // Generate particle positions, velocity and internal energies
  for (int i=0; i<Npart; i++) {
    Particle<ndim>& part = hydro->GetParticlePointer(i);

    // Set mass
    part.m = massi;

    // ---------------------------------------
    // POSITIONS
    // ---------------------------------------

    // Phi
    const FLOAT phi = 2*M_PI*std::rand()/RAND_MAX;

    // Radius, constructed from r*sigma
    FLOAT f_test, r, f;
    do {
      r = rin + (rout - rin) * std::rand()/RAND_MAX;
      f_test = std::rand()/RAND_MAX*f_max;
      f = std::pow(r/rin,-(p-1.) );
    } while (f_test >= f);

    // Set x and y
    part.r[0] = r*sin(phi);
    part.r[1] = r*cos(phi);


    // Z, using gaussian profile
    const FLOAT cs = cs0*pow(r/rin,-q);
    const FLOAT H = sqrt(2.)*cs*pow(r,1.5)/sqrt(GStar_M);
    const FLOAT z_min = -3*H;
    const FLOAT z_max = 3*H;
    const FLOAT zf_max = 1;
    FLOAT zf_test=zf_max;
    FLOAT z, zf;
    do {
      z = z_min+(z_max-z_min)*std::rand()/RAND_MAX;
      zf_test = zf_max*std::rand()/RAND_MAX;
      zf = std::exp(-pow(z/H,2));
    } while (zf_test >= zf);

    if (ndim == 3) {
      part.r[2] = z;
    }

    // ---------------------------------------
    // VELOCITIES
    // ---------------------------------------

    // Keplerian velocity
    const FLOAT vk = std::sqrt(GStar_M/r);
    const FLOAT vphi = vk*std::sqrt(1.-0.5*(H/r)*(H/r)*(1.5+p+q));
    part.v[0] = -vphi*cos(phi);
    part.v[1] = vphi*sin(phi);
    if (ndim == 3) {
      part.v[2] = 0;
    }

    // ---------------------------------------
    // DENSITY
    // ---------------------------------------

    // Rough density for smoothing length and density
    const FLOAT sigma = sig0*f/(r/rin);

    const FLOAT rhoz = zf/(H*std::sqrt(M_PI));
    if (ndim == 3) {
      part.rho = sigma*rhoz;
    }
    else {
      part.rho = sigma;
    }
    part.h = hydro->h_fac*std::pow(part.m/part.rho,1./ndim);

    // ---------------------------------------
    // THERMAL PHYSICS
    // ---------------------------------------
    part.sound = cs;
    part.u = part.sound*part.sound/gammaone;

  }


  sim->initial_h_provided = true;

  Nbody<ndim>* nbody = sim->nbody;     // Pointer to Nbody object

  // Set up the star
  StarParticle<ndim>& star = nbody->stardata[0];
  for (int k=0; k<ndim; k++) star.r[k] = (FLOAT) 0.0;
  for (int k=0; k<ndim; k++) star.v[k] = (FLOAT) 0.0;
  star.m = 1;
  star.h = rin/hydro->kernp->kernrange;

  // Set up the planet
  if (simparams->intparams["DiscIcPlanet"] == 1) {
    const FLOAT e = simparams->floatparams["DiscIcPlanetEccen"];
    const FLOAT rp = simparams->floatparams["DiscIcPlanetRadius"];
    const FLOAT i = simparams->floatparams["DiscIcPlanetIncl"];
    StarParticle<ndim>& planet = nbody->stardata[1];
    planet.r[0] = rp*(1.+e);
    planet.r[1] = 0.0;
    if (ndim==3) planet.r[2] = 0.0;
    planet.v[0] = 0.0;
    planet.v[1] = 1.0/std::sqrt(rp)*std::sqrt( (1.0-e)/(1.0+e))* \
        std::cos(i*M_PI/180.0);
    if (ndim==3) planet.v[2] = planet.v[1]*std::sin(i*M_PI/180.0)/  \
        std::cos(i*M_PI/180.0);;
    planet.m = simparams->floatparams["DiscIcPlanetMass"];
    planet.h = simparams->floatparams["DiscIcPlanetAccretionRadiusHill"]*
        rp*pow(simparams->floatparams["DiscIcPlanetMass"],1./3)/hydro->kernp->kernrange;
  }


}

template class DiscIc<1>;
template class DiscIc<2>;
template class DiscIc<3>;

