//=================================================================================================
//  BossBodenheimerIc.cpp
//  Class for generating initial conditions for Boss-Bodenheimer cloud simulations.
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
//  BossBodenheimerIc::BossBodenheimerIc
/// Set-up BossBodenheimer-type simulation initial conditions.
//=================================================================================================
template <int ndim>
BossBodenheimerIc<ndim>::BossBodenheimerIc(Simulation<ndim>* _sim, Hydrodynamics<ndim>* _hydro, FLOAT _invndim) :
  Ic<ndim>(_sim, _hydro, _invndim)
{
  // Some sanity checking to ensure correct dimensionality and (dimensionless) units are used
  if (simparams->intparams["ndim"] != 3) {
    ExceptionHandler::getIstance().raise("Boss-Bodenheimer test only runs in 3D");
  }
}



//=================================================================================================
//  BossBodenheimer::Generate
/// Set-up BossBodenheimer-type simulation initial conditions.
//=================================================================================================
template <int ndim>
void BossBodenheimerIc<ndim>::Generate(void)
{
  // Only compile for 3-dimensional case
  //-----------------------------------------------------------------------------------------------
  if (ndim == 3) {

    int i;                               // Particle counter
    int k;                               // Dimension counter
    int Nsphere;                         // Actual number of particles in sphere
    FLOAT mp;                            // Mass of one particle
    FLOAT rcentre[ndim];                 // Position of sphere centre
    FLOAT rho;                           // Fluid density
    FLOAT *r;                            // Positions of all particles
    FLOAT *v;                            // Velocities of all particles

    // Create local copies of initial conditions parameters
    int Npart      = simparams->intparams["Nhydro"];
    FLOAT amp      = simparams->floatparams["amp"];
    FLOAT angvel   = simparams->floatparams["angvel"];
    FLOAT mcloud   = simparams->floatparams["mcloud"];
    FLOAT mu_bar   = simparams->floatparams["mu_bar"];
    FLOAT press    = simparams->floatparams["press1"];
    FLOAT radius   = simparams->floatparams["radius"];
    FLOAT temp0    = simparams->floatparams["temp0"];
    FLOAT gammaone = simparams->floatparams["gamma_eos"] - 1.0;
    string particle_dist = simparams->stringparams["particle_distribution"];

    debug2("[BossBodenheimerIc::Generate]");

    // Convert any parameters to code units
    angvel /= simunits.angvel.outscale;
    mcloud /= simunits.m.outscale;
    press  /= simunits.press.outscale;
    radius /= simunits.r.outscale;
    temp0  /= simunits.temp.outscale;

    r = new FLOAT[ndim*Npart];
    v = new FLOAT[ndim*Npart];

    // Add a sphere of random particles with origin 'rcentre' and radius 'radius'
    for (k=0; k<ndim; k++) rcentre[k] = (FLOAT) 0.0;

    // Create the sphere depending on the choice of initial particle distribution
    if (particle_dist == "random") {
      Ic<ndim>::AddRandomSphere(Npart, rcentre, radius, r, sim->randnumb);
    }
    else if (particle_dist == "cubic_lattice" || particle_dist == "hexagonal_lattice") {
      Nsphere = Ic<ndim>::AddLatticeSphere(Npart, rcentre, radius, particle_dist, r, sim->randnumb);
      if (Nsphere != Npart) {
        cout << "Warning! Unable to converge to required "
             << "no. of ptcls due to lattice symmetry" << endl;
      }
      Npart = Nsphere;
    }
    else {
      string message = "Invalid particle distribution option";
      ExceptionHandler::getIstance().raise(message);
    }

    // Allocate local and main particle memory
    hydro->Nhydro = Npart;

    bool dusty_collapse =
          simparams->stringparams["dust_forces"] != "none" ;

    if (dusty_collapse)
      hydro->Nhydro *= 2 ;

    sim->AllocateParticleMemory();


    mp = mcloud / (FLOAT) Npart;
    rho = (FLOAT) 3.0*mcloud / ((FLOAT)4.0*pi*pow(radius,3));

    // Perturb positions of particles in cloud
    Ic<ndim>::AddAzimuthalDensityPerturbation(Npart, 2, amp, rcentre, r);

    // Add solid-body rotational velocity field
    Ic<ndim>::AddRotationalVelocityField(Npart, angvel, rcentre, r, v);

    // Record particle properties in main memory
    //#pragma omp parallel for default(none) shared(gammaone,Npart,mp,mu_bar,r,rho,temp0,v) private(i,k)
    for (i=0; i<Npart; i++) {
      Particle<ndim>& part = hydro->GetParticlePointer(i);
      for (k=0; k<ndim; k++) part.r[k] = r[ndim*i + k];
      for (k=0; k<ndim; k++) part.v[k] = v[ndim*i + k];
      part.m = mp;
      part.h = hydro->h_fac*pow(part.m/rho,invndim);
      part.u = temp0/gammaone/mu_bar;
      //if (hydro->gas_eos == "isothermal" || hydro->gas_eos == "barotropic")
      //  part.u = temp0/gammaone/mu_bar;
      //else
      //  part.u = press/rho/gammaone;
    }


    if (dusty_collapse){
      FLOAT d2g = simparams->floatparams["dust_mass_factor"] ;
      for (i = 0; i < Npart; ++i){
        Particle<ndim>& Pg = hydro->GetParticlePointer(i) ;
        Particle<ndim>& Pd = hydro->GetParticlePointer(i+Npart) ;
        Pd = Pg ;
        Pd.m *= d2g ;


        for (k=0; k < ndim; k++)
        Pd.r[k] += 0.01 * Pd.h ;


        Pd.h_dust = Pd.h ;
        Pd.u = 0 ;

        Pg.ptype = gas_type ;
        Pd.ptype = dust_type ;
      }
    }


    sim->initial_h_provided = true;

    delete[] v;
    delete[] r;


  }
  //-----------------------------------------------------------------------------------------------

  return;
}



template class BossBodenheimerIc<1>;
template class BossBodenheimerIc<2>;
template class BossBodenheimerIc<3>;
