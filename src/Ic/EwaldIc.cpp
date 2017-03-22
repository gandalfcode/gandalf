//=================================================================================================
//  EwaldIc.cpp
//  Class for generating initial conditions for ...
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
//  EwaldIc::EwaldIc
/// Set-up SILCC-type simulation initial conditions.
//=================================================================================================
template <int ndim>
EwaldIc<ndim>::EwaldIc(Simulation<ndim>* _sim, FLOAT _invndim) :
  Ic<ndim>(_sim, _invndim)
{
  if (simparams->intparams["ndim"] != 3) {
    ExceptionHandler::getIstance().raise("Ewald tests can only be generated in 3d");
  }
}



//=================================================================================================
//  EwaldIc::Generate
/// ...
//=================================================================================================
template <int ndim>
void EwaldIc<ndim>::Generate(void)
{
  // Only compile for 3 dimensions
  //-----------------------------------------------------------------------------------------------
  if (ndim == 3) {

    int i,k;                             // Particle and dimension counters
    int Nlattice1[3];                    // Lattice size
    int Npart;                           // No. of particles in lattice
    int periodicity;                     // orientation of the axis of symmetry
    FLOAT csound;                        // (Isothermal) sound speed
    FLOAT h0;                            // Slab scale height
    FLOAT a2inv;                         // Squared inverse scale height for cylinder
    FLOAT lambda;                        // Wavelength of perturbation
    FLOAT kwave;                         // Wave number of perturbing sound wave
    FLOAT ugas;                          // Internal energy of gas
    FLOAT volume;                        // Simulation box volume
    FLOAT *r;                            // Particle positions
    Nbody<ndim>* nbody = sim->nbody;     // Pointer to N-body object

    // Make local copies of parameters for setting up problem
    Nlattice1[0]    = simparams->intparams["Nlattice1[0]"];
    Nlattice1[1]    = simparams->intparams["Nlattice1[1]"];
    Nlattice1[2]    = simparams->intparams["Nlattice1[2]"];
    string x_boundary_lhs = simparams->stringparams["boundary_lhs[0]"];
    string x_boundary_rhs = simparams->stringparams["boundary_rhs[0]"];
    string y_boundary_lhs = simparams->stringparams["boundary_lhs[1]"];
    string y_boundary_rhs = simparams->stringparams["boundary_rhs[1]"];
    string z_boundary_lhs = simparams->stringparams["boundary_lhs[2]"];
    string z_boundary_rhs = simparams->stringparams["boundary_rhs[2]"];
    string ic       = simparams->stringparams["ic"];
    string simType  = simparams->stringparams["sim"];
    FLOAT rhofluid1 = simparams->floatparams["rhofluid1"];
    FLOAT press1    = simparams->floatparams["press1"];
    FLOAT gamma     = simparams->floatparams["gamma_eos"];
    FLOAT gammaone  = gamma - (FLOAT) 1.0;
    FLOAT amp       = simparams->floatparams["amp"];
    FLOAT temp0     = simparams->floatparams["temp0"];
    FLOAT mu_bar    = simparams->floatparams["mu_bar"];
    //FLOAT zmax      = simparams->floatparams["zmax"];

    debug2("[EwaldIc::Generate]");

    if (hydro->gas_eos == "isothermal") {
      ugas   = temp0/gammaone/mu_bar;
      press1 = gammaone*rhofluid1*ugas;
      csound = sqrt(press1/rhofluid1);
    }
    else {
      ugas   = press1/rhofluid1/gammaone;
      csound = sqrt(gamma*press1/rhofluid1);
    }

    lambda = icBox.size[0];
    volume = icBox.size[0]*icBox.size[1]*icBox.size[2];
    kwave  = twopi/lambda;
    //omegawave = twopi*csound/lambda;

    // Set particle number (depending on simulation type) and allocate memory
    Npart = Nlattice1[0]*Nlattice1[1]*Nlattice1[2];
    if (simType == "nbody") {
      nbody->Nstar = Npart;
    }
    else {
      hydro->Nhydro = Npart;
    }
    sim->AllocateParticleMemory();
    r = new FLOAT[ndim*Npart];

    // Determine orientation of the problem (variable periodicity).  Controls orientation
    // of density structure in the case of mixed boundary conditions
    // (i.e. normal to the plane or axis of the cylinder)
    periodicity = 0;
    if (x_boundary_lhs == "periodic" && x_boundary_rhs == "periodic") {
      periodicity += 1;
    }
    if (y_boundary_lhs == "periodic" && y_boundary_rhs == "periodic") {
      periodicity += 2;
    }
    if (z_boundary_lhs == "periodic" && z_boundary_rhs == "periodic") {
      periodicity += 4;
    }


    // 1D sinusoidal density perturbation
    //=============================================================================================
    if (ic == "ewaldsine" || ic == "jeans") {

      lambda = icBox.max[0] - icBox.min[0];
      kwave  = twopi/lambda;

      // Add regular distribution of SPH particles
      Ic<ndim>::AddCubicLattice(Npart, Nlattice1, icBox, false, r);

      // Add sinusoidal density perturbation to particle distribution
      Ic<ndim>::AddSinusoidalDensityPerturbation(Npart, amp, lambda, r);

      // Create star-lattice for N-body simulations.  Otherwise hydro particles
      //-------------------------------------------------------------------------------------------
      if (simType == "nbody") {

        // Add stars
        for (int i=0; i<Npart; i++) {
          for (int k=0; k<ndim; k++) nbody->stardata[i].v[k] = (FLOAT) 0.0;
          for (int k=0; k<ndim; k++) nbody->stardata[i].r[k] = r[ndim*i + k];
          nbody->stardata[i].m    = rhofluid1*volume/(FLOAT) Npart;
          nbody->stardata[i].h    = (FLOAT) 0.1 / (FLOAT) Npart;
          nbody->stardata[i].invh = (FLOAT) 1.0 / nbody->stardata[i].h;
        }

      }
      //-------------------------------------------------------------------------------------------
      else {

        // Set all particle quantities
        for (i=0; i<Npart; i++) {
          Particle<ndim>& part = hydro->GetParticlePointer(i);

          // Set positions in main array with corresponind velocity perturbation
          for (k=0; k<ndim; k++) part.r[k] = r[ndim*i + k];
          for (k=0; k<ndim; k++) part.v[k] = (FLOAT) 0.0;
          part.m = rhofluid1*volume/(FLOAT) Npart;
          part.h = hydro->h_fac*pow(part.m/rhofluid1,invndim);

          if (hydro->gas_eos == "isothermal") part.u = temp0/gammaone/mu_bar;
          else part.u = press1/rhofluid1/gammaone;
        }

      }
      //-------------------------------------------------------------------------------------------


    }
    //=============================================================================================
    else if (ic == "ewaldsine2") {

      lambda = icBox.max[0] - icBox.min[0];
      kwave = twopi/lambda;

      // Add regular distribution of SPH particles
      Ic<ndim>::AddCubicLattice(Npart, Nlattice1, icBox, false, r);

      volume = icBox.size[0]*icBox.size[1]*icBox.size[2];
      FLOAT volp = volume/(FLOAT) Npart;

      // Set all other particle quantities
      //-------------------------------------------------------------------------------------------
      for (i=0; i<Npart; i++) {
        Particle<ndim>& part = hydro->GetParticlePointer(i);

        // Set positions in main array with corresponind velocity perturbation
        for (k=0; k<ndim; k++) part.r[k] = r[ndim*i + k];
        for (k=0; k<ndim; k++) part.v[k] = 0.0;
        part.rho = rhofluid1*(1.0+amp*sin(kwave*part.r[0]));
        part.m   = part.rho*volp;
        part.h   = hydro->h_fac*pow(part.m/part.rho,invndim);

        if (hydro->gas_eos == "isothermal") part.u = temp0/gammaone/mu_bar;
        else part.u = press1/rhofluid1/gammaone;

      }
      //-------------------------------------------------------------------------------------------

    }
    //=============================================================================================
    else if (ic == "ewaldslab") {

      // it is assumed that the slab is selfgravitating, so boundary conditions for gravity
      // must be set properly
      if ((periodicity != 3) && (periodicity != 5) && (periodicity != 6)) {
        string msg = "For this test boundary conditions must be periodic in two directions";
        ExceptionHandler::getIstance().raise(msg);
      }

      h0 = csound/sqrtf(twopi*rhofluid1);

      // Add regular distribution of SPH particles
      Ic<ndim>::AddCubicLattice(Npart, Nlattice1, icBox, false, r);

      volume = icBox.size[0]*icBox.size[1]*icBox.size[2];
      FLOAT volp = volume/(FLOAT) Npart;

      // Set all other particle quantities
      //-------------------------------------------------------------------------------------------
      for (i=0; i<Npart; i++) {
        Particle<ndim>& part = hydro->GetParticlePointer(i);

        // Set positions in main array with corresponind velocity perturbation
        for (k=0; k<ndim; k++) part.r[k] = r[ndim*i + k];
        for (k=0; k<ndim; k++) part.v[k] = 0.0;

        // normal to the slab depends on chosen BCs
        switch (periodicity) {
          case 3:
            part.rho = rhofluid1/powf(cosh(part.r[2]/h0),2);
            break;
          case 5:
            part.rho = rhofluid1/powf(cosh(part.r[1]/h0),2);
            break;
          case 6:
            part.rho = rhofluid1/powf(cosh(part.r[0]/h0),2);
            break;
        }
        part.m   = part.rho*volp;
        part.h   = hydro->h_fac*pow(part.m/part.rho,invndim);

        if (hydro->gas_eos == "isothermal") part.u = temp0/gammaone/mu_bar;
        else part.u = press1/rhofluid1/gammaone;

      }

      //-------------------------------------------------------------------------------------------


    }
    //=============================================================================================
    else if (ic == "ewaldcylinder") {

      // it is assumed that cylinder is selfgravitating, so boundary conditions for gravity
      // must be set properly
      if ((periodicity != 1) && (periodicity != 2) && (periodicity != 4)) {
        string msg="For this test boundary conditions must be periodic in one direction";
        ExceptionHandler::getIstance().raise(msg);
      }

      a2inv = pi*rhofluid1*0.5/pow(csound,2);

      // Add regular distribution of SPH particles
      Ic<ndim>::AddCubicLattice(Npart, Nlattice1, icBox, false, r);

      volume = icBox.size[0]*icBox.size[1]*icBox.size[2];
      FLOAT volp = volume/(FLOAT) Npart;

      // Set all other particle quantities
      //-------------------------------------------------------------------------------------------
      for (i=0; i<Npart; i++) {
        Particle<ndim>& part = hydro->GetParticlePointer(i);

        // Set positions in main array with corresponind velocity perturbation
        for (k=0; k<ndim; k++) part.r[k] = r[ndim*i + k];
        for (k=0; k<ndim; k++) part.v[k] = 0.0;
        switch (periodicity) {
          case 1:
            part.rho = rhofluid1/pow((1.0+a2inv*(pow(part.r[1],2) + pow(part.r[2],2))),2);
            break;
          case 2:
            part.rho = rhofluid1/pow((1.0+a2inv*(pow(part.r[0],2) + pow(part.r[2],2))),2);
            break;
          case 4:
            part.rho = rhofluid1/pow((1.0+a2inv*(pow(part.r[0],2) + pow(part.r[1],2))),2);
            break;
        }
        part.m   = part.rho*volp;
        part.h   = hydro->h_fac*pow(part.m/part.rho,invndim);

        if (hydro->gas_eos == "isothermal") part.u = temp0/gammaone/mu_bar;
        else part.u = press1/rhofluid1/gammaone;

      }

      //-------------------------------------------------------------------------------------------

    }
    //=============================================================================================

    sim->initial_h_provided = true;

  }
  //-----------------------------------------------------------------------------------------------

  return;
}



template class EwaldIc<1>;
template class EwaldIc<2>;
template class EwaldIc<3>;
