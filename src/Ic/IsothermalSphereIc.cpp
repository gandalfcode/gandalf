//=================================================================================================
//  IsothermalSphereIc.cpp
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
//  IsothermalSphereIc::IsothermalSphereIc
/// Set-up SILCC-type simulation initial conditions.
//=================================================================================================
template <int ndim>
IsothermalSphereIc<ndim>::IsothermalSphereIc(Simulation<ndim>* _sim, FLOAT _invndim) :
  Ic<ndim>(_sim, _invndim)
{
}



//=================================================================================================
//  IsothermalSphereIc::Generate
/// ...
//=================================================================================================
template <int ndim>
void IsothermalSphereIc<ndim>::Generate(void)
{
  string ic = simparams->stringparams["ic"];

  //-----------------------------------------------------------------------------------------------
  if (ic == "isothermsphere") {
    int i,k;                          // Particle and dimension counters
    FLOAT rcentre[ndim];              // Position of sphere centre
    FLOAT volume;                     // Volume of sphere
    FLOAT *r;                         // Particle position vectors

    // Local copies of important parameters
    int Npart            = simparams->intparams["Nsph"];
    FLOAT mcloud         = simparams->floatparams["mcloud"];
    FLOAT radius         = simparams->floatparams["radius"];
    FLOAT rhofluid       = simparams->floatparams["rhofluid1"];
    FLOAT press          = simparams->floatparams["press1"];
    FLOAT temp0          = simparams->floatparams["temp0"];
    FLOAT mu_bar         = simparams->floatparams["mu_bar"];
    FLOAT gammaone       = simparams->floatparams["gamma_eos"] - 1.0;
    string particle_dist = simparams->stringparams["particle_distribution"];

    debug2("[Ic::IsothermSphere]");

    mcloud   /= simunits.m.outscale;
    radius   /= simunits.r.outscale;
    rhofluid /= simunits.rho.outscale;
    press    /= simunits.press.outscale;
    r = new FLOAT[ndim*Npart];

    // Add a sphere of random particles with origin 'rcentre' and radius 'radius'
    for (k=0; k<ndim; k++) rcentre[k] = (FLOAT) 0.0;
    Ic<ndim>::Addr2Sphere(Npart, r, rcentre, radius);

    hydro->Nhydro = Npart;
    sim->AllocateParticleMemory();

    if (ndim == 1) volume = 2.0*radius;
    else if (ndim == 2) volume = pi*radius*radius;
    else if (ndim == 3) volume = 4.0*onethird*pi*pow(radius,3);
    //if (mcloud > small_number && radius > small_number)
    //  rhofluid = mcloud / volume;
    rhofluid = mcloud / volume;


    // Record particle positions and initialise all other variables
  #pragma omp parallel for default(none)\
    shared(gammaone,mcloud,Npart,press,r,rhofluid,volume,temp0,mu_bar) private(i,k)
    for (i=0; i<hydro->Nhydro; i++) {
      Particle<ndim>& part = hydro->GetParticlePointer(i);
      for (k=0; k<ndim; k++) {
        part.r[k] = r[ndim*i + k];
        part.v[k] = (FLOAT) 0.0;
        part.a[k] = (FLOAT) 0.0;
      }
      if (hydro->gas_eos == "isothermal") {
        part.u = temp0/gammaone/mu_bar;
      }
      else {
        part.u = press/rhofluid/gammaone;
      }
      part.m = mcloud / (FLOAT) Npart;
      part.h = hydro->h_fac*pow(part.m/rhofluid,invndim);

      part.iorig = i;
    }

    sim->initial_h_provided = true;

    delete[] r;

  }
  //-----------------------------------------------------------------------------------------------
  else if (ic == "rotisothermsphere") {
    int i,k;                             // Particle and dimension counters
    FLOAT rcentre[ndim];                 // Position of sphere centre
    FLOAT volume;                        // Volume of sphere
    FLOAT *r;                            // Particle position vectors
    FLOAT r_perp[3];                     // Perpendicular vector
    FLOAT r_center;                      // Radius to centre
    FLOAT norm;                          // ..
    FLOAT pc = 3.08567758*pow(10,16);    // ..

    // Local copies of important parameters
    int Npart            = simparams->intparams["Nsph"];
    FLOAT mcloud         = simparams->floatparams["mcloud"];
    FLOAT radius         = simparams->floatparams["radius"];
    FLOAT rhofluid       = simparams->floatparams["rhofluid1"];
    FLOAT press          = simparams->floatparams["press1"];
    FLOAT temp0          = simparams->floatparams["temp0"];
    FLOAT mu_bar         = simparams->floatparams["mu_bar"];
    FLOAT gammaone       = simparams->floatparams["gamma_eos"] - 1.0;
    FLOAT omega          = simparams->floatparams["omega"];
    FLOAT omega_inv      = 1.0/omega;
    string particle_dist = simparams->stringparams["particle_distribution"];

    debug2("[Ic::RotIsothermSphere]");

    mcloud    /= simunits.m.outscale;
    radius    /= simunits.r.outscale;
    rhofluid  /= simunits.rho.outscale;
    press     /= simunits.press.outscale;
    omega_inv /= simunits.t.outscale;
    std::cout << omega << endl;
    omega = 1 / omega_inv;
    std::cout << omega << endl;
    r = new FLOAT[ndim*Npart];

    // Add a sphere of random particles with origin 'rcentre' and radius 'radius'
    for (k=0; k<ndim; k++) rcentre[k] = (FLOAT) 0.0;
    Ic<ndim>::Addr2Sphere(Npart, r, rcentre, radius);

    hydro->Nhydro = Npart;
    sim->AllocateParticleMemory();

    if (ndim == 1) volume = 2.0*radius;
    else if (ndim == 2) volume = pi*radius*radius;
    else if (ndim == 3) volume = 4.0*onethird*pi*pow(radius,3);
    //if (mcloud > small_number && radius > small_number)
    //  rhofluid = mcloud / volume;
    rhofluid = mcloud / volume;


    // Record particle positions and initialise all other variables
    #pragma omp parallel for default(none)\
    shared(gammaone,mcloud,Npart,press,r,rhofluid,volume,temp0,mu_bar) private(i,k)
    for (i=0; i<hydro->Nhydro; i++) {
      Particle<ndim>& part = hydro->GetParticlePointer(i);
      for (k=0; k<ndim; k++) {
        part.r[k] = r[ndim*i + k];
        part.a[k] = (FLOAT) 0.0;
      }

      if (hydro->gas_eos == "isothermal") {
        part.u = temp0/gammaone/mu_bar;
      }
      else {
        part.u = press/rhofluid/gammaone;
      }
      part.m = mcloud / (FLOAT) Npart;
      part.h = hydro->h_fac*pow(part.m/rhofluid, invndim);
    }

    if (ndim == 3) {
      for (i=0; i<hydro->Nhydro; i++) {
        Particle<ndim>& part = hydro->GetParticlePointer(i);
        part.iorig = i;

        r_perp[0] = -part.r[1];
        r_perp[1] = part.r[0];
        r_perp[2] = (FLOAT) 0.0;
        norm      = pow(pow(r_perp[0],2) + pow(r_perp[1],2),0.5);
        r_perp[0] = r_perp[0]/norm;
        r_perp[1] = r_perp[1]/norm;
        r_perp[2] = r_perp[2]/norm;
        r_center  = pow(pow(part.r[0],2) + pow(part.r[1],2) + pow(part.r[2],2),0.5);
        for (k=0; k<ndim; k++) part.v[k] = r_perp[k]*r_center*omega*pc;
      }
    }

    sim->initial_h_provided = true;

    delete[] r;

  }
  //-----------------------------------------------------------------------------------------------
  else if (ic == "turbisothermsphere") {
    int i,k;                          // Particle and dimension counters
    FLOAT rcentre[ndim];              // Position of sphere centre
    FLOAT gpecloud;                   // ..
    FLOAT keturb;                     // ..
    FLOAT mp;                         // Mass of one particle
    FLOAT rho;                        // Fluid density
    FLOAT xmin;                       // ..
    FLOAT vfactor;                    // ..
    FLOAT *r;                         // Positions of all particles
    FLOAT *v;                         // Velocities of all particles
    FLOAT dxgrid;                     // ..
    FLOAT rmax[ndim];                 // ..
    FLOAT rmin[ndim];                 // ..
    DOUBLE *vfield;                   // ..

    // Local copies of important parameters
    int Npart            = simparams->intparams["Nsph"];
    int field_type       = simparams->intparams["field_type"];
    int gridsize         = simparams->intparams["gridsize"];
    FLOAT mcloud         = simparams->floatparams["mcloud"];
    FLOAT radius         = simparams->floatparams["radius"];
    FLOAT temp0          = simparams->floatparams["temp0"];
    FLOAT mu_bar         = simparams->floatparams["mu_bar"];
    FLOAT gammaone       = simparams->floatparams["gamma_eos"] - 1.0;
    FLOAT alpha_turb     = simparams->floatparams["alpha_turb"];
    FLOAT power_turb     = simparams->floatparams["power_turb"];
    string particle_dist = simparams->stringparams["particle_distribution"];

    debug2("[Ic::TurbIsothermSphere]");

#if !defined(FFTW_TURBULENCE)
    string message = "FFTW turbulence flag not set";
    ExceptionHandler::getIstance().raise(message);
#endif

    // Convert any parameters to code units
    mcloud /= simunits.m.outscale;
    radius /= simunits.r.outscale;
    temp0  /= simunits.temp.outscale;

    // Calculate gravitational potential energy of uniform cloud
    gpecloud = 0.6*mcloud*mcloud/radius;

    r = new FLOAT[ndim*Npart];
    v = new FLOAT[ndim*Npart];

    // Add a sphere of random particles with origin 'rcentre' and radius 'radius'
    for (k=0; k<ndim; k++) rcentre[k] = (FLOAT) 0.0;
    Ic<ndim>::Addr2Sphere(Npart, r, rcentre, radius);

    // Allocate local and main particle memory
    hydro->Nhydro = Npart;
    sim->AllocateParticleMemory();
    mp = mcloud / (FLOAT) Npart;
    rho = 3.0*mcloud / (4.0*pi*pow(radius,3));

    // Record particle properties in main memory
    for (i=0; i<Npart; i++) {
      Particle<ndim>& part = hydro->GetParticlePointer(i);
      for (k=0; k<ndim; k++) part.r[k] = r[ndim*i + k];
      for (k=0; k<ndim; k++) part.v[k] = v[ndim*i + k];
      part.m = mp;
      part.h = hydro->h_fac*pow(part.m/rho,invndim);
      part.u = temp0/gammaone/mu_bar;
    }

    sim->initial_h_provided = true;


    // Generate turbulent velocity field for given power spectrum slope
    vfield = new DOUBLE[ndim*gridsize*gridsize*gridsize];
    hydro->ComputeBoundingBox(rmax, rmin, hydro->Nhydro);
    xmin = 9.9e20;
    dxgrid = 0.0;
    for (k=0; k<ndim; k++) {
      dxgrid = max(dxgrid,(rmax[k] - rmin[k])/(FLOAT) (gridsize - 1));
      xmin = min(xmin, rmin[k]);
    }

    // Generate gridded velocity field
    Ic<ndim>::GenerateTurbulentVelocityField(field_type, gridsize, power_turb, vfield);

    // Now interpolate generated field onto particle positions
    Ic<ndim>::InterpolateVelocityField(hydro->Nhydro, gridsize, xmin, dxgrid, r, vfield, v);

    // Finally, copy velocities to main SPH particle array
    for (i=0; i<hydro->Nhydro; i++) {
      Particle<ndim>& part = hydro->GetParticlePointer(i);
      for (k=0; k<ndim; k++) part.v[k] = v[ndim*i + k];
    }

    // Change to COM frame of reference
    sim->SetComFrame();

    // Calculate total kinetic energy of turbulent velocity field
    keturb = 0.0;
    for (i=0; i<hydro->Nhydro; i++) {
      Particle<ndim>& part = hydro->GetParticlePointer(i);
      keturb += part.m*DotProduct(part.v, part.v, ndim);
    }
    keturb *= 0.5;

    vfactor = sqrt(alpha_turb*gpecloud/keturb);
    cout << "Scaling factor : " << vfactor << endl;

    // Now rescale velocities to give required turbulent energy in cloud
    for (i=0; i<hydro->Nhydro; i++) {
      Particle<ndim>& part = hydro->GetParticlePointer(i);
      for (k=0; k<ndim; k++) part.v[k] *= vfactor;
    }

    delete[] v;
    delete[] r;

  }
  //-----------------------------------------------------------------------------------------------

  return;
}



template class IsothermalSphereIc<1>;
template class IsothermalSphereIc<2>;
template class IsothermalSphereIc<3>;
