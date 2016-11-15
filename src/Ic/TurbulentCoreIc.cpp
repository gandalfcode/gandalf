//=================================================================================================
//  TurbulentCoreIc.cpp
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
//  TurbulentCoreIc::TurbulentCoreIc
/// Set-up SILCC-type simulation initial conditions.
//=================================================================================================
template <int ndim>
TurbulentCoreIc<ndim>::TurbulentCoreIc(Simulation<ndim>* _sim, Hydrodynamics<ndim>* _hydro, FLOAT _invndim) :
  Ic<ndim>(_sim, _hydro, _invndim)
{
  // Some sanity checking to ensure correct dimensionality is used
  if (simparams->intparams["ndim"] != 3) {
    ExceptionHandler::getIstance().raise("Turbulent core sim only runs in 3D");
  }
  if (simparams->intparams["dimensionless"] != 0) {
    ExceptionHandler::getIstance().raise("dimensionless units not permitted");
  }
#if !defined(FFTW_TURBULENCE)
  ExceptionHandler::getIstance().raise("FFTW turbulence flag not set");
#endif
}



//=================================================================================================
//  Silcc::Generate
/// Set-up SILCC-type simulation initial conditions.
//=================================================================================================
template <int ndim>
void TurbulentCoreIc<ndim>::Generate(void)
{
  // Only compile for 3-dimensional case
  //-----------------------------------------------------------------------------------------------
  if (ndim == 3) {

    int i;                               // Particle counter
    int k;                               // Dimension counter
    int Nsphere;                         // Actual number of particles in sphere
    FLOAT gpecloud;                      // Total grav. potential energy of entire cloud
    FLOAT keturb;                        // Total turbulent kinetic energy of entire cloud
    FLOAT mp;                            // Mass of one particle
    FLOAT rcentre[ndim];                 // Position of sphere centre
    FLOAT rho;                           // Fluid density
    FLOAT xmin;                          // Minimum coordinate value
    FLOAT vfactor;                       // Velocity scaling factor (to scale correct alpha_turb)
    FLOAT *r;                            // Positions of all particles
    FLOAT *v;                            // Velocities of all particles
    FLOAT dxgrid;                        // Grid spacing
    FLOAT rmax[ndim];                    // Maximum size of bounding box
    FLOAT rmin[ndim];                    // Minimum size of bounding box
    DOUBLE *vfield;                      // Table with turbulent velocity field from

    // Create local copies of initial conditions parameters
    int field_type   = simparams->intparams["field_type"];
    int gridsize     = simparams->intparams["gridsize"];
    int Npart        = simparams->intparams["Nhydro"];
    FLOAT alpha_turb = simparams->floatparams["alpha_turb"];
    FLOAT gammaone   = simparams->floatparams["gamma_eos"] - 1.0;
    FLOAT mcloud     = simparams->floatparams["mcloud"];
    FLOAT mu_bar     = simparams->floatparams["mu_bar"];
    FLOAT power_turb = simparams->floatparams["power_turb"];
    FLOAT radius     = simparams->floatparams["radius"];
    FLOAT temp0      = simparams->floatparams["temp0"];
    string particle_dist = simparams->stringparams["particle_distribution"];

    debug2("[TurbulentCoreIc::Generate]");

    // Convert any parameters to code units
    mcloud /= simunits.m.outscale;
    radius /= simunits.r.outscale;
    temp0  /= simunits.temp.outscale;

    // Calculate gravitational potential energy of uniform density spherical cloud
    gpecloud = (FLOAT) 0.6*mcloud*mcloud/radius;

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
      assert(Nsphere <= Npart);
      if (Nsphere != Npart)
        cout << "Warning! Unable to converge to required "
             << "no. of ptcls due to lattice symmetry" << endl;
      Npart = Nsphere;
    }
    else {
      string message = "Invalid particle distribution option";
      ExceptionHandler::getIstance().raise(message);
    }

    // Allocate local and main particle memory
    hydro->Nhydro = Npart;
    sim->AllocateParticleMemory();
    mp = mcloud / (FLOAT) Npart;
    rho = (FLOAT) 3.0*mcloud / ((FLOAT) 4.0*pi*pow(radius,3));


    // Record particle properties in main memory
    for (i=0; i<hydro->Nhydro; i++) {
      Particle<ndim>& part = hydro->GetParticlePointer(i);
      for (k=0; k<ndim; k++) part.r[k] = r[ndim*i + k];
      for (k=0; k<ndim; k++) part.v[k] = 0.0;
      part.m = mp;
      part.h = hydro->h_fac*powf(mp/rho,invndim);
      part.u = temp0/gammaone/mu_bar;
      part.ptype = gas_type;
    }

    sim->initial_h_provided = true;


    // Generate turbulent velocity field for given power spectrum slope
    vfield = new DOUBLE[ndim*gridsize*gridsize*gridsize];

    // Calculate bounding box of SPH smoothing kernels
    for (k=0; k<ndim; k++) rmin[k] = big_number;
    for (k=0; k<ndim; k++) rmax[k] = -big_number;
    for (i=0; i<Npart; i++) {
      Particle<ndim>& part = hydro->GetParticlePointer(i);
      for (k=0; k<ndim; k++) rmin[k] = min(rmin[k], part.r[k] - hydro->kernrange*part.h);
      for (k=0; k<ndim; k++) rmax[k] = max(rmax[k], part.r[k] + hydro->kernrange*part.h);
    }

    xmin = (FLOAT) 9.9e20;
    dxgrid = (FLOAT) 0.0;
    for (k=0; k<ndim; k++) {
      xmin = min(xmin, rmin[k]);
      dxgrid = max(dxgrid, (rmax[k] - rmin[k])/(FLOAT) (gridsize - 1));
      //xmin = min(xmin,rmin[k]);
    }
    dxgrid = max(dxgrid, (FLOAT) 2.0*fabs(xmin)/(FLOAT) (gridsize - 1));

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
    keturb = (FLOAT) 0.0;
    for (i=0; i<hydro->Nhydro; i++) {
      Particle<ndim>& part = hydro->GetParticlePointer(i);
      keturb += part.m*DotProduct(part.v, part.v, ndim);
    }
    keturb *= (FLOAT) 0.5;

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




template class TurbulentCoreIc<1>;
template class TurbulentCoreIc<2>;
template class TurbulentCoreIc<3>;
