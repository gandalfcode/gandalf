//=================================================================================================
//  SilccIc.cpp
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
//=================================================================================================


#include <fstream>
#include <sstream>
#include "Precision.h"
#include "Debug.h"
#include "Ic.h"
using namespace std;



//=================================================================================================
//  Ic::Silcc
/// Set-up SILCC-type simulation initial conditions.
//=================================================================================================
template <int ndim>
void SilccIc<ndim>::Generate(void)
{
  // Only compile for 3-dimensional case
  //-----------------------------------------------------------------------------------------------
  if (ndim == 3) {

    int i;                               // Particle counter
    int k;                               // Dimension counter
    FLOAT box_area;                      // Area of x-y plane of simulation box
    FLOAT m_box;                         // Total gas mass in box
    FLOAT m_exp;                         // Total gas mass in exponential profile region
    FLOAT m_uniform;                     // Total gas mass in uniform density region
    FLOAT mp;                            // Mass of individual particle
    FLOAT rho_a;                         // Density at edge of exponential midplane profile
    FLOAT rho_star;                      // Stellar density at the midplane
    FLOAT z;                             // z-position of newly inserted particle

    // Create local copies of initial conditions parameters
    int Npart          = simparams->intparams["Nhydro"];
    FLOAT a_midplane   = simparams->floatparams["a_midplane"];
    FLOAT gammaone     = simparams->floatparams["gamma_eos"] - (FLOAT) 1.0;
    FLOAT h_midplane   = simparams->floatparams["h_midplane"];
    FLOAT rho_midplane = simparams->floatparams["rho_midplane"];
    FLOAT sigma_star   = simparams->floatparams["sigma_star"];
    FLOAT z_d          = simparams->floatparams["z_s"];

    debug2("[Ic::Silcc]");

    // Some sanity checking to ensure correct units are used for these ICs
    if (simparams->stringparams["routunit"] != "pc") {
      ExceptionHandler::getIstance().raise("r unit not set to pc");
    }
    if (simparams->stringparams["sigmaoutunit"] != "m_sun_pc2") {
      ExceptionHandler::getIstance().raise("sigma unit not set to m_sun_pc2");
    }

    // Convert any parameters to code units
    a_midplane   /= simunits.r.outscale;
    h_midplane   /= simunits.r.outscale;
    rho_midplane /= simunits.rho.outscale;
    sigma_star   /= simunits.sigma.outscale;
    z_d          /= simunits.r.outscale;

    // Compute total mass of particles in simulation box by integrating in the z-direction
    box_area  = simbox.boxsize[0]*simbox.boxsize[1];
    rho_star  = (FLOAT) 0.25*sigma_star/z_d;
    rho_a     = rho_midplane*exp(-a_midplane*a_midplane/h_midplane/h_midplane);
    m_exp     = (FLOAT) 0.5*sqrt(pi)*rho_midplane*h_midplane*erf(a_midplane/h_midplane)*box_area;
    m_uniform = rho_a*box_area*(simbox.boxmax[0] - a_midplane);
    m_box     = m_exp + m_uniform;

    // Allocate local and main particle memory
    hydro->Nhydro = Npart;
    sim->AllocateParticleMemory();
    mp = m_box / (FLOAT) Npart;


    // Record particle properties in main memory
    for (i=0; i<hydro->Nhydro; i++) {
      Particle<ndim>& part = hydro->GetParticlePointer(i);

      part.r[0] = simbox.boxmin[0] + simbox.boxsize[0]*sim->randnumb->floatrand();
      part.r[1] = simbox.boxmin[1] + simbox.boxsize[1]*sim->randnumb->floatrand();
      z = m_box*((FLOAT) 2.0*sim->randnumb->floatrand() - (FLOAT) 1.0);


      for (k=0; k<ndim; k++) part.v[k] = 0.0;
      part.m = mp;
      //part.h = hydro->h_fac*powf(mp/rho,invndim);
      part.itype = gas;
    }

    sim->initial_h_provided = true;

  }
  //-----------------------------------------------------------------------------------------------

  return;
}



template class SilccIc<1>;
template class SilccIc<2>;
template class SilccIc<3>;
