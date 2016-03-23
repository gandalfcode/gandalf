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
//  Silcc::Silcc
/// Set-up SILCC-type simulation initial conditions.
//=================================================================================================
template <int ndim>
SilccIc<ndim>::SilccIc(Simulation<ndim>* _sim, Hydrodynamics<ndim>* _hydro, FLOAT _invndim) :
  Ic<ndim>(_sim, _hydro, _invndim)
{
  // Create local copies of initial conditions parameters
  a_midplane   = simparams->floatparams["a_midplane"];
  h_midplane   = simparams->floatparams["h_midplane"];
  rho_midplane = simparams->floatparams["rho_midplane"];
  sigma_star   = simparams->floatparams["sigma_star"];
  z_d          = simparams->floatparams["z_d"];

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
  box_area  = simbox.boxsize[1]*simbox.boxsize[2];
  rho_star  = (FLOAT) 0.25*sigma_star/z_d;
  rho_a     = rho_midplane*exp(-a_midplane*a_midplane/h_midplane/h_midplane);
  m_exp     = (FLOAT) 0.5*sqrt(pi)*rho_midplane*h_midplane*erf(a_midplane/h_midplane)*box_area;
  m_uniform = rho_a*box_area*(simbox.boxmax[0] - a_midplane);
  m_box     = 2.0*(m_exp + m_uniform);

  // TESTING :
  //m_box = 1.2*0.5*box_area + 0.8*0.5*box_area;
  cout << "m_box : " << m_box << "   m_exp : " << m_exp << "   m_uniform : " << m_uniform << endl;
  cout << "rho_midplane : " << rho_midplane*simunits.rho.outscale << "   rho_star : " << rho_star*simunits.rho.outscale << endl;
  cout << "rho_a : " << rho_a*simunits.rho.outscale << "   " << "    box_area : " <<  box_area << endl;
  cout << "a_midplane : " << a_midplane << "   " << h_midplane << "   " << endl;

}



//=================================================================================================
//  Silcc::Generate
/// Set-up SILCC-type simulation initial conditions.
//=================================================================================================
template <int ndim>
void SilccIc<ndim>::Generate(void)
{
  // Only compile for 3-dimensional case
  //-----------------------------------------------------------------------------------------------
  //if (ndim == 3) {

    int i;                               // Particle counter
    int k;                               // Dimension counter
    FLOAT mp;                            // Mass of individual particle
    FLOAT z;                             // z-position of newly inserted particle

    // Create local copies of initial conditions parameters
    int Npart          = simparams->intparams["Nhydro"];
    FLOAT gammaone     = simparams->floatparams["gamma_eos"] - (FLOAT) 1.0;

    debug2("[SilccIc::Generate]");


    // Allocate local and main particle memory
    hydro->Nhydro = Npart;
    sim->AllocateParticleMemory();
    mp = m_box / (FLOAT) Npart;

    std::cout << "Nhydro : " << Npart << std::endl;


    // Record particle properties in main memory
    for (i=0; i<hydro->Nhydro; i++) {
      Particle<ndim>& part = hydro->GetParticlePointer(i);

      /*part.r[0] = simbox.boxmin[0] + simbox.boxsize[0]*sim->randnumb->floatrand();
      part.r[1] = simbox.boxmin[1] + simbox.boxsize[1]*sim->randnumb->floatrand();
      z = m_box*((FLOAT) 2.0*sim->randnumb->floatrand() - (FLOAT) 1.0);*/

      for (k=0; k<ndim; k++) part.r[k] = simbox.boxmin[k] + simbox.boxsize[k]*sim->randnumb->floatrand();
      for (k=0; k<ndim; k++) part.v[k] = (FLOAT) 0.0;
      part.m = mp;
      part.u = (FLOAT) 1.5;
      //part.h = hydro->h_fac*powf(mp/rho,invndim);
      part.itype = gas;

    }

  //}
  //-----------------------------------------------------------------------------------------------

  return;
}



//=================================================================================================
//  Silcc::GetValue
/// ...
//=================================================================================================
template <int ndim>
FLOAT SilccIc<ndim>::GetValue
 (std::string var,
  FLOAT r[ndim])
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

    if (fabs(r[0]) <= a_midplane) {
      return rho_midplane*exp(-r[0]*r[0]/h_midplane/h_midplane);
    }
    else {
      return rho_a;
    }
    //if (r[1] > 0.0) return 0.5;
    //else return 1.5;
    //return 1.0 + r[1]/1000.0;
  }
  else {
    std::cout << "Invalid string variable for Silcc::GetValue" << std::endl;
    return 0.0;
  }
}



template class SilccIc<1>;
template class SilccIc<2>;
template class SilccIc<3>;
