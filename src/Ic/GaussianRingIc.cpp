//=================================================================================================
//  GaussianRingIc.cpp
//  Class for generating initial conditions ...
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
//  GaussianRingIc::GaussianRingIc
/// ...
//=================================================================================================
template <int ndim>
GaussianRingIc<ndim>::GaussianRingIc(Simulation<ndim>* _sim, FLOAT _invndim) :
  Ic<ndim>(_sim, _invndim)
{
  // Some sanity checking to ensure correct dimensionality and (dimensionless) units are used
  if (simparams->intparams["ndim"] != 2) {
    ExceptionHandler::getIstance().raise("Gaussian ring only functions in 2D");
  }
  if (simparams->intparams["dimensionless"] != 1) {
    ExceptionHandler::getIstance().raise("r unit not set to pc");
  }
}



//=================================================================================================
//  GaussianRingIc::Generate
/// ...
//=================================================================================================
template <int ndim>
void GaussianRingIc<ndim>::Generate(void)
{
  // Only compile for 3-dimensional case
  //-----------------------------------------------------------------------------------------------
  if (ndim == 2) {

    debug2("[GaussianRing::Generate]");

    //int Nhydro = simparams->intparams["Nhydro"];
    int Nhydro = 13188*2; //I've hard-coded it to reproduce Murray 1996
    const FLOAT temp0     = simparams->floatparams["temp0"];
    const FLOAT mu_bar    = simparams->floatparams["mu_bar"];
    //	const FLOAT alpha = simparams->floatparams["alpha_visc"];
    const FLOAT c_s = sqrt(temp0/mu_bar);
    cout << "sound speed: " << c_s << endl;
    //	const FLOAT nu = 1./8. * alpha * c_s * 0.01;

    // Parameters of the gaussian
    const FLOAT rcentre = 0.85;
    const FLOAT width = 0.025;
    const FLOAT inner_edge = 0.80;
    const FLOAT outer_edge = 0.90;
    const int nrings = 21;

    const int Nperring = Nhydro/nrings;
    Nhydro=nrings*Nperring;

    // Allocate memory
    hydro->Nhydro = Nhydro;
    sim->nbody->Nstar=1;
    sim->AllocateParticleMemory();

    // Set up the star
    StarParticle<ndim>& star = sim->nbody->stardata[0];
    star.m=1;

    // Set up the particles
    for (int i=0; i<Nhydro; i++) {

      Particle<ndim>& part = hydro->GetParticlePointer(i);

      // Position
      const int iring = i/Nperring;
      const FLOAT ring_spacing = (outer_edge-inner_edge)/(nrings-1);
      const FLOAT r = inner_edge + iring*ring_spacing;
      const FLOAT phi_spacing = 2*M_PI/Nperring;
      const FLOAT phi = phi_spacing*(i-iring*Nperring);

      part.r[0] = r*cos(phi);
      part.r[1] = r*sin(phi);

      // Velocity
      const FLOAT vphi = 1./sqrt(r);
      const FLOAT vr=0.0;
      //const FLOAT vr = -3*nu/2/r + 6*nu/width/width * (r-rcentre);

      part.v[0] = vr*cos(phi)-vphi*sin(phi);
      part.v[1] = vr*sin(phi)+vphi*cos(phi);


      // Density (and mass)
      const FLOAT sigma = exp (-pow((r-rcentre)/width,2));
      part.m = 0.01/Nhydro*sigma;

    }

    sim->initial_h_provided = false;
  }
  //-----------------------------------------------------------------------------------------------

  return;
}



template class GaussianRingIc<1>;
template class GaussianRingIc<2>;
template class GaussianRingIc<3>;
