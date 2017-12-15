//=================================================================================================
//  PlaneParallelRadiation.cpp
//  ...
//
//  This file is part of GANDALF :
//  Graphical Astrophysics code for N-body Dynamics And Lagrangian Fluids
//  https://github.com/gandalfcode/gandalf
//  Contact : gandalfcode@gmail.com
//
//  Copyright (C) 2013  D. A. Hubber, G. Rosotti
//            (C) 2015  R. Wunsch
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


#include "Radiation.h"
#include "chealpix.h"


//=================================================================================================
//  PlaneParallelRadiation::PlaneParallelRadiation
/// Constructor for main PlaneParallelRadiation object
//=================================================================================================
template <int ndim, template<int> class ParticleType>
PlaneParallelRadiation<ndim,ParticleType>::PlaneParallelRadiation
 (Parameters *params, SmoothingKernel<ndim> *_kern, NeighbourSearch<ndim> *_neib):
  Radiation<ndim>()
{
  debug2("[PlaneParallelRadiation::PlaneParallelRadiation]");

  kern = _kern;
  neib = _neib;

  NLyC      = params->floatparams["NLyC"];
  arecomb   = params->floatparams["arecomb"];
  xmin      = params->floatparams["boxmin[0]"];

  FLOAT gammam1   = params->floatparams["gamma_eos"] - 1.0;
  FLOAT mu_ion    = params->floatparams["mu_ion"];
  FLOAT temp_ion  = params->floatparams["temp_ion"];
  uion = temp_ion/gammam1/mu_ion;

  maxIntegral = NLyC / arecomb;

  //tree = static_cast<OctTree<ndim,ParticleType,TreeCell>* > (neib->GetTree());
}



//=================================================================================================
//  PlaneParallelRadiation::~PlaneParallelRadiation
/// Destructor for PlaneParallelRadiation object
//=================================================================================================
template <int ndim, template<int> class ParticleType>
PlaneParallelRadiation<ndim,ParticleType>::~PlaneParallelRadiation()
{
}



//=================================================================================================
//  PlaneParallelRadiation::UpdateRadiationField
/// ...
//=================================================================================================
template <int ndim, template<int> class ParticleType>
void PlaneParallelRadiation<ndim,ParticleType>::UpdateRadiationField
 (int Nhydro,                          ///< [in] No. of hydro particles
  int Nnbody,                          ///< [in] No. of N-body particles
  int Nsink,                           ///< [in] No. of sink particles
  Particle<ndim> *part_gen,            ///< [in] Generic hydro particle data array
  NbodyParticle<ndim> **nbodydata,     ///< [in] N-body data array
  SinkParticle<ndim> *sinkdata)        ///< [in] Sink data array
{
  int numSteps = 0;
  int numIonised = 0;
  FLOAT dir[ndim];
  PlanarRay ray;
  ParticleType<ndim>* partdata = static_cast<ParticleType<ndim>* >(part_gen);

  debug2("[PlaneParallelRadiation::UpdateRadiationField]");
  cout << "UPDATING PLANE PARALLEL RADIATION" << std::endl;

  // Find initial position of ray based on the size of the tree root cell
  for (int k=0; k<ndim; k++) ray.r[k] = (FLOAT) 0.5;
  ray.r[0] = xmin;

  // Initial ray points in positive x-direction
  for (int k=0; k<ndim; k++) dir[k] = (FLOAT) 0.0;
  dir[0] = (FLOAT) 1.0;

  // For now, use minimum smoothing length of particles as constant step-size
  FLOAT hmin = big_number;
  FLOAT hmax = big_number;
  for (int i=0; i<Nhydro; i++) {
    hmin = min(hmin, partdata[i].h);
    hmax = max(hmax, partdata[i].h);
  }
  FLOAT step = (FLOAT) 0.2*hmin;
  FLOAT invhmaxsqd = (FLOAT) 1.0 / hmax / hmax;

  // Compute density at initial place
  FLOAT rho0 = (FLOAT) 0.0;
  for (int i=0; i<Nhydro; i++) {
    FLOAT dr[ndim];
    for (int k=0; k<ndim; k++) dr[k] = ray.r[k] - partdata[i].r[k];
    FLOAT drsqd = DotProduct(dr, dr, ndim);
    FLOAT ssqd  = drsqd*invhmaxsqd;
    rho0       += partdata[i].m*kern->w0_s2(ssqd);
  }

  // Ray-trace through the computational domain computing the ionisation integral
  //-----------------------------------------------------------------------------------------------
  do {

    for (int k=0; k<ndim; k++) ray.r[k] += step*dir[k];

    // Compute density at new point
    FLOAT rho1 = (FLOAT) 0.0;
    for (int i=0; i<Nhydro; i++) {
      FLOAT dr[ndim];
      for (int k=0; k<ndim; k++) dr[k] = ray.r[k] - partdata[i].r[k];
      FLOAT drsqd = DotProduct(dr, dr, ndim);
      FLOAT ssqd  = drsqd*invhmaxsqd;
      rho1       += partdata[i].m*kern->w0_s2(ssqd);
    }

    FLOAT dIntegral = (FLOAT) 0.5*(rho0*rho0 + rho1*rho1)*step;
    ray.rayIntegral += dIntegral;
    numSteps++;

  } while (ray.rayIntegral < maxIntegral);
  //-----------------------------------------------------------------------------------------------

  // Set ionisation fractions of particles based on position relative to ionisation front
  for (int i=0; i<Nhydro; i++) {
    if (partdata[i].r[0] > ray.r[0]) {
      partdata[i].ionstate = 0;
    }
    else {
      partdata[i].ionstate = 1;
      partdata[i].u = uion;
      numIonised++;
      //std::cout << "FOUND IONISED PARTICLE : " << i << "  " << partdata[i].r[0] << std::endl;
    }
  }

  std::cout << "NO. OF STEPS : " << numSteps << "  " << ray.rayIntegral << "   " << maxIntegral << std::endl;
  std::cout << "NO. IONISED  : " << numIonised << "    r : " << ray.r[0] << std::endl;

  return;
}



template class PlaneParallelRadiation<3, GradhSphParticle>;
template class PlaneParallelRadiation<2, GradhSphParticle>;
template class PlaneParallelRadiation<1, GradhSphParticle>;

template class PlaneParallelRadiation<3, MeshlessFVParticle>;
template class PlaneParallelRadiation<2, MeshlessFVParticle>;
template class PlaneParallelRadiation<1, MeshlessFVParticle>;
