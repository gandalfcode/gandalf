//=================================================================================================
//  SedovBlastwaveIc.cpp
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
//  SedovBlastwaveIc::SedovBlastwaveIc
/// Set-up SILCC-type simulation initial conditions.
//=================================================================================================
template <int ndim>
SedovBlastwaveIc<ndim>::SedovBlastwaveIc(Simulation<ndim>* _sim, FLOAT _invndim) :
  Ic<ndim>(_sim, _invndim)
{
  // Some sanity checking to ensure dimensionless units are used
  if (simparams->intparams["dimensionless"] == 0) {
    ExceptionHandler::getIstance().raise("dimensionless units not permitted");
  }
}



//=================================================================================================
//  Silcc::Generate
/// Set-up SILCC-type simulation initial conditions.
//=================================================================================================
template <int ndim>
void SedovBlastwaveIc<ndim>::Generate(void)
{
  int i;                                   // Particle counter
  int k;                                   // Dimension counter
  int Nbox;                                // No. of particles in box
  int Nlattice[3];                         // Lattice size
  FLOAT drmag;                             // Distance
  FLOAT drsqd;                             // Distance squared
  FLOAT mbox;                              // Total mass inside simulation box
  FLOAT umax;                              // Maximum u of all particles
  FLOAT utot;                              // Total internal energy
  FLOAT volume;                            // Volume of box
  FLOAT *r;                                // Positions of all particles

  // Create local copies of initial conditions parameters
  Nlattice[0] = simparams->intparams["Nlattice1[0]"];
  Nlattice[1] = simparams->intparams["Nlattice1[1]"];
  Nlattice[2] = simparams->intparams["Nlattice1[2]"];
  dusty_shock = simparams->stringparams["dust_forces"] != "none";
  smooth_ic   = simparams->intparams["smooth_ic"];
  rhofluid    = simparams->floatparams["rhofluid1"];
  kefrac      = simparams->floatparams["kefrac"];
  string particle_dist = simparams->stringparams["particle_distribution"];

  debug2("[SedovBlastwaveIc::Generate]");


  // Compute size and range of fluid bounding boxes
  //-----------------------------------------------------------------------------------------------
  if (ndim == 1) {
    volume = icBox.max[0] - icBox.min[0];
    Nbox = Nlattice[0];
  }
  else if (ndim == 2) {
    volume = (icBox.max[0] - icBox.min[0])*(icBox.max[1] - icBox.min[1]);
    Nbox = Nlattice[0]*Nlattice[1];
  }
  else if (ndim == 3) {
    volume = (icBox.max[0] - icBox.min[0])*
      (icBox.max[1] - icBox.min[1])*(icBox.max[2] - icBox.min[2]);
    Nbox = Nlattice[0]*Nlattice[1]*Nlattice[2];
  }
  mbox  = volume*rhofluid;
  r_hot = hydro->h_fac*hydro->kernrange*icBox.size[0]/Nlattice[0];
  std::cout << "r_hot : " << r_hot << std::endl;
  //exit(0);

  // Allocate local and main particle memory
  hydro->Nhydro = Nbox;
  if (dusty_shock) hydro->Nhydro *= 2;


  sim->AllocateParticleMemory();
  r = new FLOAT[ndim*Nbox];

  // Add a cube of random particles defined by the simulation bounding box and
  // depending on the chosen particle distribution
  if (particle_dist == "random") {
    Ic<ndim>::AddRandomBox(Nbox, icBox, r, sim->randnumb);
  }
  else if (particle_dist == "cubic_lattice") {
    Ic<ndim>::AddCubicLattice(Nbox, Nlattice, icBox, true, r);
  }
  else if (particle_dist == "hexagonal_lattice") {
    Ic<ndim>::AddHexagonalLattice(Nbox, Nlattice, icBox, true, r);
  }
  else {
    string message = "Invalid particle distribution option";
    ExceptionHandler::getIstance().raise(message);
  }

  // Record positions in main memory
  for (i=0; i<Nbox; i++) {
    Particle<ndim>& part = hydro->GetParticlePointer(i);
    for (k=0; k<ndim; k++) part.r[k] = r[ndim*i + k];
    for (k=0; k<ndim; k++) part.v[k] = (FLOAT) 0.0;
    part.m = mbox/(FLOAT) Nbox;
    part.h = hydro->h_fac*pow(part.m/rhofluid,invndim);
    part.u = small_number;
    //part.v[0] = 0.1*part.r[0] + 0.2*part.r[1] - 0.7*part.r[2];
  }

  // Set initial smoothing lengths and create initial ghost particles
  //-----------------------------------------------------------------------------------------------
  hydro->Nghost = 0;
  hydro->Ntot = hydro->Nhydro;
  for (i=0; i<hydro->Nhydro; i++) hydro->GetParticlePointer(i).flags.set(active);

  sim->initial_h_provided = true;
  sim->rebuild_tree = true;


  delete[] r;

  return;
}



//=================================================================================================
//  Silcc::SetParticleProperties
/// Sets the properties of all particles once their positions have been allocated.
//=================================================================================================
template <int ndim>
void SedovBlastwaveIc<ndim>::SetParticleProperties()
{
  int i;                                   // Particle counter
  int k;                                   // Dimension counter
  int Nbox = hydro->Nhydro;                 // No. of particles in box
  int Ncold = 0;                               // No. of cold particles
  int Nhot = 0;                                // No. of hot particles
  FLOAT drmag;                             // Distance
  FLOAT drsqd;                             // Distance squared
  FLOAT umax = (FLOAT) 0.0;                // Maximum u of all particles
  FLOAT utot = (FLOAT) 0.0;                // Total internal energy
  FLOAT ufrac = max((FLOAT) 0.0,(FLOAT) 1.0 - kefrac);  // Internal energy fraction
  int *hotlist = new int[Nbox];            // List of 'hot' particles

  std::cout << "r_hot : " << r_hot << std::endl;

  // Now calculate which particles are hot
  //-----------------------------------------------------------------------------------------------
  for (i=0; i<Nbox; i++) {
    Particle<ndim>& part = hydro->GetParticlePointer(i);
    drsqd = DotProduct(part.r,part.r,ndim);
    hotlist[i] = 0;
    Ncold++;

    if (drsqd < r_hot*r_hot) {
      hotlist[i] = 1;
      if (smooth_ic == 1) part.u = part.m*hydro->kernp->w0(hydro->kernp->kernrange*sqrt(drsqd)/r_hot);
      else part.u = part.m;
      utot += part.u;
      umax = max(umax,part.u);
      Nhot++;
    }
    else {
      hotlist[i] = 0;
      Ncold++;
    }
  }

  // Normalise the energies
  //-----------------------------------------------------------------------------------------------
  for (i=0; i<Nbox; i++) {
    Particle<ndim>& part = hydro->GetParticlePointer(i);
    if (hotlist[i] == 1) {
      drmag = sqrt(DotProduct(part.r,part.r,ndim));
      part.u = part.u/utot/part.m;
      for (k=0; k<ndim; k++) part.v[k] = sqrt(2.0*kefrac*part.u)*part.r[k]/(drmag + small_number);
      part.u = ufrac*part.u;
    }
    else {
      part.u = 1.0e-6/part.m;
    }
  }


  // Add a slightly offset dust lattice
  if (dusty_shock){
    FLOAT d2g = simparams->floatparams["dust_mass_factor"] ;
    for (int j = 0; j < Nbox; ++j){
      Particle<ndim>& pg = hydro->GetParticlePointer(j) ;
      Particle<ndim>& pd = hydro->GetParticlePointer(j+Nbox) ;
      pd = pg ;
      pg.ptype = gas_type ;
      pd.ptype = dust_type ;
      pd.r[0] += 0.01 * pd.h ;
      pd.m *= d2g ;
      pd.u = 0 ;
      pd.h_dust = pd.h ;
    }
  }

  sim->initial_h_provided = true;

  delete[] hotlist;
}



//=================================================================================================
//  SedovBlastwaveIc::GetDensity
/// Returns the value of the density at the given position.
//=================================================================================================
template <int ndim>
FLOAT SedovBlastwaveIc<ndim>::GetDensity
 (const FLOAT r[ndim],
  const int ptype) const
{
  return rhofluid;
}



//=================================================================================================
//  Silcc::GetParticleRegularizer
/// Return the regularizer based upon the density.
//=================================================================================================
template <int ndim>
Regularization::RegularizerFunction<ndim>* SedovBlastwaveIc<ndim>::GetParticleRegularizer() const {
  using Regularization::DefaultRegularizerFunction;
  return new DefaultRegularizerFunction<ndim,SedovBlastwaveIc<ndim> >(hydro->kernp, simparams, this);
}


template class SedovBlastwaveIc<1>;
template class SedovBlastwaveIc<2>;
template class SedovBlastwaveIc<3>;
