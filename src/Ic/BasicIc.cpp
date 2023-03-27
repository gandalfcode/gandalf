//=================================================================================================
//  BasicIc.cpp
//  Class for generating simple initial conditions for testing the IC regularisation algorithms.
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
//  BasicIc::BasicIc
/// Set-up BasicIc-type simulation initial conditions.
//=================================================================================================
template <int ndim>
BasicIc<ndim>::BasicIc(Simulation<ndim>* _sim, FLOAT _invndim) :
  Ic<ndim>(_sim, _invndim)
{
  if (simparams->intparams["dimensionless"] == 0) {
    ExceptionHandler::getIstance().raise("dimensionless units required");
  }

  amp = simparams->floatparams["amp"];
  rho = simparams->floatparams["rhofluid1"];
  lambda = simparams->floatparams["boxmax[0]"] - simparams->floatparams["boxmin[0]"];

}



//=================================================================================================
//  BasicIc::Generate
/// Set-up BasicIc-type simulation initial conditions.
//=================================================================================================
template <int ndim>
void BasicIc<ndim>::Generate(void)
{
  string ic = simparams->stringparams["ic"];

  debug2("[BasicIc::Generate");

  //-----------------------------------------------------------------------------------------------
  if (ic == "basic_sine") {
    FLOAT volume;
    Box<ndim> box;

    // Local copy of important parameters
    int Npart = simparams->intparams["Nhydro"];
    for (int k=0; k<ndim; k++) box.min[k] = sim->icBox.min[k];
    for (int k=0; k<ndim; k++) box.max[k] = sim->icBox.max[k];

    FLOAT *r = new FLOAT[ndim*Npart];
    Ic<ndim>::AddMonteCarloDensityField(Npart, gas_type, box, r, sim->randnumb);

    // Compute volume and number of particles inside box
    if (ndim == 1) {
      volume = icBox.max[0] - icBox.min[0];
    }
    else if (ndim == 2) {
      volume = (icBox.max[0] - icBox.min[0])*(icBox.max[1] - icBox.min[1]);
    }
    else if (ndim == 3) {
      volume = (icBox.max[0] - icBox.min[0])*
        (icBox.max[1] - icBox.min[1])*(icBox.max[2] - icBox.min[2]);
    }
    FLOAT mp = rho*volume/(FLOAT) Npart;

    // Allocate global and local memory for all particles
    hydro->Nhydro = Npart;
    sim->AllocateParticleMemory();

    // Copy positions to main array and initialise all other variables
    for (int i=0; i<hydro->Nhydro; i++) {
      Particle<ndim>& part = hydro->GetParticlePointer(i);
      for (int k=0; k<ndim; k++) part.r[k] = r[ndim*i + k];
      part.m = mp;
    }

    // Now set all other particle properties
    SetParticleProperties();

    delete[] r;
  }
  //-----------------------------------------------------------------------------------------------
  else if (ic == "sphere") {
  }
  //-----------------------------------------------------------------------------------------------

  return;
}



//=================================================================================================
//  BasicIc::GetDensity
/// Returns the value of the density at the given position.
//=================================================================================================
template <int ndim>
FLOAT BasicIc<ndim>::GetDensity
 (const FLOAT r[ndim],
  const int ptype) const
{
  return rho*(1.0 + amp*sin(twopi*r[0]/lambda));
}



//=================================================================================================
//  BasicIc::SetParticleProperties
/// Sets the properties of all particles once their positions have been allocated.
//=================================================================================================
template <int ndim>
void BasicIc<ndim>::SetParticleProperties()
{
  // Copy positions to main array and initialise all other variables
  for (int i=0; i<hydro->Nhydro; i++) {
    Particle<ndim>& part = hydro->GetParticlePointer(i);
    for (int k=0; k<ndim; k++) {
      part.v[k] = (FLOAT) 0.0;
      part.a[k] = (FLOAT) 0.0;
    }
    part.u     = 1.0;
    part.iorig = i;
    part.ptype = gas_type;
  }

  return;
}



//=================================================================================================
//  BasicIc::GetParticleRegularizer
/// Return the regularizer based upon the density.
//=================================================================================================
template <int ndim>
Regularization::RegularizerFunction<ndim>* BasicIc<ndim>::GetParticleRegularizer() const {
  using Regularization::DefaultRegularizerFunction;
  return new DefaultRegularizerFunction<ndim,BasicIc<ndim> >(hydro->kernp, simparams, this);
}



template class BasicIc<1>;
template class BasicIc<2>;
template class BasicIc<3>;
