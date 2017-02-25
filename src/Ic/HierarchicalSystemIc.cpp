//=================================================================================================
//  HierarchicalSystemIc.cpp
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
//  HierarchicalSystemIc::HierarchicalSystemIc
/// Set-up SILCC-type simulation initial conditions.
//=================================================================================================
template <int ndim>
HierarchicalSystemIc<ndim>::HierarchicalSystemIc(Simulation<ndim>* _sim, FLOAT _invndim) :
  Ic<ndim>(_sim, _invndim)
{
  // Some sanity checking to ensure correct dimensionality is used
  if (simparams->intparams["ndim"] == 1) {
    ExceptionHandler::getIstance().raise("N-body simulations do not run in 1d");
  }
}



//=================================================================================================
//  Silcc::Generate
/// Set-up SILCC-type simulation initial conditions.
//=================================================================================================
template <int ndim>
void HierarchicalSystemIc<ndim>::Generate(void)
{
  FLOAT rbinary[ndim];                // Position of binary COM
  FLOAT vbinary[ndim];                // Velocity of binary COM
  Nbody<ndim>* nbody = sim->nbody;

  string ic = simparams->stringparams["ic"];
  for (int k=0; k<ndim; k++) rbinary[k] = (FLOAT) 0.0;
  for (int k=0; k<ndim; k++) vbinary[k] = (FLOAT) 0.0;

  debug2("[HierarchicalSystemIc::Generate");

  //-----------------------------------------------------------------------------------------------
  if (ic == "binary") {

    // Binary star parameters
    FLOAT abin     = simparams->floatparams["abin"];
    FLOAT ebin     = simparams->floatparams["ebin"];
    FLOAT m1       = simparams->floatparams["m1"];
    FLOAT m2       = simparams->floatparams["m2"];
    FLOAT phirot   = simparams->floatparams["phirot"];
    FLOAT thetarot = simparams->floatparams["thetarot"];
    FLOAT psirot   = simparams->floatparams["psirot"];

    // Allocate local and main particle memory
    nbody->Nstar = 2;
    sim->AllocateParticleMemory();

    // Add binary star
    Ic<ndim>::AddBinaryStar(abin, ebin, m1, m2, (FLOAT) 0.01, (FLOAT) 0.01, phirot, thetarot, psirot,
                            (FLOAT) 0.0, rbinary, vbinary, nbody->stardata[0], nbody->stardata[1]);

  }
  //-----------------------------------------------------------------------------------------------
  else if (ic == "triple") {

    NbodyParticle<ndim> b1;              // Inner binary COM particle

    // Triple star parameters
    FLOAT abin1    = simparams->floatparams["abin"];
    FLOAT abin2    = simparams->floatparams["abin2"];
    FLOAT ebin1    = simparams->floatparams["ebin"];
    FLOAT ebin2    = simparams->floatparams["ebin2"];
    FLOAT m1       = simparams->floatparams["m1"];
    FLOAT m2       = simparams->floatparams["m2"];
    FLOAT m3       = simparams->floatparams["m3"];
    FLOAT phirot   = simparams->floatparams["phirot"];
    FLOAT psirot   = simparams->floatparams["psirot"];
    FLOAT thetarot = simparams->floatparams["thetarot"];

    // Allocate local and main particle memory
    nbody->Nstar = 3;
    sim->AllocateParticleMemory();

    // Compute main binary orbit
    Ic<ndim>::AddBinaryStar(abin1, ebin1, m1+m2, m3, (FLOAT) 0.0001, (FLOAT) 0.0001, phirot, thetarot,
                            psirot, (FLOAT) 0.0, rbinary, vbinary, b1, nbody->stardata[2]);

    // Now compute both components
    Ic<ndim>::AddBinaryStar(abin2, ebin2, m1, m2, (FLOAT) 0.0001, (FLOAT) 0.0001, phirot, thetarot,
                            psirot, (FLOAT) 0.0, b1.r, b1.v, nbody->stardata[0], nbody->stardata[1]);


  }
  //-----------------------------------------------------------------------------------------------
  else if (ic == "quadruple") {

    NbodyParticle<ndim> b1;              // Star/binary 1
    NbodyParticle<ndim> b2;              // Star/binary 2

    // Quadruple star parameters
    FLOAT abin1    = simparams->floatparams["abin"];
    FLOAT abin2    = simparams->floatparams["abin2"];
    FLOAT ebin1    = simparams->floatparams["ebin"];
    FLOAT ebin2    = simparams->floatparams["ebin2"];
    FLOAT m1       = simparams->floatparams["m1"];
    FLOAT m2       = simparams->floatparams["m2"];
    FLOAT m3       = simparams->floatparams["m3"];
    FLOAT m4       = simparams->floatparams["m4"];
    FLOAT phirot   = simparams->floatparams["phirot"];
    FLOAT psirot   = simparams->floatparams["psirot"];
    FLOAT thetarot = simparams->floatparams["thetarot"];

    // Allocate local and main particle memory
    nbody->Nstar = 4;
    sim->AllocateParticleMemory();

    // Compute main binary orbit
    Ic<ndim>::AddBinaryStar(abin1, ebin1, m1+m2, m3+m4, (FLOAT) 0.01, (FLOAT) 0.01, phirot,
                            thetarot, psirot, (FLOAT) 0.0, rbinary, vbinary, b1, b2);

    // Now compute components of both inner binaries
    Ic<ndim>::AddBinaryStar(abin2, ebin2, m1, m2, (FLOAT) 0.0001, (FLOAT) 0.0001, phirot, thetarot,
                            psirot, (FLOAT) 0.0, b1.r, b1.v, nbody->stardata[0], nbody->stardata[1]);
    Ic<ndim>::AddBinaryStar(abin2, ebin2, m3, m4, (FLOAT) 0.0001, (FLOAT) 0.0001, phirot, thetarot,
                            psirot, (FLOAT) 0.0, b2.r, b2.v, nbody->stardata[2], nbody->stardata[3]);

  }
  //-----------------------------------------------------------------------------------------------

  return;
}




template class HierarchicalSystemIc<1>;
template class HierarchicalSystemIc<2>;
template class HierarchicalSystemIc<3>;
