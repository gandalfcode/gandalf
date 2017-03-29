//=================================================================================================
//  Supernova.cpp
//  All routines for creating Supernova events via adding new particles to simulation
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


#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <math.h>
#include "Debug.h"
#include "Exception.h"
#include "Hydrodynamics.h"
#include "Ic.h"
#include "InlineFuncs.h"
#include "Nbody.h"
#include "NbodyParticle.h"
#include "Parameters.h"
#include "Precision.h"
#include "Sinks.h"
#include "StarParticle.h"
#include "Supernova.h"
using namespace std;





//=================================================================================================
//  Supernova::SupernovaInjection
/// Supernova injection
/// Pass: SNpos, Einj, R_therm_kin, Minj, Rinj, SNid, hydro
//=================================================================================================
template <int ndim>
void Supernova<ndim>::SupernovaInjection
 (const int n,                               ///< Current integer time
  const int level_step,                      ///< ..
  const int level_max,                       ///< Max. block timestep level
  const int SNid,                            ///< i.d. of supernova
  FLOAT t,                                   ///< Physical simulation time
  FLOAT SNpos[ndim],                         ///< Position of supernova
  FLOAT Einj,                                ///< Total energy injection of supernova
  FLOAT R_therm_kin,                         ///< ..
  FLOAT Minj,                                ///< Total injected mass
  FLOAT Rinj,                                ///< ..
  Hydrodynamics<ndim> *hydro,                ///< Pointer to hydrodynamics object
  NeighbourSearch<ndim> *neibsearch,         ///< Pointer to neighbour search object
  RandomNumber *randnumb)                    ///< Pointer to random number generator
{
  int i;                                     // ..
  int j;                                     // ..
  int k;                                     // ..
  int Ninject = (int) (Minj/hydro->mmean);   // No. of new hydro ptcls to inject as hot SN gas
  int Nneib;                                 // ..
  int Nneibmax = 100;                        // ..
  int nSNinject;                             // No. of ptcls accelerated by the SN (ninject + counter)
  FLOAT *pos;                                // position array for all injected particles
  FLOAT dr[ndim];                            // ..
  FLOAT vpart[ndim];                         // radial velocities for injected particles
  FLOAT rpart[ndim];                         // ..
  FLOAT drsqd;                               // distance^2 from SNpos
  FLOAT drmag;                               // ..
  FLOAT vrad_mag;                            // radial velocity magnitude [code units velocity]
  FLOAT etherm_mag;                          // thermal energy to be injected per particle [code units E]
  FLOAT uinj;                                // internal energy to be added to the particles [code units u = erg/g]
  int *neiblist = new int[Nneibmax];         // ..

#ifdef MPI_PARALLEL
  string message = "Supernova injection currently not working with MPI!";
  ExceptionHandler::getIstance().raise(message);
#endif

  // Randomly draw the new position within Rinj
  pos = new FLOAT[ndim*Ninject];
  Ic<ndim>::AddRandomSphere(Ninject, SNpos, Rinj, pos, randnumb);
  //Ninject = Ic<ndim>::AddLatticeSphere(Ninject, SNpos, Rinj, "hexagonal_lattice", pos, randnumb);

  // Do a neighbour search to find existing particles inside sphere.  If memory is not large
  // enough, then increase the size.
  Nneib = neibsearch->GetGatherNeighbourList(SNpos, Rinj, hydro->GetParticleArrayUnsafe(),
                                             hydro->Nhydro, Nneibmax, neiblist);
  while (Nneib == -1) {
    delete[] neiblist;
    Nneibmax *= 2;
    neiblist = new int[Nneibmax];
    Nneib = neibsearch->GetGatherNeighbourList(SNpos, Rinj, hydro->GetParticleArrayUnsafe(),
                                               hydro->Nhydro, Nneibmax, neiblist);
  };

  // Total number of hot SN particles = new particles + existing particles within sphere
  nSNinject = Ninject + Nneib;
#ifdef OUTPUT_ALL
  cout << "Adding " << Ninject << " new particles.  Heating " << Nneib << " other particles" << endl;
  cout << "Rinj : " << Rinj << endl;
#endif

  // Give the particles their radial velocities
  vrad_mag = sqrt(2./( (FLOAT) nSNinject)/(hydro->mmean)*Einj*1./(R_therm_kin+1.));
  etherm_mag = (1.0/(1.0 + 1.0/R_therm_kin))*Einj/( (FLOAT) nSNinject);
  uinj = etherm_mag/(hydro->mmean);
#ifdef OUTPUT_ALL
  cout << "etot : " << (1.0/(1.0 + 1.0/R_therm_kin))*Einj
       << "    etherm_mag : " << etherm_mag << "     uinj : " << uinj << endl;
#endif

  // Loop over existing particles and accelerate them
  //-----------------------------------------------------------------------------------------------
  for (j=0; j<Nneib; j++) {
    i = neiblist[j];

    Particle<ndim>& part = hydro->GetParticlePointer(i);

    // determine radial velocity vector for this particle
    for (k=0; k<ndim; k++) dr[k] = part.r[k] - SNpos[k];
    drsqd = DotProduct(dr, dr, ndim);
    drmag = sqrt(drsqd) + small_number;
    part.u += uinj;

    for (k=0;k<ndim;k++){
      part.v[k] = (part.r[k] - SNpos[k])/drmag * vrad_mag;

      // reset acceleration to 0??
      part.a[k] = (FLOAT) 0.0;
    }
  }

  // Loop over all new particles, create them one by one and initialize
  //-----------------------------------------------------------------------------------------------
  for (j=0; j<Ninject; j++){

    for (k=0; k<ndim; k++) rpart[k] = pos[ndim*j+k];

    for (k=0; k<ndim; k++) dr[k] = rpart[k] - SNpos[k];
    drsqd = DotProduct(dr, dr, ndim);
    drmag = sqrt(drsqd) + small_number;
    for (k=0; k<ndim; k++) vpart[k] = (rpart[k] - SNpos[k])/drmag * vrad_mag;

    // Get new particle i
    hydro->CreateNewParticle(gas_type, hydro->mmean, uinj, rpart, vpart, sim);

  }
  //-----------------------------------------------------------------------------------------------


  delete[] pos;
  delete[] neiblist;

  return;
}



template class Supernova<1>;
template class Supernova<2>;
template class Supernova<3>;
