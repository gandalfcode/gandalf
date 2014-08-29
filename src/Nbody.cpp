//=================================================================================================
//  Nbody.cpp
//  Contains main N-body class functions.
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


#include <algorithm>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <ostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include "Precision.h"
#include "Debug.h"
#include "InlineFuncs.h"
#include "NbodyParticle.h"
#include "StarParticle.h"
#include "SystemParticle.h"
#include "Parameters.h"
#include "SphKernel.h"
#include "Nbody.h"

using namespace std;



//=================================================================================================
//  Nbody::Nbody
/// Nbody class constructor
//=================================================================================================
template <int ndim>
Nbody<ndim>::Nbody(int nbody_softening_aux, int sub_systems_aux,
                   DOUBLE nbody_mult_aux, string KernelName, int Npec_aux):
  nbody_softening(nbody_softening_aux),
  sub_systems(sub_systems_aux),
  nbody_mult(nbody_mult_aux),
  kerntab(TabulatedKernel<ndim>(KernelName)),
  Npec(Npec_aux),
  Nstar(0),
  Nstarmax(0),
  Nsystem(0),
  Nsystemmax(0),
  Nnbody(0),
  Nnbodymax(0),
  reset_tree(0),
  allocated(false),
  Nstellartable(0)
{
}



//=================================================================================================
//  Nbody::AllocateMemory
/// Allocate all memory required for stars and N-body system particles.
//=================================================================================================
template <int ndim>
void Nbody<ndim>::AllocateMemory(int N)
{
  debug2("[Nbody::AllocateMemory]");

  if (N > Nstarmax) {
    if (allocated) DeallocateMemory();
    Nstarmax = N;
    //Nsystem = N;
    Nsystemmax = N;
    //Nnbody = N;
    Nnbodymax = Nstarmax + Nsystemmax;
    nbodydata = new struct NbodyParticle<ndim>*[Nnbodymax];
    stardata = new struct StarParticle<ndim>[Nstarmax];
    system = new struct SystemParticle<ndim>[Nsystemmax];
    allocated = true;
  }

  return;
}



//=================================================================================================
//  Nbody::DeallocateMemory
/// Deallocate all N-body memory.
//=================================================================================================
template <int ndim>
void Nbody<ndim>::DeallocateMemory(void)
{
  debug2("[Nbody::DeallocateMemory]");

  if (allocated) {
    delete[] system;
    delete[] stardata;
    delete[] nbodydata;
  }
  allocated = false;

  return;
}



//=================================================================================================
//  Nbody::LoadStellarPropertiesTable
/// Loads table from file containing properties of high-mass stars, including mass-loss rates,
/// wind speeds, total luminosities, effective surface temperatures and ionising photon flux.
//=================================================================================================
template<int ndim>
void Nbody<ndim>::LoadStellarPropertiesTable
(SimUnits *simunits)
{
  string filename = "stellar.dat";  // Stellar table filename
  string dummystring;               // Dummy string to skip unimportant lines
  ifstream infile;                  // Input stream object
  int i;                            // Table element counter

  debug2("[Nbody::LoadStellarPropertiesTable]");

  infile.open(filename.c_str());

  infile >> Nstellartable;
  stellartable = new StellarTableElement[Nstellartable];
  getline(infile,dummystring);
  getline(infile,dummystring);
  getline(infile,dummystring);
  getline(infile,dummystring);
  getline(infile,dummystring);

  for (i=0; i<Nstellartable; i++) {
    infile >> stellartable[i].mass >> stellartable[i].luminosity >> stellartable[i].NLyC
           >> stellartable[i].Teff >> stellartable[i].mdot >> stellartable[i].vwind;
    stellartable[i].mass       /= simunits->m.inscale;
    stellartable[i].luminosity /= simunits->L.inscale;
    stellartable[i].Teff       /= simunits->temp.inscale;
    stellartable[i].mdot       /= simunits->dmdt.inscale;
    stellartable[i].vwind      /= simunits->v.inscale;
  }

  infile.close();

  return;
}



//=================================================================================================
//  Nbody::UpdateStellarProperties
/// Updates the stellar properties of all stars/sinks based on their new mass after accretion.
//=================================================================================================
template<int ndim>
void Nbody<ndim>::UpdateStellarProperties(void)
{
  int i;                            // Star counter
  int jbin;                         // Stellar table bin number
  DOUBLE frac;                      // Linear interpolation weighting factor
  DOUBLE ms;                        // Mass of star

  debug2("[Nbody::UpdateStellarProperties]");

  //-----------------------------------------------------------------------------------------------
  for (i=0; i<Nstar; i++) {
    ms = nbodydata[i]->m;

    // Search through all table entries and find bin that contains mass
    for (jbin=0; jbin<Nstellartable - 1; jbin++) {
      if (ms >= stellartable[jbin].mass && ms < stellartable[jbin+1].mass) break;
    }

    // If mass exceeds last bin, simply use largest mass values for stellar properties
    // Otherwise, use linear interpolation between tabulated values
    if (jbin == Nstellartable-1) {
      jbin = Nstellartable - 2;
      frac = 1.0;
    }
    else {
      frac = (ms - stellartable[jbin].mass)/(stellartable[jbin+1].mass - stellartable[jbin].mass);
    }

    // Set all stellar properties (only NLyC for now)
    nbodydata[i]->NLyC =
      pow(10.0,(1.0 - frac)*stellartable[jbin].NLyC + frac*stellartable[jbin+1].NLyC);

  }
  //-----------------------------------------------------------------------------------------------

  return;
}




//=================================================================================================
//  Nbody::CalculateDirectGravForces
/// Calculate all star-star force contributions for active systems using
/// direct summation with unsoftened gravity.
//=================================================================================================
template <int ndim>
void Nbody<ndim>::CalculateDirectGravForces
(int N,                             ///< Number of stars
 NbodyParticle<ndim> **star)        ///< Array of stars/systems
{
  int i,j,k;                        // Star and dimension counters
  DOUBLE dr[ndim];                  // Relative position vector
  DOUBLE drdt;                      // Rate of change of distance
  DOUBLE drsqd;                     // Distance squared
  DOUBLE dv[ndim];                  // Relative velocity vector
  DOUBLE invdrmag;                  // 1 / drmag

  debug2("[Nbody::CalculateDirectGravForces]");

  // Loop over all (active) stars
  //-----------------------------------------------------------------------------------------------
  for (i=0; i<N; i++) {
    if (star[i]->active == 0) continue;

    // Sum grav. contributions for all other stars (excluding star itself)
    //---------------------------------------------------------------------------------------------
    for (j=0; j<N; j++) {
      if (i == j) continue;

      for (k=0; k<ndim; k++) dr[k] = star[j]->r[k] - star[i]->r[k];
      for (k=0; k<ndim; k++) dv[k] = star[j]->v[k] - star[i]->v[k];
      drsqd = DotProduct(dr,dr,ndim);
      invdrmag = 1.0/sqrt(drsqd);
      drdt = DotProduct(dv,dr,ndim)*invdrmag;

      star[i]->gpot += star[j]->m*invdrmag;
      for (k=0; k<ndim; k++) star[i]->a[k] += star[j]->m*dr[k]*pow(invdrmag,3);
      for (k=0; k<ndim; k++) star[i]->adot[k] +=
        star[j]->m*pow(invdrmag,3)*(dv[k] - 3.0*drdt*invdrmag*dr[k]);

    }
    //---------------------------------------------------------------------------------------------

  }
  //-----------------------------------------------------------------------------------------------

  return;
}


//=================================================================================================
//  Nbody::IntegrateInternalMotion
/// This function integrates the internal motion of a system. First integrates the internal motion
/// of its sub-systems by recursively calling their method, then integrates their COMs.
//=================================================================================================
template <int ndim>
void Nbody<ndim>::IntegrateInternalMotion
(SystemParticle<ndim>* systemi,     ///< [inout] System to integrate the internal motionv for
 int n,                             ///< [in]    Integer time
 DOUBLE timestep,                   ///< [in]    Minimum timestep value
 DOUBLE tlocal_end)                 ///< [in]    Time to integrate the internal motion for
{
  int i;                            // Star counter
  int it;                           // Iteration counter
  int k;                            // Dimension counter
  int Nchildren;                    // No. of child systems
  int nlocal_steps = 0;             // No. of locally integrated steps
  DOUBLE dt;                        // Timestep
  DOUBLE tlocal=0.0;                // Local integration time
  DOUBLE rcom[ndim];                // Position of centre-of-mass
  DOUBLE vcom[ndim];                // Velocity of centre-of-mass
  DOUBLE acom[ndim];                // Acceleration of centre-of-mass
  DOUBLE adotcom[ndim];             // Jerk of centre-of-mass
  NbodyParticle<ndim>** children;   // Child systems


  Nchildren = systemi->Nchildren;
  children = systemi->children;

  // Zero all COM summation variables
  for (k=0; k<ndim; k++) rcom[k] = 0.0;
  for (k=0; k<ndim; k++) vcom[k] = 0.0;
  for (k=0; k<ndim; k++) acom[k] = 0.0;
  for (k=0; k<ndim; k++) adotcom[k] = 0.0;

  // Make local copies of children and calculate COM properties
  //-----------------------------------------------------------------------------------------------
  for (i=0; i<Nchildren; i++) {
    for (k=0; k<ndim; k++) rcom[k] += children[i]->m*children[i]->r[k];
    for (k=0; k<ndim; k++) vcom[k] += children[i]->m*children[i]->v[k];
    for (k=0; k<ndim; k++) acom[k] += children[i]->m*children[i]->a[k];
    for (k=0; k<ndim; k++) adotcom[k] += children[i]->m*children[i]->adot[k];
  }

  // Normalise COM values
  for (k=0; k<ndim; k++) rcom[k] /= systemi->m;
  for (k=0; k<ndim; k++) vcom[k] /= systemi->m;
  for (k=0; k<ndim; k++) acom[k] /= systemi->m;
  for (k=0; k<ndim; k++) adotcom[k] /= systemi->m;


  // Now convert to COM frame
  //-----------------------------------------------------------------------------------------------
  for (i=0; i<Nchildren; i++) {
    for (k=0; k<ndim; k++) children[i]->r[k] -= rcom[k];
    for (k=0; k<ndim; k++) children[i]->r0[k] -= rcom[k];
    for (k=0; k<ndim; k++) children[i]->v[k] -= vcom[k];
    for (k=0; k<ndim; k++) children[i]->v0[k] -= vcom[k];
    for (k=0; k<ndim; k++) children[i]->a[k] -= acom[k];
    for (k=0; k<ndim; k++) children[i]->a0[k] -= acom[k];
    for (k=0; k<ndim; k++) children[i]->adot[k] -= adotcom[k];
    for (k=0; k<ndim; k++) children[i]->adot0[k] -= adotcom[k];
    for (k=0; k<ndim; k++) children[i]->a2dot[k] -= 0.0;
    for (k=0; k<ndim; k++) children[i]->a3dot[k] -= 0.0;
    children[i]->active = true;
    children[i]->nstep = 1;
    children[i]->level = 0;
  }


  // Main time integration loop
  //===============================================================================================
  do {

    // Calculate time-step
    dt = std::min(big_number, tlocal_end - tlocal);
    for (i=0; i<Nchildren; i++) {
      dt = std::min(dt, Timestep(children[i]));
    }
    tlocal += dt;
    nlocal_steps +=1;

    // Advance position and velocities
    AdvanceParticles(nlocal_steps, Nchildren, children, dt);

    // Time-symmetric iteration loop
    //---------------------------------------------------------------------------------------------
    for (it=0; it<Npec; it++) {

      // Zero all acceleration terms
      for (i=0; i<Nchildren; i++) {
        children[i]->gpot = 0.0;
        children[i]->gpe_internal = 0.0;
        for (k=0; k<ndim; k++) children[i]->a[k] = 0.0;
        for (k=0; k<ndim; k++) children[i]->adot[k] = 0.0;
        for (k=0; k<ndim; k++) children[i]->a2dot[k] = 0.0;
        for (k=0; k<ndim; k++) children[i]->a3dot[k] = 0.0;
      }

      // Calculate forces, derivatives and other terms
      CalculateDirectGravForces(Nchildren, children);

      // Apply correction terms
      CorrectionTerms(nlocal_steps, Nchildren, children, dt);
    }

    // Now loop over children and, if they are systems, integrate
    // their internal motion
    //---------------------------------------------------------------------------------------------
    for (i=0; i<Nchildren; i++) {

      if (children[i]->Ncomp > 1) {
        // The cast is needed because the function is defined only in SystemParticle, not in
        // NbodyParticle.  The safety of the cast relies on the correctness of the Ncomp value.
        IntegrateInternalMotion(static_cast<SystemParticle<ndim>* >
          (children[i]), n, timestep, dt);
      }
    }

    // Set end-of-step variables
    EndTimestep(nlocal_steps, Nchildren, children);

  } while (tlocal < tlocal_end);
  //===============================================================================================


  // Copy children back to main coordinate system
  //-----------------------------------------------------------------------------------------------
  for (i=0; i<Nchildren; i++) {
    for (k=0; k<ndim; k++) children[i]->r[k] += systemi->r[k];
    for (k=0; k<ndim; k++) children[i]->r0[k] += systemi->r[k];
    for (k=0; k<ndim; k++) children[i]->v[k] += systemi->v[k];
    for (k=0; k<ndim; k++) children[i]->v0[k] += systemi->v[k];
    for (k=0; k<ndim; k++) children[i]->a[k] += systemi->a[k];
    for (k=0; k<ndim; k++) children[i]->a0[k] += systemi->a[k];
    for (k=0; k<ndim; k++) children[i]->adot[k] += systemi->adot[k];
    for (k=0; k<ndim; k++) children[i]->adot0[k] += systemi->adot[k];
    children[i]->gpot = children[i]->gpot + systemi->gpot;
  }

  return;
}



template class Nbody<1>;
template class Nbody<2>;
template class Nbody<3>;
