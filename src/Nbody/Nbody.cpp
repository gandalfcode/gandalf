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
#include "SmoothingKernel.h"
#include "Nbody.h"
#if defined _OPENMP
#include <omp.h>
#endif
using namespace std;



//=================================================================================================
//  Nbody::Nbody
/// Nbody class constructor
//=================================================================================================
template <int ndim>
Nbody<ndim>::Nbody(int _nbody_softening, int _perturbers, int _sub_systems,
                   DOUBLE _nbody_mult, string KernelName, int _Npec):
#if defined _OPENMP
  maxNbodyOpenMp(omp_get_max_threads()*maxNbodyPerThread),
#endif
  nbody_softening(_nbody_softening),
  perturbers(_perturbers),
  sub_systems(_sub_systems),
  Npec(_Npec),
  nbody_mult(_nbody_mult),
  kerntab(TabulatedKernel<ndim>(KernelName))
{
  allocated     = false;
  Nnbody        = 0;
  Nnbodymax     = 0;
  Nstar         = 0;
  Nstarmax      = 0;
  Nsystem       = 0;
  Nsystemmax    = 0;
  Nstellartable = 0;
  reset_tree    = 0;
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
    Nstarmax   = N;
    Nsystemmax = N;
    Nnbodymax  = Nstarmax + Nsystemmax;
    nbodydata  = new NbodyParticle<ndim>*[Nnbodymax];
    stardata   = new StarParticle<ndim>[Nstarmax];
    system     = new SystemParticle<ndim>[Nsystemmax];
    allocated  = true;
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
 (SimUnits *simunits)                  ///< [in] Pointer to main units object
{
  string filename = "stellar.dat";     // Stellar table filename
  string dummystring;                  // Dummy string to skip unimportant lines
  ifstream infile;                     // Input stream object
  int i;                               // Table element counter

  debug2("[Nbody::LoadStellarPropertiesTable]");

  infile.open(filename.c_str());

  infile >> Nstellartable;
  stellartable = new StellarTableElement[Nstellartable];
  getline(infile, dummystring);
  getline(infile, dummystring);
  getline(infile, dummystring);
  getline(infile, dummystring);
  getline(infile, dummystring);

  for (i=0; i<Nstellartable; i++) {
    infile >> stellartable[i].mass >> stellartable[i].luminosity >> stellartable[i].NLyC
           >> stellartable[i].Teff >> stellartable[i].mdot >> stellartable[i].vwind;
    stellartable[i].mass       /= simunits->m.inscale;
    stellartable[i].luminosity /= simunits->L.inscale;
    //stellartable[i].NLyC       /= simunits->t.inscale;
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
  int i;                               // Star counter
  int jbin;                            // Stellar table bin number
  FLOAT frac;                          // Linear interpolation weighting factor
  FLOAT ms;                            // Mass of star

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
      pow(10.0, (1.0 - frac)*stellartable[jbin].NLyC + frac*stellartable[jbin+1].NLyC);

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
 (int N,                               ///< [in] Number of stars
  NbodyParticle<ndim> **star,          ///< [inout] Array of stars/systems
  DomainBox<ndim> &simbox,             ///< [in] Simulation domain box
  Ewald<ndim> *ewald)                  ///< [in] Ewald gravity object pointer
{
  FLOAT aperiodic[ndim];               // Ewald periodic grav. accel correction
  FLOAT dr[ndim];                      // Relative position vector
  FLOAT dr_corr[ndim];                 // Periodic corrected position vector
  FLOAT drdt;                          // Rate of change of distance
  FLOAT drsqd;                         // Distance squared
  FLOAT dv[ndim];                      // Relative velocity vector
  FLOAT invdrmag;                      // 1 / drmag
  FLOAT potperiodic;                   // Periodic correction for grav. potential

  debug2("[Nbody::CalculateDirectGravForces]");

  // Loop over all (active) stars
  //-----------------------------------------------------------------------------------------------
#pragma omp parallel for if (N > maxNbodyOpenMp) default(none) shared(ewald, N, simbox, star) \
private(aperiodic, dr, dr_corr, drdt, drsqd, dv, invdrmag, potperiodic)
  for (int i=0; i<N; i++) {
    if (not star[i]->flags.check(active)) continue;

    // Sum grav. contributions for all other stars (excluding star itself)
    //---------------------------------------------------------------------------------------------
    for (int j=0; j<N; j++) {
      if (i == j) continue;

      for (int k=0; k<ndim; k++) dr[k] = star[j]->r[k] - star[i]->r[k];
      for (int k=0; k<ndim; k++) dv[k] = star[j]->v[k] - star[i]->v[k];
      NearestPeriodicVector(simbox, dr, dr_corr);
      drsqd    = DotProduct(dr, dr, ndim);
      invdrmag = (FLOAT) 1.0/sqrt(drsqd);
      drdt     = DotProduct(dv,dr,ndim)*invdrmag;
      star[i]->gpot += star[j]->m*invdrmag;
      for (int k=0; k<ndim; k++) star[i]->a[k] += star[j]->m*dr[k]*pow(invdrmag,3);
      for (int k=0; k<ndim; k++) star[i]->adot[k] +=
        star[j]->m*pow(invdrmag,3)*(dv[k] - 3.0*drdt*invdrmag*dr[k]);

      // Add periodic gravity contribution (if activated)
      if (simbox.PeriodicGravity) {
        ewald->CalculatePeriodicCorrection(star[j]->m, dr, aperiodic, potperiodic);
        for (int k=0; k<ndim; k++) star[i]->a[k] += aperiodic[k];
        star[i]->gpot += potperiodic;
      }

    }
    //---------------------------------------------------------------------------------------------

  }
  //-----------------------------------------------------------------------------------------------

  return;
}



//=================================================================================================
//  SphIntegration::CheckBoundaries
/// Check all particles to see if any have crossed the simulation bounding box.
/// If so, then move the particles to their new location on the other side of the periodic box.
//=================================================================================================
template <int ndim>
void Nbody<ndim>::CheckBoundaries
 (int n,                               ///< [in] Integer time
  int N,                               ///< [in] No. of stars/systems
  FLOAT t,                             ///< [in] Current time
  FLOAT timestep,                      ///< [in] Smallest timestep value
  DomainBox<ndim> &simbox,             ///< [in] Main simulation domain box
  NbodyParticle<ndim> **star)          ///< [inout] Main star/system array
{
  debug2("[Nbody::CheckBoundaries]");


  // Loop over all particles and check if any lie outside the periodic box.
  // If so, then re-position with periodic wrapping.
  //===============================================================================================
#pragma omp parallel for if (N > maxNbodyOpenMp) default(none) shared(N,simbox,star)
  for (int i=0; i<N; i++) {

    // --------------------------------------------------------------------------------------------
    for (int k=0; k<ndim; k++) {

      // Check if particle has crossed LHS boundary
      //-------------------------------------------------------------------------------------------
      if (star[i]->r[k] < simbox.min[k]) {

        // Check if periodic boundary
        if (simbox.boundary_lhs[k] == periodicBoundary) {
          star[i]->r[k]  += simbox.size[k];
          star[i]->r0[k] += simbox.size[k];
        }

        // Check if wall or mirror boundary
        if (simbox.boundary_lhs[k] == mirrorBoundary || simbox.boundary_lhs[k] == wallBoundary) {
          star[i]->r[k]  = (FLOAT) 2.0*simbox.min[k] - star[i]->r[k];
          star[i]->r0[k] = (FLOAT) 2.0*simbox.min[k] - star[i]->r0[k];
          star[i]->v[k]  = -star[i]->v[k];
          star[i]->v0[k] = -star[i]->v0[k];
          star[i]->a[k]  = -star[i]->a[k];
          star[i]->a0[k] = -star[i]->a0[k];
        }

      }

      // Check if particle has crossed RHS boundary
      //-------------------------------------------------------------------------------------------
      if (star[i]->r[k] > simbox.max[k]) {

        // Check if periodic boundary
        if (simbox.boundary_rhs[k] == periodicBoundary) {
          star[i]->r[k]  -= simbox.size[k];
          star[i]->r0[k] -= simbox.size[k];
        }

        // Check if wall or mirror boundary
        if (simbox.boundary_rhs[k] == mirrorBoundary || simbox.boundary_rhs[k] == wallBoundary) {
          star[i]->r[k]  = (FLOAT) 2.0*simbox.max[k] - star[i]->r[k];
          star[i]->r0[k] = (FLOAT) 2.0*simbox.max[k] - star[i]->r0[k];
          star[i]->v[k]  = -star[i]->v[k];
          star[i]->v0[k] = -star[i]->v0[k];
          star[i]->a[k]  = -star[i]->a[k];
          star[i]->a0[k] = -star[i]->a0[k];
        }

      }


    }
    //---------------------------------------------------------------------------------------------

  }
  //===============================================================================================

  return;
}




//=================================================================================================
//  Nbody::CalculatePerturberForces
/// Calculate perturber tidal forces on all stars in a N-body sub-system.
//=================================================================================================
template <int ndim>
void Nbody<ndim>::CalculatePerturberForces
 (int N,                               ///< [in] Number of stars
  int Npert,                           ///< [in] Number of perturbing stars
  NbodyParticle<ndim> **star,          ///< [in] Array of stars/systems
  NbodyParticle<ndim> *perturber,      ///< [in] Array of perturbing stars/systems
  DomainBox<ndim> &simbox,             ///< [in] Simulation domain box
  Ewald<ndim> *ewald,                  ///< [in] Ewald gravity object pointer
  FLOAT *apert,                        ///< [out] Tidal acceleration due to perturbers
  FLOAT *adotpert)                     ///< [out] Tidal jerk due to perturbers
{
  int i,j,k;                           // Star and dimension counters
  FLOAT aperiodic[ndim];               // Ewald periodic grav. accel correction
  FLOAT dr[ndim];                      // Relative position vector
  FLOAT dr_corr[ndim];                 // Periodic corrected position vector
  FLOAT drdt;                          // Rate of change of distance
  FLOAT drsqd;                         // Distance squared
  FLOAT dv[ndim];                      // Relative velocity vector
  FLOAT invdrmag;                      // 1 / drmag
  FLOAT potperiodic;                   // Periodic correction for grav. potential
  FLOAT rcom[ndim];                    // Position of centre-of-mass
  FLOAT vcom[ndim];                    // Velocity of centre-of-mass
  FLOAT msystot = (FLOAT) 0.0;         // Total system mass

  debug2("[Nbody::CalculatePerturberForces]");

  // First, compute position and velocity of system COM
  for (k=0; k<ndim; k++) rcom[k] = (FLOAT) 0.0;
  for (k=0; k<ndim; k++) vcom[k] = (FLOAT) 0.0;
  for (i=0; i<N; i++) {
    msystot += star[i]->m;
    for (k=0; k<ndim; k++) rcom[k] += star[i]->m*star[i]->r[k];
    for (k=0; k<ndim; k++) vcom[k] += star[i]->m*star[i]->v[k];
  }
  for (k=0; k<ndim; k++) rcom[k] /= msystot;
  for (k=0; k<ndim; k++) vcom[k] /= msystot;

  // Calculate the accel. and jerk of the perturber on the system COM
  for (j=0; j<Npert; j++) {
    for (k=0; k<ndim; k++) dr[k] = rcom[k] - perturber[j].r[k];
    NearestPeriodicVector(simbox, dr, dr_corr);
    for (k=0; k<ndim; k++) dv[k] = vcom[k] - perturber[j].v[k];
    drsqd = DotProduct(dr, dr, ndim);
    invdrmag = (FLOAT) 1.0/sqrt(drsqd);
    drdt = DotProduct(dv, dr, ndim)*invdrmag;
    for (k=0; k<ndim; k++) apert[ndim*j + k] = -msystot*dr[k]*pow(invdrmag,3);
    for (k=0; k<ndim; k++) adotpert[ndim*j + k] =
      -msystot*pow(invdrmag,3)*(dv[k] - (FLOAT) 3.0*drdt*invdrmag*dr[k]);

    // Add periodic gravity contribution (if activated)
    if (simbox.PeriodicGravity) {
      ewald->CalculatePeriodicCorrection(msystot, dr, aperiodic, potperiodic);
      for (k=0; k<ndim; k++) apert[k] += aperiodic[k];
    }
  }

  // Loop over all (active) stars
  //-----------------------------------------------------------------------------------------------
  for (i=0; i<N; i++) {
    if (not star[i]->flags.check(active)) continue;

    // Sum grav. contributions for all perturbing stars.
    //---------------------------------------------------------------------------------------------
    for (j=0; j<Npert; j++) {

      for (k=0; k<ndim; k++) dr[k] = perturber[j].r[k] - star[i]->r[k];
      NearestPeriodicVector(simbox, dr, dr_corr);

      for (k=0; k<ndim; k++) dv[k] = perturber[j].v[k] - star[i]->v[k];
      drsqd = DotProduct(dr, dr, ndim);
      invdrmag = (FLOAT) 1.0/sqrt(drsqd);
      drdt = DotProduct(dv,dr,ndim)*invdrmag;

      // First, add contribution of perturber to star
      star[i]->gpe_pert += star[i]->m*perturber[j].m*invdrmag;
      star[i]->gpot += perturber[j].m*invdrmag;
      for (k=0; k<ndim; k++) star[i]->a[k] += perturber[j].m*dr[k]*pow(invdrmag,3);
      for (k=0; k<ndim; k++) star[i]->adot[k] += perturber[j].m*
        pow(invdrmag,3)*(dv[k] - (FLOAT) 3.0*drdt*invdrmag*dr[k]);

      // Next, add contribution of star to perturber
      for (k=0; k<ndim; k++) apert[ndim*j + k] -= star[i]->m*dr[k]*pow(invdrmag,3);
      for (k=0; k<ndim; k++) adotpert[ndim*j + k] -=
        star[i]->m*pow(invdrmag,3)*(dv[k] - (FLOAT) 3.0*drdt*invdrmag*dr[k]);

      // Add periodic gravity contribution (if activated)
      if (simbox.PeriodicGravity) {
        ewald->CalculatePeriodicCorrection(star[i]->m, dr, aperiodic, potperiodic);
        for (k=0; k<ndim; k++) apert[k] += aperiodic[k];
      }

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
 (SystemParticle<ndim>* systemi,       ///< [inout] System to integrate the internal motionv for
  const int n,                         ///< [in]    Integer time
  const FLOAT tstart,                  ///< [in]    Initial (local) simulation time
  const FLOAT tend,                    ///< [in]    Final (current) simulation
  DomainBox<ndim> &simbox,             ///< [in]    Simulation domain box
  Ewald<ndim> *ewald)                  ///< [in]    Ewald gravity object pointer
{
  int i;                               // Particle counter
  int it;                              // Iteration counter
  int k;                               // Dimension counter
  //int Nstar;                           // Total no. of stars
  int nsteps_local=0;                  // Local no. of steps
  FLOAT aext[ndim];                    // Acceleration due to external stars
  FLOAT adotext[ndim];                 // Jerk due to external stars
  FLOAT a2dotext[ndim];                // Snap due to external stars
  FLOAT dt;                            // Local timestep
  FLOAT gpotext;                       // Grav. potential due to external sources
  FLOAT tlocal=tstart;                 // Local time counter
  FLOAT tpert;                         // Time since beginning of step for perturbing stars
  FLOAT *apert = 0;                    // Acceleration for perturbers
  FLOAT *adotpert = 0;                 // Jerk for perturbers
  NbodyParticle<ndim>** children;      // Child systems
  NbodyParticle<ndim>* perturber = 0;  // Local array of perturber properties
  const int Nchildren = systemi->Nchildren;
  const int Npert = systemi->Npert;


  // Only integrate internal motion once COM motion has finished
  if (n - systemi->nlast != systemi->nstep) return;

  debug2("[Nbody::IntegrateInternalMotion]");

  // Allocate memory for both stars and perturbers
  //Nstar     = Nchildren + Npert;
  children  = systemi->children;

  // Record total acceleration and jerk terms in order to compute external force
  // (i.e. force due to all object outside system that are not perturbers)
  for (k=0; k<ndim; k++) aext[k] = systemi->m*systemi->a0[k];
  for (k=0; k<ndim; k++) adotext[k] = systemi->m*systemi->adot0[k];
  for (k=0; k<ndim; k++) a2dotext[k] = systemi->m*systemi->a2dot0[k];
  gpotext = systemi->m*systemi->gpot;


  // If using perturbers, record local copies and remove contribution to
  // external acceleration and jerk terms
  //-----------------------------------------------------------------------------------------------
  if (perturbers == 1 && Npert > 0) {
    perturber = new NbodyParticle<ndim>[Npert];
    apert     = new FLOAT[ndim*Npert];
    adotpert  = new FLOAT[ndim*Npert];

    // Create local copies of perturbers
    for (i=0; i<Npert; i++) {
      perturber[i].m     = systemi->perturber[i]->m;
      perturber[i].nlast = systemi->perturber[i]->nlast;
      perturber[i].tlast = systemi->perturber[i]->tlast;

      for (k=0; k<ndim; k++) perturber[i].apert[k]    = (FLOAT) 0.0;
      for (k=0; k<ndim; k++) perturber[i].adotpert[k] = (FLOAT) 0.0;
      for (k=0; k<ndim; k++) perturber[i].r0[k]       = systemi->perturber[i]->r0[k];
      for (k=0; k<ndim; k++) perturber[i].v0[k]       = systemi->perturber[i]->v0[k];
      for (k=0; k<ndim; k++) perturber[i].a0[k]       = systemi->perturber[i]->a0[k];
      for (k=0; k<ndim; k++) perturber[i].adot[k]     = systemi->perturber[i]->adot0[k];

      // Set properties of perturbers at beginning of current step
      tpert = tlocal - perturber[i].tlast;
      for (k=0; k<ndim; k++) perturber[i].r[k] = perturber[i].r0[k] + perturber[i].v0[k]*tpert +
        (FLOAT) 0.5*perturber[i].a0[k]*tpert*tpert +
        onesixth*perturber[i].adot0[k]*tpert*tpert*tpert;
      for (k=0; k<ndim; k++) perturber[i].v[k] = perturber[i].v0[k] +
        perturber[i].a0[k]*tpert + (FLOAT) 0.5*perturber[i].adot0[k]*tpert*tpert;
      for (k=0; k<ndim; k++) perturber[i].a[k] = perturber[i].a0[k] + perturber[i].adot0[k]*tpert;
      for (k=0; k<ndim; k++) perturber[i].adot[k] = perturber[i].adot0[k];
    }
  }


  // Initialise child properties
  //-----------------------------------------------------------------------------------------------
  for (i=0; i<Nchildren; i++) {
    for (k=0; k<ndim; k++) children[i]->a[k]        = (FLOAT) 0.0;
    for (k=0; k<ndim; k++) children[i]->adot[k]     = (FLOAT) 0.0;
    for (k=0; k<ndim; k++) children[i]->apert[k]    = (FLOAT) 0.0;
    for (k=0; k<ndim; k++) children[i]->adotpert[k] = (FLOAT) 0.0;
    children[i]->gpot         = (FLOAT) 0.0;
    children[i]->gpe          = (FLOAT) 0.0;
    children[i]->gpe_pert     = (FLOAT) 0.0;
    children[i]->gpe_internal = (FLOAT) 0.0;
    children[i]->flags.set(active);
  }

  if (perturbers == 1 && Npert > 0) {
    this->CalculatePerturberForces(Nchildren, Npert, children, perturber, simbox, ewald, apert, adotpert);
    for (i=0; i<Nchildren; i++) {
      for (k=0; k<ndim; k++) aext[k]     -= children[i]->m*children[i]->a[k];
      for (k=0; k<ndim; k++) adotext[k]  -= children[i]->m*children[i]->adot[k];
      for (k=0; k<ndim; k++) a2dotext[k] -= children[i]->m*children[i]->a2dot[k];
      gpotext -= children[i]->gpe_pert;
    }
  }
  for (k=0; k<ndim; k++) aext[k]     /= systemi->m;
  for (k=0; k<ndim; k++) adotext[k]  /= systemi->m;
  for (k=0; k<ndim; k++) a2dotext[k] /= systemi->m;
  gpotext /= systemi->m;

  // Calculate forces, derivatives and other terms
  CalculateDirectGravForces(Nchildren, children, simbox, ewald);

  // Add perturbing force and jerk to all child stars
  for (i=0; i<Nchildren; i++) {
    for (k=0; k<ndim; k++) children[i]->a[k]    += aext[k];
    for (k=0; k<ndim; k++) children[i]->adot[k] += adotext[k];
    for (k=0; k<ndim; k++) children[i]->a0[k]    = children[i]->a[k];
    for (k=0; k<ndim; k++) children[i]->adot0[k] = children[i]->adot[k];
  }

  // Calculate higher order derivatives
  this->CalculateAllStartupQuantities(Nchildren, children, simbox, ewald);


  // Main time integration loop
  //===============================================================================================
  do {

    // Calculate time-step
    dt = std::min(big_number, (FLOAT) 1.00000000000001*(tend - tlocal));
    for (i=0; i<Nchildren; i++) {
      children[i]->nlast = n - 1;
      children[i]->nstep = 1;
      dt = std::min(dt, (FLOAT) Timestep(children[i]));
    }
    nsteps_local += 1;
    tlocal += dt;

    // Advance position and velocities
    AdvanceParticles(n, Nchildren, tlocal, dt, children);

    // Advance positions and velocities of perturbers
    if (perturbers == 1 && Npert > 0) {
      for (i=0; i<Npert; i++) {
        tpert = tlocal - perturber[i].tlast;
        for (k=0; k<ndim; k++) perturber[i].r[k] = perturber[i].r0[k] + perturber[i].v0[k]*tpert +
          (FLOAT) 0.5*perturber[i].a0[k]*tpert*tpert +
          onesixth*perturber[i].adot0[k]*tpert*tpert*tpert;
        for (k=0; k<ndim; k++) perturber[i].v[k] = perturber[i].v0[k] +
          perturber[i].a0[k]*tpert + (FLOAT) 0.5*perturber[i].adot0[k]*tpert*tpert;
        for (k=0; k<ndim; k++) perturber[i].a[k] = perturber[i].a0[k] +
          perturber[i].adot0[k]*tpert;
      }
    }


    // Time-symmetric iteration loop
    //---------------------------------------------------------------------------------------------
    for (it=0; it<Npec; it++) {

      // Zero all acceleration terms
      for (i=0; i<Nchildren; i++) {
        children[i]->flags.set(active);
        children[i]->gpot     = gpotext;
        children[i]->gpe_pert = children[i]->m*gpotext;
        for (k=0; k<ndim; k++) children[i]->a[k]        = aext[k];
        for (k=0; k<ndim; k++) children[i]->adot[k]     = adotext[k];
        for (k=0; k<ndim; k++) children[i]->a2dot[k]    = a2dotext[k];
        for (k=0; k<ndim; k++) children[i]->apert[k]    = (FLOAT) 0.0;
        for (k=0; k<ndim; k++) children[i]->adotpert[k] = (FLOAT) 0.0;
      }

      // Calculate forces, derivatives and other terms
      CalculateDirectGravForces(Nchildren, children, simbox, ewald);

      // Add perturbation terms
      if (perturbers == 1 && Npert > 0) {
        this->CalculatePerturberForces(Nchildren, Npert, children, perturber, simbox, ewald, apert, adotpert);
      }

      // Apply correction terms
      CorrectionTerms(n, Nchildren, tlocal, dt, children);

    }
    //---------------------------------------------------------------------------------------------

    // Add perturber forces to local arrays
    //if (perturbers == 1 && Npert > 0) {
    //  for (i=0; i<Npert; i++) {
    //    for (k=0; k<ndim; k++) perturber[i].apert[k] += apert[ndim*i + k]*dt;
    //    for (k=0; k<ndim; k++) perturber[i].adotpert[k] += adotpert[ndim*i + k]*dt;
    //  }
    //}

    // Now loop over children and, if they are systems, integrate their internal motion
    for (i=0; i<Nchildren; i++) {
      if (children[i]->Ncomp > 1) {
        //cout << "Integrating internal motion again!! : " << i << "    " << children[i]->Ncomp << endl;
        // The cast is needed because the function is defined only in SystemParticle, not in
        // NbodyParticle.  The safety of the cast relies on the correctness of the Ncomp value.
        IntegrateInternalMotion(static_cast<SystemParticle<ndim>* > (children[i]),
                                n, tlocal - dt, tlocal, simbox, ewald);
      }
    }

    // Calculate correction terms on perturbing stars due to sub-systems
    //if (perturbers == 1 && Npert > 0) {
    //  this->PerturberCorrectionTerms(n, Nchildren, tlocal, dt, children);
    //  CorrectionTerms(n, Nchildren, tlocal, dt, children);
    //}

    // Correct positions of all child stars in any hierarchical systems
    this->UpdateChildStars(systemi);
    //for (i=0; i<Nchildren; i++) {
      //cout << "NCOMP? : " << i << "   " << children[i]->Ncomp << "    " << Nchildren << endl;
    //  if (children[i]->Ncomp > 1) {
    //    this->UpdateChildStars(static_cast<SystemParticle<ndim>* > (children[i]));
    //  }
    //}

    // Set end-of-step variables
    EndTimestep(n, Nchildren, tlocal, dt, children);


  } while (tlocal < tend);
  //===============================================================================================


  // Copy children back to main coordinate system
  //-----------------------------------------------------------------------------------------------
  for (i=0; i<Nchildren; i++) {
    children[i]->gpe = children[i]->m*children[i]->gpot;
  }

  systemi->gpe_internal = (FLOAT) 0.0;
  for (i=0; i<Nchildren; i++) {
    systemi->gpe_internal += (FLOAT) 0.5*children[i]->m*
      (children[i]->gpot - children[i]->gpe_pert/children[i]->m - gpotext);
  }


  // Finally, add perturbations on perturber itself to main arrays before deallocating local memory
  if (perturbers == 1 && Npert > 0) {
    //for (i=0; i<Npert; i++) {
    //  for (k=0; k<ndim; k++) systemi->perturber[i]->apert[k] += perturber[i].apert[k];
    //  for (k=0; k<ndim; k++) systemi->perturber[i]->adotpert[k] += perturber[i].adotpert[k];
    //}
    delete[] adotpert;
    delete[] apert;
    delete[] perturber;
  }

  return;
}



//=================================================================================================
//  Nbody::UpdateChildStars
/// Update/correct the positions and velocities of all child stars in sub-systems to correctly
/// represent the new COM of the system.
//=================================================================================================
template <int ndim>
void Nbody<ndim>::UpdateChildStars
 (SystemParticle<ndim>* systemi)       ///< [inout] ..
{
  int i;                               // ..
  int k;                               // ..
  int Nchildren;                       // ..
  FLOAT msystot = (FLOAT)0.0;          // ..
  FLOAT rcom[ndim];                    // ..
  FLOAT vcom[ndim];                    // ..
  FLOAT acom[ndim];                    // ..
  FLOAT adotcom[ndim];                 // ..
  FLOAT a2dotcom[ndim];                // ..
  NbodyParticle<ndim>** children;      // Child systems

  // Only correct positions at end of step
  //if (n - systemi->nlast != systemi->nstep) return;

  debug2("[Nbody::UpdateChildStars]");

  // Allocate memory for both stars and perturbers
  Nchildren = systemi->Nchildren;
  children = systemi->children;

  // First calculate old COM
  for (k=0; k<ndim; k++) rcom[k]     = (FLOAT) 0.0;
  for (k=0; k<ndim; k++) vcom[k]     = (FLOAT) 0.0;
  for (k=0; k<ndim; k++) acom[k]     = (FLOAT) 0.0;
  for (k=0; k<ndim; k++) adotcom[k]  = (FLOAT) 0.0;
  for (k=0; k<ndim; k++) a2dotcom[k] = (FLOAT) 0.0;

  for (i=0; i<Nchildren; i++) {
    msystot += children[i]->m;
    for (k=0; k<ndim; k++) rcom[k]     += children[i]->m*children[i]->r[k];
    for (k=0; k<ndim; k++) vcom[k]     += children[i]->m*children[i]->v[k];
    for (k=0; k<ndim; k++) acom[k]     += children[i]->m*children[i]->a[k];
    for (k=0; k<ndim; k++) adotcom[k]  += children[i]->m*children[i]->adot[k];
    for (k=0; k<ndim; k++) a2dotcom[k] += children[i]->m*children[i]->a2dot[k];
  }

  for (k=0; k<ndim; k++) rcom[k]     /= msystot;
  for (k=0; k<ndim; k++) vcom[k]     /= msystot;
  for (k=0; k<ndim; k++) acom[k]     /= msystot;
  for (k=0; k<ndim; k++) adotcom[k]  /= msystot;
  for (k=0; k<ndim; k++) a2dotcom[k] /= msystot;


  // Now translate positions to new COM
  for (i=0; i<Nchildren; i++) {
    for (k=0; k<ndim; k++) children[i]->r[k] += systemi->r[k] - rcom[k];
    for (k=0; k<ndim; k++) children[i]->v[k] += systemi->v[k] - vcom[k];
    for (k=0; k<ndim; k++) children[i]->a[k] += systemi->a[k] - acom[k];
    for (k=0; k<ndim; k++) children[i]->adot[k] += systemi->adot[k] - adotcom[k];
    for (k=0; k<ndim; k++) children[i]->a2dot[k] += systemi->a2dot[k] - a2dotcom[k];
    for (k=0; k<ndim; k++) children[i]->r0[k] = children[i]->r[k];
    for (k=0; k<ndim; k++) children[i]->v0[k] = children[i]->v[k];
    for (k=0; k<ndim; k++) children[i]->a0[k] = children[i]->a[k];
    for (k=0; k<ndim; k++) children[i]->adot0[k] = children[i]->adot[k];
    for (k=0; k<ndim; k++) children[i]->a2dot0[k] = children[i]->a2dot[k];
  }

  // Now update the positions of any 'grand-children' (i.e. for hierarchies)
  for (i=0; i<Nchildren; i++) {
    if (children[i]->Ncomp > 1) {
      UpdateChildStars(static_cast<SystemParticle<ndim>* > (children[i]));
    }
  }


  return;
}



template class Nbody<1>;
template class Nbody<2>;
template class Nbody<3>;
