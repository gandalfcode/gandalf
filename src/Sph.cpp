//=============================================================================
//  Sph.cpp
//  Contains important default routines for Sph class.
//
//  This file is part of GANDALF :
//  Graphical Astrophysics code for N-body Dynamics and Lagrangian Fluids
//  https://github.com/gandalfcode/gandalf
//  Contact : gandalfcode@gmail.com
//
//  Copyright (C) 2013  D. A. Hubber, G Rosotti
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
//=============================================================================


#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <math.h>
#include "Precision.h"
#include "Sph.h"
#include "SphKernel.h"
#include "SphParticle.h"
#include "Parameters.h"
#include "EOS.h"
#include "Debug.h"
#include "InlineFuncs.h"
using namespace std;

template <int ndim>
const FLOAT Sph<ndim>::invndim;



//=============================================================================
//  Sph::Sph
/// Constructor for parent SPH class.  Initialises important variables and 
/// sets important parameters using initialialisation lists.
//=============================================================================
template <int ndim>
Sph<ndim>::Sph(int hydro_forces_aux, int self_gravity_aux, 
  FLOAT alpha_visc_aux, FLOAT beta_visc_aux, FLOAT h_fac_aux, 
  FLOAT h_converge_aux, aviscenum avisc_aux, acondenum acond_aux, 
  string gas_eos_aux, string KernelName):
  hydro_forces(hydro_forces_aux),
  self_gravity(self_gravity_aux),
  alpha_visc(alpha_visc_aux),
  beta_visc(beta_visc_aux),
  h_fac(h_fac_aux),
  h_converge(h_converge_aux),
  hmin_sink(big_number),
  gas_eos(gas_eos_aux),
  kerntab(TabulatedKernel<ndim>(KernelName)),
  allocated(false),
  Nsph(0),
  Nsphmax(0),
  avisc(avisc_aux),
  acond(acond_aux)
{
}



//=============================================================================
//  Sph::AllocateMemory
/// Allocate main SPH particle array.  Currently sets the maximum memory to 
/// be 10 times the numbers of particles to allow space for ghost particles 
/// and new particle creation.
//=============================================================================
template <int ndim>
void Sph<ndim>::AllocateMemory(int N)
{
  debug2("[Sph::AllocateMemory]");

  if (N > Nsphmax) {
    if (allocated) DeallocateMemory();
    //TODO: perhaps this 10 could be made a user-provided parameter
    //(to handle the case where one doesn't want to waste memory)
    Nsphmax = 10*N;
    iorder = new int[Nsphmax];
    rsph = new FLOAT[ndim*Nsphmax];
    sphdata = new struct SphParticle<ndim>[Nsphmax];
    sphintdata = new SphIntParticle<ndim>[Nsphmax];
    for (int i=0; i<Nsphmax; i++) {
     sphintdata[i].part=&sphdata[i];
    }
    allocated = true;
  }

  return;
}



//=============================================================================
//  Sph::DeallocateMemory
/// Deallocate main array containing SPH particle data.
//=============================================================================
template <int ndim>
void Sph<ndim>::DeallocateMemory(void)
{
  debug2("[Sph::DeallocateMemory]");

  if (allocated) {
    delete[] sphintdata;
    delete[] sphdata;
    delete[] rsph;
    delete[] iorder;
  }
  allocated = false;

  return;
}



//=============================================================================
//  Sph::DeleteParticles
/// Delete selected SPH particles from the main arrays.
//=============================================================================
template <int ndim>
void Sph<ndim>::DeleteParticles
(int Ndead,                         ///< No. of 'dead' particles
 int *deadlist)                     ///< List of 'dead' particle ids
{
  int i;                            // Particle counter
  int idead = 0;                    // Aux. dead particle counter
  int ilive = 0;                    // 'Live' particle counter

  debug2("[Sph::DeleteParticles]");

  // Make sure list of dead particles is in ascending order
  InsertionSortIds(Ndead,deadlist);

  for (i=0; i<Ndead-1; i++) {
    if (deadlist[i+1] < deadlist[i]) {
      cout << "PROBLEM WITH DEADLIST : " << i << "    " 
	   << deadlist[i] << "    " << deadlist[i+1] << endl;
      exit(0);
    }
  }

  // Determine new order of particles in arrays
  for (i=0; i<Nsph; i++) {
    if (idead < Ndead) {
      if (i == deadlist[idead]) iorder[Nsph - Ndead + idead++] = i;
      else iorder[ilive++] = i;
    }
    else iorder[ilive++] = i;
  }

  // Reorder all arrays following with new order, with dead particles at end
  ReorderParticles();

  // Reduce particle counters once dead particles have been removed
  Nsph = Nsph - Ndead;
  Ntot = Nsph;

  return;
}



//=============================================================================
//  Sph::ReorderParticles
/// Delete selected SPH particles from the main arrays.
//=============================================================================
template <int ndim>
void Sph<ndim>::ReorderParticles(void)
{
  int i;                            // ..
  SphParticle<ndim> *sphdataaux;    // ..
  SphIntParticle<ndim>* sphintdataaux;

  sphdataaux = new SphParticle<ndim>[Nsph];
  sphintdataaux = new SphIntParticle<ndim>[Nsph];

  for (i=0; i<Nsph; i++) {
    sphdataaux[i] = sphdata[i];
    sphintdataaux[i] = sphintdata[i];
  }
  for (i=0; i<Nsph; i++) {
    sphdata[i] = sphdataaux[iorder[i]];
    sphintdata[i] = sphintdataaux[iorder[i]];
  }

  delete[] sphdataaux;
  delete[] sphintdataaux;

  return;
}




//=============================================================================
//  Sph::SphBoundingBox
/// Calculate the bounding box containing all SPH particles.
//=============================================================================
template <int ndim>
void Sph<ndim>::SphBoundingBox
(FLOAT rmax[ndim],                  ///< [out] Maximum extent of bounding box
 FLOAT rmin[ndim],                  ///< [out] Minimum extent of bounding box
 int Nmax)                          ///< [in] Maximum particle i.d. in loop
{
  debug2("[Sph::SphBoundingBox]");

  for (int k=0; k<ndim; k++) rmin[k] = big_number;
  for (int k=0; k<ndim; k++) rmax[k] = -big_number;

  for (int i=0; i<Nmax; i++) {
    for (int k=0; k<ndim; k++) 
      if (sphdata[i].r[k] < rmin[k]) rmin[k] = sphdata[i].r[k];
    for (int k=0; k<ndim; k++) 
      if (sphdata[i].r[k] > rmax[k]) rmax[k] = sphdata[i].r[k];
  }

  return;
}



//=============================================================================
//  Sph::InitialSmoothingLengthGuess
/// Perform initial guess of smoothing.  In the abscence of more sophisticated 
/// techniques, we guess the smoothing length assuming a uniform density 
/// medium with the same volume and total mass.
//=============================================================================
template <int ndim>
void Sph<ndim>::InitialSmoothingLengthGuess(void)
{
  FLOAT h_guess;                    // Global guess of smoothing length
  FLOAT volume;                     // Volume of global bounding box
  FLOAT rmin[ndim];                 // Min. extent of bounding box
  FLOAT rmax[ndim];                 // Max. extent of bounding box

  debug2("[Sph::InitialSmoothingLengthGuess]");

  // Calculate bounding box containing all SPH particles
  SphBoundingBox(rmax,rmin,Nsph);

  // Depending on the dimensionality, calculate the average smoothing 
  // length assuming a uniform density distribution filling the bounding box.
  // --------------------------------------------------------------------------
  if (ndim == 1) {
    Ngather = (int) (2.0*kernp->kernrange*h_fac);
    volume = rmax[0] - rmin[0];
    h_guess = (volume*(FLOAT) Ngather)/(4.0*(FLOAT) Nsph);
  }
  // --------------------------------------------------------------------------
  else if (ndim == 2) {
    Ngather = (int) (pi*pow(kernp->kernrange*h_fac,2));
    volume = (rmax[0] - rmin[0])*(rmax[1] - rmin[1]);
    h_guess = sqrtf((volume*(FLOAT) Ngather)/(4.0*(FLOAT) Nsph));
  }
  // --------------------------------------------------------------------------
  else if (ndim == 3) {
    Ngather = (int) (4.0*pi*pow(kernp->kernrange*h_fac,3)/3.0);
    volume = (rmax[0] - rmin[0])*(rmax[1] - rmin[1])*(rmax[2] - rmin[2]);
    h_guess = powf((3.0*volume*(FLOAT) Ngather)/
		   (32.0*pi*(FLOAT) Nsph),onethird);
  }
  // --------------------------------------------------------------------------

  // Set all smoothing lengths equal to average value
  for (int i=0; i<Nsph; i++) {
    sphdata[i].h = h_guess;
    sphdata[i].invh = 1.0/h_guess;
  }

  return;
}



template class Sph<1>;
template class Sph<2>;
template class Sph<3>;
