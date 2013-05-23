//=============================================================================
//  Sph.cpp
//  Contains important default routines for Sph class.
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
using namespace std;

template <int ndim>
const FLOAT Sph<ndim>::invndim;


//=============================================================================
//  Sph::Sph
/// ..
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
/// Allocate main SPH particle array.
//=============================================================================
template <int ndim>
void Sph<ndim>::AllocateMemory(int N)
{
  debug2("[Sph::AllocateMemory]");

  if (N > Nsphmax) {
    if (allocated) DeallocateMemory();
    Nsph = N;
    //TODO: perhaps this 10 could be made a user-provided parameter
    //(to handle the case where one doesn't want to waste memory)
    Nsphmax = 10*N;
    sphdata = new struct SphParticle<ndim>[Nsphmax];
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

  if (allocated) delete[] sphdata;
  allocated = false;

  return;
}



//=============================================================================
//  Sph::SphBoundingBox
/// Calculate the bounding box containing all SPH particles.
//=============================================================================
template <int ndim>
void Sph<ndim>::SphBoundingBox
(FLOAT rmax[ndim],                  ///< Maximum extent of bounding box
 FLOAT rmin[ndim],                  ///< Minimum extent of bounding box
 int Nmax)                          ///< Maximum particle i.d. in loop
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
  int Ngather;                      // No. of neighbours (move!!)
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
    Ngather = 5;
    volume = rmax[0] - rmin[0];
    h_guess = (volume*(FLOAT) Ngather)/(4.0*(FLOAT) Nsph);
  }
  // --------------------------------------------------------------------------
  else if (ndim == 2) {
    Ngather = 20;
    volume = (rmax[0] - rmin[0])*(rmax[1] - rmin[1]);
    h_guess = sqrtf((volume*(FLOAT) Ngather)/(4.0*(FLOAT) Nsph));
  }
  // --------------------------------------------------------------------------
  else if (ndim == 3) {
    Ngather = 50;
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

  cout << "hguess : " << h_guess << "    volume : " << volume << endl;

  return;
}



template class Sph<1>;
template class Sph<2>;
template class Sph<3>;
