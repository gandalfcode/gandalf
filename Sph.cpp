// ============================================================================
// Sph.cpp
// ============================================================================


#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <math.h>
#include "Sph.h"
#include "SphKernel.h"
#include "SphParticle.h"
#include "Parameters.h"
#include "EOS.h"
#include "Debug.h"
using namespace std;



// ============================================================================
// Sph::AllocateMemory
// ============================================================================
void Sph::AllocateMemory(int N)
{
  debug2("[Sph::AllocateMemory]");

  if (N > Nsphmax) {
    if (allocated) DeallocateMemory();
    Nsph = N;
    Nsphmax = 10*N;
    sphdata = new struct SphParticle[Nsphmax];
    allocated = true;
  }

  return;
}



// ============================================================================
// Sph::DeallocateMemory
// ============================================================================
void Sph::DeallocateMemory(void)
{
  debug2("[Sph::DeallocateMemory]");

  delete[] sphdata;
  allocated = false;

  return;
}



// ============================================================================
// Sph::SphBoundingBox
// Calculate the bounding box containing all SPH particles.
// ============================================================================
void Sph::SphBoundingBox(float rmax[ndimmax],float rmin[ndimmax],int Nmax)
{
  debug2("[Sph::SphBoundingBox]");

  for (int k=0; k<ndimmax; k++) rmin[k] = big_number; 
  for (int k=0; k<ndimmax; k++) rmax[k] = -big_number;

  for (int i=0; i<Nmax; i++) {
    for (int k=0; k<ndim; k++) 
      if (sphdata[i].r[k] < rmin[k]) rmin[k] = sphdata[i].r[k];
    for (int k=0; k<ndim; k++) 
      if (sphdata[i].r[k] > rmax[k]) rmax[k] = sphdata[i].r[k];
  }

#if defined(DEBUG_ALL)
  printf("rmin : %f\n",rmin[0]);
  printf("rmax : %f\n",rmax[0]);
#endif

  return;
}



// ============================================================================
// Sph::InitialSmoothingLengthGuess
// Perform initial guess of smoothing.  In the abscence of more sophisticated 
// techniques, we guess the smoothing length assuming a uniform density 
// medium with the same volume and total mass.
// ============================================================================
void Sph::InitialSmoothingLengthGuess(void)
{
  float h_guess;
  float volume;
  float rmin[ndimmax];
  float rmax[ndimmax];
  int Ngather;

  debug2("[Sph::InitialSmoothingLengthGuess]");

  // Calculate bounding box containing all SPH particles
  SphBoundingBox(rmax,rmin,Nsph);

  // Depending on the dimensionality, calculate the average smoothing 
  // length assuming a uniform density distribution filling the boudning box.
  // --------------------------------------------------------------------------
  if (ndim == 1) {
    Ngather = 5;
    volume = rmax[0] - rmin[0];
    h_guess = (volume*(float) Ngather)/(4.0f*(float) Nsph);
  }
  // --------------------------------------------------------------------------
  else if (ndim == 2) {
    Ngather = 16;
    volume = (rmax[0] - rmin[0])*(rmax[1] - rmin[1]);
    h_guess = sqrtf((volume*(float) Ngather)/(4.0f*(float) Nsph));
  }
  // --------------------------------------------------------------------------
  else if (ndim == 3) {
    Ngather = 50;
    volume = (rmax[0] - rmin[0])*(rmax[1] - rmin[1])*(rmax[2] - rmin[2]);
    h_guess = powf((3.0f*volume*(float) Ngather)/
		   (32.0f*pi*(float) Nsph),onethird);
  }
  // --------------------------------------------------------------------------

  // Set all smoothing lengths equal to average value
  for (int i=0; i<Nsph; i++) {
    sphdata[i].h = h_guess;
    sphdata[i].invh = 1.0/h_guess;
  }

  printf ("Volume : %f\n",volume);
  printf ("Smoothing length : %f\n",h_guess);

  return;
}
