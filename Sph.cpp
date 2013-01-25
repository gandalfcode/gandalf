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
using namespace std;



// ============================================================================
// Sph::AllocateMemory
// ============================================================================
void Sph::AllocateMemory(int N)
{
  printf("[Sph::AllocateMemory]");

  if (N > Nsphmax) {
    Nsph = N;
    Nsphmax = N;
    sphdata = new struct SphParticle[Nsphmax];
  }

  return;
}



// ============================================================================
// Sph::AllocateMemory
// ============================================================================

void Sph::DeallocateMemory(void)
{
  return;
}



// ============================================================================
// Sph::SphBoundingBox
// Calculate the bounding box containing all SPH particles.
// ============================================================================
void Sph::SphBoundingBox(float rmax[ndimmax],float rmin[ndimmax])
{
  for (int k=0; k<ndim; k++) rmin[k] = big_number; 
  for (int k=0; k<ndim; k++) rmax[k] = -big_number;

  for (int i=0; i<Nsph; i++) {
    for (int k=0; k<ndim; k++) 
      if (sphdata[i].r[k] < rmin[k]) rmin[k] = sphdata[i].r[k];
    for (int k=0; k<ndim; k++) 
      if (sphdata[i].r[k] > rmax[k]) rmax[k] = sphdata[i].r[k];
  }

#if defined(DEBUG_ALL)
  printf("rmin : %f %f %f\n",rmin[0],rmin[1],rmax[2]);
  printf("rmax : %f %f %f\n",rmax[0],rmax[1],rmax[2]);
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

  // Calculate bounding box containing all SPH particles
  SphBoundingBox(rmax,rmin);

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
    sphdata[i].invh = 1.0f/h_guess;
  }

  printf ("Smoothing length : %f\n",h_guess);

  return;
}


