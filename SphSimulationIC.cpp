// ============================================================================
// SphSimulationIC.cpp
// ============================================================================


#include <iostream>
#include <string>
#include <cstdio>
#include <cstring>
#include "SphSimulation.h"
#include "Sph.h"
#include "Parameters.h"
#include "Debug.h"
using namespace std;



// ============================================================================
// SphSimulation::ShockTube
// ============================================================================
void SphSimulation::ShockTube(void)
{
  int i;
  int j;
  int k;
  int Nbox1;
  int Nbox2;
  float volume;
  float *r;
  DomainBox box1;
  DomainBox box2;

  int Nlattice1[ndimmax];
  int Nlattice2[ndimmax];
  float vfluid1[ndimmax];
  float vfluid2[ndimmax];
  float rhofluid1 = 1.0;
  float rhofluid2 = 1.0;
  float press1 = 1.0;
  float press2 = 1.0;

  Nlattice1[0] = 128;
  Nlattice2[0] = 128;
  vfluid1[0] = 4.0;
  vfluid2[0] = -4.0;

  debug2("[SphSimulation::ShockTube]");

  if (ndim == 1) {
    box1.boxmin[0] = simbox.boxmin[0];
    box1.boxmax[0] = 0.0;
    box2.boxmin[0] = 0.0;
    box2.boxmax[0] = simbox.boxmax[0];
    volume = box1.boxmax[0] - box1.boxmin[0];
    Nbox1 = Nlattice1[0];
    Nbox2 = Nlattice2[0];
  }

  sph->Nsph = Nbox1 + Nbox2;
  sph->AllocateMemory(sph->Nsph);
  cout << "Allocating memory : " << sph->Nsph << endl;;


  r = new float[ndim*sph->Nsph];

  if (Nbox1 > 0) {
    AddRegularLattice(Nbox1,Nlattice1,r,box1);

    for (i=0; i<Nbox1; i++) {
      for (k=0; k<ndim; k++) sph->sphdata[i].r[k] = r[ndim*i + k];
      for (k=0; k<ndim; k++) sph->sphdata[i].v[k] = 0.0;
      sph->sphdata[i].v[0] = vfluid1[0];
      sph->sphdata[i].m = rhofluid1*volume/(float) Nbox1;
      sph->sphdata[i].u = 1.5; //temp0/gammaone/mu_bar;
    }
  }

  if (Nbox2 > 0) {
    AddRegularLattice(Nbox2,Nlattice2,r,box2);

    for (j=0; j<Nbox2; j++) {
      i = Nbox1 + j;
      for (k=0; k<ndim; k++) sph->sphdata[i].r[k] = r[ndim*j + k];
      for (k=0; k<ndim; k++) sph->sphdata[i].v[k] = 0.0;
      sph->sphdata[i].v[0] = vfluid2[0];
      sph->sphdata[i].m = rhofluid2*volume/(float) Nbox2;
      sph->sphdata[i].u = 1.5; //temp0/gammaone/mu_bar;
    }
  }

  delete[] r;

  return;
}


// ============================================================================
// SphSimulation::RandomBox
// ============================================================================
void SphSimulation::RandomBox(void)
{
  float *r;

  debug2("[SphSimulation::RandomBox]");

  sph->AllocateMemory(sph->Nsph);
  r = new float[ndim*sph->Nsph];

  // Add a cube of random particles defined by the simulation bounding box
  AddRandomBox(sph->Nsph,r,simbox);

  // Initialise all other variables
  for (int i=0; i<sph->Nsph; i++) {
    for (int k=0; k<ndim; k++) {
      sph->sphdata[i].r[k] = r[ndim*i + k];
      sph->sphdata[i].v[k] = 0.0f;
      sph->sphdata[i].a[k] = 0.0f;
    }
    sph->sphdata[i].m = 1.0f / (float) sph->Nsph;
    sph->sphdata[i].invomega = 0.5f;
    sph->sphdata[i].iorig = i;
  }

  delete[] r;

  return;
}


// ============================================================================
// SphSimulation::AddRandomBox
// ============================================================================
void SphSimulation::AddRandomBox(int Npart, float *r, DomainBox box)
{
  debug2("[SphSimulation::AddRandomBox]");

  for (int i=0; i<Npart; i++) {
    for (int k=0; k<ndim; k++) {
      r[ndim*i + k] = box.boxmin[k] + (box.boxmax[k] - box.boxmin[k])*
	(float)(rand()%RAND_MAX)/(float)RAND_MAX;
    }
  }

  return;
}



// ============================================================================
// SphSimulation::AddRegularLattice
// ============================================================================
void SphSimulation::AddRegularLattice(int Npart, int Nlattice[ndimmax], 
				      float *r, DomainBox box)
{
  int i;
  int ii;

  debug2("[SphSimulation::AddRegularLattice]");
  
  if (ndim == 1) {
    for (ii=0; ii<Nlattice[0]; ii++) {
      i = ii;
      r[i] = box.boxmin[0] + ((float)ii + 0.5)*
	(box.boxmax[0] - box.boxmin[0])/(float)Nlattice[0];
    }
  }

  return;
}



