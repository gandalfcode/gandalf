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
  float rhofluid1 = simparams.floatparams["rhofluid1"];
  float rhofluid2 = simparams.floatparams["rhofluid2"];
  float press1 = simparams.floatparams["press1"];
  float press2 = simparams.floatparams["press2"];
  float temp0 = simparams.floatparams["temp0"];
  float mu_bar = simparams.floatparams["mu_bar"];
  float gammaone = simparams.floatparams["gamma_eos"] - 1.0;
  Nlattice1[0] = simparams.intparams["Nlattice1[0]"];
  Nlattice2[0] = simparams.intparams["Nlattice2[0]"];
  vfluid1[0] = simparams.floatparams["vfluid1[0]"];
  vfluid2[0] = simparams.floatparams["vfluid2[0]"];

  debug2("[SphSimulation::ShockTube]");

  // Compute size and range of fluid bounding boxes
  // --------------------------------------------------------------------------
  if (ndim == 1) {
    box1.boxmin[0] = simbox.boxmin[0];
    box1.boxmax[0] = 0.0;
    box2.boxmin[0] = 0.0;
    box2.boxmax[0] = simbox.boxmax[0];
    volume = box1.boxmax[0] - box1.boxmin[0];
    Nbox1 = Nlattice1[0];
    Nbox2 = Nlattice2[0];
  }

  // Allocate local and main particle memory
  sph->Nsph = Nbox1 + Nbox2;
  sph->AllocateMemory(sph->Nsph);
  r = new float[ndim*sph->Nsph];
  cout << "Allocating memory : " << sph->Nsph << endl;


  // Add particles for LHS of the shocktube
  // --------------------------------------------------------------------------
  if (Nbox1 > 0) {
    AddRegularLattice(Nbox1,Nlattice1,r,box1);

    for (i=0; i<Nbox1; i++) {
      for (k=0; k<ndim; k++) sph->sphdata[i].r[k] = r[ndim*i + k];
      for (k=0; k<ndim; k++) sph->sphdata[i].v[k] = 0.0;
      sph->sphdata[i].v[0] = vfluid1[0];
      sph->sphdata[i].m = rhofluid1*volume/(float) Nbox1;
      sph->sphdata[i].u = temp0/gammaone/mu_bar;
    }
  }

  // Add particles for RHS of the shocktube
  // --------------------------------------------------------------------------
  if (Nbox2 > 0) {
    AddRegularLattice(Nbox2,Nlattice2,r,box2);

    for (j=0; j<Nbox2; j++) {
      i = Nbox1 + j;
      for (k=0; k<ndim; k++) sph->sphdata[i].r[k] = r[ndim*j + k];
      for (k=0; k<ndim; k++) sph->sphdata[i].v[k] = 0.0;
      sph->sphdata[i].v[0] = vfluid2[0];
      sph->sphdata[i].m = rhofluid2*volume/(float) Nbox2;
      sph->sphdata[i].u = temp0/gammaone/mu_bar;
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
  int jj;
  int kk;

  debug2("[SphSimulation::AddRegularLattice]");
  
  // Create lattice depending on dimensionality
  // --------------------------------------------------------------------------
  if (ndim == 1) {
    for (ii=0; ii<Nlattice[0]; ii++) {
      i = ii;
      r[i] = box.boxmin[0] + ((float)ii + 0.5)*
	(box.boxmax[0] - box.boxmin[0])/(float)Nlattice[0];
    }
  }
  // --------------------------------------------------------------------------
  else if (ndim == 2) {
    for (jj=0; jj<Nlattice[1]; jj++) {
      for (ii=0; ii<Nlattice[0]; ii++) {
	i = jj*Nlattice[0] + ii;
	r[ndim*i] = box.boxmin[0] + ((float)ii + 0.5)*
	  (box.boxmax[0] - box.boxmin[0])/(float)Nlattice[0];
	r[ndim*i + 1] = box.boxmin[1] + ((float)jj + 0.5)*
	  (box.boxmax[1] - box.boxmin[1])/(float)Nlattice[1];
      }
    }
  // --------------------------------------------------------------------------
  else if (ndim == 3) {
    for (kk=0; kk<Nlattice[1]; kk++) {
      for (jj=0; jj<Nlattice[1]; jj++) {
	for (ii=0; ii<Nlattice[0]; ii++) {
	  i = kk*Nlattice[0]*Nlattice[1] + jj*Nlattice[0] + ii;
	  r[ndim*i] = box.boxmin[0] + ((float)ii + 0.5)*
	    (box.boxmax[0] - box.boxmin[0])/(float)Nlattice[0];
	  r[ndim*i + 1] = box.boxmin[1] + ((float)jj + 0.5)*
	    (box.boxmax[1] - box.boxmin[1])/(float)Nlattice[1];
	  r[ndim*i + 2] = box.boxmin[2] + ((float)jj + 0.5)*
	    (box.boxmax[2] - box.boxmin[2])/(float)Nlattice[2];
	}
      }
    }
  }

  return;
}



