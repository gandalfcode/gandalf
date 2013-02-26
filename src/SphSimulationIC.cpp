// ============================================================================
// SphSimulationIC.cpp
// ============================================================================


#include <iostream>
#include <string>
#include <cstdio>
#include <cstring>
#include <math.h>
#include "Precision.h"
#include "Exception.h"
#include "SphSimulation.h"
#include "Sph.h"
#include "Parameters.h"
#include "InlineFuncs.h"
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
  FLOAT volume;
  FLOAT *r;
  DomainBox box1;
  DomainBox box2;
  int Nlattice1[ndimmax];
  int Nlattice2[ndimmax];
  FLOAT vfluid1[ndimmax];
  FLOAT vfluid2[ndimmax];
  FLOAT rhofluid1 = simparams.floatparams["rhofluid1"];
  FLOAT rhofluid2 = simparams.floatparams["rhofluid2"];
  FLOAT press1 = simparams.floatparams["press1"];
  FLOAT press2 = simparams.floatparams["press2"];
  FLOAT temp0 = simparams.floatparams["temp0"];
  FLOAT mu_bar = simparams.floatparams["mu_bar"];
  FLOAT gammaone = simparams.floatparams["gamma_eos"] - 1.0;
  Nlattice1[0] = simparams.intparams["Nlattice1[0]"];
  Nlattice2[0] = simparams.intparams["Nlattice2[0]"];
  vfluid1[0] = simparams.floatparams["vfluid1[0]"];
  vfluid2[0] = simparams.floatparams["vfluid2[0]"];

  debug2("[SphSimulation::ShockTube]");

  if (ndim != 1) {
    cout << "Wrong dimensionality : " << ndim << endl;
    exit(0);
  }

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
  r = new FLOAT[ndim*sph->Nsph];
  cout << "Allocating memory : " << sph->Nsph << endl;


  // Add particles for LHS of the shocktube
  // --------------------------------------------------------------------------
  if (Nbox1 > 0) {
    AddRegularLattice(Nbox1,Nlattice1,r,box1);

    for (i=0; i<Nbox1; i++) {
      for (k=0; k<ndim; k++) sph->sphdata[i].r[k] = r[ndim*i + k];
      for (k=0; k<ndim; k++) sph->sphdata[i].v[k] = 0.0;
      sph->sphdata[i].v[0] = vfluid1[0];
      sph->sphdata[i].m = rhofluid1*volume/(FLOAT) Nbox1;
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
      sph->sphdata[i].m = rhofluid2*volume/(FLOAT) Nbox2;
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
  FLOAT *r;

  debug2("[SphSimulation::RandomBox]");

  sph->AllocateMemory(sph->Nsph);
  r = new FLOAT[ndim*sph->Nsph];

  // Add a cube of random particles defined by the simulation bounding box
  AddRandomBox(sph->Nsph,r,simbox);

  // Initialise all other variables
  for (int i=0; i<sph->Nsph; i++) {
    for (int k=0; k<ndim; k++) {
      sph->sphdata[i].r[k] = r[ndim*i + k];
      sph->sphdata[i].v[k] = 0.0f;
      sph->sphdata[i].a[k] = 0.0f;
    }
    sph->sphdata[i].m = 1.0f / (FLOAT) sph->Nsph;
    sph->sphdata[i].invomega = 0.5f;
    sph->sphdata[i].iorig = i;
    sph->sphdata[i].u = 1.5;
  }

  delete[] r;

  return;
}



// ============================================================================
// SphSimulation::LatticeBox
// ============================================================================
void SphSimulation::LatticeBox(void)
{
  FLOAT *r;
  int Nlattice1[ndimmax];
  Nlattice1[0] = simparams.intparams["Nlattice1[0]"];
  Nlattice1[1] = simparams.intparams["Nlattice1[1]"];
  Nlattice1[2] = simparams.intparams["Nlattice1[1]"];

  debug2("[SphSimulation::RandomBox]");

  if (ndim == 1) sph->Nsph = Nlattice1[0];
  else if (ndim == 2) sph->Nsph = Nlattice1[0]*Nlattice1[1];
  else if (ndim == 3) sph->Nsph = Nlattice1[0]*Nlattice1[1]*Nlattice1[2];
  sph->AllocateMemory(sph->Nsph);
  r = new FLOAT[ndim*sph->Nsph];

  // Add ..
  AddHexagonalLattice(sph->Nsph,Nlattice1,r,simbox);

  // Initialise all other variables
  for (int i=0; i<sph->Nsph; i++) {
    for (int k=0; k<ndim; k++) {
      sph->sphdata[i].r[k] = r[ndim*i + k];
      sph->sphdata[i].v[k] = 0.0f;
      sph->sphdata[i].a[k] = 0.0f;
    }
    sph->sphdata[i].m = 1.0f / (FLOAT) sph->Nsph;
    sph->sphdata[i].invomega = 0.5f;
    sph->sphdata[i].iorig = i;
    sph->sphdata[i].u = 1.5;
  }

  delete[] r;

  return;
}



// ============================================================================
// SphSimulation::RandomSphere
// ============================================================================
void SphSimulation::RandomSphere(void)
{
  FLOAT *r;
  FLOAT rcentre[ndimmax];
  FLOAT radius = 1.0;

  debug2("[SphSimulation::RandomBox]");

  sph->AllocateMemory(sph->Nsph);
  r = new FLOAT[ndim*sph->Nsph];

  for (int k=0; k<ndim; k++) rcentre[k] = 0.0;

  // Add a cube of random particles defined by the simulation bounding box
  AddRandomSphere(sph->Nsph,r,rcentre,radius);

  // Initialise all other variables
  for (int i=0; i<sph->Nsph; i++) {
    for (int k=0; k<ndim; k++) {
      sph->sphdata[i].r[k] = r[ndim*i + k];
      sph->sphdata[i].v[k] = 0.0f;
      sph->sphdata[i].a[k] = 0.0f;
    }
    sph->sphdata[i].m = 1.0f / (FLOAT) sph->Nsph;
    sph->sphdata[i].invomega = 1.0;
    sph->sphdata[i].zeta = 0.0;
    sph->sphdata[i].iorig = i;
  }

  delete[] r;

  return;
}



// ============================================================================
// SphSimulation::KHI
// ============================================================================
void SphSimulation::KHI(void)
{
  int i;
  int j;
  int k;
  int Nbox1;
  int Nbox2;
  FLOAT volume;
  FLOAT *r;
  DomainBox box1;
  DomainBox box2;
  int Nlattice1[ndimmax];
  int Nlattice2[ndimmax];
  FLOAT vfluid1[ndimmax];
  FLOAT vfluid2[ndimmax];
  FLOAT rhofluid1 = simparams.floatparams["rhofluid1"];
  FLOAT rhofluid2 = simparams.floatparams["rhofluid2"];
  FLOAT press1 = simparams.floatparams["press1"];
  FLOAT press2 = simparams.floatparams["press2"];
  FLOAT temp0 = simparams.floatparams["temp0"];
  FLOAT mu_bar = simparams.floatparams["mu_bar"];
  FLOAT gammaone = simparams.floatparams["gamma_eos"] - 1.0;
  FLOAT amp = simparams.floatparams["amp"];
  FLOAT lambda = simparams.floatparams["lambda"];
  Nlattice1[0] = simparams.intparams["Nlattice1[0]"];
  Nlattice1[1] = simparams.intparams["Nlattice1[1]"];
  Nlattice2[0] = simparams.intparams["Nlattice2[0]"];
  Nlattice2[1] = simparams.intparams["Nlattice2[1]"];
  vfluid1[0] = simparams.floatparams["vfluid1[0]"];
  vfluid2[0] = simparams.floatparams["vfluid2[0]"];

  debug2("[SphSimulation::ShockTube]");

  if (ndim != 2) {
    string message = "Kelvin-Helmholtz instability only in 2D";
    ExceptionHandler::getIstance().raise(message);
  }

  // Compute size and range of fluid bounding boxes
  // --------------------------------------------------------------------------
  box1.boxmin[0] = simbox.boxmin[0];
  box1.boxmax[0] = simbox.boxmax[0];
  box1.boxmin[1] = simbox.boxmin[1];
  box1.boxmax[1] = simbox.boxmin[1] + simbox.boxhalf[1];
  box2.boxmin[0] = simbox.boxmin[0];
  box2.boxmax[0] = simbox.boxmax[0];
  box2.boxmin[1] = simbox.boxmin[1] + simbox.boxhalf[1];
  box2.boxmax[1] = simbox.boxmax[1];

  volume = (box1.boxmax[0] - box1.boxmin[0])*(box1.boxmax[1] - box1.boxmin[1]);
  Nbox1 = Nlattice1[0]*Nlattice1[1];
  Nbox2 = Nlattice2[0]*Nlattice2[1];


  // Allocate local and main particle memory
  sph->Nsph = Nbox1 + Nbox2;
  sph->AllocateMemory(sph->Nsph);
  r = new FLOAT[ndim*sph->Nsph];
  cout << "Nbox1 : " << Nbox1 << "    Nbox2 : " << Nbox2 << endl;
  cout << "Allocating memory : " << sph->Nsph << endl;


  // Add particles for LHS of the shocktube
  // --------------------------------------------------------------------------
  if (Nbox1 > 0) {
    AddRegularLattice(Nbox1,Nlattice1,r,box1);

    for (i=0; i<Nbox1; i++) {
      for (k=0; k<ndim; k++) sph->sphdata[i].r[k] = r[ndim*i + k];
      for (k=0; k<ndim; k++) sph->sphdata[i].v[k] = 0.0;
      sph->sphdata[i].r[1] -= 0.25*simbox.boxsize[1];
      if (sph->sphdata[i].r[1] < simbox.boxmin[1]) 
	sph->sphdata[i].r[1] += simbox.boxsize[1];
      sph->sphdata[i].v[0] = vfluid1[0];
      sph->sphdata[i].m = rhofluid1*volume/(FLOAT) Nbox1;
      sph->sphdata[i].u = press1/rhofluid1/gammaone;
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
      sph->sphdata[i].r[1] -= 0.25*simbox.boxsize[1];
      if (sph->sphdata[i].r[1] < simbox.boxmin[1]) 
	sph->sphdata[i].r[1] += simbox.boxsize[1];
      sph->sphdata[i].v[0] = vfluid2[0];
      sph->sphdata[i].m = rhofluid2*volume/(FLOAT) Nbox2;
      sph->sphdata[i].u = press2/rhofluid2/gammaone;
    }
  }

  // Add velocity perturbation here
  // --------------------------------------------------------------------------
  FLOAT sigmapert = 0.05/sqrt(2.0);
  for (i=0; i<sph->Nsph; i++) {
    sph->sphdata[i].v[1] = amp*sin(2.0*pi*sph->sphdata[i].r[0]/lambda)*
      (exp(-pow(sph->sphdata[i].r[1] + 0.25,2)/2.0/sigmapert/sigmapert) +  
       exp(-pow(sph->sphdata[i].r[1] - 0.25,2)/2.0/sigmapert/sigmapert));
  }

  // Set initial smoothing lengths and create initial ghost particles
  // --------------------------------------------------------------------------
  sph->Nghost = 0;
  sph->Nghostmax = sph->Nsphmax - sph->Nsph;
  sph->Ntot = sph->Nsph;
  for (int i=0; i<sph->Nsph; i++) sph->sphdata[i].active = true;
  
  sph->InitialSmoothingLengthGuess();
  sphneib->UpdateTree(sph,simparams);
  
  // Search ghost particles
  SearchGhostParticles();

  sphneib->UpdateAllSphProperties(sph);

  // Update neighbour tre
  sphneib->UpdateTree(sph,simparams);

  // Calculate all SPH properties
  sphneib->UpdateAllSphProperties(sph);
  
  for (i=0; i<sph->Nsph; i++) 
    sph->sphdata[i].u = press1/sph->sphdata[i].rho/gammaone;

  delete[] r;

  return;
}



// ============================================================================
// SphSimulation::AddRandomBox
// ============================================================================
void SphSimulation::AddRandomBox(int Npart, FLOAT *r, DomainBox box)
{
  debug2("[SphSimulation::AddRandomBox]");

  for (int i=0; i<Npart; i++) {
    for (int k=0; k<ndim; k++) {
      r[ndim*i + k] = box.boxmin[k] + (box.boxmax[k] - box.boxmin[k])*
	(FLOAT)(rand()%RAND_MAX)/(FLOAT)RAND_MAX;
    }
  }

  return;
}



// ============================================================================
// SphSimulation::AddRandomsphere
// ============================================================================
void SphSimulation::AddRandomSphere(int Npart, FLOAT *r, 
				    FLOAT *rcentre, FLOAT radius)
{
  FLOAT rad;
  FLOAT rpos[ndimmax];

  debug2("[SphSimulation::AddRandomSphere]");

  for (int i=0; i<Npart; i++) {
    do {
      for (int k=0; k<ndim; k++) 
	rpos[k] = 1.0 - 2.0*(FLOAT)(rand()%RAND_MAX)/(FLOAT)RAND_MAX;
      rad = DotProduct(rpos,rpos,ndim);
    } while (rad > radius);
    for (int k=0; k<ndim; k++) r[ndim*i + k] = rcentre[k] + rpos[k];
  }

  return;
}



// ============================================================================
// SphSimulation::AddRegularLattice
// ============================================================================
void SphSimulation::AddRegularLattice(int Npart, int Nlattice[ndimmax], 
				      FLOAT *r, DomainBox box)
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
      r[i] = box.boxmin[0] + ((FLOAT)ii + 0.5)*
	(box.boxmax[0] - box.boxmin[0])/(FLOAT)Nlattice[0];
    }
  }
  // --------------------------------------------------------------------------
  else if (ndim == 2) {
    for (jj=0; jj<Nlattice[1]; jj++) {
      for (ii=0; ii<Nlattice[0]; ii++) {
	i = jj*Nlattice[0] + ii;
	r[ndim*i] = box.boxmin[0] + ((FLOAT)ii + 0.5)*
	  (box.boxmax[0] - box.boxmin[0])/(FLOAT)Nlattice[0];
	r[ndim*i + 1] = box.boxmin[1] + ((FLOAT)jj + 0.5)*
	  (box.boxmax[1] - box.boxmin[1])/(FLOAT)Nlattice[1];
      }
    }
  }
  // --------------------------------------------------------------------------
  else if (ndim == 3) {
    for (kk=0; kk<Nlattice[1]; kk++) {
      for (jj=0; jj<Nlattice[1]; jj++) {
	for (ii=0; ii<Nlattice[0]; ii++) {
	  i = kk*Nlattice[0]*Nlattice[1] + jj*Nlattice[0] + ii;
	  r[ndim*i] = box.boxmin[0] + ((FLOAT)ii + 0.5)*
	    (box.boxmax[0] - box.boxmin[0])/(FLOAT)Nlattice[0];
	  r[ndim*i + 1] = box.boxmin[1] + ((FLOAT)jj + 0.5)*
	    (box.boxmax[1] - box.boxmin[1])/(FLOAT)Nlattice[1];
	  r[ndim*i + 2] = box.boxmin[2] + ((FLOAT)kk + 0.5)*
	    (box.boxmax[2] - box.boxmin[2])/(FLOAT)Nlattice[2];
	}
      }
    }
  }

  return;
}



// ============================================================================
// SphSimulation::AddHexagonalLattice
// Create simple hexagonal-packed lattice using A-B-A-B pattern in z-direction
// N.B. the box is scaled to fit to the x-boxsize
// ============================================================================
void SphSimulation::AddHexagonalLattice(int Npart, int Nlattice[ndimmax], 
					FLOAT *r, DomainBox box)
{
  int i;
  int ii;
  int jj;
  int kk;
  FLOAT rad;

  debug2("[SphSimulation::AddHexagonalLattice]");
  
  // Calculate 'radius' of points for simpler calculation
  rad = 0.5*(box.boxmax[0] - box.boxmin[1])/(FLOAT) Nlattice[0];

  // Create lattice depending on dimensionality
  // ==========================================================================
  if (ndim == 1) {
    for (ii=0; ii<Nlattice[0]; ii++) {
      i = ii;
      r[i] = box.boxmin[0] + rad + 2.0*(FLOAT)ii*rad;
    }
  }

  // ==========================================================================
  else if (ndim == 2) {
    for (jj=0; jj<Nlattice[1]; jj++) {
      for (ii=0; ii<Nlattice[0]; ii++) {
	i = jj*Nlattice[0] + ii;
	r[ndim*i] = box.boxmin[0] + 0.5*rad + 
	  (2.0*(FLOAT)ii + (FLOAT)(jj%2))*rad;
	r[ndim*i + 1] = box.boxmin[1] + 0.5*sqrt(3.0)*rad + 
	  (FLOAT)jj*sqrt(3.0)*rad;
      }
    }
  }

  // ==========================================================================
  else if (ndim == 3) {
    for (kk=0; kk<Nlattice[2]; kk++) {
      for (jj=0; jj<Nlattice[1]; jj++) {
	for (ii=0; ii<Nlattice[0]; ii++) {
	  i = kk*Nlattice[0]*Nlattice[1] + jj*Nlattice[0] + ii;
	  r[ndim*i] = box.boxmin[0] + 0.5*rad + 
	    (2.0*(FLOAT)ii + (FLOAT)(jj%2) + (FLOAT)((kk+1)%2))*rad;
	  r[ndim*i + 1] = box.boxmin[1] + 0.5*sqrt(3.0)*rad + 
	    (FLOAT)jj*sqrt(3.0)*rad + (FLOAT)(kk%2)/sqrt(3.0);
	  r[ndim*i + 2] = box.boxmin[1] + sqrt(6.0)*rad/3.0 + 
	    (FLOAT)kk*2.0*sqrt(6.0)*rad/3.0;
	}
      }
    }
  }

  return;
}




// ============================================================================
// SphSimulation::CutSphere
// Cut-out a sphere containing exactly'Nsphere' particles (if selected)
// ============================================================================
int SphSimulation::CutSphere(int Nsphere, int Npart,
			     FLOAT radsphere, FLOAT *r, 
			     DomainBox box, bool exact)
{
  int i;
  int k;
  int niterations = 0;
  int Ninterior = 0;
  int Niteration = 0;
  int Nitmax = 100;
  FLOAT dr[ndimmax];
  FLOAT drsqd;
  FLOAT r_low = 0.0;
  FLOAT r_high;
  FLOAT radius;
  FLOAT rcentre[ndimmax];

  debug2("[SphSimulation::CutSphere]");

  // Find centre and shortest edge-length of bounding box
  r_high = (FLOAT) big_number;
  for (k=0; k<ndim; k++) {
    rcentre[k] = 0.5*(box.boxmin[k] + box.boxmax[k]);
    r_high = min(r_high,(FLOAT)0.5*(box.boxmax[k] - box.boxmin[k]));
  }

  // Bisection iteration to determine the radius containing the desired 
  // number of particles
  // --------------------------------------------------------------------------
  do {
    radius = 0.5*(r_low + r_high);
    Ninterior = 0;

    for (i=0; i<Npart; i++) {
      for (k=0; k<ndim; k++) dr[k] = r[ndim*i + k] - rcentre[k];
      drsqd = DotProduct(dr,dr,ndim);
      if (drsqd <= radius*radius) Ninterior++;
    }
    if (Ninterior > Nsphere) r_high = radius;
    if (Ninterior < Nsphere) r_low = radius;
    if (Niteration++ > Nitmax);
  } while (Ninterior != Nsphere);


  // Now radius containing require number has been identified, record only 
  // particle inside sphere
  // --------------------------------------------------------------------------
  Ninterior = 0;
  for (i=0; i<Npart; i++) {
    for (k=0; k<ndim; k++) dr[k] = r[ndim*i + k] - rcentre[k];
    drsqd = DotProduct(dr,dr,ndim);
    if (drsqd <= radius*radius) {
      for (k=0; k<ndim; k++) r[ndim*Ninterior + k] = r[ndim*i + k];
      Ninterior++;
    }
  }

  return Ninterior;
}
