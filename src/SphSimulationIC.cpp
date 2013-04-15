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
template <int ndim>
void SphSimulation<ndim>::ShockTube(void)
{
  int i;
  int j;
  int k;
  int Nbox1;
  int Nbox2;
  FLOAT volume;
  FLOAT *r;
  DomainBox<ndim> box1;
  DomainBox<ndim> box2;
  int Nlattice1[ndim];
  int Nlattice2[ndim];
  FLOAT vfluid1[ndim];
  FLOAT vfluid2[ndim];
  FLOAT rhofluid1 = simparams->floatparams["rhofluid1"];
  FLOAT rhofluid2 = simparams->floatparams["rhofluid2"];
  FLOAT press1 = simparams->floatparams["press1"];
  FLOAT press2 = simparams->floatparams["press2"];
  FLOAT temp0 = simparams->floatparams["temp0"];
  FLOAT mu_bar = simparams->floatparams["mu_bar"];
  FLOAT gammaone = simparams->floatparams["gamma_eos"] - 1.0;
  Nlattice1[0] = simparams->intparams["Nlattice1[0]"];
  Nlattice2[0] = simparams->intparams["Nlattice2[0]"];
  vfluid1[0] = simparams->floatparams["vfluid1[0]"];
  vfluid2[0] = simparams->floatparams["vfluid2[0]"];

  debug2("[SphSimulation::ShockTube]");

  if (ndim != 1) {
    string message = "Wrong dimensionality : " + ndim;
    ExceptionHandler::getIstance().raise(message);
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
      if (sph->gas_eos == "isothermal")
    	sph->sphdata[i].u = temp0/gammaone/mu_bar;
      else
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
      sph->sphdata[i].v[0] = vfluid2[0];
      sph->sphdata[i].m = rhofluid2*volume/(FLOAT) Nbox2;
      if (sph->gas_eos == "isothermal")
	    sph->sphdata[i].u = temp0/gammaone/mu_bar;
      else
	    sph->sphdata[i].u = press2/rhofluid2/gammaone;
    }
  }

  delete[] r;

  return;
}


// ============================================================================
// SphSimulation::RandomBox
// ============================================================================
template <int ndim>
void SphSimulation<ndim>::RandomBox(void)
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
    //sph->sphdata[i].m = (1.0 + 9.0*sph->sphdata[i].r[0])*1.0f / (FLOAT) sph->Nsph;
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
template <int ndim>
void SphSimulation<ndim>::LatticeBox(void)
{
  FLOAT *r;
  int Nlattice1[ndim];
  Nlattice1[0] = simparams->intparams["Nlattice1[0]"];
  Nlattice1[1] = simparams->intparams["Nlattice1[1]"];
  Nlattice1[2] = simparams->intparams["Nlattice1[2]"];

  debug2("[SphSimulation::RandomBox]");

  if (ndim == 1) sph->Nsph = Nlattice1[0];
  else if (ndim == 2) sph->Nsph = Nlattice1[0]*Nlattice1[1];
  else if (ndim == 3) sph->Nsph = Nlattice1[0]*Nlattice1[1]*Nlattice1[2];
  sph->AllocateMemory(sph->Nsph);
  r = new FLOAT[ndim*sph->Nsph];

  // Add ..
  AddRegularLattice(sph->Nsph,Nlattice1,r,simbox);

  // Initialise all other variables
  for (int i=0; i<sph->Nsph; i++) {
    for (int k=0; k<ndim; k++) {
      sph->sphdata[i].r[k] = r[ndim*i + k];
      sph->sphdata[i].v[k] = 0.0f;
      sph->sphdata[i].a[k] = 0.0f;
    }
    //sph->sphdata[i].m = (1.0 + 9.0*sph->sphdata[i].r[0])*1.0f / (FLOAT) sph->Nsph;

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
template <int ndim>
void SphSimulation<ndim>::RandomSphere(void)
{
  FLOAT *r;
  FLOAT rcentre[ndim];
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
// SphSimulation::ContactDiscontinuity
// ============================================================================
template <int ndim>
void SphSimulation<ndim>::ContactDiscontinuity(void)
{
  DomainBox<ndim> box1;
  DomainBox<ndim> box2;
  int Nlattice1[3];
  int Nlattice2[3];
  int i;
  int j;
  int k;
  int Nbox1;
  int Nbox2;
  FLOAT volume;
  FLOAT *r;
  FLOAT vfluid1[3];
  FLOAT vfluid2[3];
  FLOAT rhofluid1 = simparams->floatparams["rhofluid1"];
  FLOAT rhofluid2 = simparams->floatparams["rhofluid2"];
  FLOAT press1 = simparams->floatparams["press1"];
  FLOAT press2 = simparams->floatparams["press2"];
  FLOAT temp0 = simparams->floatparams["temp0"];
  FLOAT mu_bar = simparams->floatparams["mu_bar"];
  FLOAT gammaone = simparams->floatparams["gamma_eos"] - 1.0;
  FLOAT amp = simparams->floatparams["amp"];
  FLOAT lambda = simparams->floatparams["lambda"];
  Nlattice1[0] = simparams->intparams["Nlattice1[0]"];
  Nlattice1[1] = simparams->intparams["Nlattice1[1]"];
  Nlattice2[0] = simparams->intparams["Nlattice2[0]"];
  Nlattice2[1] = simparams->intparams["Nlattice2[1]"];
  vfluid1[0] = simparams->floatparams["vfluid1[0]"];
  vfluid2[0] = simparams->floatparams["vfluid2[0]"];

  debug2("[SphSimulation::ContactDiscontinuity]");


  // 1d simulation
  // ==========================================================================
  if (ndim == 1) {
	box1.boxmin[0] = simbox.boxmin[0];
	box1.boxmax[0] = 0.8*simbox.boxmax[0];
	box2.boxmin[0] = 0.8*simbox.boxmax[0];
	box2.boxmax[0] = simbox.boxmax[0];
	volume = box1.boxmax[0] - box1.boxmin[0];
	Nbox1 = Nlattice1[0];
	Nbox2 = Nlattice2[0];

	// Allocate local and main particle memory
	sph->Nsph = Nbox1 + Nbox2;
	sph->AllocateMemory(sph->Nsph);
	r = new FLOAT[ndim*sph->Nsph];
	cout << "Allocating memory : " << sph->Nsph << endl;


	// --------------------------------------------------------------------
	if (Nbox1 > 0) {
	  AddRegularLattice(Nbox1,Nlattice1,r,box1);
      volume = box1.boxmax[0] - box1.boxmin[0];
      cout << "Vol1 : " << volume << "   m1 : " 
	   << rhofluid1*volume/(FLOAT) Nbox1 << "   rho : "
	   << rhofluid1*volume << "   Nbox : " << Nbox1 << endl;
      for (i=0; i<Nbox1; i++) {
		sph->sphdata[i].r[0] = r[i] - 0.4*simbox.boxsize[0];
		if (sph->sphdata[i].r[0] < simbox.boxmin[0])
		  sph->sphdata[i].r[0] += simbox.boxsize[0];
		sph->sphdata[i].v[0] = 0.0;
        sph->sphdata[i].m = rhofluid1*volume/(FLOAT) Nbox1;
        if (sph->gas_eos == "isothermal")
	      sph->sphdata[i].u = temp0/gammaone/mu_bar;
	    else
	      sph->sphdata[i].u = press1/rhofluid1/gammaone;
	  }
	}

	// --------------------------------------------------------------------
	if (Nbox2 > 0) {
	  AddRegularLattice(Nbox2,Nlattice2,r,box2);
      volume = box2.boxmax[0] - box2.boxmin[0];
      cout << "Vol2 : " << volume << "   m2 : " 
	   << rhofluid2*volume/(FLOAT) Nbox2 
	   << "   rho : " << rhofluid2*volume << "   Nbox : " << Nbox2 << endl;
	  for (j=0; j<Nbox2; j++) {
        i = Nbox1 + j;
		sph->sphdata[i].r[0] = r[j] - 0.4*simbox.boxsize[0];
		if (sph->sphdata[i].r[0] < simbox.boxmin[0])
		  sph->sphdata[i].r[0] += simbox.boxsize[0];
		sph->sphdata[i].v[0] = 0.0;
        sph->sphdata[i].m = rhofluid2*volume/(FLOAT) Nbox2;
        if (sph->gas_eos == "isothermal")
	      sph->sphdata[i].u = temp0/gammaone/mu_bar;
	    else
	      sph->sphdata[i].u = press2/rhofluid2/gammaone;
	  }
	}

  }
  // ==========================================================================
  else if (ndim == 2) {



  }
  // ==========================================================================


  // Set initial smoothing lengths and create initial ghost particles
  // --------------------------------------------------------------------------
  sph->Nghost = 0;
  sph->Nghostmax = sph->Nsphmax - sph->Nsph;
  sph->Ntot = sph->Nsph;
  for (int i=0; i<sph->Nsph; i++) sph->sphdata[i].active = true;

  sph->InitialSmoothingLengthGuess();
  sphneib->UpdateTree(sph,*simparams);

  // Search ghost particles
  SearchGhostParticles();

  sphneib->UpdateAllSphProperties(sph);

  // Update neighbour tre
  sphneib->UpdateTree(sph,*simparams);

  // Calculate all SPH properties
  sphneib->UpdateAllSphProperties(sph);

  sphneib->UpdateTree(sph,*simparams);
  sphneib->UpdateAllSphProperties(sph);

  CopySphDataToGhosts();

  //for (i=0; i<sph->Nsph; i++) {
  //  sph->sphdata[i].u = press1/sph->sphdata[i].rho/gammaone;
  // cout << "WTF?? : " << i << "   x : " << sph->sphdata[i].rho << "   " << gammaone*sph->sphdata[i].u*sph->sphdata[i].rho << endl;
  //}

  // Calculate all SPH properties
  sphneib->UpdateAllSphProperties(sph);

  delete[] r;

  return;
}



// ============================================================================
// SphSimulation::KHI
// ============================================================================
template <int ndim>
void SphSimulation<ndim>::KHI(void)
{
  int i;
  int j;
  int k;
  int Nbox1;
  int Nbox2;
  FLOAT volume;
  FLOAT *r;
  DomainBox<ndim> box1;
  DomainBox<ndim> box2;
  int Nlattice1[ndim];
  int Nlattice2[ndim];
  FLOAT vfluid1[ndim];
  FLOAT vfluid2[ndim];
  FLOAT rhofluid1 = simparams->floatparams["rhofluid1"];
  FLOAT rhofluid2 = simparams->floatparams["rhofluid2"];
  FLOAT press1 = simparams->floatparams["press1"];
  FLOAT press2 = simparams->floatparams["press2"];
  FLOAT temp0 = simparams->floatparams["temp0"];
  FLOAT mu_bar = simparams->floatparams["mu_bar"];
  FLOAT gammaone = simparams->floatparams["gamma_eos"] - 1.0;
  FLOAT amp = simparams->floatparams["amp"];
  FLOAT lambda = simparams->floatparams["lambda"];
  Nlattice1[0] = simparams->intparams["Nlattice1[0]"];
  Nlattice1[1] = simparams->intparams["Nlattice1[1]"];
  Nlattice2[0] = simparams->intparams["Nlattice2[0]"];
  Nlattice2[1] = simparams->intparams["Nlattice2[1]"];
  vfluid1[0] = simparams->floatparams["vfluid1[0]"];
  vfluid2[0] = simparams->floatparams["vfluid2[0]"];

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
  sphneib->UpdateTree(sph,*simparams);
  
  // Search ghost particles
  SearchGhostParticles();

  sphneib->UpdateAllSphProperties(sph);

  // Update neighbour tre
  sphneib->UpdateTree(sph,*simparams);

  // Calculate all SPH properties
  sphneib->UpdateAllSphProperties(sph);
  
  for (i=0; i<sph->Nsph; i++) 
    sph->sphdata[i].u = press1/sph->sphdata[i].rho/gammaone;

  delete[] r;

  return;
}



// ============================================================================
// SphSimulation::SedovBlastWave
// ============================================================================
template <int ndim>
void SphSimulation<ndim>::SedovBlastWave(void)
{
  int i;                                // ..
  int j;                                // ..
  int k;                                // ..
  int Nbox;                             // ..
  int Ncold;                            // ..
  int Nhot;                             // ..
  int Nlattice[3];                      // ..
  int *hotlist;                         // ..
  FLOAT drmag;                          // ..
  FLOAT drsqd;                          // ..
  FLOAT mbox;                           // ..
  FLOAT r_hot;                          // ..
  FLOAT ufrac;                          // ..
  FLOAT umax;                           // ..
  FLOAT utot;                           // ..
  FLOAT volume;                         // ..
  FLOAT *r;                             // ..

  // Create local copies of initial conditions parameters
  int smooth_ic = simparams->intparams["smooth_ic"];
  FLOAT rhofluid = simparams->floatparams["rhofluid1"];
  FLOAT press = simparams->floatparams["press1"];
  FLOAT mu_bar = simparams->floatparams["mu_bar"];
  FLOAT gammaone = simparams->floatparams["gamma_eos"] - 1.0;
  FLOAT kefrac = simparams->floatparams["kefrac"];
  Nlattice[0] = simparams->intparams["Nlattice1[0]"];
  Nlattice[1] = simparams->intparams["Nlattice1[1]"];
  Nlattice[2] = simparams->intparams["Nlattice1[2]"];

  debug2("[SphSimulation::SedovBlastWave]");


  // Compute size and range of fluid bounding boxes
  // --------------------------------------------------------------------------
  if (ndim == 1) {
    volume = simbox.boxmax[0] - simbox.boxmin[0];
    Nbox = Nlattice[0];
  }
  else if (ndim == 2) {
    volume = (simbox.boxmax[0] - simbox.boxmin[0])*
      (simbox.boxmax[1] - simbox.boxmin[1]);
    Nbox = Nlattice[0]*Nlattice[1];
  }
  else if (ndim == 3) {
    volume = (simbox.boxmax[0] - simbox.boxmin[0])*
      (simbox.boxmax[1] - simbox.boxmin[1])*
      (simbox.boxmax[2] - simbox.boxmin[2]);
    Nbox = Nlattice[0]*Nlattice[1]*Nlattice[2];
  }
  mbox = volume*rhofluid;
  ufrac = max((FLOAT) 0.0,(FLOAT) 1.0 - kefrac);
  Ncold = 0;
  Nhot = 0;
  r_hot = powf(powf(4.0,ndim)/(FLOAT) Nbox,sph->invndim);

  // Allocate local and main particle memory
  sph->Nsph = Nbox;
  sph->AllocateMemory(sph->Nsph);
  r = new FLOAT[ndim*sph->Nsph];
  hotlist = new int[sph->Nsph];
  cout << "Here??? : " << endl;
  cout << "mbox : " << mbox << "   r_hot : " << r_hot 
       << "   Nbox : " << Nbox << endl;

  // Add particles for LHS of the shocktube
  // --------------------------------------------------------------------------
  if (Nbox > 0) {
    AddRegularLattice(Nbox,Nlattice,r,simbox);

    for (i=0; i<Nbox; i++) {
      for (k=0; k<ndim; k++) sph->sphdata[i].r[k] = r[ndim*i + k];
      for (k=0; k<ndim; k++) sph->sphdata[i].v[k] = 0.0;
      sph->sphdata[i].m = mbox/(FLOAT) Nbox;
      sph->sphdata[i].u = small_number;
    }
  }

  // Set initial smoothing lengths and create initial ghost particles
  // --------------------------------------------------------------------------
  sph->Nghost = 0;
  sph->Nghostmax = sph->Nsphmax - sph->Nsph;
  sph->Ntot = sph->Nsph;
  for (i=0; i<sph->Nsph; i++) sph->sphdata[i].active = true;
  
  sph->InitialSmoothingLengthGuess();
  sphneib->UpdateTree(sph,*simparams);
  
  // Search ghost particles
  SearchGhostParticles();

  sphneib->UpdateAllSphProperties(sph);

  // Update neighbour tre
  sphneib->UpdateTree(sph,*simparams);

  // Calculate all SPH properties
  sphneib->UpdateAllSphProperties(sph);

  // Now calculate which particles are hot
  // --------------------------------------------------------------------------
  umax = (FLOAT) 0.0;
  utot = (FLOAT) 0.0;
  for (i=0; i<sph->Nsph; i++) {
    drsqd = DotProduct(sph->sphdata[i].r,sph->sphdata[i].r,ndim);
    if (drsqd < r_hot*r_hot) {
      hotlist[i] = 1;
      if (smooth_ic == 1)
	sph->sphdata[i].u = sph->sphdata[i].m*
	  sph->kernp->w0(sph->kernp->kernrange*sqrt(drsqd)/r_hot);
      else
	sph->sphdata[i].u = sph->sphdata[i].m;
      utot += sph->sphdata[i].u;
      umax = max(umax,sph->sphdata[i].u);
      Nhot++;
    }
    else {
      hotlist[i] = 0;
      Ncold++;
    }
  }

  cout << "utot : " << utot << endl;

  // Normalise the energies
  // --------------------------------------------------------------------------
  for (i=0; i<sph->Nsph; i++) {
    if (hotlist[i] == 1) {
      drmag = sqrt(DotProduct(sph->sphdata[i].r,sph->sphdata[i].r,ndim));
      sph->sphdata[i].u = sph->sphdata[i].u/utot/sph->sphdata[i].m;
      //,1.0e-6*umax/sph->sphdata[i].m);
      for (k=0; k<ndim; k++) sph->sphdata[i].v[k] = 
	sqrt(2.0*kefrac*sph->sphdata[i].u)*
	sph->sphdata[i].r[k]/(drmag + small_number);
      sph->sphdata[i].u = ufrac*sph->sphdata[i].u;
      //,1.0e-6*umax/sph->sphdata[i].m);
    }
    else {
      sph->sphdata[i].u = 1.0e-6/sph->sphdata[i].m;
      for (k=0; k<ndim; k++) sph->sphdata[i].v[k] = (FLOAT) 0.0;
    }
  }

  delete[] hotlist;
  delete[] r;

  return;
}




// ============================================================================
// SphSimulation::SoundWave
// ============================================================================
template <int ndim>
void SphSimulation<ndim>::SoundWave(void)
{
  int i,k;
  FLOAT csound,diff,lambda,kwave,omegawave,ugas,volume,xold,xnew;
  FLOAT *r;

  int Nlattice1[ndim];
  int Npart = simparams->intparams["Npart"];
  FLOAT rhofluid1 = simparams->floatparams["rhofluid1"];
  FLOAT press1 = simparams->floatparams["press1"];
  FLOAT gamma = simparams->floatparams["gamma_eos"];
  FLOAT gammaone = gamma - 1.0;
  FLOAT amp = simparams->floatparams["amp"];
  FLOAT temp0 = simparams->floatparams["temp0"];
  FLOAT mu_bar = simparams->floatparams["mu_bar"];
  Nlattice1[0] = simparams->intparams["Nlattice1[0]"];

  debug2("[SphSimulation::SoundWave]");

  if (ndim != 1) {
	  string message="Sound wave only available in 1D";
	  ExceptionHandler::getIstance().raise(message);
  }
  cout << "gas_eos : " << sph->gas_eos << endl;
  if (sph->gas_eos == "isothermal") {
    ugas = temp0/gammaone/mu_bar;
    press1 = gammaone*rhofluid1*ugas;
    csound = sqrt(press1/rhofluid1);
  }
  else {
    ugas = press1/rhofluid1/gammaone;
    csound = sqrt(gamma*press1/rhofluid1);
  }
  lambda = simbox.boxmax[0] - simbox.boxmin[0];
  kwave = twopi/lambda;
  omegawave = twopi*csound/lambda;

  cout << "amp : " << amp << "    rhofluid1" << rhofluid1 << endl;
  cout << "kwave : " << kwave << "    csound : " << csound << endl;

   // Allocate local and main particle memory
   sph->Nsph = Npart;
   Nlattice1[0] = Npart;
   sph->AllocateMemory(sph->Nsph);
   r = new FLOAT[ndim*sph->Nsph];
   cout << "Allocating memory : " << sph->Nsph << endl;

   AddRegularLattice(Npart,Nlattice1,r,simbox);

   for (i=0; i<Npart; i++) {
       xnew = r[ndim*i];

       // Solve iterative procedure for particle positions in sound wave
       do {
    	 //break;
    	 xold = xnew;
    	 //xnew = lambda*(xold - amp*cos(kwave*xold))/
    	 //		 (lambda - amp*cos(kwave*lambda));
    	 //xnew = (lambda - amp*cos(kwave*lambda))*(r[ndim*i]/lambda) + amp*cos(kwave*xnew);
    	 xnew = r[ndim*i] - amp*(1.0 - cos(kwave*xnew))/kwave;
    	 diff = fabs((xnew - xold)/lambda);
       } while(diff > 1.0e-6);
       if (xnew < simbox.boxmin[0] || xnew > simbox.boxmax[0]) 
	 cout << "WTF?? : " << xnew << endl;
       if (xnew > simbox.boxmax[0]) xnew -= simbox.boxsize[0];
       if (xnew < simbox.boxmin[0]) xnew += simbox.boxsize[0];
       for (k=0; k<ndim; k++) sph->sphdata[i].r[k] = xnew;
       for (k=0; k<ndim; k++) 
	 sph->sphdata[i].v[k] = csound*amp*sin(kwave*xnew);
       sph->sphdata[i].m = rhofluid1*lambda/(FLOAT) Npart;
       if (sph->gas_eos == "isothermal")
      	 sph->sphdata[i].u = temp0/gammaone/mu_bar;
       else
         sph->sphdata[i].u = press1/rhofluid1/gammaone;
     }


  return;
}



// ============================================================================
// SphSimulation::AddRandomBox
// ============================================================================
template <int ndim>
void SphSimulation<ndim>::AddRandomBox(int Npart, FLOAT *r, DomainBox<ndim> box)
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
template <int ndim>
void SphSimulation<ndim>::AddRandomSphere(int Npart, FLOAT *r,
				    FLOAT *rcentre, FLOAT radius)
{
  FLOAT rad;
  FLOAT rpos[ndim];

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
template <int ndim>
void SphSimulation<ndim>::AddRegularLattice(int Npart, int Nlattice[3],
					    FLOAT *r, DomainBox<ndim> box)
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
    for (kk=0; kk<Nlattice[2]; kk++) {
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
template <int ndim>
void SphSimulation<ndim>::AddHexagonalLattice(int Npart, int Nlattice[ndim],
					FLOAT *r, DomainBox<ndim> box)
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
template <int ndim>
int SphSimulation<ndim>::CutSphere(int Nsphere, int Npart,
			     FLOAT radsphere, FLOAT *r, 
			     DomainBox<ndim> box, bool exact)
{
  int i;
  int k;
  int niterations = 0;
  int Ninterior = 0;
  int Niteration = 0;
  int Nitmax = 100;
  FLOAT dr[ndim];
  FLOAT drsqd;
  FLOAT r_low = 0.0;
  FLOAT r_high;
  FLOAT radius;
  FLOAT rcentre[ndim];

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

