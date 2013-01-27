// ============================================================================
// SphSimulation.cpp
// ============================================================================


#include <iostream>
#include <string>
#include <cstdio>
#include <cstring>
#include "SphSimulation.h"
#include "Parameters.h"
#include "Debug.h"
using namespace std;



// ============================================================================
// SphSimulation::SphSimulation
// ============================================================================
SphSimulation::SphSimulation()
{
  n = 0;
  Nsteps = 0;
  t = 0.0f;
}



// ============================================================================
// SphSimulation::~SphSimulation
// ============================================================================
SphSimulation::~SphSimulation()
{
}



// ============================================================================
// SphSimulation::Run
// ============================================================================
void SphSimulation::Run(int Nstepsmax, double tmax)
{

  // --------------------------------------------------------------------------
  do {

    MainLoop();

  } while (t < tmax && Nsteps < Nstepsmax);
  // --------------------------------------------------------------------------

  return;
}



// ============================================================================
// SphSimulation::GenerateIC
// ============================================================================
void SphSimulation::GenerateIC(int N) 
{
  debug2("[SphSimulation::GenerateIC]\n");

  /*sph = new GradhSph;
  sph->kern = new m4;

  sph->kern->Setup(ndim);
  sph->Nsph = N;
  
  sph->AllocateMemory(N);
  srand(1);*/

  sph->RandomBox();

  sph->InitialSmoothingLengthGuess();
  //CalculateSphProperties();

  printf("Finished generating random box\n");

  return;
}



// ============================================================================
// SphSimulation::ComputeBlockTimesteps
// ============================================================================
void SphSimulation::ComputeBlockTimesteps(void)
{
  double dt;

  debug2("[SphSimulation::ComputeBlockTimesteps]\n");

  timestep = big_number;

  for (int i=0; i<sph->Nsph; i++) {
    dt = sphint->Timestep(sph->sphdata[i],sph->eos);
    if (dt < timestep) timestep = dt;
  }

  return;
}



// ============================================================================
// SphSimulation::Setup
// Main function for setting up a new SPH simulation.
// ============================================================================
void SphSimulation::Setup(void)
{
  debug1("[SphSimulation::Setup]\n");

  // Set-up all parameters and assign default values
  simparams.SetDefaultValues();

  // Read parameters files assigning any contained variables
  simparams.ReadParamsFile(paramfile);

  // Assign dimensionality variables here (for now)
#if !defined(FIXED_DIMENSIONS)
  ndim = simparams.intparams["ndim"];
  vdim = simparams.intparams["ndim"];
  bdim = simparams.intparams["ndim"];
#endif

  // Create SPH object based on chosen method in params file
  if (simparams.stringparams["sph"] == "gradh") 
    sph = new GradhSph(ndim,vdim,bdim);
  else {
    cout << "Unrecognised parameter : " << endl; exit(0);
  }

  // Create kernel object based on params file
  if (simparams.stringparams["kernel"] == "m4") {
    sph->kern = new M4Kernel(ndim);
    //sph->kern->Setup(ndim);
  }
  else {
    cout << "Unrecognised parameter : " << endl; exit(0);
  }

  // Create neighbour searching object based on chosen method in params file
  if (simparams.stringparams["neib_search"] == "bruteforce")
    sphneib = new BruteForceSearch;
  else {
    cout << "Unrecognised parameter : " << endl; exit(0);
  }

  if (simparams.stringparams["sph_integration"] == "lfkdk")
    sphint = new SphLFKDK;
  else {
    cout << "Unrecognised parameter : " << endl; exit(0);
  }


  if (simparams.stringparams["gas_eos"] == "isothermal") sph->eos = 
    new Isothermal(simparams.floatparams["temp0"],
		   simparams.floatparams["mu_bar"],
		   simparams.floatparams["gamma_eos"]);
  else {
    cout << "Unrecognised parameter : " << endl; exit(0);
  }

  sph->Nsph = simparams.intparams["Npart"];


  // Generate initial conditions
  if (simparams.stringparams["ic"] == "random_cube")
    sph->RandomBox();
  else {
    cout << "Unrecognised parameter : " << endl; exit(0);
  }



  // Set initial smoothing lengths and create initial ghost particles
  // --------------------------------------------------------------------------
  if (sph->Nsph > 0) {
    
    sph->InitialSmoothingLengthGuess();

    // Reorder particles

    // Search ghost particles

    // Update neighbour tree

  }


  // Compute all SPH particle properties (if SPH particles exist)
  // --------------------------------------------------------------------------
  if (sph->Nsph > 0) {

    // Calculate all SPH properties
    sphneib->UpdateAllSphProperties(sph,simparams);

    // Copy data to ghosts

    // Zero accelerations (perhaps)
    for (int i=0; i<sph->Nsph; i++) {
      for (int k=0; k<ndim; k++) sph->sphdata[i].a[k] = 0.0;
      sph->sphdata[i].dudt = 0.0;
    }

    // Calculate all hydro forces
    sphneib->UpdateAllSphForces(sph,simparams);

    // Calculate all gravitational forces

  }

  // Set r0,v0,a0 for initial step
  sphint->EndTimestep(n,sph->Nsph,sph->sphdata);
  

  return;
}



// ============================================================================
// SphSimulation::MainLoop
// ============================================================================
void SphSimulation::MainLoop(void)
{
  debug1("[SphSimulation::MainLoop]\n");

  // Compute timesteps for all particles
  ComputeBlockTimesteps();

  // Advance time variables
  n = n + 1;
  Nsteps = Nsteps + 1;
  t = t + timestep;

  // Advance SPH particles positions and velocities
  sphint->AdvanceParticles(sph->Nsph,sph->sphdata,timestep);

  // --------------------------------------------------------------------------
  if (sph->Nsph > 0) {
    
    // Reorder particles

    // Search ghost particles

    // Update neighbour tree

  }


  // --------------------------------------------------------------------------
  if (sph->Nsph > 0) {

    // Calculate all SPH properties
    sphneib->UpdateAllSphProperties(sph,simparams);

    // Copy data to ghosts

    // Zero accelerations (perhaps)
    for (int i=0; i<sph->Nsph; i++) 
      for (int k=0; k<ndim; k++) sph->sphdata[i].a[k] = 0.0;

    // Calculate all hydro forces
    sphneib->UpdateAllSphForces(sph,simparams);

    // Calculate all gravitational forces

  }


  // Apply correction steps
  sphint->CorrectionTerms(sph->Nsph,sph->sphdata,timestep);

  // End-of-step
  sphint->EndTimestep(n,sph->Nsph,sph->sphdata);


  return;
}


