// ============================================================================
// SphSimulation.cpp
// ============================================================================


#include <iostream>
#include <sstream>
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
  paramfile = "params.dat";
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
void SphSimulation::Run(int Nadvance)
{
  int Ntarget;

  if (Nadvance < 0) Ntarget = Nstepsmax;
  else Ntarget = Nsteps + Nadvance;

  // Continue to run simulation until we reach the required time, or 
  // exeeded the maximum allowed number of steps.
  // --------------------------------------------------------------------------
  while (t < tend && Nsteps < Ntarget) {

    MainLoop();
    Output();

  }
  // --------------------------------------------------------------------------

  return;
}



// ============================================================================
// SphSimulation::AdvanceSteps
// ============================================================================
void SphSimulation::AdvanceSteps(int Nadvance)
{

  // Advance the simulation by 'Nadvance' integer steps.
  // --------------------------------------------------------------------------
  for (int i=0; i<Nadvance; i++) {

    MainLoop();
    Output();

  }
  // --------------------------------------------------------------------------

  return;
}



// ============================================================================
// SphSimulation::Output
// ============================================================================
void SphSimulation::Output(void)
{
  string filename;
  string nostring;
  stringstream ss;

  // Output a data snapshot if reached required time
  if (t >= tsnapnext) {
    Noutsnap++;
    tsnapnext += dt_snap;
    nostring = "";
    ss << Noutsnap;
    nostring = ss.str();
    filename = run_id + '.' + nostring;
    ss.str(std::string());
    WriteSnapshotFile(filename,"column");
  }

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

  //sph->RandomBox();

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
    dt = sphint->Timestep(sph->sphdata[i]);
    if (dt < timestep) timestep = dt;
  }

  cout << "Global timestep : " << timestep << endl;

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

  map<string, int> &intparams = simparams.intparams;
  map<string, float> &floatparams = simparams.floatparams;
  map<string, string> &stringparams = simparams.stringparams;

  // Assign dimensionality variables here (for now)

#if !defined(FIXED_DIMENSIONS)
  ndim = intparams["ndim"];
  vdim = intparams["ndim"];
  bdim = intparams["ndim"];
#endif

  // Create SPH object based on chosen method in params file
  if (stringparams["sph"] == "gradh") {
    sph = new GradhSph(ndim,vdim,bdim);
    sph->alpha_visc = floatparams["alpha_visc"];
    sph->beta_visc = floatparams["beta_visc"];
  }
  else {
    cout << "Unrecognised parameter : " << endl; exit(0);
  }

  // Create kernel object based on params file
  if (stringparams["kernel"] == "m4") {
    sph->kern = new M4Kernel(ndim);
    //sph->kern->Setup(ndim);
  }
  else {
    cout << "Unrecognised parameter : " << endl; exit(0);
  }

  // Boundary condition variables
  simbox.x_boundary_lhs = stringparams["x_boundary_lhs"];
  simbox.x_boundary_rhs = stringparams["x_boundary_rhs"];
  simbox.y_boundary_lhs = stringparams["y_boundary_lhs"];
  simbox.y_boundary_rhs = stringparams["y_boundary_rhs"];
  simbox.z_boundary_lhs = stringparams["z_boundary_lhs"];
  simbox.z_boundary_rhs = stringparams["z_boundary_rhs"];
  simbox.boxmin[0] = floatparams["boxmin[0]"];
  simbox.boxmin[1] = floatparams["boxmin[1]"];
  simbox.boxmin[2] = floatparams["boxmin[2]"];
  simbox.boxmax[0] = floatparams["boxmax[0]"];
  simbox.boxmax[1] = floatparams["boxmax[1]"];
  simbox.boxmax[2] = floatparams["boxmax[2]"];
  for (int k=0; k<3; k++) 
    simbox.boxsize[k] = simbox.boxmax[k] - simbox.boxmin[k];

  // Create neighbour searching object based on chosen method in params file
  if (stringparams["neib_search"] == "bruteforce")
    sphneib = new BruteForceSearch;
  else {
    cout << "Unrecognised parameter : " << endl; exit(0);
  }

  if (stringparams["sph_integration"] == "lfkdk") {
    sphint = new SphLFKDK(floatparams["accel_mult"],
			  floatparams["courant_mult"]);
  }
  else {
    cout << "Unrecognised parameter : " << endl; exit(0);
  }


  if (stringparams["gas_eos"] == "isothermal") sph->eos =
    new Isothermal(floatparams["temp0"],
		   floatparams["mu_bar"],
		   floatparams["gamma_eos"]);
  else {
    cout << "Unrecognised parameter : " << endl; exit(0);
  }

  sph->Nsph = intparams["Npart"];


  // Generate initial conditions
  if (stringparams["ic"] == "random_cube") 
    RandomBox();
  else if (stringparams["ic"] == "shocktube") 
    ShockTube();
  else {
    cout << "Unrecognised parameter : " << endl; exit(0);
  }


  // Set time variables here (for now)
  Noutsnap = 0;
  Nstepsmax = intparams["Nstepsmax"];
  run_id = stringparams["run_id"];
  tend = floatparams["tend"];
  dt_snap = floatparams["dt_snap"];
  tsnapnext = dt_snap;


  // Set initial smoothing lengths and create initial ghost particles
  // --------------------------------------------------------------------------
  if (sph->Nsph > 0) {

    sph->Ntot = sph->Nsph;
    
    sph->InitialSmoothingLengthGuess();

    sphneib->UpdateAllSphProperties(sph,simparams);

    // Search ghost particles
    SearchGhostParticles();

    // Update neighbour tree

  }


  // Compute all SPH particle properties (if SPH particles exist)
  // --------------------------------------------------------------------------
  if (sph->Nsph > 0) {

    cout << "Ntot : " << sph->Ntot << endl;

    // Calculate all SPH properties
    sphneib->UpdateAllSphProperties(sph,simparams);

    // Copy data to ghosts
    CopyDataToGhosts();

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

  // Check all boundary conditions
  CheckBoundaries();

  // --------------------------------------------------------------------------
  if (sph->Nsph > 0) {
    
    // Reorder particles

    // Search ghost particles
    SearchGhostParticles();

    // Update neighbour tree

  }


  // --------------------------------------------------------------------------
  if (sph->Nsph > 0) {

    // Calculate all SPH properties
    sphneib->UpdateAllSphProperties(sph,simparams);

    // Copy data to ghosts
    CopyDataToGhosts();

    // Zero accelerations (perhaps)
    for (int i=0; i<sph->Nsph; i++) {
      for (int k=0; k<ndim; k++) sph->sphdata[i].a[k] = 0.0;
      sph->sphdata[i].dudt = 0.0;
    }

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


