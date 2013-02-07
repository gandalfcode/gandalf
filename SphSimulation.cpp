// ============================================================================
// SphSimulation.cpp
// Contains all main functions controlling the SPH simulation work-flow.
// ============================================================================


#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <math.h>
#include <cstdio>
#include <cstring>
#include "Exception.h"
#include "SphSimulation.h"
#include "Parameters.h"
#include "InlineFuncs.h"
#include "Debug.h"
using namespace std;



// ============================================================================
// SphSimulation::SphSimulation
// ============================================================================
SphSimulation::SphSimulation()
{
  paramfile = "adshock.dat"; //"params.dat";
  n = 0;
  Nsteps = 0;
  t = 0.0;
  setup = false;
}



// ============================================================================
// SphSimulation::~SphSimulation
// ============================================================================
SphSimulation::~SphSimulation()
{
}



// ============================================================================
// SphSimulation::Run
// Controls the simulation main loop, including exit conditions.
// If provided, will only advance the simulation by 'Nadvance' steps.
// ============================================================================
void SphSimulation::Run(int Nadvance)
{
  int Ntarget;                              // Selected integer timestep

  debug1("[SphSimulation::Run]");

  // Set integer timestep exit condition if provided as parameter.
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

  CalculateDiagnostics();
  cout << "Eerror : " << fabs(diag0.Etot - diag.Etot)/fabs(diag0.Etot) << endl;

  return;
}



// ============================================================================
// SphSimulation::Output
// Controls when regular output snapshots are written by the code.
// ============================================================================
void SphSimulation::Output(void)
{
  string filename;
  string nostring;
  stringstream ss;
  string fileform = simparams.stringparams["out_file_form"];

  debug2("[SphSimulation::Output]");

  if (Nsteps%noutputstep == 0) cout << "t : " << t << "    Nsteps : " 
				    << Nsteps << endl;

  // Output a data snapshot if reached required time
  if (t >= tsnapnext) {
    Noutsnap++;
    tsnapnext += dt_snap;
    nostring = "";
    ss << setfill('0') << setw(5) << Noutsnap;
    nostring = ss.str();
    filename = run_id + '.' + fileform + '.' + nostring;
    ss.str(std::string());
    WriteSnapshotFile(filename,"column");
  }

  return;
}



// ============================================================================
// SphSimulation::CalculateDiagnostics
// Calculates all diagnostic quantities (e.g. conserved quantities), 
// saves to the diagnostic data structure and outputs to screen.
// ============================================================================
void SphSimulation::CalculateDiagnostics(void)
{
  int k;

  debug2("[SphSimulation::CalculateDiagnostics]");

  diag.Etot = 0.0;
  diag.utot = 0.0;
  diag.ketot = 0.0;
  diag.gpetot = 0.0;
  for (k=0; k<ndim; k++) diag.mom[k] = 0.0;
  for (k=0; k<3; k++) diag.angmom[k] = 0.0;
  for (k=0; k<ndim; k++) diag.force[k] = 0.0;
  for (k=0; k<ndim; k++) diag.force_grav[k] = 0.0;

  for (int i=0; i<sph->Nsph; i++) {
    diag.ketot += sph->sphdata[i].m*
      DotProduct(sph->sphdata[i].v,sph->sphdata[i].v,ndim);
    diag.utot += sph->sphdata[i].m*sph->sphdata[i].u;
    diag.gpetot += sph->sphdata[i].m*sph->sphdata[i].gpot;
    for (k=0; k<ndim; k++) {
      diag.mom[k] += sph->sphdata[i].m*sph->sphdata[i].v[k];
      diag.force[k] += sph->sphdata[i].m*sph->sphdata[i].a[k];
      diag.force_grav[k] += sph->sphdata[i].m*sph->sphdata[i].agrav[k];
    }
  }
  diag.ketot *= 0.5;
  diag.gpetot *= 0.5;
  diag.Etot = diag.ketot + diag.utot + diag.gpetot;

  cout << "Printing out diagnostics" << endl;
  cout << "Etot       : " << diag.Etot << endl;
  cout << "utot       : " << diag.utot << endl;
  cout << "ketot      : " << diag.ketot << endl;
  cout << "gpetot     : " << diag.gpetot << endl;
  if (ndim == 1) {
    cout << "mom        : " << diag.mom[0] << endl;
    cout << "force      : " << diag.force[0] << endl;
    cout << "force_grav : " << diag.force_grav[0] << endl;
  }
  else if (ndim == 2) {
    cout << "mom        : " << diag.mom[0] << "   " << diag.mom[1] << endl;
    cout << "force      : " << diag.force[0] << "   " << diag.force[1] << endl;
    cout << "force_grav : " << diag.force_grav[0] << "   " 
	 << diag.force_grav[1] << endl;
  }
  else if (ndim == 3) {
    cout << "mom        : " << diag.mom[0] << "   " 
	 << diag.mom[1] << "   " << diag.mom[2] << endl;
    cout << "force      : " << diag.force[0] << "   " 
	 << diag.force[1] << "   " << diag.force[2] << endl;
    cout << "force_grav : " << diag.force_grav[0] << "   " 
	 << diag.force_grav[1] << "   " << diag.force_grav[2] << endl;
  }

  return;
}



// ============================================================================
// SphSimulation::GenerateIC
// Generate initial conditions for SPH simulation chosen in parameters file.
// ============================================================================
void SphSimulation::GenerateIC(void) 
{
  debug2("[SphSimulation::GenerateIC]");

  // Generate initial conditions
  if (simparams.stringparams["ic"] == "random_cube") 
    RandomBox();
  else if (simparams.stringparams["ic"] == "random_sphere") 
    RandomSphere();
  else if (simparams.stringparams["ic"] == "lattice_cube")
    LatticeBox();
  else if (simparams.stringparams["ic"] == "shocktube") 
    ShockTube();
  else if (simparams.stringparams["ic"] == "khi") 
    KHI();
  else {
    string message = "Unrecognised parameter : ic = " + simparams.stringparams["ic"];
    ExceptionHandler::getIstance().raise(message);
  }

  return;
}



// ============================================================================
// SphSimulation::ComputeBlockTimesteps
// Computes global timestep for SPH simulation.
// ============================================================================
void SphSimulation::ComputeBlockTimesteps(void)
{
  double dt;                                // Aux. timestep variable

  debug2("[SphSimulation::ComputeBlockTimesteps]");

  timestep = big_number;

  for (int i=0; i<sph->Nsph; i++) {
    dt = sphint->Timestep(sph->sphdata[i],sph->hydro_forces);
    if (dt < timestep) timestep = dt;
    if (simparams.stringparams["gas_eos"] == "energy_eqn") {
      dt = uint->Timestep(sph->sphdata[i]);
      if (dt < timestep) timestep = dt;
    }      
  }

  //cout << "Global timestep : " << timestep << "   t : " << t << endl;

  return;
}



// ============================================================================
// SphSimulation::ProcessParameters
// Process all the options chosen in the parameters file, setting various 
// simulation variables and creating important simulation objects.
// ============================================================================
void SphSimulation::ProcessParameters(void)
{
  map<string, int> &intparams = simparams.intparams;
  map<string, float> &floatparams = simparams.floatparams;
  map<string, string> &stringparams = simparams.stringparams;

  debug2("[SphSimulation::ProcessParameters]");

  // Assign dimensionality variables here (for now)
#if !defined(FIXED_DIMENSIONS)
  ndim = intparams["ndim"];
  vdim = intparams["ndim"];
  bdim = intparams["ndim"];
#endif

  // Create SPH object based on chosen method in params file
  // --------------------------------------------------------------------------
  if (stringparams["sph"] == "gradh") {
    sph = new GradhSph(ndim,vdim,bdim);
    sph->alpha_visc = floatparams["alpha_visc"];
    sph->beta_visc = floatparams["beta_visc"];
  }
  else {
    string message = "Unrecognised parameter : sph = " + simparams.stringparams["sph"];
    ExceptionHandler::getIstance().raise(message);
  }

  // Create kernel object based on params file
  if (stringparams["kernel"] == "m4") {
    sph->kern = new M4Kernel(ndim);
  }
  else if (stringparams["kernel"] == "quintic") {
    sph->kern = new QuinticKernel(ndim);
  }
  else {
    string message = "Unrecognised parameter : kernel = " + 
      simparams.stringparams["kernel"];
    ExceptionHandler::getIstance().raise(message);
  }

  // Boundary condition variables
  // --------------------------------------------------------------------------
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
  for (int k=0; k<3; k++) {
    simbox.boxsize[k] = simbox.boxmax[k] - simbox.boxmin[k];
    simbox.boxhalf[k] = 0.5*simbox.boxsize[k];
    cout << "SIMBOX : " << k << "  " << simbox.boxsize[k] << "   " 
	 << simbox.boxhalf[k] << "   " << simbox.boxmin[k] << endl;
  }

  // Create neighbour searching object based on chosen method in params file
  // --------------------------------------------------------------------------
  if (stringparams["neib_search"] == "bruteforce")
    sphneib = new BruteForceSearch(ndim);
  else if (stringparams["neib_search"] == "grid")
    sphneib = new GridSearch(ndim);
  else {
    string message = "Unrecognised parameter : neib_search = " + simparams.stringparams["neib_search"];
    ExceptionHandler::getIstance().raise(message);
  }

  if (stringparams["sph_integration"] == "lfkdk") {
    sphint = new SphLFKDK(ndim,vdim,floatparams["accel_mult"],
			  floatparams["courant_mult"]);
  }
  else {
    string message = "Unrecognised parameter : sph_integration = " + simparams.stringparams["sph_integration"];
    ExceptionHandler::getIstance().raise(message);
  }

  // Thermal physics options
  // --------------------------------------------------------------------------
  if (stringparams["gas_eos"] == "energy_eqn") {
    sph->eos = new Adiabatic(floatparams["temp0"],
			     floatparams["mu_bar"],
			     floatparams["gamma_eos"]);
    if (stringparams["energy_integration"] == "PEC") {
      uint = new EnergyPEC(floatparams["energy_mult"]);
    }
    else {
      string message = "Unrecognised parameter : energy_integration = " + simparams.stringparams["energy_integration"];
      ExceptionHandler::getIstance().raise(message);
    }
  }
  else if (stringparams["gas_eos"] == "isothermal") 
    sph->eos = new Isothermal(floatparams["temp0"],
			      floatparams["mu_bar"],
			      floatparams["gamma_eos"]);
  else {
    string message = "Unrecognised parameter : gas_eos = " + simparams.stringparams["gas_eos"];
    ExceptionHandler::getIstance().raise(message);
  }


  sph->Nsph = intparams["Npart"];
  sph->h_fac = floatparams["h_fac"];
  sph->h_converge = floatparams["h_converge"];
  sph->hydro_forces = intparams["hydro_forces"];
  sph->self_gravity = intparams["self_gravity"];
  sph->avisc = stringparams["avisc"];
  sph->acond = stringparams["acond"];
  Nstepsmax = intparams["Nstepsmax"];
  run_id = stringparams["run_id"];
  tend = floatparams["tend"];
  dt_snap = floatparams["dt_snap"];
  noutputstep = intparams["noutputstep"];

  return;
}



// ============================================================================
// SphSimulation::Setup
// Main function for setting up a new SPH simulation.
// ============================================================================
void SphSimulation::Setup(void)
{
  debug1("[SphSimulation::Setup]");

  // Set-up all parameters and assign default values
  simparams.SetDefaultValues();

  // Read parameters files assigning any contained variables
  simparams.ReadParamsFile(paramfile);

  // Process the parameters file setting up all simulation objects
  ProcessParameters();

  // Generate initial conditions for simulation
  GenerateIC();

  // Set time variables here (for now)
  Noutsnap = 0;
  tsnapnext = dt_snap;


  // Set initial smoothing lengths and create initial ghost particles
  // --------------------------------------------------------------------------
  if (sph->Nsph > 0) {

    // Set all relevant particle counters
    sph->Nghost = 0;
    sph->Nghostmax = sph->Nsphmax - sph->Nsph;
    sph->Ntot = sph->Nsph;
    for (int i=0; i<sph->Nsph; i++) sph->sphdata[i].active = true;
    
    sph->InitialSmoothingLengthGuess();
    sphneib->UpdateTree(sph,simparams);

    sphneib->UpdateAllSphProperties(sph,simparams);

    // Search ghost particles
    SearchGhostParticles();

    // Update neighbour tre
    sphneib->UpdateTree(sph,simparams);
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
      for (int k=0; k<ndim; k++) sph->sphdata[i].agrav[k] = 0.0;
      sph->sphdata[i].gpot = 0.0;
      sph->sphdata[i].dudt = 0.0;
    }

    // Calculate all hydro forces
    sphneib->UpdateAllSphForces(sph,simparams);

    // Add accelerations
    for (int i=0; i<sph->Nsph; i++) {
      for (int k=0; k<ndim; k++) 
	sph->sphdata[i].a[k] += sph->sphdata[i].agrav[k];
    }

  }

  // Set r0,v0,a0 for initial step
  sphint->EndTimestep(n,sph->Nsph,sph->sphdata);
  if (simparams.stringparams["gas_eos"] == "energy_eqn")
    uint->EndTimestep(n,sph->Nsph,sph->sphdata);
  

  CalculateDiagnostics();
  diag0 = diag;

  setup = true;

  return;
}



// ============================================================================
// SphSimulation::MainLoop
// ============================================================================
void SphSimulation::MainLoop(void)
{
  debug2("[SphSimulation::MainLoop]");

  // Compute timesteps for all particles
  ComputeBlockTimesteps();

  // Advance time variables
  n = n + 1;
  Nsteps = Nsteps + 1;
  t = t + timestep;

  // Advance SPH particles positions and velocities
  sphint->AdvanceParticles(sph->Nsph,sph->sphdata,timestep);
  if (simparams.stringparams["gas_eos"] == "energy_eqn")
    uint->EnergyIntegration(sph->Nsph,sph->sphdata,timestep);

  // Check all boundary conditions
  CheckBoundaries();

  // --------------------------------------------------------------------------
  if (sph->Nsph > 0) {
    
    // Reorder particles

    // Search ghost particles
    SearchGhostParticles();

    // Update neighbour tree
    sphneib->UpdateTree(sph,simparams);
  }


  // --------------------------------------------------------------------------
  if (sph->Nsph > 0) {

    // Calculate all SPH properties
    sphneib->UpdateAllSphProperties(sph,simparams);

    // Copy data to ghosts
    CopyDataToGhosts();

    // Zero accelerations (perhaps)
    for (int i=0; i<sph->Nsph; i++) {
      if (sph->sphdata[i].active) {
	for (int k=0; k<ndim; k++) sph->sphdata[i].a[k] = 0.0;
	for (int k=0; k<ndim; k++) sph->sphdata[i].agrav[k] = 0.0;
	sph->sphdata[i].gpot = 0.0;
	sph->sphdata[i].dudt = 0.0;
      }
    }

    // Calculate all hydro and gravitational forces
    sphneib->UpdateAllSphForces(sph,simparams);

    // Add accelerations
    for (int i=0; i<sph->Nsph; i++) {
      for (int k=0; k<ndim; k++) 
	sph->sphdata[i].a[k] += sph->sphdata[i].agrav[k];
    }

  }


  // Apply correction steps
  sphint->CorrectionTerms(sph->Nsph,sph->sphdata,timestep);
  if (simparams.stringparams["gas_eos"] == "energy_eqn")
    uint->EnergyCorrectionTerms(sph->Nsph,sph->sphdata,timestep);

  // End-of-step
  sphint->EndTimestep(n,sph->Nsph,sph->sphdata);
  if (simparams.stringparams["gas_eos"] == "energy_eqn")
    uint->EndTimestep(n,sph->Nsph,sph->sphdata);

  return;
}
