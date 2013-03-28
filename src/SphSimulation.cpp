// ============================================================================
// SphSimulation.cpp
// Contains all main functions controlling the SPH simulation work-flow.
// ============================================================================


#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <math.h>
#include <time.h>
#include <cstdio>
#include <cstring>
#include "Precision.h"
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
  paramfile = "";
  n = 0;
  nresync = 0;
  integration_step = 1;
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
  diag.Eerror = fabs(diag0.Etot - diag.Etot)/fabs(diag0.Etot);
  cout << "Eerror : " << diag.Eerror << endl;


  return;
}



// ============================================================================
// SphSimulation::InteractiveRun
// Controls the simulation main loop, including exit conditions.
// If provided, will only advance the simulation by 'Nadvance' steps.
// ============================================================================
void SphSimulation::InteractiveRun(int Nadvance)
{
  int Ntarget;                              // Selected integer timestep
  clock_t tstart = clock();
  DOUBLE tdiff = 0.0;
  DOUBLE tpython = 8.0;

  debug1("[SphSimulation::Run]");

  // Set integer timestep exit condition if provided as parameter.
  if (Nadvance < 0) Ntarget = Nstepsmax;
  else Ntarget = Nsteps + Nadvance;

  // Continue to run simulation until we reach the required time, or
  // exeeded the maximum allowed number of steps.
  // --------------------------------------------------------------------------
  while (t < tend && Nsteps < Ntarget && tdiff < tpython) {

    MainLoop();
    Output();

    tdiff = (DOUBLE) (clock() - tstart) / (DOUBLE) CLOCKS_PER_SEC;

  }
  // --------------------------------------------------------------------------

  CalculateDiagnostics();
  diag.Eerror = fabs(diag0.Etot - diag.Etot)/fabs(diag0.Etot);

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
    filename = run_id + '.' + out_file_form + '.' + nostring;
    ss.str(std::string());
    WriteSnapshotFile(filename,"column");
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
  else if (simparams.stringparams["ic"] == "soundwave")
    SoundWave();
  else if (simparams.stringparams["ic"] == "khi") 
    KHI();
  else {
    string message = "Unrecognised parameter : ic = " 
      + simparams.stringparams["ic"];
    ExceptionHandler::getIstance().raise(message);
  }

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

  Sph::aviscenum avisc;
  Sph::acondenum acond;

  debug2("[SphSimulation::ProcessParameters]");

  // Assign dimensionality variables here (for now)
#if !defined(FIXED_DIMENSIONS)
  ndim = intparams["ndim"];
  vdim = intparams["ndim"];
  bdim = intparams["ndim"];
#endif

  // Set the enum for artificial viscosity
  if (stringparams["avisc"] == "none") {
    avisc = Sph::noneav;
  }
  else if (stringparams["avisc"] == "mon97") {
    avisc = Sph::mon97;
  }
  else {
    string message = "Unrecognised parameter : avisc = " +
        simparams.stringparams["avisc"];
    ExceptionHandler::getIstance().raise(message);
  }

  // Set the enum for artificial conductivity
  if (stringparams["acond"] == "none") {
    acond = Sph::noneac;
  }
  else if (stringparams["acond"] == "wadsley2008") {
    acond = Sph::wadsley2008;
  }
  else if (stringparams["acond"] == "price2008") {
    acond = Sph::price2008;
  }
  else {
    string message = "Unrecognised parameter : acond = " +
        simparams.stringparams["acond"];
    ExceptionHandler::getIstance().raise(message);
  }

  // Create SPH object based on chosen method in params file
  // --------------------------------------------------------------------------
  if (stringparams["sph"] == "gradh") {
    string KernelName = stringparams["kernel"];
    if (stringparams["tabulatedkernel"] == "yes") {
        sph = new GradhSph<TabulatedKernel> (ndim, vdim, bdim,
        		intparams["hydro_forces"], intparams["self_gravity"],
        		floatparams["alpha_visc"], floatparams["beta_visc"],
        		floatparams["h_fac"], floatparams["h_converge"],
        		avisc, acond,
        		stringparams["gas_eos"], KernelName);
          }
    else if (stringparams["tabulatedkernel"] == "no"){
        // Depending on the kernel, instantiate a different GradSph object
        if (KernelName == "m4") {
            sph = new GradhSph<M4Kernel> (ndim, vdim, bdim,
            		intparams["hydro_forces"], intparams["self_gravity"],
            		floatparams["alpha_visc"], floatparams["beta_visc"],
            		floatparams["h_fac"], floatparams["h_converge"],
            		avisc, acond,
            		stringparams["gas_eos"], KernelName);
        }
        else if (KernelName == "quintic") {
            sph = new GradhSph<QuinticKernel> (ndim, vdim, bdim,
            		intparams["hydro_forces"], intparams["self_gravity"],
            		floatparams["alpha_visc"], floatparams["beta_visc"],
            		floatparams["h_fac"], floatparams["h_converge"],
            		avisc, acond,
            		stringparams["gas_eos"], KernelName);
	}
        else if (KernelName == "gaussian") {
            sph = new GradhSph<GaussianKernel> (ndim, vdim, bdim,
            		intparams["hydro_forces"], intparams["self_gravity"],
            		floatparams["alpha_visc"], floatparams["beta_visc"],
            		floatparams["h_fac"], floatparams["h_converge"],
            		avisc, acond,
            		stringparams["gas_eos"], KernelName);
        }
        else {
          string message = "Unrecognised parameter : kernel = " +
            simparams.stringparams["kernel"];
          ExceptionHandler::getIstance().raise(message);
        }
    }
    else {
      string message = "Invalid option for the tabulatedkernel parameter: " +
          stringparams["tabulatedkernel"];
      ExceptionHandler::getIstance().raise(message);
    }
  }
  // --------------------------------------------------------------------------
  else if (stringparams["sph"] == "sm2012") {
    string KernelName = stringparams["kernel"];
    if (stringparams["tabulatedkernel"] == "yes") {
        sph = new SM2012Sph<TabulatedKernel> (ndim, vdim, bdim,
        		intparams["hydro_forces"], intparams["self_gravity"],
        		floatparams["alpha_visc"], floatparams["beta_visc"],
        		floatparams["h_fac"], floatparams["h_converge"],
        		avisc, acond,
        		stringparams["gas_eos"], KernelName);
          }
    else if (stringparams["tabulatedkernel"] == "no"){
        // Depending on the kernel, instantiate a different GradSph object
        if (KernelName == "m4") {
            sph = new SM2012Sph<M4Kernel> (ndim, vdim, bdim,
            		intparams["hydro_forces"], intparams["self_gravity"],
            		floatparams["alpha_visc"], floatparams["beta_visc"],
            		floatparams["h_fac"], floatparams["h_converge"],
            		avisc, acond,
            		stringparams["gas_eos"], KernelName);
        }
        else if (KernelName == "quintic") {
            sph = new SM2012Sph<QuinticKernel> (ndim, vdim, bdim,
            		intparams["hydro_forces"], intparams["self_gravity"],
            		floatparams["alpha_visc"], floatparams["beta_visc"],
            		floatparams["h_fac"], floatparams["h_converge"],
            		avisc, acond,
            		stringparams["gas_eos"], KernelName);
        }
        else if (KernelName == "gaussian") {
            sph = new SM2012Sph<GaussianKernel> (ndim, vdim, bdim,
            		intparams["hydro_forces"], intparams["self_gravity"],
            		floatparams["alpha_visc"], floatparams["beta_visc"],
            		floatparams["h_fac"], floatparams["h_converge"],
            		avisc, acond,
            		stringparams["gas_eos"], KernelName);
        }
        else {
          string message = "Unrecognised parameter : kernel = " +
            simparams.stringparams["kernel"];
          ExceptionHandler::getIstance().raise(message);
        }
    }
    else {
      string message = "Invalid option for the tabulatedkernel parameter: " +
          stringparams["tabulatedkernel"];
      ExceptionHandler::getIstance().raise(message);
    }
  }
  // --------------------------------------------------------------------------
  else if (stringparams["sph"] == "godunov") {
    string KernelName = stringparams["kernel"];
    if (stringparams["tabulatedkernel"] == "yes") {
        sph = new GodunovSph<TabulatedKernel> (ndim, vdim, bdim,
        		intparams["hydro_forces"], intparams["self_gravity"],
        		floatparams["alpha_visc"], floatparams["beta_visc"],
        		floatparams["h_fac"], floatparams["h_converge"],
        		avisc, acond,
        		stringparams["gas_eos"], KernelName);
          }
    else if (stringparams["tabulatedkernel"] == "no"){
        // Depending on the kernel, instantiate a different GradSph object
        if (KernelName == "m4") {
            sph = new GodunovSph<M4Kernel> (ndim, vdim, bdim,
            		intparams["hydro_forces"], intparams["self_gravity"],
            		floatparams["alpha_visc"], floatparams["beta_visc"],
            		floatparams["h_fac"], floatparams["h_converge"],
            		avisc, acond,
            		stringparams["gas_eos"], KernelName);
        }
        else if (KernelName == "quintic") {
            sph = new GodunovSph<QuinticKernel> (ndim, vdim, bdim,
            		intparams["hydro_forces"], intparams["self_gravity"],
            		floatparams["alpha_visc"], floatparams["beta_visc"],
            		floatparams["h_fac"], floatparams["h_converge"],
            		avisc, acond,
            		stringparams["gas_eos"], KernelName);
        }
        else if (KernelName == "gaussian") {
            sph = new GodunovSph<GaussianKernel> (ndim, vdim, bdim,
            		intparams["hydro_forces"], intparams["self_gravity"],
            		floatparams["alpha_visc"], floatparams["beta_visc"],
            		floatparams["h_fac"], floatparams["h_converge"],
            		avisc, acond,
            		stringparams["gas_eos"], KernelName);
        }
        else {
          string message = "Unrecognised parameter : kernel = " +
            simparams.stringparams["kernel"];
          ExceptionHandler::getIstance().raise(message);
        }
    }
    else {
      string message = "Invalid option for the tabulatedkernel parameter: " +
          stringparams["tabulatedkernel"];
      ExceptionHandler::getIstance().raise(message);
    }
  }
  // --------------------------------------------------------------------------
  else {
    string message = "Unrecognised parameter : sph = " 
      + simparams.stringparams["sph"];
    ExceptionHandler::getIstance().raise(message);
  }


  // Create neighbour searching object based on chosen method in params file
  // --------------------------------------------------------------------------
  if (stringparams["neib_search"] == "bruteforce")
    sphneib = new BruteForceSearch(ndim);
  else if (stringparams["neib_search"] == "grid")
    sphneib = new GridSearch(ndim);
  else {
    string message = "Unrecognised parameter : neib_search = " 
      + simparams.stringparams["neib_search"];
    ExceptionHandler::getIstance().raise(message);
  }


  // Create SPH particle integration object
  // --------------------------------------------------------------------------
  if (stringparams["sph_integration"] == "lfkdk")
    sphint = new SphLeapfrogKDK(ndim, vdim, 
				floatparams["accel_mult"],
				floatparams["courant_mult"]);
  else if (stringparams["sph_integration"] == "godunov")
    sphint = new SphGodunovIntegration(ndim, vdim, 
				       floatparams["accel_mult"],
				       floatparams["courant_mult"]);
  else {
    string message = "Unrecognised parameter : sph_integration = " 
      + simparams.stringparams["sph_integration"];
    ExceptionHandler::getIstance().raise(message);
  }


  // Thermal physics object.  If energy equation is chosen, also initiate
  // the energy integration object.
  // --------------------------------------------------------------------------
  if (stringparams["gas_eos"] == "energy_eqn") {
    sph->eos = new Adiabatic(floatparams["temp0"],
			     floatparams["mu_bar"],
			     floatparams["gamma_eos"]);

    if (stringparams["energy_integration"] == "PEC")
      uint = new EnergyPEC(floatparams["energy_mult"]);
    else if (stringparams["energy_integration"] == "godunov")
      uint = new EnergyGodunovIntegration(floatparams["energy_mult"]);
    else {
      string message = "Unrecognised parameter : energy_integration = " 
	+ simparams.stringparams["energy_integration"];
      ExceptionHandler::getIstance().raise(message);
    }

  }
  else if (stringparams["gas_eos"] == "isothermal") 
    sph->eos = new Isothermal(floatparams["temp0"],
			      floatparams["mu_bar"],
			      floatparams["gamma_eos"]);
  else {
    string message = "Unrecognised parameter : gas_eos = " 
      + simparams.stringparams["gas_eos"];
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
  }


  // Set all other parameter variables
  // --------------------------------------------------------------------------
  sph->Nsph = intparams["Npart"];
  //sph->h_fac = floatparams["h_fac"];
  //sph->h_converge = floatparams["h_converge"];
  //sph->hydro_forces = intparams["hydro_forces"];
  //sph->self_gravity = intparams["self_gravity"];
  //sph->avisc = stringparams["avisc"];
  //sph->acond = stringparams["acond"];
  //sph->alpha_visc = floatparams["alpha_visc"];
  //sph->beta_visc = floatparams["beta_visc"];
  //sph->gas_eos = stringparams["gas_eos"];
  Nstepsmax = intparams["Nstepsmax"];
  run_id = stringparams["run_id"];
  out_file_form = stringparams["out_file_form"];
  tend = floatparams["tend"];
  dt_snap = floatparams["dt_snap"];
  noutputstep = intparams["noutputstep"];
  Nlevels = intparams["Nlevels"];
  sph_single_timestep = intparams["sph_single_timestep"];
  nbody_single_timestep = intparams["nbody_single_timestep"];
  sph->riemann_solver = stringparams["riemann_solver"];
  sph->slope_limiter = stringparams["slope_limiter"];
  sph->riemann_order = intparams["riemann_order"];

  return;
}



// ============================================================================
// SphSimulation::Setup
// Main function for setting up a new SPH simulation.
// ============================================================================
void SphSimulation::Setup(void)
{
  debug1("[SphSimulation::Setup]");

  // Read parameters files assigning any contained variables
  simparams.ReadParamsFile(paramfile);

  // Now set up the simulation based on chosen parameters
  SetupSimulation();

}

void SphSimulation::PreSetupForPython(void) {
  debug1("[SphSimulation::PreSetupForPython]");

  // Read the parameters
  simparams.ReadParamsFile(paramfile);

  ProcessParameters();

  sph->AllocateMemory(sph->Nsph);

}

void SphSimulation::ImportArray(double* input, int size, string quantity) {

  //First checks that the size is correct
  if (size != sph->Nsph) {
    stringstream message;
    message << "Error: the array you are passing has a size of " << size << ", but memory has been allocated for " << sph->Nsph << " particles";
    ExceptionHandler::getIstance().raise(message.str());
  }

  //Now sets the pointer to the correct value inside the particle data structure
  FLOAT SphParticle::*quantityp;
  FLOAT (SphParticle::*quantitypvec)[ndimmax];
  bool scalar;
  int index;
  if (quantity=="x") {
    quantitypvec = &SphParticle::r;
    index=0;
    scalar=false;
  }
  else if (quantity=="y") {
    if (ndim<2) {
      string message = "Error: you tried to load a y array, but you are running a 1-d simulation";
      ExceptionHandler::getIstance().raise(message);
    }
    quantitypvec = &SphParticle::r;
    index=1;
    scalar=false;
  }
  else if (quantity == "z") {
    if (ndim<3) {
      string message = "Error: you tried to load a z array, but you are running a simulation with ndim<3";
      ExceptionHandler::getIstance().raise(message);
    }
    quantitypvec = &SphParticle::r;
    index=2;
    scalar=false;
  }
  else if (quantity == "vx") {
    quantitypvec = &SphParticle::v;
    index=0;
    scalar=false;
  }
  else if (quantity == "vy") {
    if (ndim<2) {
      string message = "Error: you tried to load a vy array, but you are running a 1-d simulation";
      ExceptionHandler::getIstance().raise(message);
    }
    quantitypvec = &SphParticle::v;
    index=1;
    scalar=false;
  }
  else if (quantity == "vz") {
    if (ndim<3) {
      string message = "Error: you tried to load a vz array, but you are running a simulation with ndim<3";
      ExceptionHandler::getIstance().raise(message);
    }
    quantitypvec = &SphParticle::v;
    index=2;
    scalar=false;
  }
  else if (quantity=="rho") {
    //TODO: at the moment, if rho or h are uploaded, they will be just ignored.
    //Add some facility to use them
    quantityp = &SphParticle::rho;
    scalar=true;
  }
  else if (quantity == "h") {
    quantityp = &SphParticle::h;
    scalar=true;
  }
  else if (quantity == "u") {
    //TODO: add some facility for uploading either u, T, or cs, and compute automatically the other ones
    //depending on the EOS
    quantityp = &SphParticle::u;
    scalar=true;
  }
  else if (quantity=="m") {
    quantityp = &SphParticle::m;
    scalar=true;
  }
  else {
    string message = "Quantity " + quantity + "not recognised";
    ExceptionHandler::getIstance().raise(message);
  }

  //Finally loops over particles and set the values
  //Note that the syntax for scalar is different from the one for vectors
  if (scalar) {
    int i=0;
    for (SphParticle* particlep = sph->sphdata; particlep < sph->sphdata+size; particlep++, i++) {
      particlep->*quantityp = input[i];
    }
  }
  else {
    int i=0;
    for (SphParticle* particlep = sph->sphdata; particlep < sph->sphdata+size; particlep++, i++) {
      (particlep->*quantitypvec)[index] = input[i];
    }
  }
  return;

}

void SphSimulation::PostSetupForPython(void) {
  debug1("[SphSimulation::PostSetupForPython]");

  PostGeneration();
}

// ============================================================================
// SphSimulation::SetupSimulation
// Main function for setting up a new SPH simulation.
// ============================================================================
void SphSimulation::SetupSimulation(void)
{
  debug1("[SphSimulation::Setup]");

  // Process the parameters file setting up all simulation objects
  ProcessParameters();

  // Generate initial conditions for simulation
  GenerateIC();

  // Call a messy function that does all the rest of the initialisation
  PostGeneration();

  return;
}


//TODO: make this mess more modular (note: initial h computation
//should be done inside the neighbour search)
void SphSimulation::PostGeneration(void) {

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

    sphneib->neibcheck = false;
    sphneib->UpdateAllSphProperties(sph);

    // Search ghost particles
    SearchGhostParticles();

    // Update neighbour tre
    sphneib->UpdateTree(sph,simparams);
  }

  // Compute all SPH particle properties (if SPH particles exist)
  // --------------------------------------------------------------------------
  if (sph->Nsph > 0) {

    cout << "Ntot : " << sph->Ntot << endl;
    level_step = 1;

    // Zero accelerations (perhaps here)
    for (int i=0; i<sph->Ntot; i++) {
      for (int k=0; k<ndim; k++) sph->sphdata[i].a[k] = (FLOAT) 0.0;
      for (int k=0; k<ndim; k++) sph->sphdata[i].agrav[k] = (FLOAT) 0.0;
      sph->sphdata[i].gpot = (FLOAT) 0.0;
      sph->sphdata[i].dudt = (FLOAT) 0.0;
      sph->sphdata[i].active = true;
      sph->sphdata[i].level = level_step;
    }

    // Calculate all SPH properties
    sphneib->neibcheck = true;
    sphneib->UpdateAllSphProperties(sph);

    // Search ghost particles
    SearchGhostParticles();

    // Update neighbour tre
    sphneib->UpdateTree(sph,simparams);

    if (simparams.stringparams["sph"] == "godunov") {
      sphneib->UpdateAllSphDerivatives(sph);
      for (int i=0; i<sph->Ntot; i++) sph->sphdata[i].dt = sph->sphdata[i].h/sph->sphdata[i].sound;
    }

    CopySphDataToGhosts();
    sphneib->UpdateAllSphForces(sph);

    // Add contributions to ghost particles from original neighbours
    //CopyAccelerationFromGhosts();

    // Add accelerations
    for (int i=0; i<sph->Nsph; i++) {
      sph->sphdata[i].active = false;
      for (int k=0; k<ndim; k++)
    sph->sphdata[i].a[k] += sph->sphdata[i].agrav[k];
    }
  }

    // Set r0,v0,a0 for initial step
    sphint->EndTimestep(n,level_step,sph->Nsph,sph->sphdata);
    if (simparams.stringparams["gas_eos"] == "energy_eqn")
      uint->EndTimestep(n,level_step,sph->Nsph,sph->sphdata);

    CalculateDiagnostics();
    diag0 = diag;

    setup = true;

}

// ============================================================================
// SphSimulation::MainLoop
// ============================================================================
void SphSimulation::MainLoop(void)
{
  int i;

  debug2("[SphSimulation::MainLoop]");

  // Compute timesteps for all particles
  if (Nlevels == 1) 
    ComputeGlobalTimestep();
  else 
    ComputeBlockTimesteps();

  // For Godunov SPH, compute compressional heating rates after the timestep 
  // for each particle is known
  if (simparams.stringparams["sph"] == "godunov") {
    for (i=0; i<sph->Ntot; i++)
      sph->sphdata[i].dudt = (FLOAT) 0.0;
    //if (sph->sphdata[i].active) sph->sphdata[i].dudt = (FLOAT) 0.0;
    sphneib->UpdateAllSphDudt(sph);
  }

  // Advance time variables
  n = n + 1;
  Nsteps = Nsteps + 1;
  t = t + timestep;

  // Advance SPH particles positions and velocities
  sphint->AdvanceParticles(n,level_step,sph->Nsph,
			   sph->sphdata,(FLOAT) timestep);
  if (simparams.stringparams["gas_eos"] == "energy_eqn")
    uint->EnergyIntegration(n,level_step,sph->Nsph,
			    sph->sphdata,(FLOAT) timestep);

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

    // Zero accelerations (perhaps)
    for (i=0; i<sph->Ntot; i++) {
      if (sph->sphdata[i].active) {
	for (int k=0; k<ndim; k++) sph->sphdata[i].a[k] = (FLOAT) 0.0;
	for (int k=0; k<ndim; k++) sph->sphdata[i].agrav[k] = (FLOAT) 0.0;
	sph->sphdata[i].gpot = (FLOAT) 0.0;
	sph->sphdata[i].dudt = (FLOAT) 0.0;
      }
    }

    // Calculate all SPH properties
    sphneib->UpdateAllSphProperties(sph);

    if (simparams.stringparams["sph"] == "godunov")
      sphneib->UpdateAllSphDerivatives(sph);

    // Copy properties from original particles to ghost particles
    CopySphDataToGhosts();

    // Calculate all SPH forces
    sphneib->UpdateAllSphForces(sph);

    // Add contributions to ghost particles from original neighbours
    //CopyAccelerationFromGhosts();

    // Add accelerations
    for (i=0; i<sph->Nsph; i++) {
      for (int k=0; k<ndim; k++) 
        sph->sphdata[i].a[k] += sph->sphdata[i].agrav[k];
    }
  }

  // Apply correction steps for both particle and energy integration
  sphint->CorrectionTerms(n,level_step,sph->Nsph,
  			  sph->sphdata,(FLOAT) timestep);
  if (simparams.stringparams["gas_eos"] == "energy_eqn")
    uint->EnergyCorrectionTerms(n,level_step,sph->Nsph,
  				sph->sphdata,(FLOAT) timestep);

  // Set all end-of-step variables
  sphint->EndTimestep(n,level_step,sph->Nsph,sph->sphdata);
  if (simparams.stringparams["gas_eos"] == "energy_eqn")
    uint->EndTimestep(n,level_step,sph->Nsph,sph->sphdata);

  return;
}
