//=============================================================================
//  SphSimulation.cpp
//  Contains all main functions controlling SPH simulation work-flow.
//=============================================================================


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
#include "Nbody.h"
#include "Sph.h"
#include "RiemannSolver.h"
#include "SphSimulationIO.cpp"
#include "SphSimulationIC.cpp"
#include "SphAnalysis.cpp"
#include "SimGhostParticles.cpp"
#include "SphSimulationTimesteps.cpp"
using namespace std;


// Create template class instances of the main SphSimulation object for 
// each dimension used (1, 2 and 3)
template class SphSimulation<1>;
template class SphSimulation<2>;
template class SphSimulation<3>;

template class GodunovSimulation<1>;
template class GodunovSimulation<2>;
template class GodunovSimulation<3>;



//=============================================================================
//  SphSimulationBase
/// Creates a simulation object depending on the dimensionality.
//=============================================================================
SimulationBase* SimulationBase::SimulationFactory
(int ndim,                          ///< [in] No. of dimensions
 Parameters* params)                ///< [in] Pointer to parameters object
{
  string SimulationType;

  //Check ndim
  if (ndim<1 || ndim>3) {
    stringstream msg;
    msg << "Error: ndim must be either 1,2,3; the value " << ndim << "is not allowed!";
    ExceptionHandler::getIstance().raise(msg.str());
  }

  //Set ndim inside the parameters
  params->intparams["ndim"]=ndim;

  //Get the simulation type from the parameters
  //TODO: should the simulation type be passes as a parameter?
  SimulationType = params->stringparams["sph"];

  //Check simulation type
  if (SimulationType != "gradh" && SimulationType != "sm2012" && SimulationType != "godunov" ) {
    string msg = "Error: the simulation type " + SimulationType + " was not recognized";
    ExceptionHandler::getIstance().raise(msg);
  }


  if (ndim==1) {
    if (SimulationType=="gradh" || SimulationType=="sm2012")
      return new SphSimulation<1>(params);
    else if (SimulationType=="godunov")
      return new GodunovSimulation<1>(params);
  }
  else if (ndim==2) {
    if (SimulationType=="gradh" || SimulationType=="sm2012")
      return new SphSimulation<2>(params);
    else if (SimulationType=="godunov")
      return new GodunovSimulation<2>(params);
  }
  else if (ndim==3) {
    if (SimulationType=="gradh" || SimulationType=="sm2012")
      return new SphSimulation<3>(params);
    else if (SimulationType=="godunov")
      return new GodunovSimulation<3>(params);
  }
  return NULL;
}



//=============================================================================
//  SphSimulation::SphSimulation
/// SphSimulation constructor, initialising important simulation variables. 
// ============================================================================
SimulationBase::SimulationBase
(Parameters* params                 ///< [in] Pointer to parameters object
 )
{
  simparams = new Parameters(*params);
  paramfile = "";
  n = 0;
  nresync = 0;
  integration_step = 1;
  Nsteps = 0;
  t = 0.0;
  setup = false;
  ParametersProcessed = false;
}



//=============================================================================
//  SphSimulation::~SphSimulation
/// SphSimulation destructor
//=============================================================================
SimulationBase::~SimulationBase()
{
}



//=============================================================================
//  SphSimulationBase::SetParam
/// Accessor function for modifying a string value. Also checks that the
/// non return point has not been reached
//=============================================================================
void SimulationBase::SetParam(string key, string value) {

  //Error checking
  if (ParametersProcessed) {
    string msg = "Error: the non-return point for setting parameters has been reached!";
    ExceptionHandler::getIstance().raise(msg);
  }
  if (key=="ndim") {
    string msg = "Error: it's not possible to change the number of dimensions!";
    ExceptionHandler::getIstance().raise(msg);
  }

  simparams->SetParameter (key, value);
}



//=============================================================================
//  SphSimulationBase::SetParam
/// Accessor function for modifying an int value, wrapper around the one for string value.
/// Also checks that the non return point has not been reached
//=============================================================================
void SimulationBase::SetParam(string key, int value) {
  ostringstream convert;
  convert << value;
  SetParam (key, convert.str());
}



//=============================================================================
//  SphSimulationBase::SetParam
/// Accessor function for modifying a float value, wrapper around the one for string value.
/// Also checks that the non return point has not been reached
//=============================================================================
void SimulationBase::SetParam(string key, float value) {
  ostringstream convert;
  convert << value;
  SetParam (key, convert.str());
}



//=============================================================================
//  SphSimulationBase::GetParam
/// Accessor function for getting a parameter value
/// Wrapper around the corresponding function in Parameters
//=============================================================================
string SimulationBase::GetParam(string key) {

  return simparams->GetParameter(key);

}



//=============================================================================
//  SphSimulation::Run
/// Controls the simulation main loop, including exit conditions.  If provided
/// (optional argument), will only advance the simulation by 'Nadvance' steps.
//=============================================================================
void SimulationBase::Run
(int Nadvance                       ///< [in] Selected max no. of integer 
 )                                  ///< timesteps (Optional argument).
{
  int Ntarget;                      // Target step no before finishing 
                                    // main code integration.

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
  OutputDiagnostics();
  UpdateDiagnostics();

  return;
}


template <int ndim>
void Simulation<ndim>::UpdateDiagnostics () {
  diag.Eerror = fabs(diag0.Etot - diag.Etot)/fabs(diag0.Etot);
  cout << "Eerror : " << diag.Eerror << endl;
}


//=============================================================================
//  SphSimulation::InteractiveRun
/// Controls the simulation main loop, including exit conditions.
/// If provided, will only advance the simulation by 'Nadvance' steps.
//=============================================================================
void SimulationBase::InteractiveRun
(int Nadvance                           ///< [in] Selected max no. of integer 
 )                                      ///< timesteps (Optional argument).
{
  int Ntarget;                          // Selected integer timestep
  DOUBLE tdiff = 0.0;                   // Measured time difference
  DOUBLE tpython = 8.0;                 // Python viewer update time
  clock_t tstart = clock();             // Initial CPU clock time

  debug2("[SphSimulation::InteractiveRun]");

  // Set integer timestep exit condition if provided as parameter.
  if (Nadvance < 0) Ntarget = Nstepsmax;
  else Ntarget = Nsteps + Nadvance;

  // Continue to run simulation until we reach the required time, or
  // exeeded the maximum allowed number of steps.
  // --------------------------------------------------------------------------
  while (t < tend && Nsteps < Ntarget && tdiff < tpython) {

    MainLoop();
    Output();

    // Measure CPU clock time difference since current function was called
    tdiff = (DOUBLE) (clock() - tstart) / (DOUBLE) CLOCKS_PER_SEC;

  }
  // --------------------------------------------------------------------------

  CalculateDiagnostics();
  UpdateDiagnostics();

  return;
}



//=============================================================================
//  SphSimulation::Output
/// Controls when regular output snapshots are written by the code.
//=============================================================================
void SimulationBase::Output(void)
{
  string filename;                  // Output snapshot filename
  string nostring;                  // ???
  stringstream ss;                  // ???

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



//============================================================================
//  SphSimulation::GenerateIC
/// Generate initial conditions for SPH simulation chosen in parameters file.
//=============================================================================
template <int ndim>
void Simulation<ndim>::GenerateIC(void)
{
  debug2("[SphSimulation::GenerateIC]");

  // Generate initial conditions
  if (simparams->stringparams["ic"] == "file")
    ReadSnapshotFile(simparams->stringparams["in_file"],
		     simparams->stringparams["in_file_form"]);
  else if (simparams->stringparams["ic"] == "box")
    UniformBox();
  else if (simparams->stringparams["ic"] == "sphere")
    UniformSphere();
  else if (simparams->stringparams["ic"] == "cdiscontinuity")
    ContactDiscontinuity();
  else if (simparams->stringparams["ic"] == "noh")
    NohProblem();
  else if (simparams->stringparams["ic"] == "sedov")
    SedovBlastWave();
  else if (simparams->stringparams["ic"] == "shearflow")
    ShearFlow();
  else if (simparams->stringparams["ic"] == "shocktube")
    ShockTube();
  else if (simparams->stringparams["ic"] == "soundwave")
    SoundWave();
  else if (simparams->stringparams["ic"] == "khi")
    KHI();
  else if (simparams->stringparams["ic"] == "binary")
    BinaryStar();
  else if (simparams->stringparams["ic"] == "python")
    return;
  else {
    string message = "Unrecognised parameter : ic = " 
      + simparams->stringparams["ic"];
    ExceptionHandler::getIstance().raise(message);
  }

  // Check that the initial conditions are valid
  CheckInitialConditions();

  return;
}



//=============================================================================
//  SphSimulation::ProcessParameters
/// Process all the options chosen in the parameters file, setting various 
/// simulation variables and creating important simulation objects.
//=============================================================================
template <int ndim>
void Simulation<ndim>::ProcessParameters(void)
{
  aviscenum avisc;                  // Artificial viscosity enum
  acondenum acond;                  // Artificial conductivity enum

  // Local references to parameter variables for brevity
  map<string, int> &intparams = simparams->intparams;
  map<string, float> &floatparams = simparams->floatparams;
  map<string, string> &stringparams = simparams->stringparams;

  debug2("[SphSimulation::ProcessParameters]");

  // Sanity check for valid dimensionality
  if (ndim < 1 || ndim > 3) {
    string message = "Invalid dimensionality chosen : ndim = " + ndim;
    ExceptionHandler::getIstance().raise(message);
  }

  // Set the enum for artificial viscosity
  if (stringparams["avisc"] == "none")
    avisc = noneav;
  else if (stringparams["avisc"] == "mon97")
    avisc = mon97;
  else {
    string message = "Unrecognised parameter : avisc = " +
      simparams->stringparams["avisc"];
    ExceptionHandler::getIstance().raise(message);
  }

  // Set the enum for artificial conductivity
  if (stringparams["acond"] == "none")
    acond = noneac;
  else if (stringparams["acond"] == "wadsley2008")
    acond = wadsley2008;
  else if (stringparams["acond"] == "price2008")
    acond = price2008;
  else {
    string message = "Unrecognised parameter : acond = " +
        simparams->stringparams["acond"];
    ExceptionHandler::getIstance().raise(message);
  }

  // Create SPH object based on chosen method in params file
  // --------------------------------------------------------------------------
  if (stringparams["sph"] == "gradh") {
    string KernelName = stringparams["kernel"];
    if (intparams["tabulated_kernel"] == 1) {
      sph = new GradhSph<ndim, TabulatedKernel> 
	(intparams["hydro_forces"], intparams["self_gravity"],
	 floatparams["alpha_visc"], floatparams["beta_visc"],
	 floatparams["h_fac"], floatparams["h_converge"], 
	 avisc, acond, stringparams["gas_eos"], KernelName);
    }
    else if (intparams["tabulated_kernel"] == 0) {
      // Depending on the kernel, instantiate a different GradSph object
      if (KernelName == "m4") {
	sph = new GradhSph<ndim, M4Kernel> 
	  (intparams["hydro_forces"], intparams["self_gravity"],
	   floatparams["alpha_visc"], floatparams["beta_visc"],
	   floatparams["h_fac"], floatparams["h_converge"],
	   avisc, acond, stringparams["gas_eos"], KernelName);
      }
      else if (KernelName == "quintic") {
	sph = new GradhSph<ndim, QuinticKernel> 
	  (intparams["hydro_forces"], intparams["self_gravity"],
	   floatparams["alpha_visc"], floatparams["beta_visc"],
	   floatparams["h_fac"], floatparams["h_converge"],
	   avisc, acond, stringparams["gas_eos"], KernelName);
      }
      else if (KernelName == "gaussian") {
	sph = new GradhSph<ndim, GaussianKernel> 
	  (intparams["hydro_forces"], intparams["self_gravity"],
	   floatparams["alpha_visc"], floatparams["beta_visc"],
	   floatparams["h_fac"], floatparams["h_converge"],
	   avisc, acond, stringparams["gas_eos"], KernelName);
      }
      else {
	string message = "Unrecognised parameter : kernel = " +
	  simparams->stringparams["kernel"];
	ExceptionHandler::getIstance().raise(message);
      }
    }
    else {
      string message = "Invalid option for the tabulated_kernel parameter: " +
	stringparams["tabulated_kernel"];
      ExceptionHandler::getIstance().raise(message);
    }
  }
  // --------------------------------------------------------------------------
  else if (stringparams["sph"] == "sm2012") {
    string KernelName = stringparams["kernel"];
    if (intparams["tabulated_kernel"] == 1) {
      sph = new SM2012Sph<ndim, TabulatedKernel> 
        (intparams["hydro_forces"], intparams["self_gravity"],
	 floatparams["alpha_visc"], floatparams["beta_visc"],
	 floatparams["h_fac"], floatparams["h_converge"],
	 avisc, acond, stringparams["gas_eos"], KernelName);
    }
    else if (intparams["tabulated_kernel"] == 0){
      // Depending on the kernel, instantiate a different GradSph object
      if (KernelName == "m4") {
	sph = new SM2012Sph<ndim, M4Kernel> 
	  (intparams["hydro_forces"], intparams["self_gravity"],
	   floatparams["alpha_visc"], floatparams["beta_visc"],
	   floatparams["h_fac"], floatparams["h_converge"],
	   avisc, acond, stringparams["gas_eos"], KernelName);
      }
      else if (KernelName == "quintic") {
	sph = new SM2012Sph<ndim, QuinticKernel> 
	  (intparams["hydro_forces"], intparams["self_gravity"],
	   floatparams["alpha_visc"], floatparams["beta_visc"],
	   floatparams["h_fac"], floatparams["h_converge"],
	   avisc, acond, stringparams["gas_eos"], KernelName);
      }
      else if (KernelName == "gaussian") {
	sph = new SM2012Sph<ndim, GaussianKernel> 
	  (intparams["hydro_forces"], intparams["self_gravity"],
	   floatparams["alpha_visc"], floatparams["beta_visc"],
	   floatparams["h_fac"], floatparams["h_converge"],
	   avisc, acond, stringparams["gas_eos"], KernelName);
      }
      else {
	string message = "Unrecognised parameter : kernel = " +
	  simparams->stringparams["kernel"];
	ExceptionHandler::getIstance().raise(message);
      }
    }
    else {
      string message = "Invalid option for the tabulated_kernel parameter: " +
	stringparams["tabulated_kernel"];
      ExceptionHandler::getIstance().raise(message);
    }
  }
  // --------------------------------------------------------------------------
  else if (stringparams["sph"] == "godunov") {
    string KernelName = stringparams["kernel"];
    if (intparams["tabulated_kernel"] == 1) {
      sph = new GodunovSph<ndim, TabulatedKernel> 
	(intparams["hydro_forces"], intparams["self_gravity"],
	 floatparams["alpha_visc"], floatparams["beta_visc"],
	 floatparams["h_fac"], floatparams["h_converge"],
	 avisc, acond, stringparams["gas_eos"], KernelName);
    }
    else if (intparams["tabulated_kernel"] == 0){
      // Depending on the kernel, instantiate a different GradSph object
      if (KernelName == "m4") {
	sph = new GodunovSph<ndim, M4Kernel> 
	  (intparams["hydro_forces"], intparams["self_gravity"],
	   floatparams["alpha_visc"], floatparams["beta_visc"],
	   floatparams["h_fac"], floatparams["h_converge"],
	   avisc, acond, stringparams["gas_eos"], KernelName);
      }
      else if (KernelName == "quintic") {
	sph = new GodunovSph<ndim, QuinticKernel> 
	  (intparams["hydro_forces"], intparams["self_gravity"],
	   floatparams["alpha_visc"], floatparams["beta_visc"],
	   floatparams["h_fac"], floatparams["h_converge"],
	   avisc, acond, stringparams["gas_eos"], KernelName);
      }
      else if (KernelName == "gaussian") {
	sph = new GodunovSph<ndim, GaussianKernel> 
	  (intparams["hydro_forces"], intparams["self_gravity"],
	   floatparams["alpha_visc"], floatparams["beta_visc"],
	   floatparams["h_fac"], floatparams["h_converge"],
	   avisc, acond, stringparams["gas_eos"], KernelName);
      }
      else {
	string message = "Unrecognised parameter : kernel = " +
	  simparams->stringparams["kernel"];
	ExceptionHandler::getIstance().raise(message);
      }
    }
    else {
      string message = "Invalid option for the tabulated_kernel parameter: " +
	stringparams["tabulated_kernel"];
      ExceptionHandler::getIstance().raise(message);
    }
  }
  // --------------------------------------------------------------------------
  else {
    string message = "Unrecognised parameter : sph = " 
      + simparams->stringparams["sph"];
    ExceptionHandler::getIstance().raise(message);
  }
  // --------------------------------------------------------------------------


  // Thermal physics object.  If energy equation is chosen, also initiate
  // the energy integration object.
  // --------------------------------------------------------------------------
  string gas_eos = stringparams["gas_eos"];
  if (gas_eos == "energy_eqn") {
    sph->eos = new Adiabatic<ndim>(floatparams["temp0"],
				   floatparams["mu_bar"],
				   floatparams["gamma_eos"]);
  }
  else if (gas_eos == "isothermal")
    sph->eos = new Isothermal<ndim>(floatparams["temp0"],
				    floatparams["mu_bar"],
				    floatparams["gamma_eos"]);
  else if (gas_eos == "barotropic")
    sph->eos = new Barotropic<ndim>(floatparams["temp0"],
				    floatparams["mu_bar"],
				    floatparams["gamma_eos"],
				    floatparams["rho_bary"]);
  else {
    string message = "Unrecognised parameter : gas_eos = " + gas_eos;
    ExceptionHandler::getIstance().raise(message);
  }


  // Riemann solver object
  // --------------------------------------------------------------------------
  string riemann = stringparams["riemann_solver"];
  if (riemann == "exact") {
    sph->riemann = new ExactRiemannSolver(floatparams["gamma_eos"]);
  }
  else if (riemann == "hllc") {
    sph->riemann = new HllcRiemannSolver(floatparams["gamma_eos"]);
  }
  else {
    string message = "Unrecognised parameter : riemann_solver = "
      + riemann;
    ExceptionHandler::getIstance().raise(message);
  }


  // Create neighbour searching object based on chosen method in params file
  // --------------------------------------------------------------------------
  if (stringparams["neib_search"] == "bruteforce")
    sphneib = new BruteForceSearch<ndim>;
  else if (stringparams["neib_search"] == "grid")
    sphneib = new GridSearch<ndim>;
  else {
    string message = "Unrecognised parameter : neib_search = " 
      + simparams->stringparams["neib_search"];
    ExceptionHandler::getIstance().raise(message);
  }


  // Create SPH particle integration object
  // --------------------------------------------------------------------------
  if (stringparams["sph_integration"] == "lfkdk") {
    sphint = new SphLeapfrogKDK<ndim>(floatparams["accel_mult"],
			              floatparams["courant_mult"]);}
  else if (stringparams["sph_integration"] == "godunov")
    sphint = new SphGodunovIntegration<ndim>(floatparams["accel_mult"],
				             floatparams["courant_mult"]);
  else {
    string message = "Unrecognised parameter : sph_integration = " 
      + simparams->stringparams["sph_integration"];
    ExceptionHandler::getIstance().raise(message);
  }


  // Energy integration object
  // --------------------------------------------------------------------------
  if (stringparams["energy_integration"] == "PEC")
    uint = new EnergyPEC<ndim>(floatparams["energy_mult"]);
  else if (stringparams["energy_integration"] == "godunov")
    uint = new EnergyGodunovIntegration<ndim>(floatparams["energy_mult"]);
  else {
    string message = "Unrecognised parameter : energy_integration = "
      + simparams->stringparams["energy_integration"];
    ExceptionHandler::getIstance().raise(message);
  }

  
  // Create N-body object based on chosen method in params file
  // --------------------------------------------------------------------------
  if (stringparams["nbody"] == "lfkdk") {
    string KernelName = stringparams["kernel"];
    if (intparams["tabulated_kernel"] == 1) {
      nbody = new NbodyLeapfrogKDK<ndim, TabulatedKernel> 
	(intparams["nbody_softening"], intparams["sub_systems"],
	 floatparams["nbody_mult"], KernelName);
    }
    else if (intparams["tabulated_kernel"] == 0) {
      // Depending on the kernel, instantiate a different GradSph object
      if (KernelName == "m4") {
	nbody = new NbodyLeapfrogKDK<ndim, M4Kernel> 
	  (intparams["nbody_softening"], intparams["sub_systems"],
	   floatparams["nbody_mult"], KernelName);
      }
      else if (KernelName == "quintic") {
	nbody = new NbodyLeapfrogKDK<ndim, QuinticKernel> 
	  (intparams["nbody_softening"], intparams["sub_systems"],
	   floatparams["nbody_mult"], KernelName);
      }
      else if (KernelName == "gaussian") {
	nbody = new NbodyLeapfrogKDK<ndim, GaussianKernel> 
	  (intparams["nbody_softening"], intparams["sub_systems"],
	   floatparams["nbody_mult"], KernelName);
      }
      else {
	string message = "Unrecognised parameter : kernel = " +
	  simparams->stringparams["kernel"];
	ExceptionHandler::getIstance().raise(message);
      }
    }
    else {
      string message = "Invalid option for the tabulated_kernel parameter: " +
	stringparams["tabulated_kernel"];
      ExceptionHandler::getIstance().raise(message);
    }
  }
  // --------------------------------------------------------------------------
  else if (stringparams["nbody"] == "hermite4") {
    string KernelName = stringparams["kernel"];
    if (intparams["tabulated_kernel"] == 1) {
      nbody = new NbodyHermite4<ndim, TabulatedKernel> 
	(intparams["nbody_softening"], intparams["sub_systems"],
	 floatparams["nbody_mult"], KernelName);
    }
    else if (intparams["tabulated_kernel"] == 0) {
      // Depending on the kernel, instantiate a different GradSph object
      if (KernelName == "m4") {
	nbody = new NbodyHermite4<ndim, M4Kernel> 
	  (intparams["nbody_softening"], intparams["sub_systems"],
	   floatparams["nbody_mult"], KernelName);
      }
      else if (KernelName == "quintic") {
	nbody = new NbodyHermite4<ndim, QuinticKernel> 
	  (intparams["nbody_softening"], intparams["sub_systems"],
	   floatparams["nbody_mult"], KernelName);
      }
      else if (KernelName == "gaussian") {
	nbody = new NbodyHermite4<ndim, GaussianKernel> 
	  (intparams["nbody_softening"], intparams["sub_systems"],
	   floatparams["nbody_mult"], KernelName);
      }
      else {
	string message = "Unrecognised parameter : kernel = " +
	  simparams->stringparams["kernel"];
	ExceptionHandler::getIstance().raise(message);
      }
    }
    else {
      string message = "Invalid option for the tabulated_kernel parameter: " +
	stringparams["tabulated_kernel"];
      ExceptionHandler::getIstance().raise(message);
    }
  }
  // --------------------------------------------------------------------------
  else {
    string message = "Unrecognised parameter : nbody = " 
      + simparams->stringparams["nbody"];
    ExceptionHandler::getIstance().raise(message);
  }
  // --------------------------------------------------------------------------


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

  // Flag that we've processed all parameters already
  ParametersProcessed = true;

  return;
}


//=============================================================================
//  SphSimulation::PreSetupForPython
/// Initialisation routine called by python interface.
//=============================================================================
template <int ndim>
void Simulation<ndim>::PreSetupForPython(void)
{
  debug1("[SphSimulation::PreSetupForPython]");

  //Check that IC type is really python
  if (simparams->stringparams["ic"] != "python") {
    string msg = "Error: you should call this function only if you are using \"python\" as \"ic\" parameter";
    ExceptionHandler::getIstance().raise(msg);
  }

  if (ParametersProcessed) {
    string msg = "Error: this function has been already called!";
    ExceptionHandler::getIstance().raise(msg);
  }

  ProcessParameters();

  sph->AllocateMemory(sph->Nsph);

  return;
}



//=============================================================================
//  SphSimulation::ImportArray
/// Import an array containing particle properties from python to C++ arrays.
//=============================================================================
template <int ndim>
void Simulation<ndim>::ImportArray
(double* input,                         ///< ..
 int size,                              ///< ..
 string quantity)                       ///< ..
{
  FLOAT SphParticle<ndim>::*quantityp;
  FLOAT (SphParticle<ndim>::*quantitypvec)[ndim];
  bool scalar;
  int index;

  debug2("[SphSimulation::ImportArray]");

  //Check that PreSetup has been called
  if (! ParametersProcessed) {
    string msg = "Error: before calling ImportArray, you need to call PreSetupForPython!";
    ExceptionHandler::getIstance().raise(msg);
  }

  //Check that the size is correct
  if (size != sph->Nsph) {
    stringstream message;
    message << "Error: the array you are passing has a size of " 
	    << size << ", but memory has been allocated for " 
	    << sph->Nsph << " particles";
    ExceptionHandler::getIstance().raise(message.str());
  }

  // Now set pointer to the correct value inside the particle data structure
  // --------------------------------------------------------------------------
  if (quantity == "x") {
    quantitypvec = &SphParticle<ndim>::r;
    index = 0;
    scalar = false;
  }
  // --------------------------------------------------------------------------
  else if (quantity == "y") {
    if (ndim < 2) {
      string message = "Error: loading y-coordinate array for ndim < 2";
      ExceptionHandler::getIstance().raise(message);
    }
    quantitypvec = &SphParticle<ndim>::r;
    index = 1;
    scalar = false;
  }
  // --------------------------------------------------------------------------
  else if (quantity == "z") {
    if (ndim < 3) {
      string message = "Error: loading y-coordinate array for ndim < 3";
      ExceptionHandler::getIstance().raise(message);
    }
    quantitypvec = &SphParticle<ndim>::r;
    index = 2;
    scalar = false;
  }
  // --------------------------------------------------------------------------
  else if (quantity == "vx") {
    quantitypvec = &SphParticle<ndim>::v;
    index = 0;
    scalar = false;
  }
  // --------------------------------------------------------------------------
  else if (quantity == "vy") {
    if (ndim < 2) {
      string message = "Error: loading vy array for ndim < 2";
      ExceptionHandler::getIstance().raise(message);
    }
    quantitypvec = &SphParticle<ndim>::v;
    index = 1;
    scalar = false;
  }
  // --------------------------------------------------------------------------
  else if (quantity == "vz") {
    if (ndim < 3) {
      string message = "Error: loading vz array for ndim < 3";
      ExceptionHandler::getIstance().raise(message);
    }
    quantitypvec = &SphParticle<ndim>::v;
    index = 2;
    scalar = false;
  }
  // --------------------------------------------------------------------------
  else if (quantity=="rho") {
    //TODO: at the moment, if rho or h are uploaded, they will be just ignored.
    //Add some facility to use them
    quantityp = &SphParticle<ndim>::rho;
    scalar = true;
  }
  // --------------------------------------------------------------------------
  else if (quantity == "h") {
    quantityp = &SphParticle<ndim>::h;
    scalar = true;
  }
  // --------------------------------------------------------------------------
  else if (quantity == "u") {
    //TODO: add some facility for uploading either u, T, or cs, and compute automatically the other ones
    //depending on the EOS
    quantityp = &SphParticle<ndim>::u;
    scalar=true;
  }
  // --------------------------------------------------------------------------
  else if (quantity=="m") {
    quantityp = &SphParticle<ndim>::m;
    scalar = true;
  }
  // --------------------------------------------------------------------------
  else {
    string message = "Quantity " + quantity + "not recognised";
    ExceptionHandler::getIstance().raise(message);
  }
  // --------------------------------------------------------------------------


  // Finally loop over particles and set all values
  // (Note that the syntax for scalar is different from the one for vectors)
  // --------------------------------------------------------------------------
  if (scalar) {
    int i=0;
    for (SphParticle<ndim>* particlep = sph->sphdata; 
	 particlep < sph->sphdata+size; particlep++, i++) {
      particlep->*quantityp = input[i];
    }
  }
  else {
    int i=0;
    for (SphParticle<ndim>* particlep = sph->sphdata; 
	 particlep < sph->sphdata+size; particlep++, i++) {
      (particlep->*quantitypvec)[index] = input[i];
    }
  }

  return;
}


//=============================================================================
//  SphSimulation::SetupSimulation
/// Main function for setting up a new SPH simulation.
//=============================================================================
void SimulationBase::SetupSimulation(void)
{
  debug1("[SphSimulation::Setup]");

  if (setup) {
    string msg = "This simulation has been already set up";
    ExceptionHandler::getIstance().raise(msg);
  }


  // Process the parameters file setting up all simulation objects
  if (simparams->stringparams["ic"]=="python") {
    if (!ParametersProcessed) {
      string msg = "Error: you are attempting to setup a simulation with initial conditions generated"
          "from Python. Before setting up the simulation, you need to import the initial conditions";
      ExceptionHandler::getIstance().raise(msg);
    }
  }
  else {
    if (ParametersProcessed) {
      string msg = "The parameters of the simulation have been already processed."
          "It means that you shouldn't be calling this function, please consult the documentation.";
      ExceptionHandler::getIstance().raise(msg);
    }
    ProcessParameters();
  }

  // Generate initial conditions for simulation
  GenerateIC();

  // Call a messy function that does all the rest of the initialisation
  PostGeneration();

  return;
}



//TODO: make this mess more modular (note: initial h computation
//should be done inside the neighbour search)
//=============================================================================
//  SphSimulation::PostGeneration
/// ..
//=============================================================================
template <int ndim>
void Simulation<ndim>::PostGeneration(void)
{
  int i;                            // Particle counter
  int k;                            // Dimension counter

  debug2("[SphSimulation::PostGeneration]");

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
    sphneib->UpdateTree(sph,*simparams);

    sphneib->neibcheck = false;
    sphneib->UpdateAllSphProperties(sph);

    // Search ghost particles
    SearchGhostParticles();

    // Update neighbour tree
    sphneib->UpdateTree(sph,*simparams);
  }

  // Compute all initial SPH particle properties (if SPH particles exist)
  // --------------------------------------------------------------------------
  if (sph->Nsph > 0) {

    cout << "Ntot : " << sph->Ntot << endl;
    level_step = 1;

    // Zero accelerations (perhaps here)
    for (int i=0; i<sph->Ntot; i++) sph->sphdata[i].active = true;

    // Calculate all SPH properties
    sphneib->neibcheck = true;
    sphneib->UpdateAllSphProperties(sph);

    // Search ghost particles
    SearchGhostParticles();

    // Update neighbour tre
    sphneib->UpdateTree(sph,*simparams);
    sphneib->UpdateAllSphProperties(sph);


    if (simparams->stringparams["sph"] == "godunov") {
      sphneib->UpdateAllSphDerivatives(sph);
      //for (int i=0; i<sph->Ntot; i++) 
      //sph->sphdata[i].dt = sph->sphdata[i].h/sph->sphdata[i].sound;
    }
  }


  // Compute all initial N-body terms
  // --------------------------------------------------------------------------
  if (nbody->Nstar > 0) {

    // Zero all acceleration terms
    for (i=0; i<nbody->Nstar; i++) {

      for (k=0; k<ndim; k++) nbody->stardata[i].a[k] = 0.0;
      for (k=0; k<ndim; k++) nbody->stardata[i].adot[k] = 0.0;
      for (k=0; k<ndim; k++) nbody->stardata[i].a2dot[k] = 0.0;
      for (k=0; k<ndim; k++) nbody->stardata[i].a3dot[k] = 0.0;
      nbody->stardata[i].gpot = 0.0;
      nbody->stardata[i].active = true;
      nbody->stardata[i].level = level_step;
      nbody->stardata[i].nstep = 1;
      nbody->nbodydata[i] = &(nbody->stardata[i]);
    }

    nbody->CalculateDirectGravForces(nbody->Nnbody,nbody->nbodydata);

  }


  // Compute all initial SPH force terms
  // --------------------------------------------------------------------------
  if (sph->Nsph > 0) {

    // Zero accelerations (here for now)
    for (i=0; i<sph->Ntot; i++) {
      for (k=0; k<ndim; k++) sph->sphdata[i].a[k] = (FLOAT) 0.0;
      for (k=0; k<ndim; k++) sph->sphdata[i].agrav[k] = (FLOAT) 0.0;
      sph->sphdata[i].gpot = (FLOAT) 0.0;
      sph->sphdata[i].dudt = (FLOAT) 0.0;
      sph->sphdata[i].active = true;
      sph->sphdata[i].level = level_step;
      sph->sphdata[i].nstep = 1;
    }

    CopySphDataToGhosts();

    // Compute timesteps for all particles
    if (simparams->stringparams["sph"] == "godunov") {
      if (Nlevels == 1) 
	ComputeGlobalTimestep();
      else 
	ComputeBlockTimesteps();
      if (sph->hydro_forces == 1) sphneib->UpdateAllSphForces(sph);
      if (sph->self_gravity == 1) sphneib->UpdateAllSphGravForces(sph);
      sphneib->UpdateAllSphDudt(sph);
    }
    else {
      if (sph->hydro_forces == 1) sphneib->UpdateAllSphForces(sph);
      if (sph->self_gravity == 1) sphneib->UpdateAllSphGravForces(sph);
    }

    // Add accelerations
    for (i=0; i<sph->Nsph; i++) {
      sph->sphdata[i].active = false;
      for (k=0; k<ndim; k++)
	sph->sphdata[i].a[k] += sph->sphdata[i].agrav[k];
    }

    CopySphDataToGhosts();

  }

  // Set particle values for initial step (e.g. r0, v0, a0)
  sphint->EndTimestep(n,level_step,sph->Nsph,sph->sphdata);
  if (simparams->stringparams["gas_eos"] == "energy_eqn")
    uint->EndTimestep(n,level_step,sph->Nsph,sph->sphdata);
  nbody->EndTimestep(n,nbody->Nstar,nbody->nbodydata);
  
  CalculateDiagnostics();
  diag0 = diag;
  
  setup = true;

  return;
}



//=============================================================================
//  SphSimulation::MainLoop
/// Main SPH simulation integration loop.
//=============================================================================
template <int ndim>
void SphSimulation<ndim>::MainLoop(void)
{
  int i;                            // Particle loop counter
  int k;                            // Dimension counter

  debug2("[SphSimulation::MainLoop]");

  // Compute timesteps for all particles
  if (simparams->stringparams["sph"] != "godunov") {
    if (Nlevels == 1)
      this->ComputeGlobalTimestep();
    else 
      this->ComputeBlockTimesteps();
  }

  // For Godunov SPH, compute compressional heating rates after the timestep 
  // for each particle is known
  //if (simparams->stringparams["sph"] == "godunov") {
  //  for (i=0; i<sph->Ntot; i++)
  //    sph->sphdata[i].dudt = (FLOAT) 0.0;
  //  //if (sph->sphdata[i].active) sph->sphdata[i].dudt = (FLOAT) 0.0;
  //  sphneib->UpdateAllSphDudt(sph);
  //}

  // Advance time variables
  n = n + 1;
  Nsteps = Nsteps + 1;
  t = t + timestep;

  // Advance SPH particles positions and velocities
  sphint->AdvanceParticles(n,level_step,sph->Nsph,
			   sph->sphdata,(FLOAT) timestep);
  if (simparams->stringparams["gas_eos"] == "energy_eqn")
    uint->EnergyIntegration(n,level_step,sph->Nsph,
			    sph->sphdata,(FLOAT) timestep);
  nbody->AdvanceParticles(n,nbody->Nstar,nbody->nbodydata,timestep);

  // Check all boundary conditions
  this->CheckBoundaries();

  // --------------------------------------------------------------------------
  if (sph->Nsph > 0) {
    
    // Reorder particles

    // Search ghost particles
    this->SearchGhostParticles();

    // Update neighbour tree
    sphneib->UpdateTree(sph,*simparams);
  }


  // --------------------------------------------------------------------------
  if (sph->Nsph > 0) {

    // Calculate all SPH properties
    sphneib->UpdateAllSphProperties(sph);

  }


  // Compute N-body forces
  // --------------------------------------------------------------------------
  if (nbody->Nstar > 0) {

    // Zero all acceleration terms
    for (i=0; i<nbody->Nnbody; i++) {
      if (nbody->nbodydata[i]->active) {
	for (k=0; k<ndim; k++) nbody->nbodydata[i]->a[k] = 0.0;
	for (k=0; k<ndim; k++) nbody->nbodydata[i]->adot[k] = 0.0;
	for (k=0; k<ndim; k++) nbody->nbodydata[i]->a2dot[k] = 0.0;
	for (k=0; k<ndim; k++) nbody->nbodydata[i]->a3dot[k] = 0.0;
	nbody->nbodydata[i]->gpot = 0.0;
      }
    }

    nbody->CalculateDirectGravForces(nbody->Nstar,nbody->nbodydata);
  }


  // --------------------------------------------------------------------------
  if (sph->Nsph > 0) {

    // Compute timesteps for all particles
    if (simparams->stringparams["sph"] == "godunov") {
      if (Nlevels == 1) 
	this->ComputeGlobalTimestep();
      else 
	this->ComputeBlockTimesteps();
    }

    if (simparams->stringparams["sph"] == "godunov")
      sphneib->UpdateAllSphDerivatives(sph);

    // Copy properties from original particles to ghost particles
    this->CopySphDataToGhosts();

    // Zero accelerations (perhaps)
    for (i=0; i<sph->Ntot; i++) {
      if (sph->sphdata[i].active) {
	for (k=0; k<ndim; k++) sph->sphdata[i].a[k] = (FLOAT) 0.0;
	for (k=0; k<ndim; k++) sph->sphdata[i].agrav[k] = (FLOAT) 0.0;
	sph->sphdata[i].gpot = (FLOAT) 0.0;
	sph->sphdata[i].dudt = (FLOAT) 0.0;
      }
    }

    // Calculate all SPH forces
    if (sph->hydro_forces == 1) sphneib->UpdateAllSphForces(sph);
    if (sph->self_gravity == 1) sphneib->UpdateAllSphGravForces(sph);

    // Compute contribution to grav. accel from stars
    for (i=0; i<sph->Nsph; i++)
      if (sph->sphdata[i].active)
	sph->ComputeStarGravForces(nbody->Nnbody,nbody->nbodydata,
				   sph->sphdata[i]);

    // Add accelerations
    for (i=0; i<sph->Nsph; i++) {
      for (k=0; k<ndim; k++) 
        sph->sphdata[i].a[k] += sph->sphdata[i].agrav[k];
    }

    if (simparams->stringparams["sph"] == "godunov")
      sphneib->UpdateAllSphDudt(sph);

  }

  // Apply correction steps for both particle and energy integration
  sphint->CorrectionTerms(n,level_step,sph->Nsph,
  			  sph->sphdata,(FLOAT) timestep);
  if (simparams->stringparams["gas_eos"] == "energy_eqn")
    uint->EnergyCorrectionTerms(n,level_step,sph->Nsph,
  				sph->sphdata,(FLOAT) timestep);
  nbody->CorrectionTerms(n,nbody->Nnbody,nbody->nbodydata,timestep);

  // Set all end-of-step variables
  sphint->EndTimestep(n,level_step,sph->Nsph,sph->sphdata);
  if (simparams->stringparams["gas_eos"] == "energy_eqn")
    uint->EndTimestep(n,level_step,sph->Nsph,sph->sphdata);
  nbody->EndTimestep(n,nbody->Nnbody,nbody->nbodydata);

  return;
}

//=============================================================================
//  SphSimulation::MainLoop
/// Main SPH simulation integration loop.
//=============================================================================
template <int ndim>
void GodunovSimulation<ndim>::MainLoop(void)
{
  int i;                            // Particle loop counter
  int k;                            // Dimension counter

  debug2("[SphSimulation::MainLoop]");

  // Compute timesteps for all particles
  if (simparams->stringparams["sph"] != "godunov") {
    if (Nlevels == 1)
      this->ComputeGlobalTimestep();
    else
      this->ComputeBlockTimesteps();
  }

  // For Godunov SPH, compute compressional heating rates after the timestep
  // for each particle is known
  //if (simparams->stringparams["sph"] == "godunov") {
  //  for (i=0; i<sph->Ntot; i++)
  //    sph->sphdata[i].dudt = (FLOAT) 0.0;
  //  //if (sph->sphdata[i].active) sph->sphdata[i].dudt = (FLOAT) 0.0;
  //  sphneib->UpdateAllSphDudt(sph);
  //}

  // Advance time variables
  n = n + 1;
  Nsteps = Nsteps + 1;
  t = t + timestep;

  // Advance SPH particles positions and velocities
  sphint->AdvanceParticles(n,level_step,sph->Nsph,
               sph->sphdata,(FLOAT) timestep);
  if (simparams->stringparams["gas_eos"] == "energy_eqn")
    uint->EnergyIntegration(n,level_step,sph->Nsph,
                sph->sphdata,(FLOAT) timestep);
  nbody->AdvanceParticles(n,nbody->Nstar,nbody->nbodydata,timestep);

  // Check all boundary conditions
  this->CheckBoundaries();

  // --------------------------------------------------------------------------
  if (sph->Nsph > 0) {

    // Reorder particles

    // Search ghost particles
    this->SearchGhostParticles();

    // Update neighbour tree
    sphneib->UpdateTree(sph,*simparams);
  }


  // --------------------------------------------------------------------------
  if (sph->Nsph > 0) {

    // Calculate all SPH properties
    sphneib->UpdateAllSphProperties(sph);

  }


  // Compute N-body forces
  // --------------------------------------------------------------------------
  if (nbody->Nstar > 0) {

    // Zero all acceleration terms
    for (i=0; i<nbody->Nnbody; i++) {
      if (nbody->nbodydata[i]->active) {
    for (k=0; k<ndim; k++) nbody->nbodydata[i]->a[k] = 0.0;
    for (k=0; k<ndim; k++) nbody->nbodydata[i]->adot[k] = 0.0;
    for (k=0; k<ndim; k++) nbody->nbodydata[i]->a2dot[k] = 0.0;
    for (k=0; k<ndim; k++) nbody->nbodydata[i]->a3dot[k] = 0.0;
    nbody->nbodydata[i]->gpot = 0.0;
      }
    }

    nbody->CalculateDirectGravForces(nbody->Nstar,nbody->nbodydata);
  }


  // --------------------------------------------------------------------------
  if (sph->Nsph > 0) {

    // Compute timesteps for all particles
    if (simparams->stringparams["sph"] == "godunov") {
      if (Nlevels == 1)
    this->ComputeGlobalTimestep();
      else
    this->ComputeBlockTimesteps();
    }

    if (simparams->stringparams["sph"] == "godunov")
      sphneib->UpdateAllSphDerivatives(sph);

    // Copy properties from original particles to ghost particles
    this->CopySphDataToGhosts();

    // Zero accelerations (perhaps)
    for (i=0; i<sph->Ntot; i++) {
      if (sph->sphdata[i].active) {
    for (k=0; k<ndim; k++) sph->sphdata[i].a[k] = (FLOAT) 0.0;
    for (k=0; k<ndim; k++) sph->sphdata[i].agrav[k] = (FLOAT) 0.0;
    sph->sphdata[i].gpot = (FLOAT) 0.0;
    sph->sphdata[i].dudt = (FLOAT) 0.0;
      }
    }

    // Calculate all SPH forces
    if (sph->hydro_forces == 1) sphneib->UpdateAllSphForces(sph);
    if (sph->self_gravity == 1) sphneib->UpdateAllSphGravForces(sph);

    // Add accelerations
    for (i=0; i<sph->Nsph; i++) {
      for (k=0; k<ndim; k++)
        sph->sphdata[i].a[k] += sph->sphdata[i].agrav[k];
    }

    if (simparams->stringparams["sph"] == "godunov")
      sphneib->UpdateAllSphDudt(sph);

  }

  // Apply correction steps for both particle and energy integration
  sphint->CorrectionTerms(n,level_step,sph->Nsph,
              sph->sphdata,(FLOAT) timestep);
  if (simparams->stringparams["gas_eos"] == "energy_eqn")
    uint->EnergyCorrectionTerms(n,level_step,sph->Nsph,
                sph->sphdata,(FLOAT) timestep);
  nbody->CorrectionTerms(n,nbody->Nnbody,nbody->nbodydata,timestep);

  // Set all end-of-step variables
  sphint->EndTimestep(n,level_step,sph->Nsph,sph->sphdata);
  if (simparams->stringparams["gas_eos"] == "energy_eqn")
    uint->EndTimestep(n,level_step,sph->Nsph,sph->sphdata);
  nbody->EndTimestep(n,nbody->Nnbody,nbody->nbodydata);

  return;
}

