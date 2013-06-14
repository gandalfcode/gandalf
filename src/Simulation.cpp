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

  CalculateDiagnostics();
  OutputDiagnostics();

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
  OutputDiagnostics();
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

  if (Nsteps%noutputstep == 0) cout << "t : " << t*simunits.t.outscale << " " 
				    << simunits.t.outunit << "    Nsteps : " 
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
  else if (simparams->stringparams["ic"] == "bb")
    BossBodenheimer();
  else if (simparams->stringparams["ic"] == "box")
    UniformBox();
  else if (simparams->stringparams["ic"] == "sphere")
    UniformSphere();
  else if (simparams->stringparams["ic"] == "cdiscontinuity")
    ContactDiscontinuity();
  else if (simparams->stringparams["ic"] == "noh")
    NohProblem();
  else if (simparams->stringparams["ic"] == "plummer")
    PlummerSphere();
  else if (simparams->stringparams["ic"] == "quadruple")
    QuadrupleStar();
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

  // Set-up all output units for scaling parameters
  if (intparams["dimensionless"] == 0) {
    simunits.r.outunit = stringparams["routunit"];
    simunits.m.outunit = stringparams["moutunit"];
    simunits.t.outunit = stringparams["toutunit"];
    simunits.v.outunit = stringparams["voutunit"];
    simunits.a.outunit = stringparams["aoutunit"];
    simunits.rho.outunit = stringparams["rhooutunit"];
    simunits.press.outunit = stringparams["pressoutunit"];
    simunits.f.outunit = stringparams["foutunit"];
    simunits.E.outunit = stringparams["Eoutunit"];
    simunits.mom.outunit = stringparams["momoutunit"];
    simunits.angmom.outunit = stringparams["angmomoutunit"];
    simunits.angvel.outunit = stringparams["angveloutunit"];
    simunits.dmdt.outunit = stringparams["dmdtoutunit"];
    simunits.u.outunit = stringparams["uoutunit"];
    simunits.dudt.outunit = stringparams["dudtoutunit"];
    simunits.temp.outunit = stringparams["tempoutunit"];
    simunits.SetupUnits(simparams);
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
				    floatparams["gamma_eos"],
                                    &simunits);
  else if (gas_eos == "barotropic")
    sph->eos = new Barotropic<ndim>(floatparams["temp0"],
				    floatparams["mu_bar"],
				    floatparams["gamma_eos"],
				    floatparams["rho_bary"],
                                    &simunits);
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
  else if (stringparams["neib_search"] == "tree")
    sphneib = new BinaryTree<ndim>(intparams["Nleafmax"],
                                   floatparams["thetamaxsqd"],
                                   stringparams["gravity_mac"]);
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
  else if (stringparams["nbody"] == "hermite4ts") {
    string KernelName = stringparams["kernel"];
    if (intparams["tabulated_kernel"] == 1) {
      nbody = new NbodyHermite4TS<ndim, TabulatedKernel>
	(intparams["nbody_softening"], intparams["sub_systems"],
	 floatparams["nbody_mult"], KernelName, intparams["Npec"]);
    }
    else if (intparams["tabulated_kernel"] == 0) {
      // Depending on the kernel, instantiate a different GradSph object
      if (KernelName == "m4") {
	nbody = new NbodyHermite4TS<ndim, M4Kernel>
	  (intparams["nbody_softening"], intparams["sub_systems"],
	   floatparams["nbody_mult"], KernelName, intparams["Npec"]);
      }
      else if (KernelName == "quintic") {
	nbody = new NbodyHermite4TS<ndim, QuinticKernel>
	  (intparams["nbody_softening"], intparams["sub_systems"],
	   floatparams["nbody_mult"], KernelName, intparams["Npec"]);
      }
      else if (KernelName == "gaussian") {
	nbody = new NbodyHermite4TS<ndim, GaussianKernel>
	  (intparams["nbody_softening"], intparams["sub_systems"],
	   floatparams["nbody_mult"], KernelName, intparams["Npec"]);
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
  simbox.boxmin[0] = floatparams["boxmin[0]"]/simunits.r.outscale;
  simbox.boxmin[1] = floatparams["boxmin[1]"]/simunits.r.outscale;
  simbox.boxmin[2] = floatparams["boxmin[2]"]/simunits.r.outscale;
  simbox.boxmax[0] = floatparams["boxmax[0]"]/simunits.r.outscale;
  simbox.boxmax[1] = floatparams["boxmax[1]"]/simunits.r.outscale;
  simbox.boxmax[2] = floatparams["boxmax[2]"]/simunits.r.outscale;
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
  tend = floatparams["tend"]/simunits.t.outscale;
  dt_snap = floatparams["dt_snap"]/simunits.t.outscale;
  noutputstep = intparams["noutputstep"];
  Nlevels = intparams["Nlevels"];
  sph_single_timestep = intparams["sph_single_timestep"];
  nbody_single_timestep = intparams["nbody_single_timestep"];
  sph->riemann_solver = stringparams["riemann_solver"];
  sph->slope_limiter = stringparams["slope_limiter"];
  sph->riemann_order = intparams["riemann_order"];
  nbody->Nstar = intparams["Nstar"];

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
    string msg = "Error: the function ProcessParameters has been already called!";
    ExceptionHandler::getIstance().raise(msg);
  }

  ProcessParameters();

  sph->AllocateMemory(sph->Nsph);

  nbody->AllocateMemory(nbody->Nstar);

  return;
}

//=============================================================================
//  SphSimulation::ImportArrayNbody
/// Import an array containing nbody particle properties from python to C++ arrays.
//=============================================================================
template <int ndim>
void Simulation<ndim>::ImportArrayNbody
(double* input,
    int size,
    string quantity)
{
    FLOAT StarParticle<ndim>::*quantityp; //Pointer to scalar quantity
    FLOAT (StarParticle<ndim>::*quantitypvec)[ndim]; //Pointer to component of vector quantity
    int index; //If it's a component of a vector quantity, we need to know its index
    bool scalar; //Is the requested quantity a scalar?

    //Check that the size is correct
    if (size != nbody->Nstar) {
      stringstream message;
      message << "Error: the array you are passing has a size of "
          << size << ", but memory has been allocated for "
          << nbody->Nstar << " star particles";
      ExceptionHandler::getIstance().raise(message.str());
    }

    // Now set pointer to the correct value inside the particle data structure
    // --------------------------------------------------------------------------
    if (quantity == "x") {
      quantitypvec = &StarParticle<ndim>::r;
      index = 0;
      scalar = false;
    }
    // --------------------------------------------------------------------------
    else if (quantity == "y") {
      if (ndim < 2) {
        string message = "Error: loading y-coordinate array for ndim < 2";
        ExceptionHandler::getIstance().raise(message);
      }
      quantitypvec = &StarParticle<ndim>::r;
      index = 1;
      scalar = false;
    }
    // --------------------------------------------------------------------------
    else if (quantity == "z") {
      if (ndim < 3) {
        string message = "Error: loading y-coordinate array for ndim < 3";
        ExceptionHandler::getIstance().raise(message);
      }
      quantitypvec = &StarParticle<ndim>::r;
      index = 2;
      scalar = false;
    }
    // --------------------------------------------------------------------------
    else if (quantity == "vx") {
      quantitypvec = &StarParticle<ndim>::v;
      index = 0;
      scalar = false;
    }
    // --------------------------------------------------------------------------
    else if (quantity == "vy") {
      if (ndim < 2) {
        string message = "Error: loading vy array for ndim < 2";
        ExceptionHandler::getIstance().raise(message);
      }
      quantitypvec = &StarParticle<ndim>::v;
      index = 1;
      scalar = false;
    }
    // --------------------------------------------------------------------------
    else if (quantity == "vz") {
      if (ndim < 3) {
        string message = "Error: loading vz array for ndim < 3";
        ExceptionHandler::getIstance().raise(message);
      }
      quantitypvec = &StarParticle<ndim>::v;
      index = 2;
      scalar = false;
    }
    // --------------------------------------------------------------------------
    else if (quantity=="m") {
      quantityp = &StarParticle<ndim>::m;
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
      for (StarParticle<ndim>* particlep = nbody->stardata;
       particlep < nbody->stardata+size; particlep++, i++) {
        particlep->*quantityp = input[i];
      }
    }
    else {
      int i=0;
      for (StarParticle<ndim>* particlep = nbody->stardata;
       particlep < nbody->stardata+size; particlep++, i++) {
        (particlep->*quantitypvec)[index] = input[i];
      }
    }

    return;
}


//=============================================================================
//  SphSimulation::ImportArraySph
/// Import an array containing sph particle properties from python to C++ arrays.
//=============================================================================
template <int ndim>
void Simulation<ndim>::ImportArraySph
(double* input,
    int size,
    string quantity)
{
  FLOAT SphParticle<ndim>::*quantityp; //Pointer to scalar quantity
  FLOAT (SphParticle<ndim>::*quantitypvec)[ndim]; //Pointer to component of vector quantity
  int index; //If it's a component of a vector quantity, we need to know its index
  bool scalar; //Is the requested quantity a scalar?

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
//  SphSimulation::ImportArray
/// Import an array containing particle properties from python to C++ arrays.
// This is a wrapper around ImportArraySph and ImportArrayNbody
//=============================================================================
template <int ndim>
void Simulation<ndim>::ImportArray
(double* input,                         ///< [in] Input array
 int size,                              ///< [in] Size of the input array
 string quantity,                       ///< [in] Which quantity should be set equal to the given array
 string type)                       ///< [in] Which particle type should be assigned the array
{

  debug2("[SphSimulation::ImportArray]");

  //Check that PreSetup has been called
  if (! ParametersProcessed) {
    string msg = "Error: before calling ImportArray, you need to call PreSetupForPython!";
    ExceptionHandler::getIstance().raise(msg);
  }

  //Call the right function depending on the passed in type
  if (type=="sph") {
    //Check sph has been allocated
    if (sph==NULL) {
      string message = "Error: memory for sph was not allocated! Are you sure that this is not a nbody-only simulation?";
      ExceptionHandler::getIstance().raise(message);
    }
    ImportArraySph(input, size, quantity);

  }
  else if (type=="star") {
    if (nbody==NULL) {
      string message = "Error: memory for nbody was not allocated! Are you sure that this is not a sph-only simulation?";
      ExceptionHandler::getIstance().raise(message);
    }
    ImportArrayNbody(input, size, quantity);
  }
  else {
    string message = "Error: we did not recognize the type " + type + ", the only allowed types are \"sph\""
        " and \"nbody\"";
    ExceptionHandler::getIstance().raise(message);
  }



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



template <int ndim>
void SphSimulation<ndim>::ProcessParameters()
{
  Simulation<ndim>::ProcessParameters();
}



template <int ndim>
void GodunovSimulation<ndim>::ProcessParameters()
{
  Simulation<ndim>::ProcessParameters();
}
