//=================================================================================================
//  Simulation.cpp
//  Contains all main functions controlling the simulation work-flow.
//
//  This file is part of GANDALF :
//  Graphical Astrophysics code for N-body Dynamics And Lagrangian Fluids
//  https://github.com/gandalfcode/gandalf
//  Contact : gandalfcode@gmail.com
//
//  Copyright (C) 2013  D. A. Hubber, G. Rosotti
//
//  GANDALF is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 2 of the License, or
//  (at your option) any later version.
//
//  GANDALF is distributed in the hope that it will be useful, but
//  WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  General Public License (http://www.gnu.org/licenses) for more details.
//=================================================================================================


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
#include "Debug.h"
#include "Simulation.h"
#include "Parameters.h"
#include "InlineFuncs.h"
#include "Nbody.h"
#include "Hydrodynamics.h"
#include "Sph.h"
#include "RiemannSolver.h"
#include "SimulationIO.hpp"
#include "SimulationIC.hpp"
#include "SimAnalysis.hpp"
#include "SphSnapshot.h"
using namespace std;


// Declare invndim constant here (prevents warnings with some compilers)
template <int ndim>
const FLOAT Simulation<ndim>::invndim = 1.0/ndim;


//=================================================================================================
//  SimulationBase::SimulationFactory
/// Creates a simulation object depending on the dimensionality.
//=================================================================================================
SimulationBase* SimulationBase::SimulationFactory
 (int ndim,                            ///< [in] No. of dimensions
  string simtype,                      ///< [in] Simulation type
  Parameters* params)                  ///< [in] Pointer to parameters object
{
  debug1("[SimulationBase::SimulationFactory]");

  // Check ndim is valid
  if (ndim < 1 || ndim > 3) {
    stringstream msg;
    msg << "Error: ndim must be either 1, 2, 3; the value " << ndim << "is not allowed!";
    ExceptionHandler::getIstance().raise(msg.str());
  }

  // Check simulation type is valid
  if (simtype != "sph" && simtype != "gradhsph" && simtype != "sm2012sph" &&
      simtype != "meshlessfv" && simtype != "mfvmuscl" && simtype != "mfvrk" &&
      simtype != "nbody" ) {
    string msg = "Error: the simulation type " + simtype + " was not recognized";
    ExceptionHandler::getIstance().raise(msg);
  }


  // Set ndim and simtype inside the parameters
  params->intparams["ndim"] = ndim;
  params->stringparams["sim"] = simtype;


  // Create and return Simulation object depending on the chosen algorithm
  // and the dimensionality.
  if (ndim == 1) {
    if (simtype == "gradhsph" || simtype == "sph") {
      return new GradhSphSimulation<1>(params);
    }
    else if (simtype == "sm2012sph") {
      return new SM2012SphSimulation<1>(params);
    }
    else if (simtype == "meshlessfv" || simtype == "mfvmuscl") {
      return new MfvMusclSimulation<1>(params);
    }
    else if (simtype == "mfvrk") {
      return new MfvRungeKuttaSimulation<1>(params);
    }
    else if (simtype == "nbody") {
      return new NbodySimulation<1>(params);
    }
  }
  else if (ndim == 2) {
    if (simtype == "gradhsph" || simtype == "sph") {
      return new GradhSphSimulation<2>(params);
    }
    else if (simtype == "sm2012sph") {
      return new SM2012SphSimulation<2>(params);
    }
    else if (simtype == "meshlessfv" || simtype == "mfvmuscl") {
      return new MfvMusclSimulation<2>(params);
    }
    else if (simtype == "mfvrk") {
      return new MfvRungeKuttaSimulation<2>(params);
    }
    else if (simtype == "nbody") {
      return new NbodySimulation<2>(params);
    }
  }
  else if (ndim == 3) {
    if (simtype == "gradhsph" || simtype == "sph") {
      return new GradhSphSimulation<3>(params);
    }
    else if (simtype == "sm2012sph") {
      return new SM2012SphSimulation<3>(params);
    }
    else if (simtype == "meshlessfv" || simtype == "mfvmuscl") {
      return new MfvMusclSimulation<3>(params);
    }
    else if (simtype == "mfvrk") {
      return new MfvRungeKuttaSimulation<3>(params);
    }
    else if (simtype == "nbody") {
      return new NbodySimulation<3>(params);
    }
  }
  return NULL;
}



//=================================================================================================
//  SimulationBase::SimulationBase
/// SimulationBase constructor, initialising important simulation variables.
//=================================================================================================
SimulationBase::SimulationBase
 (Parameters* params)                ///< [in] Pointer to parameters object
{
  simparams = new Parameters(*params);
  paramfile             = "";
  integration_step      = 1;
  litesnap              = 0;
  n                     = 0;
  nlastrestart          = 0;
  nrestartstep          = 0;
  nresync               = 0;
  Nblocksteps           = 0;
  Nfullsteps            = 0;
  Nmpi                  = 1;
  Noutsnap              = 0;
  Noutlitesnap          = 0;
  Nsteps                = 0;
  rank                  = 0;
  dt_snap_wall          = 0.0;
  t                     = 0.0;
  timestep              = 0.0;
  tsnaplast             = 0.0;
  tlitesnaplast         = 0.0;
  tsnap_wallclock       = 0.0;
  ewaldGravity          = false;
  initial_h_provided    = false;
  kill_simulation       = false;
  ParametersProcessed   = false;
  periodicBoundaries    = false;
  rescale_particle_data = false;
  restart               = false;
  setup                 = false;
#if defined _OPENMP
  if (omp_get_dynamic()) {
    cout << "Warning: the dynamic adjustment of the number threads was on. "
         << "For better load-balancing, we will disable it" << endl;
  }
  omp_set_dynamic(0);
  Nthreads = omp_get_max_threads();
  assert(Nthreads > 0);
#else
  Nthreads = 1;
#endif
}



//=================================================================================================
//  SimulationBase::~SimulationBase
/// SimulationBase destructor
//=================================================================================================
SimulationBase::~SimulationBase()
{
}



//=================================================================================================
//  SimulationBase::SplashScreen
/// Write splash screen to standard output.
//=================================================================================================
void SimulationBase::SplashScreen(void)
{
  cout << "******************************************************************************" << endl;
  cout << "*                                                                            *" << endl;
  cout << "*         *****     ****    *     *   *****     ****    *      ******        *" << endl;
  cout << "*        *     *   *    *   **    *   *    *   *    *   *      *             *" << endl;
  cout << "*        *         *    *   * *   *   *    *   *    *   *      *             *" << endl;
  cout << "*        *    **   ******   *  *  *   *    *   ******   *      ******        *" << endl;
  cout << "*        *     *   *    *   *   * *   *    *   *    *   *      *             *" << endl;
  cout << "*        *     *   *    *   *    **   *    *   *    *   *      *             *" << endl;
  cout << "*         *****    *    *   *     *   *****    *    *   *****  *             *" << endl;
  cout << "*                                                                            *" << endl;
  cout << "*   Graphical Astrophysics code for N-body Dynamics And Lagrangian Fluids    *" << endl;
  cout << "*                        Version 0.4.0 - 04/09/2015                          *" << endl;
  cout << "*                                                                            *" << endl;
  cout << "*                 Original code : D. A. Hubber & G. Rosotti                  *" << endl;
  cout << "*                                                                            *" << endl;
  cout << "*              Contributions by : S. Balfour, F. Dinnbier, S. Heigl,         *" << endl;
  cout << "*                                 O. Lomax, J. Ngoumou, P. Rohde,            *" << endl;
  cout << "*                                 S. Walch, A. P. Whitworth, R. Wunsch       *" << endl;
  cout << "*                                                                            *" << endl;
  cout << "*                  https://github.com/gandalfcode/gandalf                    *" << endl;
  cout << "*                                                                            *" << endl;
  cout << "******************************************************************************" << endl;

  return;
}



//=================================================================================================
//  SimulationBase::SetParam
/// Accessor function for modifying a string value. Also checks that the
/// non return point has not been reached
//=================================================================================================
void SimulationBase::SetParam
 (string key,                          ///< [in] Parameter string name
  string value)                        ///< [in] Parameter string value
{
  // Error checking
  if (ParametersProcessed) {
    string msg = "Error: the non-return point for setting parameters has been reached!";
    ExceptionHandler::getIstance().raise(msg);
  }
  if (key == "ndim") {
    string msg = "Error: Not possible to change the number of dimensions!";
    ExceptionHandler::getIstance().raise(msg);
  }

  simparams->SetParameter(key, value);
}



//=================================================================================================
//  SphSimulationBase::SetParam
/// Accessor function for modifying an int value, wrapper around the one for
/// string value. Also checks that the non return point has not been reached
//=================================================================================================
void SimulationBase::SetParam
 (string key,                          ///< [in] Parameter string name
  int value)                           ///< [in] Parameter integer value
{
  ostringstream convert;
  convert << value;
  SetParam (key, convert.str());
}



//=================================================================================================
//  SimulationBase::SetParam
/// Accessor function for modifying a float value, wrapper around the one for
/// string value.  Also checks that the non return point has not been reached.
//=================================================================================================
void SimulationBase::SetParam
 (string key,                          ///< [in] Parameter string name
  double value)                        ///< [in] Parameter float (double) value
{
  ostringstream convert;
  convert << value;
  SetParam (key, convert.str());
}



//=================================================================================================
//  SimulationBase::GetParam
/// Accessor function for getting a parameter value
/// Wrapper around the corresponding function in Parameters
//=================================================================================================
string SimulationBase::GetParam(string key)
{
  return simparams->GetParameter(key);
}



//=================================================================================================
//  SimulationBase::GetIntAndFloatParameterKeys
/// Returns a list containing the keys of all the int and float parameters
//=================================================================================================
std::list<string>* SimulationBase::GetIntAndFloatParameterKeys()
{
  if (! keys.empty()) return &keys;

  for (std::map<string, int>::iterator it=simparams->intparams.begin() ;
       it != simparams->intparams.end(); it++) {
    keys.push_back(it->first);
  }

  for (std::map<string, double>::iterator it=simparams->floatparams.begin() ;
       it != simparams->floatparams.end(); it++) {
    keys.push_back(it->first);
  }

  return &keys;
}



//=================================================================================================
//  SimulationBase::Run
/// Controls the simulation main loop, including exit conditions.  If provided as an optional
/// argument, will only advance the simulation by 'Nadvance' steps.
//=================================================================================================
void SimulationBase::Run
 (int Nadvance)                        ///< [in] Selected max no. of timesteps (Optional argument)
{
  int Ntarget;                         // Target step no before finishing main code integration.

  debug1("[SimulationBase::Run]");

  // Set integer timestep exit condition if provided as parameter.
  if (Nadvance < 0) Ntarget = Nstepsmax;
  else Ntarget = (int) (Nsteps + Nadvance);

  // Continue to run simulation until we reach the required time, or
  // exeeded the maximum allowed number of steps.
  //-----------------------------------------------------------------------------------------------
  while (t < tend && Nsteps < Ntarget) {

    timing->StartTimingSection("RUN");

    MainLoop();
    Output();

    // Special condition to check if maximum wall-clock time has been reached.
    if (kill_simulation || timing->WallClockTime() - timing->tstart_wall > 0.95*tmax_wallclock) {
      RestartSnapshot();
      cout << "Reached maximum wall-clock time.  Killing simulation." << endl;
      break;
    }

    timing->EndTimingSection("RUN");

  }
  //-----------------------------------------------------------------------------------------------

  FinaliseSimulation();
  CalculateDiagnostics();
  OutputDiagnostics();
  UpdateDiagnostics();

  cout << "Final t : " << t*simunits.t.outscale << " " << simunits.t.outunit
       << "    Total no. of steps : " << Nsteps << endl;


  // If reached end of simulation, remove 'cont' file to prevent automatic restart
  if (t >= tend || Nsteps >= Ntarget) {
    if (remove("cont") != 0) {
      cout << "Error deleting cont file" << endl;
    }
  }

  return;
}



//=================================================================================================
//  SimulationBase::InteractiveRun
/// Controls the simulation main loop, including exit conditions.
/// If provided, will only advance the simulation by 'Nadvance' steps.
//=================================================================================================
list<SphSnapshotBase*> SimulationBase::InteractiveRun
 (int Nadvance)                        ///< [in] Max no. of integer steps (Optional argument)
{
  int Ntarget;                // Selected integer timestep
  DOUBLE tdiff = 0.0;                  // Measured time difference
  clock_t tstart = clock();            // Initial CPU clock time
  string filename;                     // Name of the output file
  list<SphSnapshotBase*> snap_list;    // List of snapshots produced while running
                                       // that will be passed back to Python

  debug2("[SimulationBase::InteractiveRun]");

  // Set integer timestep exit condition if provided as parameter.
  if (Nadvance < 0) Ntarget = Nstepsmax;
  else Ntarget = (int) (Nsteps + Nadvance);

  // Continue to run simulation until we reach the required time, or
  // exeeded the maximum allowed number of steps.
  //-----------------------------------------------------------------------------------------------
  while (t < tend && Nsteps < Ntarget && tdiff < dt_python) {

    // Evolve the simulation one step
    MainLoop();

    // Update all diagnostics (including binaries) here for now
    if (t >= tsnapnext) CalculateDiagnostics();

    // Call output routine
    filename = Output();

    // If we have written a snapshot, create a new snapshot object
    if (filename.length() != 0) {
      SphSnapshotBase* snapshot =
        SphSnapshotBase::SphSnapshotFactory(filename, this, ndims);
      snapshot->CopyDataFromSimulation();
      snap_list.push_back(snapshot);
    }

    // Measure CPU clock time difference since current function was called
    tdiff = (DOUBLE) (clock() - tstart) / (DOUBLE) CLOCKS_PER_SEC;

  }
  //-----------------------------------------------------------------------------------------------


  // Calculate and process all diagnostic quantities
  if (t >= tend || Nsteps >= Ntarget) {
    FinaliseSimulation();
    CalculateDiagnostics();
    OutputDiagnostics();
    UpdateDiagnostics();
  }

  return snap_list;
}



//=================================================================================================
//  SimulationBase::Output
/// Controls when regular output snapshots are written by the code.
//=================================================================================================
string SimulationBase::Output(void)
{
  string filename;                  // 'Lite' output snapshot filename
  string filename2;                 // Regular output snapshot filename
  string nostring;                  // String of number of snapshots
  string fileend;                   // Name of restart file
  stringstream ss;                  // Stream object for preparing filename
  ofstream outfile;                 // Stream of restart file

  debug2("[SimulationBase::Output]");
  timing->StartTimingSection("OUTPUT");


  // Output time and no of steps for root process
  if (rank == 0) {
    if (Nsteps%noutputstep == 0) {
      cout << "t : " << t*simunits.t.outscale << " " << simunits.t.outunit
           << "    dt : " << timestep*simunits.t.outscale << " "
           << simunits.t.outunit << "    Nsteps : " << Nsteps << endl;
    }
  }


  // Output a lite-data snapshot for producing movies
  //-----------------------------------------------------------------------------------------------
  if (litesnap == 1 && t >= tlitesnapnext) {

    // Prepare filename for new snapshot
    Noutlitesnap++;
    tlitesnaplast = tlitesnapnext;
    tlitesnapnext += dt_litesnap;
    nostring = "";
    ss << setfill('0') << setw(5) << Noutlitesnap;
    nostring = ss.str();
    filename = run_id + ".slite." + nostring;
    ss.str(std::string());
    WriteSnapshotFile(filename,"slite");

  }
  //-----------------------------------------------------------------------------------------------

  // Output a data snapshot if reached required time
  //-----------------------------------------------------------------------------------------------
  if (t >= tsnapnext ) {

    // Prepare filename for new snapshot
    Noutsnap++;
    tsnaplast = tsnapnext;
    tsnapnext += dt_snap;
    nostring = "";
    ss << setfill('0') << setw(5) << Noutsnap;
    nostring = ss.str();
    filename = run_id + '.' + out_file_form + '.' + nostring;
    ss.str(std::string());
    WriteSnapshotFile(filename,out_file_form);

    // Now write name and format of snapshot to file (for restarts)
    if (rank == 0) {

      fileend = "restart";
      filename2 = run_id + "." + fileend;
      outfile.open(filename2.c_str());
      outfile << out_file_form << endl;
      outfile << filename << endl;
      outfile.close();

      // Finally, calculate wall-clock time interval since last output snapshot
      if (tsnap_wallclock > 0.0) dt_snap_wall = timing->WallClockTime() - tsnap_wallclock;
      tsnap_wallclock = timing->WallClockTime();

      // If simulation is too close to maximum wall-clock time, end prematurely
      if (timing->ttot_wall > 0.95*tmax_wallclock) {
        kill_simulation = true;
      }

    }

  }
  //-----------------------------------------------------------------------------------------------


  // Output diagnostics to screen if passed sufficient number of block steps
  if (Nblocksteps%ndiagstep == 0 && n == nresync) {
    CalculateDiagnostics();
    OutputDiagnostics();
    UpdateDiagnostics();
    timing->ComputeTimingStatistics(run_id);

  }

  // Create temporary snapshot file
  if (n == nresync && Nsteps - nlastrestart >= nrestartstep) {
    RestartSnapshot();
    nlastrestart = Nsteps;
  }


  timing->EndTimingSection("OUTPUT");

  return filename;
}



//=================================================================================================
//  SimulationBase::RestartSnapshot
/// Write the restart log file (containing the last snapshot i.d.) plus a temporary snapshot
/// file for future restarting.
//=================================================================================================
void SimulationBase::RestartSnapshot(void)
{
  string filename;                     // Temporary output snapshot filename
  string filename2;                    // Restart log filename
  stringstream ss;                     // Stream object for preparing filename
  ofstream outfile;                    // Stream of restart file

  debug2("[SimulationBase::RestartSnapshot]");

  // Prepare filename for new snapshot
  filename = run_id + "." + out_file_form + ".tmp";
  ss.str(std::string());
  WriteSnapshotFile(filename,out_file_form);

  // Now write name and format of snapshot to file (for restarts)
  filename2 = run_id + ".restart";
  outfile.open(filename2.c_str());
  outfile << out_file_form << endl;
  outfile << filename << endl;
  outfile.close();

  return;
}



//=================================================================================================
//  SimulationBase::SetupSimulation
/// Main function for setting up a new simulation.
//=================================================================================================
void SimulationBase::SetupSimulation(void)
{
  debug1("[SimulationBase::Setup]");

  timing->StartTimingSection("SETUP");

  if (setup) {
    string msg = "This simulation has been already set up";
    ExceptionHandler::getIstance().raise(msg);
  }

  // Process the parameters file setting up all simulation objects
  if (simparams->stringparams["ic"] == "python") {
    if (!ParametersProcessed) {
      string msg = "Error: you are attempting to setup a simulation with initial conditions "
                   "generated from Python. Before setting up the simulation, you need to "
                   "import the initial conditions";
      ExceptionHandler::getIstance().raise(msg);
    }
  }
  else {
    if (ParametersProcessed) {
      string msg = "The parameters of the simulation have been already processed. It means that "
                   "you shouldn't be calling this function, please consult the documentation.";
      ExceptionHandler::getIstance().raise(msg);
    }
    ProcessParameters();
  }

  // Generate initial conditions for simulation on root process (for MPI jobs)
  if (rank == 0) {
    GenerateIC();
  }

  // Change to COM frame if selected
  if (simparams->intparams["com_frame"] == 1) SetComFrame();

  // Perform the rest of the initialisation, calculating all initial particle
  // quantities and setting up trees.
  PostInitialConditionsSetup();

  // Initial output before simulation begins
  Output();

  timing->EndTimingSection("SETUP");

  return;
}



//=================================================================================================
//  Simulation::ProcessNbodyParameters
/// Process all the options chosen in the parameters file for setting up
/// objects related to stars and N-body integration.
//=================================================================================================
template <int ndim>
void Simulation<ndim>::ProcessNbodyParameters(void)
{
  map<string, int> &intparams = simparams->intparams;
  map<string, double> &floatparams = simparams->floatparams;
  map<string, string> &stringparams = simparams->stringparams;

  // Create N-body object based on chosen method in params file
  //-----------------------------------------------------------------------------------------------
  if (stringparams["nbody"] == "lfkdk") {
    string KernelName = stringparams["kernel"];
    if (intparams["tabulated_kernel"] == 1) {
      nbody = new NbodyLeapfrogKDK<ndim, TabulatedKernel>
        (intparams["nbody_softening"], intparams["perturbers"], intparams["sub_systems"],
         floatparams["nbody_mult"], KernelName);
    }
    else if (intparams["tabulated_kernel"] == 0) {
      if (KernelName == "m4") {
        nbody = new NbodyLeapfrogKDK<ndim, M4Kernel>
          (intparams["nbody_softening"], intparams["perturbers"], intparams["sub_systems"],
           floatparams["nbody_mult"], KernelName);
      }
      else if (KernelName == "quintic") {
        nbody = new NbodyLeapfrogKDK<ndim, QuinticKernel>
          (intparams["nbody_softening"], intparams["perturbers"], intparams["sub_systems"],
           floatparams["nbody_mult"], KernelName);
      }
      else if (KernelName == "gaussian") {
        nbody = new NbodyLeapfrogKDK<ndim, GaussianKernel>
          (intparams["nbody_softening"], intparams["perturbers"], intparams["sub_systems"],
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
  //-----------------------------------------------------------------------------------------------
  else if (stringparams["nbody"] == "lfdkd") {
    string KernelName = stringparams["kernel"];
    if (intparams["tabulated_kernel"] == 1) {
      nbody = new NbodyLeapfrogDKD<ndim, TabulatedKernel>
        (intparams["nbody_softening"], intparams["perturbers"], intparams["sub_systems"],
         floatparams["nbody_mult"], KernelName);
    }
    else if (intparams["tabulated_kernel"] == 0) {
      if (KernelName == "m4") {
        nbody = new NbodyLeapfrogDKD<ndim, M4Kernel>
          (intparams["nbody_softening"], intparams["perturbers"], intparams["sub_systems"],
           floatparams["nbody_mult"], KernelName);
      }
      else if (KernelName == "quintic") {
        nbody = new NbodyLeapfrogDKD<ndim, QuinticKernel>
          (intparams["nbody_softening"], intparams["perturbers"], intparams["sub_systems"],
           floatparams["nbody_mult"], KernelName);
      }
      else if (KernelName == "gaussian") {
        nbody = new NbodyLeapfrogDKD<ndim, GaussianKernel>
          (intparams["nbody_softening"], intparams["perturbers"], intparams["sub_systems"],
           floatparams["nbody_mult"], KernelName);
      }
      else {
        string message = "Unrecognised parameter : kernel = " + simparams->stringparams["kernel"];
        ExceptionHandler::getIstance().raise(message);
      }
    }
    else {
      string message = "Invalid option for the tabulated_kernel parameter: " +
        stringparams["tabulated_kernel"];
      ExceptionHandler::getIstance().raise(message);
    }
    integration_step = max(integration_step,2);
  }
  //-----------------------------------------------------------------------------------------------
  else if (stringparams["nbody"] == "hermite4") {
    string KernelName = stringparams["kernel"];
    if (intparams["tabulated_kernel"] == 1) {
      nbody = new NbodyHermite4<ndim, TabulatedKernel>
        (intparams["nbody_softening"], intparams["perturbers"], intparams["sub_systems"],
         floatparams["nbody_mult"], KernelName);
    }
    else if (intparams["tabulated_kernel"] == 0) {
      if (KernelName == "m4") {
        nbody = new NbodyHermite4<ndim, M4Kernel>
          (intparams["nbody_softening"], intparams["perturbers"], intparams["sub_systems"],
           floatparams["nbody_mult"], KernelName);
      }
      else if (KernelName == "quintic") {
        nbody = new NbodyHermite4<ndim, QuinticKernel>
          (intparams["nbody_softening"], intparams["perturbers"], intparams["sub_systems"],
           floatparams["nbody_mult"], KernelName);
      }
      else if (KernelName == "gaussian") {
        nbody = new NbodyHermite4<ndim, GaussianKernel>
          (intparams["nbody_softening"], intparams["perturbers"], intparams["sub_systems"],
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
  //-----------------------------------------------------------------------------------------------
  else if (stringparams["nbody"] == "hermite4ts") {
    string KernelName = stringparams["kernel"];
    if (intparams["tabulated_kernel"] == 1) {
      nbody = new NbodyHermite4TS<ndim, TabulatedKernel>
        (intparams["nbody_softening"], intparams["perturbers"], intparams["sub_systems"],
         floatparams["nbody_mult"], KernelName, intparams["Npec"]);
    }
    else if (intparams["tabulated_kernel"] == 0) {
      if (KernelName == "m4") {
        nbody = new NbodyHermite4TS<ndim, M4Kernel>
          (intparams["nbody_softening"], intparams["perturbers"], intparams["sub_systems"],
           floatparams["nbody_mult"], KernelName, intparams["Npec"]);
      }
      else if (KernelName == "quintic") {
        nbody = new NbodyHermite4TS<ndim, QuinticKernel>
          (intparams["nbody_softening"], intparams["perturbers"], intparams["sub_systems"],
           floatparams["nbody_mult"], KernelName, intparams["Npec"]);
      }
      else if (KernelName == "gaussian") {
        nbody = new NbodyHermite4TS<ndim, GaussianKernel>
          (intparams["nbody_softening"], intparams["perturbers"], intparams["sub_systems"],
           floatparams["nbody_mult"], KernelName, intparams["Npec"]);
      }
      else {
        string message = "Unrecognised parameter : kernel = " + simparams->stringparams["kernel"];
        ExceptionHandler::getIstance().raise(message);
      }
    }
    else {
      string message = "Invalid option for the tabulated_kernel parameter: " +
        stringparams["tabulated_kernel"];
      ExceptionHandler::getIstance().raise(message);
    }
  }
  //-----------------------------------------------------------------------------------------------
  else if (stringparams["nbody"] == "hermite6ts") {
    string KernelName = stringparams["kernel"];
    if (intparams["tabulated_kernel"] == 1) {
      nbody = new NbodyHermite6TS<ndim, TabulatedKernel>
        (intparams["nbody_softening"], intparams["perturbers"], intparams["sub_systems"],
         floatparams["nbody_mult"], KernelName, intparams["Npec"]);
    }
    else if (intparams["tabulated_kernel"] == 0) {
      if (KernelName == "m4") {
        nbody = new NbodyHermite6TS<ndim, M4Kernel>
          (intparams["nbody_softening"], intparams["perturbers"], intparams["sub_systems"],
           floatparams["nbody_mult"], KernelName, intparams["Npec"]);
      }
      else if (KernelName == "quintic") {
        nbody = new NbodyHermite6TS<ndim, QuinticKernel>
          (intparams["nbody_softening"], intparams["perturbers"], intparams["sub_systems"],
           floatparams["nbody_mult"], KernelName, intparams["Npec"]);
      }
      else if (KernelName == "gaussian") {
        nbody = new NbodyHermite6TS<ndim, GaussianKernel>
          (intparams["nbody_softening"], intparams["perturbers"], intparams["sub_systems"],
           floatparams["nbody_mult"], KernelName, intparams["Npec"]);
      }
      else {
        string message = "Unrecognised parameter : kernel = " + simparams->stringparams["kernel"];
        ExceptionHandler::getIstance().raise(message);
      }
    }
    else {
      string message = "Invalid option for the tabulated_kernel parameter: " +
        stringparams["tabulated_kernel"];
      ExceptionHandler::getIstance().raise(message);
    }
  }
  //-----------------------------------------------------------------------------------------------
  else {
    string message = "Unrecognised parameter : nbody = "
      + simparams->stringparams["nbody"];
    ExceptionHandler::getIstance().raise(message);
  }
  //-----------------------------------------------------------------------------------------------


  // Create sub-system object based on chosen method in params file
  //-----------------------------------------------------------------------------------------------
  if (intparams["sub_systems"] == 1) {

    //---------------------------------------------------------------------------------------------
    if (stringparams["sub_system_integration"] == "lfkdk") {
      string KernelName = stringparams["kernel"];
      if (intparams["tabulated_kernel"] == 1) {
        subsystem = new NbodyLeapfrogKDK<ndim, TabulatedKernel>
          (intparams["nbody_softening"], intparams["perturbers"], intparams["sub_systems"],
           floatparams["subsys_mult"], KernelName);
      }
      else if (intparams["tabulated_kernel"] == 0) {
        if (KernelName == "m4") {
          subsystem = new NbodyLeapfrogKDK<ndim, M4Kernel>
            (intparams["nbody_softening"], intparams["perturbers"], intparams["sub_systems"],
             floatparams["subsys_mult"], KernelName);
        }
        else if (KernelName == "quintic") {
          subsystem = new NbodyLeapfrogKDK<ndim, QuinticKernel>
            (intparams["nbody_softening"], intparams["perturbers"], intparams["sub_systems"],
             floatparams["subsys_mult"], KernelName);
        }
        else if (KernelName == "gaussian") {
          subsystem = new NbodyLeapfrogKDK<ndim, GaussianKernel>
            (intparams["nbody_softening"], intparams["perturbers"], intparams["sub_systems"],
             floatparams["subsys_mult"], KernelName);
        }
        else {
          string message = "Unrecognised parameter : kernel = " + simparams->stringparams["kernel"];
          ExceptionHandler::getIstance().raise(message);
        }
      }
      else {
        string message = "Invalid option for the tabulated_kernel parameter: "
          + stringparams["tabulated_kernel"];
        ExceptionHandler::getIstance().raise(message);
      }
    }
    //---------------------------------------------------------------------------------------------
    else if (stringparams["sub_system_integration"] == "hermite4") {
      string KernelName = stringparams["kernel"];
      if (intparams["tabulated_kernel"] == 1) {
        subsystem = new NbodyHermite4<ndim, TabulatedKernel>
          (intparams["nbody_softening"], intparams["perturbers"], intparams["sub_systems"],
           floatparams["subsys_mult"], KernelName);
      }
      else if (intparams["tabulated_kernel"] == 0) {
        if (KernelName == "m4") {
          subsystem = new NbodyHermite4<ndim, M4Kernel>
            (intparams["nbody_softening"], intparams["perturbers"], intparams["sub_systems"],
             floatparams["subsys_mult"], KernelName);
        }
        else if (KernelName == "quintic") {
          subsystem = new NbodyHermite4<ndim, QuinticKernel>
            (intparams["nbody_softening"], intparams["perturbers"], intparams["sub_systems"],
             floatparams["subsys_mult"], KernelName);
        }
        else if (KernelName == "gaussian") {
          subsystem = new NbodyHermite4<ndim, GaussianKernel>
            (intparams["nbody_softening"], intparams["perturbers"], intparams["sub_systems"],
             floatparams["subsys_mult"], KernelName);
        }
        else {
          string message = "Unrecognised parameter : kernel = " + simparams->stringparams["kernel"];
          ExceptionHandler::getIstance().raise(message);
        }
      }
      else {
        string message = "Invalid option for the tabulated_kernel parameter: "
          + stringparams["tabulated_kernel"];
        ExceptionHandler::getIstance().raise(message);
      }
    }
    //---------------------------------------------------------------------------------------------
    else if (stringparams["sub_system_integration"] == "hermite4ts") {
      string KernelName = stringparams["kernel"];
      if (intparams["tabulated_kernel"] == 1) {
        subsystem = new NbodyHermite4TS<ndim, TabulatedKernel>
          (intparams["nbody_softening"], intparams["perturbers"], intparams["sub_systems"],
           floatparams["subsys_mult"], KernelName, intparams["Npec"]);
      }
      else if (intparams["tabulated_kernel"] == 0) {
        if (KernelName == "m4") {
          subsystem = new NbodyHermite4TS<ndim, M4Kernel>
            (intparams["nbody_softening"], intparams["perturbers"], intparams["sub_systems"],
             floatparams["subsys_mult"], KernelName, intparams["Npec"]);
        }
        else if (KernelName == "quintic") {
          subsystem = new NbodyHermite4TS<ndim, QuinticKernel>
            (intparams["nbody_softening"], intparams["perturbers"], intparams["sub_systems"],
             floatparams["subsys_mult"], KernelName, intparams["Npec"]);
        }
        else if (KernelName == "gaussian") {
          subsystem = new NbodyHermite4TS<ndim, GaussianKernel>
            (intparams["nbody_softening"], intparams["perturbers"], intparams["sub_systems"],
             floatparams["subsys_mult"], KernelName, intparams["Npec"]);
        }
        else {
          string message = "Unrecognised parameter : kernel = " + simparams->stringparams["kernel"];
          ExceptionHandler::getIstance().raise(message);
        }
      }
      else {
        string message = "Invalid option for the tabulated_kernel parameter: "
          + stringparams["tabulated_kernel"];
        ExceptionHandler::getIstance().raise(message);
      }
    }
    //---------------------------------------------------------------------------------------------
    else if (stringparams["sub_system_integration"] == "hermite6ts") {
      string KernelName = stringparams["kernel"];
      if (intparams["tabulated_kernel"] == 1) {
        subsystem = new NbodyHermite6TS<ndim, TabulatedKernel>
          (intparams["nbody_softening"], intparams["perturbers"], intparams["sub_systems"],
           floatparams["subsys_mult"], KernelName, intparams["Npec"]);
      }
      else if (intparams["tabulated_kernel"] == 0) {
        if (KernelName == "m4") {
          subsystem = new NbodyHermite6TS<ndim, M4Kernel>
            (intparams["nbody_softening"], intparams["perturbers"], intparams["sub_systems"],
             floatparams["subsys_mult"], KernelName, intparams["Npec"]);
        }
        else if (KernelName == "quintic") {
          subsystem = new NbodyHermite6TS<ndim, QuinticKernel>
            (intparams["nbody_softening"], intparams["perturbers"], intparams["sub_systems"],
             floatparams["subsys_mult"], KernelName, intparams["Npec"]);
        }
        else if (KernelName == "gaussian") {
          subsystem = new NbodyHermite6TS<ndim, GaussianKernel>
            (intparams["nbody_softening"], intparams["perturbers"], intparams["sub_systems"],
             floatparams["subsys_mult"], KernelName, intparams["Npec"]);
        }
        else {
          string message = "Unrecognised parameter : kernel = " + simparams->stringparams["kernel"];
          ExceptionHandler::getIstance().raise(message);
        }
      }
      else {
        string message = "Invalid option for the tabulated_kernel parameter: "
          + stringparams["tabulated_kernel"];
        ExceptionHandler::getIstance().raise(message);
      }
    }
    //---------------------------------------------------------------------------------------------
    else {
      string message = "Unrecognised parameter : sub_system_integration = "
        + simparams->stringparams["sub_system_integration"];
      ExceptionHandler::getIstance().raise(message);
    }
    //---------------------------------------------------------------------------------------------

  }
  //-----------------------------------------------------------------------------------------------

  return;
}



//=================================================================================================
//  Simulation::AllocateParticleMemory
/// Allocate all memory for both SPH and N-body particles.
//=================================================================================================
template <int ndim>
void Simulation<ndim>::AllocateParticleMemory(void)
{
  int N;                               // Max. no. of stars/sinks

  debug2("[Simulation::AllocateParticleMemory]");

  // Allocate N-body memory (if using N-body)
  //-----------------------------------------------------------------------------------------------
  if (nbody) {

    // If sink particles are employed, allow enough memory for new sinks
    if (sink_particles == 1) {
      N = max(nbody->Nstar, 1024);
    }
    else N = nbody->Nstar;

    // Now call all memory allocation routines
    nbody->AllocateMemory(N);
    sinks->AllocateMemory(N);
  }
  //-----------------------------------------------------------------------------------------------

  // Allocate SPH memory, if being used
  if (hydro) hydro->AllocateMemory(hydro->Nhydro);


  return;
}



//=================================================================================================
//  Simulation::DeallocateParticleMemory
/// Deallocate all particle memory
//=================================================================================================
template <int ndim>
void Simulation<ndim>::DeallocateParticleMemory(void)
{
  debug2("[Simulation::DellocateParticleMemory]");

  sinks->DeallocateMemory();
  nbody->DeallocateMemory();
  hydro->DeallocateMemory();

  return;
}



//=================================================================================================
//  Simulation::PreSetupForPython
/// Initialisation routine called by python interface.
//=================================================================================================
template <int ndim>
void Simulation<ndim>::PreSetupForPython(void)
{
  debug1("[Simulation::PreSetupForPython]");

  // Check that IC type is really python
  if (simparams->stringparams["ic"] != "python") {
    string msg = "Error: you should call this function only if you are "
      "using \"python\" as \"ic\" parameter";
    ExceptionHandler::getIstance().raise(msg);
  }

  if (ParametersProcessed) {
    string msg = "Error: ProcessParameters has been already called!";
    ExceptionHandler::getIstance().raise(msg);
  }

  // Parse all parameters and set-up all objects required for simulation
  ProcessParameters();

  // Allocate all memory for both hydro and N-body particles
  hydro->Nhydro = simparams->intparams["Nhydro"];
  if (nbody)
    nbody->Nstar = simparams->intparams["Nstar"];
  AllocateParticleMemory();

  return;
}



//=================================================================================================
//  Simulation::ImportArrayNbody
/// Import an array containing nbody particle properties from python to C++ arrays.
//=================================================================================================
template <int ndim>
void Simulation<ndim>::ImportArrayNbody
 (double* input,                       ///< [in] Array of values from python
  int size,                            ///< [in] No. of array elements
  string quantity)                     ///< [in] String id of quantity being imported
{
  FLOAT StarParticle<ndim>::*quantityp = 0;            // Pointer to scalar quantity
  FLOAT (StarParticle<ndim>::*quantitypvec)[ndim] = 0; // Pointer to component of vector quantity
  int index = 0;                                       // Component index (if quantity is vector)
  bool scalar = false;                                 // Is the requested quantity a scalar?

  // Check that the size is correct
  if (size != nbody->Nstar) {
    stringstream message;
    message << "Error: the array you are passing has a size of " << size
            << ", but memory has been allocated for " << nbody->Nstar << " star particles";
    ExceptionHandler::getIstance().raise(message.str());
  }

  // Now set pointer to the correct value inside the particle data structure
  //-----------------------------------------------------------------------------------------------
  if (quantity == "x") {
    quantitypvec = &StarParticle<ndim>::r;
    index = 0;
    scalar = false;
  }
  //-----------------------------------------------------------------------------------------------
  else if (quantity == "y") {
    if (ndim < 2) {
      string message = "Error: loading y-coordinate array for ndim < 2";
      ExceptionHandler::getIstance().raise(message);
    }
    quantitypvec = &StarParticle<ndim>::r;
    index = 1;
    scalar = false;
  }
  //-----------------------------------------------------------------------------------------------
  else if (quantity == "z") {
    if (ndim < 3) {
      string message = "Error: loading y-coordinate array for ndim < 3";
      ExceptionHandler::getIstance().raise(message);
    }
    quantitypvec = &StarParticle<ndim>::r;
    index = 2;
    scalar = false;
  }
  //-----------------------------------------------------------------------------------------------
  else if (quantity == "vx") {
    quantitypvec = &StarParticle<ndim>::v;
    index = 0;
    scalar = false;
  }
  //-----------------------------------------------------------------------------------------------
  else if (quantity == "vy") {
    if (ndim < 2) {
      string message = "Error: loading vy array for ndim < 2";
      ExceptionHandler::getIstance().raise(message);
    }
    quantitypvec = &StarParticle<ndim>::v;
    index = 1;
    scalar = false;
  }
  //-----------------------------------------------------------------------------------------------
  else if (quantity == "vz") {
    if (ndim < 3) {
      string message = "Error: loading vz array for ndim < 3";
      ExceptionHandler::getIstance().raise(message);
    }
    quantitypvec = &StarParticle<ndim>::v;
    index = 2;
    scalar = false;
  }
  //-----------------------------------------------------------------------------------------------
  else if (quantity == "m") {
    quantityp = &StarParticle<ndim>::m;
    scalar = true;
  }
  //-----------------------------------------------------------------------------------------------
  else if (quantity == "h") {
    quantityp = &StarParticle<ndim>::h;
    scalar = true;
  }
  //-----------------------------------------------------------------------------------------------
  else {
    string message = "Quantity " + quantity + "not recognised";
    ExceptionHandler::getIstance().raise(message);
  }
  //-----------------------------------------------------------------------------------------------


  // Finally loop over particles and set all values
  // (Note that the syntax for scalar is different from the one for vectors)
  //-----------------------------------------------------------------------------------------------
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



//=================================================================================================
//  Simulation::ImportArraySph
/// Import an array containing sph particle properties from python to C++ arrays.
//=================================================================================================
template <int ndim>
void Simulation<ndim>::ImportArraySph
 (double* input,                       ///< [in] Array of values imported from python
  int size,                            ///< [in] No. of elements in array
  string quantity)                     ///< [in] String id of quantity
{
  FLOAT Particle<ndim>::*quantityp = 0;             // Pointer to scalar quantity
  FLOAT (Particle<ndim>::*quantitypvec)[ndim] = 0;  // Pointer to component of vector quantity
  int index = 0;                                    // If it's a component of a vector
                                                    // quantity, we need to know its index
  bool scalar = false;                              // Is the requested quantity a scalar?

  // Check that the size is correct
  if (size != hydro->Nhydro) {
    stringstream message;
    message << "Error: the array you are passing has a size of "
            << size << ", but memory has been allocated for " << hydro->Nhydro << " particles";
    ExceptionHandler::getIstance().raise(message.str());
  }


  // Now set pointer to the correct value inside the particle data structure
  //-----------------------------------------------------------------------------------------------
  if (quantity == "x") {
    quantitypvec = &Particle<ndim>::r;
    index = 0;
    scalar = false;
  }
  //-----------------------------------------------------------------------------------------------
  else if (quantity == "y") {
    if (ndim < 2) {
      string message = "Error: loading y-coordinate array for ndim < 2";
      ExceptionHandler::getIstance().raise(message);
    }
    quantitypvec = &Particle<ndim>::r;
    index = 1;
    scalar = false;
  }
  //-----------------------------------------------------------------------------------------------
  else if (quantity == "z") {
    if (ndim < 3) {
      string message = "Error: loading y-coordinate array for ndim < 3";
      ExceptionHandler::getIstance().raise(message);
    }
    quantitypvec = &Particle<ndim>::r;
    index = 2;
    scalar = false;
  }
  //-----------------------------------------------------------------------------------------------
  else if (quantity == "vx") {
    quantitypvec = &Particle<ndim>::v;
    index = 0;
    scalar = false;
  }
  //-----------------------------------------------------------------------------------------------
  else if (quantity == "vy") {
    if (ndim < 2) {
      string message = "Error: loading vy array for ndim < 2";
      ExceptionHandler::getIstance().raise(message);
    }
    quantitypvec = &Particle<ndim>::v;
    index = 1;
    scalar = false;
  }
  //-----------------------------------------------------------------------------------------------
  else if (quantity == "vz") {
    if (ndim < 3) {
      string message = "Error: loading vz array for ndim < 3";
      ExceptionHandler::getIstance().raise(message);
    }
    quantitypvec = &Particle<ndim>::v;
    index = 2;
    scalar = false;
  }
  //-----------------------------------------------------------------------------------------------
  else if (quantity == "rho") {
    //TODO: at the moment, if rho or h are uploaded, they will be just ignored.
    //Add some facility to use them
    quantityp = &Particle<ndim>::rho;
    scalar = true;
  }
  //-----------------------------------------------------------------------------------------------
  else if (quantity == "h") {
    quantityp = &Particle<ndim>::h;
    scalar = true;
  }
  //-----------------------------------------------------------------------------------------------
  else if (quantity == "u") {
    //TODO: add some facility for uploading either u, T, or cs, and compute automatically the
    //other ones depending on the EOS
    quantityp = &Particle<ndim>::u;
    scalar=true;
  }
  //-----------------------------------------------------------------------------------------------
  else if (quantity == "m") {
    quantityp = &Particle<ndim>::m;
    scalar = true;
  }
  //-----------------------------------------------------------------------------------------------
  else {
    string message = "Quantity " + quantity + "not recognised";
    ExceptionHandler::getIstance().raise(message);
  }
  //-----------------------------------------------------------------------------------------------


  // Finally loop over particles and set all values
  // (Note that the syntax for scalar is different from the one for vectors)
  //-----------------------------------------------------------------------------------------------
  if (scalar) {
    for (int i=0; i<size; i++) {
      Particle<ndim>& part = hydro->GetParticlePointer(i);
      part.*quantityp = input[i];
    }
  }
  else {
    for (int i=0; i<size; i++) {
      Particle<ndim>& part = hydro->GetParticlePointer(i);
      (part.*quantitypvec)[index] = input[i];
    }
  }

  return;
}



//=================================================================================================
//  Simulation::ImportArray
/// Import an array containing particle properties from python to C++ arrays.
/// This is a wrapper around ImportArraySph and ImportArrayNbody
//=================================================================================================
template <int ndim>
void Simulation<ndim>::ImportArray
 (double* input,                       ///< [in] Input array
  int size,                            ///< [in] Size of the input array
  string quantity,                     ///< [in] Quantity to be set equal to the given array
  string type)                         ///< [in] Particle type that should be assigned the array
{
  debug2("[Simulation::ImportArray]");

  // Check that PreSetup has been called
  if (! ParametersProcessed) {
    string msg = "Error: before calling ImportArray, you need to call PreSetupForPython!";
    ExceptionHandler::getIstance().raise(msg);
  }

  // Call the right function depending on the passed in type
  if (type == "sph") {
    // Check sph has been allocated
    if (hydro == NULL) {
      string message = "Error: memory for sph was not allocated! Are you sure that this is not a nbody-only simulation?";
      ExceptionHandler::getIstance().raise(message);
    }
    ImportArraySph(input, size, quantity);

  }
  else if (type == "star") {
    if (nbody == NULL) {
      string message = "Error: memory for nbody was not allocated! Are you sure that this is not a sph-only simulation?";
      ExceptionHandler::getIstance().raise(message);
    }
    ImportArrayNbody(input, size, quantity);
  }
  else {
    string message = "Error: we did not recognize the type " + type +
      ", the only allowed types are \"sph\" and \"nbody\"";
    ExceptionHandler::getIstance().raise(message);
  }

  return;
}



//=================================================================================================
//  Simulation::SetComFrame
/// Move all particles (both hydro and N-body) to centre-of-mass frame.
//=================================================================================================
template<int ndim>
void Simulation<ndim>::SetComFrame(void)
{
  int i;                            // Particle counter
  int k;                            // Dimension counter

  debug2("[Simulation::SetComFrame]");

  CalculateDiagnostics();

  for (i=0; i<hydro->Nhydro; i++) {
    Particle<ndim>& part = hydro->GetParticlePointer(i);
    for (k=0; k<ndim; k++) part.r[k] -= diag.rcom[k];
    for (k=0; k<ndim; k++) part.v[k] -= diag.vcom[k];
  }

  for (i=0; i<nbody->Nstar; i++) {
    for (k=0; k<ndim; k++) nbody->stardata[i].r[k] -= diag.rcom[k];
    for (k=0; k<ndim; k++) nbody->stardata[i].v[k] -= diag.vcom[k];
  }

  CalculateDiagnostics();

  return;
}



//=================================================================================================
//  Simulation::UpdateDiagnostics
/// Update energy error value after computing diagnostic quantities.
//=================================================================================================
template <int ndim>
void Simulation<ndim>::UpdateDiagnostics(void)
{
  if (rank == 0) {
    diag.Eerror = fabs(diag0.Etot - diag.Etot)/fabs(diag0.Etot);
    cout << "Eerror : " << diag.Eerror << endl;
  }
}



template class Simulation<1>;
template class Simulation<2>;
template class Simulation<3>;
