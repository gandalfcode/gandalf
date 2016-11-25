//=================================================================================================
//  GANDALF :
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


#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include "Exception.h"
#include "Parameters.h"
#include "Simulation.h"
#include "CodeTiming.h"
#ifdef MPI_PARALLEL
#include "mpi.h"
#endif
using namespace std;



//=================================================================================================
//  main
/// Main GANDALF program for running executables from the command line.
//=================================================================================================
int main(int argc, char** argv)
{
  bool restart = false;                              // Flag restart simulation
  int rank = 0;                                      // Local copy of MPI rank
  string paramfile;                                  // Name of parameters file
  stringstream ss;                                   // Stream char to string
  ofstream outfile;                                  // Stream of temp restart file
  Parameters* params = new Parameters();             // Parameters object
  SimulationBase* sim;                               // Main simulation object
  ExceptionHandler::makeExceptionHandler(cplusplus); // Exception handler

#ifdef MPI_PARALLEL
  // Initialise all MPI processes (if activated in Makefile)
  int mpi_thread_support;
  int n_mpi_cpus;
  int required_mpi_thread_support = MPI_THREAD_SINGLE;
#ifdef _OPENMP
  required_mpi_thread_support = MPI_THREAD_FUNNELED;
#endif
  MPI_Init_thread(&argc, &argv, required_mpi_thread_support, &mpi_thread_support);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &n_mpi_cpus);

  // Determines if we have been spawned from python
  MPI_Comm parent;
  MPI_Comm_get_parent(&parent);

  // Tell exception handler to call MPI_Abort on error
  ExceptionHandler::set_mpi(1);

  if (!isPowerOfTwo(n_mpi_cpus)) {
    string message = "Error: currently in GANDALF the number of processes "
        "needs to be a power of two!";
    ExceptionHandler::getIstance().raise(message);
  }

#ifdef _OPENMP
  // Check that OpenMP and MPI can work together
  if (mpi_thread_support == MPI_THREAD_SINGLE)
#ifdef OPEN_MPI
    if (rank==0) {
      cout << "Warning: we detected that you are running with OpenMP and MPI" << endl;
      cout << "Your MPI implementation is OpenMPI, and it reported that it cannot be used "
          "with OpenMP. However we know that often OpenMPI does that even when there is "
          "no problem. We will go ahead, but please check carefully the result" << endl;
    }
#else
    ExceptionHandler::getIstance().raise("This implementation of MPI is not interoperable with OpenMP, aborting!"
        "Refer to your system administrator to know how to solve this problem");
#endif
  else {
    string message;
    if (mpi_thread_support == MPI_THREAD_FUNNELED) {
      message = "MPI_THREAD_FUNNELED";
    }
    else if (mpi_thread_support == MPI_THREAD_SERIALIZED) {
      message = "MPI_THREAD_SERIALIZED";
    }
    else if (mpi_thread_support == MPI_THREAD_MULTIPLE) {
      message = "MPI_THREAD_MULTIPLE";
    }
    cout << "The level of MPI thread support is " << message << endl;
  }
#endif
#endif


  // Parse and process all command-line arguments.
  if (argc == 3) {
    if (string(argv[1]) == "-r") {
      restart = true;
      ss << argv[2];
      ss >> paramfile;
    }
  }
  else if (argc == 2) {
    paramfile = string(argv[1]);
  }
  else {
    string message = "No parameter file specified, aborting...";
    ExceptionHandler::getIstance().raise(message);
  }

  // Create empty file (used for automatic restarts on clusters)
  outfile.open("cont");
  outfile.close();

  // Read parameters file immediately and record to file
  params->ReadParamsFile(paramfile);
#ifdef MPI_PARALLEL
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  if (rank == 0) params->RecordParametersToFile();

  // Create simulation object with required dimensionality and parameters
  sim = SimulationBase::SimulationFactory(params->intparams["ndim"],
                                          params->stringparams["sim"], params);
  sim->restart = restart;

  // Print out splash screen
  if (rank == 0) sim->SplashScreen(paramfile);

#if defined MPI_PARALLEL
  if (rank == 0) {
  cout << "Running with MPI, using " << n_mpi_cpus << " tasks" << endl;
  }
#endif
#if defined _OPENMP
  if (rank == 0) {
  cout << "Running with OPENMP, using " << omp_get_max_threads() << " threads" << endl;
  }
#if defined MPI_PARALLEL
  if (rank == 0) {
  cout << "Hybrid OpenMP/MPI parallelization currently in use, for a total of "
       << n_mpi_cpus*omp_get_max_threads() << " cores" << endl;
  }
#endif
#endif

  // Perform all set-up procedures
  sim->SetupSimulation();

  // Run entire simulation until specified end-time in parameters file.
  sim->Run();

  // Compile timing statistics from complete simulation
  sim->timing->ComputeTimingStatistics(sim->run_id);

#ifdef MPI_PARALLEL
  if (parent != MPI_COMM_NULL){
#if MPI_VERSION>=3
    MPI_Request req;
    MPI_Ibarrier(parent,&req);
    MPI_Wait(&req,MPI_STATUS_IGNORE);
#else
     MPI_Barrier(parent);
#endif
  }
  MPI_Finalize();
#endif

  // Compile timing statistics from complete simulation
  //timing.ComputeTimingStatistics(sim->run_id);

  // Finally, delete all locally created objects
  delete sim;
  delete params;


  return 0;
}
