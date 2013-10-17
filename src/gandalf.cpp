//=============================================================================
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
//=============================================================================


#include <iostream>
#include "Exception.h"
#include "Parameters.h"
#include "Simulation.h"
#ifdef MPI_PARALLEL
#include "mpi.h"
#endif
using namespace std;


//=============================================================================
//  main
/// Main GANDALF program for running executables from the command line.
//=============================================================================
int main(int argc, char** argv)
{
  SimulationBase* sim;                               // Main simulation object
  Parameters* params = new Parameters();             // Parameters object
  string paramfile;                                  // Name of parameters file
  ExceptionHandler::makeExceptionHandler(cplusplus); // Exception handler


  // Initialise all MPI processes (if activated in Makefile)
#ifdef MPI_PARALLEL
  MPI_Init(&argc,&argv);
#endif

  // Check that a valid number of arguments have been passed
  if (argc >= 2){
    paramfile = argv[1];
  }
  else {
    cout << "No parameter file specified, aborting..." << endl;
    exit(-1);
  }

  // Read parameters file immediately
  params->ReadParamsFile(paramfile);

  // Create simulation object with required dimensionality and parameters
  sim = SimulationBase::SimulationFactory(params->intparams["ndim"], params);

  // Perform all set-up procedures
  sim->SetupSimulation();

  // Run entire simulation until specified end-time in parameters file.
  sim->Run();

  return 0;
}
