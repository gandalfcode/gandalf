#include "Exception.h"
#include <iostream>
#include "Simulation.h"
using namespace std;

int main(int argc, char** argv)
{
  SimulationBase* sim;
  Parameters* params = new Parameters();
  string paramfile;
  ExceptionHandler::makeExceptionHandler(cplusplus);

  if (argc >= 2){
    paramfile = argv[1];
  }
  else {
    cout << "No parameter file specified, aborting..." << endl;
    exit(-1);
  }

  params->ReadParamsFile(paramfile);

  sim = SimulationBase::SimulationFactory(params->intparams["ndim"], params);

  sim->SetupSimulation();

  sim->Run();

  return 0;
}
