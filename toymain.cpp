#include <iostream>
#include "SphSimulation.h"
using namespace std;

int main(int argc, char** argv)
{
  SphSimulation sim;

  if (argc >= 2){
    sim.paramfile = argv[1];
  }
  else {
    cout << "No parameter file specified, aborting..." << endl;
    exit(-1);
  }

  sim.Setup();

  sim.Run();

  return 0;
}
