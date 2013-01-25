// ============================================================================
// SphSimulation.H
// ============================================================================


#ifndef _SPH_SIMULATION_H_
#define _SPH_SIMULATION_H_


#include <map>
#include <string>
#include "Parameters.h"
#include "SimUnits.h"
#include "SphKernel.h"
#include "Sph.h"
#include "SphSnapshot.h"
#include "SphNeighbourSearch.h"
#include "SphIntegration.h"
using namespace std;


// ============================================================================
// CLASS SphSimulation
// ============================================================================
class SphSimulation
{
 public:

  // Constructor and Destructor
  // --------------------------------------------------------------------------
  SphSimulation();
  ~SphSimulation();

  void GenerateIC(int);
  void Setup(void);
  void MainLoop(void);
  void ComputeBlockTimesteps(void);

#if !defined(FIXED_DIMENSIONS)
  int ndim;
  int vdim;
  int bdim;
#endif

  int n;
  int Nsteps;
  double t;
  double timestep;
  
  string paramfile;

  Parameters simparams;
  SimUnits simunits;
  Sph *sph;
  SphNeighbourSearch *sphneib;
  SphIntegration *sphint;

};



#endif
