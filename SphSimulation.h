// ============================================================================
// SphSimulation.h
// ============================================================================


#ifndef _SPH_SIMULATION_H_
#define _SPH_SIMULATION_H_


#include <map>
#include <string>
#include "Parameters.h"
#include "SimUnits.h"
#include "SphKernel.h"
#include "Sph.h"
#include "SphNeighbourSearch.h"
#include "SphIntegration.h"
#include "EnergyEquation.h"
using namespace std;


struct DomainBox {
  string x_boundary_lhs;
  string x_boundary_rhs;
  string y_boundary_lhs;
  string y_boundary_rhs;
  string z_boundary_lhs;
  string z_boundary_rhs;
  float boxmin[3];
  float boxmax[3];
  float boxsize[3];
  float boxhalf[3];
  //float rmin[ndimmax];
  //float rmax[ndimmax];
};


struct Diagnostics {
  double Eerror;
  double Etot;
  double utot;
  double ketot;
  double gpetot;
  double mom[ndimmax];
  double angmom[3];
  double force[ndimmax];
  double force_grav[ndimmax];
};



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

  // Subroutine prototypes
  // --------------------------------------------------------------------------
  void Setup(void);
  void MainLoop(void);
  void Run(int=-1);
  void Output(void);
  void GenerateIC(void);
  void ComputeBlockTimesteps(void);
  void ProcessParameters(void);
  void CalculateDiagnostics(void);

  void SearchGhostParticles(void);
  void CreateGhostParticle(int,int,float,float);
  void CopyDataToGhosts(void);
  void CheckBoundaries(void);

  // Initial conditions routines
  // --------------------------------------------------------------------------
  void RandomBox(void);
  void RandomSphere(void);
  void ShockTube(void);
  void KHI(void);

  // Initial conditions helper routines
  // --------------------------------------------------------------------------
  void AddRandomBox(int, float *, DomainBox);
  void AddRandomSphere(int, float *, float *, float);
  void AddRegularLattice(int, int *, float *, DomainBox);

  // Input-output routines
  // --------------------------------------------------------------------------
  bool ReadSnapshotFile(string,string);
  bool ReadColumnSnapshotFile(string);

  bool WriteSnapshotFile(string,string);
  bool WriteColumnSnapshotFile(string);


#if !defined(FIXED_DIMENSIONS)
  int ndim;
  int vdim;
  int bdim;
#endif

  // Integer and physical Timestep counters
  // --------------------------------------------------------------------------
  int n;                                    // Integer time counter
  int Nsteps;                               // Total no. of steps in simulation
  int Nstepsmax;                            // Max. allowed no. of steps
  int noutputstep;                          // ..
  double t;                                 // Current simulation time
  double timestep;                          // Current timestep
  double tsnapnext;                         // Time of next snapshot
  double tend;                              // End time of simulation
  double dt_snap;                           // Snapshot time interval
  int Noutsnap;                             // No. of output snapshots

  string run_id;                            // Simulation id string
  string paramfile;                         // Name of parameters file
  bool setup;                               // Flags whether the simulation has been setup

  Parameters simparams;                     // Simulation parameters object
  SimUnits simunits;                        // Simulation units object
  DomainBox simbox;                         // Simulation boundary data
  Diagnostics diag0;                        // Initial diagnostic state
  Diagnostics diag;                         // Current diagnostic state

  Sph *sph;                                 // SPH algorithm pointer
  SphNeighbourSearch *sphneib;              // SPH Neighbour scheme pointer
  SphIntegration *sphint;                   // SPH Integration scheme pointer
  EnergyEquation *uint;                     // Energy equation pointer

};


#endif
