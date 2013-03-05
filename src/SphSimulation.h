// ============================================================================
// SphSimulation.h
// ============================================================================


#ifndef _SPH_SIMULATION_H_
#define _SPH_SIMULATION_H_


#include <map>
#include <string>
#include "Precision.h"
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
  FLOAT boxmin[3];
  FLOAT boxmax[3];
  FLOAT boxsize[3];
  FLOAT boxhalf[3];
  //FLOAT rmin[ndimmax];
  //FLOAT rmax[ndimmax];
};


struct Diagnostics {
  DOUBLE Eerror;
  DOUBLE Etot;
  DOUBLE utot;
  DOUBLE ketot;
  DOUBLE gpetot;
  DOUBLE mom[ndimmax];
  DOUBLE angmom[3];
  DOUBLE force[ndimmax];
  DOUBLE force_grav[ndimmax];
};



// ============================================================================
// Clas SphSimulation
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
  void SetupSimulation(void);
  void MainLoop(void);
  void Run(int=-1);
  void InteractiveRun(int=-1);
  void Output(void);
  void GenerateIC(void);
  void ProcessParameters(void);
  void CalculateDiagnostics(void);

  void ComputeGlobalTimestep(void);
  void ComputeBlockTimesteps(void);
#if defined(VERIFY_ALL)
  void VerifyBlockTimesteps(void);
#endif

  void SearchGhostParticles(void);
  void CreateGhostParticle(int,int,FLOAT,FLOAT,FLOAT);
  void CopyAccelerationFromGhosts(void);
  void CheckBoundaries(void);

  // Initial conditions routines
  // --------------------------------------------------------------------------
  void RandomBox(void);
  void LatticeBox(void);
  void RandomSphere(void);
  void ShockTube(void);
  void KHI(void);
  void SoundWave(void);

  // Initial conditions helper routines
  // --------------------------------------------------------------------------
  void AddRandomBox(int, FLOAT *, DomainBox);
  void AddRandomSphere(int, FLOAT *, FLOAT *, FLOAT);
  void AddRegularLattice(int, int *, FLOAT *, DomainBox);
  void AddHexagonalLattice(int, int *, FLOAT *, DomainBox);
  int CutSphere(int, int, FLOAT, FLOAT *, DomainBox, bool);

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
  int sph_single_timestep;                  // ..
  int nbody_single_timestep;                // ..
  int integration_step;                     // ..
  int level_max;                            // Maximum timestep level
  int level_step;                           // ..
  int n;                                    // Integer time counter
  int Nsteps;                               // Total no. of steps in simulation
  int Nstepsmax;                            // Max. allowed no. of steps
  int Nlevels;                              // No. of timestep levels
  int noutputstep;                          // Output frequency
  int nresync;                              // ..
  DOUBLE dt_max;                            // ..
  DOUBLE t;                                 // Current simulation time
  DOUBLE timestep;                          // Current timestep
  DOUBLE tsnapnext;                         // Time of next snapshot
  DOUBLE tend;                              // End time of simulation
  DOUBLE dt_snap;                           // Snapshot time interval
  int Noutsnap;                             // No. of output snapshots

  string run_id;                            // Simulation id string
  string paramfile;                         // Name of parameters file
  string out_file_form;                     // ..
  bool setup;                               // Flag if simulation is setup

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
