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

template <int ndim>
struct DomainBox {
  string x_boundary_lhs;
  string x_boundary_rhs;
  string y_boundary_lhs;
  string y_boundary_rhs;
  string z_boundary_lhs;
  string z_boundary_rhs;
  FLOAT boxmin[ndim];
  FLOAT boxmax[ndim];
  FLOAT boxsize[ndim];
  FLOAT boxhalf[ndim];
  //FLOAT rmin[ndimmax];
  //FLOAT rmax[ndimmax];
};

template <int ndim>
struct Diagnostics {
  DOUBLE Eerror;
  DOUBLE Etot;
  DOUBLE utot;
  DOUBLE ketot;
  DOUBLE gpetot;
  DOUBLE mom[ndim];
  //TODO: I assume that the 3 here is correct...
  DOUBLE angmom[3];
  DOUBLE force[ndim];
  DOUBLE force_grav[ndim];
};



// ============================================================================
// Clas SphSimulation
// ============================================================================
class SphSimulationBase
{
 public:

  static SphSimulationBase* SphSimulationFactory (int ndim, Parameters* params);

  // Constructor and Destructor
  // --------------------------------------------------------------------------
  SphSimulationBase(Parameters* params);
  ~SphSimulationBase();

  // Subroutine prototypes
  // --------------------------------------------------------------------------
  virtual void Setup(void)=0;
  virtual void PreSetupForPython(void)=0;
  virtual void ImportArray(double* input, int size, string quantity)=0;
  virtual void PostSetupForPython(void)=0;
  virtual void SetupSimulation(void)=0;
  virtual void PostGeneration(void)=0;
  virtual void MainLoop(void)=0;
  virtual void Run(int=-1)=0;
  virtual void InteractiveRun(int=-1)=0;
  virtual void Output(void)=0;
  virtual void GenerateIC(void)=0;
  virtual void ProcessParameters(void)=0;
  virtual void CalculateDiagnostics(void)=0;

  virtual void ComputeGlobalTimestep(void)=0;
  virtual void ComputeBlockTimesteps(void)=0;
#if defined(VERIFY_ALL)
  virtual void VerifyBlockTimesteps(void)=0;
#endif

  // Ghost particle functions
  // --------------------------------------------------------------------------
  virtual void SearchGhostParticles(void)=0;
  virtual void CreateGhostParticle(int,int,FLOAT,FLOAT,FLOAT)=0;
  virtual void CopySphDataToGhosts(void)=0;
  virtual void CopyAccelerationFromGhosts(void)=0;
  virtual void CheckBoundaries(void)=0;

  // Initial conditions routines
  // --------------------------------------------------------------------------
  virtual void RandomBox(void)=0;
  virtual void LatticeBox(void)=0;
  virtual void RandomSphere(void)=0;
  virtual void ShockTube(void)=0;
  virtual void KHI(void)=0;
  virtual void SoundWave(void)=0;



  // Input-output routines
  // --------------------------------------------------------------------------
  virtual bool ReadSnapshotFile(string,string)=0;
  virtual bool ReadColumnSnapshotFile(string)=0;
  virtual bool WriteSnapshotFile(string,string)=0;
  virtual bool WriteColumnSnapshotFile(string)=0;

//#if !defined(FIXED_DIMENSIONS)
  int ndims;
//  int vdim;
//  int bdim;
//#endif


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

  Parameters* simparams;                     // Simulation parameters object
  SimUnits simunits;                        // Simulation units object

};


template <int ndim>
class SphSimulation : public SphSimulationBase {
public:
  SphSimulation(Parameters* parameters) : SphSimulationBase(parameters) {this->ndims=ndim;};


  // Initial conditions helper routines
  // --------------------------------------------------------------------------
  void AddRandomBox(int, FLOAT *, DomainBox<ndim>);
  void AddRandomSphere(int, FLOAT *, FLOAT *, FLOAT);
  void AddRegularLattice(int, int *, FLOAT *, DomainBox<ndim>);
  void AddHexagonalLattice(int, int *, FLOAT *, DomainBox<ndim>);
  int CutSphere(int, int, FLOAT, FLOAT *, DomainBox<ndim>, bool);


  // Subroutine prototypes
  // --------------------------------------------------------------------------
  void Setup(void);
  void PreSetupForPython(void);
  void ImportArray(double* input, int size, string quantity);
  void PostSetupForPython(void);
  void SetupSimulation(void);
  void PostGeneration(void);
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

  // Ghost particle functions
  // --------------------------------------------------------------------------
  void SearchGhostParticles(void);
  void CreateGhostParticle(int,int,FLOAT,FLOAT,FLOAT);
  void CopySphDataToGhosts(void);
  void CopyAccelerationFromGhosts(void);
  virtual void CheckBoundaries(void);

  // Initial conditions routines
  // --------------------------------------------------------------------------
  virtual void RandomBox(void);
  virtual void LatticeBox(void);
  virtual void RandomSphere(void);
  virtual void ShockTube(void);
  virtual void KHI(void);
  virtual void SoundWave(void);



  // Input-output routines
  // --------------------------------------------------------------------------
  virtual bool ReadSnapshotFile(string,string);
  virtual bool ReadColumnSnapshotFile(string);
  virtual bool WriteSnapshotFile(string,string);
  virtual bool WriteColumnSnapshotFile(string);



  DomainBox<ndim> simbox;                         // Simulation boundary data
  Diagnostics<ndim> diag0;                        // Initial diagnostic state
  Diagnostics<ndim> diag;                         // Current diagnostic state
  SphNeighbourSearch<ndim> *sphneib;              // SPH Neighbour scheme pointer
  SphIntegration<ndim> *sphint;                   // SPH Integration scheme pointer
  EnergyEquation<ndim> *uint;                     // Energy equation pointer
  Sph<ndim> *sph;                                 // SPH algorithm pointer

  static const int vdim=ndim;

};

#endif
