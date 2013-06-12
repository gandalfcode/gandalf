//=============================================================================
//  SphSimulation.h
//  Contains definitions for following data structures and classes:
//  - DomainBox
//  - Diagnostics
//  - SphSimulationBase
//  - SphSimulation
//=============================================================================


#ifndef _SPH_SIMULATION_H_
#define _SPH_SIMULATION_H_


#include <map>
#include <string>
#include <list>
#include "Diagnostics.h"
#include "DomainBox.h"
#include "Precision.h"
#include "Parameters.h"
#include "SimUnits.h"
#include "SphKernel.h"
#include "Sph.h"
#include "SphNeighbourSearch.h"
#include "SphIntegration.h"
#include "EnergyEquation.h"
#include "Nbody.h"
using namespace std;


class SphSnapshotBase;

//=============================================================================
//  Class SimulationBase
/// \brief  Creates a simulation object depending on the dimensionality.
/// \author D. A. Hubber, G. Rosotti
/// \date   03/04/2013
//=============================================================================
class SimulationBase
{
  // Stuff only for internal use of the class
  virtual void CalculateDiagnostics(void)=0;
  virtual void OutputDiagnostics(void)=0;
  virtual void UpdateDiagnostics(void)=0;
  virtual void GenerateIC(void)=0;
  virtual bool ReadColumnSnapshotFile(string)=0;
  virtual bool WriteColumnSnapshotFile(string)=0;


 public:

  static SimulationBase* SimulationFactory(int ndim, Parameters* params);

  // Constructor and Destructor
  // --------------------------------------------------------------------------
  SimulationBase(Parameters* params);
  ~SimulationBase();
  
  // Subroutine prototypes
  // --------------------------------------------------------------------------
  string GetParam(string key);
  void SetParam (string key, string value);
  void SetParam (string key, int value);
  void SetParam (string ket, float value);
  virtual void PreSetupForPython(void)=0;
  virtual void ImportArray(double* input, int size, string quantity, string type="sph")=0;
  void SetupSimulation(void);
  virtual void PostGeneration(void)=0;
  virtual void MainLoop(void)=0;
  void Run(int=-1);
  list<SphSnapshotBase*> InteractiveRun(int=-1);
  string Output(void);
  virtual void ProcessParameters(void)=0;

  // Input-output routines
  // --------------------------------------------------------------------------
  bool ReadSnapshotFile(string,string);
  bool WriteSnapshotFile(string,string);

  // Variables
  // --------------------------------------------------------------------------
  bool setup;                       ///< Flag if simulation is setup
  bool ParametersProcessed;         ///< Flag if params are already processed
  int integration_step;             ///< Steps per complete integration step
  int level_max;                    ///< Maximum timestep level
  int level_step;                   ///< Level of smallest timestep unit
  int n;                            ///< Integer time counter
  int nbody_single_timestep;        ///< Flag if stars use same timestep
  int ndims;                        ///< Aux. dimensionality variable. 
                                    ///< Required for python routines.
  int noutputstep;                  ///< Output frequency
  int nresync;                      ///< Integer time for resynchronisation
  int Nsteps;                       ///< Total no. of steps in simulation
  int Nstepsmax;                    ///< Max. allowed no. of steps
  int Nlevels;                      ///< No. of timestep levels
  int Noutsnap;                     ///< No. of output snapshots
  int sph_single_timestep;          ///< Flag if SPH ptcls use same step
  DOUBLE dt_max;                    ///< Value of maximum timestep level
  DOUBLE dt_snap;                   ///< Snapshot time interval
  DOUBLE t;                         ///< Current simulation time
  DOUBLE tend;                      ///< End time of simulation
  DOUBLE timestep;                  ///< Current timestep
  DOUBLE tsnapnext;                 ///< Time of next snapshot
  string out_file_form;             ///< Output snapshot file format
  string paramfile;                 ///< Name of parameters file
  string run_id;                    ///< Simulation id string

  Parameters* simparams;            ///< Simulation parameters object (pointer)
  SimUnits simunits;                ///< Simulation units object

};



#if !defined(SWIG)
//=============================================================================
//  Class SphSimulation
/// \brief  Main Sph Simulation class.
/// \author D. A. Hubber, G. Rosotti
/// \date   03/04/2013
//=============================================================================
template <int ndim>
class Simulation : public SimulationBase {

  void ImportArraySph(double* input, int size, string quantity);
  void ImportArrayNbody(double* input, int size, string quantity);

 public:
  Simulation(Parameters* parameters) : 
    SimulationBase(parameters),
    nbody(NULL),
    sph(NULL) {this->ndims=ndim;};


  // Initial conditions helper routines
  // --------------------------------------------------------------------------
  void AddAzimuthalDensityPerturbation(int, int, FLOAT, FLOAT *, FLOAT *); 
  void AddRotationalVelocityField(int, FLOAT, FLOAT *, FLOAT *, FLOAT *); 
  void AddRandomBox(int, FLOAT *, DomainBox<ndim>);
  void AddRandomSphere(int, FLOAT *, FLOAT *, FLOAT);
  void AddCubicLattice(int, int *, FLOAT *, DomainBox<ndim>, bool);
  void AddHexagonalLattice(int, int *, FLOAT *, DomainBox<ndim>, bool);
  int AddLatticeSphere(int, FLOAT *, FLOAT *, FLOAT, string);
  int CutSphere(int, int, FLOAT, FLOAT *, DomainBox<ndim>, bool);

  // Subroutine prototypes
  // --------------------------------------------------------------------------
  virtual void PreSetupForPython(void);
  virtual void ImportArray(double* input, int size, string quantity, string type="sph");
  virtual void GenerateIC(void);
  virtual void ProcessParameters(void)=0;
  virtual void CalculateDiagnostics(void);
  virtual void OutputDiagnostics(void);
  virtual void UpdateDiagnostics(void);
  virtual void ComputeGlobalTimestep(void)=0;
  virtual void ComputeBlockTimesteps(void)=0;

#if defined(VERIFY_ALL)
  void VerifyBlockTimesteps(void);
#endif

  // Ghost particle functions -> maybe move to Sph?
  // --------------------------------------------------------------------------
  void SearchGhostParticles(void);
  void CreateGhostParticle(int,int,FLOAT,FLOAT,FLOAT);
  void CopySphDataToGhosts(void);
  void CopyAccelerationFromGhosts(void);
  void CheckBoundaries(void); //should probably go inside SphIntegration

  // Initial conditions routines -> move either to Sph, either to Nbody
  // --------------------------------------------------------------------------
  void BinaryStar(void);
  void BossBodenheimer(void);
  void CheckInitialConditions(void);
  void ContactDiscontinuity(void);
  void KHI(void);
  void NohProblem(void);
  void PlummerSphere(void);
  void ShockTube(void);
  void SedovBlastWave(void);
  void ShearFlow(void);
  void SoundWave(void);
  void UniformBox(void);
  void UniformSphere(void);

  // Input-output routines
  // --------------------------------------------------------------------------
  virtual bool ReadColumnSnapshotFile(string);
  virtual bool WriteColumnSnapshotFile(string);

  // Variables
  // --------------------------------------------------------------------------
  static const int vdim=ndim;
  static const FLOAT invndim=1.0/ndim;

  DomainBox<ndim> simbox;               ///< Simulation boundary data
  Diagnostics<ndim> diag0;              ///< Initial diagnostic state
  Diagnostics<ndim> diag;               ///< Current diagnostic state
  EnergyEquation<ndim> *uint;           ///< Energy equation pointer
  Sph<ndim> *sph;                       ///< SPH algorithm pointer
  SphIntegration<ndim> *sphint;         ///< SPH Integration scheme pointer
  SphNeighbourSearch<ndim> *sphneib;    ///< SPH Neighbour scheme pointer
  Nbody<ndim> *nbody;                   ///< N-body algorithm pointer

};



template <int ndim>
class SphSimulation : public Simulation<ndim> {
  using SimulationBase::simparams;
  using Simulation<ndim>::sph;
  using Simulation<ndim>::nbody;
  using Simulation<ndim>::sphint;
  using Simulation<ndim>::uint;
  using Simulation<ndim>::sphneib;
  using Simulation<ndim>::n;
  using Simulation<ndim>::Nlevels;
  using Simulation<ndim>::Nsteps;
  using Simulation<ndim>::t;
  using Simulation<ndim>::timestep;
  using Simulation<ndim>::level_step;
  using Simulation<ndim>::Noutsnap;
  using Simulation<ndim>::tsnapnext;
  using Simulation<ndim>::dt_snap;
  using Simulation<ndim>::level_max;
  using Simulation<ndim>::integration_step;
  using Simulation<ndim>::nresync;
  using Simulation<ndim>::dt_max;
  using Simulation<ndim>::sph_single_timestep;
public:
  SphSimulation (Parameters* parameters): Simulation<ndim>(parameters) {};
  virtual void PostGeneration(void);
  virtual void MainLoop(void);
  virtual void ComputeGlobalTimestep(void);
  virtual void ComputeBlockTimesteps(void);
  virtual void ProcessParameters(void);
};



template <int ndim>
class GodunovSimulation : public Simulation<ndim> {
  using SimulationBase::simparams;
  using Simulation<ndim>::sph;
  using Simulation<ndim>::nbody;
  using Simulation<ndim>::sphint;
  using Simulation<ndim>::uint;
  using Simulation<ndim>::sphneib;
  using Simulation<ndim>::n;
  using Simulation<ndim>::Nlevels;
  using Simulation<ndim>::Nsteps;
  using Simulation<ndim>::t;
  using Simulation<ndim>::timestep;
  using Simulation<ndim>::level_step;
  using Simulation<ndim>::Noutsnap;
  using Simulation<ndim>::tsnapnext;
  using Simulation<ndim>::dt_snap;
  using Simulation<ndim>::level_max;
  using Simulation<ndim>::integration_step;
  using Simulation<ndim>::nresync;
  using Simulation<ndim>::dt_max;
  using Simulation<ndim>::sph_single_timestep;
public:
  GodunovSimulation (Parameters* parameters): Simulation<ndim>(parameters) {};
  virtual void PostGeneration(void);
  virtual void MainLoop(void);
  virtual void ComputeGlobalTimestep(void);
  virtual void ComputeBlockTimesteps(void);
  virtual void ProcessParameters(void);
};

#endif


#endif
