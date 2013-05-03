//=============================================================================
// SphSimulation.h
// Contains definitions for following data structures and classes:
// - DomainBox
// - Diagnostics
// - SphSimulationBase
// - SphSimulation
//=============================================================================


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
#include "Nbody.h"
using namespace std;


//=============================================================================
//  Structure DomainBox
/// \brief  Bounding box data structure.
/// \author D. A. Hubber, G. Rosotti
/// \date   03/04/2013
//=============================================================================
template <int ndim>
struct DomainBox {
  string x_boundary_lhs;                ///< x-dimension LHS boundary condition
  string x_boundary_rhs;                ///< x-dimension RHS boundary condition
  string y_boundary_lhs;                ///< y-dimension LHS boundary condition
  string y_boundary_rhs;                ///< y-dimension RHS boundary condition
  string z_boundary_lhs;                ///< z-dimension LHS boundary condition
  string z_boundary_rhs;                ///< z-dimension RHS boundary condition
  FLOAT boxmin[ndim];                   ///< Minimum bounding box extent
  FLOAT boxmax[ndim];                   ///< Maximum bounding box extent
  FLOAT boxsize[ndim];                  ///< Side-lengths of bounding box
  FLOAT boxhalf[ndim];                  ///< Half side-lengths of bounding box
};



//=============================================================================
//  Structure Diagnostics
/// \brief  Structure containing snapshot of current diagnostic quantities.
/// \author D. A. Hubber, G. Rosotti
/// \date   03/04/2013
//=============================================================================
template <int ndim>
struct Diagnostics {
  DOUBLE Eerror;                        ///< Total energy error
  DOUBLE Etot;                          ///< Total energy
  DOUBLE utot;                          ///< Total thermal energy
  DOUBLE ketot;                         ///< Total kinetic energy
  DOUBLE gpetot;                        ///< Total grav. potential energy
  DOUBLE mom[ndim];                     ///< Total momentum vector
  DOUBLE angmom[3];                     ///< Total angular momentum vector
  DOUBLE force[ndim];                   ///< Net force
  DOUBLE force_grav[ndim];              ///< Net gravitational force
};



//=============================================================================
//  Class SphSimulationBase
/// \brief  Creates a simulation object depending on the dimensionality.
/// \author D. A. Hubber, G. Rosotti
/// \date   03/04/2013
// ============================================================================
class SphSimulationBase
{
 public:

  static SphSimulationBase* SphSimulationFactory(int ndim, Parameters* params);

  // Constructor and Destructor
  // --------------------------------------------------------------------------
  SphSimulationBase(Parameters* params);
  ~SphSimulationBase();

  // Subroutine prototypes
  // --------------------------------------------------------------------------
  void SetParam (string key, string value);
  void SetParam (string key, int value);
  void SetParam (string ket, float value);
  virtual void PreSetupForPython(void)=0;
  virtual void ImportArray(double* input, int size, string quantity)=0;
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
  virtual void ContactDiscontinuity(void)=0;
  virtual void LatticeBox(void)=0;
  virtual void RandomSphere(void)=0;
  virtual void ShearFlow(void) = 0;
  virtual void ShockTube(void)=0;
  virtual void KHI(void)=0;
  virtual void SoundWave(void)=0;

  // Input-output routines
  // --------------------------------------------------------------------------
  virtual bool ReadSnapshotFile(string,string)=0;
  virtual bool ReadColumnSnapshotFile(string)=0;
  virtual bool WriteSnapshotFile(string,string)=0;
  virtual bool WriteColumnSnapshotFile(string)=0;

  // Variables
  // --------------------------------------------------------------------------
  bool setup;                       ///< Flag if simulation is setup
  bool ParametersProcessed;         ///< Flag if the parameters have been already processed
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



//=============================================================================
//  Class SphSimulation
/// \brief  Main Sph Simulation class.
/// \author D. A. Hubber, G. Rosotti
/// \date   03/04/2013
//=============================================================================
template <int ndim>
class SphSimulation : public SphSimulationBase {
 public:
  SphSimulation(Parameters* parameters) : 
    SphSimulationBase(parameters) {this->ndims=ndim;};


  // Initial conditions helper routines
  // --------------------------------------------------------------------------
  void AddRandomBox(int, FLOAT *, DomainBox<ndim>);
  void AddRandomSphere(int, FLOAT *, FLOAT *, FLOAT);
  void AddRegularLattice(int, int *, FLOAT *, DomainBox<ndim>);
  void AddHexagonalLattice(int, int *, FLOAT *, DomainBox<ndim>);
  int CutSphere(int, int, FLOAT, FLOAT *, DomainBox<ndim>, bool);


  // Subroutine prototypes
  // --------------------------------------------------------------------------
  void PreSetupForPython(void);
  void ImportArray(double* input, int size, string quantity);
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
  virtual void ContactDiscontinuity(void);
  virtual void LatticeBox(void);
  virtual void RandomSphere(void);
  virtual void ShockTube(void);
  virtual void KHI(void);
  virtual void SedovBlastWave(void);
  virtual void ShearFlow(void);
  virtual void SoundWave(void);

  // Input-output routines
  // --------------------------------------------------------------------------
  virtual bool ReadSnapshotFile(string,string);
  virtual bool ReadColumnSnapshotFile(string);
  virtual bool WriteSnapshotFile(string,string);
  virtual bool WriteColumnSnapshotFile(string);

  // Variables
  // --------------------------------------------------------------------------
  static const int vdim=ndim;

  DomainBox<ndim> simbox;               ///< Simulation boundary data
  Diagnostics<ndim> diag0;              ///< Initial diagnostic state
  Diagnostics<ndim> diag;               ///< Current diagnostic state
  EnergyEquation<ndim> *uint;           ///< Energy equation pointer
  Sph<ndim> *sph;                       ///< SPH algorithm pointer
  SphIntegration<ndim> *sphint;         ///< SPH Integration scheme pointer
  SphNeighbourSearch<ndim> *sphneib;    ///< SPH Neighbour scheme pointer
  Nbody<ndim> *nbody;                   ///< N-body algorithm pointer

};

#endif
