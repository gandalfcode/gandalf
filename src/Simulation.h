//=============================================================================
//  Simulation.h
//  Contains definitions for following data structures and classes:
//  - DomainBox
//  - Diagnostics
//  - SimulationBase
//  - Simulation
//  - SphSimulation
//  - GodunovSphSimulation
//
//  This file is part of GANDALF :
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


#ifndef _SIMULATION_H_
#define _SIMULATION_H_


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
#include "NbodySystemTree.h"
#include "Ghosts.h"
#include "Sinks.h"
#include "HeaderInfo.h"
using namespace std;
#ifdef MPI_PARALLEL
#include "mpi.h"
#include "MpiControl.h"
#endif


// Forward declaration of SphSnapshotBase to prevent circular dependency
class SphSnapshotBase;


//=============================================================================
//  Class SimulationBase
/// \brief  Creates a simulation object depending on the dimensionality.
/// \author D. A. Hubber, G. Rosotti
/// \date   03/04/2013
//=============================================================================
class SimulationBase
{
  // Subroutines only for internal use of the class
  virtual void CalculateDiagnostics(void)=0;
  virtual void OutputDiagnostics(void)=0;
  virtual void UpdateDiagnostics(void)=0;
  virtual void GenerateIC(void)=0;
  virtual void ReadColumnHeaderFile(ifstream& infile, HeaderInfo& info)=0;
  virtual bool ReadColumnSnapshotFile(string)=0;
  virtual bool WriteColumnSnapshotFile(string)=0;
  virtual void ReadSerenFormHeaderFile(ifstream& infile, HeaderInfo& info)=0;
  virtual bool ReadSerenFormSnapshotFile(string)=0;
  virtual bool WriteSerenFormSnapshotFile(string)=0;

  std::list<string> keys;

 public:

  static SimulationBase* SimulationFactory(int ndim, Parameters* params);

  // Constructor and Destructor
  //---------------------------------------------------------------------------
  SimulationBase(Parameters* params);
  ~SimulationBase();
  
  // Subroutine prototypes
  //---------------------------------------------------------------------------
  string GetParam(string key);
  string Output(void);
  void SetParam(string key, string value);
  void SetParam(string key, int value);
  void SetParam(string ket, float value);
  std::list<string>* GetIntAndFloatParameterKeys();
  void SetupSimulation(void);
  void SplashScreen(void);
  void Run(int=-1);
  list<SphSnapshotBase*> InteractiveRun(int=-1);

  virtual void ImportArray(double* input, int size, 
                           string quantity, string type="sph") = 0;
  virtual void MainLoop(void)=0;
  virtual void PostInitialConditionsSetup(void)=0;
  virtual void PreSetupForPython(void)=0;
  virtual void ProcessParameters(void)=0;
  virtual void SetComFrame(void)=0;


  // Input-output routines
  //---------------------------------------------------------------------------
  bool ReadSnapshotFile(string,string);
  bool WriteSnapshotFile(string,string);
  HeaderInfo ReadHeaderSnapshotFile(string filename, string format);

  // Variables
  //---------------------------------------------------------------------------
  bool setup;                       ///< Flag if simulation is setup
  bool ParametersProcessed;         ///< Flag if params are already processed
  int integration_step;             ///< Steps per complete integration step
  int level_diff_max;               ///< Max. allowed neib timestep level diff
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
  int rank;                         ///< Process i.d. (for MPI simulations)
  int sink_particles;               ///< Switch on sink particles
  int sph_single_timestep;          ///< Flag if SPH ptcls use same step
  DOUBLE dt_max;                    ///< Value of maximum timestep level
  DOUBLE dt_snap;                   ///< Snapshot time interval
  DOUBLE dt_python;                 ///< Python window update time interval
  DOUBLE t;                         ///< Current simulation time
  DOUBLE tend;                      ///< End time of simulation
  DOUBLE timestep;                  ///< Current timestep
  DOUBLE tsnapfirst;                ///< Time of first snapshot
  DOUBLE tsnapnext;                 ///< Time of next snapshot
  string out_file_form;             ///< Output snapshot file format
  string paramfile;                 ///< Name of parameters file
  string run_id;                    ///< Simulation id string

  Parameters* simparams;            ///< Simulation parameters object (pointer)
  SimUnits simunits;                ///< Simulation units object

};



#if !defined(SWIG)
//=============================================================================
//  Class Simulation
/// \brief   Main Simulation class.
/// \details Main parent Simulation class from which all other simulation 
///          objects (e.g. SphSimulation) inherit from.
/// \author  D. A. Hubber, G. Rosotti
/// \date    03/04/2013
//=============================================================================
template <int ndim>
class Simulation : public SimulationBase 
{
  void ImportArraySph(double* input, int size, string quantity);
  void ImportArrayNbody(double* input, int size, string quantity);

 public:
  Simulation(Parameters* parameters) : 
    SimulationBase(parameters),
    nbody(NULL),
    sph(NULL) {this->ndims=ndim;};


  // Memory allocation routines
  //---------------------------------------------------------------------------
  void AllocateParticleMemory(void);
  void DeallocateParticleMemory(void);


  // Initial conditions helper routines
  //---------------------------------------------------------------------------
  void AddBinaryStar(DOUBLE, DOUBLE, DOUBLE, DOUBLE, DOUBLE, DOUBLE, DOUBLE *,
                     DOUBLE *, NbodyParticle<ndim> &, NbodyParticle<ndim> &);
  void AddAzimuthalDensityPerturbation(int, int, FLOAT, FLOAT *, FLOAT *); 
  void AddRotationalVelocityField(int, FLOAT, FLOAT *, FLOAT *, FLOAT *); 
  void AddRandomBox(int, FLOAT *, DomainBox<ndim>);
  void AddRandomSphere(int, FLOAT *, FLOAT *, FLOAT);
  void AddCubicLattice(int, int *, FLOAT *, DomainBox<ndim>, bool);
  void AddHexagonalLattice(int, int *, FLOAT *, DomainBox<ndim>, bool);
  int AddLatticeSphere(int, FLOAT *, FLOAT *, FLOAT, string);
  int CutSphere(int, int, FLOAT, FLOAT *, DomainBox<ndim>, bool);

  // Subroutine prototypes
  //---------------------------------------------------------------------------
  virtual void CalculateDiagnostics(void);
  virtual void ComputeGlobalTimestep(void)=0;
  virtual void ComputeBlockTimesteps(void)=0;
  virtual void GenerateIC(void);
  virtual void ImportArray(double* input, int size, 
                           string quantity, string type="sph");
  virtual void PreSetupForPython(void);
  virtual void ProcessGodunovSphParameters(void);
  virtual void ProcessNbodyParameters(void);
  virtual void ProcessParameters(void);
  virtual void ProcessSphParameters(void);
  virtual void OutputDiagnostics(void);
  virtual void UpdateDiagnostics(void);
  virtual void SetComFrame(void);

#if defined(VERIFY_ALL)
  void VerifyBlockTimesteps(void);
#endif

  // Initial conditions routines -> move either to Sph, either to Nbody
  //---------------------------------------------------------------------------
  void BinaryAccretion(void);
  void BinaryStar(void);
  void BossBodenheimer(void);
  void CheckInitialConditions(void);
  void ContactDiscontinuity(void);
  void KHI(void);
  void NohProblem(void);
  void PlummerSphere(void);
  void QuadrupleStar(void);
  void ShockTube(void);
  void SedovBlastWave(void);
  void ShearFlow(void);
  void SoundWave(void);
  void TripleStar(void);
  void UniformBox(void);
  void UniformSphere(void);

  // Input-output routines
  //---------------------------------------------------------------------------
  virtual void ReadColumnHeaderFile(ifstream& infile, HeaderInfo& info);
  virtual bool ReadColumnSnapshotFile(string);
  virtual bool WriteColumnSnapshotFile(string);
  virtual void ReadSerenFormHeaderFile(ifstream& infile, HeaderInfo& info);
  virtual bool ReadSerenFormSnapshotFile(string);
  virtual bool WriteSerenFormSnapshotFile(string);

  // Variables
  //---------------------------------------------------------------------------
  static const int vdim=ndim;
  static const FLOAT invndim=1.0/ndim;

  DomainBox<ndim> simbox;               ///< Simulation boundary data
  Diagnostics<ndim> diag0;              ///< Initial diagnostic state
  Diagnostics<ndim> diag;               ///< Current diagnostic state
  EnergyEquation<ndim> *uint;           ///< Energy equation pointer
  Ghosts<ndim> ghosts;                  ///< Ghost particle object
  Nbody<ndim> *nbody;                   ///< N-body algorithm pointer
  Nbody<ndim> *subsystem;               ///< N-body object for sub-systems
  NbodySystemTree<ndim> nbodytree;      ///< N-body tree to create sub-systems
  Sinks<ndim> sinks;                    ///< Sink particle object
  Sph<ndim> *sph;                       ///< SPH algorithm pointer
  SphIntegration<ndim> *sphint;         ///< SPH Integration scheme pointer
  SphNeighbourSearch<ndim> *sphneib;    ///< SPH Neighbour scheme pointer
#ifdef MPI_PARALLEL
  MpiControl<ndim> mpicontrol;          ///< MPI control object
#endif

};



//=============================================================================
//  Class SphSimulation
/// \brief   Main SphSimulation class.
/// \details Main SphSimulation class definition, inherited from Simulation, 
///          which controls the main program flow for SPH simulations.
/// \author  D. A. Hubber, G. Rosotti
/// \date    03/04/2013
//=============================================================================
template <int ndim>
class SphSimulation : public Simulation<ndim> 
{
  using SimulationBase::simparams;
  using Simulation<ndim>::sph;
  using Simulation<ndim>::nbody;
  using Simulation<ndim>::sinks;
  using Simulation<ndim>::subsystem;
  using Simulation<ndim>::nbodytree;
  using Simulation<ndim>::sphint;
  using Simulation<ndim>::uint;
  using Simulation<ndim>::sphneib;
  using Simulation<ndim>::ghosts;
  using Simulation<ndim>::simbox;
  using Simulation<ndim>::simunits;
  using Simulation<ndim>::Nstepsmax;
  using Simulation<ndim>::run_id;
  using Simulation<ndim>::out_file_form;
  using Simulation<ndim>::tend;
  using Simulation<ndim>::noutputstep;
  using Simulation<ndim>::nbody_single_timestep;
  using Simulation<ndim>::ParametersProcessed;
  using Simulation<ndim>::n;
  using Simulation<ndim>::Nlevels;
  using Simulation<ndim>::Nsteps;
  using Simulation<ndim>::t;
  using Simulation<ndim>::timestep;
  using Simulation<ndim>::level_step;
  using Simulation<ndim>::Noutsnap;
  using Simulation<ndim>::tsnapnext;
  using Simulation<ndim>::dt_snap;
  using Simulation<ndim>::dt_python;
  using Simulation<ndim>::level_diff_max;
  using Simulation<ndim>::level_max;
  using Simulation<ndim>::integration_step;
  using Simulation<ndim>::nresync;
  using Simulation<ndim>::dt_max;
  using Simulation<ndim>::sph_single_timestep;
  using Simulation<ndim>::sink_particles;

public:

  SphSimulation (Parameters* parameters): Simulation<ndim>(parameters) {};
  virtual void PostInitialConditionsSetup(void);
  virtual void MainLoop(void);
  virtual void ComputeGlobalTimestep(void);
  virtual void ComputeBlockTimesteps(void);
  //virtual void ProcessParameters(void);

};



//=============================================================================
//  Class GodunovSphSimulation
/// \brief   Main GodunovSphSimulation class.
/// \details Main GodunovSphSimulation class definition, inherited from 
///          Simulation, which controls the main program flow for Godunov 
///          SPH simulations (Inutsuka 2002).
/// \author  D. A. Hubber, G. Rosotti
/// \date    03/04/2013
//=============================================================================
template <int ndim>
class GodunovSphSimulation : public Simulation<ndim> 
{
  using SimulationBase::simparams;
  using Simulation<ndim>::sph;
  using Simulation<ndim>::nbody;
  using Simulation<ndim>::sinks;
  using Simulation<ndim>::subsystem;
  using Simulation<ndim>::nbodytree;
  using Simulation<ndim>::sphint;
  using Simulation<ndim>::uint;
  using Simulation<ndim>::sphneib;
  using Simulation<ndim>::ghosts;
  using Simulation<ndim>::simbox;
  using Simulation<ndim>::simunits;
  using Simulation<ndim>::Nstepsmax;
  using Simulation<ndim>::run_id;
  using Simulation<ndim>::out_file_form;
  using Simulation<ndim>::tend;
  using Simulation<ndim>::noutputstep;
  using Simulation<ndim>::nbody_single_timestep;
  using Simulation<ndim>::ParametersProcessed;
  using Simulation<ndim>::n;
  using Simulation<ndim>::Nlevels;
  using Simulation<ndim>::Nsteps;
  using Simulation<ndim>::t;
  using Simulation<ndim>::timestep;
  using Simulation<ndim>::level_step;
  using Simulation<ndim>::Noutsnap;
  using Simulation<ndim>::tsnapnext;
  using Simulation<ndim>::dt_snap;
  using Simulation<ndim>::dt_python;
  using Simulation<ndim>::level_diff_max;
  using Simulation<ndim>::level_max;
  using Simulation<ndim>::integration_step;
  using Simulation<ndim>::nresync;
  using Simulation<ndim>::dt_max;
  using Simulation<ndim>::sph_single_timestep;
  using Simulation<ndim>::sink_particles;

public:

  GodunovSphSimulation (Parameters* parameters): 
    Simulation<ndim>(parameters) {};
  virtual void PostInitialConditionsSetup(void);
  virtual void MainLoop(void);
  virtual void ComputeGlobalTimestep(void);
  virtual void ComputeBlockTimesteps(void);
  //virtual void ProcessParameters(void);
};



//=============================================================================
//  Class NbodySimulation
/// \brief   Main class for running N-body only simulations.
/// \details Main NbodySimulation class definition, inherited from Simulation, 
///          which controls the main program flow for N-body only simulations.
/// \author  D. A. Hubber, G. Rosotti
/// \date    03/04/2013
//=============================================================================
template <int ndim>
class NbodySimulation : public Simulation<ndim> 
{
  using SimulationBase::simparams;
  using Simulation<ndim>::sph;
  using Simulation<ndim>::nbody;
  using Simulation<ndim>::subsystem;
  using Simulation<ndim>::nbodytree;
  using Simulation<ndim>::ghosts;
  using Simulation<ndim>::simbox;
  using Simulation<ndim>::simunits;
  using Simulation<ndim>::Nstepsmax;
  using Simulation<ndim>::run_id;
  using Simulation<ndim>::out_file_form;
  using Simulation<ndim>::tend;
  using Simulation<ndim>::noutputstep;
  using Simulation<ndim>::nbody_single_timestep;
  using Simulation<ndim>::ParametersProcessed;
  using Simulation<ndim>::n;
  using Simulation<ndim>::Nlevels;
  using Simulation<ndim>::Nsteps;
  using Simulation<ndim>::t;
  using Simulation<ndim>::timestep;
  using Simulation<ndim>::level_step;
  using Simulation<ndim>::Noutsnap;
  using Simulation<ndim>::tsnapnext;
  using Simulation<ndim>::dt_snap;
  using Simulation<ndim>::dt_python;
  using Simulation<ndim>::level_max;
  using Simulation<ndim>::integration_step;
  using Simulation<ndim>::nresync;
  using Simulation<ndim>::dt_max;

public:

  NbodySimulation (Parameters* parameters): Simulation<ndim>(parameters) {};
  virtual void PostInitialConditionsSetup(void);
  virtual void MainLoop(void);
  virtual void ComputeGlobalTimestep(void);
  virtual void ComputeBlockTimesteps(void);
  //virtual void ProcessParameters(void);
};

#endif


#endif
