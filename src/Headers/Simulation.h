//=================================================================================================
//  Simulation.h
//  Contains definitions for following simulation classes:
//  - SimulationBase
//  - Simulation
//  - SphSimulation
//  - GradhSphSimulation
//  - SM2012SphSimulation
//  - MeshlessFVSimulation
//  - NbodySimulation
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
//=================================================================================================


#ifndef _SIMULATION_H_
#define _SIMULATION_H_


#include <map>
#include <string>
#include <list>
#include "CodeTiming.h"
#include "Diagnostics.h"
#include "DomainBox.h"
#include "Dust.h"
#include "Ewald.h"
#include "ExternalPotential.h"
#include "Hydrodynamics.h"
#include "MeshlessFV.h"
#include "MfvNeighbourSearch.h"
#include "Precision.h"
#include "Parameters.h"
#include "Radiation.h"
#include "RandomNumber.h"
#include "SimUnits.h"
#include "SmoothingKernel.h"
#include "Sph.h"
#include "SphNeighbourSearch.h"
#include "SphIntegration.h"
#include "TreeRay.h"
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


//=================================================================================================
//  Class SimulationBase
/// \brief  Creates a simulation object depending on the dimensionality.
/// \author D. A. Hubber, G. Rosotti
/// \date   03/04/2013
//=================================================================================================
class SimulationBase
{
  // Subroutines only for internal use of the class
  virtual void CalculateDiagnostics(void)=0;
  virtual void OutputDiagnostics(void)=0;
  virtual void OutputTestDiagnostics(void)=0;
  virtual void RecordDiagnostics(void)=0;
  virtual void UpdateDiagnostics(void)=0;
  virtual void GenerateIC(void)=0;
  virtual void ReadColumnHeaderFile(ifstream& infile, HeaderInfo& info)=0;
  virtual bool ReadColumnSnapshotFile(string)=0;
  virtual bool WriteColumnSnapshotFile(string)=0;
  virtual void ReadSerenFormHeaderFile(ifstream& infile, HeaderInfo& info)=0;
  virtual bool ReadSerenFormSnapshotFile(string)=0;
  virtual bool WriteSerenFormSnapshotFile(string)=0;
  virtual void ReadSerenUnformHeaderFile(ifstream& infile, HeaderInfo& info)=0;
  virtual bool ReadSerenUnformSnapshotFile(string)=0;
  virtual bool WriteSerenUnformSnapshotFile(string)=0;
  virtual bool WriteSerenLiteSnapshotFile(string)=0;

  std::list<string> keys;

 public:

  static SimulationBase* SimulationFactory(int ndim, string simtype, Parameters* params);


  // Constructor and Destructor
  //-----------------------------------------------------------------------------------------------
  SimulationBase(Parameters* params);
  virtual ~SimulationBase();

  // Subroutine prototypes
  //-----------------------------------------------------------------------------------------------
  list<string>* GetIntAndFloatParameterKeys();
  string GetParam(string key);
  list<SphSnapshotBase*> InteractiveRun(int=-1);
  string Output(void);
  void RestartSnapshot(void);
  void Run(int=-1);
  void SetParam(string key, string value);
  void SetParam(string key, int value);
  void SetParam(string ket, double value);
  void SetupSimulation(void);
  void SplashScreen(void);


  virtual void ImportArray(double* input, int size, string quantity, string type="sph") = 0;
  virtual void MainLoop(void)=0;
  virtual void PostInitialConditionsSetup(void)=0;
  virtual void PreSetupForPython(void)=0;
  virtual void ProcessParameters(void)=0;
  virtual void SetComFrame(void)=0;
  virtual void FinaliseSimulation(void) {};


  // Input-output routines
  //-----------------------------------------------------------------------------------------------
  bool ReadSnapshotFile(string,string);
  bool WriteSnapshotFile(string,string);
  HeaderInfo ReadHeaderSnapshotFile(string filename, string format);


  // Variables
  //-----------------------------------------------------------------------------------------------
  bool ewaldGravity;                   ///< Flag if periodic graivty is being used
  bool extra_sink_output;              ///< Create extra output files for sink particles
  bool initial_h_provided;             ///< Have initial h values been calculated?
  bool kill_simulation;                ///< Kill simulation flag
  bool ParametersProcessed;            ///< Flag if params are already processed
  bool periodicBoundaries;             ///< Flag if periodic boundaries are being used
  bool rebuild_tree;                   ///< Flag to rebuild neighbour tree
  bool recomputeRadiation;             ///< Flag to recompute radiation field from all sources
  bool rescale_particle_data;          ///< Flag to scale data to code units
  bool restart;                        ///< Flag to restart from last snapshot
  bool setup;                          ///< Flag if simulation is setup
  int integration_step;                ///< Steps per complete integration step
  int litesnap;                        ///< Activate lite snapshots (for movies)
  int nbody_single_timestep;           ///< Flag if stars use same timestep
  int ndims;                           ///< Aux. dimensionality variable.
                                       ///< Required for python routines.
  int Nmpi;                            ///< No. of MPI processes
  int Noutsnap;                        ///< No. of output snapshots
  int Noutlitesnap;                    ///< No. of lite output snapshots
  int Nthreads;                        ///< Max no. of (OpenMP) threads
  int pruning_level_min;               ///< Min. level of pruned trees for MPI
  int pruning_level_max;               ///< Max. level of pruned trees for MPI
  int rank;                            ///< Process i.d. (for MPI simulations)
  int sink_particles;                  ///< Switch on sink particles
  int sph_single_timestep;             ///< Flag if SPH ptcls use same step
  int level_diff_max;                  ///< Max. allowed neib timestep level diff
  int level_max;                       ///< Maximum timestep level
  int level_step;                      ///< Level of smallest timestep unit
  int n;                               ///< Integer time counter
  int ndiagstep;                       ///< Diagnostic output frequency (in units
                                       ///< of full block timestep steps)
  int nlastrestart;                    ///< Integer time of last restart snapshot
  int noutputstep;                     ///< Output frequency
  int nradstep;                        ///< Integer time between computing radiation field
  int nrestartstep;                    ///< Integer time between creating temp restart files
  int nresync;                         ///< Integer time for resynchronisation
  int nsystembuildstep;                ///< Integer time between rebuilding N-body system tree
  int ntreebuildstep;                  ///< Integer time between rebuilding tree
  int ntreestockstep;                  ///< Integer time between restocking tree
  int Nblocksteps;                     ///< No. of full block timestep steps
  int Nlevels;                         ///< No. of timestep levels
  int Nfullsteps;                      ///< No. of full steps in simulation
  int Nstepsmax;                       ///< Max. allowed no. of steps
  int Nsteps;                          ///< Total no. of steps in simulation
  DOUBLE dt_litesnap;                  ///< Lite-snapshot time interval
  DOUBLE dt_max;                       ///< Value of maximum timestep level
  DOUBLE dt_min_nbody;                 ///< Minimum timestep of all N-body particles
  DOUBLE dt_min_hydro;                 ///< Minimum timestep of all SPH particles
  DOUBLE dt_python;                    ///< Python window update time interval
  DOUBLE dt_snap;                      ///< Snapshot time interval
  DOUBLE dt_snap_wall;                 ///< Wallclock time between snapshots
  DOUBLE t;                            ///< Current simulation time
  DOUBLE tmax_wallclock;               ///< Max. wall-clock time for sim (in s)
  DOUBLE tend;                         ///< End time of simulation
  DOUBLE timestep;                     ///< Current timestep
  DOUBLE tsnapfirst;                   ///< Time of first snapshot
  DOUBLE tsnaplast;                    ///< (Expected) time of last snaphsot
  DOUBLE tlitesnaplast;                ///< Time of last lite-snapshot
  DOUBLE tlitesnapnext;                ///< Time of next lite-snapshot
  DOUBLE tsnapnext;                    ///< Time of next snapshot
  DOUBLE tsnap_wallclock;              ///< Wallclock time of last snapshot
  string out_file_form;                ///< Output snapshot file format
  string paramfile;                    ///< Name of parameters file
  string run_id;                       ///< Simulation id string

  Parameters* simparams;               ///< Simulation parameters object (pointer)
  SimUnits simunits;                   ///< Simulation units object
  CodeTiming *timing;                  ///< Simulation timing object (pointer)

};



#if !defined(SWIG)
//=================================================================================================
//  Class Simulation
/// \brief   Main Simulation class.
/// \details Main parent Simulation class from which all other simulation
///          objects (e.g. SphSimulation) inherit from.
/// \author  D. A. Hubber, G. Rosotti
/// \date    03/04/2013
//=================================================================================================
template <int ndim>
class Simulation : public SimulationBase
{
  void ImportArraySph(double* input, int size, string quantity);
  void ImportArrayNbody(double* input, int size, string quantity);

 public:
  Simulation(Parameters* parameters) :
    SimulationBase(parameters),
    hydro(NULL),
    nbody(NULL)
    {this->ndims=ndim;};


  // Memory allocation routines
  //-----------------------------------------------------------------------------------------------
  void AllocateParticleMemory(void);
  void DeallocateParticleMemory(void);


  // Subroutine prototypes
  //-----------------------------------------------------------------------------------------------
  virtual void CalculateDiagnostics(void);
  virtual void ComputeGlobalTimestep(void)=0;
  virtual void ComputeBlockTimesteps(void)=0;
  virtual void GenerateIC(void);
  virtual void ImportArray(double* input, int size, string quantity, string type="sph");
  virtual void OutputDiagnostics(void);
  virtual void OutputTestDiagnostics(void);
  virtual void PreSetupForPython(void);
  virtual void ProcessNbodyParameters(void);
  virtual void ProcessParameters(void)=0;
  virtual void RecordDiagnostics(void);
  virtual void RegulariseParticleDistribution(const int) {};
  virtual void SetComFrame(void);
  virtual void SmoothParticleQuantity(const int, FLOAT *) {};
  virtual void UpdateDiagnostics(void);


  // Input-output routines
  //-----------------------------------------------------------------------------------------------
  virtual void ReadColumnHeaderFile(ifstream& infile, HeaderInfo& info);
  virtual bool ReadColumnSnapshotFile(string);
  virtual bool WriteColumnSnapshotFile(string);
  virtual void ReadSerenFormHeaderFile(ifstream& infile, HeaderInfo& info);
  virtual bool ReadSerenFormSnapshotFile(string);
  virtual bool WriteSerenFormSnapshotFile(string);
  virtual void ReadSerenUnformHeaderFile(ifstream& infile, HeaderInfo& info);
  virtual bool ReadSerenUnformSnapshotFile(string);
  virtual bool WriteSerenUnformSnapshotFile(string);
  virtual bool WriteSerenLiteSnapshotFile(string);
  virtual void ConvertToCodeUnits(void);


  // Variables
  //-----------------------------------------------------------------------------------------------
  static const int vdim=ndim;          ///< Velocity vector dimensionality (same as ndim)
  static const FLOAT invndim;          ///< Local copy of 1/ndim

  DomainBox<ndim> simbox;              ///< Simulation boundary data
  Diagnostics<ndim> diag0;             ///< Initial diagnostic state
  Diagnostics<ndim> diag;              ///< Current diagnostic state
  EnergyEquation<ndim> *uint;          ///< Energy equation pointer
  ExternalPotential<ndim> *extpot;     ///< Pointer to external potential object
  Ewald<ndim> *ewald;                  ///< Ewald periodic gravity object
  Ghosts<ndim>* LocalGhosts;           ///< Periodic ghost particle object
  Hydrodynamics<ndim> *hydro;          ///< Hydrodynamics algorithm pointer
  Nbody<ndim> *nbody;                  ///< N-body algorithm pointer
  Nbody<ndim> *subsystem;              ///< N-body object for sub-systems
  NbodySystemTree<ndim> nbodytree;     ///< N-body tree to create sub-systems
  Radiation<ndim> *radiation;          ///< Radiation field object
  RandomNumber *randnumb;              ///< Random number object (pointer)
  Sinks<ndim> *sinks;                  ///< Sink particle object
  SphIntegration<ndim> *sphint;        ///< SPH Integration scheme pointer
  SphNeighbourSearch<ndim> *sphneib;   ///< SPH Neighbour scheme pointer
#ifdef MPI_PARALLEL
  MpiControl<ndim>* mpicontrol;        ///< MPI control object
  Ghosts<ndim>* MpiGhosts;             ///< MPI ghost particle object
#endif

};


// Declare invndim constant here (prevents warnings with some compilers)
//template <int ndim>
//const FLOAT Simulation<ndim>::invndim = 1.0/ndim;



//=================================================================================================
//  Class SphSimulation
/// \brief   Main SphSimulation class.
/// \details Main SphSimulation class definition, inherited from Simulation,
///          which controls the main program flow for SPH simulations.
/// \author  D. A. Hubber, G. Rosotti
/// \date    03/04/2013
//=================================================================================================
template <int ndim>
class SphSimulation : public Simulation<ndim>
{
 public:
  using Simulation<ndim>::ewaldGravity;
  using Simulation<ndim>::extra_sink_output;
  using Simulation<ndim>::periodicBoundaries;
  using Simulation<ndim>::Nmpi;
  using Simulation<ndim>::pruning_level_min;
  using Simulation<ndim>::pruning_level_max;
  using Simulation<ndim>::restart;
  using Simulation<ndim>::simparams;
  using Simulation<ndim>::timing;
  using Simulation<ndim>::extpot;
  using Simulation<ndim>::ewald;
  using Simulation<ndim>::kill_simulation;
  using Simulation<ndim>::hydro;
  using Simulation<ndim>::nbody;
  using Simulation<ndim>::sinks;
  using Simulation<ndim>::subsystem;
  using Simulation<ndim>::nbodytree;
  using Simulation<ndim>::sphint;
  using Simulation<ndim>::uint;
  using Simulation<ndim>::litesnap;
  using Simulation<ndim>::LocalGhosts;
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
  using Simulation<ndim>::Nblocksteps;
  using Simulation<ndim>::Nfullsteps;
  using Simulation<ndim>::Nlevels;
  using Simulation<ndim>::Nsteps;
  using Simulation<ndim>::t;
  using Simulation<ndim>::timestep;
  using Simulation<ndim>::level_step;
  using Simulation<ndim>::Noutsnap;
  using Simulation<ndim>::tlitesnapnext;
  using Simulation<ndim>::tsnapnext;
  using Simulation<ndim>::dt_litesnap;
  using Simulation<ndim>::dt_min_nbody;
  using Simulation<ndim>::dt_min_hydro;
  using Simulation<ndim>::dt_snap;
  using Simulation<ndim>::dt_python;
  using Simulation<ndim>::level_diff_max;
  using Simulation<ndim>::level_max;
  using Simulation<ndim>::integration_step;
  using Simulation<ndim>::nresync;
  using Simulation<ndim>::dt_max;
  using Simulation<ndim>::sph_single_timestep;
  using Simulation<ndim>::sink_particles;
  using Simulation<ndim>::randnumb;
  using Simulation<ndim>::rank;
  using Simulation<ndim>::rebuild_tree;
  using Simulation<ndim>::recomputeRadiation;
  using Simulation<ndim>::ndiagstep;
  using Simulation<ndim>::nradstep;
  using Simulation<ndim>::nrestartstep;
  using Simulation<ndim>::ntreebuildstep;
  using Simulation<ndim>::ntreestockstep;
  using Simulation<ndim>::tmax_wallclock;
  using Simulation<ndim>::sphneib;
  using Simulation<ndim>::radiation;
#ifdef MPI_PARALLEL
  using Simulation<ndim>::mpicontrol;
  using Simulation<ndim>::MpiGhosts;
#endif

  SphSimulation (Parameters* parameters)
  : Simulation<ndim>(parameters), sph(NULL), sphdust(NULL)
  {};
  virtual ~SphSimulation() {} ;

  virtual void ProcessSphParameters(void)=0;
  virtual void PostInitialConditionsSetup(void);
  virtual void MainLoop(void);
  virtual void ComputeGlobalTimestep(void);
  virtual void ComputeBlockTimesteps(void);
  virtual void ProcessParameters(void);
  virtual void WriteExtraSinkOutput(void);
  virtual void RegulariseParticleDistribution(const int);
  virtual void SmoothParticleQuantity(const int, FLOAT *);

  Sph<ndim> *sph;                      ///< SPH algorithm pointer
  DustBase<ndim>* sphdust ;               ///< Dust algorithm pointer
};



//=================================================================================================
//  Class GradhSphSimulation
/// \brief   grad-h SPH simulation class.
/// \details grad-h SPH Simulation class definition, inherited from SphSimulation, which controls
///          the main program flow for grad-h SPH simulations.
/// \author  D. A. Hubber, G. Rosotti
/// \date    17/03/2014
//=================================================================================================
template <int ndim>
class GradhSphSimulation: public SphSimulation<ndim>
{
 public:

  using Simulation<ndim>::Nmpi;
  using Simulation<ndim>::restart;
  using Simulation<ndim>::simparams;
  using Simulation<ndim>::timing;
  using Simulation<ndim>::extpot;
  using Simulation<ndim>::nbody;
  using Simulation<ndim>::subsystem;
  using Simulation<ndim>::nbodytree;
  using Simulation<ndim>::sphint;
  using Simulation<ndim>::uint;
  using Simulation<ndim>::LocalGhosts;
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
  using Simulation<ndim>::Nblocksteps;
  using Simulation<ndim>::Nlevels;
  using Simulation<ndim>::Nsteps;
  using Simulation<ndim>::randnumb;
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
  using Simulation<ndim>::rank;
  using Simulation<ndim>::rebuild_tree;
  using Simulation<ndim>::ndiagstep;
  using Simulation<ndim>::ntreebuildstep;
  using Simulation<ndim>::ntreestockstep;
  using SphSimulation<ndim>::radiation;
  using SphSimulation<ndim>::sinks;
  using SphSimulation<ndim>::sph;
  using SphSimulation<ndim>::sphneib;
  using SphSimulation<ndim>::tmax_wallclock;
  using SphSimulation<ndim>::sphdust ;
#ifdef MPI_PARALLEL
  using Simulation<ndim>::mpicontrol;
  using Simulation<ndim>::MpiGhosts;
#endif

  GradhSphSimulation (Parameters* parameters): SphSimulation<ndim>(parameters) {};
  virtual ~GradhSphSimulation() { } ;
  virtual void ProcessSphParameters(void);

};



//=================================================================================================
//  Class SM2012SphSimulation
/// \brief   Saitoh & Makino (2012) SPH simulation class.
/// \details Saitoh & Makino(2012) SPH Simulation class definition, inherited from SphSimulation,
///          which controls the main program flow for grad-h SPH simulations.
/// \author  D. A. Hubber, G. Rosotti
/// \date    17/03/2014
//=================================================================================================
template <int ndim>
class SM2012SphSimulation: public SphSimulation<ndim>
{
 public:

  using Simulation<ndim>::Nmpi;
  using Simulation<ndim>::restart;
  using Simulation<ndim>::simparams;
  using Simulation<ndim>::timing;
  using Simulation<ndim>::extpot;
  using Simulation<ndim>::nbody;
  using Simulation<ndim>::subsystem;
  using Simulation<ndim>::nbodytree;
  using Simulation<ndim>::sphint;
  using Simulation<ndim>::uint;
  using Simulation<ndim>::LocalGhosts;
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
  using Simulation<ndim>::Nblocksteps;
  using Simulation<ndim>::Nfullsteps;
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
  using Simulation<ndim>::rank;
  using Simulation<ndim>::rebuild_tree;
  using Simulation<ndim>::ndiagstep;
  using Simulation<ndim>::ntreebuildstep;
  using Simulation<ndim>::ntreestockstep;
  using SphSimulation<ndim>::sinks;
  using SphSimulation<ndim>::sph;
  using SphSimulation<ndim>::sphneib;
  using SphSimulation<ndim>::tmax_wallclock;
#ifdef MPI_PARALLEL
  using Simulation<ndim>::mpicontrol;
  using Simulation<ndim>::MpiGhosts;
#endif

  SM2012SphSimulation (Parameters* parameters): SphSimulation<ndim>(parameters) {};
  virtual ~SM2012SphSimulation() {} ;
  virtual void ProcessSphParameters(void);

};



//=================================================================================================
//  Class MeshlessFVSimulation
/// \brief   Main MeshlessFV Simulation class.
/// \details ...
/// \author  D. A. Hubber, J. Ngoumou, S. Heigl
/// \date    27/04/2015
//=================================================================================================
template <int ndim>
class MeshlessFVSimulation : public Simulation<ndim>
{
 public:
  using Simulation<ndim>::ewaldGravity;
  using Simulation<ndim>::extra_sink_output;
  using Simulation<ndim>::periodicBoundaries;
  using Simulation<ndim>::Nmpi;
  using Simulation<ndim>::pruning_level_min;
  using Simulation<ndim>::pruning_level_max;
  using Simulation<ndim>::restart;
  using Simulation<ndim>::simparams;
  using Simulation<ndim>::timing;
  using Simulation<ndim>::extpot;
  using Simulation<ndim>::ewald;
  using Simulation<ndim>::kill_simulation;
  using Simulation<ndim>::hydro;
  using Simulation<ndim>::nbody;
  using Simulation<ndim>::sinks;
  using Simulation<ndim>::subsystem;
  using Simulation<ndim>::nbodytree;
  using Simulation<ndim>::sphint;
  using Simulation<ndim>::uint;
  using Simulation<ndim>::litesnap;
  using Simulation<ndim>::LocalGhosts;
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
  using Simulation<ndim>::Nblocksteps;
  using Simulation<ndim>::Nfullsteps;
  using Simulation<ndim>::Nlevels;
  using Simulation<ndim>::Nsteps;
  using Simulation<ndim>::t;
  using Simulation<ndim>::timestep;
  using Simulation<ndim>::level_step;
  using Simulation<ndim>::Noutsnap;
  using Simulation<ndim>::tlitesnapnext;
  using Simulation<ndim>::tsnapnext;
  using Simulation<ndim>::dt_litesnap;
  using Simulation<ndim>::dt_min_nbody;
  using Simulation<ndim>::dt_min_hydro;
  using Simulation<ndim>::dt_snap;
  using Simulation<ndim>::dt_python;
  using Simulation<ndim>::level_diff_max;
  using Simulation<ndim>::level_max;
  using Simulation<ndim>::integration_step;
  using Simulation<ndim>::nresync;
  using Simulation<ndim>::dt_max;
  using Simulation<ndim>::sph_single_timestep;
  using Simulation<ndim>::sink_particles;
  using Simulation<ndim>::randnumb;
  using Simulation<ndim>::rank;
  using Simulation<ndim>::rebuild_tree;
  using Simulation<ndim>::ndiagstep;
  using Simulation<ndim>::nrestartstep;
  using Simulation<ndim>::ntreebuildstep;
  using Simulation<ndim>::ntreestockstep;
  using Simulation<ndim>::tmax_wallclock;
  using Simulation<ndim>::radiation;
#ifdef MPI_PARALLEL
  using Simulation<ndim>::mpicontrol;
  using Simulation<ndim>::MpiGhosts;
#endif

  MeshlessFVSimulation (Parameters* parameters): Simulation<ndim>(parameters) {};
  virtual void PostInitialConditionsSetup(void);
  virtual void MainLoop(void) = 0;
  virtual void ComputeGlobalTimestep(void);
  virtual void ComputeBlockTimesteps(void);
  virtual void ProcessParameters(void);
  virtual void WriteExtraSinkOutput(void);
  virtual void FinaliseSimulation(void);

  MeshlessFV<ndim> *mfv;                       ///< Meshless FV hydrodynamics algorithm pointer
  MeshlessFVNeighbourSearch<ndim> *mfvneib;    ///< Meshless FV neighbour search object pointer

};



//=================================================================================================
//  Class MfvMusclSimulation
/// \brief   Meshless Finite-volume simulation class using TVD Runge-Kutta integration
/// \details Meshless Finite-volume simulation class using TVD Runge-Kutta integration
/// \author  D. A. Hubber, J. Ngoumou
/// \date    27/04/2015
//=================================================================================================
template <int ndim>
class MfvMusclSimulation : public MeshlessFVSimulation<ndim>
{
 public:
  using Simulation<ndim>::ewaldGravity;
  using Simulation<ndim>::extra_sink_output;
  using Simulation<ndim>::periodicBoundaries;
  using Simulation<ndim>::Nmpi;
  using Simulation<ndim>::pruning_level_min;
  using Simulation<ndim>::pruning_level_max;
  using Simulation<ndim>::restart;
  using Simulation<ndim>::simparams;
  using Simulation<ndim>::timing;
  using Simulation<ndim>::extpot;
  using Simulation<ndim>::ewald;
  using Simulation<ndim>::kill_simulation;
  using Simulation<ndim>::hydro;
  using Simulation<ndim>::nbody;
  using Simulation<ndim>::sinks;
  using Simulation<ndim>::subsystem;
  using Simulation<ndim>::nbodytree;
  using Simulation<ndim>::sphint;
  using Simulation<ndim>::uint;
  using Simulation<ndim>::litesnap;
  using Simulation<ndim>::LocalGhosts;
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
  using Simulation<ndim>::Nblocksteps;
  using Simulation<ndim>::Nfullsteps;
  using Simulation<ndim>::Nlevels;
  using Simulation<ndim>::Nsteps;
  using Simulation<ndim>::t;
  using Simulation<ndim>::timestep;
  using Simulation<ndim>::level_step;
  using Simulation<ndim>::Noutsnap;
  using Simulation<ndim>::tlitesnapnext;
  using Simulation<ndim>::tsnapnext;
  using Simulation<ndim>::dt_litesnap;
  using Simulation<ndim>::dt_min_nbody;
  using Simulation<ndim>::dt_min_hydro;
  using Simulation<ndim>::dt_snap;
  using Simulation<ndim>::dt_python;
  using Simulation<ndim>::level_diff_max;
  using Simulation<ndim>::level_max;
  using Simulation<ndim>::integration_step;
  using Simulation<ndim>::nresync;
  using Simulation<ndim>::dt_max;
  using Simulation<ndim>::sph_single_timestep;
  using Simulation<ndim>::sink_particles;
  using Simulation<ndim>::randnumb;
  using Simulation<ndim>::rank;
  using Simulation<ndim>::rebuild_tree;
  using Simulation<ndim>::ndiagstep;
  using Simulation<ndim>::nrestartstep;
  using Simulation<ndim>::ntreebuildstep;
  using Simulation<ndim>::ntreestockstep;
  using Simulation<ndim>::tmax_wallclock;
  using Simulation<ndim>::radiation;
  using MeshlessFVSimulation<ndim>::mfv;
  using MeshlessFVSimulation<ndim>::mfvneib;
#ifdef MPI_PARALLEL
  using Simulation<ndim>::mpicontrol;
  using Simulation<ndim>::MpiGhosts;
#endif

  MfvMusclSimulation (Parameters* parameters) : MeshlessFVSimulation<ndim>(parameters) {};
  //virtual void PostInitialConditionsSetup(void);
  virtual void MainLoop(void);
  //virtual void ComputeGlobalTimestep(void);
  //virtual void ComputeBlockTimesteps(void);
  //virtual void ProcessParameters(void);
  //virtual void WriteExtraSinkOutput(void);

  //MeshlessFV<ndim> *mfv;                       ///< Meshless FV hydrodynamics algorithm pointer
  //MeshlessFVNeighbourSearch<ndim> *mfvneib;    ///< ..

};



//=================================================================================================
//  Class MfvRungeKutta
/// \brief   Meshless Finite-volume simulation class using TVD Runge-Kutta integration
/// \details Meshless Finite-volume simulation class using TVD Runge-Kutta integration
/// \author  D. A. Hubber, J. Ngoumou, S. Heigl
/// \date    27/04/2015
//=================================================================================================
template <int ndim>
class MfvRungeKuttaSimulation : public MeshlessFVSimulation<ndim>
{
 public:
  using Simulation<ndim>::ewaldGravity;
  using Simulation<ndim>::extra_sink_output;
  using Simulation<ndim>::periodicBoundaries;
  using Simulation<ndim>::Nmpi;
  using Simulation<ndim>::pruning_level_min;
  using Simulation<ndim>::pruning_level_max;
  using Simulation<ndim>::restart;
  using Simulation<ndim>::simparams;
  using Simulation<ndim>::timing;
  using Simulation<ndim>::extpot;
  using Simulation<ndim>::ewald;
  using Simulation<ndim>::kill_simulation;
  using Simulation<ndim>::hydro;
  //using Simulation<ndim>::mfv;
  //using Simulation<ndim>::mfvneib;
  using Simulation<ndim>::nbody;
  using Simulation<ndim>::sinks;
  using Simulation<ndim>::subsystem;
  using Simulation<ndim>::nbodytree;
  using Simulation<ndim>::sphint;
  using Simulation<ndim>::uint;
  using Simulation<ndim>::litesnap;
  using Simulation<ndim>::LocalGhosts;
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
  using Simulation<ndim>::Nblocksteps;
  using Simulation<ndim>::Nfullsteps;
  using Simulation<ndim>::Nlevels;
  using Simulation<ndim>::Nsteps;
  using Simulation<ndim>::t;
  using Simulation<ndim>::timestep;
  using Simulation<ndim>::level_step;
  using Simulation<ndim>::Noutsnap;
  using Simulation<ndim>::tlitesnapnext;
  using Simulation<ndim>::tsnapnext;
  using Simulation<ndim>::dt_litesnap;
  using Simulation<ndim>::dt_min_nbody;
  using Simulation<ndim>::dt_min_hydro;
  using Simulation<ndim>::dt_snap;
  using Simulation<ndim>::dt_python;
  using Simulation<ndim>::level_diff_max;
  using Simulation<ndim>::level_max;
  using Simulation<ndim>::integration_step;
  using Simulation<ndim>::nresync;
  using Simulation<ndim>::dt_max;
  using Simulation<ndim>::sph_single_timestep;
  using Simulation<ndim>::sink_particles;
  using Simulation<ndim>::randnumb;
  using Simulation<ndim>::rank;
  using Simulation<ndim>::rebuild_tree;
  using Simulation<ndim>::ndiagstep;
  using Simulation<ndim>::nrestartstep;
  using Simulation<ndim>::ntreebuildstep;
  using Simulation<ndim>::ntreestockstep;
  using Simulation<ndim>::tmax_wallclock;
  using Simulation<ndim>::radiation;
  using MeshlessFVSimulation<ndim>::mfv;
  using MeshlessFVSimulation<ndim>::mfvneib;
#ifdef MPI_PARALLEL
  using Simulation<ndim>::mpicontrol;
  using Simulation<ndim>::MpiGhosts;
#endif

  MfvRungeKuttaSimulation (Parameters* parameters) : MeshlessFVSimulation<ndim>(parameters) {};
  //virtual void PostInitialConditionsSetup(void);
  virtual void MainLoop(void);
  //virtual void ComputeGlobalTimestep(void);
  //virtual void ComputeBlockTimesteps(void);
  //virtual void ProcessParameters(void);
  //virtual void WriteExtraSinkOutput(void);

  //MeshlessFV<ndim> *mfv;                       ///< Meshless FV hydrodynamics algorithm pointer
  //MeshlessFVNeighbourSearch<ndim> *mfvneib;    ///< ..

};



//=================================================================================================
//  Class NbodySimulation
/// \brief   Main class for running N-body only simulations.
/// \details Main NbodySimulation class definition, inherited from Simulation,
///          which controls the main program flow for N-body only simulations.
/// \author  D. A. Hubber, G. Rosotti
/// \date    03/04/2013
//=================================================================================================
template <int ndim>
class NbodySimulation : public Simulation<ndim>
{
  using Simulation<ndim>::Nmpi;
  using Simulation<ndim>::restart;
  using Simulation<ndim>::simparams;
  using Simulation<ndim>::timing;
  using Simulation<ndim>::extpot;
  using Simulation<ndim>::nbody;
  using Simulation<ndim>::subsystem;
  using Simulation<ndim>::nbodytree;
  using Simulation<ndim>::hydro;
  using Simulation<ndim>::LocalGhosts;
  using Simulation<ndim>::simbox;
  using Simulation<ndim>::simunits;
  using Simulation<ndim>::Nstepsmax;
  using Simulation<ndim>::run_id;
  using Simulation<ndim>::out_file_form;
  using Simulation<ndim>::tend;
  using Simulation<ndim>::noutputstep;
  using Simulation<ndim>::nsystembuildstep;
  using Simulation<ndim>::nbody_single_timestep;
  using Simulation<ndim>::ParametersProcessed;
  using Simulation<ndim>::n;
  using Simulation<ndim>::Nblocksteps;
  using Simulation<ndim>::Nfullsteps;
  using Simulation<ndim>::Nlevels;
  using Simulation<ndim>::Nsteps;
  using Simulation<ndim>::ndiagstep;
  using Simulation<ndim>::randnumb;
  using Simulation<ndim>::tmax_wallclock;
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
  using Simulation<ndim>::nrestartstep;
  using Simulation<ndim>::dt_max;
  using Simulation<ndim>::sinks;

public:

  NbodySimulation (Parameters* parameters): Simulation<ndim>(parameters) {};
  virtual void PostInitialConditionsSetup(void);
  virtual void MainLoop(void);
  virtual void ComputeGlobalTimestep(void);
  virtual void ComputeBlockTimesteps(void);
  virtual void ProcessParameters(void);

};

#endif


#endif
