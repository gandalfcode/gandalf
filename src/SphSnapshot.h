//=============================================================================
//  SphSnapshot.h
//  Contains definitions for SphSnapshot class
//=============================================================================


#ifndef _SPH_SNAPSHOT_H_
#define _SPH_SNAPSHOT_H_

#include <ctime>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include "Precision.h"
#include "Sph.h"
#include "Simulation.h"
#include "UnitInfo.h"
#include "BinaryOrbit.h"
using namespace std;


//=============================================================================
/// Class SphSnapshotBase
/// \brief   Definition for Sph Snapshot (Base) for recording simulation data.
/// \details Snapshot class contains a copy of the simulation data, 
///          either from main memory or from a loaded file.
/// \author  D. A. Hubber, G. Rosotti
/// \date    03/04/2013
//=============================================================================
class SphSnapshotBase
{
  void AllocateBufferMemoryBinary();
  void AllocateBufferMemorySph();
  void AllocateBufferMemoryStar();
  void DeallocateBufferMemoryBinary();
  void DeallocateBufferMemorySph();
  void DeallocateBufferMemoryStar();

protected:
  int nneededbinary;        ///< No. of variables needed to store binary orbit
  int nneededsph;           ///< No. of variables needed to store for sph ptcl
  int nneededstar;          ///< No. of variables needed to store for star ptcl
  vector<string> _species;

 public:

  static SphSnapshotBase* SphSnapshotFactory(string filename, SimulationBase* sim, int ndim);


  SphSnapshotBase(SimUnits*, string="");
  virtual ~SphSnapshotBase();


  // Snapshot function prototypes
  // --------------------------------------------------------------------------
  void AllocateBufferMemory(void);
  void DeallocateBufferMemory(void);
  int CalculateMemoryUsage(void);
  int CalculatePredictedMemoryUsage(void);
  virtual void CopyDataFromSimulation()=0;
  UnitInfo ExtractArray(string, string, float** out_array, int* size_array,
		    float& scaling_factor, string RequestedUnit);
  virtual void ReadSnapshot(string)=0;
  int GetNTypes() {return _species.size(); };
  string GetSpecies(int ispecies) { return _species.at(ispecies); };
  string GetRealType(string);
  int GetNparticlesType(string species);

  // All variables
  // --------------------------------------------------------------------------
  bool allocated;                   ///< Is snapshot memory allocated?
  bool allocatedbinary;             ///< Is SPH particle memory allocated?
  bool allocatedsph;                ///< Is SPH particle memory allocated?
  bool allocatedstar;               ///< Is star particle memory allocated?
  bool computedbinary;              ///< Are binary properties computed?
  bool computedsph;                 ///< Are additional SPH values computed?
  bool computednbody;               ///< Are additional star values computed?
  int LastUsed;                     ///< ??
  int nallocatedbinary;             ///< No. of floats allocated for SPH
  int nallocatedsph;                ///< No. of floats allocated for SPH
  int nallocatedstar;               ///< No. of floats allocated for stars
  int ndim;                         ///< Local copy of ndim
  int Nbinary;                      ///< No. of binary stars
  int Norbit;                       ///< No. of orbits in memory
  int Norbitmax;                    ///< Max. no. of orbits
  int Nquadruple;                   ///< No. of quadruple systems
  int Nsph;                         ///< No. of SPH particles
  int Nsphmax;                      ///< Max. no. of SPH particles
  int Nstar;                        ///< No. of star particles
  int Nstarmax;                     ///< Max. no. of star particles
  int Ntriple;                      ///< No. of triple systems
  DOUBLE t;                         ///< Simulation time of snapshot

  string filename;                  ///< Filename of snapshot
  string fileform;                  ///< File format of snapshot
  //string unitname;                  ///< Aux. unit string
  string label;                     ///< Aux. latex label

  SimUnits* units;                  ///< Pointer to units object


  // Pointers for allocating memory required for storing all important
  // snapshot data
  // --------------------------------------------------------------------------
  float *x;                         ///< x-position for SPH particles
  float *y;                         ///< y-position for SPH particles
  float *z;                         ///< z-position for SPH particles
  float *vx;                        ///< x-velocity for SPH particles
  float *vy;                        ///< y-velocity for SPH particles
  float *vz;                        ///< z-velocity for SPH particles
  float *ax;                        ///< x-acceleration for SPH particles
  float *ay;                        ///< y-acceleration for SPH particles
  float *az;                        ///< z-acceleration for SPH particles
  float *m;                         ///< Masses for SPH particles
  float *h;                         ///< Smoothing lengths for SPH particles
  float *rho;                       ///< Density for SPH particles
  float *u;                         ///< Specific int. energy for SPH particles
  float *dudt;                      ///< Heating/cooling rate for SPH particles

  float *xstar;                     /// x-position for star particles
  float *ystar;                     /// y-position for star particles
  float *zstar;                     /// z-position for star particles
  float *vxstar;                    /// x-velocity for star particles
  float *vystar;                    /// y-velocity for star particles
  float *vzstar;                    /// z-velocity for star particles
  float *axstar;                    /// x-acceleration for star particles
  float *aystar;                    /// y-acceleration for star particles
  float *azstar;                    /// z-acceleration for star particles
  float *mstar;                     /// Masses for star particles
  float *hstar;                     /// Smoothing length for star particles

  float *ecc;                       /// Binary orbital eccentricity
  float *mbin;                      /// Binary total mass 
  float *period;                    /// Binary orbital period
  float *qbin;                      /// Binary mass ratio
  float *sma;                       /// Binary orbital semi-major axis

};



//=============================================================================
/// Class SphSnapshot
/// \brief   Definition for Sph Snapshot for recording simulation data.
/// \details Snapshot class contains a copy of the simulation data, 
///          either from main memory or from a loaded file.
/// \author  D. A. Hubber, G. Rosotti
/// \date    03/04/2013
//=============================================================================
template <int ndims>
class SphSnapshot : public SphSnapshotBase {
public:
  SphSnapshot (string, SimulationBase* );
  ~SphSnapshot() {};
  void CopyDataFromSimulation();
  void ReadSnapshot(string);

  Simulation<ndims>* simulation;
};
#endif
