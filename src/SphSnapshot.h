// ============================================================================
// SphSnapshot.h
// ============================================================================


#ifndef _SPH_SNAPSHOT_H_
#define _SPH_SNAPSHOT_H_

#include <ctime>
#include <iostream>
#include <map>
#include <string>
#include "Precision.h"
#include "Sph.h"
#include "SphSimulation.h"
#include "UnitInfo.h"
using namespace std;


// ============================================================================
// Class SphSnapshot
// Snapshot class contains a copy of the simulation data, either from main
// memory or from a loaded file.
// ============================================================================
class SphSnapshotBase
{
 public:

  static SphSnapshotBase* SphSnapshotFactory(string filename, SimulationBase* sim, int ndim);


  SphSnapshotBase(string="");
  ~SphSnapshotBase();


  // Snapshot function prototypes
  // --------------------------------------------------------------------------
  void AllocateBufferMemory(void);
  void DeallocateBufferMemory(void);
  int CalculateMemoryUsage(void);
  virtual void CopyDataFromSimulation()=0;
  UnitInfo ExtractArray(string, float** out_array, int* size_array,
		    float& scaling_factor, string RequestedUnit);
  virtual void ReadSnapshot(string)=0;


  // All variables
  // --------------------------------------------------------------------------
  bool allocated;                           // Is snapshot memory allocated?
  int LastUsed;                             // ??
  int nallocated;                           // No. of floats allocated
  int ndim;                                 // Local copy of ndim
  int Nsph;                                 // No. of SPH particles
  int Nmax;                                 // Max. no. of SPH particles
  DOUBLE t;                                 // Simulation time of snapshot

  string filename;                          // Filename of snapshot
  string fileform;                          // File format of snapshot
//  string unitname;                          // Aux. unit string
  string label;                             // Aux. latex label

  SimUnits* units;                          // Pointer to units object


  // Pointers for allocating memory required for storing all important
  // snapshot data
  // --------------------------------------------------------------------------
  float *x;
  float *y;
  float *z;
  float *vx;
  float *vy;
  float *vz;
  float *ax;
  float *ay;
  float *az;
  float *m;
  float *h;
  float *rho;
  float *u;
  float *dudt;

};

template <int ndims>
class SphSnapshot : public SphSnapshotBase {
public:
  SphSnapshot (string, SimulationBase* );
  void CopyDataFromSimulation();
  void ReadSnapshot(string);

  Simulation<ndims>* simulation;
};
#endif
