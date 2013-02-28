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
using namespace std;


// ============================================================================
// Class SphSnapshot
// ============================================================================
class SphSnapshot
{
 public:

  SphSnapshot(string="");
  ~SphSnapshot();

  void AllocateBufferMemory(void);
  void DeallocateBufferMemory(void);
  int CalculateMemoryUsage(void);
  void CopyDataFromSimulation(int,int,SphParticle*);
  void ExtractArray(string, float** out_array, int* size_array, float& scaling_factor, string RequestedUnit);
  void ReadSnapshot(string, SphSimulation *);

  bool allocated;
  int nallocated;
  int ndim;
  int Nsph;
  int Nmax;
  DOUBLE t;
  std::string filename;
  std::string fileform;
  int LastUsed;

  SimUnits* units;

  string unitname;
  string label;


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


#endif
