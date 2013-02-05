// ============================================================================
// SphSnapshot.h
// ============================================================================


#ifndef _SPH_SNAPSHOT_H_
#define _SPH_SNAPSHOT_H_


#include <iostream>
#include <map>
#include <string>
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
  void ExtractArray(std::string, float** out_array, int* size_array);
  void ReadSnapshot(string, SphSimulation *);

  bool allocated;
  int nallocated;
  int ndim;
  int Nsph;
  float t;
  std::string filename;
  std::string fileform;


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

};


#endif
