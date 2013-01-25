// ============================================================================
// SNAPSHOT.H
// ============================================================================


#ifndef _SphSnapshot_H_
#define _SphSnapshot_H_


#include <iostream>
#include <map>
#include <string>
#include "Sph.h"
using namespace std;


// ============================================================================
// CLASS SNAPSHOT
// ============================================================================
class SphSnapshot
{
 public:

  SphSnapshot();
  ~SphSnapshot();

  void AllocateBufferMemory(void);
  void DeallocateBufferMemory(void);
  int CalculateMemoryUsage(void);
  void CopyDataFromSimulation(int,int,SphParticle*);
  void ExtractArray(std::string, float** out_array, int* size_array);

  bool allocated;
  int nallocated;
  int ndim;
  int Nsph;
  float t;
  std::string filename;
  std::string fileform;

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
