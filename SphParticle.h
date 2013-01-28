// ============================================================================
// SphParticle.h
// ============================================================================


#ifndef _SPH_PARTICLE_H_
#define _SPH_PARTICLE_H_


#include "Dimensions.h"


struct SphParticle {
  //SphParticle();
  int iorig;
  int itype;
  float r[ndimmax];
  float v[vdimmax];
  float a[vdimmax];
  float r0[ndimmax];
  float v0[vdimmax];
  float a0[vdimmax];
  float m;
  float h;
  float invh;
  float rho;
  float invrho;
  float hfactor;
  float pfactor;
  float div_v;
  float u;
  float dudt;
  float invomega;
  float zeta;
  float sound;
};

/*
SphParticle::SphParticle()
{
  iorig = -1;
  itype = -1;
  for (int k=0; k<ndimmax; k++) r[k] = 0.0;
  for (int k=0; k<ndimmax; k++) v[k] = 0.0;
  for (int k=0; k<ndimmax; k++) a[k] = 0.0;
  for (int k=0; k<ndimmax; k++) r0[k] = 0.0;
  for (int k=0; k<ndimmax; k++) v0[k] = 0.0;
  for (int k=0; k<ndimmax; k++) a0[k] = 0.0;
  m = 0;
  h = 0;
  invh = 0.0;
  rho = 0.0;
  invrho = 0.0;
  hfactor = 0.0;
  u = 0.0;
  dudt = 0.0;
  invomega = 0.0;
  zeta = 0.0;
} 
*/

#endif
