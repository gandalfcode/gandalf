// ============================================================================
// SphParticle.h
// Main SPH particle data structure
// ============================================================================


#ifndef _SPH_PARTICLE_H_
#define _SPH_PARTICLE_H_


#include "Dimensions.h"
#include "Precision.h"


struct SphParticle {
  bool active;
  //bool active_grav;
  int iorig;
  int itype;
  int level;
  FLOAT r[ndimmax];
  FLOAT v[vdimmax];
  FLOAT a[vdimmax];
  FLOAT r0[vdimmax];
  FLOAT v0[vdimmax];
  FLOAT a0[vdimmax];
  FLOAT agrav[ndimmax];
  FLOAT u;
  FLOAT u0;
  FLOAT dudt;
  FLOAT dudt0;
  FLOAT m;
  FLOAT h;
  FLOAT invh;
  FLOAT hfactor;
  FLOAT pfactor;
  FLOAT rho;
  FLOAT invrho;
  FLOAT div_v;
  FLOAT invomega;
  FLOAT zeta;
  FLOAT sound;
  FLOAT gpot;
  DOUBLE dt;

  SphParticle()
  {
    active = false;
    iorig = -1;
    itype = -1;
    level = 0;
    for (int k=0; k<ndimmax; k++) r[k] = (FLOAT) 0.0;
    for (int k=0; k<ndimmax; k++) v[k] = (FLOAT) 0.0;
    for (int k=0; k<ndimmax; k++) a[k] = (FLOAT) 0.0;
    for (int k=0; k<ndimmax; k++) r0[k] = (FLOAT) 0.0;
    for (int k=0; k<ndimmax; k++) v0[k] = (FLOAT) 0.0;
    for (int k=0; k<ndimmax; k++) a0[k] = (FLOAT) 0.0;
    for (int k=0; k<ndimmax; k++) agrav[k] = (FLOAT) 0.0;
    u = (FLOAT) 0.0;
    u0 = (FLOAT) 0.0;
    dudt = (FLOAT) 0.0;
    dudt0 = (FLOAT) 0.0;
    m = (FLOAT) 0.0;
    h = (FLOAT) 0.0;
    invh = (FLOAT) 0.0;
    hfactor = (FLOAT) 0.0;
    pfactor = (FLOAT) 0.0;
    rho = (FLOAT) 0.0;
    invrho = (FLOAT) 0.0;
    invomega = (FLOAT) 0.0;
    zeta = (FLOAT) 0.0;
    sound = (FLOAT) 0.0;
    dt = (FLOAT) 0.0;
  } 

};


#endif
