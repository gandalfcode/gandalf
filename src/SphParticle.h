// ============================================================================
// SphParticle.h
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
  FLOAT a0[vdimmax];
  FLOAT agrav[ndimmax];
  FLOAT u;
  FLOAT u0;
  FLOAT dudt;
  FLOAT dudt0;
  FLOAT m;
  FLOAT h;
  FLOAT invh;
  FLOAT rho;
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
    for (int k=0; k<ndimmax; k++) r[k] = 0.0;
    for (int k=0; k<ndimmax; k++) v[k] = 0.0;
    for (int k=0; k<ndimmax; k++) a[k] = 0.0;
    for (int k=0; k<ndimmax; k++) a0[k] = 0.0;
    for (int k=0; k<ndimmax; k++) agrav[k] = 0.0;
    m = 0;
    h = 0;
    invh = 0.0;
    rho = 0.0;
    u = 0.0;
    dudt = 0.0;
    invomega = 0.0;
    zeta = 0.0;
    dt = 0.0;
  } 

};


#endif
