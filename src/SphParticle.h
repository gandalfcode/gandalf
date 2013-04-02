// ============================================================================
// SphParticle.h
// Main SPH particle data structure
// ============================================================================


#ifndef _SPH_PARTICLE_H_
#define _SPH_PARTICLE_H_


#include "Dimensions.h"
#include "Precision.h"
#include "Constants.h"

template <int ndim>
struct SphParticle {
    bool active;
    //bool active_grav;
    int iorig;
    int itype;
    int level;
    FLOAT r[ndim];
    FLOAT v[ndim];
    FLOAT a[ndim];
    FLOAT r0[ndim];
    FLOAT v0[ndim];
    FLOAT a0[ndim];
    FLOAT agrav[ndim];
    FLOAT u;
    FLOAT u0;
    FLOAT dudt;
    FLOAT dudt0;
    FLOAT m;
    FLOAT h;
    FLOAT invh;
    FLOAT hfactor;
    FLOAT rho;
    FLOAT invrho;
    FLOAT press;
    FLOAT pfactor;
    FLOAT div_v;
    FLOAT invomega;
    FLOAT zeta;
    FLOAT q;
    FLOAT invq;
    FLOAT sound;
    FLOAT gpot;
    DOUBLE dt;
    FLOAT gradP[ndim];
    FLOAT gradrho[ndim];
    FLOAT gradv[ndim][ndim];
  	FLOAT rhomax;
  	FLOAT rhomin;
  	FLOAT pressmax;
  	FLOAT pressmin;
  	FLOAT vmax[ndim];
  	FLOAT vmin[ndim];


    SphParticle()
    {
      active = false;
      iorig = -1;
      itype = -1;
      level = 0;
      for (int k=0; k<ndim; k++) r[k] = (FLOAT) 0.0;
      for (int k=0; k<ndim; k++) v[k] = (FLOAT) 0.0;
      for (int k=0; k<ndim; k++) a[k] = (FLOAT) 0.0;
      for (int k=0; k<ndim; k++) r0[k] = (FLOAT) 0.0;
      for (int k=0; k<ndim; k++) v0[k] = (FLOAT) 0.0;
      for (int k=0; k<ndim; k++) a0[k] = (FLOAT) 0.0;
      for (int k=0; k<ndim; k++) agrav[k] = (FLOAT) 0.0;
      u = (FLOAT) 0.0;
      u0 = (FLOAT) 0.0;
      dudt = (FLOAT) 0.0;
      dudt0 = (FLOAT) 0.0;
      m = (FLOAT) 0.0;
      h = (FLOAT) 0.0;
      invh = (FLOAT) 0.0;
      hfactor = (FLOAT) 0.0;
      rho = (FLOAT) 0.0;
      invrho = (FLOAT) 0.0;
      press = (FLOAT) 0.0;
      pfactor = (FLOAT) 0.0;
      invomega = (FLOAT) 0.0;
      zeta = (FLOAT) 0.0;
      q = (FLOAT) 0.0;
      invq = (FLOAT) q;
      sound = (FLOAT) 0.0;
      for (int k=0; k<ndim; k++) gradP[k] = (FLOAT) 0.0;
      for (int k=0; k<ndim; k++) gradrho[k] = (FLOAT) 0.0;
      for (int k=0; k<ndim; k++)
        for (int kk=0; kk<ndim; kk++) gradv[k][kk] = (FLOAT) 0.0;
      dt = (DOUBLE) 0.0;
      rhomin = (FLOAT) 0.0;
      rhomax = (FLOAT) 0.0;
      pressmin = (FLOAT) 0.0;
      pressmax = (FLOAT) 0.0;
      for (int k=0; k<ndimmax; k++) vmax[k] = 0.0;
      for (int k=0; k<ndimmax; k++) vmin[k] = 0.0;
    }

  };


#endif
