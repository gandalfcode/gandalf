// ============================================================================
// StarParticle.h
// ============================================================================


#ifndef _STAR_PARTICLE_H_
#define _STAR_PARTICLE_H_


#include "Dimensions.h"
#include "Precision.h"


struct StarParticle {
  bool active;
  int ilevel;
  DOUBLE r[ndimmax];
  DOUBLE v[ndimmax];
  DOUBLE a[ndimmax];
  DOUBLE adot[ndimmax];
  DOUBLE r0[ndimmax];
  DOUBLE v0[ndimmax];
  DOUBLE a0[ndimmax];
  DOUBLE adot0[ndimmax];
  DOUBLE m;
  DOUBLE h;
  DOUBLE invh;
  DOUBLE hfactor;
  DOUBLE gpot;

  StarParticle()
  {
    active = false;
    for (int k=0; k<ndimmax; k++) r[k] = 0.0;
    for (int k=0; k<ndimmax; k++) v[k] = 0.0;
    for (int k=0; k<ndimmax; k++) a[k] = 0.0;
    for (int k=0; k<ndimmax; k++) r0[k] = 0.0;
    for (int k=0; k<ndimmax; k++) v0[k] = 0.0;
    for (int k=0; k<ndimmax; k++) a0[k] = 0.0;
    for (int k=0; k<ndimmax; k++) agrav[k] = 0.0;
    m = 0;
    h = 0;
    invh = 0.0;
    hfactor = 0.0;
  } 

};


#endif
