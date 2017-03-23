//=================================================================================================
//  BondiAccretionIc.cpp
//  Class for generating initial conditions for ...
//
//  This file is part of GANDALF :
//  Graphical Astrophysics code for N-body Dynamics And Lagrangian Fluids
//  https://github.com/gandalfcode/gandalf
//  Contact : gandalfcode@gmail.com
//
//  Copyright (C) 2013  D. A. Hubber, G. Rosotti
//
//  GANDALF is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 2 of the License, or
//  (at your option) any later version.
//
//  GANDALF is distributed in the hope that it will be useful, but
//  WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  General Public License (http://www.gnu.org/licenses) for more details.
//=================================================================================================


#include <fstream>
#include <sstream>
#include "Precision.h"
#include "Debug.h"
#include "Ic.h"
using namespace std;



//=================================================================================================
//  BondiAccretionIc::BondiAccretionIc
/// Set-up SILCC-type simulation initial conditions.
//=================================================================================================
template <int ndim>
BondiAccretionIc<ndim>::BondiAccretionIc(Simulation<ndim>* _sim, FLOAT _invndim) :
  Ic<ndim>(_sim, _invndim)
{
}



//=================================================================================================
//  BondiAccretionIc::ComputeBondiSolution
/// Compute the numerical solution to the Bondi accretion problem.
/// Translated from F90 subroutine written by A. P. Whitworth.
//=================================================================================================
template <int ndim>
void BondiAccretionIc<ndim>::ComputeBondiSolution
 (int Ntable,                          ///< ..
  FLOAT *w,                            ///< ..
  FLOAT *x,                            ///< ..
  FLOAT *y,                            ///< ..
  FLOAT *z)                            ///< ..
{
  const int isonic = 2*Ntable/3;       // Position in table for sonic point
  const FLOAT asympt = exp(1.5);       // ..
  FLOAT epsilon = 0.001;               // ..

  int i;                               // ..
  FLOAT disc,dx,hdx;                   // ..
  FLOAT f1,f2,f3,f4,f5;                // ..
  FLOAT g1,g2,g3,g4,g5;                // ..
  //FLOAT w0,x0,y0,z0;                   // ..
  FLOAT x1,x2,x4;                      // ..
  // FLOAT x3;
  FLOAT w1,w2,w3,w4;                   // ..
  FLOAT z1;                            // ..
  FLOAT *dlnx,*lw,*lx,*ly,*lz;         // ..

  debug2("[Ic::ComputeBondiSolution]");

  dlnx = new FLOAT[Ntable];
  lw   = new FLOAT[Ntable];
  lx   = new FLOAT[Ntable];
  ly   = new FLOAT[Ntable];
  lz   = new FLOAT[Ntable];

  // Set values for the sonic point directly
  x[isonic]  = (FLOAT) 1.0;
  lx[isonic] = (FLOAT) 0.0;
  w[isonic]  = (FLOAT) 1.0;
  lw[isonic] = (FLOAT) 0.0;
  y[isonic]  = asympt;
  ly[isonic] = log10(y[isonic]);
  z[isonic]  = (FLOAT) 2.4102434440;
  lz[isonic] = log10(z[isonic]);


  f5 = (FLOAT) 0.0;
  g5 = (FLOAT) 0.0;
  x1 = (FLOAT) 1.0 + epsilon;
  w1 = (FLOAT) 1.0 - epsilon;
  z1 = (FLOAT) 2.4102434440 + asympt*epsilon;
  i = isonic + 1;

  //-----------------------------------------------------------------------------------------------
  do {
    disc = (FLOAT) 100.0*log10(x1) + (FLOAT) isonic - (FLOAT) i;
    if (disc > (FLOAT) 0.0) {
      lx[i]   = (FLOAT) 0.01*(FLOAT)(i - isonic);
      x[i]    = powf((FLOAT) 10.0,lx[i]);
      dlnx[i] = (x[i] - x1)/x1;
      w[i]    = w1 + f5*(x[i] - x1);
      y[i]    = asympt/(x[i]*x[i]*w[i]);
      z[i]    = z1 + g5*(x[i] - x1);
      lw[i]   = log10(w[i]);
      ly[i]   = log10(y[i]);
      lz[i]   = log10(z[i]);
      i++;
    }
    dx = x1*epsilon;
    hdx = (FLOAT) 0.5*dx;
    f1 = (FLOAT) 2.0*(((FLOAT) 1.0/x1) - ((FLOAT) 1.0/(x1*x1)))/(w1 - ((FLOAT) 1.0/w1));
    g1 = asympt/w1;
    x2 = x1 + hdx;
    w2 = w1 + f1*hdx;
    f2 = (FLOAT) 2.0*(((FLOAT) 1.0/x2) - ((FLOAT) 1.0/(x2*x2)))/(w2 - ((FLOAT) 1.0/w2));
    g2 = asympt/w2;
    w3 = w1 + f2*hdx;
    f3 = (FLOAT) 2.0*(((FLOAT) 1.0/x2) - ((FLOAT) 1.0/(x2*x2)))/(w3 - ((FLOAT) 1.0/w3));
    g3 = asympt/w3;
    x4 = x1 + dx;
    w4 = w1 + f3*dx;
    f4 = (FLOAT) 2.0*(((FLOAT) 1.0/x4) - ((FLOAT) 1.0/(x4*x4)))/(w4 - ((FLOAT) 1.0/w4));
    g4 = asympt/w4;
    x1 = x4;
    f5 = (f1 + (FLOAT) 2.0*f2 + (FLOAT) 2.0*f3 + f4)/(FLOAT) 6.0;
    w1 = w1 + f5*dx;
    g5 = (g1 + (FLOAT) 2.0*g2 + (FLOAT) 2.0*g3 + g4)/(FLOAT) 6.0;
    z1 = z1 + g5*dx;
  } while (i < Ntable);


  epsilon = -epsilon;
  x1 = (FLOAT) 1.0 + epsilon;
  w1 = (FLOAT) 1.0 - epsilon;
  z1 = (FLOAT) 2.4102434440 + asympt*epsilon;
  i = isonic - 1;

  do {
    disc = (FLOAT) i - (FLOAT) 100.0*log10(x1) - (FLOAT) isonic;
    if (disc > (FLOAT) 0.0) {
      lx[i]   = (FLOAT) 0.01*(FLOAT) (i - isonic);
      x[i]    = powf((FLOAT) 10.0,lx[i]);
      dlnx[i] = (x[i] - x1)/x1;
      w[i]    = w1 + f5*(x[i] - x1);
      lw[i]   = log10(w[i]);
      y[i]    = asympt/(x[i]*x[i]*w[i]);
      ly[i]   = log10(y[i]);
      z[i]    = z1 + g5*(x[i] - x1);
      lz[i]   = -(FLOAT) 8.0;
      if (z[i] > (FLOAT) 0.1E-07) lz[i] = log10(z[i]);
      i--;
    }
    dx = x1*epsilon;
    hdx = (FLOAT) 0.5*dx;
    f1 = (FLOAT) 2.0*(((FLOAT) 1.0/x1) - ((FLOAT) 1.0/(x1*x1)))/(w1 - ((FLOAT) 1.0/w1));
    g1 = asympt/w1;
    x2 = x1 + hdx;
    w2 = w1 + f1*hdx;
    f2 = (FLOAT) 2.0*(((FLOAT) 1.0/x2) - ((FLOAT) 1.0/(x2*x2)))/(w2 - ((FLOAT) 1.0/w2));
    g2 = asympt/w2;
    w3 = w1 + f2*hdx;
    f3 = (FLOAT) 2.0*(((FLOAT) 1.0/x2) - ((FLOAT) 1.0/(x2*x2)))/(w3 - ((FLOAT) 1.0/w3));
    g3 = asympt/w3;
    x4 = x1 + dx;
    w4 = w1 + f3*dx;
    f4 = (FLOAT) 2.0*(((FLOAT) 1.0/x4) - ((FLOAT) 1.0/(x4*x4)))/(w4 - ((FLOAT) 1.0/w4));
    g4 = asympt/w4;
    x1 = x4;
    f5 = (f1 + (FLOAT) 2.0*f2 + (FLOAT) 2.0*f3 + f4)/(FLOAT) 6.0;
    w1 = w1 + f5*dx;
    g5 = (g1 + (FLOAT) 2.0*g2 + (FLOAT) 2.0*g3 + g4)/(FLOAT) 6.0;
    z1 = z1 + g5*dx;
  } while (i >= 0);


  // Write to file (for checking)
  ofstream outfile;
  string filename = "BONDI_SOLUTION.dat";
  outfile.open(filename.c_str());
  for (i=0; i<Ntable; i++) {
    outfile << x[i] << "    " << z[i] << "    " << y[i] << "    " << w[i] << endl;
  }
  outfile.close();


  // Free up all locally allocated arrays
  delete[] lz;
  delete[] ly;
  delete[] lx;
  delete[] lw;
  delete[] dlnx;

  return;
}



//=================================================================================================
//  BondiAccretionIc::Generate
/// ...
//=================================================================================================
template <int ndim>
void BondiAccretionIc<ndim>::Generate(void)
{
  if (ndim == 3) {
    int i;                            // Particle counter
    int k;                            // Dimension counter
    int Nsphere;                      // Actual number of particles in sphere
    FLOAT dr[ndim];                   // Relative position vector
    FLOAT drsqd;                      // Distance squared
    FLOAT mp;                         // Mass of one particle
    FLOAT racc;                       // Bondi accretion radius
    FLOAT rcentre[ndim];              // Position of sphere centre
    FLOAT rsonic;                     // Sonic radius
    FLOAT *r;                         // Positions of all particles
    FLOAT *v;                         // Velocities of all particles
    FLOAT *w,*x,*y,*z;                // Arrays for numerical Bondi solution

    Nbody<ndim>* nbody = sim->nbody;
    Sinks<ndim>* sinks = sim->sinks;

    // Create local copies of initial conditions parameters
    int Npart     = simparams->intparams["Nhydro"];
    FLOAT temp0   = simparams->floatparams["temp0"];
    FLOAT mu_bar  = simparams->floatparams["mu_bar"];
    FLOAT mcloud  = simparams->floatparams["mcloud"];
    FLOAT msink   = simparams->floatparams["m1"];
    FLOAT rsink   = simparams->floatparams["sink_radius"];
    FLOAT rhogas  = simparams->floatparams["rhofluid1"];
    string particle_dist = simparams->stringparams["particle_distribution"];

    debug2("[Ic::BondiAccretion]");

    // Convert any parameters to code units
    temp0  /= simunits.temp.outscale;
    mcloud /= simunits.m.outscale;
    msink  /= simunits.m.outscale;
    FLOAT radius = (FLOAT) 1.0;

    FLOAT asound = sqrtf(temp0/mu_bar);
    racc = (FLOAT) 2.0*msink/asound/asound;
    rsonic = (FLOAT) 0.5*msink/asound/asound;
    FLOAT dmdt = exp((FLOAT) 1.5)*pi*msink*msink*rhogas/powf(asound,3);

    r = new FLOAT[ndim*Npart];
    v = new FLOAT[ndim*Npart];

    // Add a sphere of random particles with origin 'rcentre' and radius 'radius'
    for (k=0; k<ndim; k++) rcentre[k] = (FLOAT) 0.0;

    // Create the sphere depending on the choice of initial particle distribution
    if (particle_dist == "random") {
      Ic<ndim>::AddRandomSphere(Npart, rcentre, radius, r, sim->randnumb);
    }
    else if (particle_dist == "cubic_lattice" || particle_dist == "hexagonal_lattice") {
      Nsphere = Ic<ndim>::AddLatticeSphere(Npart, rcentre, radius, particle_dist, r, sim->randnumb);
      if (Nsphere != Npart)
        cout << "Warning! Unable to converge to required "
             << "no. of ptcls due to lattice symmetry" << endl;
      Npart = Nsphere;
    }
    else {
      string message = "Invalid particle distribution option";
      ExceptionHandler::getIstance().raise(message);
    }

    // Allocate local and main particle memory
    hydro->Nhydro = Npart;
    sim->AllocateParticleMemory();
    //mp = mcloud / (FLOAT) Npart;
    //mp = 4.0*pi*powf(rsonic,3)*rhogas*mcloud / (FLOAT) Npart;

    // Find Bondi solution
    const int Ntablemax = 4000;
    w = new FLOAT[Ntablemax];
    x = new FLOAT[Ntablemax];
    y = new FLOAT[Ntablemax];
    z = new FLOAT[Ntablemax];
    ComputeBondiSolution(Ntablemax, w, x, y, z);

    // First check that mass of cloud is not greater than that provided by the table
    if (mcloud > z[Ntablemax-1]) {
      ExceptionHandler::getIstance().raise("Ic::BondiAccretion problem; Cloud mass too big");
    }

    // Find radius of cloud
    int iradius;
    for (i=0; i<Ntablemax; i++) {
      iradius = i;
      if (z[i] >= mcloud) break;
    }
    FLOAT remainder = (mcloud - z[iradius-1])/(z[iradius] - z[iradius - 1]);
    radius = x[iradius - 1] + remainder*(x[iradius] - x[iradius - 1]);
    cout << "Cloud mass : " << mcloud        << "    iradius   : " << iradius << endl;
    cout << "Radius     : " << radius*rsonic << "    remainder : " << remainder << endl;
    cout << "rsonic     : " << rsonic        << "    racc      : " << racc << endl;
    cout << "temp0      : " << temp0         << "    asound    : " << asound << endl;
    cout << "dmdt       : " << dmdt << endl;


    // Now stretch particles in radial direction and rescale masses to reproduce the required
    // density distribution provided by the numiercal solution.
    //-----------------------------------------------------------------------------------------------
    //mp = 4.0_PR*pi*pow(rsonic,3)*rhogas*mcloud / (FLOAT) Npart;
    //mp = 0.5_PR*PI*rhogas*mcloud*msink**3 / (agas**6) / real(pgas,PR)
    //mp = mcloud/(FLOAT) Npart;
    mp = 4.0*pi*powf(rsonic,3)*rhogas*mcloud / (FLOAT) Npart;
    cout << "mp : " << mp << "    "
         << 0.5*pi*rhogas*mcloud*powf(msink,3)/powf(asound,6)/(FLOAT) Npart << endl;


    for (i=0; i<Npart; i++) {
      Particle<ndim>& part = hydro->GetParticlePointer(i);
      for (k=0; k<ndim; k++) dr[k] = r[ndim*i + k] - rcentre[k];
      FLOAT drmag = sqrtf(DotProduct(dr, dr, ndim)) + small_number;
      FLOAT mint = mcloud*powf(drmag,3);
      for (int j=0; j<Ntablemax; j++) {
        iradius = j;
        if (z[j] >= mint) break;
      }
      remainder = (mint - z[iradius-1])/(z[iradius] - z[iradius - 1]);
      FLOAT radp = x[iradius - 1] + remainder*(x[iradius] - x[iradius - 1]);
      FLOAT vradp = w[iradius - 1] + remainder*(w[iradius] - w[iradius - 1]);
      for (k=0; k<ndim; k++) r[ndim*i + k] = rsonic*r[ndim*i + k]*radp/drmag;
      for (k=0; k<ndim; k++) part.r[k] = r[ndim*i + k];
      for (k=0; k<ndim; k++) part.v[k] = -asound*vradp*dr[k]/drmag;
      part.m = mp;
    }


    // Now add star/sink particle
    rsink *= rsonic;
    nbody->Nstar = 1;
    sinks->Nsink = 1;
    for (k=0; k<ndim; k++) nbody->stardata[0].r[k] = 0.0;
    for (k=0; k<ndim; k++) nbody->stardata[0].v[k] = 0.0;
    nbody->stardata[0].m      = msink;
    nbody->stardata[0].radius = rsink;
    nbody->stardata[0].h      = nbody->kernp->invkernrange*nbody->stardata[0].radius;
    nbody->stardata[0].invh   = 1.0/nbody->stardata[0].h;
    sinks->sink[0].star        = &(nbody->stardata[0]);
    sinks->sink[0].radius      = rsink;
    sinks->sink[0].racc        = rsink;
    sinks->sink[0].mmax        = 0.0;
    sinks->sink[0].menc        = 0.0;


    // Find total mass inside sink and set to mmax
    for (i=0; i<Npart; i++) {
      Particle<ndim>& part = hydro->GetParticlePointer(i);
      for (k=0; k<ndim; k++) dr[k] = part.r[k] - rcentre[k];
      drsqd = DotProduct(dr, dr, ndim);
      if (drsqd > rsink*rsink) continue;
      sinks->sink[0].mmax += part.m;
    }

    cout << "mmax : " << sinks->sink[0].mmax << "    " << sinks->sink[0].mmax/mp << endl;
    cout << "rsink : " << rsink << "     rsink/rsonic : " << rsink/rsonic << endl;

    sim->initial_h_provided = false;

    delete[] v;
    delete[] r;
  }

  return;
}



template class BondiAccretionIc<1>;
template class BondiAccretionIc<2>;
template class BondiAccretionIc<3>;
