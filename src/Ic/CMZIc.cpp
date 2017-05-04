//=================================================================================================
//  CMZIc.cpp
//  Class for generating initial conditions for turbulent cloud on CMZ orbit
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


#include "Precision.h"
#include "Debug.h"
#include "Ic.h"



//=================================================================================================
//  CMZIc::CMZIc
/// Constructor
//=================================================================================================
template <int ndim>
CMZIc<ndim>::CMZIc(Simulation<ndim>* _sim, FLOAT _invndim) :
  Ic<ndim>(_sim, _invndim)
{
}



//=================================================================================================
//  CMZIc::Generate
/// Set-up turbulent cloud on CMZ orbit by reading in cloud from SPHNG dump
//=================================================================================================
template <int ndim>
void CMZIc<ndim>::Generate(void)
{
  int i;                               // Particle counter
  int k;                               // Dimension counter
  int Nsphere;                         // Actual number of particles in sphere
  FLOAT gpecloud;                      // Total grav. potential energy of entire cloud
  FLOAT keturb;                        // Total turbulent kinetic energy of entire cloud
  FLOAT mp;                            // Mass of one particle
  FLOAT rcentre[ndim];                 // Position of sphere centre
  FLOAT rmax;                          // Current max size of cloud
  FLOAT rp;                            // Particle radial position
  FLOAT rscale;                        // Length scaling factor
  FLOAT vfactor;                       // Velocity scale factor
  FLOAT *rx;                           // x Positions of all particles
  FLOAT *ry;                           // y Positions of all particles
  FLOAT *rz;                           // z Positions of all particles
  FLOAT *vx;                           // x Velocities of all particles
  FLOAT *vy;                           // y Velocities of all particles
  FLOAT *vz;                           // z Velocities of all particles
  FLOAT *h;                            // smoothing lengths of all particles
  FLOAT *m;                            // masses of all particles
  FLOAT sigma3d;                       // 3D velocity dispersion
  FLOAT kewant;                        // desired kinetic energy
  FLOAT vxmean;                        // RMS x velocity
  FLOAT vymean;                        // RMS y velocity
  FLOAT vzmean;                        // RMS z velocity
  FLOAT vmean;                         // Mean of the RMS components
  FLOAT vfactorx;                      // correction factor in x dirn
  FLOAT vfactory;                      // correction factor in y dirn
  FLOAT vfactorz;                      // correction factor in z dirn

  // Create local copies of initial conditions parameters
  int field_type   = simparams->intparams["field_type"];
  int gridsize     = simparams->intparams["gridsize"];
  int Npart        = simparams->intparams["Nhydro"];
  FLOAT gammaone   = simparams->floatparams["gamma_eos"] - 1.0;
  FLOAT alpha_turb = simparams->floatparams["alpha_turb"];
  FLOAT mcloud     = simparams->floatparams["mcloud"];
  FLOAT mu_bar     = simparams->floatparams["mu_bar"];
  FLOAT radius     = simparams->floatparams["radius"];
  FLOAT temp0      = simparams->floatparams["temp0"];
  FLOAT cloudposx  = simparams->floatparams["cloudposx"];
  FLOAT cloudposy  = simparams->floatparams["cloudposy"];
  FLOAT cloudposz  = simparams->floatparams["cloudposz"];
  FLOAT cloudvelx  = simparams->floatparams["cloudvelx"];
  FLOAT cloudvely  = simparams->floatparams["cloudvely"];
  FLOAT cloudvelz  = simparams->floatparams["cloudvelz"];

  debug2("[CMZIc::Generate]");

  // Convert any parameters to code units
  mcloud /= simunits.m.outscale;
  radius /= simunits.r.outscale;
  temp0  /= simunits.temp.outscale;
  cloudposx /= simunits.r.outscale;
  cloudposy /= simunits.r.outscale;
  cloudposz /= simunits.r.outscale;
  cloudvelx /= simunits.v.outscale;
  cloudvely /= simunits.v.outscale;
  cloudvelz /= simunits.v.outscale;

  rx = new FLOAT[Npart];
  ry = new FLOAT[Npart];
  rz = new FLOAT[Npart];
  vx = new FLOAT[Npart];
  vy = new FLOAT[Npart];
  vz = new FLOAT[Npart];
  h = new FLOAT[Npart];
  m = new FLOAT[Npart];

  // Read in particle properties

  cout << "Using SPH_NG ics" << endl;

  ifstream inFile("sanity.out");
  if (inFile.fail()){
       cerr << "Unable to open CMZ SPHNG IC file." << endl;
  }
  cout <<"[CMZ_SPHNG] opening SPHNG IC file"<<endl;

  int il=0;
  for (il=0;il<Npart;il++) {

    inFile >> rx[il] >> ry[il] >> rz[il] >> m[il] >> h[il] >> vx[il] >> vy[il] >> vz[il];

    }

  inFile.close();

  cout <<"[CMZ_SPHNG] closing SPHNG IC file"<<endl;

  // Allocate local and main particle memory
  hydro->Nhydro = Npart;
  sim->AllocateParticleMemory();

  // Determine current cloud radius

  rmax=-1.00;
  for (i=0; i<hydro->Nhydro; i++) {
      rp=(rx[i]*rx[i]+ry[i]*ry[i]+rz[i]*rz[i]);
      rmax=max(rp,rmax);

  }

  rmax=pow(rmax,0.5);

  cout << "current rmax: " << rmax << endl;

  rscale=radius/rmax;

  cout << "length scaling factor: " << rscale << endl;

  // Determine required particle mass

  mp = mcloud / (FLOAT) Npart;

  // Estimate GPE

  gpecloud = (FLOAT) 1.75*mcloud*mcloud/radius;

  sigma3d = sqrt(0.6*alpha_turb*mcloud/radius);

  kewant=0.5*mcloud*sigma3d*sigma3d;

  // Calculate total kinetic energy of turbulent velocity field
  // and compare means in x, y, z direction
  keturb = (FLOAT) 0.0;
  vxmean = (FLOAT) 0.0;
  vymean = (FLOAT) 0.0;
  vzmean = (FLOAT) 0.0;

  for (i=0; i<hydro->Nhydro; i++) {
    keturb += mp*(vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i]);
    vxmean += (vx[i]*vx[i]);
    vymean += (vy[i]*vy[i]);
    vzmean += (vz[i]*vz[i]);
  }
  keturb *= (FLOAT) 0.5;
  vxmean = sqrt(vxmean/(FLOAT) Npart);
  vymean = sqrt(vymean/(FLOAT) Npart);
  vzmean = sqrt(vzmean/(FLOAT) Npart);

  vmean=(vxmean+vymean+vzmean)/3.;
  vfactorx=vmean/vxmean;
  vfactory=vmean/vymean;
  vfactorz=vmean/vzmean;

  // force velocity field to be isotropic

  keturb = (FLOAT) 0.0;
  for (i=0; i<hydro->Nhydro; i++) {
    vx[i] *= vfactorx;
    vy[i] *= vfactory;
    vz[i] *= vfactorz;
    keturb += mp*(vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i]);
  }

  keturb *= (FLOAT) 0.5;

  vfactor = sqrt(kewant/keturb);
  cout << "kewant,keturb : " << keturb <<" "<< kewant << endl;

  // Record particle properties in main memory
  for (i=0; i<hydro->Nhydro; i++) {
      Particle<ndim>& part = hydro->GetParticlePointer(i);
      part.r[0] = rx[i]*rscale;
      part.r[1] = ry[i]*rscale;
      part.r[2] = rz[i]*rscale;
      part.v[0] = vx[i]*vfactor;
      part.v[1] = vy[i]*vfactor;
      part.v[2] = vz[i]*vfactor;
      part.m = mp;
      part.h = h[i]*rscale;
      part.u = temp0/gammaone/mu_bar;
      part.ptype = gas_type;
  }

  keturb = (FLOAT) 0.0;
  for (i=0; i<hydro->Nhydro; i++) {
      Particle<ndim>& part = hydro->GetParticlePointer(i);
      keturb += mp*(part.v[0]*part.v[0]+part.v[1]*part.v[1]+part.v[2]*part.v[2]);
  }
  keturb *= (FLOAT) 0.5;

  vfactor = sqrt(kewant/keturb);
  cout << "kewant,keturb : " << keturb <<" "<< kewant << endl;


  // Change to COM frame of reference
  sim->SetComFrame();

  sim->initial_h_provided = true;

  vxmean = (FLOAT) 0.0;
  vymean = (FLOAT) 0.0;
  vzmean = (FLOAT) 0.0;

  for (i=0; i<hydro->Nhydro; i++) {
    Particle<ndim>& part = hydro->GetParticlePointer(i);
    vxmean += (part.v[0]*part.v[0]);
    vymean += (part.v[1]*part.v[1]);
    vzmean += (part.v[2]*part.v[2]);
  }

  vxmean = sqrt(vxmean/(FLOAT) Npart);
  vymean = sqrt(vymean/(FLOAT) Npart);
  vzmean = sqrt(vzmean/(FLOAT) Npart);

  cout << "mean vx, vy, vz: " << vxmean << " " << vymean << " " << vzmean << endl;

  // Finally, move cloud to required position and bulk velocity
  for (i=0; i<hydro->Nhydro; i++) {
      Particle<ndim>& part = hydro->GetParticlePointer(i);
      part.r[0]+=cloudposx;
      part.r[1]+=cloudposy;
      part.r[2]+=cloudposz;
      part.v[0]+=cloudvelx;
      part.v[1]+=cloudvely;
      part.v[2]+=cloudvelz;
  }

  delete[] vx;
  delete[] vy;
  delete[] vz;
  delete[] rx;
  delete[] ry;
  delete[] rz;
  delete[] h;
  delete[] m;

  return;
}



template class CMZIc<1>;
template class CMZIc<2>;
template class CMZIc<3>;
