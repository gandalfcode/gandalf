//=================================================================================================
//  TestSinks.cpp
//  ..
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


#include <math.h>
#include "Constants.h"
#include "Hydrodynamics.h"
#include "Nbody.h"
#include "Parameters.h"
#include "RandomNumber.h"
#include "Sinks.h"
#include "Sph.h"
#include "SphNeighbourSearch.h"
#include "gtest/gtest.h"


//=================================================================================================
//  class SinkTest
//=================================================================================================
class SinkTest : public testing::Test
{
public:

  void SetUp(void);
  void TearDown(void);
  void CalculateCom(FLOAT &, FLOAT *, FLOAT *, FLOAT *, Sph<3> *, Nbody<3> *);

  const int hydro_forces    = 1;
  const int nbody_softening = 1;
  const int self_gravity    = 1;
  const int sub_systems     = 0;
  const FLOAT alpha_visc    = 1.0;
  const FLOAT beta_visc     = 1.0;
  const FLOAT h_fac         = 1.2;
  const FLOAT h_converge    = 0.01;
  const FLOAT nbody_mult    = 0.1;
  const aviscenum avisc     = mon97;
  const acondenum acond     = noac;
  const tdaviscenum tdavisc = notdav;
  const string gas_eos      = "energy_eqn";
  const string kernelName   = "m4";

  static const unsigned long rseed = 1000;
  static const int Npart = 1024;

  CodeTiming timing;
  DomainBox<3> simbox;
  GradhSph<3,M4Kernel> *sph;
  NbodyLeapfrogKDK<3,M4Kernel> *nbody;
  StarParticle<3> *star;
  Sinks<3> sinks;
  Parameters params;
  RandomNumber *randnumb;
  GradhSphBruteForce<3,GradhSphParticle> *sphneib;
  SimUnits simunits;

};



//=================================================================================================
//  SinkTest::Setup
//=================================================================================================
void SinkTest::SetUp(void)
{
  // Create all objects for test
  sph = new GradhSph<3,M4Kernel>(hydro_forces, self_gravity, alpha_visc, beta_visc, h_fac,
                                 h_converge, avisc, acond, tdavisc, gas_eos, kernelName);
  nbody = new NbodyLeapfrogKDK<3,M4Kernel>(nbody_softening, sub_systems, nbody_mult, kernelName);
  sphneib = new GradhSphBruteForce<3,GradhSphParticle>(sph->kernrange, &simbox, sph->kernp, &timing);
  randnumb = new XorshiftRand(rseed);
  sph->eos = new Isothermal<3>(params.floatparams["temp0"], params.floatparams["mu_bar"],
                               params.floatparams["gamma_eos"], &simunits);

  // Allocate all memory
  sph->AllocateMemory(Npart);
  nbody->AllocateMemory(Npart);
  sinks.AllocateMemory(Npart);
  sph->Nhydro = Npart;
  sph->Ntot = Npart;

  // Set important sink parameters
  sinks.create_sinks = 1;
  sph->create_sinks = 1;
  sinks.smooth_accretion = 0;
  sinks.timing = &timing;


  // Create some simple particle configuration (random population of a cube)
  for (int i=0; i<Npart; i++) {
    SphParticle<3> &part = sph->GetSphParticlePointer(i);
    for (int k=0; k<3; k++) {
      part.r[k]  = 1.0 - 2.0*randnumb->floatrand();
      part.r0[k] = part.r[k];
      part.a[k]  = 0.0;
      part.a0[k] = 0.0;
    }
    part.m      = 1.0 / (FLOAT) Npart;
    part.h      = 0.1;
    part.flags.set(active);
    part.ptype  = gas_type;
    part.flags.unset(potmin);
    part.gpot   = 0.0;
    part.nstep  = 1;
  }

  for (int k=0; k<3; k++) {
    simbox.boxmin[k]  = -1.0;
    simbox.boxmax[k]  = 1.0;
    simbox.boxsize[k] = 2.0;
    simbox.boxhalf[k] = 1.0;
  }

  sph->InitialSmoothingLengthGuess();

  // Calculate smoothing lengths of all particles
  sph->InitialSmoothingLengthGuess();
  sphneib->UpdateAllSphProperties(Npart,Npart,sph->GetSphParticleArray(),sph,nbody);
  sphneib->UpdateAllSphProperties(Npart,Npart,sph->GetSphParticleArray(),sph,nbody);

  // Calculate forces of all particles
  for (int i=0; i<Npart; i++) {
    SphParticle<3> &part = sph->GetSphParticlePointer(i);
    part.gpot = (FLOAT) 0.0;
    for (int k=0; k<3; k++) part.a[k] = (FLOAT) 0.0;
  }
  sphneib->UpdateAllSphForces(Npart,Npart,sph->GetSphParticleArray(),sph,nbody);

  // Calculate smoothing lengths of all particles
  sphneib->UpdateAllSphProperties(Npart,Npart,sph->GetSphParticleArray(),sph,nbody);

  // Calculate forces of all particles
  for (int i=0; i<Npart; i++) {
    SphParticle<3> &part = sph->GetSphParticlePointer(i);
    part.gpot = (FLOAT) 0.0;
    for (int k=0; k<3; k++) part.a[k] = (FLOAT) 0.0;
  }
  sphneib->UpdateAllSphForces(Npart,Npart,sph->GetSphParticleArray(),sph,nbody);

  return;
}



//=================================================================================================
//  SinkTest::TearDown
//=================================================================================================
void SinkTest::TearDown(void)
{
  // Delete all objects created for test
  delete randnumb;
  delete nbody;
  delete sph;

  return;
}



//=================================================================================================
//  SinkTest::CalculateCom
//=================================================================================================
void SinkTest::CalculateCom(FLOAT &mtot, FLOAT rcom[3], FLOAT vcom[3], FLOAT acom[3],
                            Sph<3> *sph, Nbody<3> *nbody)
{
  mtot = 0.0;
  for (int k=0; k<3; k++) {
    rcom[k] = 0.0;
    vcom[k] = 0.0;
    acom[k] = 0.0;
  }

  for (int i=0; i<sph->Nhydro; i++) {
    SphParticle<3> &part = sph->GetSphParticlePointer(i);
    if (part.flags.is_dead()) continue;
    for (int k=0; k<3; k++) {
      rcom[k] += part.m*part.r[k];
      vcom[k] += part.m*part.v[k];
      acom[k] += part.m*part.a[k];
    }
    mtot += part.m;
  }

  for (int i=0; i<nbody->Nnbody; i++) {
    NbodyParticle<3> *star = nbody->nbodydata[i];
    for (int k=0; k<3; k++) {
      rcom[k] += star->m*star->r[k];
      vcom[k] += star->m*star->v[k];
      acom[k] += star->m*star->a[k];
    }
    mtot += star->m;
  }

  for (int k=0; k<3; k++) {
    rcom[k] /= mtot;
    vcom[k] /= mtot;
    acom[k] /= mtot;
  }

  cout << "mtot : " << mtot << endl;
  cout << "acom : " << acom[0] << "    " << acom[1] << "    " << acom[2] << endl;

  return;
}



//=================================================================================================
//  CreateNewSinkTest
//=================================================================================================
TEST_F(SinkTest, CreateNewSinkTest)
{
  FLOAT mtot1, mtot2;
  FLOAT acom1[3],rcom1[3],vcom1[3];
  FLOAT acom2[3],rcom2[3],vcom2[3];

  // Calculate COM before sink is created
  CalculateCom(mtot1, rcom1, vcom1, acom1, sph, nbody);

  // Now create new sink particle
  sinks.rho_sink = 0.0;
  sinks.SearchForNewSinkParticles(0.0, 0.0, sph, nbody);

  cout << "Nsink : " << sinks.Nsink << endl;

  // Calculate COM after sink is created
  CalculateCom(mtot2, rcom2, vcom2, acom2, sph, nbody);

  EXPECT_FLOAT_EQ(mtot1,mtot2);
  EXPECT_FLOAT_EQ(rcom1[0],rcom2[0]);
  EXPECT_FLOAT_EQ(rcom1[1],rcom2[1]);
  EXPECT_FLOAT_EQ(rcom1[2],rcom2[2]);
  EXPECT_FLOAT_EQ(vcom1[0],vcom2[0]);
  EXPECT_FLOAT_EQ(vcom1[1],vcom2[1]);
  EXPECT_FLOAT_EQ(vcom1[2],vcom2[2]);
  EXPECT_FLOAT_EQ(acom1[0],acom2[0]);
  EXPECT_FLOAT_EQ(acom1[1],acom2[1]);
  EXPECT_FLOAT_EQ(acom1[2],acom2[2]);


  // Accrete mass to sinks
  sinks.AccreteMassToSinks(sph, nbody, 0, 0.0);
  cout << "Nhydro : " << sph->Nhydro << " particle(s)" << endl;
  sph->DeleteDeadParticles();

  // Calculate COM after sink accretes
  CalculateCom(mtot2, rcom2, vcom2, acom2, sph, nbody);

  EXPECT_FLOAT_EQ(mtot1,mtot2);
  EXPECT_FLOAT_EQ(rcom1[0],rcom2[0]);
  EXPECT_FLOAT_EQ(rcom1[1],rcom2[1]);
  EXPECT_FLOAT_EQ(rcom1[2],rcom2[2]);
  EXPECT_FLOAT_EQ(vcom1[0],vcom2[0]);
  EXPECT_FLOAT_EQ(vcom1[1],vcom2[1]);
  EXPECT_FLOAT_EQ(vcom1[2],vcom2[2]);
  EXPECT_FLOAT_EQ(acom1[0],acom2[0]);
  EXPECT_FLOAT_EQ(acom1[1],acom2[1]);
  EXPECT_FLOAT_EQ(acom1[2],acom2[2]);

}
