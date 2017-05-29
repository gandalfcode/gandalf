//=================================================================================================
//  TestTree.cpp
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
#include "Parameters.h"
#include "RandomNumber.h"
#include "Sph.h"
#include "KDTree.h"
#include "gtest/gtest.h"


//=================================================================================================
//  class TreeTest
//=================================================================================================
class TreeTest : public testing::Test
{
public:

  void SetUp(void);
  void TearDown(void);

  static const int Npart = 1000;
  static const unsigned long rseed = 1000;
  static const int Nleafmax = 6;
  static const FLOAT thetamaxsqd = 0.1;
  static const FLOAT kernrange = 2.0;
  static const FLOAT macerror = 0.0001;
  string gravity_mac = "geometric";
  string multipole = "monopole";

  SphParticle<3> *partdata;
  Parameters params;
  KDTree<3,SphParticle,KDTreeCell> *kdtree;
  RandomNumber *randnumb;
  DomainBox<3> simbox;

};



//=================================================================================================
//  TreeTest::Setup
//=================================================================================================
void TreeTest::SetUp(void)
{
  // Create all objects for test
  partdata = new SphParticle<3>[Npart];
  randnumb = new XorshiftRand(rseed);
  kdtree = new KDTree<3,SphParticle,KDTreeCell>(Nleafmax, thetamaxsqd, kernrange,
                                                macerror, gravity_mac, multipole);

  // Create some simple particle configuration (random population of a cube)
  for (int i=0; i<Npart; i++) {
    for (int k=0; k<3; k++) partdata[i].r[k] = 1.0 - 2.0*randnumb->floatrand();
    partdata[i].m = 1.0 / (FLOAT) Npart;
    partdata[i].h = 0.2*randnumb->floatrand();
  }
  for (int k=0; k<3; k++) {
    simbox.boxmin[k] = -1.0;
    simbox.boxmax[k] = 1.0;
    simbox.boxsize[k] = 2.0;
    simbox.boxhalf[k] = 1.0;
  }

  // Now build the tree using the particle configuration
  kdtree->Ntot       = Npart;
  kdtree->Ntotmaxold = 0;
  kdtree->Ntotmax    = Npart;
  kdtree->ifirst     = 0;
  kdtree->ilast      = Npart - 1;
  kdtree->BuildTree(0, Npart-1, Npart, Npart, partdata, 0.0);
  kdtree->StockTree(kdtree->celldata[0], partdata);

  return;
}



//=================================================================================================
//  TreeTest::TearDown
//=================================================================================================
void TreeTest::TearDown(void)
{
  // Delete all objects created for test
  delete kdtree;
  delete randnumb;
  delete[] partdata;

  return;
}



//=================================================================================================
//  StructureTest
//=================================================================================================
TEST_F(TreeTest, StructureTest)
{
  bool overlap_flag;                   // ..
  int c;                               // Cell counter
  int cc;                              // Aux. cell counter
  int i;                               // Particle counter
  int j;                               // Aux. particle counter
  int k;                               // ..
  int l;                               // Tree level
  int activecount;                     // ..
  int Nactivecount=0;                  // Counter for total no. of active ptcls
  int leafcount=0;                     // ..
  int Ncount=0;                        // Total particle counter
  int *ccount;                         // Array for counting cells
  int *lcount;                         // Array for counting ptcls on each level
  int *pcount;                         // Array for counting particles in tree
  KDTreeCell<3> cell;                  // Local copy of tree cell

  ccount = new int[kdtree->Ncellmax];
  pcount = new int[Npart];
  lcount = new int[kdtree->ltot+1];

  for (i=0; i<Npart; i++) pcount[i] = 0;
  for (c=0; c<kdtree->Ncellmax; c++) ccount[c] = 0;
  for (l=0; l<kdtree->ltot; l++) lcount[l] = 0;
  Ncount = 0;
  Nactivecount = 0;

  // Count how many times we enter a cell in a full tree walk
  c = 0;
  while (c < kdtree->Ncell) {
    ccount[c]++;
    if (kdtree->celldata[c].c1 != -1) c = kdtree->celldata[c].c1;
    else c = kdtree->celldata[c].cnext;
  }

  // Now check we enter all cells once and once only
  for (c=0; c<kdtree->Ncell; c++) {
    ASSERT_EQ(ccount[c], 1);
  }

  // Verify linked lists are valid for all levels of tree
  //-----------------------------------------------------------------------------------------------
  for (l=0; l<=kdtree->ltot; l++) {
    for (i=0; i<Npart; i++) pcount[i] = 0;

    for (c=0; c<kdtree->Ncell; c++) {
      cell = kdtree->celldata[c];

      // Check that particles are not in linked lists more than once
      if (cell.level == l) {
        i = cell.ifirst;
        while (i != -1) {
          pcount[i]++;
          if (i == cell.ilast) break;
          i = kdtree->inext[i];
        }

      }
    }

    // Check particles are included in the tree once and once only
    for (i=0; i<Npart; i++) {
      ASSERT_EQ(pcount[i], 1);
    }

  }
  for (i=0; i<Npart; i++) pcount[i] = 0;


  // Loop over all cells in tree
  //-----------------------------------------------------------------------------------------------
  for (c=0; c<kdtree->Ncell; c++) {
    cell = kdtree->celldata[c];
    activecount = 0;
    leafcount = 0;

    // Add particles from current level
    lcount[cell.level] += cell.N;

    // Check that particles are not in linked lists more than once
    if (cell.level == kdtree->ltot) {
      i = cell.ifirst;
      while (i != -1) {
        pcount[i]++;
        leafcount++;
        Ncount++;
        if (partdata[i].flags.check(active)) activecount++;
        if (partdata[i].flags.check(active)) Nactivecount++;
        ASSERT_LE(partdata[i].h, cell.hmax);
        for (k=0; k<3; k++) ASSERT_LE(partdata[i].r[k], cell.bbmax[k]);
        for (k=0; k<3; k++) ASSERT_GE(partdata[i].r[k], cell.bbmin[k]);
        if (i == cell.ilast) break;
        i = kdtree->inext[i];
      }
      //ASSERT_LE(leafcount,Nleafmax);
      ASSERT_LE(activecount, leafcount);
      for (k=0; k<3; k++) ASSERT_LE(cell.rcell[k], cell.bbmax[k]);
      for (k=0; k<3; k++) ASSERT_GE(cell.rcell[k], cell.bbmin[k]);
    }

    // Check that bounding boxes of cells on each level do not overlap each other
    for (cc=0; cc<kdtree->Ncell; cc++) {
      overlap_flag = false;
      if (c != cc && kdtree->celldata[cc].level == cell.level) {
        if (cell.bbmin[0] < kdtree->celldata[cc].bbmax[0] &&
            cell.bbmax[0] > kdtree->celldata[cc].bbmin[0] &&
            cell.bbmin[1] < kdtree->celldata[cc].bbmax[1] &&
            cell.bbmax[1] > kdtree->celldata[cc].bbmin[1] &&
            cell.bbmin[2] < kdtree->celldata[cc].bbmax[2] &&
            cell.bbmax[2] > kdtree->celldata[cc].bbmin[2]) {
          overlap_flag = true;
        }
      }
      ASSERT_EQ(false, overlap_flag);
    }
  }
  //-----------------------------------------------------------------------------------------------


  delete[] lcount;
  delete[] pcount;
  delete[] ccount;

}



//=================================================================================================
//  ComTest
//=================================================================================================
TEST_F(TreeTest, ComTest)
{
  FLOAT mtot = 0.0;
  FLOAT rcom[3],vcom[3];

  for (int k=0; k<3; k++) rcom[k] = 0.0;
  for (int k=0; k<3; k++) vcom[k] = 0.0;

  for (int i=0; i<Npart; i++) {
    for (int k=0; k<3; k++) rcom[k] += partdata[i].m*partdata[i].r[k];
    for (int k=0; k<3; k++) vcom[k] += partdata[i].m*partdata[i].v[k];
    mtot += partdata[i].m;
  }
  for (int k=0; k<3; k++) rcom[k] /= mtot;
  for (int k=0; k<3; k++) vcom[k] /= mtot;

  EXPECT_FLOAT_EQ(mtot, kdtree->celldata[0].m);
  EXPECT_FLOAT_EQ(rcom[0], kdtree->celldata[0].r[0]);
  EXPECT_FLOAT_EQ(rcom[1], kdtree->celldata[0].r[1]);
  EXPECT_FLOAT_EQ(rcom[2], kdtree->celldata[0].r[2]);
  EXPECT_FLOAT_EQ(vcom[0], kdtree->celldata[0].v[0]);
  EXPECT_FLOAT_EQ(vcom[1], kdtree->celldata[0].v[1]);
  EXPECT_FLOAT_EQ(vcom[2], kdtree->celldata[0].v[2]);
}



//=================================================================================================
//  GatherTest
//=================================================================================================
TEST_F(TreeTest, GatherTest)
{
  int i;
  int j;
  int k;
  int Nneib;
  int Nneibtree = 0;
  int *neiblist;
  FLOAT dr[3];
  FLOAT drsqd;
  FLOAT hrange;

  neiblist = new int[Npart];

  for (i=0; i<Npart; i++) {
    Nneib = 0;
    hrange = 2.0*partdata[i].h;
    for (j=0; j<Npart; j++) {
      for (k=0; k<3; k++) dr[k] = partdata[i].r[k] - partdata[j].r[k];
      drsqd = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
      if (drsqd < hrange*hrange) Nneib++;
    }
    Nneibtree = 0;
    Nneibtree = kdtree->ComputeGatherNeighbourList
      (partdata, partdata[i].r, hrange, Npart, Nneibtree, neiblist);

    ASSERT_EQ(Nneib, Nneibtree);
  }

}



//=================================================================================================
//  PeriodicWalkTest
//=================================================================================================
TEST_F(TreeTest, PeriodicWalkTest)
{
  /*for (int cc=0; cc<kdtree->Ncell; cc++) {
    bool overlap_flag = false;
    if (c != cc && kdtree->celldata[cc].level == cell.level) {
      if (cell.bbmin[0] < kdtree->celldata[cc].bbmax[0] &&
          cell.bbmax[0] > kdtree->celldata[cc].bbmin[0] &&
          cell.bbmin[1] < kdtree->celldata[cc].bbmax[1] &&
          cell.bbmax[1] > kdtree->celldata[cc].bbmin[1] &&
          cell.bbmin[2] < kdtree->celldata[cc].bbmax[2] &&
          cell.bbmax[2] > kdtree->celldata[cc].bbmin[2]) {
          overlap_flag = true;
      }
    }
    ASSERT_EQ(false, overlap_flag);
  }
  bool okflag = kdtree->ComputePeriodicGravityInteractionList
    (cell, sphdata, simbox, macfactor, Nneibmax, Ngravcellmax, Nneib, Nhydroneib,
     Ndirect, Ngravcell, neiblist, sphlist, directlist, gravcell, neibpart);*/
}
