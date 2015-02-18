//=================================================================================================
//  TestScaling.cpp
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
#include "SimUnits.h"
#include "gtest/gtest.h"


//=================================================================================================
//  class ScalingTest
//=================================================================================================
class ScalingTest : public testing::Test
{
public:

  void SetUp(void);
  void TearDown(void);

  SimUnits units;
  Parameters params;

};



//=================================================================================================
//  ScalingTest::Setup
//=================================================================================================
void ScalingTest::SetUp(void)
{
  units.SetupUnits(&params);

  return;
}



//=================================================================================================
//  ScalingTest::TearDown
//=================================================================================================
void ScalingTest::TearDown(void)
{
  return;
}



TEST_F(ScalingTest, rTest)
{
  EXPECT_DOUBLE_EQ(units.r.inscale*units.r.inSI,units.r.outscale*units.r.outSI);
  EXPECT_DOUBLE_EQ(units.r.OutputScale("pc"),units.r.outscale);
  EXPECT_DOUBLE_EQ(units.r.OutputScale("m"),units.r.outscale*units.r.outSI);
  EXPECT_DOUBLE_EQ(units.r.OutputScale("cm"),units.r.outscale*units.r.outcgs);
  EXPECT_DOUBLE_EQ(units.r.OutputScale("au"),r_pc/r_au);
}



TEST_F(ScalingTest, mTest)
{
  EXPECT_DOUBLE_EQ(units.m.inscale*units.m.inSI,units.m.outscale*units.m.outSI);
  EXPECT_DOUBLE_EQ(units.m.OutputScale("m_sun"),units.m.outscale);
  EXPECT_DOUBLE_EQ(units.m.OutputScale("kg"),units.m.outscale*units.m.outSI);
  EXPECT_DOUBLE_EQ(units.m.OutputScale("g"),units.m.outscale*units.m.outcgs);
  EXPECT_DOUBLE_EQ(units.m.OutputScale("m_earth"),m_sun/m_earth);
}



TEST_F(ScalingTest, tTest)
{
  EXPECT_DOUBLE_EQ(units.t.inscale*units.t.inSI,units.t.outscale*units.t.outSI);
  EXPECT_DOUBLE_EQ(units.t.OutputScale("myr"),units.t.outscale);
  EXPECT_DOUBLE_EQ(units.t.OutputScale("s"),units.t.outscale*units.t.outSI);
  EXPECT_DOUBLE_EQ(units.t.OutputScale("s"),units.t.outscale*units.t.outcgs);
  EXPECT_DOUBLE_EQ(units.t.OutputScale("day")/units.t.OutputScale("myr"),myr/day);
}



TEST_F(ScalingTest, vTest)
{
  EXPECT_DOUBLE_EQ(units.v.inscale*units.v.inSI,units.v.outscale*units.v.outSI);
  EXPECT_DOUBLE_EQ(units.v.OutputScale("km_s"),units.v.outscale);
  EXPECT_DOUBLE_EQ(units.v.OutputScale("m_s"),units.v.outscale*units.v.outSI);
  EXPECT_DOUBLE_EQ(units.v.OutputScale("cm_s"),units.v.outscale*units.v.outcgs);
  EXPECT_DOUBLE_EQ(units.v.outscale*units.v.outSI,units.r.outscale*units.r.outSI/
                   (units.t.outscale*units.t.outSI));
}



TEST_F(ScalingTest, rhoTest)
{
  EXPECT_DOUBLE_EQ(units.rho.inscale*units.rho.inSI,units.rho.outscale*units.rho.outSI);
  EXPECT_DOUBLE_EQ(units.rho.OutputScale("g_cm3"),units.rho.outscale);
  EXPECT_DOUBLE_EQ(units.rho.OutputScale("kg_m3"),units.rho.outscale*units.rho.outSI);
  EXPECT_DOUBLE_EQ(units.rho.OutputScale("g_cm3"),units.rho.outscale*units.rho.outcgs);
  EXPECT_DOUBLE_EQ(units.rho.outscale*units.rho.outSI,units.m.outscale*units.m.outSI/
                   pow(units.r.outscale*units.r.outSI,3));
}



TEST_F(ScalingTest, ETest)
{
  EXPECT_DOUBLE_EQ(units.E.inscale*units.E.inSI,units.E.outscale*units.E.outSI);
  EXPECT_DOUBLE_EQ(units.E.OutputScale("J"),units.E.outscale);
  EXPECT_DOUBLE_EQ(units.E.OutputScale("J"),units.E.outscale*units.E.outSI);
  EXPECT_DOUBLE_EQ(units.E.OutputScale("erg"),units.E.outscale*units.E.outcgs);
  EXPECT_DOUBLE_EQ(units.E.outscale*units.E.outSI,units.m.outscale*units.m.outSI*
                   pow(units.r.outscale*units.r.outSI,2)/pow(units.t.outscale*units.t.outSI,2));
}



TEST_F(ScalingTest, uTest)
{
  EXPECT_DOUBLE_EQ(units.u.inscale*units.u.inSI,units.u.outscale*units.u.outSI);
  EXPECT_DOUBLE_EQ(units.u.OutputScale("J_kg"),units.u.outscale);
  EXPECT_DOUBLE_EQ(units.u.OutputScale("J_kg"),units.u.outscale*units.u.outSI);
  EXPECT_DOUBLE_EQ(units.u.OutputScale("erg_g"),units.u.outscale*units.u.outcgs);
  EXPECT_DOUBLE_EQ(units.u.outscale*units.u.outSI,pow(units.r.outscale*units.r.outSI,2)/
                   pow(units.t.outscale*units.t.outSI,2));
}



TEST_F(ScalingTest, tempTest)
{
  EXPECT_DOUBLE_EQ(units.temp.inscale*units.temp.inSI,units.temp.outscale*units.temp.outSI);
  EXPECT_DOUBLE_EQ(units.temp.OutputScale("K"),units.temp.outscale);
  EXPECT_DOUBLE_EQ(units.temp.OutputScale("K"),units.temp.outscale*units.temp.outSI);
  EXPECT_DOUBLE_EQ(units.temp.OutputScale("K"),units.temp.outscale*units.temp.outcgs);
}
