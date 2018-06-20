//=================================================================================================
//  TestMhdRiemann.cpp
//  Unit tests for MHD HLLD Riemann solver
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
//  Class MhdRiemannTest
/// \brief  Unit test class for MHD HLLD Riemann solver.
/// \author S. Ganguly
/// \date   20/06/2018
//=================================================================================================
class MhdRiemannTest : public testing::Test
{
public:

  void SetUp();
  void TearDown();

};



//=================================================================================================
//  MhdRiemannTest::Setup
//=================================================================================================
void MhdRiemannTest::SetUp()
{
}



//=================================================================================================
//  MhdRiemannTest::TearDown
//=================================================================================================
void MhdRiemannTest::TearDown()
{
}



TEST_F(MhdRiemannTest, DummyTest)
{
}
