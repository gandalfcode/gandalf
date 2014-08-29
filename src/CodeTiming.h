//=============================================================================
//  CodeTiming.h
//  Contains class for controlling internal timing routines.
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
//=============================================================================


#ifndef _CODE_TIMING_H_
#define _CODE_TIMING_H_


#include <map>
#include <string>
#include <time.h>
#include <iostream>
#include <iomanip>
#include <sys/time.h>
#include "Precision.h"
using namespace std;



//=============================================================================
//  Class CodeTiming
/// \brief   Internal timing control class
/// \details Internal timing control class
/// \author  D. A. Hubber
/// \date    12/02/2014
//=============================================================================
class CodeTiming
{
 public:

  // Constructor and destructor
  //---------------------------------------------------------------------------
  CodeTiming();
  ~CodeTiming();
  
  
  // Other functions
  //---------------------------------------------------------------------------
  void StartTimingSection(string,int);
  void EndTimingSection(string);
  void ComputeTimingStatistics(string);
  double WallClockTime(void);
  
  
  // CodeTiming class variables
  //---------------------------------------------------------------------------
  static const int Nblockmax=256;   ///< Max. no. of code timing blocks
  int Nblock;                       ///< No. of timing blocks
  int Nlevelmax;                    ///< Max. no. of timing levels
  double tstart_wall;               ///< ..
  double tend_wall;                 ///< ..
  DOUBLE ttot;                      ///< ..
  DOUBLE ttot_wall;                 ///< ..
  clock_t tstart;                   ///< ..
  clock_t tend;                     ///< ..
  map<string,int> blockno;          ///< ..


  struct TimingBlock {
    bool timing_flag;
    int timing_level;
    int Ncalled;
    double tstart_wall;
    double tend_wall;
    double ttot_wall;
    double tfraction_wall;
    DOUBLE ttot;
    DOUBLE tfraction;
    clock_t tstart;
    clock_t tend;
    string block_name;
    
    TimingBlock()
    {
      timing_flag = false;
      timing_level = 0;
      Ncalled = 0;
      ttot = 0.0;
      ttot_wall = 0.0;
      tfraction = 0.0;
      tfraction_wall = 0.0;
      block_name = "";
    }
    
  };

  TimingBlock block[Nblockmax];
    
};
#endif
